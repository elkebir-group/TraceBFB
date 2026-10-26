// Command-line front end for the nested BFB distance (nbfb_simple.hpp).
//
//   nbfb -i INSTANCE.txt
//       One source per line, the target on the last line, in the format of
//       bfbdist's read_instance. Prints foldback_t (with --foldback-first),
//       t_cycles, types, operations and optimal_histories lines, then the
//       distance on the last line, as an integer or "inf".
//
//   nbfb --sources SOURCES.txt --targets TARGETS.txt [-j N] [-o OUT.tsv]
//       One shared source set and one target per line. Prints one TSV row
//       per target, with one optimal history (types, ops) and the number of
//       optimal histories (optimal_histories).
//
//   Nestedness is the foldable-prefix definition: every string of the
//   history is t[1,psi] y with y or the reverse complement of y a
//   substring of t.
//
//   --blocks (default) | --general   (either mode)
//       The sources are every contiguous block of the source lines and the
//       reverse complement of every such block (a file listing every block
//       denotes the same set), or the source lines exactly as listed.
//   --foldback-first (default) | --no-foldback-first   (either mode)
//       Minimize first the number of t-BFB cycles whose join is
//       fold-back-shaped (the segment before the join and the first segment
//       of the appended source lie on the same chromosome in opposite
//       orientations, with breakpoints at most --foldback-span BP apart,
//       default 50000, or 'any'), then t-BFB cycles, then all cycles. Such
//       cycles print as 'T'. --no-foldback-first minimizes t-BFB cycles,
//       then all cycles.
//   --all-optimal | --max-histories N   (either mode)
//       The number of optimal histories (same strings and cycle types = same
//       history) is always reported. These options also list the optimal
//       histories: all of them, or at most N. Single mode: one "optimal:"
//       line per history (types, tab, operations) before the distance.
//       Batch mode: an all_ops column with the histories separated by '|'.

#include "nbfb.hpp"
#include "nbfb_simple.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <unordered_map>
#include <vector>

namespace {

// -------------------------------------------------------------- parsing

// Maps every distinct oriented segment to a signed identifier. Two
// segments get identifiers of the same magnitude exactly when one is the
// reverse complement of the other, so reverse complementing a sequence is
// reversing the order and negating every identifier.
class SegmentInterner {
  public:
    int id(const std::string& token) {
        std::string chromosome;
        long long start = 0;
        long long end = 0;
        bool complement = false;
        parse(token, chromosome, start, end, complement);
        const std::string forward = canonical(chromosome, start, end, complement);
        auto found = table_.find(forward);
        if (found != table_.end()) return found->second;
        const std::string reverse =
            canonical(chromosome, end, start, !complement);
        auto other = table_.find(reverse);
        if (other != table_.end()) {
            const int value = -other->second;
            table_.emplace(forward, value);
            return value;
        }
        const int value = ++next_;
        table_.emplace(forward, value);
        table_.emplace(reverse, -value);
        geometry_.push_back({chromosome, start, end, complement});
        return value;
    }

    // Fold-back-shaped joins u -> v between all interned segments: same
    // chromosome, opposite orientations, and the coordinate where u is left
    // within span bp of the coordinate where v is entered (span < 0: any).
    nbfb::JoinTable join_table(long long span) const {
        nbfb::JoinTable joins;
        joins.max_id = next_;
        const std::size_t width = 2 * static_cast<std::size_t>(next_) + 1;
        joins.foldback.assign(width * width, 0);
        for (int u = -next_; u <= next_; ++u) {
            if (u == 0) continue;
            const Geometry left = oriented(u);
            for (int v = -next_; v <= next_; ++v) {
                if (v == 0) continue;
                const Geometry right = oriented(v);
                if (left.chromosome != right.chromosome ||
                    left.complement == right.complement) {
                    continue;
                }
                const long long gap = std::llabs(left.end - right.start);
                if (span >= 0 && gap > span) continue;
                joins.foldback[static_cast<std::size_t>(u + next_) * width +
                               static_cast<std::size_t>(v + next_)] = 1;
            }
        }
        return joins;
    }

    std::vector<int> sequence(const std::string& line) {
        std::vector<int> result;
        std::size_t begin = 0;
        while (true) {
            const std::size_t stop = line.find(';', begin);
            const std::string token = line.substr(
                begin, stop == std::string::npos ? std::string::npos
                                                 : stop - begin);
            result.push_back(id(token));
            if (stop == std::string::npos) break;
            begin = stop + 1;
        }
        return result;
    }

  private:
    static std::string trim(const std::string& text) {
        std::size_t begin = 0;
        std::size_t stop = text.size();
        while (begin < stop && std::isspace(static_cast<unsigned char>(text[begin]))) ++begin;
        while (stop > begin && std::isspace(static_cast<unsigned char>(text[stop - 1]))) --stop;
        return text.substr(begin, stop - begin);
    }

    static void parse(const std::string& raw, std::string& chromosome,
                      long long& start, long long& end, bool& complement) {
        const std::string text = trim(raw);
        const std::size_t colon = text.find(':');
        const std::size_t open = text.find('[');
        if (colon == std::string::npos || open == std::string::npos ||
            open + 3 != text.size() || text[open + 2] != ']' ||
            (text[open + 1] != 'R' && text[open + 1] != 'C')) {
            throw std::invalid_argument("invalid segment '" + text + "'");
        }
        const std::size_t dash = text.find('-', colon + 1);
        if (dash == std::string::npos || dash > open) {
            throw std::invalid_argument("invalid segment '" + text + "'");
        }
        chromosome = text.substr(0, colon);
        if (chromosome.empty() ||
            chromosome.find_first_of(";: \t") != std::string::npos) {
            throw std::invalid_argument("invalid segment '" + text + "'");
        }
        start = digits(text, colon + 1, dash, text);
        end = digits(text, dash + 1, open, text);
        complement = text[open + 1] == 'C';
    }

    static long long digits(const std::string& text, std::size_t begin,
                            std::size_t stop, const std::string& whole) {
        if (begin >= stop) {
            throw std::invalid_argument("invalid segment '" + whole + "'");
        }
        long long value = 0;
        for (std::size_t i = begin; i < stop; ++i) {
            if (text[i] < '0' || text[i] > '9') {
                throw std::invalid_argument("invalid segment '" + whole + "'");
            }
            value = value * 10 + (text[i] - '0');
        }
        return value;
    }

    static std::string canonical(const std::string& chromosome, long long start,
                                 long long end, bool complement) {
        return chromosome + ':' + std::to_string(start) + '-' +
               std::to_string(end) + (complement ? "[C]" : "[R]");
    }

    struct Geometry {
        std::string chromosome;
        long long start;
        long long end;
        bool complement;
    };

    Geometry oriented(int id) const {
        Geometry g = geometry_[static_cast<std::size_t>(std::abs(id)) - 1];
        if (id < 0) {
            std::swap(g.start, g.end);
            g.complement = !g.complement;
        }
        return g;
    }

    std::unordered_map<std::string, int> table_;
    std::vector<Geometry> geometry_;  // orientation of each positive id
    int next_ = 0;

};

std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot read '" + path + "'");
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

bool blank(const std::string& line) {
    for (char c : line) {
        if (!std::isspace(static_cast<unsigned char>(c))) return false;
    }
    return true;
}

// -------------------------------------------------------------- output

std::string types_of(const std::vector<nbfb::Operation>& history) {
    std::string text;
    text.reserve(history.size());
    for (const nbfb::Operation& operation : history) {
        text.push_back(operation.foldback ? 'T' : operation.kind);
    }
    return text;
}

std::string ops_of(const std::vector<nbfb::Operation>& history) {
    std::string text;
    for (const nbfb::Operation& operation : history) {
        if (!text.empty()) text.push_back(';');
        if (operation.kind == 't') {
            text += "t,";
            text += std::to_string(operation.source);
        } else {
            text.push_back(operation.kind);
            text.push_back(',');
            text += std::to_string(operation.x);
            text.push_back(',');
            text += std::to_string(operation.y);
        }
    }
    return text;
}

// -------------------------------------------------------------- modes

// Set by --blocks / --general and --all-optimal / --max-histories.
bool blocks = true;
std::size_t max_histories = 0;

nbfb::Result solve(const std::vector<std::vector<int>>& sources,
                   const std::vector<int>& target,
                   nbfb::MemoryBudget& budget, const nbfb::JoinTable* joins) {
    return nbfb::nbfb_simple_dist(sources, target, budget, joins, blocks,
                                  max_histories);
}

int run_single(const std::string& path,
               nbfb::MemoryBudget& budget, bool foldback_first,
               long long foldback_span) {
    const std::vector<std::string> lines = read_lines(path);
    bool bad = lines.size() < 2;
    for (const std::string& line : lines) {
        if (blank(line)) bad = true;
    }
    if (bad) {
        std::cerr << "Error: input file must contain one or more sources "
                     "followed by the target, one nonempty sequence per line\n";
        return 2;
    }
    SegmentInterner interner;
    std::vector<std::vector<int>> sources;
    for (std::size_t i = 0; i + 1 < lines.size(); ++i) {
        sources.push_back(interner.sequence(lines[i]));
    }
    const std::vector<int> target = interner.sequence(lines.back());

    nbfb::JoinTable joins;
    if (foldback_first) joins = interner.join_table(foldback_span);
    const nbfb::Result result = solve(
        sources, target, budget, foldback_first ? &joins : nullptr);
    if (foldback_first) {
        std::cout << "foldback_t: "
                  << (result.reachable ? std::to_string(result.foldback_t)
                                       : std::string("inf"))
                  << "\n";
    }
    if (result.reachable) {
        std::cout << "t_cycles: " << result.t_cycles << "\n";
        std::cout << "types: " << types_of(result.history) << "\n";
        std::cout << "operations: " << ops_of(result.history) << "\n";
    } else {
        std::cout << "t_cycles: inf\n";
        std::cout << "types: \n";
        std::cout << "operations: \n";
    }
    std::cout << "optimal_histories: "
              << (result.reachable ? result.optimal_count : std::string("0"))
              << "\n";
    for (const std::vector<nbfb::Operation>& history : result.optimal_histories) {
        std::cout << "optimal: " << types_of(history) << '\t' << ops_of(history)
                  << "\n";
    }
    if (result.reachable) {
        std::cout << result.total << std::endl;
    } else {
        std::cout << "inf" << std::endl;
    }
    return 0;
}

struct Row {
    std::size_t index = 0;
    std::size_t tokens = 0;
    bool ok = false;
    bool reachable = false;
    int foldback_t = 0;
    int t_cycles = 0;
    int total = 0;
    std::string types;
    std::string ops;
    std::string optimal_count;
    std::string all_ops;
    double seconds = 0.0;
    std::string error;
};

int run_batch(const std::string& sources_path, const std::string& targets_path,
              int threads, const std::string& out_path,
              nbfb::MemoryBudget& budget, bool foldback_first,
              long long foldback_span) {
    SegmentInterner interner;
    std::vector<std::vector<int>> sources;
    for (const std::string& line : read_lines(sources_path)) {
        if (blank(line)) continue;
        sources.push_back(interner.sequence(line));
    }
    if (sources.empty()) {
        std::cerr << "Error: no source in '" << sources_path << "'\n";
        return 2;
    }
    std::vector<std::string> target_lines;
    for (const std::string& line : read_lines(targets_path)) {
        if (blank(line)) continue;
        target_lines.push_back(line);
    }

    // A target introducing a new segment extends the interner, so the
    // targets are interned before the worker threads start.
    std::vector<std::vector<int>> targets;
    targets.reserve(target_lines.size());
    for (const std::string& line : target_lines) {
        targets.push_back(interner.sequence(line));
    }
    nbfb::JoinTable joins;
    if (foldback_first) joins = interner.join_table(foldback_span);
    const nbfb::JoinTable* joins_ptr = foldback_first ? &joins : nullptr;

    std::vector<Row> rows(targets.size());
    std::atomic<std::size_t> next(0);
    auto worker = [&]() {
        while (true) {
            const std::size_t i = next.fetch_add(1);
            if (i >= targets.size()) return;
            Row& row = rows[i];
            row.index = i;
            row.tokens = targets[i].size();
            const auto started = std::chrono::steady_clock::now();
            try {
                const nbfb::Result result =
                    solve(sources, targets[i], budget, joins_ptr);
                row.ok = true;
                row.reachable = result.reachable;
                row.foldback_t = result.foldback_t;
                row.t_cycles = result.t_cycles;
                row.total = result.total;
                row.types = types_of(result.history);
                row.ops = ops_of(result.history);
                row.optimal_count = result.reachable ? result.optimal_count : "0";
                for (const auto& history : result.optimal_histories) {
                    if (!row.all_ops.empty()) row.all_ops.push_back('|');
                    row.all_ops += ops_of(history);
                }
            } catch (const std::exception& error) {
                row.ok = false;
                row.error = error.what();
            }
            row.seconds = std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - started)
                              .count();
        }
    };

    if (threads <= 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (int i = 0; i < threads; ++i) pool.emplace_back(worker);
        for (std::thread& item : pool) item.join();
    }

    std::ostream* out = &std::cout;
    std::ofstream file;
    if (!out_path.empty()) {
        file.open(out_path);
        if (!file) {
            std::cerr << "Error: cannot write '" << out_path << "'\n";
            return 2;
        }
        out = &file;
    }
    *out << "index\ttokens\t" << (foldback_first ? "foldback_t\t" : "")
         << "t_cycles\ttotal\ttypes\tseconds\tops\toptimal_histories";
    if (max_histories > 0) *out << "\tall_ops";
    *out << "\n";
    int failures = 0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const Row& row = rows[i];
        *out << i << '\t' << row.tokens << '\t';
        if (foldback_first) {
            if (!row.ok) *out << "error\t";
            else if (!row.reachable) *out << "inf\t";
            else *out << row.foldback_t << '\t';
        }
        if (!row.ok) {
            ++failures;
            *out << "error\terror\t" << row.error << '\t';
        } else if (!row.reachable) {
            *out << "inf\tinf\t\t";
        } else {
            *out << row.t_cycles << '\t' << row.total << '\t' << row.types
                 << '\t';
        }
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "%.4f", row.seconds);
        *out << buffer;
        *out << '\t' << row.ops;
        *out << '\t' << (row.ok ? row.optimal_count : "error");
        if (max_histories > 0) *out << '\t' << row.all_ops;
        *out << "\n";
    }
    out->flush();
    return failures == 0 ? 0 : 1;
}

void usage() {
    std::cerr
        << "usage: nbfb -i INSTANCE.txt [--max-mem-gb G]\n"
           "       nbfb --sources SOURCES.txt --targets TARGETS.txt "
           "[-j N] [-o OUT.tsv] [--max-mem-gb G]\n"
           "       either form: [--blocks | --general] "
           "[--foldback-first | --no-foldback-first] [--foldback-span BP|any] "
           "[--all-optimal | --max-histories N]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string instance;
    std::string sources;
    std::string targets;
    std::string out;
    int threads = 1;
    double max_mem_gb = 2.0;
    bool foldback_first = true;
    long long foldback_span = 50000;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << name << " needs a value\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (argument == "-i" || argument == "--input") {
            instance = value("-i");
        } else if (argument == "--sources") {
            sources = value("--sources");
        } else if (argument == "--targets") {
            targets = value("--targets");
        } else if (argument == "-o" || argument == "--out") {
            out = value("-o");
        } else if (argument == "-j" || argument == "--threads") {
            threads = std::atoi(value("-j").c_str());
        } else if (argument == "--max-mem-gb") {
            max_mem_gb = std::atof(value("--max-mem-gb").c_str());
        } else if (argument == "--foldback-first") {
            foldback_first = true;
        } else if (argument == "--no-foldback-first") {
            foldback_first = false;
        } else if (argument == "--foldback-span") {
            const std::string span = value("--foldback-span");
            if (span == "any") {
                foldback_span = -1;
            } else {
                char* end = nullptr;
                foldback_span = std::strtoll(span.c_str(), &end, 10);
                if (span.empty() || *end != '\0' || foldback_span < 0) {
                    std::cerr << "Error: --foldback-span must be a "
                                 "non-negative integer or 'any'\n";
                    return 2;
                }
            }
        } else if (argument == "--blocks") {
            blocks = true;
        } else if (argument == "--general") {
            blocks = false;
        } else if (argument == "--all-optimal") {
            max_histories = static_cast<std::size_t>(-1);
        } else if (argument == "--max-histories") {
            const std::string text = value("--max-histories");
            char* end = nullptr;
            const unsigned long long n = std::strtoull(text.c_str(), &end, 10);
            if (text.empty() || *end != '\0' || text[0] == '-' || text[0] == '+') {
                std::cerr << "Error: --max-histories must be a non-negative integer\n";
                return 2;
            }
            max_histories = static_cast<std::size_t>(n);
        } else if (argument == "-h" || argument == "--help") {
            usage();
            return 0;
        } else {
            std::cerr << "Error: unknown argument '" << argument << "'\n";
            usage();
            return 2;
        }
    }

    if (max_mem_gb <= 0.0) {
        std::cerr << "Error: --max-mem-gb must be positive\n";
        return 2;
    }
    if (threads < 1) threads = 1;
    nbfb::MemoryBudget budget(
        static_cast<std::size_t>(max_mem_gb * 1024.0 * 1024.0 * 1024.0));

    try {
        if (!instance.empty()) {
            if (!sources.empty() || !targets.empty()) {
                std::cerr << "Error: -i cannot be combined with --sources or "
                             "--targets\n";
                return 2;
            }
            return run_single(instance, budget, foldback_first,
                              foldback_span);
        }
        if (sources.empty() || targets.empty()) {
            usage();
            return 2;
        }
        return run_batch(sources, targets, threads, out, budget,
                         foldback_first, foldback_span);
    } catch (const nbfb::MemoryLimitExceeded& error) {
        std::cerr << "Error: " << error.what()
                  << " (raise --max-mem-gb to continue)\n";
        return 3;
    } catch (const std::bad_alloc&) {
        std::cerr << "Error: out of memory (the process limit was reached "
                     "before --max-mem-gb)\n";
        return 3;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 2;
    }
}
