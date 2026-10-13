// Command-line front end for the nested BFB distance.
//
//   nbfb -i INSTANCE.txt [--history]
//       One source per line, the target on the last line, in the format of
//       bfb_dist's read_instance. Prints the distance on the last line of
//       stdout, as an integer or "inf".
//
//   nbfb --sources SOURCES.txt --targets TARGETS.txt [-j N] [-o OUT.tsv]
//       One shared source set and one target per line. Prints one TSV row
//       per target.
//
//   --foldback-first [--foldback-span BP|any]   (experimental, either mode)
//       Minimize first the number of t-BFB cycles whose join is
//       fold-back-shaped: the segment before the join and the first segment
//       of the appended source lie on the same chromosome in opposite
//       orientations, with breakpoints at most BP apart (default 50000).
//       Adds a foldback_t line / column; such cycles print as 'T'.
//   --foldback-total   (experimental, implies --foldback-first)
//       Also report the objective (fold-back-shaped t-BFB cycles, total
//       cycles, t-BFB cycles): fb_total, fb_total_t and fb_total_types.
//   --definition prefix|lcp   (either mode)
//       Use the direct search over strings of nbfb_state.hpp instead of
//       the index-tuple search. "prefix" is the nested class of
//       2026-10-02 (every string is t[1,psi] y, y empty, an inverted window
//       of t, or a source); "lcp" is the earlier LCP-monotone class that
//       the index-tuple search solves, for differential testing.
//   --blocks   (either mode)
//       The sources are all contiguous substrings of the given source
//       lines and their reverse complements (nbfb_blocks.hpp): the
//       index-tuple search for the prefix definition of nestedness.

#include "nbfb_fast.hpp"
#include "nbfb_state.hpp"
#include "nbfb_blocks.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
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

// Set by --definition; null selects the index-tuple search.
const nbfb::Definition* state_definition = nullptr;
// Set by --blocks: the sources are all substrings of the given lines.
bool use_blocks = false;

nbfb::Result solve(const nbfb::PreparedSources& sources,
                   const std::vector<int>& target, bool want_history,
                   nbfb::MemoryBudget& budget, const nbfb::JoinTable* joins) {
    if (use_blocks) {
        return nbfb::nbfb_blocks_dist(sources, target, want_history, budget,
                                      joins);
    }
    if (state_definition != nullptr) {
        return nbfb::nbfb_state_dist(sources, target, want_history, budget,
                                     *state_definition, joins);
    }
    return nbfb::nbfb_fast_dist(sources, target, want_history, budget, joins);
}

int run_single(const std::string& path, bool want_history,
               nbfb::MemoryBudget& budget, bool foldback_first,
               long long foldback_span, bool foldback_total) {
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

    nbfb::PreparedSources prepared(std::move(sources));
    nbfb::JoinTable joins;
    if (foldback_first) joins = interner.join_table(foldback_span);
    const nbfb::Result result = solve(
        prepared, target, want_history, budget,
        foldback_first ? &joins : nullptr);
    if (foldback_first) {
        std::cout << "foldback_t: "
                  << (result.reachable ? std::to_string(result.foldback_t)
                                       : std::string("inf"))
                  << "\n";
    }
    if (foldback_total) {
        if (result.reachable) {
            std::cout << "fb_total: " << result.fb_total << "\n";
            std::cout << "fb_total_t: " << result.fb_total_t << "\n";
            if (want_history) {
                std::cout << "fb_total_types: "
                          << types_of(result.fb_total_history) << "\n";
                std::cout << "fb_total_operations: "
                          << ops_of(result.fb_total_history) << "\n";
            }
        } else {
            std::cout << "fb_total: inf\nfb_total_t: inf\n";
            if (want_history) {
                std::cout << "fb_total_types: \nfb_total_operations: \n";
            }
        }
    }
    if (want_history) {
        if (result.reachable) {
            std::cout << "t_cycles: " << result.t_cycles << "\n";
            std::cout << "types: " << types_of(result.history) << "\n";
            std::cout << "operations: " << ops_of(result.history) << "\n";
        } else {
            std::cout << "t_cycles: inf\n";
            std::cout << "types: \n";
            std::cout << "operations: \n";
        }
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
    int fb_total = 0;
    int fb_total_t = 0;
    std::string fb_total_types;
    std::string fb_total_ops;
    double seconds = 0.0;
    std::string error;
};

int run_batch(const std::string& sources_path, const std::string& targets_path,
              int threads, const std::string& out_path, bool want_history,
              nbfb::MemoryBudget& budget, bool foldback_first,
              long long foldback_span, bool foldback_total) {
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

    // The interner is the only table shared across targets that depends on
    // the sources alone, together with the prepared sources themselves; a
    // target introducing a new segment extends the interner, so the
    // targets are interned before the worker threads start.
    std::vector<std::vector<int>> targets;
    targets.reserve(target_lines.size());
    for (const std::string& line : target_lines) {
        targets.push_back(interner.sequence(line));
    }
    const nbfb::PreparedSources prepared(std::move(sources));
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
                    solve(prepared, targets[i], true, budget, joins_ptr);
                row.ok = true;
                row.reachable = result.reachable;
                row.foldback_t = result.foldback_t;
                row.t_cycles = result.t_cycles;
                row.total = result.total;
                row.types = types_of(result.history);
                if (want_history) row.ops = ops_of(result.history);
                row.fb_total = result.fb_total;
                row.fb_total_t = result.fb_total_t;
                row.fb_total_types = types_of(result.fb_total_history);
                if (want_history) {
                    row.fb_total_ops = ops_of(result.fb_total_history);
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
         << "t_cycles\ttotal\ttypes\tseconds";
    if (want_history) *out << "\tops";
    if (foldback_total) {
        *out << "\tfb_total\tfb_total_t\tfb_total_types";
        if (want_history) *out << "\tfb_total_ops";
    }
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
        if (want_history) *out << '\t' << row.ops;
        if (foldback_total) {
            if (!row.ok) *out << "\terror\terror\t";
            else if (!row.reachable) *out << "\tinf\tinf\t";
            else *out << '\t' << row.fb_total << '\t' << row.fb_total_t
                      << '\t' << row.fb_total_types;
            if (want_history) *out << '\t' << row.fb_total_ops;
        }
        *out << "\n";
    }
    out->flush();
    return failures == 0 ? 0 : 1;
}

void usage() {
    std::cerr
        << "usage: nbfb -i INSTANCE.txt [--history] [--max-mem-gb G]\n"
           "       nbfb --sources SOURCES.txt --targets TARGETS.txt "
           "[-j N] [-o OUT.tsv] [--history] [--max-mem-gb G]\n"
           "       either form: [--foldback-first [--foldback-span BP|any]] "
           "[--foldback-total]\n"
           "       either form: [--definition prefix|lcp] [--blocks]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string instance;
    std::string sources;
    std::string targets;
    std::string out;
    int threads = 1;
    bool want_history = false;
    double max_mem_gb = 2.0;
    bool foldback_first = false;
    bool foldback_total = false;
    long long foldback_span = 50000;
    nbfb::Definition definition = nbfb::Definition::Prefix;

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
        } else if (argument == "--history") {
            want_history = true;
        } else if (argument == "--max-mem-gb") {
            max_mem_gb = std::atof(value("--max-mem-gb").c_str());
        } else if (argument == "--foldback-first") {
            foldback_first = true;
        } else if (argument == "--foldback-total") {
            foldback_first = true;
            foldback_total = true;
        } else if (argument == "--foldback-span") {
            const std::string span = value("--foldback-span");
            foldback_span = span == "any" ? -1 : std::atoll(span.c_str());
        } else if (argument == "--definition") {
            const std::string name = value("--definition");
            if (name == "prefix") {
                definition = nbfb::Definition::Prefix;
            } else if (name == "lcp") {
                definition = nbfb::Definition::Lcp;
            } else {
                std::cerr << "Error: --definition must be prefix or lcp\n";
                return 2;
            }
            state_definition = &definition;
        } else if (argument == "--blocks") {
            use_blocks = true;
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
    if (use_blocks && (state_definition != nullptr || foldback_total)) {
        std::cerr << "Error: --blocks cannot be combined with --definition "
                     "or --foldback-total\n";
        return 2;
    }
    if (state_definition != nullptr && foldback_total) {
        std::cerr << "Error: --foldback-total needs the index-tuple search\n";
        return 2;
    }
    nbfb::MemoryBudget budget(
        static_cast<std::size_t>(max_mem_gb * 1024.0 * 1024.0 * 1024.0));

    try {
        if (!instance.empty()) {
            if (!sources.empty() || !targets.empty()) {
                std::cerr << "Error: -i cannot be combined with --sources or "
                             "--targets\n";
                return 2;
            }
            return run_single(instance, want_history, budget, foldback_first,
                              foldback_span, foldback_total);
        }
        if (sources.empty() || targets.empty()) {
            usage();
            return 2;
        }
        return run_batch(sources, targets, threads, out, want_history, budget,
                         foldback_first, foldback_span, foldback_total);
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
