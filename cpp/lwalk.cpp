// lwalk: exact longest walks of a breakpoint graph.
//
// A walk starts at a loose end, alternates segment traversals and bonds (reference adjacencies
// and junctions), and uses every segment at most its copy number and every bond at most its
// multiplicity.  A walk's length is its number of segment traversals, and its letters are the
// segment id when the segment is entered at its start (node 2i) and -id (in the older format the
// lower-cased name, or -name when it has no lower case) when entered at its end (node 2i+1).  GRAPH is the breakpoint graph format of docs/file_formats.md
// (start labels X.L / X.R), or the older segments/bonds/starts format.
//
// Dynamic program over residual states.  A state is (start, current exit node, remaining bond
// multiplicities); the remaining segment copy numbers follow from these, because every segment
// traversal except the first is entered through a bond and every one except the current is left
// through a bond.  For each state the program stores the maximum number of further traversals and
// the number of distinct continuations reaching it (128-bit).  Parallel bonds between the same two
// nodes are merged, so distinct continuations are distinct letter strings.
//
//   lwalk GRAPH.txt [--starts LABEL,LABEL] [--max-mem-gb 4] [--max-foldbacks]   counts only
//       --max-foldbacks: maximize the number of fold-back junction traversals first and the number
//       of traversals only as a tie-break, so that no walk wins by leaving a fold-back to another
//       molecule
//   lwalk GRAPH.txt ... --enumerate -o PREFIX [--dedup-rc] [--max-walks N] [--chunk N]  all longest walks
//   lwalk GRAPH.txt ... --sample N [--seed S] -o PREFIX [--dedup-rc]      uniform sample
//
// PREFIX.walks.txt gets one walk per line as space-separated letters, PREFIX.targets.txt the same
// walk as ';'-separated tokens (the target format of cpp/nbfb --targets).  With --dedup-rc a walk
// whose reverse complement is also an admissible walk (it ends at an allowed start) is written
// only in the lexicographically smaller of its two readings.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

typedef unsigned __int128 u128;
static const int KW = 4;  // key words: 256 bits

static std::string u128str(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v) { s.push_back(char('0' + int(v % 10))); v /= 10; }
    std::reverse(s.begin(), s.end());
    return s;
}

struct Graph {
    int S = 0, M = 0;
    std::vector<std::string> name, fwd, rev;
    std::vector<int> cap;
    std::vector<int> ba, bb, bm, bfb;
    std::vector<std::vector<std::pair<int, int>>> adj;  // node -> (bond, other node)
    std::vector<int> starts;
    std::vector<std::string> start_label;
    std::string header;
    bool signed_ids = false;  // graph format: reverse traversals print as -id
};

// The reverse of a name without coordinates: lower case, or -name when lower-casing changes
// nothing (s0, 12).
static std::string rc_name(const std::string& n) {
    std::string l = n;
    for (auto& c : l) c = char(std::tolower(c));
    return l != n ? l : "-" + n;
}

static std::string rc_token(const std::string& t) {
    // chrom:a-b[R] -> chrom:b-a[C]; a token without coordinates (a letter) as rc_name
    size_t colon = t.rfind(':');
    if (colon == std::string::npos || t.find('[', colon) == std::string::npos) return rc_name(t);
    size_t dash = t.find('-', colon), open = t.find('[', dash);
    std::string chrom = t.substr(0, colon), a = t.substr(colon + 1, dash - colon - 1),
                b = t.substr(dash + 1, open - dash - 1);
    char strand = t[open + 1] == 'R' ? 'C' : 'R';
    return chrom + ":" + b + "-" + a + "[" + strand + "]";
}

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> f;
    size_t begin = 0;
    while (true) {
        size_t stop = line.find('\t', begin);
        f.push_back(line.substr(begin, stop == std::string::npos ? std::string::npos : stop - begin));
        if (stop == std::string::npos) break;
        begin = stop + 1;
    }
    return f;
}

static Graph read_graph(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    Graph g;
    std::string line, word;
    std::map<std::pair<int, int>, int> merged, foldback;
    std::vector<std::pair<int, int>> order;
    auto add_bond = [&](int a, int b, int m, int fb) {
        auto key = std::make_pair(std::min(a, b), std::max(a, b));
        if (!merged.count(key)) order.push_back(key);
        merged[key] += m;
        foldback[key] = foldback.count(key) ? std::max(foldback[key], fb) : fb;
    };
    std::vector<std::string> lines;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    bool tsv = false;
    for (const std::string& l : lines) {
        if (l.empty() || l[0] == '#') continue;
        tsv = l.rfind("H\t", 0) == 0;
        break;
    }
    if (tsv) {
        g.signed_ids = true;
        // The breakpoint graph format of docs/file_formats.md. Segment i has node 2i at X.L
        // and 2i+1 at X.R; a free end is a start labeled with its segment end.
        std::map<std::string, int> index;
        std::vector<std::vector<std::string>> rest;
        for (const std::string& l : lines) {
            if (l.empty() || l[0] == '#') continue;
            std::vector<std::string> f = split_tabs(l);
            if (f[0] == "S") {
                if (f.size() < 7) throw std::runtime_error("bad S record: " + l);
                index[f[1]] = g.S++;
                g.name.push_back(f[1]);
                g.cap.push_back(std::stoi(f[6]));
                std::string tk = f[2] == "*" ? f[1]
                                 : f[5] == "+" ? f[2] + ":" + f[3] + "-" + f[4] + "[R]"
                                               : f[2] + ":" + f[4] + "-" + f[3] + "[C]";
                g.fwd.push_back(tk);
                g.rev.push_back(f[2] == "*" ? f[1] : rc_token(tk));
            } else {
                rest.push_back(f);
            }
        }
        auto node = [&](const std::string& end) {
            size_t dot = end.rfind('.');
            auto it = index.find(end.substr(0, dot));
            if (dot == std::string::npos || it == index.end())
                throw std::runtime_error("unknown segment end " + end);
            return 2 * it->second + (end.substr(dot + 1) == "R" ? 1 : 0);
        };
        for (const auto& f : rest) {
            const size_t need = f[0] == "R" ? 4 : f[0] == "J" ? 6 : f[0] == "F" ? 3 : 0;
            if (f.size() < need) throw std::runtime_error("bad " + f[0] + " record: too few fields");
            if (f[0] == "R") add_bond(node(f[1]), node(f[2]), std::stoi(f[3]), 0);
            else if (f[0] == "J") add_bond(node(f[2]), node(f[3]), std::stoi(f[4]), std::stoi(f[5]));
            else if (f[0] == "F") { g.starts.push_back(node(f[1])); g.start_label.push_back(f[1]); }
        }
    } else {
        // The older segments/bonds/starts format.
        size_t k = 0;
        auto next = [&](std::string& out) {
            if (k >= lines.size()) return false;
            out = lines[k++];
            return true;
        };
        while (next(line)) {
            if (line.empty()) continue;
            if (line[0] == '#') { g.header = line; continue; }
            std::istringstream ls(line);
            ls >> word;
            int n = 0;
            ls >> n;
            if (word == "segments") {
                g.S = n;
                for (int i = 0; i < n; ++i) {
                    next(line);
                    std::istringstream s(line);
                    std::string nm, tk;
                    int c;
                    s >> nm >> c >> tk;
                    g.name.push_back(nm); g.cap.push_back(c); g.fwd.push_back(tk); g.rev.push_back(rc_token(tk));
                }
            } else if (word == "bonds") {
                for (int i = 0; i < n; ++i) {
                    next(line);
                    std::istringstream s(line);
                    int a, b, m, fb = 0;
                    std::string nm;
                    s >> a >> b >> m >> nm >> fb;
                    add_bond(a, b, m, fb);
                }
            } else if (word == "starts") {
                for (int i = 0; i < n; ++i) {
                    next(line);
                    std::istringstream s(line);
                    int v; std::string lab;
                    s >> v >> lab;
                    g.starts.push_back(v); g.start_label.push_back(lab);
                }
            }
        }
    }
    g.adj.assign(2 * g.S, {});
    for (auto& key : order) {
        int b = g.M++;
        g.ba.push_back(key.first); g.bb.push_back(key.second); g.bm.push_back(merged[key]);
        g.bfb.push_back(foldback.count(key) ? foldback[key] : 0);
        g.adj[key.first].push_back({b, key.second});
        if (key.first != key.second) g.adj[key.second].push_back({b, key.first});
    }
    return g;
}

// ------------------------------------------------------------------ memo table (open addressing)
struct Entry {
    u128 count;
    uint64_t key[KW];
    int32_t best;  // -1 = empty slot
    int32_t pad;
};

struct Table {
    std::vector<Entry> slots;
    size_t used = 0, limit_bytes;
    explicit Table(size_t limit) : limit_bytes(limit) { grow(1 << 16); }
    static uint64_t mix(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    static uint64_t hash(const uint64_t* k) {
        uint64_t h = 0;
        for (int i = 0; i < KW; ++i) h = mix(h ^ k[i]);
        return h;
    }
    void grow(size_t n) {
        size_t bytes = n * sizeof(Entry) + slots.size() * sizeof(Entry);
        if (bytes > limit_bytes)
            throw std::runtime_error("memo table would exceed --max-mem-gb (" + std::to_string(used) +
                                     " states stored)");
        std::vector<Entry> old;
        old.swap(slots);
        Entry empty{};
        empty.best = -1;
        slots.assign(n, empty);
        used = 0;
        for (auto& e : old)
            if (e.best >= 0) insert(e.key, e.best, e.count);
    }
    Entry* find(const uint64_t* k) {
        size_t mask = slots.size() - 1, i = hash(k) & mask;
        while (true) {
            Entry& e = slots[i];
            if (e.best < 0) return nullptr;
            if (std::memcmp(e.key, k, sizeof(e.key)) == 0) return &e;
            i = (i + 1) & mask;
        }
    }
    void insert(const uint64_t* k, int best, u128 count) {
        if ((used + 1) * 10 > slots.size() * 7) grow(slots.size() * 2);
        size_t mask = slots.size() - 1, i = hash(k) & mask;
        while (slots[i].best >= 0) i = (i + 1) & mask;
        std::memcpy(slots[i].key, k, sizeof(slots[i].key));
        slots[i].best = best;
        slots[i].count = count;
        ++used;
    }
};

// The value maximized at each state: with --max-foldbacks, the number of fold-back-shaped bonds
// traversed comes first and the number of traversals breaks ties; otherwise only the traversals
// count.  Both fit in one int32 as (fold-backs << 12) | traversals.
static bool MAX_FB = false;
static inline int value(int fb, int len) { return MAX_FB ? (fb << 12) | len : len; }
static inline int val_fb(int v) { return MAX_FB ? (v >> 12) : 0; }
static inline int val_len(int v) { return MAX_FB ? (v & 0xfff) : v; }

// ------------------------------------------------------------------ solver
struct Solver {
    const Graph& g;
    std::vector<int> brem, crem;
    std::vector<int> off, width;
    uint64_t key[KW];
    Table table;
    Solver(const Graph& graph, size_t limit) : g(graph), table(limit) {
        brem = g.bm;
        crem = g.cap;
        // A bond needs no field in the key when its residual follows from the others.  Every
        // traversal of segment s touches both of its nodes once, so the bonds used at one node of s,
        // less those at the other, are fixed by the start and current nodes.  Each segment thus
        // settles one bond it touches unevenly once all its other such bonds are settled
        // (reference adjacencies, chain by chain from the window ends).
        std::vector<std::vector<std::pair<int, int>>> touch(g.S);  // segment -> (bond, coefficient)
        for (int b = 0; b < g.M; ++b) {
            std::map<int, int> c;
            for (int x : {g.ba[b], g.bb[b]}) c[x >> 1] += (x & 1) ? -1 : 1;
            for (auto& sc : c) if (sc.second) touch[sc.first].push_back({b, sc.second});
        }
        std::vector<char> settled(g.M, 0), used(g.S, 0);
        for (bool changed = true; changed;) {
            changed = false;
            for (int s = 0; s < g.S; ++s) {
                if (used[s]) continue;
                int open = -1, n = 0;
                for (auto& bc : touch[s]) if (!settled[bc.first]) { open = bc.first; ++n; }
                if (n == 1) { settled[open] = 1; used[s] = 1; changed = true; }
            }
        }
        int bit = 16;  // bits 0..7 start node, 8..15 current node
        for (int b = 0; b < g.M; ++b) {
            if (settled[b]) { off.push_back(-1); width.push_back(0); continue; }
            int w = 1;
            while ((1 << w) <= g.bm[b]) ++w;
            if (bit % 64 + w > 64) bit += 64 - bit % 64;  // a field never straddles two words
            off.push_back(bit); width.push_back(w);
            bit += w;
        }
        if (bit > 64 * KW) throw std::runtime_error("state does not fit in 256 bits");
        if (2 * g.S > 255) throw std::runtime_error("too many segments (node ids must fit in 8 bits)");
        std::memset(key, 0, sizeof(key));
        for (int b = 0; b < g.M; ++b) add_bond(b, g.bm[b]);
    }
    void add_bond(int b, int delta) {  // signed change of the packed residual of bond b
        int o = off[b];
        if (o < 0) return;
        uint64_t v = uint64_t(delta < 0 ? -delta : delta) << (o % 64);
        if (delta < 0) key[o / 64] -= v; else key[o / 64] += v;
    }
    void set_nodes(int start, int cur) {
        key[0] = (key[0] & ~uint64_t(0xffff)) | uint64_t(start) | (uint64_t(cur) << 8);
    }
    int start = 0;

    // maximum further traversals and number of distinct continuations from exit node w
    std::pair<int, u128> solve(int w) {
        set_nodes(start, w);
        if (Entry* e = table.find(key)) return {e->best, e->count};
        uint64_t here[KW];
        std::memcpy(here, key, sizeof(key));
        int best = 0;
        u128 count = 1;
        for (auto& bx : g.adj[w]) {
            int b = bx.first, x = bx.second, s = x >> 1;
            if (!brem[b] || !crem[s]) continue;
            --brem[b]; --crem[s]; add_bond(b, -1);
            auto r = solve(x ^ 1);
            ++brem[b]; ++crem[s]; add_bond(b, +1);
            int v = r.first + value(g.bfb[b], 1);
            if (v > best) { best = v; count = r.second; }
            else if (v == best) count += r.second;
        }
        std::memcpy(key, here, sizeof(key));
        table.insert(key, best, count);
        return {best, count};
    }
    Entry* state_after(int b, int x) {  // lookup the state reached by taking bond b into node x
        --brem[b]; --crem[x >> 1]; add_bond(b, -1);
        set_nodes(start, x ^ 1);
        Entry* e = table.find(key);
        ++brem[b]; ++crem[x >> 1]; add_bond(b, +1);
        return e;
    }
};

static std::string letter(const Graph& g, int x) {
    const std::string& n = g.name[x >> 1];
    if ((x & 1) == 0) return n;
    return g.signed_ids ? "-" + n : rc_name(n);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: lwalk GRAPH.txt [--starts L1,L2] [--max-mem-gb G] "
                     "[--enumerate | --sample N [--seed S]] [-o PREFIX] [--dedup-rc] [--max-walks N]\n";
        return 2;
    }
    std::string path = argv[1], prefix, starts_arg;
    double mem_gb = 4;
    bool enumerate = false, dedup = false;
    long long sample = 0, max_walks = 5000000, chunk = 0;
    unsigned long long seed = 1;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() { if (i + 1 >= argc) { std::cerr << a << " needs a value\n"; std::exit(2); } return std::string(argv[++i]); };
        if (a == "--starts") starts_arg = val();
        else if (a == "--max-mem-gb") mem_gb = std::atof(val().c_str());
        else if (a == "--enumerate") enumerate = true;
        else if (a == "--sample") sample = std::atoll(val().c_str());
        else if (a == "--seed") seed = std::strtoull(val().c_str(), nullptr, 10);
        else if (a == "-o") prefix = val();
        else if (a == "--dedup-rc") dedup = true;
        else if (a == "--max-walks") max_walks = std::atoll(val().c_str());
        else if (a == "--chunk") chunk = std::atoll(val().c_str());
        else if (a == "--max-foldbacks") MAX_FB = true;
        else { std::cerr << "unknown argument " << a << "\n"; return 2; }
    }
    try {
        Graph g = read_graph(path);
        std::vector<int> starts;
        std::vector<std::string> labels;
        for (size_t i = 0; i < g.starts.size(); ++i) {
            bool ok = starts_arg.empty();
            std::stringstream ss(starts_arg);
            std::string item;
            while (std::getline(ss, item, ',')) ok = ok || item == g.start_label[i];
            if (ok) { starts.push_back(g.starts[i]); labels.push_back(g.start_label[i]); }
        }
        if (starts.empty())
            throw std::runtime_error(starts_arg.empty() ? "the graph has no free end"
                                                        : "no free end matches --starts " + starts_arg);
        std::set<int> start_set(starts.begin(), starts.end());
        Solver sv(g, size_t(mem_gb * 1073741824.0));
        auto t0 = std::chrono::steady_clock::now();
        int L = 0;
        u128 total = 0;
        std::vector<std::pair<int, u128>> per(starts.size(), {0, 0});
        long total_trav = 0;
        for (int c : g.cap) total_trav += c;
        if (MAX_FB && total_trav > 0xfff)
            throw std::runtime_error("--max-foldbacks supports at most 4095 traversals, the graph has " +
                                     std::to_string(total_trav));
        for (size_t i = 0; i < starts.size(); ++i) {
            int s0 = starts[i];
            if (!sv.crem[s0 >> 1]) continue;
            sv.start = s0;
            --sv.crem[s0 >> 1];
            auto r = sv.solve(s0 ^ 1);
            ++sv.crem[s0 >> 1];
            const int v0 = value(0, 1) + r.first;   // the start segment: one traversal, no bond
            per[i] = {v0, r.second};
            if (v0 > L) { L = v0; total = 0; }
            if (v0 == L) total += r.second;
        }
        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << g.header << "\n";
        std::cout << (MAX_FB ? "best walk: " : "longest walk: ") << val_len(L) << " of " << total_trav
                  << " traversals" << (MAX_FB ? ", " + std::to_string(val_fb(L)) + " fold-back traversals" : "")
                  << "; distinct best walks (each reading from an allowed start): " << u128str(total) << "\n";
        for (size_t i = 0; i < starts.size(); ++i)
            std::cout << "  start " << labels[i] << ": " << val_len(per[i].first) << " traversals"
                      << (MAX_FB ? ", " + std::to_string(val_fb(per[i].first)) + " fold-backs" : "")
                      << ", walks of that value " << u128str(per[i].second)
                      << (per[i].first == L ? "" : "  (worse than the best)") << "\n";
        std::cout << "states " << sv.table.used << ", table " << (sv.table.slots.size() * sizeof(Entry) >> 20)
                  << " MiB, " << secs << " s\n" << std::flush;
        if (!enumerate && sample == 0) return 0;
        if (prefix.empty()) { std::cerr << "-o PREFIX is required to write walks\n"; return 2; }
        if (total == 0) { std::cerr << "no walk from the selected starts\n"; return 3; }
        if (enumerate && total > u128(max_walks)) {
            std::cerr << "refusing to enumerate " << u128str(total) << " walks (> --max-walks " << max_walks << ")\n";
            return 4;
        }
        // with --chunk N the targets go to PREFIX.targets.0000.txt, .0001, ... (N walks each)
        std::ofstream fw(prefix + ".walks.txt"), ft;
        long long written = 0, skipped = 0;
        auto open_targets = [&]() {
            if (ft.is_open()) ft.close();
            char buf[32];
            std::snprintf(buf, sizeof(buf), ".%04lld", chunk ? written / chunk : 0LL);
            ft.open(prefix + ".targets" + (chunk ? std::string(buf) : std::string()) + ".txt");
        };
        open_targets();
        std::set<std::string> seen;
        auto emit = [&](const std::vector<int>& path) -> bool {
            std::string w, t, r;
            for (size_t k = 0; k < path.size(); ++k) {
                if (k) { w += ' '; t += ';'; }
                w += letter(g, path[k]);
                t += (path[k] & 1) ? g.rev[path[k] >> 1] : g.fwd[path[k] >> 1];
            }
            if (dedup && start_set.count(path.back() ^ 1)) {
                for (size_t k = path.size(); k-- > 0;) {
                    if (k + 1 < path.size()) r += ' ';
                    r += letter(g, path[k] ^ 1);
                }
                if (r < w) { ++skipped; return false; }
            }
            if (sample) {
                if (seen.count(w)) return false;
                seen.insert(w);
            }
            if (chunk && written > 0 && written % chunk == 0) open_targets();
            fw << w << "\n";
            ft << t << "\n";
            ++written;
            return true;
        };
        std::vector<int> pathv;
        if (enumerate) {
            // depth-first over optimal moves only
            std::function<void(int, int)> rec = [&](int w, int left) {
                if (left == 0) { emit(pathv); return; }
                for (auto& bx : g.adj[w]) {
                    int b = bx.first, x = bx.second, s = x >> 1;
                    if (!sv.brem[b] || !sv.crem[s]) continue;
                    Entry* e = sv.state_after(b, x);
                    const int step = value(g.bfb[b], 1);
                    if (!e || step + e->best != left) continue;
                    --sv.brem[b]; --sv.crem[s]; sv.add_bond(b, -1);
                    pathv.push_back(x);
                    rec(x ^ 1, left - step);
                    pathv.pop_back();
                    ++sv.brem[b]; ++sv.crem[s]; sv.add_bond(b, +1);
                }
            };
            for (size_t i = 0; i < starts.size(); ++i) {
                if (per[i].first != L) continue;
                int s0 = starts[i];
                sv.start = s0;
                --sv.crem[s0 >> 1];
                pathv.assign(1, s0);
                rec(s0 ^ 1, L - value(0, 1));
                ++sv.crem[s0 >> 1];
            }
        } else {
            std::mt19937_64 rng(seed);
            auto rnd = [&](u128 n) { u128 v = (u128(rng()) << 64) | rng(); return v % n; };
            long long attempts = 0;
            while (written < sample && attempts < 50 * sample) {
                ++attempts;
                u128 r = rnd(total);
                size_t i = 0;
                for (; i < starts.size(); ++i) {
                    if (per[i].first != L) continue;
                    if (r < per[i].second) break;
                    r -= per[i].second;
                }
                int s0 = starts[i];
                sv.start = s0;
                --sv.crem[s0 >> 1];
                pathv.assign(1, s0);
                int w = s0 ^ 1, left = L - value(0, 1);
                std::vector<std::pair<int, int>> taken;
                while (left > 0) {
                    std::vector<std::pair<std::pair<int, int>, u128>> opts;
                    u128 sum = 0;
                    for (auto& bx : g.adj[w]) {
                        int b = bx.first, x = bx.second;
                        if (!sv.brem[b] || !sv.crem[x >> 1]) continue;
                        Entry* e = sv.state_after(b, x);
                        if (!e || value(g.bfb[b], 1) + e->best != left) continue;
                        opts.push_back({bx, e->count});
                        sum += e->count;
                    }
                    u128 pick = rnd(sum);
                    size_t k = 0;
                    while (pick >= opts[k].second) { pick -= opts[k].second; ++k; }
                    int b = opts[k].first.first, x = opts[k].first.second;
                    --sv.brem[b]; --sv.crem[x >> 1]; sv.add_bond(b, -1);
                    taken.push_back({b, x});
                    pathv.push_back(x);
                    w = x ^ 1;
                    left -= value(g.bfb[b], 1);
                }
                emit(pathv);
                for (auto& bx : taken) { ++sv.brem[bx.first]; ++sv.crem[bx.second >> 1]; sv.add_bond(bx.first, +1); }
                ++sv.crem[s0 >> 1];
            }
        }
        std::cout << "wrote " << written << " walks to " << prefix << ".{walks,targets}.txt"
                  << (dedup ? " (" + std::to_string(skipped) + " reverse-complement readings skipped)" : "") << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }
    return 0;
}
