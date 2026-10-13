#include "nbfb.hpp"

#include <string>
#include <utility>

namespace nbfb {

// ---------------------------------------------------------------- budget

MemoryBudget::MemoryBudget(std::size_t limit_bytes)
    : limit_(limit_bytes), used_(0), peak_(0) {}

void MemoryBudget::acquire(std::size_t bytes) {
    std::size_t now = used_.fetch_add(bytes, std::memory_order_relaxed) + bytes;
    if (now > limit_) {
        used_.fetch_sub(bytes, std::memory_order_relaxed);
        throw MemoryLimitExceeded(
            "memory limit exceeded: the search needs more than " +
            std::to_string(limit_) + " bytes");
    }
    std::size_t seen = peak_.load(std::memory_order_relaxed);
    while (now > seen &&
           !peak_.compare_exchange_weak(seen, now, std::memory_order_relaxed)) {
    }
}

void MemoryBudget::release(std::size_t bytes) {
    used_.fetch_sub(bytes, std::memory_order_relaxed);
}

PreparedSources::PreparedSources(std::vector<std::vector<int>> sources)
    : forward(std::move(sources)) {
    reverse.reserve(forward.size());
    for (const std::vector<int>& source : forward) {
        std::vector<int> other(source.rbegin(), source.rend());
        for (int& value : other) value = -value;
        reverse.push_back(std::move(other));
    }
}

}  // namespace nbfb
