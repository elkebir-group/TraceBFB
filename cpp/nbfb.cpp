#include "nbfb.hpp"

#include <string>

namespace nbfb {

// ---------------------------------------------------------------- budget

MemoryBudget::MemoryBudget(std::size_t limit_bytes)
    : limit_(limit_bytes), used_(0) {}

void MemoryBudget::acquire(std::size_t bytes) {
    std::size_t now = used_.fetch_add(bytes, std::memory_order_relaxed) + bytes;
    if (now > limit_) {
        used_.fetch_sub(bytes, std::memory_order_relaxed);
        throw MemoryLimitExceeded(
            "memory limit exceeded: the search needs more than " +
            std::to_string(limit_) + " bytes");
    }
}

void MemoryBudget::release(std::size_t bytes) {
    used_.fetch_sub(bytes, std::memory_order_relaxed);
}

}  // namespace nbfb
