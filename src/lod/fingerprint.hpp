#pragma once
#include <bit>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include "alod/lod.hpp"
namespace alod {
class FingerprintBuilder {
public:
    void add_u64(std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash_ ^= static_cast<unsigned char>((value >> (8 * byte)) & 0xffU);
            hash_ *= 1099511628211ULL;
        }
    }

    void add_i64(std::int64_t value) {
        add_u64(static_cast<std::uint64_t>(value));
    }

    void add_double(double value) {
        // Positive and negative zero are numerically indistinguishable in all
        // estimator operations, so give them one canonical representation.
        if (value == 0.0) value = 0.0;
        add_u64(std::bit_cast<std::uint64_t>(value));
    }

    void add_complex(const Complex &value) {
        add_double(value.real());
        add_double(value.imag());
    }

    void add_string(const std::string &value) {
        add_u64(value.size());
        for (unsigned char byte : value) {
            hash_ ^= byte;
            hash_ *= 1099511628211ULL;
        }
    }

    std::string finish() const {
        std::ostringstream encoded;
        encoded << "fnv1a64:" << std::hex << std::setw(16)
                << std::setfill('0') << hash_;
        return encoded.str();
    }

private:
    std::uint64_t hash_ = 14695981039346656037ULL;
};

}
