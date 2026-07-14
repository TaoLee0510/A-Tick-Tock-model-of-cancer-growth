#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>

namespace atcg3d {

inline std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

inline std::uint64_t rng_word(std::uint64_t seed,
                              std::uint64_t uid,
                              std::uint64_t event_type,
                              std::uint64_t event_sequence,
                              std::uint64_t draw = 0) {
    std::uint64_t value = splitmix64(seed);
    value ^= splitmix64(uid + 0x632be59bd9b4e019ULL);
    value ^= splitmix64(event_type + 0x8cb92baa3f3d8dd7ULL);
    value ^= splitmix64(event_sequence + 0x58f38dedf4c2f89bULL);
    value ^= splitmix64(draw + 0xa0761d6478bd642fULL);
    return splitmix64(value);
}

inline double rng_unit(std::uint64_t seed,
                       std::uint64_t uid,
                       std::uint64_t event_type,
                       std::uint64_t event_sequence,
                       std::uint64_t draw = 0) {
    constexpr double scale = 1.0 / static_cast<double>(1ULL << 53U);
    return static_cast<double>(rng_word(seed, uid, event_type, event_sequence, draw) >> 11U) * scale;
}

inline std::uint64_t rng_bounded(std::uint64_t seed,
                                 std::uint64_t uid,
                                 std::uint64_t event_type,
                                 std::uint64_t event_sequence,
                                 std::uint64_t upper_exclusive,
                                 std::uint64_t draw = 0) {
    if (upper_exclusive == 0) {
        return 0;
    }
    return rng_word(seed, uid, event_type, event_sequence, draw) % upper_exclusive;
}

template <class RandomIt>
void deterministic_shuffle(RandomIt first,
                           RandomIt last,
                           std::uint64_t seed,
                           std::uint64_t uid,
                           std::uint64_t event_type,
                           std::uint64_t event_sequence) {
    using difference_type = typename std::iterator_traits<RandomIt>::difference_type;
    const difference_type count = last - first;
    for (difference_type i = count - 1; i > 0; --i) {
        const auto j = static_cast<difference_type>(
            rng_bounded(seed, uid, event_type, event_sequence, static_cast<std::uint64_t>(i + 1),
                        static_cast<std::uint64_t>(count - i)));
        std::iter_swap(first + i, first + j);
    }
}

}  // namespace atcg3d
