#ifndef RECMUONALG_UTILS_DETECTORTYPE_HPP_
#define RECMUONALG_UTILS_DETECTORTYPE_HPP_

enum class DetectorType {
    UNKNOWN = 0,
    CD = 1 << 0,
    WP = 1 << 1,
    TT = 1 << 2
};

inline constexpr DetectorType operator|(DetectorType lhs, DetectorType rhs) noexcept {
    return static_cast<DetectorType>(static_cast<int>(lhs) | static_cast<int>(rhs));
}

inline constexpr DetectorType operator&(DetectorType lhs, DetectorType rhs) noexcept {
    return static_cast<DetectorType>(static_cast<int>(lhs) & static_cast<int>(rhs));
}

inline constexpr DetectorType operator^(DetectorType lhs, DetectorType rhs) noexcept {
    return static_cast<DetectorType>(static_cast<int>(lhs) ^ static_cast<int>(rhs));
}

inline constexpr DetectorType operator~(DetectorType rhs) noexcept {
    return static_cast<DetectorType>(~static_cast<int>(rhs));
}

inline constexpr DetectorType& operator|=(DetectorType& lhs, DetectorType rhs) noexcept {
    lhs = lhs | rhs;
    return lhs;
}

inline constexpr DetectorType& operator&=(DetectorType& lhs, DetectorType rhs) noexcept {
    lhs = lhs & rhs;
    return lhs;
}

inline constexpr DetectorType& operator^=(DetectorType& lhs, DetectorType rhs) noexcept {
    lhs = lhs ^ rhs;
    return lhs;
}

#endif // RECMUONALG_UTILS_DETECTORTYPE_HPP_