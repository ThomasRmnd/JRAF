#ifndef TIMESTAMP_HPP_
#define TIMESTAMP_HPP_

#include <ctime>
#include <iomanip>
#include <iostream>

class timestamp {

public:

    timestamp() : sec{0l}, nsec{0} {}
    timestamp(time_t sec_, int nsec_) : sec{sec_}, nsec{nsec_} {
        normalize();
    }

    ~timestamp() = default;

    time_t sec;
    int nsec;

private:

    void normalize() {
        while (nsec >= 1000000000) {
            nsec -= 1000000000;
            ++sec;
        }
        while (nsec < 0) {
            nsec += 1000000000;
            --sec;
        }
    }

};

inline timestamp operator+(const timestamp& lhs, const timestamp& rhs) {
    return timestamp{lhs.sec + rhs.sec, lhs.nsec + rhs.nsec};
}

inline timestamp operator-(const timestamp& lhs, const timestamp& rhs) {
    return timestamp{lhs.sec - rhs.sec, lhs.nsec - rhs.nsec};
}

inline bool operator==(const timestamp& lhs, const timestamp& rhs) {
    return lhs.sec == rhs.sec && lhs.nsec == rhs.nsec;
}

inline bool operator!=(const timestamp& lhs, const timestamp& rhs) {
    return lhs.sec != rhs.sec || lhs.nsec != rhs.nsec;
}

inline bool operator<(const timestamp& lhs, const timestamp& rhs) {
    return (lhs.sec < rhs.sec) || (lhs.sec == rhs.sec && lhs.nsec < rhs.nsec);
}

inline bool operator<=(const timestamp& lhs, const timestamp& rhs) {
    return (lhs.sec < rhs.sec) || (lhs.sec == rhs.sec && lhs.nsec <= rhs.nsec);
}

inline bool operator>(const timestamp& lhs, const timestamp& rhs) {
    return (lhs.sec > rhs.sec) || (lhs.sec == rhs.sec && lhs.nsec > rhs.nsec);
}

inline bool operator>=(const timestamp& lhs, const timestamp& rhs) {
    return (lhs.sec > rhs.sec) || (lhs.sec == rhs.sec && lhs.nsec >= rhs.nsec);
}

inline std::ostream& operator<<(std::ostream& os, const timestamp& ts) {
    std::tm* tm_ptr = std::gmtime(&ts.sec);
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", tm_ptr);
    os << buffer << '.' << std::setw(9) << std::setfill('0') << ts.nsec;
    return os;
}

inline double timestamp_to_double(const timestamp& ts) {
    return static_cast<double>(ts.sec) + static_cast<double>(ts.nsec) * 1.0e-9;
}

#endif // TIMESTAMP_HPP_