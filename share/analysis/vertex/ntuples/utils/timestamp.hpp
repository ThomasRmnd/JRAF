#ifndef UTILS_TIMESTAMP_HPP_
#define UTILS_TIMESTAMP_HPP_

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

class timestamp {

public:

    timestamp() : sec{0l}, nsec{0} {}
    timestamp(time_t sec_, int nsec_) : sec{sec_}, nsec{nsec_} {
        normalize();
    }

    timestamp(const timestamp& other) = default;
    timestamp(timestamp&& other) = default;

    ~timestamp() = default;

    timestamp& operator=(const timestamp& other) = default;
    timestamp& operator=(timestamp&& other) = default;

    timestamp& operator+=(const timestamp& other) {
        sec += other.sec;
        nsec += other.nsec;
        normalize();
        return *this;
    }

    timestamp& operator-=(const timestamp& other) {
        sec -= other.sec;
        nsec -= other.nsec;
        normalize();
        return *this;
    }

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
    
    std::ios::fmtflags f(os.flags());
    char fill = os.fill();
    
    os << buffer << '.' << std::setw(9) << std::setfill('0') << ts.nsec;
    
    os.flags(f);
    os.fill(fill);
    
    return os;
}

inline double timestamp_to_double(const timestamp& ts) {
    return static_cast<double>(ts.sec) + static_cast<double>(ts.nsec) * 1.0e-9;
}

inline std::string timestamp_to_string(const timestamp& ts) {
    std::ostringstream oss;
    oss << ts;
    return oss.str();
}

#endif // UTILS_TIMESTAMP_HPP_