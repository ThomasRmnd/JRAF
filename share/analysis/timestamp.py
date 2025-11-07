from datetime import datetime, timezone 

class TimeStamp:
    __slots__ = ('sec', 'nsec')

    def __init__(self, sec: int, nsec: int):
        self.sec = int(sec)
        self.nsec = int(nsec)
        self.normalize()

    def normalize(self):
        while (self.nsec < 0):
            self.nsec += 1_000_000_000
            self.sec -= 1
        while (self.nsec >= 1_000_000_000):
            self.nsec -= 1_000_000_000
            self.sec += 1

    def to_sec(self):
        return self.sec + self.nsec * 1e-9

    def to_nsec(self):
        return self.sec * 1_000_000_000 + self.nsec

    def __add__(self, other):
        return TimeStamp(self.sec + other.sec, self.nsec + other.nsec)

    def __sub__(self, other):
        return TimeStamp(self.sec - other.sec, self.nsec - other.nsec)

    def __lt__(self, other):
        return (self.sec, self.nsec) < (other.sec, other.nsec)

    def __le__(self, other):
        return (self.sec, self.nsec) <= (other.sec, other.nsec)

    def __gt__(self, other):
        return (self.sec, self.nsec) > (other.sec, other.nsec)

    def __ge__(self, other):
        return (self.sec, self.nsec) >= (other.sec, other.nsec)

    def __eq__(self, other):
        return (self.sec, self.nsec) == (other.sec, other.nsec)

    def __repr__(self):
        dt = datetime.fromtimestamp(self.sec, tz=timezone.utc)
        return dt.strftime(f"%Y-%m-%d %H:%M:%S.{self.nsec:09d} UTC")