"""Helpers for the language revision: exact replacement and span replacement on the manuscript."""
import io

BASE = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability"
PATH = BASE + r"\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)


class Doc:
    def __init__(self):
        raw = io.open(PATH, encoding="utf-8", newline="").read()
        self.crlf = (CR + LF) in raw
        self.s = raw.replace(CR + LF, LF)

    def rep(self, old, new, count=1):
        n = self.s.count(old)
        assert n == count, (n, old[:90])
        self.s = self.s.replace(old, new)

    def span(self, start, end, new, include_end=True):
        """Replace from the unique ``start`` marker to the first ``end`` marker after it (inclusive)."""
        assert self.s.count(start) == 1, (self.s.count(start), start[:90])
        i = self.s.index(start)
        j = self.s.index(end, i)
        j += len(end) if include_end else 0
        self.s = self.s[:i] + new + self.s[j:]

    def save(self):
        io.open(PATH, "w", encoding="utf-8", newline="").write(self.s.replace(LF, CR + LF) if self.crlf else self.s)
