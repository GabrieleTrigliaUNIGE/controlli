"""Transfer function and frequency sampling."""

import numpy as np


class TransferFunction:
    """Transfer function in time-constant form."""

    def __init__(self, K=1.0, g=0, real_zeros=None, real_poles=None,
                 complex_zeros=None, complex_poles=None):
        self.K = float(K)
        self.g = int(g)
        self.real_zeros = list(real_zeros or [])         # time constants for (1+s*tau)
        self.real_poles = list(real_poles or [])          # time constants
        self.complex_zeros = list(complex_zeros or [])    # natural frequency and damping ratio pairs
        self.complex_poles = list(complex_poles or [])    # natural frequency and damping ratio pairs

    # ---------- exact frequency response for the Nyquist plot -----------------
    def response(self, w):
        w = np.asarray(w, dtype=float)
        s = 1j * w

        num = self.K * np.ones_like(s, dtype=complex)
        if self.g > 0:
            num = num * s ** self.g
        for tau in self.real_zeros:
            num = num * (1 + s * tau)
        for wn, zeta in self.complex_zeros:
            num = num * (1 + 2 * zeta * s / wn + (s / wn) ** 2)

        den = np.ones_like(s, dtype=complex)
        if self.g < 0:
            den = den * s ** (-self.g)
        for tau in self.real_poles:
            den = den * (1 + s * tau)
        for wn, zeta in self.complex_poles:
            den = den * (1 + 2 * zeta * s / wn + (s / wn) ** 2)

        return num / den

    # ---------- break frequencies for automatic frequency range selection ---
    def breakpoints(self):
        bp = []
        for tau in self.real_zeros + self.real_poles:
            if tau != 0:
                bp.append(1.0 / abs(tau))
        for wn, _ in self.complex_zeros + self.complex_poles:
            if wn != 0:
                bp.append(wn)
        return bp

    def default_freq_range(self):
        bp = self.breakpoints()
        if bp:
            wmin = min(bp) / 100.0
            wmax = max(bp) * 100.0
        else:
            wmin, wmax = 1e-2, 1e2
        if np.log10(wmax / wmin) < 3:
            center = np.sqrt(wmin * wmax)
            span = 10 ** 1.5
            wmin, wmax = center / span, center * span
        return wmin, wmax

    # ---------- first-order factor contribution ----------
    @staticmethod
    def _real_factor(w, tau, sign, asymptotic):
        """Use sign +1 for a zero and -1 for a pole."""
        wc = 1.0 / abs(tau)
        if not asymptotic:
            mag = sign * 20 * np.log10(np.sqrt(1 + (w * tau) ** 2))
            phase = sign * np.degrees(np.arctan(w * tau))
        else:
            mag = sign * np.where(w < wc, 0.0, 20 * np.log10(w / wc))
            w1, w2 = wc / 10.0, wc * 10.0
            ramp_end = sign * np.sign(tau) * 90.0
            phase = np.where(
                w <= w1, 0.0,
                np.where(w >= w2, ramp_end,
                         sign * np.sign(tau) * 45.0 * np.log10(w / w1))
            )
        return mag, phase

    # ---------- second-order factor contribution ------
    @staticmethod
    def _complex_factor(w, wn, zeta, sign, asymptotic):
        """Use sign +1 for a zero and -1 for a pole. The asymptotic phase assumes positive damping."""
        if not asymptotic:
            re = 1 - (w / wn) ** 2
            im = 2 * zeta * w / wn
            mag = sign * 20 * np.log10(np.sqrt(re ** 2 + im ** 2))
            phase = sign * np.degrees(np.arctan2(im, re))
        else:
            mag = sign * np.where(w < wn, 0.0, 40 * np.log10(w / wn))
            w1, w2 = wn / 10.0, wn * 10.0
            ramp_end = sign * 180.0
            phase = np.where(
                w <= w1, 0.0,
                np.where(w >= w2, ramp_end,
                         sign * 90.0 * np.log10(w / w1))
            )
        return mag, phase

    # ---------- total magnitude (dB) and phase (degrees) --------------------------
    def bode(self, w, asymptotic=False):
        w = np.asarray(w, dtype=float)
        mag = np.full_like(w, 20 * np.log10(abs(self.K))) + 20 * self.g * np.log10(w)
        phase = np.full_like(w, 180.0 if self.K < 0 else 0.0) + self.g * 90.0

        for tau in self.real_zeros:
            m, p = self._real_factor(w, tau, +1, asymptotic)
            mag, phase = mag + m, phase + p
        for tau in self.real_poles:
            m, p = self._real_factor(w, tau, -1, asymptotic)
            mag, phase = mag + m, phase + p
        for wn, zeta in self.complex_zeros:
            m, p = self._complex_factor(w, wn, zeta, +1, asymptotic)
            mag, phase = mag + m, phase + p
        for wn, zeta in self.complex_poles:
            m, p = self._complex_factor(w, wn, zeta, -1, asymptotic)
            mag, phase = mag + m, phase + p

        return mag, phase

    # ---------- text representation ----------------------------------
    def to_string(self):
        def fmt(x):
            return f"{x:g}"

        num_terms, den_terms = [], []
        if self.g > 0:
            num_terms.append(f"s^{self.g}" if self.g != 1 else "s")
        elif self.g < 0:
            den_terms.append(f"s^{-self.g}" if self.g != -1 else "s")
        for tau in self.real_zeros:
            num_terms.append(f"(1+{fmt(tau)}s)")
        for wn, zeta in self.complex_zeros:
            num_terms.append(f"(1+{fmt(2*zeta/wn)}s+{fmt(1/wn**2)}s^2)")
        for tau in self.real_poles:
            den_terms.append(f"(1+{fmt(tau)}s)")
        for wn, zeta in self.complex_poles:
            den_terms.append(f"(1+{fmt(2*zeta/wn)}s+{fmt(1/wn**2)}s^2)")

        num = " * ".join(num_terms) if num_terms else "1"
        den = " * ".join(den_terms) if den_terms else "1"
        return f"G(s) = {fmt(self.K)} * {num} / [{den}]"


def make_freq_vector(tf, wmin=None, wmax=None, n=1000):
    if wmin is None or wmax is None:
        wmin_d, wmax_d = tf.default_freq_range()
        wmin = wmin if wmin is not None else wmin_d
        wmax = wmax if wmax is not None else wmax_d
    return np.logspace(np.log10(wmin), np.log10(wmax), n)
