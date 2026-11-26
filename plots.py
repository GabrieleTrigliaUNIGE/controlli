"""Draw Bode and Nyquist plots on Matplotlib figures."""

import numpy as np


def plot_bode(tf, w, fig):
    ax_m, ax_p = fig.axes

    mag_e, ph_e = tf.bode(w, asymptotic=False)
    mag_a, ph_a = tf.bode(w, asymptotic=True)

    ax_m.semilogx(w, mag_a, '--', color='tab:orange', lw=2, label='Asymptotic')
    ax_m.semilogx(w, mag_e, '-', color='tab:blue', lw=1.3, label='Exact')
    ax_m.axhline(0, color='k', lw=0.6)
    ax_m.set_ylabel('Magnitude |G(jw)| [dB]')
    ax_m.grid(True, which='both', ls=':', alpha=0.6)
    ax_m.legend(loc='best', fontsize=9)
    ax_m.set_title(tf.to_string(), fontsize=10)

    ax_p.semilogx(w, ph_a, '--', color='tab:orange', lw=2, label='Asymptotic')
    ax_p.semilogx(w, ph_e, '-', color='tab:blue', lw=1.3, label='Exact')
    ax_p.set_ylabel('Phase of G(jw) [deg]')
    ax_p.set_xlabel('w [rad/s]')
    ax_p.grid(True, which='both', ls=':', alpha=0.6)
    ymin, ymax = ax_p.get_ylim()
    ticks = np.arange(np.floor(ymin / 45) * 45, np.ceil(ymax / 45) * 45 + 1, 45)
    ax_p.set_yticks(ticks)
    ax_p.legend(loc='best', fontsize=9)

    fig.tight_layout()


def plot_nyquist(tf, w, fig):
    ax = fig.axes[0]
    G = tf.response(w)

    ax.plot(G.real, G.imag, '-', color='tab:blue', lw=1.5, label='w: 0 -> +inf')
    ax.plot(G.real, -G.imag, '--', color='tab:blue', lw=1.0, alpha=0.6,
            label='w: -inf -> 0 (mirrored)')
    ax.plot(-1, 0, 'rx', ms=9, mew=2, label='Critical point (-1, 0)')
    ax.axhline(0, color='k', lw=0.6)
    ax.axvline(0, color='k', lw=0.6)
    ax.set_xlabel('Re[G(jw)]')
    ax.set_ylabel('Im[G(jw)]')
    ax.set_title(tf.to_string(), fontsize=10)
    ax.grid(True, ls=':', alpha=0.6)

    finite = np.isfinite(G.real) & np.isfinite(G.imag)
    if finite.any():
        re, im = G.real[finite], G.imag[finite]
        radius = np.hypot(re, im)
        # Exclude distant tails so the curve near the critical point stays visible.
        typical = np.median(radius)
        cap = max(typical * 15, 3.0)
        core = radius <= cap
        if core.sum() < 10:
            core = np.ones_like(radius, dtype=bool)
        re_c, im_c = re[core], im[core]

        xlo = min(re_c.min(), -1.3)
        xhi = max(re_c.max(), 0.3)
        center_x = (xlo + xhi) / 2
        # Include both mirrored branches and leave room around the critical point.
        half_span = 1.3 * max((xhi - xlo) / 2, np.max(np.abs(im_c)), 0.3)
        ax.set_xlim(center_x - half_span, center_x + half_span)
        ax.set_ylim(-half_span, half_span)
        ax.set_aspect('equal', adjustable='datalim')

    ax.legend(loc='best', fontsize=9)

    mid = len(w) // 2
    ax.annotate('', xy=(G.real[mid + 1], G.imag[mid + 1]),
                 xytext=(G.real[mid], G.imag[mid]),
                 arrowprops=dict(arrowstyle='->', color='tab:blue', lw=1.5))

    fig.tight_layout()
