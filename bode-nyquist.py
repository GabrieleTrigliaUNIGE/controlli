#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diagrammi di Bode (asintotici/ideali) e di Nyquist
per funzioni di trasferimento espresse in forma a costanti di tempo:

    G(s) = K * s^g * Prod_i(1 + s*tau_zi) * Prod_k(1 + 2*zeta_zk*s/wn_zk + (s/wn_zk)^2)
           -------------------------------------------------------------------------
                 Prod_j(1 + s*tau_pj) * Prod_l(1 + 2*zeta_pl*s/wn_pl + (s/wn_pl)^2)

dove g e' il "tipo" del sistema:
    g < 0  ->  |g| poli nell'origine
    g > 0  ->   g  zeri nell'origine

Per il diagramma di Bode viene mostrato sia il tracciato ASINTOTICO (le
classiche rette a tratti, l'approssimazione "ideale" usata nello studio
dei sistemi) sia quello ESATTO, per potersi confrontare.
Per il diagramma di Nyquist viene mostrato il tracciato esatto di G(jw).

In ogni finestra dei grafici e' possibile zoomare con la rotellina del
mouse (centrato sul punto sotto il cursore), con i bottoni "Zoom +/-",
oppure con lo strumento a riquadro della toolbar di matplotlib; il
bottone "Reset vista" e l'icona "Home" della toolbar ripristinano la
visualizzazione iniziale.

Esecuzione:
    python3 bode_nyquist_gui.py

Dipendenze:
    pip install numpy matplotlib
(tkinter e' incluso in Python; su alcune distribuzioni Linux va
installato separatamente con: sudo apt install python3-tk)
"""

try:
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except ImportError as e:
    raise SystemExit(f"Dipendenza mancante: {e}\nInstallare con: pip install numpy matplotlib")

try:
    import tkinter as tk
    from tkinter import ttk, messagebox
except ImportError:
    raise SystemExit("Questo script richiede tkinter (incluso in Python; "
                      "su Linux potrebbe servire: sudo apt install python3-tk)")


# ============================================================================
# 1. MODELLO DELLA FUNZIONE DI TRASFERIMENTO
# ============================================================================

class TransferFunction:
    """Funzione di trasferimento in forma a costanti di tempo."""

    def __init__(self, K=1.0, g=0, real_zeros=None, real_poles=None,
                 complex_zeros=None, complex_poles=None):
        self.K = float(K)
        self.g = int(g)
        self.real_zeros = list(real_zeros or [])         # lista di tau -> (1+s*tau)
        self.real_poles = list(real_poles or [])          # lista di tau
        self.complex_zeros = list(complex_zeros or [])    # lista di (wn, zeta)
        self.complex_poles = list(complex_poles or [])    # lista di (wn, zeta)

    # ---------- valore esatto di G(jw), usato per Nyquist -----------------
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

    # ---------- punti notevoli, per stimare un range di frequenze utile ---
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

    # ---------- contributo di un singolo fattore reale (1+s*tau) ----------
    @staticmethod
    def _real_factor(w, tau, sign, asymptotic):
        """sign = +1 per zero, -1 per polo."""
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

    # ---------- contributo di una coppia complessa di secondo ordine ------
    @staticmethod
    def _complex_factor(w, wn, zeta, sign, asymptotic):
        """sign = +1 per zero, -1 per polo. Si assume zeta > 0 (stabile)."""
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

    # ---------- modulo (dB) e fase (gradi) totali --------------------------
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

    # ---------- rappresentazione testuale ----------------------------------
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


# ============================================================================
# 2. FUNZIONI DI PLOT
# ============================================================================

def make_freq_vector(tf, wmin=None, wmax=None, n=1000):
    if wmin is None or wmax is None:
        wmin_d, wmax_d = tf.default_freq_range()
        wmin = wmin if wmin is not None else wmin_d
        wmax = wmax if wmax is not None else wmax_d
    return np.logspace(np.log10(wmin), np.log10(wmax), n)


def plot_bode(tf, w, fig):
    ax_m, ax_p = fig.axes

    mag_e, ph_e = tf.bode(w, asymptotic=False)
    mag_a, ph_a = tf.bode(w, asymptotic=True)

    ax_m.semilogx(w, mag_a, '--', color='tab:orange', lw=2, label='Asintotico (ideale)')
    ax_m.semilogx(w, mag_e, '-', color='tab:blue', lw=1.3, label='Esatto')
    ax_m.axhline(0, color='k', lw=0.6)
    ax_m.set_ylabel('Modulo |G(jw)| [dB]')
    ax_m.grid(True, which='both', ls=':', alpha=0.6)
    ax_m.legend(loc='best', fontsize=9)
    ax_m.set_title(tf.to_string(), fontsize=10)

    ax_p.semilogx(w, ph_a, '--', color='tab:orange', lw=2, label='Asintotico (ideale)')
    ax_p.semilogx(w, ph_e, '-', color='tab:blue', lw=1.3, label='Esatto')
    ax_p.set_ylabel('Fase G(jw) [grad]')
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
            label='w: -inf -> 0 (speculare)')
    ax.plot(-1, 0, 'rx', ms=9, mew=2, label='punto (-1, 0)')
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
        # Per i sistemi con poli/zeri nell'origine il modulo tende a
        # infinito per w->0 (o w->inf): quella coda andrebbe esclusa dal
        # calcolo dei limiti degli assi, altrimenti schiaccia la parte
        # interessante della curva (zona vicino a -1).
        typical = np.median(radius)
        cap = max(typical * 15, 3.0)
        core = radius <= cap
        if core.sum() < 10:
            core = np.ones_like(radius, dtype=bool)
        re_c, im_c = re[core], im[core]

        xlo, xhi = re_c.min(), re_c.max()
        ylo, yhi = im_c.min(), im_c.max()
        xlo, xhi = min(xlo, -1.3), max(xhi, 0.3)
        ylo, yhi = min(ylo, -0.3), max(yhi, 0.3)
        padx = 0.15 * (xhi - xlo + 1e-9)
        pady = 0.15 * (yhi - ylo + 1e-9)
        ax.set_xlim(xlo - padx, xhi + padx)
        ax.set_ylim(ylo - pady, yhi + pady)

        xspan = (xhi + padx) - (xlo - padx)
        yspan = (yhi + pady) - (ylo - pady)
        if xspan > 0 and yspan > 0 and max(xspan, yspan) / min(xspan, yspan) <= 8:
            ax.set_aspect('equal', adjustable='box')

    ax.legend(loc='best', fontsize=9)

    mid = len(w) // 2
    ax.annotate('', xy=(G.real[mid + 1], G.imag[mid + 1]),
                 xytext=(G.real[mid], G.imag[mid]),
                 arrowprops=dict(arrowstyle='->', color='tab:blue', lw=1.5))

    fig.tight_layout()


# ============================================================================
# 2bis. ZOOM SUI GRAFICI (rotellina del mouse + pulsanti +/-)
# ============================================================================

def _zoom_range(lo, hi, center, factor, scale):
    """Restringe/allarga l'intervallo [lo, hi] di un fattore 'factor'
    attorno a 'center', rispettando la scala (lineare o logaritmica)."""
    if scale == 'log':
        return center * (lo / center) ** factor, center * (hi / center) ** factor
    return center - (center - lo) * factor, center + (hi - center) * factor


def zoom_axes(axes, factor, xc=None, active_ax=None, yc=None):
    """Applica lo zoom a tutti gli assi di una figura.

    - factor < 1  ->  zoom IN (l'intervallo si restringe)
    - factor > 1  ->  zoom OUT (l'intervallo si allarga)
    - xc: centro x da usare per tutti gli assi (utile per gli assi x
      condivisi, come nel Bode); se None si usa il centro di ciascun asse.
    - active_ax/yc: se forniti, l'asse 'active_ax' viene zoomato in y
      attorno a 'yc' (es. la posizione del mouse) invece che al suo centro.

    I limiti originali di TUTTI gli assi vengono letti prima di
    modificarne alcuno, per evitare un doppio zoom sugli assi x
    condivisi (sharex) tra i due sottografici del Bode.
    """
    snapshot = [(ax.get_xlim(), ax.get_ylim(), ax.get_xscale(), ax.get_yscale())
                for ax in axes]
    for ax, (xlim, ylim, xscale, yscale) in zip(axes, snapshot):
        xlo, xhi = xlim
        ylo, yhi = ylim
        cx = xc if xc is not None else (
            np.sqrt(xlo * xhi) if xscale == 'log' else (xlo + xhi) / 2)
        if active_ax is not None and ax is active_ax and yc is not None:
            cy = yc
        else:
            cy = np.sqrt(ylo * yhi) if yscale == 'log' else (ylo + yhi) / 2
        ax.set_xlim(*_zoom_range(xlo, xhi, cx, factor, xscale))
        ax.set_ylim(*_zoom_range(ylo, yhi, cy, factor, yscale))


def zoom_figure(fig, canvas, factor):
    """Zoom centrato sulla vista corrente (per i pulsanti +/-)."""
    zoom_axes(fig.axes, factor)
    canvas.draw_idle()


def reset_view(fig, canvas, original_limits):
    """Ripristina i limiti degli assi salvati al momento della creazione."""
    for ax, (xlim, ylim) in zip(fig.axes, original_limits):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    canvas.draw_idle()


def on_scroll(event, fig, canvas):
    """Zoom con la rotellina del mouse, centrato sul punto del cursore."""
    ax = event.inaxes
    if ax is None or event.xdata is None or event.ydata is None:
        return
    factor = (1 / 1.3) if event.button == 'up' else 1.3
    zoom_axes(fig.axes, factor, xc=event.xdata, active_ax=ax, yc=event.ydata)
    canvas.draw_idle()


# ============================================================================
# 3. COMPONENTI GUI RIUTILIZZABILI
# ============================================================================

class RealListFrame(ttk.LabelFrame):
    """Lista di costanti di tempo reali tau, per fattori (1+s*tau)."""

    def __init__(self, master, title):
        super().__init__(master, text=title, padding=8)
        self.values = []

        row = ttk.Frame(self)
        row.pack(fill='x')
        ttk.Label(row, text="tau [s]:").pack(side='left')
        self.entry = ttk.Entry(row, width=10)
        self.entry.pack(side='left', padx=4)
        ttk.Button(row, text="Aggiungi", command=self.add).pack(side='left', padx=2)
        ttk.Button(row, text="Rimuovi sel.", command=self.remove).pack(side='left', padx=2)

        self.listbox = tk.Listbox(self, height=4)
        self.listbox.pack(fill='x', pady=4)

    def add_value(self, val):
        self.values.append(val)
        self.listbox.insert('end', f"tau = {val:g} s")

    def add(self):
        txt = self.entry.get().strip().replace(',', '.')
        if not txt:
            return
        try:
            val = float(txt)
        except ValueError:
            messagebox.showerror("Errore", f"Valore non valido: {txt}")
            return
        if val == 0:
            messagebox.showerror("Errore", "tau non puo' essere 0.")
            return
        self.add_value(val)
        self.entry.delete(0, 'end')

    def remove(self):
        sel = self.listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        self.listbox.delete(idx)
        self.values.pop(idx)

    def clear(self):
        self.values = []
        self.listbox.delete(0, 'end')


class ComplexListFrame(ttk.LabelFrame):
    """Lista di coppie (wn, zeta), per fattori (1+2*zeta*s/wn+(s/wn)^2)."""

    def __init__(self, master, title):
        super().__init__(master, text=title, padding=8)
        self.values = []

        row = ttk.Frame(self)
        row.pack(fill='x')
        ttk.Label(row, text="wn [rad/s]:").pack(side='left')
        self.entry_wn = ttk.Entry(row, width=8)
        self.entry_wn.pack(side='left', padx=4)
        ttk.Label(row, text="zeta:").pack(side='left')
        self.entry_zeta = ttk.Entry(row, width=8)
        self.entry_zeta.pack(side='left', padx=4)
        ttk.Button(row, text="Aggiungi", command=self.add).pack(side='left', padx=2)
        ttk.Button(row, text="Rimuovi sel.", command=self.remove).pack(side='left', padx=2)

        self.listbox = tk.Listbox(self, height=4)
        self.listbox.pack(fill='x', pady=4)

    def add_value(self, wn, zeta):
        self.values.append((wn, zeta))
        self.listbox.insert('end', f"wn = {wn:g} rad/s,  zeta = {zeta:g}")

    def add(self):
        try:
            wn = float(self.entry_wn.get().strip().replace(',', '.'))
            zeta = float(self.entry_zeta.get().strip().replace(',', '.'))
        except ValueError:
            messagebox.showerror("Errore", "wn e zeta devono essere numeri.")
            return
        if wn <= 0:
            messagebox.showerror("Errore", "wn deve essere positivo.")
            return
        self.add_value(wn, zeta)
        self.entry_wn.delete(0, 'end')
        self.entry_zeta.delete(0, 'end')

    def remove(self):
        sel = self.listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        self.listbox.delete(idx)
        self.values.pop(idx)

    def clear(self):
        self.values = []
        self.listbox.delete(0, 'end')


# ============================================================================
# 4. APPLICAZIONE PRINCIPALE
# ============================================================================

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Diagrammi di Bode e Nyquist - Funzione di trasferimento")
        self.geometry("600x760")
        self._build_ui()

    def _set_spin(self, value):
        self.spin_g.delete(0, 'end')
        self.spin_g.insert(0, str(value))

    def _build_ui(self):
        pad = dict(padx=8, pady=6)

        top = ttk.LabelFrame(self, text="Parametri generali", padding=8)
        top.pack(fill='x', **pad)
        ttk.Label(top, text="Guadagno K:").grid(row=0, column=0, sticky='w')
        self.entry_K = ttk.Entry(top, width=10)
        self.entry_K.insert(0, "1")
        self.entry_K.grid(row=0, column=1, padx=6)

        ttk.Label(top, text="Tipo g (poli/zeri nell'origine):").grid(row=0, column=2, sticky='w')
        self.spin_g = ttk.Spinbox(top, from_=-5, to=5, width=5)
        self.spin_g.grid(row=0, column=3, padx=6)
        self._set_spin(0)

        self.frame_rz = RealListFrame(self, "Zeri reali  (1 + s*tau)")
        self.frame_rz.pack(fill='x', **pad)
        self.frame_rp = RealListFrame(self, "Poli reali  (1 + s*tau)")
        self.frame_rp.pack(fill='x', **pad)
        self.frame_cz = ComplexListFrame(self, "Zeri complessi  (1 + 2*zeta*s/wn + s^2/wn^2)")
        self.frame_cz.pack(fill='x', **pad)
        self.frame_cp = ComplexListFrame(self, "Poli complessi  (1 + 2*zeta*s/wn + s^2/wn^2)")
        self.frame_cp.pack(fill='x', **pad)

        freq = ttk.LabelFrame(self, text="Intervallo di frequenza", padding=8)
        freq.pack(fill='x', **pad)
        self.auto_range = tk.BooleanVar(value=True)
        ttk.Checkbutton(freq, text="Automatico", variable=self.auto_range,
                         command=self._toggle_range).grid(row=0, column=0, sticky='w')
        ttk.Label(freq, text="w min [rad/s]:").grid(row=1, column=0, sticky='w')
        self.entry_wmin = ttk.Entry(freq, width=10, state='disabled')
        self.entry_wmin.grid(row=1, column=1, padx=6)
        ttk.Label(freq, text="w max [rad/s]:").grid(row=1, column=2, sticky='w')
        self.entry_wmax = ttk.Entry(freq, width=10, state='disabled')
        self.entry_wmax.grid(row=1, column=3, padx=6)

        btns = ttk.Frame(self)
        btns.pack(fill='x', **pad)
        ttk.Button(btns, text="Carica esempio", command=self.load_example).pack(side='left', padx=4)
        ttk.Button(btns, text="Genera diagrammi", command=self.generate).pack(side='left', padx=4)
        ttk.Button(btns, text="Esci", command=self.destroy).pack(side='right', padx=4)

    def _toggle_range(self):
        state = 'disabled' if self.auto_range.get() else 'normal'
        self.entry_wmin.config(state=state)
        self.entry_wmax.config(state=state)

    def load_example(self):
        self.entry_K.delete(0, 'end')
        self.entry_K.insert(0, "100")
        self._set_spin(-2)
        self.frame_rz.clear()
        self.frame_rz.add_value(0.1)
        self.frame_rp.clear()
        self.frame_cz.clear()
        self.frame_cp.clear()
        self.frame_cp.add_value(50.0, 0.2)
        messagebox.showinfo(
            "Esempio caricato",
            "G(s) = 100*(1+0.1s) / [s^2*(1+0.4s/50+s^2/2500)]"
        )

    def _build_tf(self):
        try:
            K = float(self.entry_K.get().strip().replace(',', '.'))
            g = int(self.spin_g.get())
        except ValueError:
            messagebox.showerror("Errore", "Controllare i valori di K e g.")
            return None
        if K == 0:
            messagebox.showerror("Errore", "K non puo' essere 0.")
            return None
        return TransferFunction(
            K=K, g=g,
            real_zeros=self.frame_rz.values,
            real_poles=self.frame_rp.values,
            complex_zeros=self.frame_cz.values,
            complex_poles=self.frame_cp.values,
        )

    def _build_freq_vector(self, tf):
        if self.auto_range.get():
            return make_freq_vector(tf)
        try:
            wmin = float(self.entry_wmin.get().strip().replace(',', '.'))
            wmax = float(self.entry_wmax.get().strip().replace(',', '.'))
        except ValueError:
            messagebox.showerror("Errore", "Intervallo di frequenza non valido.")
            return None
        if wmin <= 0 or wmax <= wmin:
            messagebox.showerror("Errore", "Deve essere 0 < w_min < w_max.")
            return None
        return make_freq_vector(tf, wmin, wmax)

    def generate(self):
        tf = self._build_tf()
        if tf is None:
            return
        w = self._build_freq_vector(tf)
        if w is None:
            return

        self._open_bode_window(tf, w)
        self._open_nyquist_window(tf, w)

    def _open_bode_window(self, tf, w):
        win = tk.Toplevel(self)
        win.title("Diagrammi di Bode")
        fig = Figure(figsize=(8, 7))
        ax_m = fig.add_subplot(211)
        fig.add_subplot(212, sharex=ax_m)
        plot_bode(tf, w, fig)
        self._embed_figure(win, fig)

    def _open_nyquist_window(self, tf, w):
        win = tk.Toplevel(self)
        win.title("Diagramma di Nyquist")
        fig = Figure(figsize=(6.5, 6.5))
        fig.add_subplot(111)
        plot_nyquist(tf, w, fig)
        self._embed_figure(win, fig)

    @staticmethod
    def _embed_figure(win, fig):
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().pack(side='top', fill='both', expand=True)

        controls = ttk.Frame(win)
        controls.pack(side='bottom', fill='x')

        toolbar = NavigationToolbar2Tk(canvas, controls, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side='left')

        # --- memorizza i limiti iniziali di ciascun asse, per "Reset vista"
        original_limits = {ax: (ax.get_xlim(), ax.get_ylim()) for ax in fig.axes}

        def zoom_axis(get_lim, set_lim, scale, factor, center=None):
            lo, hi = get_lim()
            if center is None:
                center = np.sqrt(lo * hi) if scale == 'log' else (lo + hi) / 2
            if scale == 'log':
                set_lim(center * (lo / center) ** factor, center * (hi / center) ** factor)
            else:
                set_lim(center + (lo - center) * factor, center + (hi - center) * factor)

        def zoom_all(factor):
            """Zoom centrato sul punto medio corrente di ciascun asse (bottoni +/-)."""
            for ax in fig.axes:
                zoom_axis(ax.get_xlim, ax.set_xlim, ax.get_xscale(), factor)
                zoom_axis(ax.get_ylim, ax.set_ylim, ax.get_yscale(), factor)
            canvas.draw_idle()

        def reset_view():
            for ax, (xlim, ylim) in original_limits.items():
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
            canvas.draw_idle()

        def on_scroll(event):
            """Zoom con la rotellina del mouse, centrato sul cursore."""
            ax = event.inaxes
            if ax is None or event.xdata is None or event.ydata is None:
                return
            factor = (1 / 1.2) if event.button == 'up' else 1.2
            zoom_axis(ax.get_xlim, ax.set_xlim, ax.get_xscale(), factor, center=event.xdata)
            zoom_axis(ax.get_ylim, ax.set_ylim, ax.get_yscale(), factor, center=event.ydata)
            canvas.draw_idle()

        zoom_frame = ttk.Frame(controls)
        zoom_frame.pack(side='right', padx=6)
        ttk.Button(zoom_frame, text="Zoom +", width=8,
                   command=lambda: zoom_all(1 / 1.3)).pack(side='left', padx=2)
        ttk.Button(zoom_frame, text="Zoom -", width=8,
                   command=lambda: zoom_all(1.3)).pack(side='left', padx=2)
        ttk.Button(zoom_frame, text="Reset vista", width=10,
                   command=reset_view).pack(side='left', padx=2)

        canvas.mpl_connect('scroll_event', on_scroll)
        win.protocol("WM_DELETE_WINDOW", win.destroy)


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
