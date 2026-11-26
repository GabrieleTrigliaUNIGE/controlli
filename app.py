"""Main window and parameter input."""

import tkinter as tk
from tkinter import ttk, messagebox

from widgets import ComplexListFrame, RealListFrame
from plot_windows import open_bode_window, open_nyquist_window
from transfer_function import TransferFunction, make_freq_vector


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Bode and Nyquist Plots - Transfer Function")
        self.geometry("600x760")
        self._build_ui()

    def _set_spin(self, value):
        self.spin_g.delete(0, 'end')
        self.spin_g.insert(0, str(value))

    def _build_ui(self):
        pad = dict(padx=8, pady=6)

        top = ttk.LabelFrame(self, text="General parameters", padding=8)
        top.pack(fill='x', **pad)
        ttk.Label(top, text="Gain K:").grid(row=0, column=0, sticky='w')
        self.entry_K = ttk.Entry(top, width=10)
        self.entry_K.insert(0, "1")
        self.entry_K.grid(row=0, column=1, padx=6)

        ttk.Label(top, text="Exponent g in s^g:").grid(row=0, column=2, sticky='w')
        self.spin_g = ttk.Spinbox(top, from_=-5, to=5, width=5)
        self.spin_g.grid(row=0, column=3, padx=6)
        self._set_spin(0)

        self.frame_rz = RealListFrame(self, "Real zeros  (1 + s*tau)")
        self.frame_rz.pack(fill='x', **pad)
        self.frame_rp = RealListFrame(self, "Real poles  (1 + s*tau)")
        self.frame_rp.pack(fill='x', **pad)
        self.frame_cz = ComplexListFrame(self, "Complex zeros  (1 + 2*zeta*s/wn + s^2/wn^2)")
        self.frame_cz.pack(fill='x', **pad)
        self.frame_cp = ComplexListFrame(self, "Complex poles  (1 + 2*zeta*s/wn + s^2/wn^2)")
        self.frame_cp.pack(fill='x', **pad)

        freq = ttk.LabelFrame(self, text="Frequency range", padding=8)
        freq.pack(fill='x', **pad)
        self.auto_range = tk.BooleanVar(value=True)
        ttk.Checkbutton(freq, text="Automatic", variable=self.auto_range,
                         command=self._toggle_range).grid(row=0, column=0, sticky='w')
        ttk.Label(freq, text="w min [rad/s]:").grid(row=1, column=0, sticky='w')
        self.entry_wmin = ttk.Entry(freq, width=10, state='disabled')
        self.entry_wmin.grid(row=1, column=1, padx=6)
        ttk.Label(freq, text="w max [rad/s]:").grid(row=1, column=2, sticky='w')
        self.entry_wmax = ttk.Entry(freq, width=10, state='disabled')
        self.entry_wmax.grid(row=1, column=3, padx=6)

        btns = ttk.Frame(self)
        btns.pack(fill='x', **pad)
        ttk.Button(btns, text="Load example", command=self.load_example).pack(side='left', padx=4)
        ttk.Button(btns, text="Generate plots", command=self.generate).pack(side='left', padx=4)
        ttk.Button(btns, text="Exit", command=self.destroy).pack(side='right', padx=4)

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
            "Example loaded",
            "G(s) = 100*(1+0.1s) / [s^2*(1+0.4s/50+s^2/2500)]"
        )

    def _build_tf(self):
        try:
            K = float(self.entry_K.get().strip().replace(',', '.'))
            g = int(self.spin_g.get())
        except ValueError:
            messagebox.showerror("Error", "Check the values of K and g.")
            return None
        if K == 0:
            messagebox.showerror("Error", "K must be nonzero.")
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
            messagebox.showerror("Error", "Invalid frequency range.")
            return None
        if wmin <= 0 or wmax <= wmin:
            messagebox.showerror("Error", "The range must satisfy 0 < w_min < w_max.")
            return None
        return make_freq_vector(tf, wmin, wmax)

    def generate(self):
        tf = self._build_tf()
        if tf is None:
            return
        w = self._build_freq_vector(tf)
        if w is None:
            return

        open_bode_window(self, tf, w)
        open_nyquist_window(self, tf, w)
