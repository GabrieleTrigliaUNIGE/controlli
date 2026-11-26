"""Main window and parameter input."""

import tkinter as tk
from tkinter import ttk, messagebox

from widgets import ComplexListFrame, RealListFrame
from plot_windows import open_bode_window, open_nyquist_window
from transfer_function import TransferFunction, make_freq_vector


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Bode & Nyquist")
        self._configure_style()
        self._build_ui()
        # Measure the active theme and font before choosing the window size.
        self.update_idletasks()
        self.minsize(self.winfo_reqwidth(), self.winfo_reqheight())

    def _configure_style(self):
        style = ttk.Style(self)
        style.theme_use("clam")
        self.configure(background="#f4f7fa")
        style.configure(".", background="#f4f7fa", foreground="#20354a")
        style.configure("TButton", padding=(9, 5))
        style.configure("Primary.TButton", background="#17665f", foreground="white")
        style.map("Primary.TButton", background=[("pressed", "#104640"),
                                                  ("active", "#1d7870")])
        style.configure("TNotebook", borderwidth=0)
        style.configure("TNotebook.Tab", padding=(14, 6))
        style.map("TNotebook.Tab", foreground=[("selected", "#125b55")],
                  background=[("selected", "#e1efed")])
        style.configure("TLabelframe", bordercolor="#cbd7df", relief="solid")
        style.configure("TLabelframe.Label", foreground="#284761")
        self.option_add("*Listbox.relief", "flat")
        self.option_add("*Listbox.highlightThickness", 1)
        self.option_add("*Listbox.highlightBackground", "#cbd7df")
        self.option_add("*Listbox.selectBackground", "#17665f")

    def _set_spin(self, value):
        self.spin_g.delete(0, 'end')
        self.spin_g.insert(0, str(value))

    def _build_ui(self):
        # Reserve the action bar before allocating space to the editor.
        actions = ttk.Frame(self, padding=(16, 10))
        actions.pack(side="bottom", fill="x")
        ttk.Button(actions, text="Load example", command=self.load_example).pack(side="left")
        ttk.Button(actions, text="Generate plots", style="Primary.TButton",
                   command=self.generate).pack(side="right")
        ttk.Button(actions, text="Exit", command=self.destroy).pack(side="right", padx=8)
        ttk.Separator(self).pack(side="bottom", fill="x")

        tabs = ttk.Notebook(self)
        tabs.pack(fill="both", expand=True, padx=16, pady=(14, 8))
        settings = ttk.Frame(tabs, padding=10)
        zeros = ttk.Frame(tabs, padding=10)
        poles = ttk.Frame(tabs, padding=10)
        for page, title in ((settings, "Parameters"), (zeros, "Zeros"), (poles, "Poles")):
            tabs.add(page, text=title)
        tabs.enable_traversal()

        top = ttk.LabelFrame(settings, text="General parameters", padding=10)
        top.pack(fill="x", pady=(0, 8))
        ttk.Label(top, text="Gain K:").grid(row=0, column=0, sticky='w')
        self.entry_K = ttk.Entry(top, width=10)
        self.entry_K.insert(0, "1")
        self.entry_K.grid(row=0, column=1, padx=6)

        ttk.Label(top, text="Exponent g in s^g:").grid(row=0, column=2, sticky='w')
        self.spin_g = ttk.Spinbox(top, from_=-5, to=5, width=5)
        self.spin_g.grid(row=0, column=3, padx=6)
        self._set_spin(0)

        self.frame_rz = RealListFrame(zeros, "Real zeros  (1 + s*tau)")
        self.frame_rz.pack(fill='x', pady=(0, 6))
        self.frame_rp = RealListFrame(poles, "Real poles  (1 + s*tau)")
        self.frame_rp.pack(fill='x', pady=(0, 6))
        self.frame_cz = ComplexListFrame(zeros, "Complex zeros  (1 + 2*zeta*s/wn + s^2/wn^2)")
        self.frame_cz.pack(fill='x', pady=(0, 6))
        self.frame_cp = ComplexListFrame(poles, "Complex poles  (1 + 2*zeta*s/wn + s^2/wn^2)")
        self.frame_cp.pack(fill='x', pady=(0, 6))

        freq = ttk.LabelFrame(settings, text="Frequency range", padding=10)
        freq.pack(fill='x', pady=(0, 6))
        self.auto_range = tk.BooleanVar(value=True)
        ttk.Checkbutton(freq, text="Automatic", variable=self.auto_range,
                         command=self._toggle_range).grid(row=0, column=0, sticky='w')
        ttk.Label(freq, text="w min [rad/s]:").grid(row=1, column=0, sticky='w')
        self.entry_wmin = ttk.Entry(freq, width=10, state='disabled')
        self.entry_wmin.grid(row=1, column=1, padx=6)
        ttk.Label(freq, text="w max [rad/s]:").grid(row=1, column=2, sticky='w')
        self.entry_wmax = ttk.Entry(freq, width=10, state='disabled')
        self.entry_wmax.grid(row=1, column=3, padx=6)

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
