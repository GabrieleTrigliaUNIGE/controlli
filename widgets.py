"""Input widgets for poles and zeros."""

import tkinter as tk
from tkinter import ttk, messagebox


class RealListFrame(ttk.LabelFrame):
    """Edit time constants for real first-order factors."""

    def __init__(self, master, title):
        super().__init__(master, text=title, padding=8)
        self.values = []

        row = ttk.Frame(self)
        row.pack(fill='x')
        ttk.Label(row, text="tau [s]:").pack(side='left')
        self.entry = ttk.Entry(row, width=10)
        self.entry.pack(side='left', padx=4)
        ttk.Button(row, text="Add", command=self.add).pack(side='left', padx=2)
        ttk.Button(row, text="Remove selected", command=self.remove).pack(side='left', padx=2)

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
            messagebox.showerror("Error", f"Invalid value: {txt}")
            return
        if val == 0:
            messagebox.showerror("Error", "tau must be nonzero.")
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
    """Edit natural frequency and damping ratio pairs for second-order factors."""

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
        ttk.Button(row, text="Add", command=self.add).pack(side='left', padx=2)
        ttk.Button(row, text="Remove selected", command=self.remove).pack(side='left', padx=2)

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
            messagebox.showerror("Error", "wn and zeta must be numbers.")
            return
        if wn <= 0:
            messagebox.showerror("Error", "wn must be positive.")
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
