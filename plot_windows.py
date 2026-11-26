"""Plot windows and navigation controls."""

import tkinter as tk
from tkinter import ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from plots import plot_bode, plot_nyquist


def open_bode_window(parent, tf, w):
    win = tk.Toplevel(parent)
    win.title("Bode Plots")
    fig = Figure(figsize=(8, 7))
    ax_m = fig.add_subplot(211)
    fig.add_subplot(212, sharex=ax_m)
    plot_bode(tf, w, fig)
    _embed_figure(win, fig)

def open_nyquist_window(parent, tf, w):
    win = tk.Toplevel(parent)
    win.title("Nyquist Plot")
    fig = Figure(figsize=(6.5, 6.5))
    fig.add_subplot(111)
    plot_nyquist(tf, w, fig)
    _embed_figure(win, fig)


def _embed_figure(win, fig):
    canvas = FigureCanvasTkAgg(fig, master=win)
    canvas.draw()
    canvas.get_tk_widget().pack(side='top', fill='both', expand=True)

    controls = ttk.Frame(win)
    controls.pack(side='bottom', fill='x')

    toolbar = NavigationToolbar2Tk(canvas, controls, pack_toolbar=False)
    toolbar.update()
    toolbar.pack(side='left')

    # --- save the initial axis limits for the reset button
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
        """Zoom around the current midpoint of each axis."""
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
        """Zoom around the pointer when the mouse wheel turns."""
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
    ttk.Button(zoom_frame, text="Reset view", width=10,
               command=reset_view).pack(side='left', padx=2)

    canvas.mpl_connect('scroll_event', on_scroll)
    win.protocol("WM_DELETE_WINDOW", win.destroy)
