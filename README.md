# Bode & Nyquist

A Python desktop tool for studying continuous-time transfer functions. Enter the gain, poles, and zeros in factored form to compare exact and asymptotic Bode plots and explore the Nyquist curve around the critical point (−1, 0).

## Run on Ubuntu

Install Python 3, Tkinter, NumPy, and Matplotlib:

```bash
sudo apt update
sudo apt install python3 python3-tk python3-numpy python3-matplotlib
```

Download and extract the source ZIP, or clone the repository:

```bash
git clone https://github.com/GabrieleTrigliaUNIGE/controlli.git
cd controlli
```

From the directory containing the Python files, start the application:

```bash
python3 bode_nyquist_gui.py
```

## Plot a transfer function

Choose **Load example** and then **Generate plots** for a first run. To enter your own system, use the three tabs:

| Tab | Inputs |
| --- | --- |
| Parameters | Gain K, exponent g in sᵍ, and automatic or manual frequency range in rad/s. Use g = −1 for one integrator. |
| Zeros | Numerator factors. Enter a time constant τ in seconds for a first-order factor, or natural frequency ωₙ in rad/s and damping ratio ζ for a second-order factor. |
| Poles | Denominator factors, using the same conventions as zeros. |

For example, to plot **G(s) = 10(1 + s) / [s(1 + 0.3s)(1 + 0.01s)]**, start with an empty form. Set K to 10 and g to −1. Add one real zero with τ = 1 and two real poles with τ = 0.3 and τ = 0.01. Select **Generate plots**.

The Bode window displays magnitude in dB and phase in degrees. The dashed lines show the asymptotic approximation. The Nyquist window displays the sampled positive-frequency response and its mirrored branch. For poles on the imaginary axis, the plot omits the contour detours required for a full Nyquist stability analysis.

Use the mouse wheel to zoom around the pointer, or choose **Reset view** to restore the initial limits.

## Source files

| File | Responsibility |
| --- | --- |
| `bode_nyquist_gui.py` | Application entry point |
| `app.py` | Main window, parameter input, and built-in example |
| `widgets.py` | Pole and zero editors |
| `transfer_function.py` | Frequency response, Bode approximations, and frequency sampling |
| `plots.py` | Bode and Nyquist drawing routines |
| `plot_windows.py` | Plot windows, navigation toolbar, and zoom controls |
