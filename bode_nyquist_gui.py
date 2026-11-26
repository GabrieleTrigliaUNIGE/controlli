#!/usr/bin/env python3
"""Launch the Bode and Nyquist plot application."""

try:
    from app import App
except ModuleNotFoundError as exc:
    raise SystemExit(f"Missing module: {exc.name}") from exc


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
