import control as ct
import matplotlib.pyplot as plt
import numpy as np

def plot_risposta_frequenza():
    print("--- Generatore Diagrammi di Bode e Nyquist ---")
    print("Inserisci i coefficienti dei polinomi separati da spazio.")
    print("Esempio: per s^2 + 2s + 10 inserisci -> 1 2 10")
    print("-" * 45)

    try:
        # 1. Input dei coefficienti
        num_str = input("Inserisci coefficienti NUMERATORE: ")
        den_str = input("Inserisci coefficienti DENOMINATORE: ")

        # Conversione stringa -> lista di float
        num = [float(x) for x in num_str.split()]
        den = [float(x) for x in den_str.split()]

        # 2. Creazione della Funzione di Trasferimento
        sys = ct.tf(num, den)
        
        print("\nFunzione di Trasferimento rilevata:")
        print(sys)

        # 3. Diagramma di Bode
        # Crea una nuova figura per Bode
        plt.figure(figsize=(10, 7))
        # plt.title("Diagramma di Bode")
        # La funzione bode_plot calcola e plotta. dB=True converte il modulo in decibel.
        # deg=True converte la fase in gradi.
        ct.bode_plot(sys, dB=True, deg=True, display_margins=True, grid=True)
        plt.suptitle("")

        # 4. Diagramma di Nyquist
        # Crea una nuova figura separata per Nyquist
        plt.figure(figsize=(8, 8))
        ct.nyquist_plot(sys, arrows=True)
        plt.title("Diagramma di Nyquist")
        # --- MODIFICA PER RIMUOVERE L'OVERLAP ---
        # 1. Recuperiamo gli assi correnti
        ax = plt.gca()
        
        # 2. Impostiamo il titolo direttamente sugli assi (sovrascrive quello di default "Nyquist plot for...")
        ax.set_title("Diagramma di Nyquist")
        
        # 3. Per sicurezza, rimuoviamo l'eventuale titolo superiore (suptitle) se la libreria lo ha creato
        plt.suptitle("")
        plt.grid(True)
        
        # Aggiungo un cerchio unitario per riferimento (utile per stabilità)
        theta = np.linspace(0, 2*np.pi, 100)
        plt.plot(np.cos(theta), np.sin(theta), 'r--', label='Cerchio Unitario (-1, j0)')
        plt.legend()

        # Mostra i grafici
        print("\nGenerazione grafici in corso...")
        # plt.show()

        # ... codice precedente ...

        # Invece di plt.show() usa:
        print("\nSalvataggio grafici in corso...")
        plt.figure(1).savefig('diagramma_bode.png')
        plt.figure(2).savefig('diagramma_nyquist.png')
        print("Grafici salvati come 'diagramma_bode.png' e 'diagramma_nyquist.png'")
        
        # plt.show()  <-- Commenta questa riga se non hai un display server

    except ValueError:
        print("\nERRORE: Assicurati di inserire solo numeri separati da spazi.")
    except Exception as e:
        print(f"\nSi è verificato un errore: {e}")

if __name__ == "__main__":
    plot_risposta_frequenza()