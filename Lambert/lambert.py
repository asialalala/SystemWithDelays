import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

def solve_dde_lambert(a, b, tau, t_max, n_branches=5):
    """
    Wyznacza przybliżone rozwiązanie DDE: x'(t) = ax(t) + bx(t-tau)
    za pomocą n_branches gałęzi funkcji Lamberta.
    """
    t = np.linspace(0, t_max, 500)
    x_t = np.zeros_like(t, dtype=complex)
    
    # Obliczamy s_k dla każdej gałęzi k
    # Formuła: s_k = a + (1/tau) * W_k(b * tau * exp(-a * tau))
    for k in range(-n_branches, n_branches + 1):
        # Argument funkcji W Lamberta
        arg = b * tau * np.exp(-a * tau)
        
        # Wyznaczenie s_k dla gałęzi k
        s_k = a + (1/tau) * lambertw(arg, k=k)
        
        # Zakładamy dla uproszczenia stały warunek początkowy (C_k = 1)
        # W rzeczywistości C_k zależą od historii układu x(t) dla t < 0
        x_t += np.exp(s_k * t)

    return t, x_t.real

# Parametry układu
a = -0.5
b = -1.0
tau = 1.0  # Opóźnienie

t, response = solve_dde_lambert(a, b, tau, t_max=20)

# Wizualizacja
plt.figure(figsize=(10, 5))
plt.plot(t, response, label=f'Suma {11} gałęzi $W_k$', color='blue')
plt.title(f"Odpowiedź układu z opóźnieniem $\\tau={tau}$")
plt.xlabel("Czas [t]")
plt.ylabel("x(t)")
plt.grid(True)
plt.legend()
plt.show()