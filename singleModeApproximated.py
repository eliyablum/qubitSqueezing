# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 10:44:11 2024

@author: Eliya
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt

#%% Physical parameters (all frequencies in GHz and time in ns)

Wrabi = 40e-3 * 2 * np.pi        # Rabi frequency [GHz]
chi = -0.05e-3 * 2 * np.pi       # Dispersive shift [GHz]
chi_a = -0.06e-3 * 2 * np.pi     # Dispersive shift for mode a [GHz]
chi_b = -0.04e-3 * 2 * np.pi     # Dispersive shift for mode b [GHz]
deltaW = 1.6e-3 * 2 * np.pi      # Detuning of sidebands from Rabi frequency [GHz]
Wsb1 = Wrabi + deltaW            # Upper sideband frequency [GHz]
Wsb2 = Wrabi - deltaW            # Lower sideband frequency [GHz]
gamma = 0.16e-3 * 2 * np.pi      # Coupling rate [GHz]

# Complex amplitude for displacement
abar = -gamma * 2 / chi * np.exp(1j * np.pi / 4)  # Optimal displacement amplitude

#%% Time parameters
tf = 20000                                        # Final time in ns
times = np.linspace(0, tf, int(tf + 1))           # Time array from 0 to tf [ns]

#%% Basis states
g = basis(2, 1)                                   # Ground state |g⟩
e = basis(2, 0)                                   # Excited state |e⟩

#%% Hilbert space truncation
N = 400                                           # Photon number cutoff
a = destroy(N)                                    # Annihilation operator for mode a
b = destroy(N)                                    # Annihilation operator for mode b (unused here)
vac = basis(N, 0)                                 # Vacuum state of cavity

#%% Initial states
init_displaced = tensor(g, vac).unit()            # Initial state: |g⟩ ⊗ |0⟩
init_RWA = tensor(g + e, vac)                     # Superposition init state (not normalized)
init_ideal = tensor(g + e, vac)                   # Ideal evolution initial state (not normalized)

#%% Quadrature operators (for squeezing observables)
X = (a.dag() + a) / np.sqrt(2)                    # Position quadrature
P = 1j * (a.dag() - a) / np.sqrt(2)               # Momentum quadrature

# Symmetric/antisymmetric combinations for potential two-mode squeezing
Xp = (tensor(qeye(2), X, qeye(N)) + tensor(qeye(2), qeye(N), X)) / np.sqrt(2)
Xm = (tensor(qeye(2), X, qeye(N)) - tensor(qeye(2), qeye(N), X)) / np.sqrt(2)
Pp = (tensor(qeye(2), P, qeye(N)) + tensor(qeye(2), qeye(N), P)) / np.sqrt(2)
Pm = (tensor(qeye(2), P, qeye(N)) - tensor(qeye(2), qeye(N), P)) / np.sqrt(2)

# Observables list for expectation value calculations
e_ops = [Xp, Pp, Xm, Pm, 
         Xp**2, Xm**2, Pp**2, Pm**2]

#%% Hamiltonian definitions

is_high_order = False  # Flag to include higher-order corrections

# Leading-order sideband-induced squeezing Hamiltonian (effective model)
H = chi**2 / (8 * deltaW) * (
    abar.conjugate()**2 * tensor(sigmaz(), a * a) +
    abar**2 * tensor(sigmaz(), a.dag() * a.dag())
)

# Optional higher-order correction terms (nonlinear squeezing contributions)
H3 = 1j * chi**3 / (8 * deltaW**2) * (
    abar * abar.conjugate()**2 * tensor(sigmam(), a) -
    abar**2 * abar.conjugate() * tensor(sigmap(), a.dag()) +
    abar ** 3 * tensor(sigmam(), a.dag()**3) -
    abar.conjugate()**3 * tensor(sigmap(), a**3) +
    abar * abar.conjugate()**2 * tensor(sigmam(), a.dag() * a**2) -
    abar**2 * abar.conjugate() * tensor(sigmap(), a.dag()**2 * a)
)

if is_high_order:
    H += H3

#%% Observables for measurement
e_ops = [tensor(sig, qeye(N)) for sig in [sigmax(), sigmay(), sigmaz()]]  # Qubit Pauli operators
e_ops.append(tensor(qeye(2), a.dag() * a))                # Cavity photon number
e_ops.append(tensor(qeye(2), basis(N, N-1).proj()))       # Projection on Fock state N-1
e_ops.append(tensor(sigmaz(), qeye(N)))                   # σ_z

# Projectors on qubit ground and excited states
e_ops.append(tensor(g.unit().proj(), qeye(N)))
e_ops.append(tensor(e.unit().proj(), qeye(N)))

# Cavity quadrature moments
e_ops.append(tensor(qeye(2), X**2))
e_ops.append(tensor(qeye(2), P**2))
e_ops.append(tensor(qeye(2), P * X))
e_ops.append(tensor(qeye(2), X * P))

# Conditional quadrature moments when qubit is in ground state
e_ops.append(tensor(g.unit().proj(), X**2))
e_ops.append(tensor(g.unit().proj(), P**2))
e_ops.append(tensor(g.unit().proj(), P * X))
e_ops.append(tensor(g.unit().proj(), X * P))

# Conditional quadrature moments when qubit is in excited state
e_ops.append(tensor(e.unit().proj(), X**2))
e_ops.append(tensor(e.unit().proj(), P**2))
e_ops.append(tensor(e.unit().proj(), P * X))
e_ops.append(tensor(e.unit().proj(), X * P))

#%% Simulation (unitary evolution)
res_ideal = mesolve(
    H,
    init_ideal,
    times,
    e_ops=e_ops,
    c_ops=None,  # No dissipation
    options=Options(store_states=False, max_step=1),
    progress_bar=True
)

#%% Squeezing calculation (in dB) conditioned on qubit state

# Ground state population and moments
g_res = res_ideal.expect[-14]
XX_g_res = res_ideal.expect[-8]
YY_g_res = res_ideal.expect[-7]

# Excited state population and moments
e_res = res_ideal.expect[-13]
XX_e_res = res_ideal.expect[-4]
YY_e_res = res_ideal.expect[-3]

# Compute squeezing in dB for X and P quadratures
# Compared to vacuum variance of 0.5 (→ factor of 2 here)
Sq_g_dB_XX = -10 * np.log10(XX_g_res / g_res * 2)
Sq_g_dB_YY = -10 * np.log10(YY_g_res / g_res * 2)
Sq_e_dB_XX = -10 * np.log10(XX_e_res / e_res * 2)
Sq_e_dB_YY = -10 * np.log10(YY_e_res / e_res * 2)

#%% Plot squeezing vs time

plt.figure()
plt.plot(times, Sq_g_dB_XX, 'b-', label='XXg')  # Ground state, X quadrature
plt.plot(times, Sq_g_dB_YY, 'r-', label='YYg')  # Ground state, P quadrature
plt.plot(times, Sq_e_dB_XX, 'b--', label='XXe') # Excited state, X quadrature
plt.plot(times, Sq_e_dB_YY, 'r--', label='YYe') # Excited state, P quadrature
plt.ylabel('Squeezed amplitude [dB]')
plt.xlabel('Time [ns]')
plt.legend()

#%% Plot qubit expectation values (Bloch vector components)

fig, axs = plt.subplots(3, 1, sharex=True)
labels = ['X', 'Y', 'Z']

for ax, ex, lb in zip(axs, res_ideal.expect[:3], labels):
    plt.sca(ax)
    plt.plot(times, ex, label='Ideal')
    plt.ylabel(lb)
    plt.legend()
