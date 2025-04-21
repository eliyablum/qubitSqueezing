# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 10:44:11 2024

@author: Eliya
"""

from qutip import *                    # QuTiP: Quantum Toolbox in Python
import numpy as np                    # NumPy for numerical operations
import matplotlib.pyplot as plt       # Matplotlib for plotting

#%% PARAMETERS
Wrabi = 40e-3 * 2 * np.pi             # Rabi frequency [GHz]
chi = -0.05e-3 * 2 * np.pi            # Dispersive coupling strength [GHz]
deltaW = 1.6e-3 * 2 * np.pi           # Sideband frequency detuning [GHz]
Wsb1 = Wrabi + deltaW                 # Upper sideband frequency
Wsb2 = Wrabi - deltaW                 # Lower sideband frequency
gamma = 0.16e-3 * 2 * np.pi           # Effective rate [GHz]

abar = 1j * np.sqrt(-1j) * gamma * 2 / chi  # Complex amplitude (displacement)

tf = 20000                                  # Final simulation time [ns]
times = np.linspace(0, tf, int(tf + 1))     # Time array
N = 400                                     # Truncation for Fock space (cavity dimension)

#%% STATES AND OPERATORS
g = basis(2, 1)                  # Qubit ground state |g>
e = basis(2, 0)                  # Qubit excited state |e>

a = destroy(N)                   # Annihilation operator for the cavity
vac = basis(N, 0)                # Vacuum state of the cavity

init_displaced = tensor(g + e, vac).unit()  # Initial state: superposition of |g> and |e> with vacuum

# Define quadratures X and P (cavity operators)
X = (a.dag() + a) / np.sqrt(2)   # Position-like quadrature
P = 1j * (a.dag() - a) / np.sqrt(2)  # Momentum-like quadrature

# Expectation value operators for measurement
e_ops = [tensor(sig, qeye(N)) for sig in [sigmax(), sigmay(), sigmaz()]]  # Qubit Pauli operators
e_ops.append(tensor(qeye(2), a.dag() * a))                                # Photon number
e_ops.append(tensor(qeye(2), basis(N, N - 1).proj()))                     # Projector onto highest Fock state
e_ops.append(tensor(sigmaz(), qeye(N)))                                   # σ_z operator
e_ops.append(tensor(g.proj(), qeye(N)))                                   # Projector onto |g>
e_ops.append(tensor(e.proj(), qeye(N)))                                   # Projector onto |e>

# Quadrature variances and correlations
e_ops.append(tensor(qeye(2), X**2))      # ⟨X²⟩
e_ops.append(tensor(qeye(2), P**2))      # ⟨P²⟩
e_ops.append(tensor(qeye(2), P * X))     # ⟨PX⟩
e_ops.append(tensor(qeye(2), X * P))     # ⟨XP⟩

#%% HAMILTONIAN

# Bare system Hamiltonian
H = Wrabi / 2 * tensor(sigmax(), qeye(N)) + chi / 2 * tensor(sigmaz(), a.dag() * a)

# Time-dependent interaction terms
interaction = -chi / 2 * abar.conjugate() * (np.sin(times * Wsb1) + 1j * np.cos(times * Wsb2))
time_dependent_ss = chi / 2 * np.abs(abar)**2 * np.sin(2 * Wrabi * times) * np.sin(2 * deltaW * times)

# Full time-dependent Hamiltonian as a list for mesolve
H_tot = [
    H,
    [tensor(sigmaz(), a), interaction],
    [tensor(sigmaz(), a.dag()), interaction.conjugate()],
    [tensor(sigmaz(), qeye(N)), time_dependent_ss]
]

#%% SIMULATION

# Run the master equation solver (no dissipation here, c_ops = None)
res_displaced = mesolve(
    H_tot, 
    init_displaced, 
    times, 
    e_ops=e_ops, 
    c_ops=None, 
    options=Options(store_states=False), 
    progress_bar=True
)

#%% PLOT THE RESULTS

# Plot qubit Bloch sphere components over time
fig, axs = plt.subplots(3, 1, sharex=True)
labels = ['X', 'Y', 'Z']

for ax, ex, lb in zip(axs, res_displaced.expect[:3], labels):
    plt.sca(ax)
    plt.plot(times, ex, label='Full')
    plt.ylabel(lb)
    plt.legend()

# Compute and plot quadrature squeezing in dB
XX_res = res_displaced.expect[-4]     # ⟨X²⟩
YY_res = res_displaced.expect[-3]     # ⟨P²⟩

# Squeezing in decibels: -10 * log10(2⟨X²⟩) etc.
Sq_dB_XX = -10 * np.log10(XX_res * 2)
Sq_dB_YY = -10 * np.log10(YY_res * 2)

plt.figure()
plt.plot(times, Sq_dB_XX, 'b-', label='XX full')
plt.plot(times, Sq_dB_YY, 'r-', label='YY full')
plt.ylabel('Squeezed amplitude [dB]')
plt.xlabel('Time [ns]')
plt.legend()
