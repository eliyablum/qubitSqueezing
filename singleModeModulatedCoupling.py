# -*- coding: utf-8 -*-
"""
Created on Wed Jun 5 09:32:18 2024
@author: Shay

All frequencies are in GHz, and time is in ns.
Gamma is a coupling rate.
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt

#%% Define physical parameters
Wrabi = 40e-3 * 2 * np.pi             # Rabi frequency [GHz]
chi = -0.05e-3 * 2 * np.pi            # Dispersive coupling [GHz]
deltaW = 1.6e-3 * 2 * np.pi           # Sideband detuning [GHz]
Wsb1 = Wrabi + deltaW                 # Upper sideband frequency [GHz]
Wsb2 = Wrabi - deltaW                 # Lower sideband frequency [GHz]
gamma = 0.16e-3 * 2 * np.pi           # Coupling rate [GHz]

# Displacement amplitude (used for squeezing), complex-valued
abar = 1j * np.sqrt(-1j) * gamma * 2 / chi

# Simulation time in ns
tf = 30000
times = np.linspace(0, tf, int(tf + 1))  # Time array [ns]

#%% Define quantum states and operators

g = basis(2, 1)              # Ground state |g⟩
e = basis(2, 0)              # Excited state |e⟩

N = 1000                    # Cavity photon number cutoff
a = destroy(N)              # Annihilation operator
vac = basis(N, 0)           # Vacuum state

# Initial state: superposition (|g⟩ + |e⟩) ⊗ |0⟩
init = tensor(g + e, vac).unit()

# Cavity quadrature operators
X = (a.dag() + a) / np.sqrt(2)       # Position quadrature
P = 1j * (a.dag() - a) / np.sqrt(2)  # Momentum quadrature

#%% Measurement observables (expectation values)
e_ops = [tensor(sig, qeye(N)) for sig in [sigmax(), sigmay(), sigmaz()]]  # Qubit Pauli operators
e_ops.append(tensor(qeye(2), a.dag() * a))                 # Photon number ⟨a†a⟩
e_ops.append(tensor(qeye(2), basis(N, N - 1).proj()))      # Projection on Fock state |N-1⟩
e_ops.append(tensor(sigmaz(), qeye(N)))                    # σ_z ⊗ I

# Projectors on qubit states
e_ops.append(tensor(g.proj(), qeye(N)))
e_ops.append(tensor(e.proj(), qeye(N)))

# Cavity quadrature moments
e_ops.append(tensor(qeye(2), X ** 2))
e_ops.append(tensor(qeye(2), P ** 2))
e_ops.append(tensor(qeye(2), P * X))
e_ops.append(tensor(qeye(2), X * P))

#%% Time-dependent Hamiltonian (beyond RWA)

# Static part of Hamiltonian (here zero, but placeholder for generality)
H = tensor(sigmax(), qeye(N)) * 0

# Complex time-dependent coefficients:
interaction1 = chi / 2 * abar * 1j * (np.exp(1j * deltaW * times) + np.exp(-1j * deltaW * times)) / 2
interaction2 = -chi / 2 * abar * (np.exp(1j * deltaW * times) - np.exp(-1j * deltaW * times)) / (2j)

# Corresponding operators
interaction_op1 = tensor(sigmap(), a.dag())     # σ⁺ ⊗ a†
interaction_op2 = tensor(sigmam(), a.dag())     # σ⁻ ⊗ a†

# Full Hamiltonian list: [H0, [op1, coeff1(t)], ...]
H_tot = [
    H,
    [interaction_op1, interaction1],
    [interaction_op2, interaction2],
    [interaction_op1.dag(), interaction1.conjugate()],
    [interaction_op2.dag(), interaction2.conjugate()]
]

#%% Additional observables (conditioned squeezing moments)

qb_state1 = g.unit()
qb_state2 = e.unit()

# Add projectors conditioned on qubit being in |g⟩ or |e⟩ for quadrature observables
e_ops.extend([
    tensor(qb_state1.proj(), X ** 2),
    tensor(qb_state1.proj(), P ** 2),
    tensor(qb_state1.proj(), P * X),
    tensor(qb_state1.proj(), X * P),

    tensor(qb_state2.proj(), X ** 2),
    tensor(qb_state2.proj(), P ** 2),
    tensor(qb_state2.proj(), P * X),
    tensor(qb_state2.proj(), X * P),
])

#%% Solve the time-dependent Schrödinger equation
res = mesolve(
    H_tot,
    init,
    times,
    e_ops=e_ops,
    c_ops=False,  # Unitary evolution (no dissipation)
    options=Options(store_states=False, store_final_state=True),
    progress_bar=True
)

#%% Extract conditional squeezing values for |g⟩ and |e⟩

g_res = res.expect[-14]         # ⟨|g⟩⟨g|⟩
XX_g_res = res.expect[-8]       # ⟨X²⟩ conditioned on |g⟩
YY_g_res = res.expect[-7]       # ⟨P²⟩ conditioned on |g⟩
e_res = res.expect[-13]         # ⟨|e⟩⟨e|⟩
XX_e_res = res.expect[-4]       # ⟨X²⟩ conditioned on |e⟩
YY_e_res = res.expect[-3]       # ⟨P²⟩ conditioned on |e⟩

# Convert to squeezing in dB:
# Reference variance is 0.5 → normalize by ×2
Sq_g_dB_XX = -10 * np.log10(XX_g_res / g_res * 2)
Sq_g_dB_YY = -10 * np.log10(YY_g_res / g_res * 2)
Sq_e_dB_XX = -10 * np.log10(XX_e_res / e_res * 2)
Sq_e_dB_YY = -10 * np.log10(YY_e_res / e_res * 2)

#%% Plot squeezing vs time

plt.figure()
plt.plot(times, Sq_g_dB_XX, 'b-', label='XXg RWA')  # Ground state, X quadrature
plt.plot(times, Sq_g_dB_YY, 'r-', label='YYg RWA')  # Ground state, P quadrature
plt.plot(times, Sq_e_dB_XX, 'b--', label='XXe RWA') # Excited state, X quadrature
plt.plot(times, Sq_e_dB_YY, 'r--', label='YYe RWA') # Excited state, P quadrature
plt.ylabel('Squeezed amplitude [dB]')
plt.xlabel('Time [ns]')
plt.legend()
