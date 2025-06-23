# Qubit‑Squeezing

Companion code for the paper\
**"Single‑ and two‑mode squeezing by modulated coupling to a Rabi‑driven qubit"**\
E. Blumenthal1, N. Gutman2, I. Kaminer2, S. Hacohen‑Gourgy1 (2025)

## Abstract

We present a novel method for generating conditional single‑ and two‑mode squeezing using a Rabi‑driven qubit dispersively coupled to one or two harmonic oscillators. This enables universal control over bosonic modes, with predicted intra‑cavity squeezing of **13 dB** (single‑mode) and **12 dB** (two‑mode).

---

## Repository layout

```
├── singleModeApproximated.py        # Effective squeezing simulation (Fig. 2a/3)
├── singleModeDisplacedFrame.py      # Full driven model in the displaced frame simulation (Fig. 2 and Fig. 3)
├── singleModeModulatedCoupling.py   # Effective modulated coupling simulation (Fig. 2a/3)
├── singleModeSqueezing.nb           # Mathematica notebook deriving single-mode Hamiltonian
├── twoModeSqueezing.nb              # Mathematica notebook deriving two-mode Hamiltonian
└── README.md
```



## Notebooks

The `.nb` files are Mathematica notebooks that generate the high‑resolution Wigner‑function plots shown in the paper. If you do not have Mathematica you can still open them with *Wolfram Player*.

---

## How to cite

If you use this code, please cite **the paper**:

```bibtex
@article{Blumenthal2025Squeezing,
  title   = {Single- and two-mode squeezing by modulated coupling to a Rabi driven qubit},
  author  = {E. Blumenthal and N. Gutman and I. Kaminer and S. Hacohen-Gourgy},
  journal = {arXiv preprint arXiv:XXXXX},
  year    = {2025}
}


---

## License

This project is licensed under the **MIT License** – see [`LICENSE`](LICENSE) for details.

## Contact

For questions, open an [issue](https://github.com/<your‑username>/qubitSqueezing/issues) or email [**eliya.blumenthal@technion.ac.il**](mailto\:eliya.blumenthal@technion.ac.il).

---

### Acknowledgements

This work was done in the Technion – Israel Institute of Technology.\
The code benefitted from QuTiP community resources.

