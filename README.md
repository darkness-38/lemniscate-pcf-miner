# lemniscate-pcf-miner

> **High-Performance Polynomial Continued Fraction (PCF) Search Engine for Gauss's Lemniscate and Elliptic Constants**

[![Rust](https://img.shields.io/badge/rust-1.80%2B-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/tests-passing-brightgreen.svg)]()

`lemniscate-pcf-miner` is a high-throughput, multi-threaded mathematical discovery engine written in Rust. Applying the methodology of the **Ramanujan Machine**, it systematically explores, discovers, and proves polynomial continued fraction (PCF) representations for **Gauss's Lemniscate Constant ($\varpi \approx 2.62205755...$)**, the **Gauss AGM Constant ($\pi/\varpi \approx 1.19814023...$)**, and their associated Möbius transformations ($\frac{m_1 S + m_2}{m_3 S + m_4}$).

---

## 🌟 Key Mathematical Discoveries

1. **Tail Evaluation Architecture ($K = \frac{b_1}{a_1 + \dots}$):**  
   Evaluates generalized continued fractions directly from $n=1$, eliminating arbitrary integer shifts induced by $a_0$.
2. **Twin Duality Symmetry ($2/\varpi \longleftrightarrow \pi/\varpi$):**  
   In the Euler-Pfaff family $a_n = 6n, \; b_n(k) = -(2n - 1)(4n - (2k + 1))$:
   - **Odd $k$ ($k=1, 3$):** Converges strictly to the Lemniscate family: $K(1) = \frac{2}{\varpi} - 1$, $K(3) = \frac{4}{\varpi} - 1$.
   - **Even $k$ ($k=2, 4$):** Converges strictly to the Gauss AGM family: $K(2) = \frac{\pi}{\varpi} - 1$, $K(4) = \frac{3\pi}{2\varpi} - 1$.
3. **Equivalence Reduction from 4th Degree to Linear Roots:**  
   Of **103 high-precision identities** discovered, exactly **52 are independent canonical root identities**. All degree-4 polynomial configurations ($b_n \sim n^4$) were proven to be polynomial equivalence transformations ($c_n = n+1$) of linear roots ($a_n' = (n+1)a_n, \; b_n' = n(n+1)b_n$).
4. **Record Convergence Rate:**  
   Achieved up to **82 decimal digits of precision in 150 steps** ($\Delta < 10^{-82}$) in the Gauss-Pfaff family ($a_n = 4n - 3, \; b_n = -n(2n - 5) \implies \frac{3\varpi + 6}{\varpi + 4}$).

---

## 🚀 Performance & Architecture

Optimized for **Intel Core i3-N305** (8 Gracemont E-Cores, AVX2, 6MB L3 cache):
- **Zero-Allocation Inner Loop:** Register-only stack primitives during Euler-Wallis recurrence.
- **Cache-Friendly Chunking:** Lock-free Rayon work-stealing with batching.
- **Multi-Tier Floating Scaling:** Multipliers preventing underflow/overflow beyond $10^{\pm 80}$.
- **Throughput:** Scans **$31.07 \times 10^6$ combinations/second**, evaluating $1.81 \times 10^7$ candidates in **582 milliseconds**.
- **Arbitrary-Precision Verification:** 150-step exact integer arithmetic (`IBig`) + 160-digit decimal division (`DBig`) via `dashu`.

---

## 📦 Project Structure

```text
lemniscate-pcf-miner/
├── Cargo.toml
├── README.md                              # This overview
├── PAPER_DRAFT.md                         # Academic pre-print draft & proofs
├── discovered_lemniscate_identities.txt   # Complete log of 103 verified identities
└── src/
    ├── lib.rs                             # Library entry point
    ├── main.rs                            # Parallel search engine CLI
    ├── pcf.rs                             # PCF templates, O(1) Darboux filters, Wallis recurrence
    ├── verify.rs                          # Target Möbius table & 160-digit dashu verification
    ├── analyze.rs                         # Equivalence reduction & shift analysis module
    └── bin/
        └── analyze.rs                     # Analysis tool generating PAPER_DRAFT.md
```

---

## 🛠️ Quickstart

### Prerequisites
- [Rust toolchain](https://www.rust-lang.org/) (version 1.80 or newer recommended)

### Build and Run Tests
```bash
cargo test
```

### Run Full Search
Runs the 18.1-million combination search across all Mod-4 hypergeometric families with default 7 threads:
```bash
cargo run --release
```

### Run Identity Analysis & Paper Generator
Classifies the 103 discovered identities, performs equivalence reduction, runs parametric shift analysis, and updates `PAPER_DRAFT.md`:
```bash
cargo run --release --bin analyze
```

---

## 📊 Discovered Formulas Highlights

| Constant | $a_n$ | $b_n$ | Closed Form Target | Error ($\Delta_{150}$) |
|:---|:---|:---|:---|:---:|
| **$\varpi$** | $4n - 3$ | $-n(2n - 5)$ | $\dfrac{3\varpi + 6}{\varpi + 4}$ | **$5.2 \times 10^{-82}$** |
| **$\varpi$** | $4n + 1$ | $-n(2n - 1)$ | $\dfrac{2}{\varpi} - 1$ | **$1.1 \times 10^{-81}$** |
| **$\varpi$** | $4n + 5$ | $-(n+2)(2n - 1)$ | $\varpi - 3$ | **$3.8 \times 10^{-81}$** |
| **$\varpi$** | $6n$ | $-(2n - 1)(4n - 3)$ | $\dfrac{2}{\varpi} - 1$ | **$4.8 \times 10^{-48}$** |
| **$\pi/\varpi$** | $4n - 1$ | $-n(2n - 3)$ | $\dfrac{\pi/\varpi}{\pi/\varpi + 2}$ | **$3.4 \times 10^{-52}$** |
| **$\pi/\varpi$** | $6n$ | $-(2n - 1)(4n - 5)$ | $\dfrac{\pi}{\varpi} - 1$ | **$7.6 \times 10^{-51}$** |
| **$\pi/\varpi$** | $4n + 3$ | $-n(2n + 1)$ | $\dfrac{-3(\pi/\varpi) + 4}{\pi/\varpi - 2}$ | **$5.4 \times 10^{-51}$** |

For the full catalog of identities and mathematical proofs, refer to [PAPER_DRAFT.md](PAPER_DRAFT.md).

---

## 📜 Citation & License

This project is licensed under the [MIT License](LICENSE).

```bibtex
@misc{lemniscate_pcf_miner_2026,
  author = {Betül Daşdemir},
  title = {Polynomial Continued Fractions for Gauss's Lemniscate and Elliptic Constants: Twin Symmetries and Parametric Shift Families},
  year = {2026},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/darkness-38/lemniscate-pcf-miner}}
}
```
