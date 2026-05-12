# Rough Volatility Models: Pricing & Calibration

### Overview

This project implements three **rough volatility models** for option pricing: Hurst exponent estimation, Lifted Heston (Carr-Madan / FFT) and Lifted Bergomi (Monte Carlo) motivated by the empirical finding that realized volatility behaves like a fractional Brownian motion with $H \approx 0.1$ (Gatheral et al., 2018).

## 1 - Hurst Exponent Estimation
**[Full Doc](docs/hurst_estimation.md) | [Notebook](notebooks/hurst_estimation.ipynb)**

The Hurst exponent $H$ controls process roughness through increment scaling: $\mathbb{E}[|X_{t+\Delta} - X_t|^2] \sim \Delta^{2H}$

We simulate three processes over $[0,1]$ ($m = 10^5$ steps) and estimate $H$ via log-moment regression across moments $q \in \{0.5, 1, 1.5, 2\}$ and lags $\Delta = 1,\ldots,10$:

$$\log(m(q,\Delta)) \approx \zeta_q \log(\Delta) + C_q \implies \hat{H} = \frac{\zeta_q}{q}$$

<div align="center">

| Process | Simulation method | Expected $H$ |
|:-------:|:-----------------:|:------------:|
| Brownian Motion | Euler | 0.5 |
| Fractional Brownian Motion | Cholesky | 0.1 |
| Lifted Heston variance $V^n$ | Implicit-Explicit Euler | ? |
</div>

**Results:** The Lifted Heston variance path is visually indistinguishable from an fBM, yet its estimated Hurst exponent sits at $\hat{H} \approx 0.47$ at the variance level (a semimartingale effect) before dropping toward $\hat{H} \approx 0.2$ at finer timescales on volatility. In all three cases, the linearity of $\zeta_q$ in $q$ confirms **monofractality**.

## 2 - Implied Volatility in the Lifted Heston Model
**[Full Doc](docs/lifted_heston_iv.md) | [Notebook](notebooks/lifted_heston_iv.ipynb)**

The affine structure yields a closed-form characteristic function through Riccati ODEs, enabling Fourier pricing through the **Carr-Madan formula**:

$$C_0 = \frac{e^{-r T - \alpha_2 \log K}}{\pi} \int_0^{\infty} \text{Re} \left( \frac{\Phi_T(u - (\alpha_2+1)i)}{(\alpha_2 + iu)(\alpha_2 + 1 + iu)} e^{-i\log(K)u} \right) du$$

We study smile convergence as $n \to \infty$ across two regimes:

<div align="center">

| Maturity | Convergence |
|:--------:|:-----------:|
| $T = 1$ |  Fast: $n \approx 10$ sufficient |
| $T = 1/26$ | Slow: $n > 20$ needed |
</div>

**Observations:** At long maturity, the smile converges quickly with $n \approx 10$ factors already sufficient. However, at short maturity, convergence is notably slower (due to the rough volatility).

## 3 - Implied Volatility in the Lifted Bergomi Model
**[Full Doc](docs/lifted_bergomi_iv.md) | [Notebook](notebooks/lifted_bergomi_iv.ipynb)**

The Lifted Bergomi replaces the affine structure with an exponential variance inherited from Rough Bergomi:

$$V_t^n = \xi_0(t) \exp\left(\eta\sqrt{2H}\Gamma\left(H+\tfrac{1}{2}\right)\int_0^t K^n(t-s)dW_s - H\eta^2\Gamma^2\left(H+\tfrac{1}{2}\right)\int_0^t (K^n(s))^2ds\right)$$

Without a tractable characteristic function, pricing is done via **Monte Carlo**: $C_0 = e^{-rT}\mathbb{E}[(S_T^n - K)^+]$

**Observations:** The Lifted Bergomi produces a steeper and more realistic skew than the Lifted Heston, though this comes at the cost of Monte Carlo noise in the wings and a significantly higher computational burden at short maturities.

## Summary

<div align="center">

| Categories | Lifted Heston | Lifted Bergomi |
|:---:|:---:|:---:|
| **Pricing** | Carr-Madan (FFT) | Monte Carlo |
| **Characteristic function** | Closed-form | Not available |
| **Convergence in $n$** | Fast (long mat) / Slow (short mat) | Moderate |
| **Skew** | Moderate | Steep |
| **Short maturity** | Tractable | Expensive |
</div>

## Perspectives: VIX Options Pricing

A natural extension would be **VIX options pricing** under each model. In the Lifted Heston, the affine structure extends to forward variance, making a Fourier approach feasible. In the Lifted Bergomi, we could reuse the Monte Carlo engine from Part 3 by simulating forward variance paths and computing $\text{VIX}^2$ as a realized average but at the cost of higher variance and computational load in the wings.

## References

- Abi Jaber, E. (2019). *Lifting the Heston model*. Quantitative Finance.
- El Euch, O., Rosenbaum, M. (2019). *The characteristic function of rough Heston models*. Mathematical Finance.
- Gatheral, J., Jaisson, T., Rosenbaum, M. (2018). *Volatility is rough*. Quantitative Finance.
- Bayer, C., Friz, P., Gatheral, J. (2016). *Pricing under rough volatility*. Quantitative Finance.
- Carr, P., Madan, D. (1999). *Option valuation using the fast Fourier transform*. Journal of Computational Finance.
