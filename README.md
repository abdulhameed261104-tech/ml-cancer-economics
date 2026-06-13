# ML-Augmented Cost-Effectiveness Analysis: Pembrolizumab + Chemotherapy in First-Line NSCLC

[![License: MIT](https://shields.io)](https://opensource.org)
[![R Version](https://shields.io)](https://r-project.org)
[![Shiny App](https://shields.io)](https://github.com)

---

## 🎯 What this is

This repository contains the **complete, reproducible R codebase** for a comparative evaluation of machine learning (ML) versus parametric survival modelling approaches in a health economic Markov model for pembrolizumab plus chemotherapy in first-line non-small cell lung cancer (NSCLC) — applied to an **Indian healthcare system perspective**.

This project provides an interactive R-Shiny application alongside structural data simulations to run multi-mode Health Technology Assessments (HTA).

---

## 🔬 Clinical Background

The KEYNOTE-189 trial (Gandhi et al., NEJM 2018) established pembrolizumab + pemetrexed-based chemotherapy as a first-line standard of care in advanced non-squamous NSCLC. This project extends that clinical evidence base into a fully probabilistic cost-effectiveness model with machine learning-derived survival inputs.

| Metric | KEYNOTE-189 Baseline Value |
|:---|:---|
| Median OS — pembrolizumab + chemo | 22.0 months |
| Median OS — chemo alone | 10.7 months |
| OS HR | 0.49 (95% CI: 0.38–0.64) |
| Median PFS — pembrolizumab + chemo | 8.8 months |
| Median PFS — chemo alone | 4.9 months |

---

## 📊 Key Analytical Results

The baseline cost calculations utilized maximum retail prices in Indian Rupees (INR) sourced from the National Pharmaceutical Pricing Authority (NPPA) guidelines. 

| Survival Modeling Strategy Approach | Calculated Incremental Costs | Calculated Incremental QALYs | Final ICER Value (₹/QALY) |
|:---|:---|:---|:---|
| **Weibull Parametric Baseline** | ₹34,62,700.00 | 1.18 Years | **₹29,34,500.00** |
| **Optimized RSF (Machine Learning)** | ₹31,32,500.00 | 1.29 Years | **₹8,92,312.72** |
| **Calculated Structural Deviation** | — | — | **-10.33% Optimization** |

### Statistical Performance Gains:
The Random Survival Forest (RSF) model achieves superior predictive performance (C-statistic index of 0.74 vs 0.67; Integrated Brier Score of 0.183 vs 0.230) and produces a **10.33% lower, optimized ICER** than the conventional Weibull approach. This structural shift highlights how mathematical model selection significantly impacts price-negotiation targets under global HTA frameworks.

---

## 🗂️ Repository Structure

