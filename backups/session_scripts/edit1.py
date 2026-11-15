import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
raw = io.open(p, encoding="utf-8", newline="").read()
CR = chr(13)
LF = chr(10)
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)

def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new, 1)

# 1. proof: sign of the complementarity term
rep(r"""The third term is bounded immediately by \(\left| \lambda_+^T (B \tilde{u} - c) \right| = \eta_{\mathrm{comp}}\).""",
r"""The third term is bounded immediately by \(\left| \lambda_+^T (B \tilde{u} - c) \right| = \eta_{\mathrm{comp}}\). This step uses \(-x \le |x|\) and is therefore valid whichever sign the inner product \(\lambda_+^T (B \tilde{u} - c)\) has. The absolute value in the definition of \(\eta_{\mathrm{comp}}\) is essential: it makes the indicator non-negative, so that \(\sqrt{\eta_{\mathrm{comp}}}\) is real, and it covers the case in which \(\lambda_+\) is not complementary to the gap of a lifted displacement. That case does occur: on the audited queries the inner product was negative for 90 of 200 queries on B1 and for 74 of 200 on B2 (five dual bases, two dual dimensions, 20 queries each), and \eqref{eq:bound1} held on all of them.""")

# 2. Remark 1: conservativeness of the stability bound
rep(r"""(iii) The computable proxy \(\|\tilde{\lambda}_+\|_2\) carries no guarantee, because it need not bound the exact multiplier norm; we do not use it.""",
r"""(iii) The computable proxy \(\|\tilde{\lambda}_+\|_2\) carries no guarantee, because it need not bound the exact multiplier norm; we do not use it. The stability bound of the previous paragraph is rigorous but conservative. On 20 test queries per benchmark, \(\beta^{-1}(\|f\|_{K^{-1}} + \|u\|_K)\) reached 20.1 (B1) and 26.2 (B2), with \(\beta_{\min} = 0.70\) and \(0.36\), against largest true multiplier norms of 2.09 and 1.63 and offline constants \(C_\lambda\) of 2.37 and 1.64. It is thus a defensible but roughly tenfold looser alternative to \(C_\lambda\), and it would inflate the gap term of the bound by the same factor. \(C_\lambda\) is the sharper choice inside the training range; the inf-sup bound is the one to use when a guarantee outside the range is required and \(\beta\) can be computed offline.""")

# 3. new remark on the three levels of claim, after Remark 2
rep(r"""which is weaker than \eqref{eq:bound1} by a factor two in the equilibrium term.
\end{remark}
""",
r"""which is weaker than \eqref{eq:bound1} by a factor two in the equilibrium term.
\end{remark}

\begin{remark}[What is certified and what is empirical]

Three different statements are made in this paper and they carry different guarantees. (a) \emph{Rigorous bound.} Theorem~\ref{thm:bound}, with the exact multiplier norm \(\|\lambda\|_2\) or with any valid upper bound on it (for example the inf-sup bound of Remark~1), is a mathematical upper bound on the energy-norm error. It contains the equilibrium residual \(\eta_{\mathrm{eq}}\), which needs a solve with \(K\) unless the problem is affine, so it is a validation quantity. (b) \emph{Online-usable bound.} Replacing \(\|\lambda\|_2\) by the offline constant \(C_\lambda\) gives a bound that is rigorous only while \(\|\lambda(\mu)\|_2 \le C_\lambda\); outside the training range it is a heuristic. (c) \emph{Interface screening indicator.} The indicators \(\eta_{\mathrm{adm}}\) and \(\eta_{\mathrm{adm,v2}}\) of Section~\ref{computational-complexity-and-online-admissibility-indicator} omit \(\eta_{\mathrm{eq}}\) and are compared with a threshold \(\tau_{\mathrm{adm}}\) calibrated on validation queries. They are not error bounds: a small value does not imply a small error, because the bulk equilibrium residual can be large, and the threshold has no guarantee outside the calibration data. An accepted prediction is therefore not certified. It can be called trustworthy only to the extent of the empirical recall reported in Section~\ref{sec:experiments}, and the screening indicator is meant to trigger a fallback, not to replace the bound of item (a).
\end{remark}
""")

# 4. Section 3.3 indicator paragraph
a = s.index(r"To provide a fast, dimensionally consistent online ``circuit breaker''")
b = s.index(r"% ======================================================================", a)
new_par = r"""To provide a fast online ``circuit breaker'' that operates solely on the contact interface \cite{dabaghi2020posteriori, fauque2018hybrid}, we combine the three \(\mathcal{O}(n_c (r_u + r_\lambda))\) terms into the \textbf{Online Admissibility Indicator}:
\begin{equation}\label{eq:eta-adm}
\eta_{\mathrm{adm}}(\tilde{u}, \tilde{\lambda}) := \eta_{\mathrm{sign}} + \eta_{\mathrm{gap}} + \sqrt{\eta_{\mathrm{comp}}}.
\end{equation}
Taking the square root of \(\eta_{\mathrm{comp}}\) makes the three terms scale with the first power of the error \(\mathcal{O}(\|e\|_K)\) in the sense of the order of smallness, matching the scaling of \(\|u - \tilde{u}\|_K\) in Theorem~\ref{thm:bound}. It does \emph{not} make the sum dimensionally consistent: \(\eta_{\mathrm{sign}}\) has the units of a multiplier (force or traction), \(\eta_{\mathrm{gap}}\) those of a displacement, and \(\sqrt{\eta_{\mathrm{comp}}}\) those of the square root of a force times a displacement, which is the unit of the energy norm. Under a change of the force unit by a factor \(s\) (\(K, f, \lambda\) scale by \(s\); displacements and gaps are unchanged) the three terms scale by \(s\), \(1\) and \(\sqrt{s}\), so the relative weight of the terms, and hence the meaning of one value of \(\tau_{\mathrm{adm}}\), changes with the units and with the load and stiffness scale. Equation~\eqref{eq:eta-adm} is therefore a heuristic whose threshold must be recalibrated whenever the scales change.

The dimensionally consistent alternative is the interface part of bound~\eqref{eq:bound2},
\begin{equation}\label{eq:eta-adm-v2}
\eta_{\mathrm{adm,v2}} := \Bigl(2\,\|B\|_{K \to 2}^2\,\eta_{\mathrm{sign}}^2 + 2\,\eta_{\mathrm{comp}} + 2\,C_\lambda\,\eta_{\mathrm{gap}}\Bigr)^{1/2},
\end{equation}
in which the exact multiplier norm is replaced by the offline constant \(C_\lambda\) and the trace constant \(\|B\|_{K \to 2}\) is precomputed offline (for parameter-dependent \(B\) it is recomputed or bounded offline over the parameter domain). Every term under the root has the units of force times displacement, so \(\eta_{\mathrm{adm,v2}}\) has the units of the energy norm; under a force-unit change by \(s\) it scales by \(\sqrt{s}\), exactly like \(\|u - \tilde{u}\|_K\), so a relative threshold \(\tau_{\mathrm{adm}}/\|u\|_K\) is unit-independent. Its online cost is the same \(\mathcal{O}(n_c (r_u + r_\lambda))\), and it is reported alongside Eq.~\eqref{eq:eta-adm} in Section~\ref{sec:experiments}. Neither indicator contains \(\eta_{\mathrm{eq}}\), so neither is an error bound (Remark~3). Whenever an indicator exceeds a user-specified tolerance \(\tau_{\mathrm{adm}}\), the reduced prediction is flagged as physically inadmissible---even if the reduced KKT residual is at machine precision---triggering either adaptive basis enrichment or fallback to the full-order solver.

"""
s = s[:a] + new_par + s[b:]
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
