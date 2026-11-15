import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)

def rep(old, new):
    global s
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new, 1)

# abstract: the interface indicators are an empirical screen, not a bound
rep(r"""Together, the three interface indicators form an $\mathcal{O}(n_c (r_u + r_\lambda))$ online admissibility circuit breaker ($\eta_{\mathrm{adm}}$) independent of the bulk mesh dimension $N$.""",
    r"""The three interface indicators, which omit the equilibrium residual and therefore are not an error bound, are combined into an empirical $\mathcal{O}(n_c (r_u + r_\lambda))$ online admissibility screen ($\eta_{\mathrm{adm}}$, and a dimensionally consistent variant $\eta_{\mathrm{adm,v2}}$) independent of the bulk mesh dimension $N$, to be used with a validation-calibrated threshold.""")

# contributions: separate the three claims
rep(r"""acting as an online circuit breaker.""",
    r"""acting as an empirical online circuit breaker. We keep three claims apart: the extension of the full-order residual decomposition to arbitrary lifted pairs (the setting), the energy-norm bound itself (a rigorous result, with the exact multiplier norm or a valid upper bound on it), and the interface screen (an empirical detector with a calibrated threshold, not a certificate).""")

# alpha in the cost model
rep(r"""(or fill-in factor \(\alpha \in [1.5, 2.0]\) under sparse Cholesky / algebraic multigrid preconditioning \cite{franceschini2022block})""",
    r"""(or, in three dimensions, a cost growing as \(N^\alpha\), where the exponent \(\alpha\) lumps the fill-in of the sparse Cholesky factor and the ordering quality; \(\alpha \approx 2\) is the classical value for nested-dissection orderings in 3D, and smaller values apply to well-ordered 2D problems \cite{franceschini2022block})""")

io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
