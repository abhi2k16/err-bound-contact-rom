import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)

anchor = r"\hypertarget{offline-and-online-procedure}{%"
assert s.count(anchor) == 1

block = r"""\hypertarget{cost-model-assumptions}{%
\subsubsection{Assumptions and Domain of Validity}\label{cost-model-assumptions}}

Equations~\eqref{eq:cfom}--\eqref{eq:crom} are an operation-count model, and the break-even sizes they give hold only under the assumptions below. We separate what the model predicts (floating-point operations) from what is measured (single-thread wall-clock time of our Python implementation).

\begin{enumerate}
\tightlist
\item \textbf{Same contact solver on both sides.} The full-order and the reduced problem are both solved by the same dense active-set LCP method without warm starts, in one thread, for frictionless unilateral contact. A different complementarity solver (semismooth Newton, interior point, warm-started active sets) changes the iteration counts $k_{\mathrm{FOM}}$ and $k_{\mathrm{ROM}}$ and the per-iteration costs.
\item \textbf{Direct factorization in the full-order model.} The first two terms of Eq.~\eqref{eq:cfom} assume that $K$ is factored directly and that $B K^{-1} B^T$ is formed explicitly with $n_c$ triangular solves. The implementation counts the banded cost $c_{\mathrm{fact}} N b_w^2$ for every benchmark, with the bandwidth $b_w$ of the assembled matrix and unfitted constants $c_{\mathrm{fact}} = c_{\mathrm{solve}} = 1$ and $c_{\mathrm{proj}} = 2$, whereas the measured full-order solver is a sparse LU with a fill-reducing ordering. The $N^\alpha$ form quoted above for three dimensions was not used. The break-even sizes are therefore sensitive to these constants, and we did not study that sensitivity.
\item \textbf{Iterative full-order solvers are not modelled.} A preconditioned Krylov or multigrid solver costs about $n_{\mathrm{it}}\,(2\,\mathrm{nnz}(K) + C_{\mathrm{prec}})$ per right-hand side, which grows roughly linearly in $N$, against $N^{3/2}$ (two dimensions) or $N^{2}$ (three dimensions) for a direct factorization. Forming $B K^{-1} B^T$ then needs $n_c$ solves unless a primal--dual or Schur-complement-free method is used. Either way the first two terms of Eq.~\eqref{eq:cfom} would be replaced and the break-even size would move to larger $N$. We did not test an iterative full-order solver, so our break-even sizes apply to direct solvers.
\item \textbf{Memory is not modelled.} The model counts operations, not memory traffic or storage. A sparse factor has fill-in of order $N \log N$ (two dimensions) to $N^{4/3}$ (three dimensions) with nested dissection, whereas the reduced model stores dense $r_u \times r_u$, $r_\lambda \times r_\lambda$ and $n_c \times r_\lambda$ blocks, and $N \times r_u$ if the full field is lifted. The memory advantage of the reduced model at large $N$ is not quantified here.
\item \textbf{Projection in the varying-stiffness regime.} For $K = K(\mu)$ the reduced operators are re-projected from the full-order matrices for every query, at the cost $c_{\mathrm{proj}}\bigl(\mathrm{nnz}(K) r_u + \mathrm{nnz}(B)(r_u + r_\lambda)\bigr)$. With an affine decomposition or hyper-reduction this term is smaller and the break-even size would move down.
\item \textbf{Operation counts are not run times.} The model omits call overhead, memory traffic and the lower efficiency of small dense and of sparse kernels. In our implementation a reduced solve takes about $0.2$--$0.3$\,ms even when it needs only a few thousand operations, so the measured speedup exceeds the predicted one at small $N$ and is far below it at large $N$ (Table~\ref{tab:timing}). Only the location of the break-even size, not the magnitude of the speedup, should be read from the model, and both are specific to this code and machine.
\end{enumerate}

"""
s = s.replace(anchor, block + anchor, 1)
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
