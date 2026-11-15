import io
base = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability"
CR, LF = chr(13), chr(10)

# ------------------------------------------------------------------------------------------------ bibliography
BIB = r"""
@article{veeser2001efficient,
  author = {Veeser, A.},
  title = {Efficient and reliable a posteriori error estimators for elliptic obstacle problems},
  journal = {SIAM Journal on Numerical Analysis},
  volume = {39},
  number = {1},
  pages = {146--167},
  year = {2001},
  doi = {10.1137/S0036142900370812}
}

@article{braess2007convergence,
  author = {Braess, D. and Carstensen, C. and Hoppe, R. H. W.},
  title = {Convergence analysis of a conforming adaptive finite element method for an obstacle problem},
  journal = {Numerische Mathematik},
  volume = {107},
  number = {3},
  pages = {455--471},
  year = {2007},
  doi = {10.1007/s00211-007-0098-6}
}

@article{wohlmuth2007posteriori,
  author = {Wohlmuth, B. I.},
  title = {An a posteriori error estimator for two-body contact problems on non-matching meshes},
  journal = {Journal of Scientific Computing},
  volume = {33},
  number = {1},
  pages = {25--45},
  year = {2007},
  doi = {10.1007/s10915-007-9139-7}
}

@article{weiss2009posteriori,
  author = {Weiss, A. and Wohlmuth, B. I.},
  title = {A posteriori error estimator and error control for contact problems},
  journal = {Mathematics of Computation},
  volume = {78},
  number = {267},
  pages = {1237--1267},
  year = {2009},
  doi = {10.1090/S0025-5718-09-02235-2}
}

@article{bader2016empirical,
  author = {Bader, E. and Zhang, Z. and Veroy, K.},
  title = {An empirical interpolation approach to reduced basis approximations for variational inequalities},
  journal = {Mathematical and Computer Modelling of Dynamical Systems},
  volume = {22},
  number = {4},
  pages = {345--361},
  year = {2016},
  doi = {10.1080/13873954.2016.1198388}
}

@article{zhang2016slack,
  author = {Zhang, Z. and Bader, E. and Veroy, K.},
  title = {A slack approach to reduced-basis approximation and error estimation for variational inequalities},
  journal = {Comptes Rendus Math{\'e}matique},
  volume = {354},
  number = {3},
  pages = {283--289},
  year = {2016},
  doi = {10.1016/j.crma.2015.10.024}
}

@article{gerner2012certified,
  author = {Gerner, A.-L. and Veroy, K.},
  title = {Certified reduced basis methods for parametrized saddle point problems},
  journal = {SIAM Journal on Scientific Computing},
  volume = {34},
  number = {5},
  pages = {A2812--A2836},
  year = {2012},
  doi = {10.1137/110854084}
}

@article{ballarin2015supremizer,
  author = {Ballarin, F. and Manzoni, A. and Quarteroni, A. and Rozza, G.},
  title = {Supremizer stabilization of {POD}--{G}alerkin approximation of parametrized steady incompressible {N}avier--{S}tokes equations},
  journal = {International Journal for Numerical Methods in Engineering},
  volume = {102},
  number = {5},
  pages = {1136--1161},
  year = {2015},
  doi = {10.1002/nme.4772}
}

@article{rozza2007stability,
  author = {Rozza, G. and Veroy, K.},
  title = {On the stability of the reduced basis method for {S}tokes equations in parametrized domains},
  journal = {Computer Methods in Applied Mechanics and Engineering},
  volume = {196},
  number = {7},
  pages = {1244--1260},
  year = {2007},
  doi = {10.1016/j.cma.2006.09.005}
}

@article{benner2015survey,
  author = {Benner, P. and Gugercin, S. and Willcox, K.},
  title = {A survey of projection-based model reduction methods for parametric dynamical systems},
  journal = {SIAM Review},
  volume = {57},
  number = {4},
  pages = {483--531},
  year = {2015},
  doi = {10.1137/130932715}
}

@article{chaturantabut2010nonlinear,
  author = {Chaturantabut, S. and Sorensen, D. C.},
  title = {Nonlinear model reduction via discrete empirical interpolation},
  journal = {SIAM Journal on Scientific Computing},
  volume = {32},
  number = {5},
  pages = {2737--2764},
  year = {2010},
  doi = {10.1137/090766498}
}

@article{hernandez2017dimensional,
  author = {Hern{\'a}ndez, J. A. and Caicedo, M. A. and Ferrer, A.},
  title = {Dimensional hyper-reduction of nonlinear finite element models via empirical cubature},
  journal = {Computer Methods in Applied Mechanics and Engineering},
  volume = {313},
  pages = {687--722},
  year = {2017},
  doi = {10.1016/j.cma.2016.10.022}
}

@article{farhat2014dimensional,
  author = {Farhat, C. and Avery, P. and Chapman, T. and Cortial, J.},
  title = {Dimensional reduction of nonlinear finite element dynamic models with finite rotations and energy-based mesh sampling and weighting for computational efficiency},
  journal = {International Journal for Numerical Methods in Engineering},
  volume = {98},
  number = {9},
  pages = {625--662},
  year = {2014},
  doi = {10.1002/nme.4668}
}

@article{lee1999learning,
  author = {Lee, D. D. and Seung, H. S.},
  title = {Learning the parts of objects by non-negative matrix factorization},
  journal = {Nature},
  volume = {401},
  number = {6755},
  pages = {788--791},
  year = {1999},
  doi = {10.1038/44565}
}

@article{hild2010stabilized,
  author = {Hild, P. and Renard, Y.},
  title = {A stabilized {L}agrange multiplier method for the finite element approximation of contact problems in elastostatics},
  journal = {Numerische Mathematik},
  volume = {115},
  number = {1},
  pages = {101--129},
  year = {2010},
  doi = {10.1007/s00211-009-0273-z}
}

@article{wohlmuth2011variationally,
  author = {Wohlmuth, B.},
  title = {Variationally consistent discretization schemes and numerical algorithms for contact problems},
  journal = {Acta Numerica},
  volume = {20},
  pages = {569--734},
  year = {2011},
  doi = {10.1017/S0962492911000079}
}

@article{hintermuller2002primal,
  author = {Hinterm{\"u}ller, M. and Ito, K. and Kunisch, K.},
  title = {The primal-dual active set strategy as a semismooth {N}ewton method},
  journal = {SIAM Journal on Optimization},
  volume = {13},
  number = {3},
  pages = {865--888},
  year = {2002},
  doi = {10.1137/S1052623401383558}
}

@article{kornhuber2001adaptive,
  author = {Kornhuber, R. and Krause, R.},
  title = {Adaptive multigrid methods for {S}ignorini's problem in linear elasticity},
  journal = {Computing and Visualization in Science},
  volume = {4},
  number = {1},
  pages = {9--20},
  year = {2001},
  doi = {10.1007/s007910100052}
}

@article{wohlmuth2003monotone,
  author = {Wohlmuth, B. I. and Krause, R. H.},
  title = {Monotone multigrid methods on nonmatching grids for nonlinear multibody contact problems},
  journal = {SIAM Journal on Scientific Computing},
  volume = {25},
  number = {1},
  pages = {324--347},
  year = {2003},
  doi = {10.1137/S1064827502405318}
}

@article{george1973nested,
  author = {George, A.},
  title = {Nested dissection of a regular finite element mesh},
  journal = {SIAM Journal on Numerical Analysis},
  volume = {10},
  number = {2},
  pages = {345--363},
  year = {1973},
  doi = {10.1137/0710032}
}

@article{lipton1979generalized,
  author = {Lipton, R. J. and Rose, D. J. and Tarjan, R. E.},
  title = {Generalized nested dissection},
  journal = {SIAM Journal on Numerical Analysis},
  volume = {16},
  number = {2},
  pages = {346--358},
  year = {1979},
  doi = {10.1137/0716027}
}

@article{davis2016survey,
  author = {Davis, T. A. and Rajamanickam, S. and Sid-Lakhdar, W. M.},
  title = {A survey of direct methods for sparse linear systems},
  journal = {Acta Numerica},
  volume = {25},
  pages = {383--566},
  year = {2016},
  doi = {10.1017/S0962492916000076}
}

@article{demmel1999supernodal,
  author = {Demmel, J. W. and Eisenstat, S. C. and Gilbert, J. R. and Li, X. S. and Liu, J. W. H.},
  title = {A supernodal approach to sparse partial pivoting},
  journal = {SIAM Journal on Matrix Analysis and Applications},
  volume = {20},
  number = {3},
  pages = {720--755},
  year = {1999},
  doi = {10.1137/S0895479895291765}
}

@article{efron1979bootstrap,
  author = {Efron, B.},
  title = {Bootstrap methods: another look at the jackknife},
  journal = {The Annals of Statistics},
  volume = {7},
  number = {1},
  pages = {1--26},
  year = {1979},
  doi = {10.1214/aos/1176344552}
}

@article{fawcett2006introduction,
  author = {Fawcett, T.},
  title = {An introduction to {ROC} analysis},
  journal = {Pattern Recognition Letters},
  volume = {27},
  number = {8},
  pages = {861--874},
  year = {2006},
  doi = {10.1016/j.patrec.2005.10.010}
}

@article{kollepara2022limitations,
  author = {Kollepara, K. S. and Navarro-Jim{\'e}nez, J. M. and Le Guennec, Y. and Silva, L. and Aguado, J. V.},
  title = {On the limitations of low-rank approximations in contact mechanics problems},
  journal = {International Journal for Numerical Methods in Engineering},
  volume = {124},
  number = {1},
  pages = {217--234},
  year = {2022},
  doi = {10.1002/nme.7118}
}
"""
bp = base + r"\references.bib"
braw = io.open(bp, encoding="utf-8", newline="").read()
bcrlf = (CR + LF) in braw
b = braw.replace(CR + LF, LF).rstrip(LF) + LF
for key in [l.split("{")[1].rstrip(",") for l in BIB.splitlines() if l.startswith("@article")]:
    assert ("{" + key + ",") not in b, key
b += BIB
io.open(bp, "w", encoding="utf-8", newline="").write(b.replace(LF, CR + LF) if bcrlf else b)

# ------------------------------------------------------------------------------------------------ manuscript
p = base + r"\manuscript_master_v2.tex"
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)


def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:90])
    s = s.replace(old, new, 1)


rep(r"in design optimization, uncertainty quantification, and real-time control \cite{balajewicz2016projection, benaceur2020reduced}",
    r"in design optimization, uncertainty quantification, and real-time control \cite{benner2015survey, balajewicz2016projection, benaceur2020reduced}")
rep(r"alongside empirical cubature and parametric contact reduction strategies \cite{koller2021hyper, aresdeparga2023petrov, kollepara2024sparse, pawar2025efficient}",
    r"alongside empirical cubature and parametric contact reduction strategies \cite{farhat2014dimensional, hernandez2017dimensional, koller2021hyper, aresdeparga2023petrov, kollepara2024sparse, pawar2025efficient}")
rep(r"\cite{coorevits2000posteriori, hild2007residual, dabaghi2020posteriori, banz2024generic, fontana2024guaranteed, banz2026posteriori}. Within the reduced-basis framework",
    r"\cite{veeser2001efficient, coorevits2000posteriori, hild2007residual, braess2007convergence, wohlmuth2007posteriori, weiss2009posteriori, dabaghi2020posteriori, banz2024generic, fontana2024guaranteed, banz2026posteriori}. Within the reduced-basis framework")
rep(r"for reduced basis variational inequalities \cite{haasdonk2012reduced}",
    r"for reduced basis variational inequalities \cite{haasdonk2012reduced, bader2016empirical, zhang2016slack}")
rep(r"\cite{coorevits2000posteriori, wohlmuth2000mortar, franceschini2022block})",
    r"\cite{coorevits2000posteriori, wohlmuth2000mortar, wohlmuth2011variationally, hild2010stabilized, franceschini2022block})")
rep(r"on standard non-redundant contact interfaces \cite{franceschini2022block}",
    r"on standard non-redundant contact interfaces \cite{wohlmuth2011variationally, franceschini2022block}")
rep(r"r_{\text{rep}}\) \cite{leberre2022condition, rozza2013reduced}. Simultaneously",
    r"r_{\text{rep}}\) \cite{leberre2022condition, rozza2013reduced, rozza2007stability, gerner2012certified, ballarin2015supremizer}. Simultaneously")
rep(r"non-negative matrix factorization (NMF) \cite{balajewicz2016projection}, cone-projected greedy algorithms",
    r"non-negative matrix factorization (NMF) \cite{lee1999learning, balajewicz2016projection}, cone-projected greedy algorithms")
rep(r"state the CPG, gIS and RID constructions as implemented in our experiments",
    r"state the CPG, gIS and RID constructions as implemented in our experiments (the RID construction uses the discrete empirical interpolation method \cite{chaturantabut2010nonlinear})")
rep(r"smaller values apply to well-ordered 2D problems \cite{franceschini2022block}",
    r"smaller values apply to well-ordered 2D problems \cite{george1973nested, lipton1979generalized, davis2016survey, franceschini2022block}")
rep(r"are factored with SuperLU and the contact problems", r"are factored with SuperLU \cite{demmel1999supernodal} and the contact problems")
rep(r"(semismooth Newton, interior point, warm-started active sets)", r"(semismooth Newton \cite{hintermuller2002primal}, interior point, warm-started active sets)")
rep(r"A preconditioned Krylov or multigrid solver costs", r"A preconditioned Krylov or multigrid solver \cite{kornhuber2001adaptive, wohlmuth2003monotone} costs")
rep(r"Uncertainty is a $95\%$ percentile bootstrap interval resampling whole seeds",
    r"Uncertainty is a $95\%$ percentile bootstrap interval \cite{efron1979bootstrap} resampling whole seeds")
rep(r"Figure~\ref{fig:tradeoff} plots the false-negative rate against the share of queries accepted.",
    r"Figure~\ref{fig:tradeoff} plots the false-negative rate against the share of queries accepted, a trade-off curve of the receiver-operating-characteristic type \cite{fawcett2006introduction}.")
rep(r"provide no per-query error check; the authors of the greedy index selection",
    r"provide no per-query error check, and limitations of low-rank approximations in contact problems have been reported \cite{kollepara2022limitations}; the authors of the greedy index selection")
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
