"""Entry point. Pins BLAS/OpenMP to one thread *before* numpy is imported (E4 protocol: single-thread timings)."""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

from .run import main  # noqa: E402  (must come after the environment is set)

main()
