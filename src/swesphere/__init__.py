"""swesphere: shallow-water equations on the sphere for data-assimilation testbeds.

    from swesphere import presets, climatology
    model, x0 = presets.two_jets()
    S = climatology.build_climatology(model, x0, spinup_days=40, n_snapshots=200, every_days=2)
"""
from .model import SWEModel
from .grid import SphereGrid
from .forcing import ZonalRelaxation
from . import presets, climatology, initial_conditions, integrators, dynamics, diagnostics

# Taken from the installed package metadata, so it cannot drift away from
# pyproject.toml the way a hand-written string does. The fallback covers
# running from a source tree that was never installed.
try:                                        # pragma: no cover
    from importlib.metadata import PackageNotFoundError, version as _version

    try:
        __version__ = _version("swesphere")
    except PackageNotFoundError:
        __version__ = "0.0.0+unknown"
except ImportError:                         # pragma: no cover
    __version__ = "0.0.0+unknown"
__all__ = ["SWEModel", "SphereGrid", "ZonalRelaxation", "presets", "climatology", "initial_conditions", "integrators", "dynamics", "diagnostics", "__version__"]
