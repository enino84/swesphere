"""swesphere: shallow-water equations on the sphere for data-assimilation testbeds.

    from swesphere import presets, climatology
    model, x0 = presets.two_jets()
    S = climatology.build_climatology(model, x0, spinup_days=40, n_snapshots=200, every_days=2)
"""
from .model import SWEModel
from .grid import SphereGrid
from .forcing import ZonalRelaxation
from . import presets, climatology, initial_conditions, integrators, dynamics

__version__ = "0.1.0"
__all__ = ["SWEModel", "SphereGrid", "ZonalRelaxation", "presets", "climatology", "initial_conditions", "integrators", "dynamics", "__version__"]
