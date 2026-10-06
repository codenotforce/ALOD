"""Compatibility fields shared by new runs and checkpoint validation."""

# These are migration defaults, not today's production defaults. In particular,
# an old checkpoint without an enrichment policy must retain archived adjoint
# tests rather than silently switching to the current kernel-lift formulation.
COMPATIBILITY_DEFAULTS = {
    "checkpoint_interval_cycles": 0,
    "checkpoint_retention": "all",
    "keep_final": False,
    "exact_target": -1.0,
    "exact_scope": "nominal",
    "audit_mode": "full",
    "ell_absolute_threshold": -1.0,
    "wavenumber": 16,
    "ell_ratio_mode": "raw",
    "ell_threshold": 0.0,
    "enrichment_tests": "adjoint",
}


def with_compatibility_defaults(config):
    return {**COMPATIBILITY_DEFAULTS, **config}
