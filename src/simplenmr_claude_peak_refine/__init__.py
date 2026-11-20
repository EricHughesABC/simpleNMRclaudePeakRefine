"""simpleNMRclaudePeakRefine - automatic refinement and harmonisation of NMR peak picking.

Refines the 1D and 2D peak lists exported from Mnova (simpleNMR JSON) so that every
spectrum agrees on one proton and one carbon shift list, ready for simpleNMR.
Developed by Eric Hughes (Durham University); code written by Claude (Anthropic)
under his direction.
"""
from .refine import VERSION as __version__, Params, refine, LabelError  # noqa: F401
from .check import consistency  # noqa: F401
