# Backward-compatible re-exports: CRISPRessoPlot and _upsetplot moved to
# CRISPResso2.plots but can still be imported from the top-level package.
# We can remove this once we update CRISPRessoPro's importing
from CRISPResso2.plots import CRISPRessoPlot  # noqa: F401
from CRISPResso2.plots import _upsetplot  # noqa: F401
