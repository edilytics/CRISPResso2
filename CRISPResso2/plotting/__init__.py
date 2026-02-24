"""CRISPResso2 plotting data contract.

This subpackage provides PlotContext, a read-only dataclass that
CRISPRessoPro and plugins use to access analysis results for
generating custom plots.
"""

from CRISPResso2.plotting.plot_context import PlotContext

__all__ = ['PlotContext']
