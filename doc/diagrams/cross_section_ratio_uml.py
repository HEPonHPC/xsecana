import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False)
g.add_node(
    0,
    label=uml_formatter('CrossSectionRatio<HistType>',
                        methods=['+ SetNTargetsNumerator(ntargets : double) : void',
                                 '+ SetNTargetsDenominator(ntargets : double) : void',
                                 '+ Eval(data : const HistType) : HistType'],

                        members=['- fNumerator : CrossSection<HistType>',
                                 '- fDenominator : CrossSection<HistType>',
                                 ]),
    shape='record',
    fontsize=20,
)

g.layout(prog='dot')
g.draw('cross_section_ratio_uml.pdf', prog='dot')
