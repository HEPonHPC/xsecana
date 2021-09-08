import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False)
g.add_node(
    0,
    label=uml_formatter('CrossSection<HistType>',
                        methods=['+ SetNTargets(ntargets : double) : void',
                                 '+ Eval(data : const HistType) : HistType'],

                        members=['- fEfficiency : IEfficiency<HistType>*',
                                 '- fSignalEstimator : ISignalEstimator<HistType>*',
                                 '- fFlux : IFlux*',
                                 '- fUnfold : IUnfold*']),
    shape='record',
    fontsize=20,
)

g.layout(prog='dot')
g.draw('cross_section_uml.svg', prog='dot')
