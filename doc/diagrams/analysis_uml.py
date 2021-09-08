import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False)
g.add_node(
    0,
    label=uml_formatter('Analysis<HistType>',
                        methods=['+ FractionalUncertainty(syst_name : string) : HistType',
                                 '+ AbsoluteUncertainty(syst_name : string) : HistType',
                                 '+ TotalFractionalUncertainty() : std::pair<HistType, HistType>',
                                 '+ TotalAbsoluteUncertainty() : std::pair<HistType, HistType>'],
                        members=['- fNominalMeasurement : IMeasurement<HistType>*',
                                 '- fShiftedMeasurments : Systematic<IMeasurement*>',
                                 '- fData : const HistType']),
    shape='record',
    fontsize=20,
)

g.layout(prog='dot')
g.draw('analysis_uml.svg', prog='dot')
