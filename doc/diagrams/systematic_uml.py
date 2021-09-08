import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False)
g.add_node(
    0,
    label=uml_formatter('Systematic<T>',
                        methods=['+ Invoke(f : class F, args : class ... Args) \l: Systematic<std::invoke_result<F, T, Args>>',
                                 '+ Eval(data : const HistType) : HistType',
                                 '+ NSigmaShift(nsigma : double, nominal : T, args ... : Args) \l: HistType',
                                 '- BinSigma(nsigma : double, \luniverses : std::vector<ScalarType>, \lnominal & ScalarType) \l: ScalarType'],

                        members=['- fContainer : std::vector<T>',
                                 '- fName : std::string'
                                 '- fType : SystType_t']),
    shape='record',
    fontsize=20,
)

g.layout(prog='dot')
g.draw('systematic_uml.svg', prog='dot')
