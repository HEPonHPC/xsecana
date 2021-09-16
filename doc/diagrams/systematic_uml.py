import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False)
g.add_node(
    0,
    label=uml_formatter('Systematic<T>',
                        methods=['+ Eval(args... : class ... Args) : HistType'],
                        members=['- fContainer : std::vector<T>',
                                 '- fName : std::string'
                                 '- fType : SystType_t']),
    shape='record',
    fontsize=20,
)

g.add_node(1,
           label=uml_formatter('Systematic<IEfficiency>',
                               methods='+ Eval() : HistType')

)

g.layout(prog='dot')
g.draw('systematic_uml.pdf', prog='dot')
