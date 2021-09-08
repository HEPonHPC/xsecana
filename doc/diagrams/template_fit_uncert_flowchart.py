import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False, directed=True)

g.add_node('SSE', label='Signal Estimators', shape='box3d')
g.add_node('NSE', label='Signal Estimator', shape='box')
g.add_node('Data', shape='box')



g.layout(prog='dot')
g.draw('template_fit_uncert_flowchart.svg', prog='dot')