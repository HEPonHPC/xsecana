import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False, directed=True)

g.add_node('SM', label='Systematically-shifted \\nMeasurements', shape='box3d')
g.add_node('NM', label='Nominal Measurement', shape='box')
g.add_node('Data', shape='box')

g.add_node('SE', label='Apply(Measurement::Eval)', shape='parallelogram')
g.add_edge('SM', 'SE')
g.add_edge('Data', 'SE')

g.add_node('NE', label=' Measurement::Eval', shape='parallelogram')
g.add_edge('NM', 'NE')
g.add_edge('Data', 'NE')

g.add_node('SH', label='Systematically-shifted \\nHistograms', shape='box3d')
g.add_node('NH', label='Nominal Histogram', shape='box')
g.add_edge('NE', 'NH')
g.add_edge('SE', 'SH')

g.add_node('sub', label='Subtract', shape='parallelogram')
g.add_edge('SH', 'sub')
g.add_edge('NH', 'sub')

g.add_node('div', label='Divide', shape='parallelogram')
g.add_edge('sub', 'div')
g.add_edge('NH', 'div')
g.add_node('FU', label='Fractional Uncertainty', shape='box')
g.add_edge('div', 'FU')

g.layout(prog='dot')
g.draw('simple_quad_sum_flowchart.pdf', prog='dot')

