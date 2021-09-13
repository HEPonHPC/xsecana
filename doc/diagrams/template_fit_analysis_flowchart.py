import pygraphviz as pgv
from format_utils import *

g = pgv.AGraph(strict=False, directed=True)
g.graph_attr['compound'] = True



g.add_node('Fill', label='Histogram Fill Service', shape='cylinder')
g.add_node('I', label='Interfacing Layer', shape='parallelogram')
g.add_edge('Fill', 'I')

g.add_node('Data', shape='box')
g.add_edge('I', 'Data')

g.add_node(
    'NM',
    label='{<TF>TemplateFitSignalEstimator | <F>SimpleFlux | <E>SimpleEfficiency | <U>IterativeUnfold}',
    shape='record',
)
g.add_edge('I', 'NM', lhead='cluster_Nominal')
g.add_subgraph(['NM'], label='Nominal', name='cluster_Nominal')

g.add_node(
    'SM',
    label='{<TF>TemplateFitSignalEstimator | <F>SimpleFlux | <E>SimpleEfficiency | <U>IterativeUnfold}',
    shape='record',
)
g.add_subgraph(['SM'], label='Systematically Shifted', name='cluster_Shifted')
g.add_edge('I', 'SM', lhead='cluster_Shifted')

g.add_node('cov', label='Calculate \\nCovariance Matrix')
g.add_edge('NM', 'cov', tailport='TF')
g.add_edge('SM', 'cov', tailport='TF')

g.add_node('Fit')
g.add_edge('cov', 'Fit')
g.add_edge('Data', 'Fit')
g.add_edge('NM', 'Fit', tailport='TF')

g.add_node('SE', label='Signal Estimation', shape='box')
g.add_node('SEU', label='Signal Estimation \\nUncertainty', shape='box')
g.add_edge('Fit', 'SE')
g.add_edge('Fit', 'SEU')

g.add_node('NH', label='Nominal Hist')
g.add_node('SH', label='Shifted Hists')
g.add_edge('NM', 'NH', tailport='F')
g.add_edge('NM', 'NH', tailport='E')
g.add_edge('NM', 'NH', tailport='U')

g.add_edge('SM', 'SH', tailport='F')
g.add_edge('SM', 'SH', tailport='E')
g.add_edge('SM', 'SH', tailport='U')

g.add_node('QS', label='Quadrature Sum')
g.add_edge('div2', 'QS')
g.add_edge('SEU', 'QS')

g.add_node('div1', label='Divide', shape='diamond')
g.add_edge('NH', 'div1')
g.add_edge('SE', 'div1')

g.add_node('div2', label='Divide', shape='diamond')
g.add_edge('NH', 'div2')
g.add_edge('SH', 'div2')

g.add_node('R', label='Result + Uncertainty', shape='box')
g.add_edge('QS', 'R')
g.add_edge('div1', 'R')

g.add_subgraph(
    ['SH', 'cov', 'NH', 'Fit', 'SEU', 'SE', 'QS', 'div1', 'div2'],
    name='cluster_propagator',
    label='TemplateFitUncertaintyPropagator',
    style='rounded'
)

g.layout(prog='dot')
g.draw('template_fit_analysis_flowchart.pdf', prog='dot')