import pygraphviz as pgv

G = pgv.AGraph(strict=False, directed=True)

G.add_node('Histogram Fill Service', shape='cylinder')
G.add_node('Interfacing Layer', shape='parallelogram')
G.add_edge('Histogram Fill Service', 'Interfacing Layer')

G.add_node('Nominal Measurement', shape='box')
G.add_node('Systematically-shifted Measurements', shape='box3d')
G.add_node('Data1', shape='box')
G.add_node('Data2', shape='box')

G.add_subgraph(
    ['Nominal Measurement',
     'Systematically-shifted Measurements'],
    name='cluster_Analysis',
    label='Analysis'
)

G.add_edge('Interfacing Layer', 'Nominal Measurement')
G.add_edge('Interfacing Layer', 'Systematically-shifted Measurements')
G.add_edge('Interfacing Layer', 'Data1')
G.add_edge('Interfacing Layer', 'Data2')
#-------------------------

G.add_node('Uncertainty Propagation')

G.add_edge('Nominal Measurement', 'Uncertainty Propagation')
G.add_edge('Systematically-shifted Measurements', 'Uncertainty Propagation')
G.add_edge('Data1', 'Uncertainty Propagation')
G.add_edge('Data2', 'Uncertainty Propagation')
#--------------------------

G.add_node('Result + Uncertainty', shape='box')
G.add_edge('Uncertainty Propagation', 'Result + Uncertainty')

G.layout(prog='dot')
G.draw('analysis_flowchart_ratio.svg', prog='dot')

