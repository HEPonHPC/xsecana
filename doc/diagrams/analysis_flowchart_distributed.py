import pygraphviz as pgv

G = pgv.AGraph(strict=False, directed=True)

G.add_node('Histogram Fill Service', shape='cylinder')
G.add_node('Interfacing Layer', shape='parallelogram')
G.add_edge('Histogram Fill Service', 'Interfacing Layer')

G.add_node('Measurements', shape='box3d')
G.add_edge('Interfacing Layer', 'Measurements')

G.add_node('Histograms', shape='box3d')
G.add_edge('Interfacing Layer', 'Histograms')
#--------

G.add_subgraph(['Measurements'], label='Analysis', name='cluster_Analysis')
G.add_subgraph(['Histograms'], label='Data', name='cluster_data')


#-------------------------

G.add_node('Uncertainty Propagation')
G.add_edge('Measurements', 'Uncertainty Propagation', label=' Reduce')
G.add_edge('Histograms', 'Uncertainty Propagation', label=' Reduce')
#--------------------------

G.add_node('Result + Uncertainty', shape='box')
G.add_edge('Uncertainty Propagation', 'Result + Uncertainty')

G.layout(prog='dot')
G.draw('analysis_flowchart_distributed.pdf', prog='dot')

