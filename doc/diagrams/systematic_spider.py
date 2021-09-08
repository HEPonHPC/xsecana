import pygraphviz as pgv

g = pgv.AGraph()

g.add_node('Systematic<T>')
g.graph_attr['overlap'] = False

g.add_node('1-sided\\n(1 T)', shape='box3d')
g.add_edge('Systematic<T>', '1-sided\\n(1 T)')

g.add_node('2-sided\\n(2 T)', shape='box3d')
g.add_edge('Systematic<T>', '2-sided\\n(2 T)')

g.add_node('Multiverse\\n(N T)', shape='box3d')
g.add_edge('Systematic<T>', 'Multiverse\\n(N T)')

g.layout(prog='neato')
g.draw('systematic_spider.svg', prog='neato')