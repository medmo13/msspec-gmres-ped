from __future__ import print_function
import json
import tempfile
# from pprint import pprint
import pydot
from .ilayoutmanager import ILayoutManager

# def pydot_sample():
#     import pydot  # import pydot or you're not going to get anywhere my friend :D
#
#     # first you create a new graph, you do that with pydot.Dot()
#     graph = pydot.Dot(graph_type='graph')
#
#     # the idea here is not to cover how to represent the hierarchical data
#     # but rather how to graph it, so I'm not going to work on some fancy
#     # recursive function to traverse a multidimensional array...
#     # I'm going to hardcode stuff... sorry if that offends you
#
#     # let's add the relationship between the king and vassals
#     for i in range(3):
#         # we can get right into action by "drawing" edges between the nodes in our graph
#         # we do not need to CREATE nodes, but if you want to give them some custom style
#         # then I would recomend you to do so... let's cover that later
#         # the pydot.Edge() constructor receives two parameters, a source node and a destination
#         # node, they are just strings like you can see
#         edge = pydot.Edge("king", "lord%d" % i)
#         # and we obviosuly need to add the edge to our graph
#         graph.add_edge(edge)
#
#     # now let us add some vassals
#     vassal_num = 0
#     for i in range(3):
#         # we create new edges, now between our previous lords and the new vassals
#         # let us create two vassals for each lord
#         for j in range(2):
#             edge = pydot.Edge("lord%d" % i, "vassal%d" % vassal_num)
#             graph.add_edge(edge)
#             vassal_num += 1
#
#     # ok, we are set, let's save our graph into a file
#     graph.write_png('/tmp/example1.png') # .write_png('example1_graph.png')
#
#     # and we are done!

# def graphviz_sample():
#     from graphviz import Digraph
#
#     dot = Digraph(comment='The Round Table')
#     dot.node('A', 'King Arthur')
#     dot.node('B', 'Sir Bedevere the Wise')
#     dot.node('L', 'Sir Lancelot the Brave')
#
#     dot.edges(['AB', 'AL'])
#     dot.edge('B', 'L', constraint='false')
#
#     # print(dot.source)
#     # // The Round Table
#     # digraph {
#     #     A [label="King Arthur"]
#     #     B [label="Sir Bedevere the Wise"]
#     #     L [label="Sir Lancelot the Brave"]
#     #         A -> B
#     #         A -> L
#     #         B -> L [constraint=false]
#     # }
#
#     dot.render('/tmp/round-table.png', view=True)


class DotLayoutManager(ILayoutManager):
    """A layout manager that uses graphviz's dot
    """

    def __init__(self):
        super(DotLayoutManager, self).__init__()
        self._method = 'plain'  # at the moment we use the plain method instead of the more robust json method because json method is too recent to be supported on all platforms (on windows, graphviz 2.38 doesn't yet supports json output format)
        self._tmp_file_path = tempfile.mktemp(suffix='.%s' % self._method, prefix='msspec_layout_')
        print('DotLayoutManager.__init__ : self._tmp_file_path = %s' % self._tmp_file_path)

    def compute_operators_position(self, data_flow):
        """
        :param DataFlow data_flow:
        """
        # pydot_sample()
        graph = pydot.Dot(graph_type='graph')
        
        # node [shape = record,height=.1];
        
        # attr_to_label = {}
        
        graph.set_node_defaults(shape='record')
        for op in data_flow.operators:
            node = pydot.Node(op.name)
            label = '"'
            plug_index = 0
            for plug in op.get_input_plugs() + op.get_output_plugs():
                if plug.is_pluggable:
                    if len(label) > 1:
                        label += '|'
                    # label += "<%s> %s" % (plug.name, plug.name)
                    label += "<%s> " % (plug.name)
                    plug_index += 1
            label += '"'
            node.set('label', label)
            graph.add_node(node)
            
        for wire in data_flow.wires:
            #: :type wire: Wire
            src_plug = wire.input_plug
            dst_plug = wire.output_plug
            
            src_op = src_plug.operator
            dst_op = dst_plug.operator

            # print("src_op.name = %s" % src_op.name)
            # print("dst_op.name = %s" % dst_op.name)
            edge = pydot.Edge('"%s":%s' % (src_op.name, src_plug.name), '"%s":%s' % (dst_op.name, dst_plug.name))
            graph.add_edge(edge)
            
        # print(graph.to_string())
        # graph.write('/tmp/toto.pdf', format='pdf')
        func = {'plain': self._graph_to_pos_plain, 'json': self._graph_to_pos_json}[self._method]
        return func(graph, data_flow)
    
    def _graph_to_pos_plain(self, graph, data_flow):
        """extracts operator positions from the pydot graph using dot's plain output format
        """
        class GraphField(object):
            X_SIZE = 2
            Y_SIZE = 3
        
        class NodeField(object):
            NAME = 1
            X = 2
            Y = 3

        graph.write(self._tmp_file_path, format='plain')
        
        # example output
        # graph 1 1.9444 2.5417
        # node "ase.lattice.bulk_1" 0.375 2.2847 0.75 0.51389 "<output_cluster> " solid record black lightgrey
        # node "ase.repeat_2" 0.56944 1.2708 0.75 0.51389 "<input_cluster> |<output_cluster> " solid record black lightgrey
        # node "ase.lattice.bulk_3" 1.5694 1.2708 0.75 0.51389 "<output_cluster> " solid record black lightgrey
        # node "ase.build.add_adsorbate_4" 1.3056 0.25694 0.83333 0.51389 "<slab> |<adsorbate> |<output_cluster> " solid record black lightgrey
        # edge "ase.lattice.bulk_1" "ase.repeat_2" 4 0.375 2.0301 0.375 1.8816 0.375 1.6906 0.375 1.5208 solid black
        # edge "ase.repeat_2" "ase.build.add_adsorbate_4" 4 0.76389 1.0208 0.76389 0.76408 1.0278 0.76369 1.0278 0.50694 solid black
        # edge "ase.lattice.bulk_3" "ase.build.add_adsorbate_4" 4 1.4344 1.0167 1.3695 0.87155 1.3056 0.68373 1.3056 0.50694 solid black
        # stop
        
        positions = {}  #: :type positions: dict(str, (float,float))

        scale = 100.0
        y_max = None
        with open(self._tmp_file_path) as data_file:
            for line in data_file.readlines():
                # print(line)
                tokens = line.split(' ')
                if tokens[0] == 'graph':
                    y_max = float(tokens[GraphField.Y_SIZE])
                elif tokens[0] == 'node':
                    op_name = tokens[NodeField.NAME][1:-1]
                    x = float(tokens[NodeField.X])
                    y = float(tokens[NodeField.Y])
                    # print( op_name, x, y)
                    positions[data_flow.find_operator(op_name)] = ((y_max - y) * scale, x * scale)  # graphviz' dot lays th flow of arcs in the down direction by default, and we want it on the right
                elif tokens[0] == 'edge':
                    break  # for efficiency reason, stop parsing as we're not interested in the rest of the file
        return positions
    
    def _graph_to_pos_json(self, graph, data_flow):
        """extracts operator positions from the pydot graph using dot's json output format
        """
        graph.write('toto.json', format='json')
        
        with open('toto.json') as data_file:
            data = json.load(data_file)
        # pprint(data)

        y_max = float(data["bb"].split(",")[3])

        positions = {}  #: :type positions: dict(str, (float,float))
        if "objects" in data:
            for op_node in data["objects"]:
                op_name = op_node["name"]
                # print(op_node["name"])
                # print(op_node["pos"])
                coords = op_node["pos"].split(',')
                x = float(coords[0])
                y = float(coords[1])
                # print( op_name, x, y)
                positions[data_flow.find_operator(op_name)] = (y_max - y, x)  # graphviz' dot lays th flow of arcs in the down direction by default, and we want it on the right
        
        return positions
