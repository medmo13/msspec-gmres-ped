import msspecgui
from msspecgui import scenegraph2d
# import msspecgui.datafloweditor as datafloweditor
import msspecgui.dataflow as dataflow
import widget


class WireWidget(widget.Widget):
    '''
    The widget representing a wire (a link connectiong two plugs)
    '''
    
    def __init__(self, data_flow_view):
        """
        :param data_flow_view: the dataflowview to which this operator belongs
        :type data_flow_view: datafloweditor.DataflowView
        
        """
        super(WireWidget, self).__init__(data_flow_view)
        # print('WireWidget constructor')
        self.main_shape_node = None
        self._source_plug_widget = None  #: :type self._source_plug_widget: msspecgui.datafloweditor.PlugWidget
        self._dest_plug_widget = None
        self.main_shape_node = scenegraph2d.Line()
        self.main_shape_node.fill = scenegraph2d.Color(0, 255, 255)
    
    @property
    def wire(self):
        return self.dest_plug_widget.plug.incoming_wire
    
    @property
    def source_plug_widget(self):
        return self._source_plug_widget

    @source_plug_widget.setter
    def source_plug_widget(self, plug_widget):
        self._source_plug_widget = plug_widget
        self._set_from_pos(plug_widget.centre_pos)
    
    @property
    def dest_plug_widget(self):
        return self._dest_plug_widget

    @dest_plug_widget.setter
    def dest_plug_widget(self, plug_widget):
        self._dest_plug_widget = plug_widget
        self._set_to_pos(plug_widget.centre_pos)
    
    def _set_from_pos(self, pos):
        (self.main_shape_node.x1, self.main_shape_node.y1) = pos
    
    def _set_to_pos(self, pos):
        (self.main_shape_node.x2, self.main_shape_node.y2) = pos

    def render_to_scene_graph(self, scenegraph_group):
        """
        :param scenegraph_group: the group node that contains the drawing of this element
        :type scenegraph_group: scenegraph.Group
        """
        self.scenegraph_group = scenegraph_group
        scenegraph_group.add_child(self.main_shape_node)
        
    def update_position(self):
        self._set_from_pos(self.source_plug_widget.centre_pos)
        self._set_to_pos(self.dest_plug_widget.centre_pos)
        
    def remove_from_scene_graph(self):
        parent = self.main_shape_node.parent
        parent.remove_child(self.main_shape_node)
        
    def is_construction_wire(self):
        """ Tells if this widget represents a wire being built interactively by the user
        """
        return self.source_plug_widget is None or self.dest_plug_widget is None
        
    def set_pointer_pos(self, pos):
        """
        :type pos: a tuple (x,y)
        """
        assert(self.is_construction_wire())  # this call only makes sens if this wire widget represents a wire in construction
        if self.source_plug_widget is None:
            self._set_from_pos(pos)
        else:
            self._set_to_pos(pos)
        
    def is_valid_final_plug_widget(self, plug_widget):
        """
        checks whether plug_widget is a valid plug to complete this connection
        
        :type plug_widget: datafloweditor.PlugWidget
        """
        assert(isinstance(plug_widget, msspecgui.datafloweditor.PlugWidget))
        assert(self.is_construction_wire())
        if self.source_plug_widget is None:
            return plug_widget.plug.is_output_plug()
        else:
            return plug_widget.plug.is_input_plug()
        
    def set_final_plug_widget(self, plug_widget):
        """
        completes this wire with the given plug_widget
        
        :type plug_widget: datafloweditor.PlugWidget
        """
        assert(isinstance(plug_widget, msspecgui.datafloweditor.PlugWidget))
        assert(self.is_construction_wire())
        assert(self.is_valid_final_plug_widget(plug_widget))
        if self.source_plug_widget is None:
            self.source_plug_widget = plug_widget
        else:
            self.dest_plug_widget = plug_widget
