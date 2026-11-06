from msspecgui import scenegraph2d
from .widget import Widget
# import dataflow


class PlugWidget(Widget):
    '''
    The widget representing a plug (input or output) of an operator in the dataflow

    '''
    
    DATA_IS_UNAVAILABLE_COLOR = scenegraph2d.Color(128, 64, 64)
    DATA_IS_AVAILABLE_COLOR = scenegraph2d.Color(64, 128, 64)
    HIGHLITED_COLOR = scenegraph2d.Color(128, 255, 128)
    
    def __init__(self, plug, operator_widget, data_flow_view):
        """
        :param plug: the dataflow plug that this widget graphically represents
        :type plug: dataflow.Plug
        @type plug: dataflow.Plug
        :param operator_widget: the operator widget that this plug widget is attached to
        :type operator_widget: datafloweditor.OperatorWidget
        :param data_flow_view: the dataflowview to which this operator belongs
        :type data_flow_view: datafloweditor.DataflowView
        
        """
        super(PlugWidget, self).__init__(data_flow_view)
        self.operator_widget = operator_widget
        self._plug = plug
        self.main_shape_node = None
    
    @property
    def plug(self):
        """
        :rtype: dataflow.plug.Plug
        """
        # plug = dataflow.plug.Plug
        return self._plug
    
    @property
    def centre_pos(self):
        """ Returns the position of the centre of this plug widget
        """
        # print('PlugWidget.centre_pos : local_to_parent_matrix = %s' % str(self.scenegraph_group.local_to_parent_matrix()))
        # defs = []
        # print('PlugWidget.centre_pos : self.scenegraph_group = %s' % str(self.scenegraph_group._xml(defs)))
        # defs = []
        # print('PlugWidget.centre_pos : self.scenegraph_group.parent = %s' % str(self.scenegraph_group.parent._xml(defs)))
        
        return self.scenegraph_group.local_to_world_matrix().project(x=0, y=0)
    
    def is_connectable(self):
        if self._plug.is_connected():
            return False
        if not self._plug.value_is_available():
            return False
        return True
    
    def update_appearance(self, mouse_is_above):
        color = {False: self.DATA_IS_UNAVAILABLE_COLOR,
                 True: self.DATA_IS_AVAILABLE_COLOR}[self._plug.value_is_available()]
        if mouse_is_above and self.is_connectable():
            color = self.HIGHLITED_COLOR
        
        self.main_shape_node.fill = color
    
    def render_to_scene_graph(self, scenegraph_group):
        """
        :param scenegraph_group: the group node that contains the drawing of this element
        :type scenegraph_group: scenegraph.Group
        """

        self.scenegraph_group = scenegraph_group
        circle_node = scenegraph2d.Circle()
        circle_node.cx = 0.0
        circle_node.cy = 0.0
        circle_node.r = 7.0
        self.main_shape_node = circle_node
        scenegraph_group.add_child(circle_node)
        self.update_appearance(mouse_is_above=False)

    def on_hover(self, is_entering):
        self.update_appearance(mouse_is_above=is_entering)
