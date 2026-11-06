# import wx.lib.wxcairo
# import math
import msspecgui
from msspecgui import scenegraph2d
from .widget import Widget
# from msspecgui import datafloweditor
# from numpy.distutils.misc_util import cxx_ext_match


class OperatorWidget(Widget):
    '''
    The widget representing a node in the dataflow
    '''
    _g_last_id = 0
    
    def __init__(self, operator, data_flow_view):
        """
        :param operator: the dataflow operator that this widget graphically represents
        :type operator: msspecgui.dataflow.Operator
        :param data_flow_view: the dataflowview to which this operator belongs
        :type data_flow_view: msspecgui.datafloweditor.DataflowView
        
        """
        super(OperatorWidget, self).__init__(data_flow_view)
        self._operator = operator
        self._id = OperatorWidget.get_new_id()
        self.is_selected = False
        self.widget_background = None
        self.main_shape_node = scenegraph2d.Group()
        self.plug_to_widget = {}  # an associative array that gives the widget associated to each plug name
        
        for side in ['input', 'output']:
            plugs = {'input': self._operator.get_input_plugs(),
                     'output': self._operator.get_output_plugs()}[side]
            # plug_index = 0
            for plug in plugs:
                if plug.is_pluggable:
                    plug_widget = msspecgui.datafloweditor.PlugWidget(plug, self, data_flow_view)
                    self.plug_to_widget[plug.name] = plug_widget
                # plug_index += 1
    
    @classmethod
    def get_new_id(cls):
        OperatorWidget._g_last_id += 1
        return OperatorWidget._g_last_id
    
    def get_id(self):
        return self._id
    
    @property
    def operator(self):
        return self._operator

    def update_appearance(self):
        for plug_widget in self.plug_to_widget.itervalues():
            plug_widget.update_appearance(mouse_is_above=False)  # FIXME : handle the case where the mouse is over the plug widget

    def get_plug_widget(self, plug):
        """
        Returns the plug widget associated with the given plug
        :type plug: Plug
        :rtype: PlugWidget
        """
        plug_widget = self.plug_to_widget[plug.name]
        return plug_widget
    
    def set_selected_state(self, is_selected):
        self.is_selected = is_selected
        self.widget_background.fill = {False: scenegraph2d.Color(128, 128, 128), True: scenegraph2d.Color(192, 192, 128)}[self.is_selected]
    
    def set_position(self, x, y):
        self.main_shape_node.parent.transform = [msspecgui.scenegraph2d.Translate(x, y)]
        assert self.operator is not None
        for wire in self.operator.wires:
            if wire in self._data_flow_view.wire_to_widget:  # while load the dataflow, wires are not guaranteed to have a widget yet
                wire_widget = self._data_flow_view.wire_to_widget[wire]
                wire_widget.update_position()
    
    def render_to_scene_graph(self, scenegraph_group):
        """
        :param scenegraph_group: the group node that contains the drawing of this element
        :type scenegraph_group: scenegraph.Group
        """
        rect = scenegraph2d.Rectangle()
        rect.width = 70.0
        rect.height = 70.0
        rect.x = -35.0
        rect.y = -35.0
        rect.fill = scenegraph2d.Color(128, 128, 128)
        self.widget_background = rect
        self.main_shape_node.add_child(rect)
        title = "%s (%d)" % (self._operator.creator.get_operator_type_id().split('.')[-1], self.get_id())
        # import my.source.module
        # c = my.source.module.Color()
        title_text = scenegraph2d.Text(title)
        title_text.fill = scenegraph2d.Color.black
        title_text.text_anchor = "middle"
        self.main_shape_node.add_child(title_text)

        for side in ['input', 'output']:
            (cx, plugs) = {'input': (-25.0, self._operator.get_input_plugs()),
                           'output': (25.0, self._operator.get_output_plugs())}[side]
            plug_index = 0
            for p in plugs:
                if p.is_pluggable:
                    plug_group = scenegraph2d.Group()
                    plug_group.x = cx
                    plug_group.y = -25.0 + 10.0 * plug_index
                    self.main_shape_node.add_child(plug_group)
    
                    plug_widget = self.plug_to_widget[p.name]
                    plug_widget.render_to_scene_graph(plug_group)
                    plug_index += 1
        self.scenegraph_group = scenegraph_group
        self.scenegraph_group.add_child(self.main_shape_node)

    def remove_from_scene_graph(self):
        parent = self.main_shape_node.parent
        parent.remove_child(self.main_shape_node)
