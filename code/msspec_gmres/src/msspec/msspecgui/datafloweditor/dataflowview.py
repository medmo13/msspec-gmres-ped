"""
a graphical editor for dataflow
"""
from __future__ import print_function
# import wx
import wx.lib.wxcairo
import cairo
# from wx.lib.scrolledpanel import ScrolledPanel
# import scenegraph2d.xml
import msspecgui.scenegraph2d.cairo

# import dataflow
import msspecgui.dataflow as dataflow
# import msspecgui.datafloweditor as datafloweditor

from .plugwidget import PlugWidget
from .wirewidget import WireWidget
from msspecgui.datafloweditor.dotlayoutmanager import DotLayoutManager
# class DataflowView(wx.ScrolledWindow):


class DataflowView(wx.Panel):
    '''
    a dataflow Graphical User Interface
    
    :type self.wire_being_created_widget: WireWidget
    '''

    # see https://wiki.wxwidgets.org/Scrolling

    class DataflowEventsHandler(dataflow.DataFlow.IDataFlowEventsHandler):
        """
        handles dataflow events
        """
        def __init__(self, data_flow_view):
            """
            :param msspec.datafloweditor.DataflowView data_flow_view:
            """
            super(DataflowView.DataflowEventsHandler, self).__init__()
            self.data_flow_view = data_flow_view
            
        def on_added_operator(self, operator):
            super(DataflowView.DataflowEventsHandler, self).on_added_operator(operator)
            # a new operator has just been added to the dataflow
            # create a new widget to graphically manipulate the added operator
            self.data_flow_view.add_operator_widget(msspecgui.datafloweditor.OperatorWidget(operator, self.data_flow_view))

        def on_deleted_operator(self, operator):
            super(DataflowView.DataflowEventsHandler, self).on_deleted_operator(operator)
            self.data_flow_view.delete_widget_for_operator(operator)
            self.data_flow_view.update_appearance()  # possibly update the appearance of the plugs that have now become available
            
        def on_added_wire(self, wire):
            super(DataflowView.DataflowEventsHandler, self).on_added_wire(wire)
            # create a new widget to graphically represent and manipulate the added wire
            wire_widget = self.data_flow_view.add_wire_widget(wire)  # @UnusedVariable
            self.data_flow_view.update_appearance()  # possibly update the appearance of the plugs that have now become available

        def on_deleted_wire(self, wire):
            super(DataflowView.DataflowEventsHandler, self).on_deleted_wire(wire)
            self.data_flow_view.delete_widget_for_wire(wire)
            self.data_flow_view.update_appearance()  # possibly update the appearance of the plugs that have now become available

    def __init__(self, parent, data_flow, log):
        """
        :type parent: wx.Window
        :type data_flow: DataFlow
        """
        wx.Panel.__init__(self, parent, -1)
        self.scenegraph_group_to_widget = {}  # this array associates a datafloweditor.Widget instance to a scenegraph group that represents this widget
        self.layout_is_enabled = True

        self.log = log
        self.layout_manager = DotLayoutManager()
        self.scene = None  # the 2d scenegraph that describes the graphical aspect of the dataflow
        self.scene = msspecgui.scenegraph2d.Group()
        # with open('/Users/graffy/ownCloud/ipr/msspec/rectangle.svg') as f:
        #     self.scene = scenegraph2d.xml.parse(f.read())
        self.operator_to_widget = {}
        self.wire_to_widget = {}

        self.dataflow = data_flow

        # self.selected_widget = None
        self.hovered_widget = None  # the widget hovered on by the mouse pointer
        self.plug_being_connected = None  # when the user creates a wire, memorizes the first selected plug
        self.wire_being_created_widget = None  # while the use creates a wire, the widget that is used to represent it
        self.is_left_down = False

        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.Bind(wx.EVT_LEFT_DOWN, self.on_left_down)
        self.Bind(wx.EVT_LEFT_UP, self.on_left_up)
        self.Bind(wx.EVT_MOTION, self.on_move)
        self.Bind(wx.EVT_CONTEXT_MENU, self.on_context_menu)

        self.data_flow_view_updater = DataflowView.DataflowEventsHandler(self)
        self.layout_is_enabled = False  # temporarily disable the layout manager for efficiency (prevent complete layout computation for each operator in the data_flow)
        data_flow.add_dataflow_events_handler(self.data_flow_view_updater)
        # initialize the widgets to reflect the state of the cluster flow
        for node in self.dataflow.operators:
            self.data_flow_view_updater.on_added_operator(node)

        for wire in self.dataflow.wires:
            self.data_flow_view_updater.on_added_wire(wire)
        self.layout_is_enabled = True
        self.update_operators_position()

    def on_paint(self, evt):  # IGNORE:unused-argument
        """handler for the wx.EVT_PAINT event
        """
        if self.IsDoubleBuffered():
            display_context = wx.PaintDC(self)
        else:
            display_context = wx.BufferedPaintDC(self)
        display_context.SetBackground(wx.Brush('white'))
        display_context.Clear()
        
        self.render(display_context)

    def update_appearance(self):
        for operator_widget in self.operator_to_widget.itervalues():
            operator_widget.update_appearance()

    def update_operators_position(self):
        if self.layout_is_enabled:
            op_pos = self.layout_manager.compute_operators_position(self.dataflow)
            dodgy_offset = (20.0, 20.0)  # TODO: make it better
            dodgy_scale = 1.5
            for op, pos in op_pos.iteritems():
                x, y = pos
                print("update_operators_position : %f %f" % (x, y))
                if op in self.operator_to_widget:  # while loading the dataflow, there are operators that don't have widgets yet
                    op_widget = self.operator_to_widget[op]
                    op_widget.set_position(x * dodgy_scale + dodgy_offset[0], y * dodgy_scale + dodgy_offset[1])

    def add_operator_widget(self, operator_widget):
        """
        :type operator_widget: an instance of OperatorWidget
        """
        self.operator_to_widget[operator_widget.operator] = operator_widget
        widget_root = msspecgui.scenegraph2d.Group()
        self.scene.add_child(widget_root)
        # widget_root.transform = [msspecgui.scenegraph2d.Translate(operator_widget.get_id() * 100.0 + 50.0, 60.0)]
        operator_widget.render_to_scene_graph(widget_root)
        self.update_operators_position()
        # self.UpdateWindowUI(wx.UPDATE_UI_RECURSE)
        # print("self.UpdateWindowUI has executed")
        # self.Refresh()

    def delete_widget_for_operator(self, operator):
        """
        :param Operator operator: the operator for which we want to delete the widget
        """
        op_widget = self.get_operator_widget(operator)
        op_widget.remove_from_scene_graph()
        
        self.operator_to_widget.pop(operator)

    def add_wire_widget(self, wire):
        """
        creates a wire widget and adds it to this dataflow view
        :param wire: the wire represented by the widget to create
        :type wire: Wire
        :rtype: WireWidget
        """
        wire_widget = WireWidget(self)
        wire_widget.source_plug_widget = self.get_plug_widget(wire.input_plug)
        wire_widget.dest_plug_widget = self.get_plug_widget(wire.output_plug)
        self.wire_to_widget[wire] = wire_widget
        widget_root = msspecgui.scenegraph2d.Group()
        self.scene.add_child(widget_root)
        wire_widget.render_to_scene_graph(widget_root)
        self.update_operators_position()
        return wire_widget

    def delete_widget_for_wire(self, wire):
        """
        :param Wire wire: the wire for which we want to delete the widget
        """
        wire_widget = self.get_wire_widget(wire)
        wire_widget.remove_from_scene_graph()
        # widget_root = wire_widget.parent()
        # self.scene.remove_child(widget_root)
        
        self.wire_to_widget.pop(wire)
        self.update_operators_position()

    def get_wire_widget(self, wire):
        """
        Returns the widget associated with the given wire
        :param Wire wire:
        :return WireWidget:
        """
        return self.wire_to_widget[wire]

    def get_operator_widget(self, operator):
        """
        Returns the operator widget associated with the given operator
        :param Operator operator:
        :return OperatorWidget:
        """
        return self.operator_to_widget[operator]

    def get_plug_widget(self, plug):
        """
        Returns the plug widget associated with the given plug
        :type plug: Plug
        :rtype: PlugWidget
        """
        operator_widget = self.get_operator_widget(plug.operator)
        plug_widget = operator_widget.get_plug_widget(plug)
        return plug_widget

    def render(self, display_context):
        """
        renders this dataflow view in the given drawing context
        
        :param display_context: the drawing context in which to render the dataflow
        :type display_context: a wx context
        """
        #cairo_context = wx.lib.wxcairo.ContextFromDC(display_context)
        w, h = self.GetClientSize()
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
        cairo_context = cairo.Context(surface)

        msspecgui.scenegraph2d.cairo.render_scene(self.scene, cairo_context)

        bitmap = wx.lib.wxcairo.BitmapFromImageSurface(surface)
        display_context.DrawBitmap(bitmap, 0, 0)
        
    def get_pointed_widget(self, pointer_position):
        """
        returns the widget that is at the position pointer_position
        """
        hits = self.scene.pick(pointer_position[0], pointer_position[1])
        if len(hits) > 0:
            svg_path, local_pos = hits[-1]  # @UnusedVariable
            # print("hit svg node stack : %s" % str(svg_path))
            for svg_node in reversed(svg_path):
                if svg_node in self.scenegraph_group_to_widget:
                    return self.scenegraph_group_to_widget[svg_node]
        
        # for widget in self.operator_to_widget.itervalues():
        #     if widget.get_bounding_box().Contains(pointer_position):
        #         return widget
        return None

#     def select_widget(self, widget):
#         if self.selected_widget is not None:
#             self.selected_widget.set_selected_state(False)
#         self.selected_widget = widget
#         self.selected_widget.set_selected_state(True)

    def on_left_down(self, event):
        pos = event.GetPositionTuple()
        display_context = wx.ClientDC(self)
        display_context.DrawCircle(pos[0], pos[1], 5)
        widget = self.get_pointed_widget(pos)
        if widget is not None:
            # self.select_widget(widget)
            if isinstance(widget, PlugWidget):
                plug_widget = widget
                if plug_widget.is_connectable():
                    self.plug_being_connected = plug_widget
                    wire_widget = WireWidget(self)
                    if plug_widget.plug.is_source():
                        wire_widget.source_plug_widget = plug_widget
                    else:
                        wire_widget.dest_plug_widget = plug_widget
                    self.wire_being_created_widget = wire_widget
                    self.wire_being_created_widget.render_to_scene_graph(self.scene)
        self.is_left_down = True
        self.Refresh()

    def on_left_up(self, event):
        self.is_left_down = False
        self.plug_being_connected = None

        if self.wire_being_created_widget is not None:
            pos = event.GetPositionTuple()
            pointed_widget = self.get_pointed_widget(pos)
            if pointed_widget is not None and isinstance(pointed_widget, PlugWidget):
                if self.wire_being_created_widget.is_valid_final_plug_widget(pointed_widget):
                    self.wire_being_created_widget.set_final_plug_widget(pointed_widget)
                    self.dataflow.create_wire(self.wire_being_created_widget.source_plug_widget.plug, self.wire_being_created_widget.dest_plug_widget.plug)
            self.wire_being_created_widget.remove_from_scene_graph()
            self.wire_being_created_widget = None
        self.Refresh()

    def on_move(self, event):
        pos = event.GetPositionTuple()
        if self.wire_being_created_widget is not None:
            self.wire_being_created_widget.set_pointer_pos(pos)
        widget = self.get_pointed_widget(pos)
        # print('widget at (%d,%d) : %s\n' % (pos[0], pos[1], widget))
        if self.hovered_widget is not None:
            if self.hovered_widget != widget:
                self.hovered_widget.on_hover(False)
                self.hovered_widget = None
        if widget is not None:
            if self.hovered_widget is None:
                widget.on_hover(True)
                self.hovered_widget = widget
        
        if self.is_left_down:
            display_context = wx.ClientDC(self)
            display_context.DrawCircle(pos[0], pos[1], 3)
        self.Refresh()

    def on_context_menu(self, event):
        pos = self.ScreenToClient(event.GetPosition())

        for wire_widget in self.wire_to_widget.itervalues():
            if wire_widget.get_bounding_box(border=5).Contains(pos):
                self.on_wire_context_menu(event, wire_widget.wire)
                return

        for operator_widget in self.operator_to_widget.itervalues():
            if operator_widget.get_bounding_box().Contains(pos):
                self.on_operator_context_menu(event, operator_widget.operator)
                return
        
        self.on_background_context_menu(event)

    def on_operator_context_menu(self, event, operator):
        '''
        called whenever the user right-clicks in an operator of the dataflow
        '''
        pass  # possibly implemented in derived classes

    def on_wire_context_menu(self, event, wire):
        '''
        called whenever the user right-clicks in an operator of the dataflow
        :param msspec.dataflow.Wire wire: the wire on which the context menu is supposed to act
        '''
        pass  # possibly implemented in derived classes

    def on_background_context_menu(self, event):
        '''
        called whenever the user right-clicks in the background of the dataflow
        '''
        pass  # possibly implemented in derived classes
