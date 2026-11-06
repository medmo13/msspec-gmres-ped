import wx
from wx.lib.agw import floatspin as fs
from msspecgui.dataflow import IOperatorCreator


# bug 778 on https://www.brainwy.com/tracker/PyDev/#
# autocomplete (code completion) doesn't work for list items inside self members
# Hi, The code below shows a case of code completion that doesn't work although it works on very close cases (see below)
#
# class Toto(object):
#     def dummy(self):
#         pass
#
#
# class Titi(object):
#     def __init__(self, l, s):
#         """
#         :param list(Toto) l:
#         :param str s:
#         """
#         self._l = l
#         self._s = s
#
#         s.capitalize()  # <-- autocomplete works on s
#         self._s.capitalize()  # <-- autocomplete works on self._s
#
#         l.pop()  # <-- autocomplete works on l
#         for p1 in l:
#             p1.dummy()  # <-- autocomplete works on p1
#
#         self._l.pop()  # <-- autocomplete works on self._l
#         for p2 in self._l:
#             p2.dummy()  # <-- autocomplete doesn't work on p2

class PlugsGuiFrame(wx.Dialog):
    """ The frame containing the graphical user interface for the given plugs
    """

    def __init__(self, parent, window_title, plugs):
        """
        :param wx.Widget parent: window that contains this window
        :param str window_title: window title
        :param list(msspecgui.dataflow.Plug) plugs: the plugs that this gui allows the user to modify
        :param msspecgui.dataflow.Plug plug: the plugs that this gui allows the user to modify
        """
        super(PlugsGuiFrame, self).__init__(parent, title=window_title)
        self._plugs = plugs
        self._widgets = {}  #: :type self._widgets: dict(str, wx.Control)
        self.initui()
        self.Centre()
        self.Show()

    def initui(self):
        """
        the constructor of the interface
        """
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.panel = wx.Panel(parent=self, id=1)
        main_params_sizer = wx.FlexGridSizer(0, 2, 10, 10)
        # print("PlugsGuiFrame.initui len(self._plugs) = %d" % len(list(self._plugs)))
        
        for plug in self._plugs:  #: :type plug: msspecgui.dataflow.Plug
            if not plug.is_pluggable:
                label = wx.StaticText(self.panel, label=plug.name)
                main_params_sizer.Add(label)
                widget = None
                # print("plug.data_type.get_type_id() = %s" % plug.data_type.get_type_id())
                if plug.data_type.get_type_id() == 'int':
                    widget = wx.SpinCtrl(self.panel)
                    widget.SetValue(plug.get_value())
                elif plug.data_type.get_type_id() == 'float':
                    widget = fs.FloatSpin(self.panel,
                                          value=plug.get_value(), size=(180, -1),
                                          min_val=0.1, max_val=10.00,
                                          increment=0.01,
                                          style=fs.FS_CENTRE)
                    widget.SetFormat("%f")
                    widget.SetDigits(6)
                else:
                    raise NotImplementedError
                main_params_sizer.Add(widget)
                self._widgets[plug.name] = widget

        ok_cancel_reset_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # reset_btn = wx.Button(self.panel, 2, label="Reset")
        # reset_btn.Bind(wx.EVT_BUTTON, self.on_reset)
        cancel_btn = wx.Button(self.panel, wx.ID_CANCEL, label="Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.on_close, id=wx.ID_CANCEL)
        self.ok_btn = wx.Button(self.panel, wx.ID_OK, label="OK")
        # self.ok_btn.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        self.ok_btn.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        
        # ok_cancel_reset_sizer.Add(reset_btn, flag=wx.ALL, border=5)
        ok_cancel_reset_sizer.Add(cancel_btn)
        ok_cancel_reset_sizer.Add(self.ok_btn, flag=wx.LEFT, border=5)

        self.vbox.Add(main_params_sizer, proportion=2, flag=wx.ALL | wx.EXPAND, border=15)
        
        self.vbox.Add(ok_cancel_reset_sizer, flag=wx.ALIGN_RIGHT | wx.ALL, border=5)
        self.panel.SetSizer(self.vbox)
        self.panel.Fit()
        self.Fit()

    def on_reset(self, event):
        pass

    def on_close(self, event):
        event.Skip()  # propagate the event so that the dialog closes

    def on_ok(self, event):
        # print("PlugsGuiFrame.on_ok len(self._plugs) = %d" % len(list(self._plugs)))
        for plug in self._plugs:  #: :type plug: msspecgui.dataflow.Plug
            # print("plug.data_type.get_type_id() = %s" % plug.data_type.get_type_id())
            if not plug.is_pluggable:
                new_value = None
                widget = self._widgets[plug.name]
                if plug.data_type.get_type_id() == 'int':
                    new_value = widget.GetValue()
                    # print("PlugsGuiFrame.on_ok : new_value (int) = %d" % new_value)
                elif plug.data_type.get_type_id() == 'float':
                    new_value = widget.GetValue()
                    # print("PlugsGuiFrame.on_ok : new_value (float) = %f" % new_value)
                plug.set_value(new_value)
        event.Skip()  # propagate the event so that the dialog closes


class OperatorGui(IOperatorCreator.IAction):
    """a graphic user interface to allow the user to modify an operator's parameters
    
    """
    def get_name(self):
        return 'properties via gui'
    
    def execute_on_operator(self, operator):
        """
        :param dataflow.operator.Operator operator: the operator related to this action
        """
        dialog = PlugsGuiFrame(None, window_title='operator %s' % operator.name, plugs=operator.get_input_plugs())
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            operator.data_flow.on_modified_operator(operator)
        dialog.Destroy()
