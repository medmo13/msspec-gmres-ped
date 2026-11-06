from __future__ import print_function
import abc
import functools
import wx.lib.agw.aui as aui
import wx


class IViewCreator(object):

    @abc.abstractproperty
    def view_type_name(self):
        """
        :return str:
        """
        pass

    @abc.abstractmethod
    def create_view(self, parent):
        """
        :param wx.Window parent: the wx.Window that owns the view
        :return wx.Panel:
        """
        pass


class ViewFactory(object):

    def __init__(self):
        self._view_creators = {}  #: :type self._view_creators: dict(str, IViewCreator)

    def register_view_creator(self, view_creator):
        """
        :param IViewCreator view_creator:
        """
        self._view_creators[view_creator.view_type_name] = view_creator

    def get_view_type_names(self):
        return list(self._view_creators.iterkeys())

    def create_view(self, view_type_name, container_window):
        """
        :param str view_type_name: the type of view to create
        :param wx.Window container_window: the parent wx window
        """
        return self._view_creators[view_type_name].create_view(container_window)


class View(wx.Panel):
    """
    panel who will be used to create the interface
    """

    def __init__(self, view_manager, view_id):
        """
        :param ViewManager view_manager: the manager that manages this view
        :param int id: a number uniquely identifying this view in its view manager
        """
        wx.Panel.__init__(self, view_manager, -1)
        self.view_manager = view_manager
        self.id = view_id

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(vbox)
        self.panel = None

    @property
    def view_factory(self):
        return self.view_manager.view_factory

    def set_view_type(self, view_type):
        """
        :param str view_type:
        """
        self.on_view_set(None, view_type)

    def on_view_set(self, evt, view_type_name):
        """
        called when the menu item that assigns a view type to this view is clicked by the user
        
        :param str view_type_name: the type of view that the user just chose for this view
        """

        pane = self.view_manager.aui_manager.GetPaneByWidget(self)

        index = self.view_manager.get_free_pane_index(view_type_name)

        pane_name = view_type_name + ("" if index == 0 else " " + str(index))
        pane.caption = pane_name

        if self.panel is not None:
            self.panel.Destroy()
            self.panel = None

        self.panel = self.view_factory.create_view(view_type_name, self)
        self.GetSizer().Add(self.panel, proportion=1, flag=wx.EXPAND)  # the new panel is the only item in the vbox sizer that is allowed to change its size, so any value other than 0 will do

#        self.Update()
#        self.Refresh()
#        self.view_manager.Update()
        # self.view_manager.Refresh()
        self.view_manager.aui_manager.Update()
        self.view_manager.update_view_menu()
        # print self.selected_crystal_str


class ViewManager(wx.Panel):
    """a user interface that contains views
    """

    def __init__(self, parent, view_factory):
        """
        :param ViewFactory view_factory: the factory that is used to create views
        """
        super(ViewManager, self).__init__(parent)
        self.view_factory = view_factory

        self.parent = parent

        # manages the panes associated with it for a particular wxFrame
        self.aui_manager = aui.AuiManager(self)
        self.last_created_view_id = 0
        self.views = []
        self.view_set_menu = None  # :param wxMenu self.view_set_menu: the menu View->Set that allows the user to assign a view type to a view

        view = self.add_view(False, False)  # @UnusedVariable the 1st view is special : it's not closable nor movable (to avoid holes)

        self.aui_manager.Update()
        self.update_view_menu()

    def pane_exists(self, panes, pane_name_prefix, pane_index):
        for pane_info in panes:
            if pane_info.caption == "%s %d" % (pane_name_prefix, pane_index):
                return True
        return False

    def get_free_pane_index(self, pane_name_prefix):
        """
        finds a non-used integer to uniquely identify a pane that we want to name "<pane_name_prefix> <i>"
        
        :param str pane_name_prefix: the pane name prefix (usually the view type name)
        """
        index = 1
        while self.pane_exists(self.aui_manager.GetAllPanes(), pane_name_prefix, index):
            index += 1
        return index

    def add_view(self, movable=True, enable_close_button=True):
        """
        :param bool movable: if True, the view is movable
        :param bool enable_close_button: if True a close button is added to the view
        :return View:
        """
        self.last_created_view_id += 1
        view = View(self, self.last_created_view_id)

        index = self.get_free_pane_index("Empty")

        caption = 'Empty %d' % (index)
        pane_info = aui.AuiPaneInfo().Movable(movable).Floatable(False).CloseButton(enable_close_button).Caption(caption).DestroyOnClose(True)
        if view.id % 2 == 1:
            pane_info = pane_info.Center()
        else:
            pane_info = pane_info.Right()
        pane_info = pane_info.Layer(1)
        self.aui_manager.AddPane(view, pane_info)

        self.aui_manager.Bind(aui.EVT_AUI_PANE_CLOSE, self.on_close_pane)

        self.aui_manager.OnFloatingPaneResized(view, None)
        self.aui_manager.Update()
        self.views.append(view)
        self.update_view_menu()
        return view

    def on_close_pane(self, event):
        """
        :param wx.Event event:
        """
        pane = event.GetPane()
        self.views.remove(pane.window)
        items = self.view_set_menu.GetMenuItems()
        for i in items:
            if i.GetText() == pane.caption:
                self.view_set_menu.DestroyItem(i)

    def update_view_menu(self):
        """
        builds and updates the view menu to reflect the views in the view manager
        """
        def populate_view_set_submenu(view_set_submenu, views, aui_manager, view_type_names, window_receiving_events):
            """
                populate the View/Set menu
                
                :param wxMenu view_set_submenu: the menu View/Set
                :param list(View) views: the list of existing views in the view manager
                :param wxAuiManager aui_manager: the aui manager that manages the views (which are widgets)
                :param list(str) view_type_names: the list of available view types
                :param wxWindow window_receiving_events: the window that receives the menu events
            """
            for view in views:
                view_item = wx.Menu()
                view_set_submenu.AppendSubMenu(view_item, aui_manager.GetPaneByWidget(view).caption)
                for view_type_name in view_type_names:
                    menu_item_id = wx.NewId()
                    view_type_item = wx.MenuItem(view_item, menu_item_id, view_type_name)
                    view_item.AppendItem(view_type_item)
                    window_receiving_events.Bind(wx.EVT_MENU, functools.partial(view.on_view_set, view_type_name=view_type_name), view_type_item)
        
        # create the menu bar if it doesn't exist yet
        if self.parent.GetMenuBar():
            menu_bar = self.parent.GetMenuBar()
        else:
            menu_bar = wx.MenuBar()

        # delete the View menu if it already exists, as we are going to re-create it
        i = 0
        while i < menu_bar.GetMenuCount():
            if menu_bar.GetMenuLabelText(i) == "View":
                menu_bar.Remove(i)
            i += 1

        # create the view menu
        menu_view = wx.Menu()
        self.view_set_menu = wx.Menu()

        # create the menu item that allows the user to add a view to the view manager
        add_view_menu_item = wx.MenuItem(menu_view, wx.NewId(), "Add")
        menu_view.AppendItem(add_view_menu_item)
        # create the sub menu that allows the user to assign a view type to a view
        menu_view.AppendSubMenu(self.view_set_menu, "Set")
        menu_bar.Append(menu_view, "View")
        self.parent.Bind(wx.EVT_MENU, self.add_view, add_view_menu_item)
        self.parent.SetMenuBar(menu_bar)

        populate_view_set_submenu(self.view_set_menu, self.views, self.aui_manager, self.view_factory.get_view_type_names(), self.parent)

    """
    def close_win(self, event):
        self.Close()

    def maximize(self, event):
        self.Maximize(True)

    def minimize(self, event):
        self.Maximize(False)
    """
