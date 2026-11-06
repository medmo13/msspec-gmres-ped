'''
Created on Jan 4, 2016

@author: graffy
'''
import wx
from .filemenu import FileMenu
# import clustermenu
from msspecgui.msspec.gui.viewmanager import ViewFactory
from msspecgui.msspec.gui.viewmanager import ViewManager
from msspecgui.msspec.gui.clustereditor import ClusterEditor
from msspecgui.msspec.gui.clusterview import ClusterView


class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=(300, 300))
        # self.c  ontrol = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.CreateStatusBar()  # A Statusbar in the bottom of the window

        self.view_manager = None
        self.set_workspace(None)

        self._file_menu = FileMenu(self)
        # self.m_clusterMenu = clustermenu.ClusterMenu(self)
 
        # Creating the menubar.
        menu_bar = wx.MenuBar()
        menu_bar.Append(self._file_menu, "&File")  # Adding the "fileMenu" to the MenuBar
        # menu_bar.Append(self.m_clusterMenu,"&Cluster") # Adding the "fileMenu" to the MenuBar
        self.SetMenuBar(menu_bar)  # Adding the MenuBar to the Frame content.

        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.Show(True)  # Show the frame.

    def set_workspace(self, workspace):
        """
        :type workspace: Workspace
        """
        self.m_workspace = workspace
        window_title = 'MsSpecGui'
        if workspace is not None:
            window_title += ' - ' + workspace.workspace_path
        self.SetTitle(window_title)

        if workspace is None:
            if self.view_manager is not None:
                self.view_manager.Destroy()
                self.view_manager = None
        else:
            view_factory = ViewFactory()
            view_factory.register_view_creator(ClusterEditor.Creator(self.get_workspace()))
            view_factory.register_view_creator(ClusterView.Creator(self.get_workspace().get_cluster_flow()))
            self.view_manager = ViewManager(self, view_factory)
            self.view_manager.views[0].set_view_type(ClusterEditor.Creator.VIEW_TYPE_NAME)
            view2 = self.view_manager.add_view()
            view2.set_view_type(ClusterView.Creator.VIEW_TYPE_NAME)
        self.Hide()  # without this call, I couldn't make the workspace window updated (at least on osx); as a result, I had to drag the corner of the window to force a refresh
        self.Show()

    def get_workspace(self):
        return self.m_workspace

    def close_workspace(self, ask_confirmation=True):
        if self.m_workspace is not None:
            dlg = wx.MessageDialog(self, message="Do you want to save your workspace before closing it ?",
                                   caption="Save current workspace ?",
                                   style=wx.YES_NO | wx.CANCEL | wx.YES_DEFAULT | wx.ICON_QUESTION)
            clicked_button_id = dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            if clicked_button_id == wx.ID_YES:
                self.save_workspace()
            elif clicked_button_id == wx.ID_CANCEL:
                return False
        self.set_workspace(None)
        
    def save_workspace(self):
        self.m_workspace.save()

    def on_close(self, event):
        print("MainWindow.on_close")
        self.close_workspace(ask_confirmation=True)
        event.Skip()  # propagate the event so that the dialog closes
