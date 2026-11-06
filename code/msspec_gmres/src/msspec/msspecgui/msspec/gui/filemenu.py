'''
Created on Jan 4, 2016

@author: graffy
'''

import wx
import msspecgui.msspec.workspace


class FileMenu(wx.Menu):
    def __init__(self, main_frame):
        wx.Menu.__init__(self)
        self._main_frame = main_frame

        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        self._new_workspace_menu_item = self.Append(wx.ID_ANY, "N&ew workspace", " create a new workspace")
        self._open_workspace_menu_item = self.Append(wx.ID_ANY, "O&pen workspace", " load an existing workspace")
        self.AppendSeparator()
        self._save_workspace_menu_item = self.Append(wx.ID_ANY, "S&ave workspace", " save the current workspace")
        self._close_workspace_menu_item = self.Append(wx.ID_ANY, "C&lose workspace", " close the current workspace")
        self.AppendSeparator()
        self._about_menu_item = self.Append(wx.ID_ABOUT, "&About", " Information about this program")
        self.AppendSeparator()
        self._exit_menu_item = self.Append(wx.ID_EXIT, "E&xit", " Terminate the program")
    
        self._main_frame.Bind(wx.EVT_MENU, self.on_about, self._about_menu_item)
        self._main_frame.Bind(wx.EVT_MENU, self.on_exit, self._exit_menu_item)

        self._main_frame.Bind(wx.EVT_MENU, self.on_new_workspace, self._new_workspace_menu_item)
        self._main_frame.Bind(wx.EVT_MENU, self.on_open_workspace, self._open_workspace_menu_item)
        self._main_frame.Bind(wx.EVT_MENU, self.on_save_workspace, self._save_workspace_menu_item)
        self._main_frame.Bind(wx.EVT_MENU, self.on_close_workspace, self._close_workspace_menu_item)
        
        self.update()
    
    def update(self):
        self._close_workspace_menu_item.Enable(self._main_frame.m_workspace is not None)
        self._save_workspace_menu_item.Enable(self._main_frame.m_workspace is not None)
 
    def on_new_workspace(self, event):
        if self.on_close_workspace(event=None):
            default_path = ''
            dlg = wx.DirDialog(self._main_frame, "Choose location to store this new workspace", default_path)
            if dlg.ShowModal() == wx.ID_OK:
                workspace_path = dlg.GetPath()
                self._main_frame.set_workspace(msspecgui.msspec.workspace.Workspace(workspace_path, is_new_workspace=True))
               
                # f = open(, 'r')
                # self.control.SetValue(f.read())
                # f.close()
            dlg.Destroy()
        self.update()

    def on_open_workspace(self, event):
        if self.on_close_workspace(event=None):
            default_path = ''
            dlg = wx.DirDialog(self._main_frame, "Choose the workspace location", default_path)
            if dlg.ShowModal() == wx.ID_OK:
                workspace_path = dlg.GetPath()
                workspace = msspecgui.msspec.workspace.Workspace(workspace_path, is_new_workspace=False)
                self._main_frame.set_workspace(workspace)
                
                # f = open(, 'r')
                # self.control.SetValue(f.read())
                # f.close()
            dlg.Destroy()
        self.update()
    
    def on_save_workspace(self, event):
        print('FileMenu.on_save_workspace')
        self._main_frame.save_workspace()
        self.update()

    def on_close_workspace(self, event):
        self._main_frame.close_workspace()
        self.update()
        return True

    def on_about(self, event):
        # A message dialog box with an OK button. wx.OK is a standard ID in wxWidgets.
        dlg = wx.MessageDialog(self._main_frame, "A multiple-scattering package for spectroscopies using electrons to probe materials", "MsSpec 1.1", wx.OK)
        dlg.ShowModal()  # Show it
        dlg.Destroy()  # finally destroy it when finished.

    def on_exit(self, event):
        self.on_close_workspace(event=event)
        self._main_frame.Close(True)  # Close the frame.
