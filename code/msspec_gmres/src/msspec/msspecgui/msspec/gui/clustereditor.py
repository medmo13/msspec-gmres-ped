'''
Created on Jan 18, 2016

@author: graffy
'''

# standard libraries
import os
import functools
# import math
# import sys
# from fileinput import close

import ase.io
# 3rd party libraries
import wx
# import cairo
# import wx.lib.wxcairo

# local libraries
import msspecgui.datafloweditor as datafloweditor
import msspecgui.dataflow as dataflow
from msspecgui.dataflow import IDataType
from msspecgui.msspec.gui.bulk_interface import BulkFrame
from msspecgui.msspec.gui.viewmanager import IViewCreator
from msspecgui.datafloweditor import OperatorGui

# based on example https://github.com/wxWidgets/wxPython/blob/master/demo/Cairo.py


class BackgroundContextMenu(wx.Menu):
    '''
    the context menu that the user sees when right-clicking on the background of the dataflow
    
    '''
    def __init__(self, panel):
        """
        :type panel: msspecgui.msspec.gui.ClusterEditor
        """
        wx.Menu.__init__(self)
        self.m_panel = panel

        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        # self.add_cluster_modifier_menu_item = self.Append(wx.ID_ANY,"Add Cluster Modifier"," adds a cluster modifier in the cluster editor")

        self.add_cluster_modifier_menu = wx.Menu()
        self.add_cluster_modifier_menu_item = self.AppendMenu(wx.ID_ANY, "Add Cluster Modifier...", self.add_cluster_modifier_menu)

        # self.m_panel.Bind(wx.EVT_MENU, self.on_add_cluster_modifier, self.add_cluster_modifier_menu_item)
        
        data_flow = self.m_panel.get_cluster_flow()
        
        for creator in data_flow.get_operator_creators():
            operator_type_id = creator.get_operator_type_id()
            
            menu_item = self.add_cluster_modifier_menu.Append(wx.ID_ANY, operator_type_id, " adds a cluster modifier in the cluster editor")
            # self.m_panel.Bind(wx.EVT_MENU, functools.partial( creator.create_operator, dflow = data_flow), menu_item)
            self.m_panel.Bind(wx.EVT_MENU, functools.partial(self.on_add_cluster_modifier, operator_type_id=operator_type_id), menu_item)

    def on_add_cluster_modifier(self, event, operator_type_id):
        """Callback that is invoked when the user chooses to add a specific cluster modifier using the context menu
        
        :param operator_type_id: the type of cluster modifier to add (eg 'ipr.msspec.cutsphere')
        """
        self.m_panel.create_operator(operator_type_id)


class OperatorContextMenu(wx.Menu):
    '''
    the context menu that the user sees when right-clicking on an operator of the dataflow
    
    '''
    def __init__(self, panel, operator):
        """
        :type panel: msspecgui.msspec.gui.ClusterEditor
        :type operator: dataflow.operator.Operator
        """
        wx.Menu.__init__(self)
        self.m_panel = panel

        menu_item = self.Append(wx.ID_ANY, 'delete', "deletes this operator")
        self.m_panel.Bind(wx.EVT_MENU, functools.partial(self.on_delete_operator, operator=operator), menu_item)

        for action in operator.creator.get_actions():
            
            menu_item = self.Append(wx.ID_ANY, action.get_name(), " performs this action on the operator")
            # self.m_panel.Bind(wx.EVT_MENU, functools.partial( creator.create_operator, dflow = data_flow), menu_item)
            self.m_panel.Bind(wx.EVT_MENU, functools.partial(self.on_perform_operator_action, action=action, operator=operator), menu_item)
            
    def on_perform_operator_action(self, event, action, operator):
        """Callback that is invoked when the user chooses to add a specific cluster modifier using the context menu
        
        :param operator_type_id: the type of cluster modifier to add (eg 'ipr.msspec.cutsphere')
        :type action: dataflow.ioperatorcreator.IOperatorCreator.IAction
        :type operator: dataflow.operator.Operator
        """
        action.execute_on_operator(operator)

    def on_delete_operator(self, event, operator):
        """Callback that is invoked when the user chooses to delete an operator
        
        :param dataflow.Operator operator: the operator that needs to be deleted
        """
        data_flow = operator.data_flow
        data_flow.delete_operator(operator)


class WireContextMenu(wx.Menu):
    '''
    the context menu that the user sees when right-clicking on a wire of the dataflow
    
    '''
    def __init__(self, panel, wire):
        """
        :param msspecgui.msspec.gui.ClusterEditor panel:
        :param dataflow.Wire wire: the wire for which this context menu is created
        """
        wx.Menu.__init__(self)
        self.m_panel = panel
        
        #         # create a menu item for each wire action
        #         for wire_action in self.m_panel.get_cluster_flow().wire_actions:
        #             menu_item = self.Append(wx.ID_ANY, data_action.get_name(), " performs this action on the operator")
        #             self.m_panel.Bind(wx.EVT_MENU, functools.partial(self.on_perform_data_action, data_action=data_action, data=wire.input_plug.get_value()), menu_item)
        menu_item = self.Append(wx.ID_ANY, 'delete', "deletes this wire")
        self.m_panel.Bind(wx.EVT_MENU, functools.partial(self.on_delete_wire, wire=wire), menu_item)

        for data_action in wire.data_type.get_actions():
            menu_item = self.Append(wx.ID_ANY, data_action.get_name(), " performs this action on the operator")
            self.m_panel.Bind(wx.EVT_MENU, functools.partial(self.on_perform_data_action, data_action=data_action, data=wire.input_plug.get_value()), menu_item)
            
    #     def on_perform_wire_action(self, event, wire_action, wire):
    #         """Callback that is invoked when the user chooses to add a specific cluster modifier using the context menu
    #
    #         :param dataflow.Wire.IAction wire_action: the action to perform
    #         :param dataflow.Wire wire: the wire on which to perform the action
    #         """
    #         wire_action.execute_on_wire(wire)

    def on_perform_data_action(self, event, data_action, data):
        """Callback that is invoked when the user chooses to add a specific data action using the context menu
        
        :param dataflow.IDataType.IAction data_action: the action to perform on the given data
        :param data: the data on which to perform the action
        """
        assert isinstance(data, data_action.datatype.get_python_class())
        data_action.execute_on_data(data)

    def on_delete_wire(self, event, wire):
        """Callback that is invoked when the user chooses to delete a wire using the context menu
        
        :param dataflow.Wire wire: the wire that needs to be deleted
        """
        data_flow = wire.data_flow
        data_flow.delete_wire(wire)


def opj(path):
    """Convert paths to the platform-specific separator"""
    platform_path = apply(os.path.join, tuple(path.split('/')))
    # HACK: on Linux, a leading / gets lost...
    if path.startswith('/'):
        platform_path = '/' + platform_path
    return platform_path


class BulkGui(dataflow.IOperatorCreator.IAction):
    
    def get_name(self):
        return 'properties via gui'
    
    def execute_on_operator(self, operator):
        """
        :type operator: dataflow.operator.Operator
        """
        dialog = BulkFrame(None, title="create a single cell", bulk_operator=operator)
        print("execute_on_operator : before MainLoop")
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            print("execute_on_operator : signaling operator %d as modified" % operator.id)
            operator.data_flow.on_modified_operator(operator)
            print("OK")
        else:
            print("Cancel")
        dialog.Destroy()


class ExportCluster(IDataType.IAction):
    
    def __init__(self, cluster_flow):
        """
        :param msspecgui.cluster.ClusterFlow cluster_flow:
        """
        super(ExportCluster, self).__init__(cluster_flow.get_data_type('physics.atomscluster'))
    
    def get_name(self):
        return 'export cluster'
    
    def execute_on_data(self, data):
        """
        :param data: the data on which this action needs to be performed
        """
        assert isinstance(data, self.datatype.get_python_class())
        dlg = wx.FileDialog(None, message="Choose a file to store this cluster (the format of the file is defined by the file's extension)", defaultDir='', defaultFile='', style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            cluster_file_path = dlg.GetPath()
            ase.io.write(cluster_file_path, data)
        dlg.Destroy()


class ClusterEditor(datafloweditor.DataflowView):

    class Creator(IViewCreator):
        VIEW_TYPE_NAME = 'cluster editor'
        
        def __init__(self, workspace):
            """
            :param msspecgui.msspec.Workspace workspace: the workspace that is associated with the cluster editors created by this creator
            """
            self._workspace = workspace
        
        @property
        def view_type_name(self):
            """
            :return str:
            """
            return self.VIEW_TYPE_NAME
        
        def create_view(self, parent):
            """
            :param wx.Window parent: the wx.Window that owns the view
            :return wx.Panel:
            """
            return ClusterEditor(parent, self._workspace.get_cluster_flow())

    def __init__(self, parent, cluster_flow):
        """
        :param cluster_flow: the dataflow that this editor manipulates
        :type cluster_flow: msspec.cluster.clusterflow.ClusterFlow
        
        """
        super(ClusterEditor, self).__init__(parent, cluster_flow, -1)
            
        for operator_creator in cluster_flow.get_operator_creators():
            if operator_creator.get_operator_type_id() == 'ase.lattice.bulk':
                operator_creator.register_operator_action(BulkGui())
            else:
                operator_creator.register_operator_action(OperatorGui())

        datatype = cluster_flow.get_data_type('physics.atomscluster')
        datatype.register_datatype_action(ExportCluster(cluster_flow))
        
    def create_operator(self, operator_type_id):
        operator = self.get_cluster_flow().create_operator(operator_type_id)
        self.get_cluster_flow().add_operator(operator)

    def get_cluster_flow(self):
        """
        :rtype: msspec.cluster.clusterflow.ClusterFlow
        """
        return self.dataflow

    def on_background_context_menu(self, event):
        '''
        called whenever the user right-clicks in the background of the dataflow
        '''
        pos = event.GetPosition()
        # print(pos)
        pos = self.ScreenToClient(pos)
        self.PopupMenu(BackgroundContextMenu(self), pos)

    def on_operator_context_menu(self, event, operator):
        '''
        called whenever the user right-clicks in an operator of the dataflow
        
        :type event: wx.Event
        :type operator: dataflow.operator.Operator
        '''
        pos = event.GetPosition()
        # print(pos)
        pos = self.ScreenToClient(pos)
        self.PopupMenu(OperatorContextMenu(self, operator), pos)

    def on_wire_context_menu(self, event, wire):
        '''
        called whenever the user right-clicks in a wire of the dataflow
        
        :type event: wx.Event
        :param dataflow.Wire wire:
        '''
        pos = event.GetPosition()
        # print(pos)
        pos = self.ScreenToClient(pos)
        self.PopupMenu(WireContextMenu(self, wire), pos)
