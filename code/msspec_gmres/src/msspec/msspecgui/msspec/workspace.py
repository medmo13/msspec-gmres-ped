'''
Created on Dec 23, 2015

@author: graffy
'''

import msspecgui.msspec.cluster.clusterflow
# import sys
import msspecgui.dataflowxmlserializer


class Workspace(object):
    '''
    classdocs
    '''

    SETTINGS_FILE_NAME = 'settings.xml'

    def __init__(self, workspace_path, is_new_workspace):
        '''Constructor
        
        :param workspace_path: the file path to the workspace directory
        :type workspace_path: str
        :param is_new_workspace: if True, the workspace is created, otherwise the workspace is loaded from it path
        
        '''
        self._workspace_path = workspace_path
        if not is_new_workspace:
            self.load()
        else:
            self._cluster_flow = msspecgui.msspec.cluster.clusterflow.ClusterFlow()
            debugging = False
            if debugging is True:
                serializer = msspecgui.dataflowxmlserializer.DataflowSerializer()
                serializer.save_dataflow(self._cluster_flow, '%s/%s' % (self._workspace_path, self.SETTINGS_FILE_NAME))
    
    def load(self):
        '''Loads the workspace
        '''
        print('Workspace.load')
        serializer = msspecgui.dataflowxmlserializer.DataflowSerializer()
        self._cluster_flow = msspecgui.msspec.cluster.clusterflow.ClusterFlow()
        serializer.load_dataflow('%s/%s' % (self._workspace_path, self.SETTINGS_FILE_NAME), self._cluster_flow)
    
    def save(self):
        serializer = msspecgui.dataflowxmlserializer.DataflowSerializer()
        serializer.save_dataflow(self._cluster_flow, '%s/%s' % (self._workspace_path, self.SETTINGS_FILE_NAME))
    
    def get_cluster_flow(self):
        return self._cluster_flow
    
    @property
    def workspace_path(self):
        return self._workspace_path
