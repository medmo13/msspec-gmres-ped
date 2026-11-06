'''
Created on Jan 18, 2016

@author: graffy
'''
from msspecgui import dataflow
from msspecgui.dataflow.datatypes import StringDataType
from msspecgui.dataflow.datatypes import FloatDataType
from msspecgui.dataflow.datatypes import BoolDataType
from msspecgui.dataflow.datatypes import IntDataType
from .datatypes import ClusterDataType
from .generators import Bulk
from .generators import Repeat
from .generators import AddAdsorbate


class ClusterFlow(dataflow.DataFlow):
    '''
    a dataflow of cluster generator nodes
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(ClusterFlow, self).__init__()
        self.register_operator_creator(Bulk.Creator())
        self.register_operator_creator(Repeat.Creator())
        self.register_operator_creator(AddAdsorbate.Creator())
        self.register_data_type(ClusterDataType())
        self.register_data_type(StringDataType())
        self.register_data_type(FloatDataType())
        self.register_data_type(BoolDataType())
        self.register_data_type(IntDataType())
