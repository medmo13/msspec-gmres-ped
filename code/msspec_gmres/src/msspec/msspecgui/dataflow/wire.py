'''
Created on Jul 1, 2016

@author: graffy
'''


class Wire(object):
    """
    a wire connects an input plug to an output plug
    """
    
    def __init__(self, input_plug, output_plug):
        """
        :type input_plug: dataflow.Plug
        :type output_plug: dataflow.Plug
        """
        self.input_plug = input_plug
        self.output_plug = output_plug

    @property
    def data_type(self):
        """
        :return IDataType: the datatype of this wire
        """
        return self.input_plug.data_type

    @property
    def data_flow(self):
        """
        :return msspecgui.dataflow.DataFlow
        """
        return self.input_plug.operator.data_flow
