# import plug
# from .dflow import DataFlow
from .plug import Plug


class Operator(object):
    '''
    an abstract operator that generates outputs from its inputs
    '''

    def __init__(self, data_flow, creator):
        '''
        Constructor
        
        :param data_flow: the data_flow that will contain this operator
        :type data_flow: DataFlow
        :param creator: an instance of IOperatorCreator, the creator of this object (consider it as an object that stores informations that are common to all operators of this type).
        :type creator: dataflow.ioperatorcreator.IOperatorCreator
        '''
        self._data_flow = data_flow
        self._creator = creator
        self._id = data_flow.get_new_operator_id()
        self._input_plugs = {}
        self._output_plugs = {}
        for attr in self._creator.get_input_attributes():
            data_type = self._data_flow.get_data_type(attr.type_name)
            assert data_type is not None
            p = Plug(attr, data_type, Plug.PlugSide.INPUT, self)
            self._input_plugs[p.name] = p
        for attr in self._creator.get_output_attributes():
            data_type = self._data_flow.get_data_type(attr.type_name)
            assert data_type is not None
            p = Plug(attr, data_type, Plug.PlugSide.OUTPUT, self)
            self._output_plugs[p.name] = p

    @property
    def name(self):
        return '%s_%d' % (self._creator.get_operator_type_id(), self._id)

    def get_input_plugs(self):
        """
        :rtype: list(msspecgui.dataflow.Plug)
        """
        return list(self._input_plugs.itervalues())

    def get_output_plugs(self):
        """
        :rtype: list(msspecgui.dataflow.Plug)
        """
        return list(self._output_plugs.itervalues())

    def get_plug(self, plug_name):
        """
        :type plug_name: str
        :rtype: msspecgui.dataflow.Plug
        """
        if plug_name in self._input_plugs:
            return self._input_plugs[plug_name]
        if plug_name in self._output_plugs:
            return self._output_plugs[plug_name]

    @property
    def wires(self):
        """
        :return list(msspecgui.dataflow.Wire)
        """
        wires = []
        for p in self._input_plugs.itervalues():
            if p.is_connected():
                wires.append(p.incoming_wire)
        for p in self._output_plugs.itervalues():
            if p.is_connected():
                for wire in p.outgoing_wires:
                    wires.append(wire)
        return wires

    def set_dirty(self):
        """
        indicates that the output values of this operator are obsolete and need to be recomputed
        """
        print("setting operator %d as dirty" % self.id)
        for plug in self.get_output_plugs():
            plug.set_dirty()

    def update(self):
        """
        compute this operator's output values from its input values
        
        this method is supposed to be overriden in derived classes
        """
        assert False  # should be implemented in derived classes

    def all_inputs_are_available(self):
        for input_plug in self.get_input_plugs():
            if not input_plug.value_is_available():
                return False
        return True

    @property
    def id(self):
        """
        :rtype: int
        """
        return self._id

    @id.setter
    def id(self, operator_id):
        """changes the unique identifier of this operator (use at your own risk)
        
        :param int operator_id: the new unique identifier for this operator
        """
        self._id = operator_id
    
    @property
    def creator(self):
        return self._creator
    
    @property
    def data_flow(self):
        return self._data_flow
