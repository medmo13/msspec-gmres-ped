# from __builtin__ import False


class Plug(object):
    """
    A Plug represents an input or output attribute on a dataflow operator
    
    :param __operator: the operator that contains this plug
    :type __operator: dataflow.Operator
    """
    
    class PlugSide(object):
        INPUT = 1
        OUTPUT = 2

    def __init__(self, attr, data_type, plug_side, operator):
        """Constructor
        
        :param attr: an instance of Attribute
        :type attr: msspec.dataflow.Attribute
        :param data_type: the type of the attribute
        :type data_type: IDataType
        :param plug_side: tells if this plug is an input plug or an output plug
        :type plug_side: Plug.PlugSide enum
        :param operator: the operator that contains this plug
        :type operator: msspec.dataflow.Operator
        """
        self._attr = attr
        self._data_type = data_type
        self._plug_side = plug_side
        self._operator = operator
        self._outgoing_wires = []
        self._incoming_wire = None
        self._value = None
        
    @property
    def name(self):
        """
        :rtype: str
        """
        return self._attr.name
    
    @property
    def id(self):
        """
        :rtype: str
        """
        return '%s.%s' % (self.operator.name, self.name)
        
    @property
    def data_type(self):
        return self._data_type
        
    @property
    def operator(self):
        """
        :rtype: msspec.dataflow.Operator
        """
        return self._operator
    
    @property
    def is_pluggable(self):
        return self._attr.is_pluggable
    
    def get_value(self):
        print("getting value of plug %s" % self.id)
        assert self.value_is_available()
        if self.is_dirty():
            self._update_value()
            assert not self.is_dirty()
        return self._value

    def set_value(self, value):
        """
        :type value: any type of data that is compatible with this plug's datatype
        """
        assert isinstance(value, self._data_type.get_python_class())
        if self.is_pluggable:
            assert self.is_dirty()  # we expect to only have to use this method when this plug is dirty
        self._value = value
    
    def value_is_available(self):
        """
        indicates if the value of this plug can be computed
        """
        if self.is_destination():
            if not self.is_pluggable:
                return True  # if this plug is not pluggable, then its value is available since it doesn't depend on the result of other operators
            # print('value_is_available: self.is_connected() = %s' % str(self.is_connected()))
            if self.is_connected():
                src_plug = self.get_incoming_plug()
                return src_plug.value_is_available()
            else:
                return False
        else:
            # this is an output plug; its value is available if all the operators's input plugs are available
            return self.operator.all_inputs_are_available()
        
    def set_dirty(self):
        """
        indicates that the value of this plug is obsolete and needs to be recomputed
        """
        print("setting plug %s as dirty" % self.id)
        self._value = None
        if self.is_output_plug():
            if self.is_connected():
                for connected_plug in self.get_outgoing_plugs():
                    connected_plug.set_dirty()
        else:
            # this is an operator's input plug. Propagate the dirtiness to its out plugs
            self.operator.set_dirty()
            
    def is_dirty(self):
        """
        indicates if the value needs to be computed
        """
        assert self.value_is_available()
        if not self.is_pluggable:
            return False  # a non-pluggable plug is never dirty, it always has a valid value
        return self._value is None
    
    def _update_value(self):
        print("updating value of plug %s" % self.id)
        assert self.value_is_available()
        assert self.is_dirty()
        if self.is_destination():
            if self.is_connected():
                self._value = self.get_incoming_plug().get_value()
            else:
                assert False  # if we are in this case, this would mean that the value is not available
        else:
            print("updating operator %s" % self._operator.name)
            self._operator.update()
        assert not self.is_dirty()
    
    def is_connected(self):
        if self.is_source():
            return len(self._outgoing_wires) != 0
        else:
            return self._incoming_wire is not None
    
    def is_source(self):
        return self._plug_side == Plug.PlugSide.OUTPUT

    def is_destination(self):
        return self._plug_side == Plug.PlugSide.INPUT

    def is_output_plug(self):
        return self._plug_side == Plug.PlugSide.OUTPUT

    def is_input_plug(self):
        return self._plug_side == Plug.PlugSide.INPUT

    def connect_to(self, other_plug):
        """
        :type other_plug: dataflow.Plug
        :rtype: dataflow.Wire
        """
        return self._operator.data_flow.create_wire(self, other_plug)

    def add_outgoing_wire(self, wire):
        """
        :type wire: dataflow.Wire
        """
        assert self.is_source()
        self._outgoing_wires.append(wire)
    
    @property
    def incoming_wire(self):
        assert self.is_destination()
        assert self.is_connected()
        return self._incoming_wire
    
    @incoming_wire.setter
    def incoming_wire(self, wire):
        """
        :type wire: dataflow.Wire
        """
        assert self.is_destination()
        self._incoming_wire = wire
        self.operator.set_dirty()  # one of the operator's input plugs has changed, therefore its output plugs values are obsolete

    @property
    def outgoing_wires(self):
        """
        :return list(dataflow.Wire):
        """
        assert self.is_source()
        assert self.is_connected()
        return self._outgoing_wires

    def detach_wire(self, wire):
        """
        :param msspec.dataflow.Wire wire:
        """
        assert self.is_connected()
        if self.is_destination():
            assert self._incoming_wire == wire
            self.incoming_wire = None
            self.set_dirty()  # the data on this plug is no longer available
        else:
            assert self.is_source()
            self._outgoing_wires.remove(wire)

    def get_incoming_plug(self):
        assert self.is_destination()
        assert self.is_connected()
        return self.incoming_wire.input_plug  # (scenegraph_group.setter seen as a method, see https://github.com/PyCQA/pylint/issues/870) pylint: disable=no-member

    def get_outgoing_plugs(self):
        """
        :rtype: list(dataflow.Plug)
        """
        assert self.is_source()
        assert self.is_connected()
        outgoing_plugs = []
        for wire in self._outgoing_wires:
            outgoing_plugs.append(wire.output_plug)
        return outgoing_plugs
