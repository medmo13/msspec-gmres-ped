# from __future__ import absolute_import
# from .plug import Plug
from .wire import Wire


class DataFlow(object):
    ''' a flow of operators, each of them having inputs and outputs that are connected together
    '''
    
    class IDataFlowEventsHandler(object):
        """an abstract class that listens and responds to events affecting a flow of operators
        
        """
        def on_added_operator(self, operator):
            """this method is invoked whenever an operator is added to the flow
            
            :param dataflow.Operator operator: the operator that has just been added to the flow
            """
            pass

        def on_deleted_operator(self, operator):
            """this method is invoked whenever an operator is deleted from the data flow
            
            :param dataflow.Operator operator: the operator that has just been deleted
            """
            pass

        def on_modified_operator(self, operator):
            """this method is invoked whenever an operator is modified in the flow
            
            :param operator: the operator that has just been modified
            :type operator: dataflow.Operator
            """
            pass
        
        def on_added_wire(self, wire):
            """this method is invoked whenever a wire is added to the flow
            
            :param wire: the wire that has just been added to the flow
            :type wire: dataflow.Wire
            """
            pass

        def on_deleted_wire(self, wire):
            """this method is invoked whenever a wire is deleted from the flow
            
            :param dataflow.Wire wire: the wire that has just been deleted
            """
            pass

    def __init__(self):
        self._last_created_operator_id = 0  # each operator (node) has a uniqueid in the flow
        self._operators = {}  #: :type self._operators: dict[int, msspec.dataflow.Operator]
        self._registered_creators = {}  # This factory lists and creates cluster creator instances
        self._registered_data_types = {}  #: :type self._registered_data_types: dict(str, msspecgui.dataflow.IDataType) ; the list of known data types (for links between nodes) in this dataflow
        self._dataflow_events_handlers = []
        self._wires = []
    
    @property
    def last_created_operator_id(self):
        return self._last_created_operator_id
    
    @last_created_operator_id.setter
    def last_created_operator_id(self, operator_id):
        """initializes the operator id generator
        
        useful for restoring the dataflow state, eg after loading
        
        :param int operator_id:
        """
        self._last_created_operator_id = operator_id

    def add_dataflow_events_handler(self, dataflow_events_handler):
        """adds a dataflow events handler to the list of dataflow events handler
        
        :type dataflow_events_handler: IDataFlowEventsHandler
        """
        self._dataflow_events_handlers.append(dataflow_events_handler)

    def remove_dataflow_events_handler(self, dataflow_events_handler):
        """removes a dataflow events handler to the list of dataflow events handler
        
        :type dataflow_events_handler: IDataFlowEventsHandler
        """
        self._dataflow_events_handlers.remove(dataflow_events_handler)
    
    def add_operator(self, operator):
        """
        :type operator: msspec.dataflow.Operator
        """
        self._operators[operator.id] = operator
        for dataflow_events_handler in self._dataflow_events_handlers:
            dataflow_events_handler.on_added_operator(operator)

    def delete_operator(self, operator):
        """
        :param msspec.dataflow.Operator operator:
        """
        data_flow = operator.data_flow
        
        # delete incoming wires
        for plug in operator.get_input_plugs():
            if plug.is_connected():
                wire = plug.incoming_wire
                data_flow.delete_wire(wire)

        # delete outgoing wires
        for plug in operator.get_output_plugs():  #: :type plug: Plug
            if plug.is_connected():
                for wire in plug.outgoing_wires:
                    data_flow.delete_wire(wire)
        
        self._operators.pop(operator.id)
        for dataflow_events_handler in self._dataflow_events_handlers:
            dataflow_events_handler.on_deleted_operator(operator)

    def on_modified_operator(self, operator):
        """tells the dataflow that the given operator has been modified (eg its paraeters have been changed)
        
        :type operator: Operator
        """
        # print("on_modified_operator : operator %d" % operator.id)
        operator.set_dirty()
        for deh in self._dataflow_events_handlers:
            deh.on_modified_operator(operator)

    def get_operator(self, operator_id):
        """
        :param int operator_id:
        :rtype: msspec.dataflow.Operator
        """
        return self._operators[operator_id]

    def find_operator(self, operator_name):
        """
        :param str operator_name:
        :rtype: msspec.dataflow.Operator
        """
        for op in self.operators:
            if op.name == operator_name:
                return op
        return None

    @property
    def operators(self):
        """
        :rtype: list(Operator)
        """
        return list(self._operators.itervalues())
    
    def get_new_operator_id(self):
        """
        :rtype: int
        """
        self._last_created_operator_id += 1
        return self._last_created_operator_id
    
    def register_operator_creator(self, operator_creator):
        '''
        
        :param operator_creator:
        :type operator_creator: derived from IOperatorCreator
        '''
        self._registered_creators[operator_creator.get_operator_type_id()] = operator_creator

    def register_data_type(self, data_type):
        '''
        
        :param data_type: the type of this argument is expected to be derived from IDataType
        :type data_type: derived from IDataType
        '''
        self._registered_data_types[data_type.get_type_id()] = data_type

    def get_data_type(self, data_type_id):
        """
        :return msspecggui.dataflow.IDataType:
        """
        return self._registered_data_types[data_type_id]

    def get_operator_creators(self):
        """
        :rtype: list(dataflow.ioperatorcreator.IOperatorCreator)
        """
        return list(self._registered_creators.itervalues())

    @property
    def wires(self):
        """
        :rtype: list(msspec.dataflow.Wire)
        """
        return self._wires

    def create_operator(self, operator_type_id):
        """ creates an operator of the given type and adds it to this flow
        
        :param operator_type_id: the type of operator to create. This type must be one of the allowed types for this flow (in other words, this type must have been registered to this flow)
        """
        creator = self._registered_creators[operator_type_id]
        return creator.create_operator(self)
    
    def create_wire(self, input_plug, output_plug):
        """
        :param Plug input_plug: the input plug
        :param Plug output_plug:
        :return dataflow.Wire:
        """
        wire = Wire(input_plug, output_plug)
        input_plug.add_outgoing_wire(wire)
        output_plug.incoming_wire = wire
        self._wires.append(wire)
        
        for dataflow_events_handler in self._dataflow_events_handlers:
            dataflow_events_handler.on_added_wire(wire)

        return wire

    def delete_wire(self, wire):
        """
        :param Wire wire: the wire that needs to be destroyed
        """
        wire.input_plug.detach_wire(wire)
        wire.output_plug.detach_wire(wire)

        self._wires.remove(wire)
        
        for dataflow_events_handler in self._dataflow_events_handlers:
            dataflow_events_handler.on_deleted_wire(wire)

        return wire
