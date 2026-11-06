
from .attribute import Attribute
import abc


class IOperatorCreator(object):
    '''abstract base class that allows registration of operator creators
    '''
    
    class IAction(object):
        """an abstract class that represents an action that can be performed on an operator type
        
        for example, the ase.lattice.bulk operator's gui (property settings) is an action that can be performed on an operator of type ase.lattice.bulk
        """
        def __init__(self):
            pass
        
        @abc.abstractmethod
        def get_name(self):
            """ returns a name that uniquely identifies the action amongst actions associated with a given IOperatorCreator
            """
            assert False
        
        @abc.abstractmethod
        def execute_on_operator(self, operator):
            """executes this action on the given operator
            
            :param operator: the operator on which this action is performed
            :type operator: dataflow.Operator
            """
            assert False
    
    def __init__(self):
        self._input_attrs = {}
        self._output_attrs = {}
        self._actions = {}
    
    @abc.abstractmethod
    def get_operator_type_id(self):  # pylint: disable=no-self-use
        '''returns the unique id of the type of the operator that this creator creates (eg ipr.msspec.simpleclustercreator)
        '''
        assert False  # this method should be implemented in a derived class
    
    @abc.abstractmethod
    def create_operator(self, dflow):  # pylint: disable=no-self-use
        '''creates an application specific operator in the given Dataflow dflow
        
        :type dflow: DataFlow
        '''
        assert False  # this method should be implemented in a derived class

    def register_operator_action(self, action):
        """
        
        :type action: dataflow.ioperatorcreator.IOperatorCreator.IAction
        """
        self._actions[action.get_name()] = action

    def get_actions(self):
        """
        :rtype: list(dataflow.ioperatorcreator.IOperatorCreator.IAction)
        """
        return list(self._actions.itervalues())

    def add_input_attribute(self, attr_name, attr_type_name, attr_desc, is_pluggable=True):
        """Declares a new input attribute for the related operator
        
        :param attr_name: the name of the attribute
        :type attr_name: str
        :param attr_type_name: the name of the type of the attribute (eg 'math.position3')
        :type attr_type_name: str
        :param attr_desc: the description of the attribute
        :type attr_desc: str
        
        :rtype: Attribute
        """
        attr = Attribute(attr_name, attr_type_name, attr_desc, is_pluggable=is_pluggable)
        self._input_attrs[attr_name] = attr
        return attr

    def add_output_attribute(self, attr_name, attr_type_name, attr_desc, is_pluggable=True):
        """Declares a new output attribute for the related operator
        
        :param attr_name: the name of the attribute
        :type attr_name: str
        :param attr_type_name: the name of the type of the attribute (eg 'math.position3')
        :type attr_type_name: str
        :param attr_desc: the description of the attribute
        :type attr_desc: str
        
        :rtype: Attribute
        """
        attr = Attribute(attr_name, attr_type_name, attr_desc, is_pluggable=is_pluggable)
        self._output_attrs[attr_name] = attr
        return attr

    def get_input_attributes(self):
        return self._input_attrs.itervalues()
            
    def get_output_attributes(self):
        return self._output_attrs.itervalues()
