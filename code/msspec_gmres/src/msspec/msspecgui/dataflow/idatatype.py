import abc


class IDataType(object):
    
    class IAction(object):
        """an abstract class that represents an action that can be performed on a data of the given data type
        
        for example, the export cluster is an action (with a gui that allows the user to select a file path) that can be performed on a data of type physics.atomscluster. This mechanism is used in the dataflow editor : when the user right-clicks on a wire, a list of actions related to the wire's datatype is offered to the user.
        """
        def __init__(self, datatype):
            """
            :param msspecgui.dataflow.IDataType datatype: the datatype on which this action is expected to perform
            """
            self.datatype = datatype
        
        @abc.abstractmethod
        def get_name(self):
            """ returns a name that uniquely identifies the action amongst actions associated with a given IDataType
            """
            assert False
        
        @abc.abstractmethod
        def execute_on_data(self, data):
            """executes this action on the given data
            
            :param data: the data on which this action is performed
            """
            assert False
    
    def __init__(self):
        self._actions = {}  #: :type self._actions: dict(str, IDataType.IAction)
    
    def get_type_id(self):  # pylint: disable=no-self-use
        """
        Returns a string iuniquely identifying this data type
        :rtype: str
        """
        assert False  # this method is expected to be defined in a derived class

    def get_python_class(self):  # pylint: disable=no-self-use
        """
        returns the python class associated with this data type
        """
        assert False  # this method is expected to be defined in a derived class

    def register_datatype_action(self, action):
        """
        
        :type action: IDataType.IAction
        """
        assert action.datatype == self
        self._actions[action.get_name()] = action

    def get_actions(self):
        """
        :rtype: list(IDataType.IAction)
        """
        return list(self._actions.itervalues())
