from msspecgui.dataflow import IDataType
from ase import Atoms


class ClusterDataType(IDataType):
    def __init__(self):
        IDataType.__init__(self)
        
    def get_type_id(self):
        """
        Returns a string uniquely identifying this data type
        """
        return 'physics.atomscluster'

    def get_python_class(self):  # pylint: disable=no-self-use
        """
        see IDataType.get_python_class
        """
        return Atoms
