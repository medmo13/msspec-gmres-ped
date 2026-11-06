from .idatatype import IDataType


class StringDataType(IDataType):
    def __init__(self):
        IDataType.__init__(self)
        
    def get_type_id(self):
        """
        Returns a string uniquely identifying this data type
        """
        return 'string'

    def get_python_class(self):  # pylint: disable=no-self-use
        """
        see IDataType.get_python_class
        """
        return str


class FloatDataType(IDataType):
    def __init__(self):
        IDataType.__init__(self)
        
    def get_type_id(self):
        """
        Returns a string uniquely identifying this data type
        """
        return 'float'

    def get_python_class(self):  # pylint: disable=no-self-use
        """
        see IDataType.get_python_class
        """
        return float


class BoolDataType(IDataType):
    def __init__(self):
        IDataType.__init__(self)
        
    def get_type_id(self):
        """
        Returns a string uniquely identifying this data type
        """
        return 'bool'

    def get_python_class(self):  # pylint: disable=no-self-use
        """
        see IDataType.get_python_class
        """
        return bool


class IntDataType(IDataType):
    def __init__(self):
        IDataType.__init__(self)
        
    def get_type_id(self):
        """
        Returns a string uniquely identifying this data type
        """
        return 'int'

    def get_python_class(self):  # pylint: disable=no-self-use
        """
        see IDataType.get_python_class
        """
        return int
