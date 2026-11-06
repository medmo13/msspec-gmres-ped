

class Attribute(object):
    def __init__(self, attr_name, attr_type_name, attr_desc, is_pluggable=True):
        """
        :param attr_name: the name of this attribute, eg 'file_name'
        :type attr_name: str
        :param attr_type_name: the id of the type of this attribute (eg 'string'). Note that this type of attribute should have been registered in the dataflow first
        :type attr_type_name: str
        :param attr_desc: a string describing the role of this attribute
        :type attr_desc: str
        """
        self._attr_name = attr_name
        self._attr_type_name = attr_type_name
        self._attr_desc = attr_desc
        self.is_pluggable = is_pluggable  # indicates that this attribute can be plugged to another one using a wire. If not, then the value of this attribute is set using other means (eg graphical user interface)
    
    @property
    def name(self):
        return self._attr_name
    
    @property
    def type_name(self):
        return self._attr_type_name
