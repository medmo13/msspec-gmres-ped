import abc
from xml.dom import minidom

# from msspecgui.dataflow import Operator
from msspecgui.dataflow import IDataflowSerializer
from msspecgui.dataflow.datatypes import StringDataType
from msspecgui.dataflow.datatypes import FloatDataType
from msspecgui.dataflow.datatypes import BoolDataType
from msspecgui.dataflow.datatypes import IntDataType


class IDataTypeSerializer(object):
    
    @abc.abstractmethod
    def get_data_type(self):
        """
        :return msspecgui.dataflow.IDatatype: the id of the type this seralize deals with
        """
        return None

    @abc.abstractmethod
    def value_to_xml(self, value, xml_node):
        """
        :param minidom.Element xml_node: the container node of the serialized value
        """
        pass
    
    @abc.abstractmethod
    def xml_to_value(self, xml_node):
        """
        :param minidom.Element xml_node: the container node of the serialized value
        """
        pass
    
    
class StringSerializer(IDataTypeSerializer):
    
    def get_data_type(self):
        return StringDataType()

    def value_to_xml(self, value, xml_node):
        xml_node.setAttribute('value', '%s' % value)
    
    def xml_to_value(self, xml_node):
        return xml_node.GetAttribute('value')


class FloatSerializer(IDataTypeSerializer):
    
    def get_data_type(self):
        return FloatDataType()

    def value_to_xml(self, value, xml_node):
        xml_node.setAttribute('value', '%f' % value)
    
    def xml_to_value(self, xml_node):
        return float(xml_node.GetAttribute('value'))


class BoolSerializer(IDataTypeSerializer):
    
    def get_data_type(self):
        return BoolDataType()

    def value_to_xml(self, value, xml_node):
        xml_node.setAttribute('value', {False: 'false', True: 'true'}[value])
    
    def xml_to_value(self, xml_node):
        value_as_str = xml_node.GetAttribute('value')
        return {'false': False, 'true': True}[value_as_str]


class IntSerializer(IDataTypeSerializer):
    
    def get_data_type(self):
        return IntDataType()

    def value_to_xml(self, value, xml_node):
        xml_node.setAttribute('value', '%d' % value)
    
    def xml_to_value(self, xml_node):
        return int(xml_node.GetAttribute('value'))


class DataflowSerializer(IDataflowSerializer):
    
    FORMAT_VERSION = 1
    
    def __init__(self):
        super(DataflowSerializer, self).__init__()
        self._data_type_serializers = {}  #: :type self._data_type_serializers: dict[str, IDataTypeSerializer]
        self._register_data_type_serializer(StringSerializer())
        self._register_data_type_serializer(FloatSerializer())
        self._register_data_type_serializer(BoolSerializer())
        self._register_data_type_serializer(IntSerializer())
    
    def _register_data_type_serializer(self, data_type_serializer):
        """
        :param IDataTypeSerializer data_type_serializer: the datatype serialize that needs to be registered
        """
        self._data_type_serializers[data_type_serializer.get_data_type().get_type_id()] = data_type_serializer
    
    def _get_datatype_serializer(self, data_type_id):
        """
        :param str data_type_id:
        """
        return self._data_type_serializers[data_type_id]
    
    def _operator_as_xml(self, operator, xml_doc):
        """creates the xml representation of the given operator
        
        :param msspecgui.dataflow.Operator operator: the operator
        :param minidom.Document xml_doc: the xml document that be used to create the xml node
        :return minidom.Element: the xml representation of the operator
        """
        op_xml_node = xml_doc.createElement('operator')
        op_xml_node.setAttribute('nodeId', '%d' % operator.id)
        op_xml_node.setAttribute('operatorTypeId', '%s' % operator.creator.get_operator_type_id())
        for plug in operator.get_input_plugs():
            if not plug.is_pluggable:
                plug_node = xml_doc.createElement(plug.name)
                data_type_id = plug.data_type.get_type_id()
                # print('data_type_id = %s' % data_type_id)
                data_type_serializer = self._get_datatype_serializer(data_type_id)
                data_type_serializer.value_to_xml(plug.get_value(), plug_node)
                op_xml_node.appendChild(plug_node)
        return op_xml_node

    def _create_operator_from_xml(self, op_xml_node, dataflow):
        """
        :param minidom.Element op_xml_node: the xml node describing the operator and its data
        :param msspec.dataflow.IDataFlow dataflow: the dataflow that contains the resulting operator
        
        :return msspec.dataflow.Operator: the created operator
        """
        operator_type_id = str(op_xml_node.getAttribute('operatorTypeId'))
        operator = dataflow.create_operator(operator_type_id)
        operator.id = int(op_xml_node.getAttribute('nodeId'))
        dataflow.add_operator(operator)
        return operator
    
    def save_dataflow(self, dataflow, file_path):
        """
        :type dataflow: msspecgui.dataflow.DataFlow

        """
        xml_doc = minidom.Document()
        root_xml_node = xml_doc.createElement('dataflow')
        xml_doc.appendChild(root_xml_node)

        # store a format version so that if the format needs changing, then it will be possible to detect and support old file formats
        format_version_xml_node = xml_doc.createElement('formatVersion')
        format_version_xml_node.setAttribute('value', '%d' % self.FORMAT_VERSION)
        root_xml_node.appendChild(format_version_xml_node)

        last_create_op_id_xml_node = xml_doc.createElement('lastCreatedOperatorId')
        last_create_op_id_xml_node.setAttribute('value', '%d' % dataflow.last_created_operator_id)
        root_xml_node.appendChild(last_create_op_id_xml_node)

        operators_xml_node = xml_doc.createElement('operators')
        root_xml_node.appendChild(operators_xml_node)

        for op in dataflow.operators:
            op_xml_node = self._operator_as_xml(op, xml_doc)
            operators_xml_node.appendChild(op_xml_node)
        
        wires_xml_node = xml_doc.createElement('wires')
        root_xml_node.appendChild(wires_xml_node)
        
        for wire in dataflow.wires:  #: :type wire: msspec.dataflow.Wire
            wire_xml_node = xml_doc.createElement('wire')
            
            wire_xml_node.setAttribute('fromOperator', '%d' % wire.input_plug.operator.id)
            wire_xml_node.setAttribute('fromAttr', wire.input_plug.name)
            wire_xml_node.setAttribute('toOperator', '%d' % wire.output_plug.operator.id)
            wire_xml_node.setAttribute('toAttr', wire.output_plug.name)
            wires_xml_node.appendChild(wire_xml_node)
            
        print('save_dataflow : saving to %s\n' % file_path)
        with open(file_path, 'w') as f:
            f.write(xml_doc.toprettyxml())
        
    def load_dataflow(self, file_path, dataflow):
        """
        :param msspecgui.dataflow.DataFlow dataflow: an empty dataflow that will be filled

        """
        xml_doc = minidom.parse(file_path)
        root_xml_node = xml_doc.documentElement

        last_create_op_id_xml_node = root_xml_node.getElementsByTagName('lastCreatedOperatorId')[0]
        last_created_operator_id = int(last_create_op_id_xml_node.getAttribute('value'))
        dataflow.last_created_operator_id = last_created_operator_id
        
        operators_xml_node = root_xml_node.getElementsByTagName('operators')[0]
        for op_xml_node in operators_xml_node.getElementsByTagName('operator'):
            # print('load_dataflow : creating operator')
            self._create_operator_from_xml(op_xml_node, dataflow)
            
        wires_xml_node = root_xml_node.getElementsByTagName('wires')[0]
        for wire_xml_node in wires_xml_node.getElementsByTagName('wire'):
            from_op_id = int(wire_xml_node.getAttribute('fromOperator'))
            from_attr = wire_xml_node.getAttribute('fromAttr')
            to_op_id = int(wire_xml_node.getAttribute('toOperator'))
            to_attr = wire_xml_node.getAttribute('toAttr')
            
            from_plug = dataflow.get_operator(from_op_id).get_plug(from_attr)
            to_plug = dataflow.get_operator(to_op_id).get_plug(to_attr)
            
            # print('load_dataflow : creating wire')
            dataflow.create_wire(from_plug, to_plug)
            
        return dataflow
