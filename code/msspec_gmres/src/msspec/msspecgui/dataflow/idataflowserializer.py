import abc


class IDataflowSerializer(object):
    
    @abc.abstractmethod
    def save_dataflow(self, dataflow, filepath):
        """ saves the given dataflow at the given file location
        
        :param dataflow: the dataflow to save
        :type dataflow: msspecgui.dataflow.DataFlow
        :param file_path: the path of the file that will store this dataflow in serialized form
        :type file_path: str
        """
        pass

    @abc.abstractmethod
    def load_dataflow(self, file_path, dataflow):
        """loads the dataflow from the given serialized dataflow file

        :param str file_path: the xml file describng the dataflow to load
        :param msspecgui.dataflow.DataFlow dataflow: an empty dataflow that will be filled

        :rtype: DataFlow
        """
        return None
