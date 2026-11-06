import abc


class ILayoutManager(object):
    
    @abc.abstractmethod
    def compute_operators_position(self, data_flow):
        """
        :param DataFlow data_flow: the data flow that we want to compute the layout for
        :return dict(msspecgui.dataflow.Operator, (int, int)): a (x, y) position for each operator
        """
        return {}
