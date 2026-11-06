from msspecgui import dataflow
from ase.build import add_adsorbate


class AddAdsorbate(dataflow.Operator):
    """
    an operator for ase.build.add_absorbate
    """

    class Creator(dataflow.IOperatorCreator):
        def __init__(self):
            dataflow.IOperatorCreator.__init__(self)
            self.add_input_attribute('slab', 'physics.atomscluster', 'the surface onto which the adsorbate should be added')
            self.add_input_attribute('adsorbate', 'physics.atomscluster', 'the adsorbate')
            self.add_input_attribute('height', 'float', 'height above the surface.', is_pluggable=False)

            self.add_output_attribute('output_cluster', 'physics.atomscluster', 'the resulting cluster')
            
        def get_operator_type_id(self):
            return 'ase.build.add_adsorbate'
        
        def create_operator(self, dflow):
            return AddAdsorbate(dflow, self)

    def __init__(self, data_flow, creator):
        '''
        Constructor
        
        :param data_flow: the dataflow that will contain this operator
        '''
        dataflow.Operator.__init__(self, data_flow, creator)
        self.get_plug('height').set_value(0.0)

    def update(self):
        """
        see dataflow.Operator.update
        """
        slab = self.get_plug('slab').get_value()
        adsorbate = self.get_plug('adsorbate').get_value()
        height = self.get_plug('height').get_value()
        output_cluster = slab.copy()
        add_adsorbate(output_cluster, adsorbate, height)
        self.get_plug('output_cluster').set_value(output_cluster)
