from msspecgui import dataflow


class Repeat(dataflow.Operator):
    """
    an operator that repeats the input cluser a given number of times along each axis
    """
    # an operator that removes the atoms of the cluster that are outside the given sphere

    class Creator(dataflow.IOperatorCreator):
        def __init__(self):
            dataflow.IOperatorCreator.__init__(self)
            self.add_input_attribute('input_cluster', 'physics.atomscluster', 'the cluster from which we want to keep the atoms that are inside the given sphere')
            self.add_input_attribute('repeat_x', 'int', 'the number of repetitions along x axis', is_pluggable=False)
            self.add_input_attribute('repeat_y', 'int', 'the number of repetitions along y axis', is_pluggable=False)
            self.add_input_attribute('repeat_z', 'int', 'the number of repetitions along z axis', is_pluggable=False)

            # self.add_input_attribute('inputSphere.center', 'position3', 'the position of the center of the sphere')
            # self.add_input_attribute('inputSphere.radius', 'float', 'the radius of the center of the sphere')
            self.add_output_attribute('output_cluster', 'physics.atomscluster', 'the resulting cluster')
            
        def get_operator_type_id(self):
            return 'ase.repeat'
        
        def create_operator(self, dflow):
            return Repeat(dflow, self)

    def __init__(self, data_flow, creator):
        '''
        Constructor
        
        :param data_flow: the dataflow that will contain this operator
        '''
        dataflow.Operator.__init__(self, data_flow, creator)
        self.get_plug('repeat_x').set_value(2)
        self.get_plug('repeat_y').set_value(3)
        self.get_plug('repeat_z').set_value(4)

    def update(self):
        """
        see dataflow.Operator.update
        """
        input_cluster = self.get_plug('input_cluster').get_value()
        rx = self.get_plug('repeat_x').get_value()
        ry = self.get_plug('repeat_y').get_value()
        rz = self.get_plug('repeat_z').get_value()
        output_cluster = input_cluster.repeat((rx, ry, rz))
        self.get_plug('output_cluster').set_value(output_cluster)
