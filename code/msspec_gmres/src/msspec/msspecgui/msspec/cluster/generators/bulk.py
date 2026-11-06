from ase.lattice import bulk

from msspecgui import dataflow


class Bulk(dataflow.Operator):
    '''
    a cluster creator that performs what ase.lattice.bulk performs
    '''

    class Creator(dataflow.IOperatorCreator):
        def __init__(self):
            dataflow.IOperatorCreator.__init__(self)
            self.add_output_attribute('output_cluster', 'physics.atomscluster', 'the resulting cluster')
            
            # the following attributes are set via a graphical user interface, hence not pluggable
            self.add_input_attribute('atoms', 'string', 'the atoms that constitute this cluster (eg "MgO")', is_pluggable=False)
            self.add_input_attribute('structure', 'string', 'the lattice structure (eg "rocksalt")', is_pluggable=False)
            self.add_input_attribute('a', 'float', 'the lattice constant a', is_pluggable=False)
            self.add_input_attribute('c', 'float', 'the lattice constant c', is_pluggable=False)
            self.add_input_attribute('u', 'float', 'the lattice constant u', is_pluggable=False)
            self.add_input_attribute('is_cubic', 'bool', 'is the lattice cubic', is_pluggable=False)
            self.add_input_attribute('is_orthorombic', 'bool', 'is the lattice orthorhombic', is_pluggable=False)
        
        def get_operator_type_id(self):
            return 'ase.lattice.bulk'
        
        def create_operator(self, dflow):
            return Bulk(dflow, self)

    def __init__(self, data_flow, creator):
        '''
        Constructor
        
        :param data_flow: the dataflow that will contain this operator
        '''
        dataflow.Operator.__init__(self, data_flow, creator)
        self.atoms = 'MgO'
        self.structure = 'rocksalt'
        self.a = 4.219
        self.c = 0.0
        self.u = 0.0
        self.is_cubic = True
        self.is_orthorombic = False

    @property
    def atoms(self):
        return self.get_plug('atoms').get_value()

    @atoms.setter
    def atoms(self, atoms):
        """
        :type atoms: str
        """
        self.get_plug('atoms').set_value(atoms)
    
    @property
    def structure(self):
        return self.get_plug('structure').get_value()

    @structure.setter
    def structure(self, structure):
        """
        :type structure: str
        """
        self.get_plug('structure').set_value(structure)
        
    @property
    def a(self):
        return self.get_plug('a').get_value()

    @a.setter
    def a(self, a):
        """
        :type a: float
        """
        self.get_plug('a').set_value(a)

    @property
    def c(self):
        return self.get_plug('c').get_value()

    @c.setter
    def c(self, c):
        """
        :type c: float
        """
        self.get_plug('c').set_value(c)

    @property
    def u(self):
        return self.get_plug('u').get_value()

    @u.setter
    def u(self, u):
        """
        :type u: float
        """
        self.get_plug('u').set_value(u)

    @property
    def is_cubic(self):
        return self.get_plug('is_cubic').get_value()

    @is_cubic.setter
    def is_cubic(self, is_cubic):
        """
        :type is_cubic: bool
        """
        self.get_plug('is_cubic').set_value(is_cubic)

    @property
    def is_orthorombic(self):
        return self.get_plug('is_orthorombic').get_value()

    @is_orthorombic.setter
    def is_orthorombic(self, is_orthorombic):
        """
        :type is_orthorombic: bool
        """
        self.get_plug('is_orthorombic').set_value(is_orthorombic)
    
    def update(self):
        """
        see dataflow.Operator.update
        """
        # >>> ase.lattice.bulk('MgO', crystalstructure='rocksalt', a=4.219, c=0.1, u=0.1, cubic=True, orthorhombic=False)
        # Atoms(symbols='MgOMgOMgOMgO', positions=..., cell=[4.219, 4.219, 4.219], pbc=[True, True, True])
        print('self.atoms = %s (type = %s)' % (self.atoms, type(self.atoms)))
        cluster = bulk(self.atoms, crystalstructure=self.structure, a=self.a, c=self.c, u=self.u, cubic=self.is_cubic, orthorhombic=self.is_orthorombic)
        self.get_plug('output_cluster').set_value(cluster)
