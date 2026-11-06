#!/usr/bin/python
# -*- coding: utf-8 -*-
import wx
from wx.lib.agw import floatspin as fs
"""
bulk interface that the user use to create a single cell crystal
"""


class StructureDef(object):
    # enabled widgets
    class EnabledWidgets(object):
        covera_flg =  0b00001
        u_flg =       0b00010
        c_flg =       0b00100
        orth_flg =    0b01000
        cubic_flg =   0b10000
    
    def __init__(self, ase_name, long_name, min_atoms, max_atoms, default_atoms, enabled_widgets, default_a, default_c=0.1):
        """
        :param ase_name: crystal structure identifier used by ase
        :param long_name: crystal long name
        :param min_atoms: minimal number of atoms allowed in crystal
        :param max_atoms: minimal number of atoms allowed in crystal
        """
        self.long_name = long_name
        self.ase_name = ase_name
        self.min_atoms = min_atoms
        self.max_atoms = max_atoms
        self.default_atoms = default_atoms
        self.enabled_widgets = enabled_widgets
        self.default_a = default_a
        self.default_c = default_c


class StructureDefs(object):
    
    def __init__(self):
        ew = StructureDef.EnabledWidgets
        self._long_name_to_structure_def = {}
        self._ase_name_to_structure_def = {}
        self._add_structure_def(StructureDef('sc',             'Simple Cubic',         -1,  1, 'Po',   ew.cubic_flg | ew.orth_flg,                         3.80))
        self._add_structure_def(StructureDef('fcc',            'Face Centered Cubic',  -1,  1, 'Cu',   ew.cubic_flg | ew.orth_flg,                         3.615))
        self._add_structure_def(StructureDef('bcc',            'Body Centered Cubic',  -1,  1, 'Fe',   ew.cubic_flg | ew.orth_flg,                         2.866))
        self._add_structure_def(StructureDef('hcp',            'Hexagonal Close-Pack',  2, -1, 'Mg',   ew.cubic_flg | ew.c_flg | ew.covera_flg,            3.209, 5.210))
        self._add_structure_def(StructureDef('diamond',        'Diamond',               2, -1, 'C',    ew.cubic_flg | ew.orth_flg,                         3.5668))
        self._add_structure_def(StructureDef('zincblende',     'Zincblende',            2, -1, 'ZnS',  ew.cubic_flg | ew.orth_flg,                         5.420))
        self._add_structure_def(StructureDef('rocksalt',       'Rocksalt',              2, -1, 'NaCl', ew.cubic_flg | ew.orth_flg,                         5.630))
        self._add_structure_def(StructureDef('cesiumchloride', 'Cesiumchloride',       -1,  2, 'CsCl', ew.cubic_flg,                                       4.110))
        self._add_structure_def(StructureDef('fluorite',       'Fluorite',             -1,  1, 'CaFF', ew.orth_flg | ew.c_flg | ew.covera_flg,             5.463))
        self._add_structure_def(StructureDef('wurtzite',       'Wurtzite',             -1,  1, 'ZnS',  ew.cubic_flg | ew.c_flg | ew.u_flg | ew.covera_flg, 3.820, 6.260))
        
    def _add_structure_def(self, structure_def):
        """
        :type structure_def: StructureDef
        """
        self._long_name_to_structure_def[structure_def.long_name] = structure_def
        self._ase_name_to_structure_def[structure_def.ase_name] = structure_def

    def get_structure_def(self, long_name=None, ase_name=None):
        """
        :rtype: StructureDef
        """
        assert (long_name is None) ^ (ase_name is None)
        if long_name is not None:
            return self._long_name_to_structure_def[long_name]
        if ase_name is not None:
            return self._ase_name_to_structure_def[ase_name]


class BulkFrame(wx.Dialog):
    """
    :param structure_defs: list of crystal structures
    :type structure_defs: dictionary
    :param selected_crystal_str: the structure that will be
    used to create the cluster
    """

    structure_defs = StructureDefs()

    @classmethod
    def get_structure_defs(cls):
        """
        """
        return cls.structure_defs

    @classmethod
    def get_structure_id(cls, structure_name):
        """
        :type structure_name: str
        :rtype: str
        """
        return cls.structure_defs.get_structure_def(long_name=structure_name).ase_name
        
    @classmethod
    def get_structure_name(cls, structure_id):
        """
        :type structure_id: str
        :rtype: str
        """
        return cls.structure_defs.get_structure_def(ase_name=structure_id).long_name

    def get_crystal_def(self, long_name):
        """
        :param long_name : the structure long name (eg 'Face Centered Cubic'
        :type  long_name : str
        :rtype: StructureDef
        """

        return BulkFrame.get_structure_defs().get_structure_def(long_name=long_name)

    selected_crystal_str = ""

    def __init__(self, parent, title, bulk_operator):
        """
        :param parent: window that contains this window
        :type parent: wx.Widget
        :param bulk_operator: the bulk_operator that this gui allows the user to modify
        :type bulk_operator: msspec.cluster.generators.Bulk
        """
        self._bulk_operator = bulk_operator
        super(BulkFrame, self).__init__(parent, title=title)
        self.initui()
        self.Centre()
        self.Show()

    # --------------------------------------------------------------------------

    def initui(self):
        """
        the constructor of the interface
        """
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.panel = wx.Panel(parent=self, id=1)
        main_params_sizer = wx.GridBagSizer()

        # all crystal structures that can be selected to create the cluster
        self.crystal_structures = ['Simple Cubic',
                                   'Face Centered Cubic',
                                   'Body Centered Cubic',
                                   'Hexagonal Close-Pack',
                                   'Diamond',
                                   'Zincblende',
                                   'Rocksalt',
                                   'Cesiumchloride',
                                   'Wurtzite']

        # crystal structure
        crystalstructure_txt = wx.StaticText(self.panel, label="Crystal structure : ")
        main_params_sizer.Add(crystalstructure_txt, pos=(1, 0), flag=wx.LEFT | wx.TOP, border=10)
        # list of crystal structures
        self.crystr_cbx = wx.ComboBox(self.panel, choices=self.crystal_structures)
        main_params_sizer.Add(self.crystr_cbx, pos=(1, 1), border=10)
        crystal_structure_name = BulkFrame.get_structure_name(self._bulk_operator.structure)
        self.crystr_cbx.Select(self.crystal_structures.index(crystal_structure_name))
        # a ToolTip associated to the combobox
        self.crystr_cbx.SetToolTip(wx.ToolTip("choose a crystal structure"))
        self.crystr_cbx.Bind(wx.EVT_COMBOBOX, self.on_crystal_select)

        # the atoms that constitute the crystal
        self.atoms_txt = wx.StaticText(self.panel, label="Atoms : ")
        main_params_sizer.Add(self.atoms_txt, pos=(3, 0), flag=wx.LEFT, border=10)
        # to insert the name of the cluster
        self.atoms_txc = wx.TextCtrl(self.panel)
        main_params_sizer.Add(self.atoms_txc, pos=(3, 1), flag=wx.EXPAND, border=10)
        # default value is 'MgO'
        self.atoms_txc.SetLabel(self._bulk_operator.atoms)
        self.atoms_txc.SetToolTip(wx.ToolTip("the atoms that form the crystal (eg MgO or NaCl)"))
        self.atoms_txc.Bind(wx.EVT_TEXT, self.alert_crystal_name)
        self.numatoms_txt = wx.StaticText(self.panel, label="this crystal is composed of two atoms")
        main_params_sizer.Add(self.numatoms_txt, pos=(4, 1), flag=wx.LEFT, border=10)

        # is cubic
        self.is_cubic_txt = wx.StaticText(self.panel, label="Cubic : ")
        main_params_sizer.Add(self.is_cubic_txt, pos=(5, 0), flag=wx.TOP | wx.LEFT, border=10)
        self.is_cubic_chk = wx.CheckBox(self.panel)
        main_params_sizer.Add(self.is_cubic_chk, pos=(5, 1), flag=wx.TOP | wx.EXPAND, border=10)
        self.is_cubic_chk.SetValue(self._bulk_operator.is_cubic)

        # is othorhombic
        self.is_orth_txt = wx.StaticText(self.panel, label="Orthorhombic : ")
        main_params_sizer.Add(self.is_orth_txt, pos=(7, 0), flag=wx.TOP | wx.LEFT, border=10)
        self.is_orth_chk = wx.CheckBox(self.panel)
        main_params_sizer.Add(self.is_orth_chk, pos=(7, 1), flag=wx.TOP | wx.EXPAND, border=10)
        self.is_orth_chk.SetValue(self._bulk_operator.is_orthorombic)

        # lattice parameters
        lattice_params_sizer = wx.BoxSizer(wx.VERTICAL)
        self.lattice_params_stbox = wx.StaticBox(self.panel, label='Lattice Parameters')
        lattice_params_subsizer = wx.StaticBoxSizer(self.lattice_params_stbox, wx.VERTICAL)
        
        lattice_params_subsizer.AddSpacer(3)

        # lattice parameter a
        a_sizer = wx.BoxSizer(wx.HORIZONTAL)
        a_txt = wx.StaticText(self.panel, label=" a  = ")
        self.a_fspin = fs.FloatSpin(self.panel,
                                    value=self._bulk_operator.a, size=(180, -1),
                                    min_val=0.1, max_val=10.00,
                                    increment=0.01,
                                    style=fs.FS_CENTRE)
        self.a_fspin.SetFormat("%f")
        self.a_fspin.SetDigits(6)
        self.a_fspin.Bind(fs.EVT_FLOATSPIN, self.alert_a)
        self.a_fspin.Bind(wx.EVT_SPINCTRL, self.on_a_changed)
        a_sizer.Add(a_txt)
        a_sizer.Add(self.a_fspin)

        lattice_params_subsizer.AddSpacer(10)

        # lattice parameter c
        c_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.c_txt = wx.StaticText(self.panel, 6, label=" c  = ")
        self.c_txt.Disable()
        self.c_fspin = fs.FloatSpin(self.panel,
                                    value=self._bulk_operator.c, size=(160, -1),
                                    min_val=0.1, max_val=10.00, increment=0.01)
        self.c_fspin.SetFormat("%f")
        self.c_fspin.SetDigits(6)
        self.c_fspin.Bind(wx.EVT_SPINCTRL, self.on_c_changed)
        self.c_fspin.Disable()
        c_sizer.Add(self.c_txt)
        c_sizer.Add(self.c_fspin)

        # lattice parameter u
        u_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.u_txt = wx.StaticText(self.panel, 7, label=" u  = ")
        self.u_txt.Disable()
        self.u_fspin = fs.FloatSpin(self.panel,
                                    value=self._bulk_operator.u, size=(160, -1),
                                    min_val=0.1, max_val=10.00, increment=0.01)
        self.u_fspin.SetFormat("%f")
        self.u_fspin.SetDigits(6)
        self.u_fspin.Disable()
        u_sizer.Add(self.u_txt)
        u_sizer.Add(self.u_fspin)

        # lattice parameter c/a
        covera_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.covera_txt = wx.StaticText(self.panel, label="c/a ratio : ")
        self.covera_txt.Disable()
        self.covera_fspin = fs.FloatSpin(self.panel,
                                         value='0.0', size=(160, -1),
                                         min_val=0.1, max_val=10.00, increment=0.01)
        self.covera_fspin.SetFormat("%f")
        self.covera_fspin.SetDigits(6)
        self.covera_fspin.Bind(wx.EVT_SPINCTRL, self.on_covera_changed)
        self.covera_fspin.Disable()
        covera_sizer.Add(self.covera_txt)
        covera_sizer.Add(self.covera_fspin)

        lattice_params_subsizer.Add(a_sizer)
        lattice_params_subsizer.AddSpacer(10)
        lattice_params_subsizer.Add(c_sizer)
        lattice_params_subsizer.AddSpacer(10)
        lattice_params_subsizer.Add(u_sizer)
        lattice_params_subsizer.AddSpacer(10)
        lattice_params_subsizer.Add(covera_sizer)
        lattice_params_subsizer.AddSpacer(10)

        lattice_params_sizer.Add(lattice_params_subsizer, flag=wx.ALL | wx.EXPAND, border=5)

        ok_cancel_reset_sizer = wx.BoxSizer(wx.HORIZONTAL)
        reset_btn = wx.Button(self.panel, 2, label="Reset")
        reset_btn.Bind(wx.EVT_BUTTON, self.on_reset_crystal_params)
        cancel_btn = wx.Button(self.panel, wx.ID_CANCEL, label="Cancel")
        cancel_btn.Bind(wx.EVT_BUTTON, self.on_close, id=wx.ID_CANCEL)
        self.ok_btn = wx.Button(self.panel, wx.ID_OK, label="OK")
        self.ok_btn.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        ok_cancel_reset_sizer.Add(reset_btn)
        ok_cancel_reset_sizer.Add(cancel_btn)
        ok_cancel_reset_sizer.Add(self.ok_btn)

        self.vbox.Add(main_params_sizer)
        self.vbox.Add(lattice_params_sizer, proportion=1, flag=wx.LEFT | wx.EXPAND, border=10)
        self.vbox.Add(ok_cancel_reset_sizer, flag=wx.ALIGN_RIGHT, border=5)
        self.panel.SetSizer(self.vbox)
        self.panel.Fit()
        self.Fit()

    def on_reset_crystal_params(self, event):
        self.on_crystal_select(event)

    def on_close(self, event):
        event.Skip()  # propagate the event so that the dialog closes

    def on_ok(self, event):
        a_val = str(self.a_fspin.GetValue())
        c_val = str(self.c_fspin.GetValue())
        u_val = str(self.u_fspin.GetValue())
        covera_val = str(self.covera_fspin.GetValue())

        print('bulk =  (%s,' % self.atoms_txc.GetValue() +
              ' %s,' % self.crystr_cbx.GetValue() + ' a= %s,' % a_val +
              ' c= %s,' % c_val +
              ' u= %s,' % u_val +
              ' c/a ratio= %s,' % covera_val +
              ' cubic= %s,' % str(self.is_cubic_chk.GetValue()) +
              ' orthorombic= %s)' % str(self.is_orth_chk.GetValue()))
        
        self._bulk_operator.atoms = self.atoms_txc.GetValue().encode('ascii', 'ignore')
        self._bulk_operator.structure = BulkFrame.get_structure_id(self.crystr_cbx.GetValue())
        self._bulk_operator.a = float(a_val)
        self._bulk_operator.c = float(c_val)
        self._bulk_operator.u = float(u_val)
        self._bulk_operator.is_cubic = self.is_cubic_chk.GetValue()
        self._bulk_operator.is_orthorombic = self.is_orth_chk.GetValue()
        event.Skip()  # propagate the event so that the dialog closes

    def alert_a(self, event):
        if self.a_fspin.GetValue() == 0.00:
            wx.MessageBox("Please enter the parameter a",
                          style=wx.OK | wx.ICON_MASK)

    def alert_crystal_name(self, event):
        if self.atoms_txc.GetValue() == '':
            wx.MessageBox("Please enter the Atoms that\ncomposed the crystal ",
                          style=wx.CANCEL | wx.ICON_MASK)
            
    def on_a_changed(self, event):
        """
        change values of c and covera when the user insert the value of the parameter a
        """
        self.covera_fspin.SetValue(self.c_fspin.GetValue() / self.a_fspin.GetValue())

    def on_c_changed(self, event):
        """
        change values of a and covera when the user insert the value of the parameter c
        """
        self.covera_fspin.SetValue(self.c_fspin.GetValue() / self.a_fspin.GetValue())

    """
    def change_covera_val(self):
        print ("---_cov_---")
        # if self.c_fspin.GetValue() != 0.00000:
        self.covera_fspin.SetValue(
            self.c_fspin.GetValue() / self.a_fspin.GetValue())
        print ("---***_cov_***---")
    """
    
    def on_covera_changed(self, event):
        """
            change values of c and a when the user
            insert the value of the parameter covera
        """
        self.change_covera_val()
        
    def change_covera_val(self):
        self.c_fspin.SetValue(self.covera_fspin.GetValue() * self.a_fspin.GetValue())

    def get_atoms_constraint_msg(self, num_atoms_min, num_atoms_max):
        """
        return the constraint of atoms number
        """
        message = ''
        if num_atoms_min != -1:
            message += 'minimum %d atoms' % num_atoms_min
        if num_atoms_max != -1:
            if len(message) > 0:
                message += ', '
            message += 'maximum %d atoms' % num_atoms_max
        return message
    
    def update_widgets(self, enabled_widgets):
        """
        manage the widgets
        """
        ew = StructureDef.EnabledWidgets

        if enabled_widgets & ew.cubic_flg:
            self.is_cubic_txt.Enable()
            self.is_cubic_chk.Enable()
            self.is_cubic_chk.SetValue(False)
        else:
            self.is_cubic_chk.SetValue(False)
            self.is_cubic_txt.Disable()
            self.is_cubic_chk.Disable()

        if enabled_widgets & ew.orth_flg:
            self.is_orth_txt.Enable()
            self.is_orth_chk.Enable()
            self.is_orth_chk.SetValue(False)
        else:
            self.is_orth_chk.SetValue(False)
            self.is_orth_txt.Disable()
            self.is_orth_chk.Disable()

        if enabled_widgets & ew.c_flg:
            self.c_txt.Enable()
            self.c_fspin.Enable()
        else:
            self.c_txt.Disable()
            self.c_fspin.Disable()

        if enabled_widgets & ew.u_flg:
            self.u_txt.Enable()
            self.u_fspin.Enable()
        else:
            self.u_txt.Disable()
            self.u_fspin.Disable()

        if enabled_widgets & ew.covera_flg:
            self.covera_txt.Enable()
            self.covera_fspin.Enable()
        else:
            self.covera_txt.Disable()
            self.covera_fspin.Disable()

    def on_crystal_select(self, event):
        """
        when the user choose one structure this function will be triggered
        """

        self.selected_crystal_str = self.crystr_cbx.GetValue()
        # print self.selected_crystal_str

        crystal_def = self.get_crystal_def(self.selected_crystal_str)
        # print ('values : %s' % crystal_def)
        self.numatoms_txt.SetLabel(str(self.get_atoms_constraint_msg(crystal_def.min_atoms, crystal_def.max_atoms)))
        self.update_widgets(crystal_def.enabled_widgets)
        self.atoms_txc.SetLabel(crystal_def.default_atoms)
        self.a_fspin.SetValue(crystal_def.default_a)
        self.c_fspin.SetValue(crystal_def.default_c)
        self.covera_fspin.SetValue(crystal_def.default_c / crystal_def.default_a)

        self.a_fspin.Bind(fs.EVT_FLOATSPIN, self.alert_a)
        self.lattice_params_stbox.Refresh()
        self.vbox.RecalcSizes()
        self.panel.SetSizer(self.vbox)


if __name__ == '__main__':
    app = wx.App()
    BulkFrame(None, title="create a single cell", bulk_operator=None)
    app.MainLoop()
