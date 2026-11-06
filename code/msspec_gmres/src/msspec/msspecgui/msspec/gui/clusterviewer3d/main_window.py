# coding: utf-8

import wx
from GLPanel import GLPanel
from msspec.utils import hemispherical_cluster
from ase.build import bulk


class GLFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)
        self.pnl = GLPanel(self)

    def set_atoms(self, atoms):
        self.pnl.set_atoms(atoms)


if __name__ == "__main__":

    lattice = bulk('MgO', crystalstructure='rocksalt', a=4.21, cubic=True)
    cluster = hemispherical_cluster(lattice, 0, 0, diameter=20.)


    app = wx.App(False)
    frame = GLFrame(None, -1, 'GL Window', size=(640, 480))
    frame.set_atoms(cluster)
    frame.Show()

    app.MainLoop()
    app.Destroy()