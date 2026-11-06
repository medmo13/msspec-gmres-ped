# -*- encoding: utf-8 -*-
# vim: set fdm=indent ts=2 sw=2 sts=2 et tw=80 cc=+1 mouse=a nu : #

# import sys
# sys.path.append('/home/sylvain/dev/msspec/trunk/src/python/MsSpecGui/')

import wx
from clusterviewer import ClusterViewer

from ase.lattice import bulk


class MyFrame(wx.Frame):
    def __init__(self, **kwargs):
        p = {}
        p.setdefault('parent', None)
        p.setdefault('title', 'Test Frame')
        p.setdefault('size', (640, 480))
        p.update(kwargs)

        wx.Frame.__init__(self, **p)

        self.Window = ClusterViewer(self)
        self.Show()

# display an MgO cell
MgO = bulk('MgO', crystalstructure='rocksalt', orthorhombic=True, a=4.21)
MgO = MgO.repeat((20, 20, 30))

app = wx.App(False)
frame = MyFrame(title='MsSpec Viewer')
frame.Window.light_mode_threshold = 2
frame.Window.set_atoms(MgO, rescale=True)

frame.Show()

app.MainLoop()
app.Destroy()
