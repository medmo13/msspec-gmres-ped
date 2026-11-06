'''
Created on Dec 23, 2015

@author: graffy

how to install ase on osx macports
    [OK]graffy@pr079234:~/ownCloud/ipr/msspec/toto[12:50:42]>export ASE_TAGS=https://svn.fysik.dtu.dk/projects/ase/tags/
    [OK]graffy@pr079234:~/ownCloud/ipr/msspec/toto[12:50:57]>svn co -r 4567 $ASE_TAGS/3.9.1 ase-3.9.1
    [OK]graffy@pr079234:~[12:52:41]>cd ownCloud/ipr/msspec/toto/ase-3.9.1
    sudo /opt/local/bin/python setup.py install
'''
import sys
# import re
# sudo port install py27-wxpython-3.0
# import os

sys.path.append('../../')
import wx  # @IgnorePep8
import ase  # @IgnorePep8
# import msspec
import msspecgui.msspec.gui.mainwindow as mainwindow  # @IgnorePep8
# import cairosvg


def view_atoms(atoms):
    '''Displays the atoms.
    
    This is the equivalent of ase.visualize.view(molecule), but this function doesn't make use of the external python program ase-3.9.1/tools/ase-gui. It was built by looking at what ase.visualize.view does.
     
    '''
    import ase.gui.gui  # IGNORE:redefined-outer-name
    # import ase.gui.ag
    import ase.io
    # import ase.gui.gtkexcepthook
    # from ase.visualize import view
    # from ase.io import write
    # from ase.gui.ag import main
    # http://wiki.wxpython.org/Getting%20Started
    
    import ase.gui.images
    file_format = 'traj'
    file_path = '/tmp/toto.%s' % file_format
    # ase.io.write(file_path, molecule, file_format)
    ase.io.write(file_path, images=atoms, format=file_format)

    images = ase.gui.images.Images()
    images.read([file_path], ase.io.string2index(':'))
    gui = ase.gui.gui.GUI(images, '', 1, False)
    gui.run(None)
    

def test_ase():
    d = 1.10
    molecule = ase.Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])
    view_atoms(molecule)


def main():
    print('hello, world')
    app = wx.App(False)  # Create a new app, don't redirect stdout/stderr to a window.
    app.SetAppName('MsSpec')
    # app.SetAppDisplayName('MsSpec')
    frame = mainwindow.MainWindow(None, 'MsSpec')  # A Frame is a top-level window. @UnusedVariable

    # test_ase()
    app.MainLoop()
    

if __name__ == '__main__':
    main()
