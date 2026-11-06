#!/usr/bin/env python
# coding: utf-8
#
# Copyright © 2016-2020 - Rennes Physics Institute
#
# This file is part of msspec.
#
# msspec is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# msspec is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this msspec.  If not, see <http://www.gnu.org/licenses/>.
#
# Source file  : src/msspec/iodata_wx.py
# Last modified: Mon, 27 Sep 2021 17:49:48 +0200
# Committed by : sylvain tricot <sylvain.tricot@univ-rennes1.fr>


"""
Module iodata
=============

This module contains all classes useful to manipulate, store and display
data results.

The :py:class:`Data` and :py:class:`DataSet` are the two enduser classes
important to manipulate the data.
Here is an example of how to store values in a Data object:

.. code-block:: python

    from msspec.iodata import Data
    import numpy as np


    # Let's create first some dumy data
    X = np.arange(0, 20)
    Y = X**2

    # Create a Data object. You need to give a title as an argument
    data = Data('all my data')
    # and append a new DataSet with its title
    dset = data.add_dset('Dataset 0')

    # To feed the DataSet with columns, use the add_columns method
    # and provide as many keywords as you like. Each key being the
    # column name and each value being an array holding the column
    # data.
    dset.add_columns(x=X, y=Y, z=X+2, w=Y**3)
    # you can provide parameters with their values with keywords as well
    dset.add_parameter(name='truc', group='main', value='3.14', unit='eV')

    # To plot these data, you need to add a 'view' with its title
    view = dset.add_view('my view')
    # You then need to select which columns you which to plot and
    # and under wich conditions (with the 'which' keyword)
    view.select('x', 'y', where="z<10", legend=r"z = 0")
    view.select('x', 'y', where="z>10", legend=r"z = 1")

    # To pop up the graphical window
    data.view()

"""


import os
import sys
from distutils.version import LooseVersion
from distutils.version import StrictVersion
from io import StringIO
from datetime import datetime

import ase.io
import h5py
import numpy as np
import wx.grid
from lxml import etree
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
#from matplotlib.backends.backend_wxagg import FigureCanvasWx as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
from matplotlib.figure import Figure
from terminaltables import AsciiTable

import msspec
from msspec.msspecgui.msspec.gui.clusterviewer import ClusterViewer
from msspec.misc import LOGGER


def cols2matrix(x, y, z, nx=88*1+1, ny=360*1+1):
    # mix the values of existing theta and new theta and return the
    # unique values
    newx = np.linspace(np.min(x), np.max(x), nx)
    newy = np.linspace(np.min(y), np.max(y), ny)
    ux = np.unique(np.append(x, newx))
    uy = np.unique(np.append(y, newy))

    # create an empty matrix to hold the results
    zz = np.empty((len(ux), len(uy)))
    zz[:] = np.nan

    for p in zip(x, y, z):
        i = np.argwhere(ux == p[0])
        j = np.argwhere(uy == p[1])
        zz[i, j] = p[2]

    for i in range(len(ux)):
        #ok, = np.where(-np.isnan(zz[i,:]))
        ok, = np.where(~np.isnan(zz[i, :]))
        if len(ok) > 0:
            xp = uy[ok]
            fp = zz[i, ok]
            zz[i,:] = np.interp(uy, xp, fp)

    for i in range(len(uy)):
        #ok, = np.where(-np.isnan(zz[:,i]))
        ok, = np.where(~np.isnan(zz[:, i]))
        if len(ok) > 0:
            xp = ux[ok]
            fp = zz[ok, i]
            zz[:,i] = np.interp(ux, xp, fp)

    return ux, uy, zz


def is_title_valid(title):
    """ Ensure the string does not contain special characters:
        /\:*?"<>|
    """
    special_chars = ('/', '\\', ':', '*', '?', '\"', '<', '>', '|')
    for char in special_chars:
        if title.find(char) > -1:
            return False
    return True



class _DataPoint(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def __getattr__(self, name):
        if name in list(self.keys()):
            return self[name]
        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(
                self.__class__.__name__, name))

class DataSet(object):
    """
    This class can create an object to hold column-oriented data.

    :param title: The text used to entitled the dataset
    :type title: str
    :param notes: Some comments to add to the data
    :type notes: str

    """
    def __init__(self, title, notes=""):
        assert is_title_valid(title), '\/:*?"<>| are not allowed in the string'
        #self._col_names = []
        #self._col_arrays = []
        self.__dict__['_col_names'] = []
        self.__dict__['_col_arrays'] = []

        self.title = title
        self.notes = notes
        self._views = []
        self._parameters = []
        self.attributes = {}


        self._defaults = {'bool': False, 'str': '', 'int': 0, 'float': 0.,
                          'complex': complex(0)}
        self._formats = {bool: '{:s}', str: '{:s}', int: '{:<20d}',
                         float: '{:<20.10e}', complex: 's'}
        self._formats = ((np.integer, '{:<20d}'),
                         (np.floating, '{:<20.10e}'),
                         (np.complex, '({0.real:<.10e} {0.imag:<.10e}j)'),
                         (np.bool, '{:s}'),
                         (str, '{:s}'))


    def _empty_array(self, val):
        if isinstance(val, str):
            t = 'S256'
        else:
            t = np.dtype(type(val))

        if isinstance(val, bool):
            default = self._defaults['bool']
        elif isinstance(val, str):
            default = self._defaults['str']
        elif isinstance(val, int):
            default = self._defaults['int']
        elif isinstance(val, float):
            default = self._defaults['float']
        elif isinstance(val, complex):
            default = self._defaults['complex']
        else:
            raise TypeError('Not a supported type')

        return np.array([default]*len(self), dtype=t)

    def add_row(self, **kwargs):
        """Add a row of data into the dataset.

        :param kwargs: Each keyword is a column name. The number of keywords (columns) must be coherent with the
            number of existing columns. If no column are defined yet, they will be created.

        """
        for k, v in list(kwargs.items()):
            if k not in self._col_names:
                self._col_names.append(k)
                self._col_arrays.append(self._empty_array(v))
        for k, v in list(kwargs.items()):
            i = self._col_names.index(k)
            arr = self._col_arrays[i]
            arr = np.append(arr, v)
            self._col_arrays[i] = arr

    def add_columns(self, **kwargs):
        """
        Add columns to the dataset.

        You can provide as many columns as you want to this function. This
        function can be called several times on the same dataset but each time
        with different column names. Column names are given as keywords.

        :Example:

        >>> from iodata import DataSet
        >>> dset = DataSet('My Dataset', notes="Just an example")
        >>> xdata = range(10)
        >>> ydata = [i**2 for i in xdata]
        >>> dset.add_columns(x=xdata, y=ydata)
        >>> print dset
        >>> +-------+
        >>> | x  y  |
        >>> +-------+
        >>> | 0  0  |
        >>> | 1  1  |
        >>> | 2  4  |
        >>> | 3  9  |
        >>> | 4  16 |
        >>> | 5  25 |
        >>> | 6  36 |
        >>> | 7  49 |
        >>> | 8  64 |
        >>> | 9  81 |
        >>> +-------+

        """
        for k, vv in list(kwargs.items()):
            assert k not in self._col_names, ("'{}' column already exists"
                                              "".format(k))
            #if len(self) > 0:
            #    assert len(vv) == len(self), (
            #        'Too many values in the column (max = {})'.format(
            #            len(self)))
        for k, vv in list(kwargs.items()):
            arr = np.array(vv)
            self._col_names.append(k)
            self._col_arrays.append(arr)

    def delete_rows(self, itemspec):
        """
        Delete the rows specified with itemspec.

        """
        for i in range(len(self._col_names)):
            self._col_arrays[i] = np.delete(self._col_arrays[i], itemspec)

    def delete_columns(self, *tags):
        """
        Removes all columns name passed as arguments

        :param tags: column names.
        :type tags: str

        """
        for tag in tags:
            i = self._col_names.index(tag)
            self._col_names.pop(i)
            self._col_arrays.pop(i)

    def columns(self):
        """
        Get all the column names.

        :return: List of column names.
        :rtype: List of str

        """
        return self._col_names

    def add_view(self, name, overwrite=False, **plotopts):
        """
        Creates a new view named *name* with specied plot options.

        :param name: name of the view.
        :type name: str
        :param plotopts: list of keywords for configuring the plots.
        :return: a view.
        :rtype: :py:class:`iodata._DataSetView`
        """
        if overwrite:
            self.delete_view(name)

        if isinstance(name, str):
            v = _DataSetView(self, name, **plotopts)
        else:
            v = name
            v.dataset = self
        self._views.append(v)
        return v

    def delete_view(self, name):
        view_titles = [_.title for _ in self._views]
        try:
            i = view_titles.index(name)
            self._views.pop(i)
        except:
            pass

    def views(self):
        """Returns all the defined views in the dataset.

        :return: A list of view
        :rtype: List of :py:class:`iodata._DataSetView`
        """
        return self._views

    def add_parameter(self, **kwargs):
        """Add a parameter to store with the dataset.

        :param kwargs: list of keywords with str values.

        These keywords are:
            * name: the name of the parameter.
            * group: the name of a group it belongs to.
            * value: the value of the parameter.
            * unit: the unit of the parameter.

        For example:

        .. code-block:: python

            from iodata import DataSet

            mydset = DataSet("Experiment")
            mydset.add_parameter(name='Spectrometer', group='misc', value='Omicron', unit='')

        """
        self._parameters.append(kwargs)

    def parameters(self):
        """
        Returns the list of defined parameters.

        :return: all parameters defined in the :py:class:`iodata.DataSet` object.
        :rtype: List of dict
        """
        return self._parameters

    def get_parameter(self, group=None, name=None):
        """Retrieves all parameters for a given name and group.

        * If *name* is given and *group* is None, returns all parameters with such a *name* in all groups.
        * If *group* is given and *name* is None, returns all parameters in such a *group*
        * If both *name* and *group* are None. Returns all parameters (equivalent to
            :py:func:`iodata.DataSet.parameters`).

        :param group: The group name or None.
        :type group: str
        :param name: The parameter's name or None.
        :type name: str
        :return: A list of parameters.
        :rtype: List of dict
        """
        p = []
        for _ in self._parameters:
            if _['group'] == group or group == None:
                if _['name'] == name or name == None:
                    p.append(_)
        return p[0] if len(p) == 1 else p

    def get_cluster(self):
        """Get all the atoms in the cluster.

        :return: The cluster
        :rtype: :py:class:`ase.Atoms`
        """
        s = StringIO()
        s.write(self.get_parameter(group='Cluster', name='cluster')['value'])
        return ase.io.read(s, format='xyz')


    def select(self, *args, **kwargs):
        condition = kwargs.get('where', 'True')
        indices = []

    def export_views(self, folder):
        for view in self.views():
            f = view.get_figure()
            fname = os.path.join(folder, view.title) + '.png'
            f.savefig(fname)


    def export(self, filename="", mode="w"):
        """Export the DataSet to the given *filename*.

        :param filename: The name of the file.
        :type filename: str

        .. warning::

            Not yet implemented
        """

        rule = '#' * 80 + '\n'

        def header():
            s = '# PARAMETERS:\n'
            groups = []
            for p in self.parameters():
                g = p['group']
                if g not in groups:
                    groups.append(g)
            parameters = {}
            for group in groups:
                parameters[group] = self.get_parameter(group=group)
            for k, v in parameters.items():
                if k == 'Cluster':
                    continue
                s += "# {}:\n".format(k)
                if not(isinstance(v, list)):
                    v = [v,]
                for p in v:
                    s += "#       {} = {} {}\n".format(p['name'], p['value'], p['unit'])
            return s

        colnames = self.columns()
        with open(filename, mode) as fd:
            # write the date and time of export
            now = datetime.now()
            fd.write("# Data exported on {}\n".format(now))
            fd.write(rule)

            # Append notes
            fd.write("# NOTES:\n")
            for line in self.notes.split('\n'):
                fd.write("#    {}\n".format(line))
            fd.write(rule)

            # Append parameters
            fd.write(header())
            fd.write(rule)

            # Append the data
            fd.write("# DATA:\n")
            fd.write("# " + ("{:<20s}" * len(colnames)).format(*colnames
                                                               ) + "\n")
            for i in range(len(self)):
                row = self[i]
                for key in row.columns():
                    value = row[key][0]
                    fmt = '{:s}'
                    #print value
                    for t, f in self._formats:
                        if isinstance(value, t):
                            fmt = f
                            break
                    #fd.write('  ')
                    fd.write(fmt.format(value))
                    #fd.write(str(value) + ', ')
                fd.write('\n')

    def __getitem__(self, itemspec):
        if isinstance(itemspec, str):
            return getattr(self, itemspec)
        title = 'untitled'
        new = DataSet(title)

        new._col_names = self.columns()
        for arr in self._col_arrays:
            new._col_arrays.append(np.asarray(arr)[itemspec].flatten())

        return new

    def __setstate__(self, state):
        self.__dict__ = state

    def __getstate__(self):
        return self.__dict__

    def __getattr__(self, name):
        if name in self._col_names:
            i = self._col_names.index(name)
            return self._col_arrays[i]
        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(
                self.__class__.__name__, name))

    def __setattr__(self, name, value):
        if name in self._col_names:
            i = self._col_names.index(name)
            self._col_arrays[i] = value
        else:
            self.__dict__[name] = value

    def __iter__(self):
        for i in range(len(self)):
            _ = {k: arr[i] for k, arr in zip(self._col_names,
                                          self._col_arrays)}
            point = _DataPoint(_)
            yield point

    def __len__(self):
        try:
            #length = len(self._col_arrays[0])
            length = 0
            for array in self._col_arrays:
                length = max(length, len(array))
        except IndexError:
            length = 0
        return length

    def __str__(self):
        max_len = 10
        max_col = 10
        ncols = min(max_col, len(self._col_arrays))
        table_data = [self._col_names[:ncols]]
        table_data[0].insert(0, "")

        all_indices = np.arange(0, len(self))
        indices = all_indices
        if len(self) > max_len:
            indices = list(range(int(max_len/2))) + list(range(int(-max_len/2), 0))

        _i = 0
        for i in indices:
            if i < _i:
                row = ['...' for _ in range(ncols + 1)]
                table_data.append(row)
            row = [str(all_indices[i]),]
            for j in range(ncols):
                arr = self._col_arrays[j]
                row.append(str(arr[i]))
            if len(self._col_names) > max_col:
                row.append('...')
            table_data.append(row)
            _i = i

        table = AsciiTable(table_data)
        table.outer_border = True
        table.title = self.title
        table.inner_column_border = False
        return table.table

    def __repr__(self):
        s = "<{}('{}')>".format(self.__class__.__name__, self.title)
        return s

class Data(object):
    """Creates a new Data object to store DataSets.

    :param title: The title of the Data object.
    :type str:

    """
    def __init__(self, title=''):
        assert is_title_valid(title), '\/:*?"<>| are not allowed in the string'
        self.title = title
        self._datasets = []
        self._dirty = False

    def add_dset(self, title, overwrite=False):
        """Adds a new DataSet in the Data object.

        :param title: The name of the DataSet.
        :type title: str
        :param overwrite: Tells whether to re-create the dataset if it exists.
        :type overwrite: bool
        :return: The newly created DataSet.
        :rtype: :py:class:`iodata.DataSet`
        """
        if overwrite:
            try:
                self.delete_dset(title)
            except Exception as err:
                pass

        titles = [d.title for d in self._datasets]
        if not title in titles:
            dset = DataSet(title)
            self._datasets.append(dset)
            self._dirty = True
            return dset
        else:
            raise NameError('A Dataset with that name already exists!')

    def delete_dset(self, title):
        """Removes a DataSet from the Data object.

        :param title: The DataSet name to be removed.
        :type title: str

        """
        titles = [d.title for d in self._datasets]
        i = titles.index(title)
        self._datasets.pop(i)
        self._dirty = True

    def get_last_dset(self):
        """Get the last DataSet of the Data object.

        :return: The lastly created DataSet in the Data object
        :rtype: :py:class:`iodata.DataSet`
        """
        return self._datasets[-1]

    def is_dirty(self):
        """Wether the Data object needs to be saved.

        :return: A boolean value to indicate if Data has changed since last dump to hard drive.
        :rtype: bool
        """
        return self._dirty


    def save(self, filename, append=False):
        """Saves the current Data to the hard drive.

        The Data, all its content along with parameters, defined views... are saved to the hard drive in the HDF5
        file format. Please see `hdfgroup <https://support.hdfgroup.org/HDF5/>`_ for more details about HDF5.

        :param filename: The name of the file to create or to append to.
        :type filename: str
        :param append: Wether to create a neww file or to append to an existing one.
        :type append: bool

        """
        mode = 'a' if append else 'w'
        titles = [d.title for d in self._datasets]
        with h5py.File(filename, mode) as fd:
            if append:
                try:
                    data_grp = fd['DATA']
                    meta_grp = fd['MsSpec viewer metainfo']
                except Exception as err:
                    fd.close()
                    self.save(filename, append=False)
                    return
            else:
                data_grp = fd.create_group('DATA')
                meta_grp = fd.create_group('MsSpec viewer metainfo')

            data_grp.attrs['title'] = self.title
            for dset in self._datasets:
                if dset.title in data_grp:
                    LOGGER.warning('dataset \"{}\" already exists in file \"{}\", not overwritting'.format(
                        dset.title, os.path.abspath(filename)))
                    continue
                grp = data_grp.create_group(dset.title)
                grp.attrs['notes'] = dset.notes
                for col_name in dset.columns():
                    data = dset[col_name]
                    grp.create_dataset(col_name, data=data)

            meta_grp.attrs['version'] = msspec.__version__

            root = etree.Element('metainfo')
            # xmlize views
            for dset in self._datasets:
                views_node = etree.SubElement(root, 'views', dataset=dset.title)
                for view in dset.views():
                    view_el = etree.fromstring(view.to_xml())
                    views_node.append(view_el)

            # xmlize parameters
            for dset in self._datasets:
                param_node = etree.SubElement(root, 'parameters', dataset=dset.title)
                for p in dset.parameters():
                    child = etree.SubElement(param_node, 'parameter')
                    for k, v in list(p.items()):
                        child.attrib[k] = v
            xml_str = etree.tostring(root, pretty_print=False)
            try:
                del meta_grp['info']
            except:
                pass
            finally:
                meta_grp.create_dataset('info', data=np.array((xml_str,)).view('S1'))
        self._dirty = False
        LOGGER.info('Data saved in {}'.format(os.path.abspath(filename)))

    def export(self, folder, overwrite=False):
        os.makedirs(folder, exist_ok=overwrite)
        for dset in self._datasets:
            dset_name = dset.title.replace(' ', '_')
            p = os.path.join(folder, dset_name)
            os.makedirs(p, exist_ok=overwrite)
            fname = os.path.join(p, dset_name) + '.txt'
            dset.export(fname)
            dset.export_views(p)

    @staticmethod
    def load(filename):
        """Loads an HDF5 file from the disc.

        :param filename: The path to the file to laod.
        :type filename: str
        :return: A Data object.
        :rtype: :py:class:`iodata.Data`
        """
        output = Data()
        with h5py.File(filename, 'r') as fd:
            parameters = {}
            views = {}

            output.title = fd['DATA'].attrs['title']
            for dset_name in fd['DATA'] :
                parameters[dset_name] = []
                views[dset_name] = []
                dset = output.add_dset(dset_name)
                dset.notes = fd['DATA'][dset_name].attrs['notes']
                for h5dset in fd['DATA'][dset_name]:
                    dset.add_columns(**{h5dset: fd['DATA'][dset_name][h5dset].value})

            try:
                vfile = LooseVersion(fd['MsSpec viewer metainfo'].attrs['version'])
                if vfile > LooseVersion(msspec.__version__):
                    raise NameError('File was saved with a more recent format')
                xml = fd['MsSpec viewer metainfo']['info'].value.tostring()
                root = etree.fromstring(xml)
                for elt0 in root.iter('parameters'):
                    dset_name = elt0.attrib['dataset']
                    for elt1 in elt0.iter('parameter'):
                        parameters[dset_name].append(elt1.attrib)

                for elt0 in root.iter('views'):
                    dset_name = elt0.attrib['dataset']
                    for elt1 in elt0.iter('view'):
                        view = _DataSetView(None, "")
                        view.from_xml(etree.tostring(elt1))
                        views[dset_name].append(view)

            except Exception as err:
                print(err)


            for dset in output:
                for v in views[dset.title]:
                    dset.add_view(v)
                for p in parameters[dset.title]:
                    dset.add_parameter(**p)

        output._dirty = False
        return output

    def __iter__(self):
        for dset in self._datasets:
            yield dset

    def __getitem__(self, key):
        try:
            titles = [d.title for d in self._datasets]
            i = titles.index(key)
        except ValueError:
            i = key
        return self._datasets[i]

    def __len__(self):
        return len(self._datasets)

    def __str__(self):
        s = str([dset.title for dset in self._datasets])
        return s

    def __repr__(self):
        s = "<Data('{}')>".format(self.title)
        return s

    def view(self):
        """Pops up a grphical window to show all the defined views of the Data object.

        """
        app = wx.App(False)
        app.SetAppName('MsSpec Data Viewer')
        frame = _DataWindow(self)
        frame.Show(True)
        app.MainLoop()


class _DataSetView(object):
    def __init__(self, dset, name, **plotopts):
        assert is_title_valid(name), '\/:*?"<>| are not allowed in the string'
        self.title = name
        self.dataset = dset
        self._plotopts = dict(
            title='No title',
            xlabel='', ylabel='', grid=True, legend=[], colorbar=False,
            projection='rectilinear', xlim=[None, None], ylim=[None, None],
            scale='linear',
            marker=None, autoscale=False)
        self._plotopts.update(plotopts)
        self._selection_tags = []
        self._selection_conditions = []

    def set_plot_options(self, **kwargs):
        self._plotopts.update(kwargs)

    def select(self, *args, **kwargs):
        condition = kwargs.get('where', 'True')
        legend = kwargs.get('legend', '')
        index = kwargs.get('index', None)
        if index is None:
            self._selection_conditions.append(condition)
            self._selection_tags.append(args)
            self._plotopts['legend'].append(legend)
        else:
            self._selection_conditions[index] = condition
            self._selection_tags[index] = args
            self._plotopts['legend'][index] = legend

    def tags(self):
        return self._selection_tags

    def get_data(self):
        data = []
        for condition, tags in zip(self._selection_conditions,
                                   self._selection_tags):
            indices = []
            # replace all occurence of tags
            for tag in self.dataset.columns():
                condition = condition.replace(tag, "p['{}']".format(tag))

            for i, p in enumerate(self.dataset):
                if eval(condition):
                    indices.append(i)

            values = []
            for tag in tags:
                values.append(getattr(self.dataset[indices], tag))

            data.append(values)
        return data

    def get_figure(self):
        opts = self._plotopts

        figure = Figure()
        axes = None
        proj = opts['projection']
        scale = opts['scale']
        if proj == 'rectilinear':
            axes = figure.add_subplot(111, projection='rectilinear')
        elif proj in ('polar', 'ortho', 'stereo'):
            axes = figure.add_subplot(111, projection='polar')

        for values, label in zip(self.get_data(), opts['legend']):
            # if we have only one column to plot, select a bar graph
            if np.shape(values)[0] == 1:
                xvalues = list(range(len(values[0])))
                axes.bar(xvalues, values[0], label=label,
                          picker=5)
                axes.set_xticks(xvalues)
            else:
                if proj in ('ortho', 'stereo'):
                    theta, phi, Xsec = cols2matrix(*values)
                    theta_ticks = np.arange(0, 91, 15)
                    if proj == 'ortho':
                        R = np.sin(np.radians(theta))
                        R_ticks = np.sin(np.radians(theta_ticks))
                    elif proj == 'stereo':
                        R = 2 * np.tan(np.radians(theta/2.))
                        R_ticks = 2 * np.tan(np.radians(theta_ticks/2.))
                    #R = np.tan(np.radians(theta/2.))
                    X, Y = np.meshgrid(np.radians(phi), R)
                    im = axes.pcolormesh(X, Y, Xsec)
                    axes.set_yticks(R_ticks)
                    axes.set_yticklabels(theta_ticks)

                    figure.colorbar(im)

                elif proj == 'polar':
                    values[0] = np.radians(values[0])
                    axes.plot(*values, label=label, picker=5,
                              marker=opts['marker'])
                else:
                    if scale == 'semilogx':
                        pltcmd = axes.semilogx
                    elif scale == 'semilogy':
                        pltcmd = axes.semilogy
                    elif scale == 'log':
                        pltcmd = axes.loglog
                    else:
                        pltcmd = axes.plot
                    pltcmd(*values, label=label, picker=5,
                           marker=opts['marker'])
            axes.grid(opts['grid'])
            axes.set_title(opts['title'])
            axes.set_xlabel(opts['xlabel'])
            axes.set_ylabel(opts['ylabel'])
            axes.set_xlim(*opts['xlim'])
            axes.set_ylim(*opts['ylim'])
            if label:
                axes.legend()
        axes.autoscale(enable=opts['autoscale'])

        canvas = FigureCanvasAgg(figure)
        return figure

    def serialize(self):
        data = {
            'name': self.title,
            'selection_conditions': self._selection_conditions,
            'selection_tags': self._selection_tags,
            'plotopts': self._plotopts
            }
        root = etree.Element('root')

        return data

    def to_xml(self):
        plotopts = self._plotopts.copy()
        legends = plotopts.pop('legend')

        root = etree.Element('view', name=self.title)
        for key, value in list(plotopts.items()):
            root.attrib[key] = str(value)
        #root.attrib['dataset_name'] = self.dataset.title

        for tags, cond, legend in zip(self._selection_tags,
                                      self._selection_conditions,
                                      legends):
            curve = etree.SubElement(root, 'curve')
            curve.attrib['legend'] = legend
            curve.attrib['condition'] = cond
            axes = etree.SubElement(curve, 'axes')
            for tag in tags:
                variable = etree.SubElement(axes, 'axis', name=tag)


        return etree.tostring(root, pretty_print=False)

    def from_xml(self, xmlstr):
        root = etree.fromstring(xmlstr)
        self.title = root.attrib['name']
        #self._plotopts['title'] = root.attrib['title']
        #self._plotopts['xlabel'] = root.attrib['xlabel']
        # self._plotopts['ylabel'] = root.attrib['ylabel']
        # self._plotopts['grid'] = bool(root.attrib['grid'])
        # self._plotopts['colorbar'] = bool(root.attrib['colorbar'])
        # self._plotopts['projection'] = root.attrib['projection']
        # self._plotopts['marker'] = root.attrib['marker']
        for key in list(self._plotopts.keys()):
            try:
                self._plotopts[key] = eval(root.attrib.get(key))
            except:
                self._plotopts[key] = root.attrib.get(key)



        legends = []
        conditions = []
        tags = []
        for curve in root.iter("curve"):
            legends.append(curve.attrib['legend'])
            conditions.append(curve.attrib['condition'])
            variables = []
            for var in curve.iter('axis'):
                variables.append(var.attrib['name'])
            tags.append(tuple(variables))

        self._selection_conditions = conditions
        self._selection_tags = tags
        self._plotopts['legend'] = legends

    def __repr__(self):
        s = "<{}('{}')>".format(self.__class__.__name__, self.title)
        return s

    def __str__(self):
        try:
            dset_title = self.dataset.title
        except AttributeError:
            dset_title = "unknown"
        s = '{}:\n'.format(self.__class__.__name__)
        s += '\tname    : %s\n' % self.title
        s += '\tdataset : %s\n' % dset_title
        s += '\ttags    : %s\n' % str(self._selection_tags)
        s += '\tconditions : %s\n' % str(self._selection_conditions)
        return s

class _GridWindow(wx.Frame):
    def __init__(self, dset, parent=None):
        title = 'Data: ' + dset.title
        wx.Frame.__init__(self, parent, title=title, size=(640, 480))
        self.create_grid(dset)

    def create_grid(self, dset):
        grid = wx.grid.Grid(self, -1)
        grid.CreateGrid(len(dset), len(dset.columns()))
        for ic, c in enumerate(dset.columns()):
            grid.SetColLabelValue(ic, c)
            for iv, v in enumerate(dset[c]):
                grid.SetCellValue(iv, ic, str(v))

class _ParametersWindow(wx.Frame):
    def __init__(self, dset, parent=None):
        title = 'Parameters: ' + dset.title
        wx.Frame.__init__(self, parent, title=title, size=(400, 480))
        self.create_tree(dset)

    def create_tree(self, dset):
        datatree = {}
        for p in dset.parameters():
            is_hidden = p.get('hidden', "False")
            if is_hidden == "True":
                continue
            group = datatree.get(p['group'], [])
            #strval = str(p['value'] * p['unit'] if p['unit'] else p['value'])
            #group.append("{:s} = {:s}".format(p['name'], strval))
            group.append("{} = {} {}".format(p['name'], p['value'], p['unit']))
            datatree[p['group']] = group

        tree = wx.TreeCtrl(self, -1)
        root = tree.AddRoot('Parameters')

        for key in list(datatree.keys()):
            item0 = tree.AppendItem(root, key)
            for item in datatree[key]:
                tree.AppendItem(item0, item)
        tree.ExpandAll()
        tree.SelectItem(root)

class _DataWindow(wx.Frame):
    def __init__(self, data):
        assert isinstance(data, (Data, DataSet))

        if isinstance(data, DataSet):
            dset = data
            data = Data()
            data.first = dset
        self.data = data
        self._filename = None
        self._current_dset = None

        wx.Frame.__init__(self, None, title="", size=(640, 480))

        self.Bind(wx.EVT_CLOSE, self.on_close)

        # Populate the menu bar
        self.create_menu()

        # Create the status bar
        statusbar = wx.StatusBar(self, -1)
        statusbar.SetFieldsCount(3)
        statusbar.SetStatusWidths([-2, -1, -1])
        self.SetStatusBar(statusbar)

        # Add the notebook to hold all graphs
        self.notebooks = {}
        sizer = wx.BoxSizer(wx.VERTICAL)
        #sizer.Add(self.notebook)
        self.SetSizer(sizer)

        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.on_page_changed)

        self.create_notebooks()

        self.update_title()

    def create_notebooks(self):
        for key in list(self.notebooks.keys()):
            nb = self.notebooks.pop(key)
            nb.Destroy()

        for dset in self.data:
            nb = wx.Notebook(self, -1)
            self.notebooks[dset.title] = nb
            #self.GetSizer().Add(nb, 1, wx.ALL|wx.EXPAND)
            self.GetSizer().Add(nb, proportion=1, flag=wx.ALL|wx.EXPAND)
            for view in dset.views():
                self.create_page(nb, view)

        self.create_menu()

        self.show_dataset(self.data[0].title)


    def create_menu(self):
        menubar = wx.MenuBar()
        menu1 = wx.Menu()
        menu1.Append(110, "Open\tCtrl+O")
        menu1.Append(120, "Save\tCtrl+S")
        menu1.Append(130, "Save as...")
        menu1.Append(140, "Export\tCtrl+E")
        menu1.AppendSeparator()
        menu1.Append(199, "Close\tCtrl+Q")

        menu2 = wx.Menu()
        for i, dset in enumerate(self.data):
            menu_id = 201 + i
            menu2.AppendRadioItem(menu_id, dset.title)
            self.Bind(wx.EVT_MENU, self.on_menu_dataset, id=menu_id)

        self.Bind(wx.EVT_MENU, self.on_open, id=110)
        self.Bind(wx.EVT_MENU, self.on_save, id=120)
        self.Bind(wx.EVT_MENU, self.on_saveas, id=130)
        self.Bind(wx.EVT_MENU, self.on_export, id=140)
        self.Bind(wx.EVT_MENU, self.on_close, id=199)


        menu3 = wx.Menu()
        menu3.Append(301, "Data")
        menu3.Append(302, "Cluster")
        menu3.Append(303, "Parameters")

        self.Bind(wx.EVT_MENU, self.on_viewdata, id=301)
        self.Bind(wx.EVT_MENU, self.on_viewcluster, id=302)
        self.Bind(wx.EVT_MENU, self.on_viewparameters, id=303)

        menubar.Append(menu1, "&File")
        menubar.Append(menu2, "&Datasets")
        menubar.Append(menu3, "&View")
        self.SetMenuBar(menubar)

    def on_open(self, event):
        if self.data.is_dirty():
            mbx = wx.MessageDialog(self, ('Displayed data is unsaved. Do '
                                          'you wish to save before opening'
                                          'another file ?'),
                                   'Warning: Unsaved data',
                                   wx.YES_NO | wx.ICON_WARNING)
            if mbx.ShowModal() == wx.ID_YES:
                self.on_saveas(wx.Event())
            mbx.Destroy()

        wildcard = "HDF5 files (*.hdf5)|*.hdf5"
        dlg = wx.FileDialog(
            self, message="Open a file...", defaultDir=os.getcwd(),
            defaultFile="", wildcard=wildcard, style=wx.FD_OPEN
            )

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self._filename = path
            self.data = Data.load(path)
            self.create_notebooks()
        dlg.Destroy()
        self.update_title()

    def on_save(self, event):
        if self._filename:
            if self.data.is_dirty():
                self.data.save(self._filename)
        else:
            self.on_saveas(event)

    def on_saveas(self, event):
        overwrite = True
        wildcard = "HDF5 files (*.hdf5)|*.hdf5|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Save file as ...", defaultDir=os.getcwd(),
            defaultFile='{}.hdf5'.format(self.data.title.replace(' ','_')),
            wildcard=wildcard, style=wx.FD_SAVE)
        dlg.SetFilterIndex(0)

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if os.path.exists(path):
                mbx = wx.MessageDialog(self, ('This file already exists. '
                                              'Do you wish to overwrite it ?'),
                                       'Warning: File exists',
                                       wx.YES_NO | wx.ICON_WARNING)
                if mbx.ShowModal() == wx.ID_NO:
                    overwrite = False
                mbx.Destroy()
            if overwrite:
                self.data.save(path)
                self._filename = path
        dlg.Destroy()
        self.update_title()

    def on_export(self, event):
        overwrite = True
        dlg = wx.DirDialog(
            self, message="Export data...", defaultPath=os.getcwd(),
            style=wx.DD_DEFAULT_STYLE)

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if os.listdir(path):
                mbx = wx.MessageDialog(self,
                                       ('This folder is not empty. '
                                        'Exporting tour data here may '
                                        'overwrite its content. Do you wish '
                                        'to continue ?'),
                                       'Warning: Folder is not empty',
                                       wx.YES_NO | wx.ICON_WARNING)
                if mbx.ShowModal() == wx.ID_NO:
                    overwrite = False
                mbx.Destroy()
            self.data.export(path, overwrite)
        dlg.Destroy()

    def on_viewdata(self, event):
        dset = self.data[self._current_dset]
        frame = _GridWindow(dset, parent=self)
        frame.Show()

    def on_viewcluster(self, event):
        win = wx.Frame(None, size=wx.Size(480, 340))
        cluster_viewer = ClusterViewer(win, size=wx.Size(480, 340))

        dset = self.data[self._current_dset]
        s = StringIO()
        s.write(dset.get_parameter(group='Cluster', name='cluster')['value'])
        atoms = ase.io.read(s, format='xyz')
        cluster_viewer.set_atoms(atoms, rescale=True, center=True)
        cluster_viewer.rotate_atoms(0., 180.)
        cluster_viewer.rotate_atoms(-45., -45.)
        #cluster_viewer.show_emitter(True)
        win.Show()

    def on_viewparameters(self, event):
        dset = self.data[self._current_dset]
        frame = _ParametersWindow(dset, parent=self)
        frame.Show()

    def on_close(self, event):
        if self.data.is_dirty():
            mbx = wx.MessageDialog(self, ('Displayed data is unsaved. Do you '
                                          'really want to quit ?'),
                                   'Warning: Unsaved data',
                                   wx.YES_NO | wx.ICON_WARNING)
            if mbx.ShowModal() == wx.ID_NO:
                mbx.Destroy()
                return
        self.Destroy()


    def on_menu_dataset(self, event):
        menu_id = event.GetId()
        dset_name = self.GetMenuBar().FindItemById(menu_id).GetText()
        self.show_dataset(dset_name)


    def show_dataset(self, name):
        for nb in list(self.notebooks.values()):
            nb.Hide()
        self.notebooks[name].Show()
        self.Layout()
        self.update_statusbar()
        self._current_dset = name

    def create_page(self, nb, view):
        # Get the matplotlib figure
        figure = view.get_figure()

        # Create a panel
        p = wx.Panel(nb, -1)

        # Create a matplotlib canvas for the figure
        canvas = FigureCanvas(p, -1, figure)
        sizer = wx.BoxSizer(wx.VERTICAL)

        toolbar = NavigationToolbar2WxAgg(canvas)
        toolbar.Realize()

        #sizer.Add(toolbar, 0, wx.ALL|wx.EXPAND)
        sizer.Add(toolbar, proportion=0, flag=wx.ALL|wx.EXPAND)
        toolbar.update()

        #sizer.Add(canvas, 5, wx.ALL|wx.EXPAND)
        sizer.Add(canvas, proportion=1, flag=wx.ALL|wx.EXPAND)

        p.SetSizer(sizer)
        p.Fit()
        p.Show()

        # MPL events
        figure.canvas.mpl_connect('motion_notify_event', self.on_mpl_motion)
        figure.canvas.mpl_connect('pick_event', self.on_mpl_pick)

        nb.AddPage(p, view.title)
        canvas.draw()


    def OLDcreate_page(self, nb, view):
        opts = view._plotopts
        p = wx.Panel(nb, -1)

        figure = Figure()

        axes = None
        proj = opts['projection']
        scale = opts['scale']
        if proj == 'rectilinear':
            axes = figure.add_subplot(111, projection='rectilinear')
        elif proj in ('polar', 'ortho', 'stereo'):
            axes = figure.add_subplot(111, projection='polar')

        canvas = FigureCanvas(p, -1, figure)
        sizer = wx.BoxSizer(wx.VERTICAL)

        toolbar = NavigationToolbar2WxAgg(canvas)
        toolbar.Realize()

        sizer.Add(toolbar, 0, wx.ALL|wx.EXPAND)
        toolbar.update()

        sizer.Add(canvas, 5, wx.ALL|wx.EXPAND)

        p.SetSizer(sizer)
        p.Fit()
        p.Show()


        for values, label in zip(view.get_data(), opts['legend']):
            # if we have only one column to plot, select a bar graph
            if np.shape(values)[0] == 1:
                xvalues = list(range(len(values[0])))
                axes.bar(xvalues, values[0], label=label,
                          picker=5)
                axes.set_xticks(xvalues)
            else:
                if proj in ('ortho', 'stereo'):
                    theta, phi, Xsec = cols2matrix(*values)
                    theta_ticks = np.arange(0, 91, 15)
                    if proj == 'ortho':
                        R = np.sin(np.radians(theta))
                        R_ticks = np.sin(np.radians(theta_ticks))
                    elif proj == 'stereo':
                        R = 2 * np.tan(np.radians(theta/2.))
                        R_ticks = 2 * np.tan(np.radians(theta_ticks/2.))
                    #R = np.tan(np.radians(theta/2.))
                    X, Y = np.meshgrid(np.radians(phi), R)
                    im = axes.pcolormesh(X, Y, Xsec)
                    axes.set_yticks(R_ticks)
                    axes.set_yticklabels(theta_ticks)

                    figure.colorbar(im)

                elif proj == 'polar':
                    values[0] = np.radians(values[0])
                    axes.plot(*values, label=label, picker=5,
                              marker=opts['marker'])
                else:
                    if scale == 'semilogx':
                        pltcmd = axes.semilogx
                    elif scale == 'semilogy':
                        pltcmd = axes.semilogy
                    elif scale == 'log':
                        pltcmd = axes.loglog
                    else:
                        pltcmd = axes.plot
                    pltcmd(*values, label=label, picker=5,
                           marker=opts['marker'])
            axes.grid(opts['grid'])
            axes.set_title(opts['title'])
            axes.set_xlabel(opts['xlabel'])
            axes.set_ylabel(opts['ylabel'])
            axes.set_xlim(*opts['xlim'])
            axes.set_ylim(*opts['ylim'])
            if label:
                axes.legend()
        axes.autoscale(enable=opts['autoscale'])


        # MPL events
        figure.canvas.mpl_connect('motion_notify_event', self.on_mpl_motion)
        figure.canvas.mpl_connect('pick_event', self.on_mpl_pick)

        nb.AddPage(p, view.title)


    def update_statusbar(self):
        sb = self.GetStatusBar()
        menu_id = self.GetMenuBar().FindMenu('Datasets')
        menu = self.GetMenuBar().GetMenu(menu_id)
        for item in menu.GetMenuItems():
            if item.IsChecked():
                sb.SetStatusText("%s" % item.GetText(), 1)
                break

    def update_title(self):
        title = "MsSpec Data Viewer"
        if self.data.title:
            title += ": " + self.data.title
        if self._filename:
            title += " [" + os.path.basename(self._filename) + "]"
        self.SetTitle(title)

    def on_mpl_motion(self, event):
        sb = self.GetStatusBar()
        try:
            txt = "[{:.3f}, {:.3f}]".format(event.xdata, event.ydata)
            sb.SetStatusText(txt, 2)
        except Exception:
            pass

    def on_mpl_pick(self, event):
        print(event.artist)

    def on_page_changed(self, event):
        self.update_statusbar()











if __name__ == "__main__":
    if True:
        data = Data('all my data')
        dset = data.add_dset('Dataset 0')
        X = np.arange(0, 20)
        Y = X**2

        dset.add_columns(x=X, y=Y, z=X+2, w=Y**3)
        dset.add_parameter(name='truc', group='main', value='3.14', unit='eV')
        dset.add_parameter(name='machin', group='main', value='abc', unit='')

    #     Z = [0,1]
    #
    #     for z in Z:
    #         for x, y in zip(X, Y):
    #             dset.add_row(x=x, y=y, z=z, random=np.random.rand())
    #
    #
        view = dset.add_view('my view', autoscale=True)
        view.select('x', 'y', where="z<10", legend=r"z = 0")
        view.select('x', 'y', where="z>10", legend=r"z = 1")
        print(dset.get_parameter(group='main'))
        constraint = lambda a, b: (a > 10 and a < 15) and b > 0
        indices = list(map(constraint, dset.x, dset.w))
        print(dset.y[indices])

        data.view()
    #import sys
    #data = Data.load(sys.argv[1])
    #data.view()
