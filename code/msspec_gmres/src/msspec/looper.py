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
# Source file  : src/msspec/looper.py
# Last modified: Wed, 26 Feb 2025 11:15:54 +0100
# Committed by : Sylvain Tricot <sylvain.tricot@univ-rennes.fr>

from collections import OrderedDict
from functools import partial
import itertools
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import time


logger = logging.getLogger(__name__)


class Variable:
    def __init__(self, name, doc=""):
        self.name = name
        self.doc = doc

    def __repr__(self):
        return "<Variable(\'{}\')>".format(self.name)

class Sweep:
    def __init__(self, key, comments="", unit=None,
                 values=None,
                 start=None, stop=None, step=None,
                                   num=10, scale='lin', base=10,
                 default=None,
                 folded=False,
                 unzip=False,
                 group=None,
                 accumulations=1,
                 linked_to=None):
        self.key = key    # The variable's name
        self.comments = comments
        self.default = default
        self.folded = folded
        self.unzip = unzip
        self.linked_to = linked_to
        self.group = None
        # First use case: values are specidied
        if values is not None:
            self.values = values
        else:
            assert start is not None and stop is not None
            self.start = start
            self.stop = stop
            if step is not None:
                self.step = step
                self.values = np.arange(start, stop, step)
            else:
                self.num = num
                if scale == 'lin':
                    self.values = np.linspace(start, stop, num)
                elif scale == 'log':
                    self.base = base
                    self.values = np.logspace(start, stop, num, base=base)

    def __getitem__(self, index):
        return (self.key, self.values[index])

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return "<{}(\'{}\')>".format(self.__class__.__name__, self.key)

    @property
    def unfolded(self):
        return not(self.folded)


class SweepRange:
    def __init__(self, *sweeps):
        self.sweeps = sweeps
        self.index = 0

        # First check that sweeps that are linked to another on are all included
        links = {}
        for sweep in sweeps:
            if sweep.linked_to is not None:
                if sweep.linked_to not in sweeps:
                    raise NameError(("The sweep \'{}\' is linked to \'{}\' "
                                     "but is not included in the loop!").format(
                                        sweep.key, sweep.linked_to.key))
                # add the linked sweep to the list of links
                if sweep.linked_to not in links.keys():
                    links[sweep.linked_to] = [sweep,]
                else:
                    links[sweep.linked_to].append(sweep)

        # The cumulative product of lengths
        lengths = []
        for sweep in sweeps:
            if sweep.linked_to is None:
                lengths.append(len(sweep) if sweep.unfolded else 1)
            else:
                lengths.append(1)
        cn = np.cumprod(lengths)
        # Get the total number of combinations is the last one
        ntot = cn[-1]

        self.links = links
        self.cn = cn
        self.ntot = ntot

    def __len__(self):
        return self.ntot

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index < self.ntot:
            i = self.index
            self.index += 1

            row = OrderedDict({k: None for k in self.columns})
            for isweep, sweep in enumerate(self.sweeps):
                # If sweep is linked to another one, do nothing as its value
                # will be added by its parent
                if sweep.linked_to is not None:
                    continue

                # Compute the index
                if sweep.folded:
                    key, value = sweep.key, sweep.values
                    row[key] = value
                else:
                    idx = int(i/(self.ntot/self.cn[isweep])) % len(sweep)

                    # If this sweep has links, add also all values from its
                    # children
                    children = self.links.get(sweep, [])
                    for s in [sweep,] + children:
                        key, value = s[idx]
                        row[key] = value
                row['sweep_index'] = i
            return row
        else:
            raise StopIteration

    @property
    def columns(self):
        cols  = ['sweep_index']
        cols += [sweep.key for sweep in self.sweeps]
        return cols

    @property
    def values(self):
        return [list(row.values()) for row in self]

    @property
    def items(self):
        items = []
        for row in self:
            items.append(row)
        return items


class Looper:
    def __init__(self):
        self.data = None

    @property
    def pipeline(self):
        try:
            return self._pipeline
        except AttributeError as error:
            return None

    @pipeline.setter
    def pipeline(self, value):
        self._pipeline = value

    def _wrapper(self, x):
        logger.debug("Pipeline called with {}".format(x))
        return self.pipeline(**x)

    def run(self, *sweeps, ncpu=1, **kwargs):
        logger.info("Loop starts...")
        # prepare the list of inputs
        sr = SweepRange(*sweeps)
        items = sr.items

        data = []

        t0 = time.time()

        if ncpu == 1:
            # serial processing...
            logger.info("serial processing...")
            for i, values in enumerate(items):
                values.update(kwargs)
                result = self._wrapper(values)
                data.append(result)
        else:
            # Parallel processing...
            chunksize = 1 #int(nsets/ncpu)
            [values.update(kwargs) for values in items]
            logger.info(("Parallel processing over {:d} cpu (chunksize={:d})..."
                         "").format(ncpu, chunksize))
            t0 = time.time()

            pool = mp.Pool(processes=ncpu)
            data = pool.map(self._wrapper, items, chunksize=chunksize)
            pool.close()
            pool.join()

        t1 = time.time()
        dt = t1 - t0
        logger.info(("Processed {:d} sets of inputs in {:.3f} seconds"
                        "").format(len(sr), dt))

        # Create the DataFrame
        dfdata = []
        columns = sr.columns + list(kwargs.keys()) + ['output',]

        for i in range(len(sr)):
            row = list(items[i].values())
            row.append(data[i])
            dfdata.append(row)


        df = pd.DataFrame(dfdata, columns=columns)
        df = df.drop(columns=['sweep_index'])

        self.data = df

        #return df
        # return a list of [x0, y0, x1, y1,...xn, yn] and a list
        # of corresponding dict of parameters {'keyA': [val0,...valn],
        # 'keyB': [val0,...valn], ...}

        # all_xy = []
        # for irow, row in df.iterrows():
            # all_xy.append(row.output[0])
            # all_xy.append(row.output[1])
        # parameters = df.to_dict()
        # parameters.pop('output')

        return self.data #all_xy, parameters






if __name__ == "__main__":
    import numpy as np
    import time
    import logging


    logging.basicConfig(level=logging.DEBUG)



    def bar(**kwargs):
        i = kwargs.get('sweep_index')
        return np.linspace(0,i,10)

    theta = Sweep(key='theta', comments="The polar angle",
                  start=-70, stop=70, num=3,
                  unit='degree',
                  folded=True)

    phi = Sweep('phi', comments="The azimutal angle",
                values=[0, 45],
                unit='degree',
                folded=True)

    emitter_plane = Sweep('emitter_plane', comments="The emitter plane",
                          start=0, stop=3, step=1)

    emitter = Sweep(key='emitter', values=('Ti', 'Sr'))

    levels = Sweep(key='level', values=('2p', '3d'), linked_to=emitter)

    lmaxt = Sweep(key='lmaxt', values=(25, 29, 30))#, linked_to=emitter)

    uij = Sweep(key='uij', values=(0.01, 0.02, 0.03))#, linked_to=lmaxt)

    sweeps = [theta, phi, emitter_plane, emitter, levels, lmaxt, uij]


    looper = Looper()
    looper.pipeline = bar
    other_kws = {'un':1, 'deux':2}
    data = looper.run(emitter, emitter_plane, uij, theta, levels, ncpu=4, **other_kws)

    # Print the dataframe
    print(data)

    # Accessing the parameters and ouput values for a given sweep (e.g the last one)
    print(looper.data.iloc[-1])

    # Post-process the output values. For example here, the output is a 1D-array,
    # make the sum of sweeps with 'Sr' emitter
    X = np.array([ x for x in data[data.emitter == 'Sr'].output]).sum(axis=0)
    print(X)
