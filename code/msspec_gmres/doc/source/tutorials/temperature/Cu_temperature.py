# coding: utf8

from ase.build import bulk
import numpy as np
from msspec.calculator import MSSPEC, XRaySource
from msspec.utils import hemispherical_cluster, get_atom_index


def create_clusters(nplanes=3):
    copper  = bulk('Cu', a=3.6, cubic=True)
    clusters = []
    for emitter_plane in range(nplanes):
        cluster = hemispherical_cluster(copper,
                                        emitter_plane=emitter_plane,
                                        planes=emitter_plane+2,
                                        shape='cylindrical')
        cluster.absorber = get_atom_index(cluster, 0, 0, 0)
        cluster.info.update({
            'emitter_plane': emitter_plane,
        })
        clusters.append(cluster)
    return clusters

def compute(clusters, all_theta=[45., 83.],
            all_T=np.arange(300., 1000., 400.)):
    data = None
    for ic, cluster in enumerate(clusters):
        # Retrieve info from cluster object
        plane   = cluster.info['emitter_plane']

        calc = MSSPEC(spectroscopy='PED', algorithm='expansion')
        calc.source_parameters.energy = XRaySource.AL_KALPHA
        calc.muffintin_parameters.interstitial_potential = 14.1

        calc.calculation_parameters.vibrational_damping = 'averaged_tl'
        calc.calculation_parameters.use_debye_model = True
        calc.calculation_parameters.debye_temperature = 343
        calc.calculation_parameters.vibration_scaling = 1.2

        calc.detector_parameters.average_sampling = 'high'
        calc.detector_parameters.angular_acceptance = 5.7

        for atom in cluster:
            atom.set('forward_angle', 30)
            atom.set('backward_angle', 30)
        filters = ['forward_scattering', 'backward_scattering']
        calc.calculation_parameters.path_filtering = filters

        calc.calculation_parameters.RA_cutoff = 2

        for T in all_T:
            for theta in all_theta:
                calc.calculation_parameters.temperature = T
                calc.calculation_parameters.scattering_order = min(1 + plane, 3)
                calc.set_atoms(cluster)
                data = calc.get_phi_scan(level='2p3/2', theta=theta,
                                         phi=np.linspace(0, 100),
                                         kinetic_energy=560, data=data)
                dset = data[-1]
                dset.title = "plane #{:d}, T={:f}K, theta={:f}°".format(plane,
                                                                        T,
                                                                        theta)

                dset.add_parameter(group='misc', name='plane', value=plane, unit='')
                dset.add_parameter(group='misc', name='T', value=T, unit='')
                dset.add_parameter(group='misc', name='theta', value=theta, unit='')
    return data


def analysis(data):
    all_plane = []
    all_T     = []
    all_theta = []
    for dset in data:
        plane = dset.get_parameter('misc', 'plane')['value']
        T     = dset.get_parameter('misc', 'T')['value']
        theta = dset.get_parameter('misc', 'theta')['value']
        cs    = dset.cross_section.copy()
        phi   = dset.phi.copy()

        if plane not in all_plane: all_plane.append(plane)
        if T not in all_T: all_T.append(T)
        if theta not in all_theta: all_theta.append(theta)

    def get_anisotropy(theta, T):
        cs = None
        for dset in data:
            try:
                _T     = dset.get_parameter('misc', 'T')['value']
                _theta = dset.get_parameter('misc', 'theta')['value']
                _cs    = dset.cross_section.copy()
                phi    = dset.phi.copy()
            except:
                continue
            if _theta == theta and _T == T:
                try:
                    cs += _cs
                except:
                    cs = _cs
        Imax = np.max(cs)
        Imin = np.min(cs)
        return (Imax - Imin)/Imax

    # create a substrate dataset for each T and theta
    anisotropy_dset = data.add_dset("all")
    anisotropy_view = anisotropy_dset.add_view('Anisotropies',
                         title='Relative anisotropies for Cu(2p)',
                         marker='o',
                         xlabel='T (K)',
                         ylabel=r'$\frac{\Delta I / I_{max}(T)}{\Delta I_{300}'
                                r'/ I_{max}(300)} (\%)$')

    for theta in all_theta:
        for T in all_T:
            A  = get_anisotropy(theta, T)
            A0 = get_anisotropy(theta, np.min(all_T))

            anisotropy_dset.add_row(temperature=T, theta=theta, anisotropy=A/A0)

        anisotropy_view.select('temperature', 'anisotropy',
                               where='theta == {:f}'.format(theta),
                               legend=r'$\theta = {:.0f} \degree$'.format(theta))
    return data


clusters = create_clusters()
data     = compute(clusters)
data     = analysis(data)
data.view()
