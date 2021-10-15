from collections import defaultdict, OrderedDict
from functools import partial

import astropy.units as u
import numpy as np
from astropy.table.table import Table
from astropy.units.quantity import Quantity
from roentgen.absorption.material import Compound, MassAttenuationCoefficient, Material

from stix.fits.products.quicklook import ENERGY_CHANNELS

MIL_SI = 0.0254 * u.mm

COMPONENTS = OrderedDict([
    ('front_window', [('solarblack', 0.005 * u.mm), ('be', 2*u.mm)]),
    ('rear_window', [('be', 1*u.mm)]),
    ('grid_covers', [('kapton', 4 * 2 * MIL_SI )]),
    ('dem', [('kapton', 2 * 3 * MIL_SI )]),
    ('attenuator', [('al', 0.6 * u.mm)]),
    ('mli', [('al', 1000 * u.angstrom), ('kapton', 3 * MIL_SI), ('al', 40 * 1000 * u.angstrom),
            ('mylar', 20 * 0.25 * MIL_SI), ('pet', 21 * 0.005 * u.mm), ('kapton', 3 * MIL_SI),
            ('al', 1000 * u.angstrom)]),
    ('calibration_foil', [('al', 4 * 1000 * u.angstrom), ('kapton', 4 * 2 * MIL_SI)]),
    ('dead_layer', [('te_o2', 392 * u.nm)]),
])

MATERIALS = OrderedDict([
    ('al', ({'Al': 1.0}, 2.7 * u.g/u.cm**3)),
    ('be', ({'Be': 1.0}, 1.85 * u.g/u.cm**3)),
    ('kapton', ({'H': 0.026362, 'C': 0.691133, 'N': 0.073270, 'O': 0.209235},
               1.43 * u.g / u.cm ** 3)),
    ('mylar', ({'H': 0.041959, 'C': 0.625017, 'O': 0.333025}, 1.38 * u.g / u.cm ** 3)),
    ('pet', ({'H': 0.041960, 'C': 0.625016, 'O': 0.333024}, 1.370 * u.g / u.cm**3)),
    ('solarblack_oxygen', ({'H': 0.002, 'O': 0.415, 'Ca': 0.396, 'P': 0.187}, 3.2 * u.g/u.cm**3)),
    ('solarblack_carbon', ({'C': 0.301, 'Ca': 0.503, 'P': 0.195}, 3.2 * u.g/u.cm**3)),
    ('te_o2', ({'Te': 0.7995088158691722, 'O': 0.20049124678825841}, 5.670 * u.g/u.cm**3))]
)


class Response:
    """
    Calculate the response of the
    """
    def __init__(self, solarblack='solarblack_carbon'):
        if solarblack not in ['solarblack_oxygen', 'solarblack_carbon']:
            raise ValueError('solarblack must be either solarblack_oxygen or solarblack_carbon.')

        self.solarblack = solarblack
        self.materials = MATERIALS
        self.components_ = COMPONENTS
        self.components = dict()
        self.energies = [ENERGY_CHANNELS[i]['e_lower'] for i in range(1, 32)] * u.keV

        for name, layers in COMPONENTS.items():
            parts = []
            for material, thickness in layers:
                if material == 'solarblack':
                    material = self.solarblack
                mass_frac, den = MATERIALS[material]
                parts.append(self.create_material(name=material, fractional_masses=mass_frac,
                                                  thickness=thickness, density=den))
            self.components[name] = Compound(parts)

    def get_transmission(self, energies=None, attenuator=False):
        base_comps = [self.components[name] for name in ['front_window', 'rear_window', 'dem',
                                                         'mli', 'calibration_foil',
                                                         'dead_layer']]

        if energies is None:
            energies = self.energies

        if attenuator:
            base_comps.append(self.components['attenuator'])

        base = Compound(base_comps)
        base_trans = base.transmission(energies[:-1] + 0.5 * np.diff(energies))
        fine = Compound(base_comps + [self.components['grid_covers']])
        fine_trans = fine.transmission(energies[:-1] + 0.5 * np.diff(energies))

        fine_grids = np.array([11, 13, 18, 12, 19, 17]) - 1
        transmission = Table()
        # transmission['sci_channel'] = range(1, 31)
        for i in range(33):
            name = f'det-{i}'
            if np.isin(i, fine_grids):
                transmission[name] = fine_trans
            else:
                transmission[name] = base_trans
        return transmission

    def get_transmission_by_compnent(self):
        return self.components

    def get_transmission_by_material(self):
        material_thickness = dict()
        for name, layers in COMPONENTS.items():
            for material_name, thickness in layers:
                if material_name == 'solarblack':
                    material_name = self.solarblack
                if material_name in material_thickness.keys():
                    material_thickness[material_name] += thickness.to('mm')
                else:
                    material_thickness[material_name] = thickness.to('mm')
        res = {}
        for name, thickness in material_thickness.items():
            frac_mass, density = self.materials[name]
            mat = self.create_material(name=name, fractional_masses=frac_mass,
                                       density=density, thickness=thickness)
            res[name] = mat

        return res

    @classmethod
    def create_material(cls, name=None, fractional_masses=None, thickness=None, density=None):
        material = Material('h', thickness, density)
        material.name = name
        # probbably don't need this
        material.density = density

        def func(fractional_masses, e):
            return sum([MassAttenuationCoefficient(element).func(e) * frac_mass
                        for element, frac_mass in fractional_masses.items()])

        material_func = partial(func, fractional_masses)

        material.mass_attenuation_coefficient.func = material_func
        return material

    def _bin_average(self, obj, n_points=1001):
        out = []
        for i in range(self.energies.size-1):
            ee = np.linspace(self.energies[i], self.energies[i+1], n_points)
            vals = obj.transmission(ee)
            avg = np.average(vals)
            out.append(avg)
        return np.hstack(out)
