# This file is part of the ctmicro add-on package to Cantera.
# See LICENSE in the top-level directory

import csv
import cantera as ct
import numpy as np
from ._ctmicro import ChannelFlow, ProfileBase


__all__ = ['ChannelFlame']


class ChannelFlame(ct.FlameBase):
    """A flame in a channel."""
    __slots__ = ('inlet', 'outlet', 'flame')

    def __init__(self, gas, profile=None):
        """
        A domain of type IdealGasFlow named 'flame' will be created to represent
        the flame and set to free flow. The three domains comprising the stack
        are stored as ``self.inlet``, ``self.flame``, and ``self.outlet``.

        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on a fixed grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.
        """
        assert isinstance(profile, ProfileBase), 'invalid profile'

        self.inlet = ct.Inlet1D(name='inlet', phase=gas)
        self.outlet = ct.Outlet1D(name='outlet', phase=gas)
        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            self.flame = ChannelFlow(gas, name='channel')
            # self.flame.set_free_flow()
            self.flame.set_wall_profile(profile)

        z_up, z_dn = profile.z_range
        grid = np.linspace(z_up, z_dn, 10)

        super(ChannelFlame, self).__init__((self.inlet, self.flame, self.outlet),
                                           gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.inlet.T = gas.T
        self.inlet.X = gas.X

#     def __init__(self, gas, grid=None, **kwargs):
#         """
#         :param gas:
#             `Solution` (using the IdealGas thermodynamic model) used to
#             evaluate all gas properties and reaction rates.
#         :param grid:
#             Array of initial grid points

#         A domain of class `ChannelFlow` named ``channel`` will
#         be created to represent the channel. The three domains comprising the
#         stack are stored as ``self.inlet``, ``self.flame``, and
#         ``self.outlet``.
#         """
#         self.inlet = ct.Inlet1D(name='inlet', phase=gas)  # Channel
#         self.inlet.T = gas.T
#         self.outlet = ct.Outlet1D(name='outlet', phase=gas)  # Channel
#         self.flame = ChannelFlow(gas, name='channel')

#         super().__init__((self.inlet, self.flame, self.outlet), gas, grid)

    def set_initial_guess(self, locs=[0.0, 0.3, 0.5, 1.0]):
        """
        Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the inlet gas
        composition.

        :param locs:
            A list of four locations to define the temperature and mass fraction profiles.
            Profiles rise linearly between the second and third location.
            Locations are given as a fraction of the entire domain
        """
        super(ChannelFlame, self).set_initial_guess()
        self.gas.TPY = self.inlet.T, self.P, self.inlet.Y

        if not self.inlet.mdot:
            # nonzero initial guess increases likelihood of convergence
            self.inlet.mdot = 1.0 * self.gas.density

        Y0 = self.inlet.Y
        u0 = self.inlet.mdot / self.gas.density
        T0 = self.inlet.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.inlet.mdot / self.gas.density

        self.set_profile('u', locs, [u0, u0, u1, u1])
        self.set_profile('T', locs, [T0, T0, Teq, Teq])

        # Pick the location of the fixed temperature point, using an existing
        # point if a reasonable choice exists
        T = self.T
        Tmid = 0.75 * T0 + 0.25 * Teq
        i = np.flatnonzero(T < Tmid)[-1]  # last point less than Tmid
        if Tmid - T[i] < 0.2 * (Tmid - T0):
            self.set_fixed_temperature(T[i])
        elif T[i + 1] - Tmid < 0.2 * (Teq - Tmid):
            self.set_fixed_temperature(T[i + 1])
        else:
            self.set_fixed_temperature(Tmid)

        for n in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(n),
                             locs, [Y0[n], Y0[n], Yeq[n], Yeq[n]])

#     def set_initial_guess(self, lw=.6):
#         """
#         Set the initial guess for the solution. The adiabatic channel
#         temperature and equilibrium composition are computed for the inlet
#         gas composition. The temperature profile rises linearly in the first
#         20% of the channel to Tad, then is flat. The mass fraction profiles are
#         set similarly.
#         """
#         super().set_initial_guess()

#         Y0 = self.inlet.Y
#         u0 = self.inlet.mdot / self.gas.density
#         T0 = self.inlet.T

#         locs = [.0, lw, lw + .1, 1.0]
#         tw = self.flame.get_wall_profile()
#         T1 = (1. - lw) * tw[0] + lw * tw[-1]
#         T2 = tw[-1]

#         self.gas.TPY = T1, self.P, Y0

#         # get adiabatic channel temperature and composition
#         self.gas.equilibrate('HP')
#         Teq = self.gas.T
#         Yeq = self.gas.Y
#         u1 = self.inlet.mdot / self.gas.density
#         print(T1)
#         print(Teq)
#         # raise(TypeError('stop'))

#         # raise(TypeError('stop'))
#         print(locs)
#         print([T0, T1, Teq, T2])
#         self.set_profile('u', locs, [u0, u0 * T1 / T0, u1, u1 * T2 / Teq])
#         self.set_profile('T', locs, [T0, T1, Teq, T2])
#         for n in range(self.gas.n_species):
#             self.set_profile(self.gas.species_name(n),
#                              locs, [Y0[n], Y0[n], Yeq[n], Yeq[n]])

    def solve(self, loglevel=1, refine_grid=True, auto=False):
        """
        Solve the problem.

        :param loglevel:
            integer flag controlling the amount of diagnostic output. Zero
            suppresses all output, and 5 produces very verbose output.
        :param refine_grid:
            if True, enable grid refinement.
        :param auto: if True, sequentially execute the different solution stages
            and attempt to automatically recover from errors. Attempts to first
            solve on the initial grid with energy enabled. If that does not
            succeed, a fixed-temperature solution will be tried followed by
            enabling the energy equation, and then with grid refinement enabled.
            If non-default tolerances have been specified or multicomponent
            transport is enabled, an additional solution using these options
            will be calculated.
        """
        if not auto:
            return super(ChannelFlame, self).solve(loglevel, refine_grid, auto)

        # Use a callback function to check that the domain is actually wide
        # enough to contain the flame after each steady-state solve. If the user
        # provided a callback, store this so it can called in addition to our
        # callback, and restored at the end.
        original_callback = self._steady_callback

        class DomainTooNarrow(Exception):
            pass

        def check_width(t):
            T = self.T
            x = self.grid
            mRef = (T[-1] - T[0]) / (x[-1] - x[0])
            mLeft = (T[1] - T[0]) / (x[1] - x[0]) / mRef
            mRight = (T[-3] - T[-1]) / (x[-3] - x[-1]) / mRef

            # The domain is considered too narrow if gradient at the left or
            # right edge is significant, compared to the average gradient across
            # the domain.
            if mLeft > 0.02 or mRight > 0.02:
                raise DomainTooNarrow()

            if original_callback:
                return original_callback(t)
            else:
                return 0.0

        self.set_steady_callback(check_width)

        for _ in range(12):
            try:
                return super(ChannelFlame, self).solve(loglevel, refine_grid, auto)
            except DomainTooNarrow:
                self.flame.grid *= 2
                if loglevel > 0:
                    print('Expanding domain to accomodate flame thickness. '
                          'New width: {} m'.format(
                              self.flame.grid[-1] - self.flame.grid[0]))
                if refine_grid:
                    self.refine(loglevel)

        self.set_steady_callback(original_callback)

    #@property
    def Tw(self):
        """ Wall temperature [K] at this point. """
        # note that we're using V's allocated memory
        return self.flame.get_wall_profile()  # profile(self.flame,'V')

    def write_csv(self, filename, species='X', quiet=True):
        """
        Write the velocity, temperature, density, and species profiles
        to a CSV file.

        :param filename:
            Output file name
        :param species:
            Attribute to use obtaining species profiles, e.g. ``X`` for
            mole fractions or ``Y`` for mass fractions.
        """

        z = self.grid
        T = self.T
        u = self.u
        Tw = self.flame.get_wall_profile()

        csvfile = open(filename, 'w')
        writer = csv.writer(csvfile)
        writer.writerow(['z (m)', 'u (m/s)', 'Tw (K)',
                         'T (K)', 'rho (kg/m3)'] + self.gas.species_names)
        for n in range(self.flame.n_points):
            self.set_gas_state(n)
            writer.writerow([z[n], u[n], Tw[n], T[n], self.gas.density] +
                            list(getattr(self.gas, species)))
        csvfile.close()
        if not quiet:
            print("Solution saved to '{0}'.".format(filename))
