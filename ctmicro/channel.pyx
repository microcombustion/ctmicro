# This file is part of the ctmicro add-on package to Cantera.
# See LICENSE in the top-level directory

__all__ += ['ChannelFlow']


cdef class ChannelFlow(ct._FlowBase):
    """
    new flow object
    """

    def __cinit__(self, ct._SolutionBase thermo, *args, **kwargs):
        gas = ct.getIdealGasPhase(thermo)
        self.channel = new CxxChannelFlow(gas, thermo.n_species, 2)
        self.flow = <ct.CxxStFlow * >(self.channel)

    property Nu:
        """ Nusselt Number """

        def __get__(self):
            return self.channel.getNusselt()

        def __set__(self, double nu):
            self.channel.setNusselt(nu)

    property diameter:
        """ scale  """

        def __get__(self):
            return self.channel.getDiameter()

        def __set__(self, double dia):
            self.channel.setDiameter(dia)

    def set_wall_profile(self, ProfileBase profile):
        """ pass call to cdef'd function """
        self.cxx_set_wall_profile(profile)

    cdef cxx_set_wall_profile(self, ProfileBase profile):
        """ set wall temperature profile """
        self.profile = profile
        self.channel.setWallProfile(profile._profile)

    def get_wall_profile(self):
        """ get wall temperature """
        cdef np.ndarray[np.double_t, ndim = 1] data = np.empty(self.n_points)
        for j, z in enumerate(self.grid):
            data[j] = self.profile.T(z)
        return data

    property n_states:
        def __get__(self):
            return self.channel.nStates()

    def get_state(self):
        cdef np.ndarray[np.double_t, ndim=1] data = np.zeros(self.n_states)

        self.channel.getState(&data[0])
        return data.reshape(self.n_points, self.n_components).T

    def get_residual(self):
        cdef np.ndarray[np.double_t, ndim=1] data = np.zeros(self.n_states)

        self.channel.getResidual(&data[0])
        return data.reshape(self.n_points, self.n_components).T

    def get_t_mask(self):
        cdef np.ndarray[np.double_t, ndim=1] data = np.zeros(self.n_states)

        self.channel.getTMask(&data[0])
        return data.reshape(self.n_points, self.n_components).T.astype(np.int)
