# This file is part of the ctmicro add-on package to Cantera.
# See LICENSE in the top-level directory

from numpy import polyval

__all__ += ['ConstProfile', 'LinearProfile', 'ErfProfile',
            'InterpolatedProfile', 'PolynomialProfile']


cdef class ProfileBase(object):
    """
    Base type for wall temperature profile objects
    """
    profile_type = ""

    def __cinit__(self):
        self._profile = CxxNewProfile(stringify((self.profile_type)))
        self.profile = self._profile.get()

    def __init__(self, (double, double) z_range, (double, double) T_range):
        self.profile.setup(z_range[0], z_range[1], T_range[0], T_range[1])

    property z_range:
        """ Position range (upstream, downstream) """

        def __get__(self):
            return (self.profile.getZup(), self.profile.getZdn())

    property T_range:
        """ Temperature range (upstream, downstream) """

        def __get__(self):
            return (self.profile.getTup(), self.profile.getTdn())

    def T(self, double z):
        """ Calculate temperature at location z """

        return self.profile.calcT(z)


cdef class ConstProfile(ProfileBase):
    """
    A profile with constant wall temperature
    """
    profile_type = "constant"

    def __init__(self, (double, double) z_range=(0., 1.), double T_wall=300.):
        super().__init__(z_range, (T_wall, T_wall))


cdef class LinearProfile(ProfileBase):
    """
    A linear wall temperature profile
    """
    profile_type = "linear"

    def __init__(self, (double, double) z_range=(0., 1.), (double, double) T_range=(300., 1500.)):
        super().__init__(z_range, T_range)


cdef class ErfProfile(ProfileBase):
    """
    A wall temperature profile based on an error function
    """
    profile_type = "error-function"

    def __cinit__(self):
        self.cxx = <CxxErfProfile *>(self.profile)

    def __init__(self, (double, double) z_range=(-1., 1.), (double, double) T_range=(300., 1500.),
                 double mu=0., double sigma=1.):
        super().__init__(z_range, T_range)
        self.mu = mu
        self.sigma = sigma

    property mu:
        """ Error function mid point """

        def __get__(self):
            return self.cxx.mu()

        def __set__(self, double mu):
            self.cxx.setMu(mu)

    property sigma:
        """ Error function width """

        def __get__(self):
            return self.cxx.sigma()

        def __set__(self, double sigma):
            self.cxx.setSigma(sigma)


cdef class InterpolatedProfile(ProfileBase):
    """
    A wall temperature profile based on interpolated nodes
    """
    profile_type = "interpolated"

    def __cinit__(self):
        self.cxx = <CxxInterpolatedProfile *>(self.profile)

    def __init__(self, zfixed=[0., .2, .5, 1.], tfixed=[300., 400, 1200, 1500.]):
        cdef vector[double] x, y
        for z in zfixed:
            x.push_back(z)
        for t in tfixed:
            y.push_back(t)
        super().__init__((zfixed[0], zfixed[-1]), (tfixed[0], tfixed[-1]))
        self.cxx.setValues(x, y)


cdef class PolynomialProfile(ProfileBase):
    """
    A wall temperature profile with a polynomial progression
    """
    profile_type = "polynomial"

    def __cinit__(self):
        self.cxx = <CxxPolynomialProfile *>(self.profile)

    def __init__(self, z_range=(0., 1.), coeffs=[1200., 0., 300.]):
        """ C++ polyN functions use reverted indexing compared to polyval """
        cdef vector[double] x
        for c in coeffs[-1::-1]:
            x.push_back(c)
        cdef size_t n = len(coeffs)
        assert n > 2 and n < 7, 'invalid polynomial order'
        T_range = polyval(coeffs, z_range)
        super().__init__(z_range, tuple(T_range))
        self.cxx.setValues(n, x)
