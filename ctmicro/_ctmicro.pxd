# This file is part of the ctmicro add-on package to Cantera.
# See LICENSE in the top-level directory
# cython: embedsignature=True
# cython: language_level=3
# distutils: language = c++
'''
_ctmixcro.pxd

This file is part of the ctmicro add-on package to Cantera.
See License.txt in the top-level directory
'''

#import numpy as np
# cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as cbool

cdef extern from "<memory>":
    cppclass shared_ptr "std::shared_ptr" [T]:
        T* get()
        void reset(T*)

cdef string stringify(x) except *


# cimport _cantera
cimport cantera as ct

cdef extern from "cantera/cython/funcWrapper.h":
    cdef int translate_exception()


cdef extern from "ChannelFlow.h":
    cdef cppclass CxxChannelFlow "CanteraApp::ChannelFlow":
        CxxChannelFlow(ct.CxxIdealGasPhase *, int, int)
        void setNusselt(double)
        double getNusselt()
        void setDiameter(double)
        double getDiameter()
        void setWallProfile(shared_ptr[CxxProfile]) except +translate_exception
        double Twall(int)
        size_t nStates()
        void getState(double *)
        void getResidual(double *)
        void getTMask(double *)

cdef class ChannelFlow(ct._FlowBase):
    cdef CxxChannelFlow * channel
    cdef ProfileBase profile
    cdef cxx_set_wall_profile(self, ProfileBase)

cdef extern from "Profile.h":

    cdef shared_ptr[CxxProfile] CxxNewProfile "CanteraApp::newProfile" (string) except +translate_exception

    cdef cppclass CxxProfile "CanteraApp::Profile":
        void setup(double, double, double, double)
        double getZup()
        double getZdn()
        double getTup()
        double getTdn()
        double calcT(double)

    cdef cppclass CxxErfProfile "CanteraApp::ErfProfile":
        double mu()
        void setMu(double)
        double sigma()
        void setSigma(double)

    cdef cppclass CxxInterpolatedProfile "CanteraApp::InterpolatedProfile":
        void setValues(vector[double]&, vector[double]&)

    cdef cppclass CxxPolynomialProfile "CanteraApp::PolynomialProfile":
        void setValues(size_t, vector[double]&)


cdef class ProfileBase:
    cdef shared_ptr[CxxProfile] _profile
    cdef CxxProfile* profile


cdef class ConstProfile(ProfileBase):
    pass


cdef class LinearProfile(ProfileBase):
    pass


cdef class ErfProfile(ProfileBase):
    cdef CxxErfProfile* cxx


cdef class InterpolatedProfile(ProfileBase):
    cdef CxxInterpolatedProfile* cxx


cdef class PolynomialProfile(ProfileBase):
    cdef CxxPolynomialProfile* cxx
