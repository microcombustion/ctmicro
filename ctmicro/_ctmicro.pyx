# This file is part of the ctmicro add-on package to Cantera.
# See LICENSE in the top-level directory

import numpy as np
cimport numpy as np

import cantera as ct
cimport cantera as ct

from cython.operator cimport dereference as deref

from cantera import interrupts

__all__ = []

cdef string stringify(x) except *:
    """ Converts Python strings to std::string. """
    if isinstance(x, bytes):
        return string(<bytes>x)
    else:
        tmp = bytes(x.encode())
        return string(tmp)

include "profile.pyx"
include "channel.pyx"
