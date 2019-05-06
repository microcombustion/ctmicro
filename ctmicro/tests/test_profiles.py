#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""unit tests for wall profile definitions
"""

#import logging
import numpy as np
#import os
import unittest
#import importlib

import ctmicro as cm

class Test_Profiles(unittest.TestCase):

    def assertNear(self, a, b, rtol=1e-8, atol=1e-12, msg=None):
        cmp = 2 * abs(a - b)/(abs(a) + abs(b) + 2 * atol / rtol)
        if cmp > rtol:
            message = ('AssertNear: %.14g - %.14g = %.14g\n' % (a, b, a-b) +
                       'Relative error of %10e exceeds rtol = %10e' % (cmp, rtol))
            if msg:
                message = msg + '\n' + message
            self.fail(message)

    def test_constant0(self):
        p = cm.ConstProfile()
        self.assertEqual(len(p.z_range), 2)
        self.assertEqual(len(p.T_range), 2)

    def test_constant1(self):
        p = cm.ConstProfile((-1., 1.), 300.)
        self.assertEqual(p.z_range, (-1, 1))
        self.assertEqual(p.T_range, (300, 300))
        self.assertNear(p.T(0), 300.)

    def test_linear1(self):
        p = cm.ConstProfile()
        self.assertEqual(len(p.z_range), 2)
        self.assertEqual(len(p.T_range), 2)

    def test_linear2(self):
        p = cm.LinearProfile(z_range=(-1., 1.), T_range=(300., 900.))
        self.assertEqual(p.z_range, (-1, 1))
        self.assertEqual(p.T_range, (300, 900))
        self.assertNear(p.T(0), 600.)

    def test_erf1(self):
        p = cm.ErfProfile()
        self.assertEqual(len(p.z_range), 2)
        self.assertEqual(len(p.T_range), 2)

    def test_erf2(self):
        mu = .05
        sigma = .1
        p = cm.ErfProfile(z_range=(-.1, .1), T_range=(300., 1500.), mu=mu, sigma=sigma)
        self.assertEqual(p.z_range, (-.1, .1))
        self.assertEqual(p.T_range, (300., 1500.))
        self.assertEqual(p.mu, mu)
        self.assertEqual(p.sigma, sigma)
        self.assertNear(p.T(mu), 900.)

    def test_interpolated1(self):
        p = cm.InterpolatedProfile()
        self.assertEqual(len(p.z_range), 2)
        self.assertEqual(len(p.T_range), 2)

    def test_interpolated2(self):
        p = cm.InterpolatedProfile(zfixed=[0, .3, .5, 1.], tfixed=[300., 400., 1200., 1500.])
        self.assertEqual(p.z_range, (0, 1))
        self.assertEqual(p.T_range, (300., 1500.))
        self.assertNear(p.T(.3), 400)
        self.assertNear(p.T(.4), 800)

    def test_polynomial1(self):
        p = cm.PolynomialProfile()
        self.assertEqual(len(p.z_range), 2)
        self.assertEqual(len(p.T_range), 2)

    def test_polynomial2(self):
        p = cm.PolynomialProfile(z_range=(0.,1.), coeffs=[1200., 0., 300.])
        self.assertEqual(p.z_range, (0, 1))
        self.assertNear(p.T(0), 300)
        self.assertNear(p.T(.5), 600)
        self.assertNear(p.T(1), 1500)
        self.assertEqual(p.T_range, (300., 1500.))
