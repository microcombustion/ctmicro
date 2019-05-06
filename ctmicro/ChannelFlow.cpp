//! @file ChannelFlow.cpp

// This file is part of the ctmicro add-on package to Cantera.
// See LICENSE in the top-level directory

// the code below is partially based on StFlow.cpp (Cantera 2.5.0a2)

#include <stdlib.h>
#include <time.h>

#include "ChannelFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"

#include <cstdio>

// using namespace ctml;
using namespace std;

using namespace Cantera;

namespace CanteraApp {

ChannelFlow::ChannelFlow(IdealGasPhase *ph, size_t nsp, size_t points)
  : StFlow(ph, nsp, points), m_dia(1.), m_nu(4.) {
  m_dovisc = false;
  m_type = cAxisymmetricStagnationFlow; // not a freely propagating flame
  setID("flame");
}

void ChannelFlow::evalFlowResidual(double *x, double *rsd, int *tmask, double rdt,
                                   size_t jmin, size_t jmax) {

  //----------------------------------------------------
  // evaluate the residual equations at all required
  // grid points
  //----------------------------------------------------

  for (size_t j = jmin; j <= jmax; j++) {
    //----------------------------------------------
    //         left boundary
    //----------------------------------------------

    if (j == 0) {
      // these may be modified by a boundary object

      // Continuity. This propagates information right-to-left, since
      // rho_u at point 0 is dependent on rho_u at point 1, but not on
      // mdot from the inlet.
      rsd[index(c_offset_U, 0)] = -(rho_u(x, 1) - rho_u(x, 0)) / m_dz[0] -
                                  (density(1) * V(x, 1) + density(0) * V(x, 0));

      // the inlet (or other) object connected to this one will modify
      // these equations by subtracting its values for V, T, and mdot. As
      // a result, these residual equations will force the solution
      // variables to the values for the boundary object
      rsd[index(c_offset_V, 0)] = V(x, 0);
      rsd[index(c_offset_L, 0)] = -rho_u(x, 0);

      // The default boundary condition for species is zero flux. However,
      // the boundary object may modify this.
      double sum = 0.0;
      for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x, k, 0);
        rsd[index(c_offset_Y + k, 0)] =
            -(m_flux(k, 0) + rho_u(x, 0) * Y(x, k, 0));
      }
      rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

      // set residual of poisson's equ to zero
      rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
    } else if (j == m_points - 1) {
      evalRightBoundary(x, rsd, tmask, rdt);
      // set residual of poisson's equ to zero
      rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
    } else { // interior points
      evalContinuity(j, x, rsd, tmask, rdt);
      // set residual of poisson's equ to zero
      rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

      //------------------------------------------------
      //    Radial momentum equation
      //
      //    \rho dV/dt + \rho u dV/dz + \rho V^2
      //       = d(\mu dV/dz)/dz - lambda
      //-------------------------------------------------
      rsd[index(c_offset_V, j)] =
          (shear(x, j) - lambda(x, j) - rho_u(x, j) * dVdz(x, j) -
           m_rho[j] * V(x, j) * V(x, j)) /
              m_rho[j] -
          rdt * (V(x, j) - V_prev(j));
      tmask[index(c_offset_V, j)] = true;

      //-------------------------------------------------
      //    Species equations
      //
      //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
      //   = M_k\omega_k
      //-------------------------------------------------
      getWdot(x, j);
      for (size_t k = 0; k < m_nsp; k++) {
        double convec = rho_u(x, j) * dYdz(x, k, j);
        double diffus =
            2.0 * (m_flux(k, j) - m_flux(k, j - 1)) / (z(j + 1) - z(j - 1));
        rsd[index(c_offset_Y + k, j)] =
            (m_wt[k] * (wdot(k, j)) - convec - diffus) / m_rho[j] -
            rdt * (Y(x, k, j) - Y_prev(k, j));
        tmask[index(c_offset_Y + k, j)] = true;
      }

      rsd[index(c_offset_L, j)] = lambda(x, j) - lambda(x, j - 1);
      tmask[index(c_offset_L, j)] = false;
    }
  }
}

void ChannelFlow::evalEnergyResidual(double *x, double *rsd, int *tmask, double rdt,
                                     size_t jmin, size_t jmax) {

  // coefficients
  double Tw;
  double htcoeff = m_nu / m_dia;
  htcoeff *= 4. / m_dia;

  //----------------------------------------------------
  // evaluate the residual equations at all required
  // grid points
  //----------------------------------------------------

  for (size_t j = jmin; j <= jmax; j++) {
    //----------------------------------------------
    //         left boundary
    //----------------------------------------------

    if (j == 0) {
      // these may be modified by a boundary object

      if (doEnergy(0)) {
        rsd[index(c_offset_T, 0)] = T(x, 0);
      } else {
        rsd[index(c_offset_T, 0)] = T(x, 0) - T_fixed(0);
      }
    } else if (j == m_points - 1) {
    } else { // interior points
      // evalContinuity(j, x, rsd, tmask, rdt);

      //-----------------------------------------------
      //    energy equation
      //
      //    \rho c_p dT/dt + \rho c_p u dT/dz
      //    = d(k dT/dz)/dz
      //      - sum_k(\omega_k h_k_ref)
      //      - sum_k(J_k c_p_k / M_k) dT/dz
      //-----------------------------------------------
      if (m_do_energy[j]) {
        setGas(x, j);

        const vector_fp &h_RT = m_thermo->enthalpy_RT_ref();
        const vector_fp &cp_R = m_thermo->cp_R_ref();
        double dtdzj = dTdz(x, j);

        // heat release term and interspecies fluxes
        m_q_reaction[j] = 0.;
        m_q_inter[j] = 0.;
        for (size_t k = 0; k < m_nsp; k++) {
          double flxk = 0.5 * (m_flux(k, j - 1) + m_flux(k, j));
          m_q_reaction[j] += wdot(k, j) * h_RT[k];
          m_q_inter[j] += flxk * cp_R[k] / m_wt[k];
        }
        m_q_reaction[j] *= GasConstant * T(x, j);
        m_q_inter[j] *= GasConstant * dtdzj;

        // convection and conduction
        m_q_flow[j] = m_cp[j] * rho_u(x, j) * dtdzj;
        m_q_cond[j] = divHeatFlux(x, j);

        // wall heat transfer
        Tw = m_profile->calcT(z(j));
        m_q_wall[j] = htcoeff * m_tcon[j] * (T(x, j) - Tw);

        // sum up terms
        rsd[index(c_offset_T, j)] =
          - m_q_cond[j] - m_q_flow[j] - m_q_reaction[j] - m_q_inter[j] - m_q_wall[j];

        rsd[index(c_offset_T, j)] /= (m_rho[j] * m_cp[j]);
        rsd[index(c_offset_T, j)] -= rdt * (T(x, j) - T_prev(j));
        tmask[index(c_offset_T, j)] = true;
      } else {
        // residual equations if the energy equation is disabled
        rsd[index(c_offset_T, j)] = T(x, j) - m_profile->calcT(z(j)); //T_fixed(j);
        tmask[index(c_offset_T, j)] = false;
      }

    }
  }
}

void ChannelFlow::resize(size_t components, size_t points) {
  // call parent
  StFlow::resize(components, points);

  m_q_flow.resize(m_points, 0.);
  m_q_cond.resize(m_points, 0.);
  m_q_inter.resize(m_points, 0.);
  m_q_reaction.resize(m_points, 0.);
  m_q_wall.resize(m_points, 0.);

  m_state.resize(nStates(), 0.);
  m_residual.resize(nStates(), 0.);
  m_tmask.resize(nStates(), 0);
}

} // namespace CanteraApp
