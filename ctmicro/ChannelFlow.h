//! @file ChannelFlow.h

// This file is part of the ctmicro add-on package to Cantera.
// See LICENSE in the top-level directory

// the code below is partially based on StFlow.h (Cantera 2.5.0a2)

#ifndef CT_CHANNELFLOW_H
#define CT_CHANNELFLOW_H

#include "cantera/oneD/StFlow.h"
#include "Profile.h"

namespace CanteraApp {

using namespace Cantera;

//-----------------------------------------------------------
//  Class ChannelFlow
//-----------------------------------------------------------

/**
 * A class for non-adiabatic flames with heat losses.
 */
class ChannelFlow : public StFlow {
public:
  ChannelFlow(IdealGasPhase *ph = 0, size_t nsp = 1, size_t points = 1);

  double Twall(size_t j) const { return m_profile->calcT(z(j)); }

  virtual void resize(size_t components, size_t points);

  void setWallProfile(shared_ptr<Profile> profile) {
    m_profile = profile;
  }

  // bool getTube() { return m_bl; }
  // void setTube(bool tube) { m_bl = tube; }

  //! Get Nusselt number
  double getNusselt() const { return m_nu; }

  //! Set Nusselt number
  void setNusselt(double nu) { m_nu = nu; }

  //! Get tube diameter
  double getDiameter() const { return m_dia; }

  //! Set tube diameter
  virtual void setDiameter(double dia) { m_dia = dia; }

  size_t nStates() const { return m_nv*m_points; }

  virtual void getState(double *x) {
    for (size_t i=0; i<nStates(); i++) {
      x[i] = m_state[i];
    }
  }

  virtual void getResidual(double *rsd) {
    for (size_t i=0; i<nStates(); i++) {
      rsd[i] = m_residual[i];
    }
  }

  virtual void getTMask(double *tmask) {
    for (size_t i=0; i<nStates(); i++) {
      tmask[i] = m_tmask[i];
    }
  }

protected:
  //! Evaluate the residual function. This function is called in eval
  //! after updateProperties is called.
  virtual void evalResidual(double *x, double *rsd, int *tmask, double rdt,
                            size_t jmin, size_t jmax) {
    evalFlowResidual(x, rsd, tmask, rdt, jmin, jmax);
    evalEnergyResidual(x, rsd, tmask, rdt, jmin, jmax);
    bufferData(x, rsd, tmask);
  }

  virtual void bufferData(double *x, double *rsd, int *tmask) {

    for (size_t i=0; i<m_nv*m_points; i++) {
      m_state[i] = x[i];
      m_residual[i] = rsd[i];
      m_tmask[i] = tmask[i];
    }
  }

  virtual void evalFlowResidual(double *x, double *rsd, int *tmask, double rdt,
                                size_t jmin, size_t jmax);

  virtual void evalEnergyResidual(double *x, double *rsd, int *tmask, double rdt,
                                  size_t jmin, size_t jmax);

  double m_dia; //!< Tube diameter
  double m_nu; //!< Nusselt number

  vector_fp m_q_flow;
  vector_fp m_q_cond;
  vector_fp m_q_inter;
  vector_fp m_q_reaction;
  vector_fp m_q_wall;

  vector_fp m_state;
  vector_fp m_residual;
  vector_int m_tmask;

  shared_ptr<Profile> m_profile; //!< Wall temperature profile

  //bool m_bl;
};

} // namespace CanteraApp
#endif
