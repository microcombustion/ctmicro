//! @file Profile.cpp

// This file is part of the ctmicro add-on package to Cantera.
// See LICENSE in the top-level directory

#include "Profile.h"
#include <math.h>
#include "cantera/numerics/funcs.h"
#include "cantera/base/utilities.h"

namespace CanteraApp {

  using namespace std;
  using namespace Cantera;

  ProfileFactory* ProfileFactory::s_factory = 0;
  mutex ProfileFactory::profile_mutex;

  Profile::Profile()
    : m_Zup(0.), m_Zdn(1.), m_Tup(300.), m_Tdn(300.) {}

  void Profile::setup(double Zup, double Zdn, double Tup, double Tdn) {
    m_Zup = Zup;
    m_Zdn = Zdn;
    m_Tup = Tup;
    m_Tdn = Tdn;
    m_Tmid = .5 * (m_Tup + m_Tdn);
    m_DeltaT = .5 * (m_Tdn - m_Tup);
  }

  ErfProfile::ErfProfile()
    : Profile(), m_mu(0.), m_sigma(1.) {}

  double ErfProfile::calcT(double z) {
    double xi = (z - m_mu) / m_sigma;
    return m_Tmid + m_DeltaT * erf(xi);
  }

  InterpolatedProfile::InterpolatedProfile() : Profile() {}

  void InterpolatedProfile::setValues(vector_fp &zfixed,
                                      vector_fp &tfixed) {
    m_wallZfix = zfixed;
    m_wallTfix = tfixed;
  }

  double InterpolatedProfile::calcT(double z) {
    return linearInterp(z, m_wallZfix, m_wallTfix);
  }

  PolynomialProfile::PolynomialProfile()
    : Profile(), m_order(0) {}

  void PolynomialProfile::setValues(size_t order,
                                    vector_fp &coeffs) {
    m_order = order;
    m_coeffs = coeffs;
  }

  double PolynomialProfile::calcT(double z) {
    if (m_order == 6) {
      return poly6(z, m_coeffs.data());
    } else if (m_order == 5) {
      return poly5(z, m_coeffs.data());
    } else if (m_order == 4) {
      return poly4(z, m_coeffs.data());
    } else if (m_order == 3) {
      return poly3(z, m_coeffs.data());
    }
    return -1.;
  }

  ProfileFactory::ProfileFactory()
  {
    reg("linear", []() { return new Profile(); });
    addAlias("linear", "constant");
    reg("error-function", []() { return new ErfProfile(); });
    reg("interpolated", []() { return new InterpolatedProfile(); });
    reg("polynomial", []() { return new PolynomialProfile(); });
  }

  //! Create a Profile object of the specified type
  shared_ptr<Profile> newProfile(const string& model) {
    return shared_ptr<Profile>(ProfileFactory::factory()->create(model));
  }

} // namespace CanteraApp
