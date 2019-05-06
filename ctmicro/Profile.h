//! @file Profile.h

// This file is part of the ctmicro add-on package to Cantera.
// See LICENSE in the top-level directory

#ifndef CT_PROFILE_H
#define CT_PROFILE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/FactoryBase.h"

namespace CanteraApp {

  using namespace Cantera;

  //! Base class for wall profile definitions
  /*!
   * By default, the profile is linear
   */

  class Profile {

  public:
    Profile(); //!< Default constructor

    //! Set upstream/downstream locations and temperatures
    void setup(double Zup, double Zdn, double Tup, double Tdn);

    double getZup() { return m_Zup; } //!< Get upstream location
    double getZdn() { return m_Zdn; } //!< Get downstream location
    double getTup() { return m_Tup; } //!< Get upstream temperature
    double getTdn() { return m_Tdn; } //!< Get downstream temperature

    //! Calculate temperature at position z
    virtual double calcT(double z) {
      if (z < m_Zup) {
        return m_Tup;
      } else if (z > m_Zdn) {
        return m_Tdn;
      } else {
        return m_Tup + (z - m_Zup) / (m_Zdn - m_Zup) * (m_Tdn - m_Tup);
      }
    }

  protected:
    double m_Zup; //!< Upstream location
    double m_Zdn; //!< Downstream location
    double m_Tup; //!< Upstream temperature
    double m_Tdn; //!< Downstream temperature
    double m_Tmid; //!< Average of m_Tup and m_Tdn
    double m_DeltaT; //!< Change from m_Tup to m_Tdn
  };

  //! Derived class ErfProfile
  /*!
   * Class defines error-function shaped wall temperature profile
   */

  class ErfProfile : public Profile {

  public:
    ErfProfile(); //!< Default constructor.

    virtual double calcT(double z);

    //! Get error function mid point
    double mu() const { return m_mu; }

    //! Set error function mid point
    void setMu(const double mu) { m_mu = mu; }

    //! Get error function width
    double sigma() const { return m_sigma; }

    //! Set error function width
    void setSigma(const double sigma) { m_sigma = sigma; }

  protected:
    double m_mu; //!< Error function mid point
    double m_sigma; //!< Error function width
  };

  //! Derived class InterpolatedProfile
  /*!
   * Class defines wall temperature profile based on interpolated segments
   */

  class InterpolatedProfile : public Profile {

  public:
    InterpolatedProfile(); //!< Default constructor.

    //! Set parameters, i.e. vectors for locations and temperatures
    void setValues(vector_fp &zfixed, vector_fp &tfixed);

    virtual double calcT(double z);

  protected:
    vector_fp m_wallZfix; //!< Vector of locations
    vector_fp m_wallTfix; //!< Vector of temperatures
  };

  //! Derived class PolynomialProfile
  /*!
   * Class defines polynomial wall temperature profile. The polynomial order
   * can range from 3 to 6.
   */

  class PolynomialProfile : public Profile {

  public:
    PolynomialProfile(); //!< Default constructor.

    //! Set parameters, i.e. polynomial order and coefficients
    void setValues(size_t order, vector_fp &coeffs);

    virtual double calcT(double z);

  protected:
    size_t m_order; //! Polynomial order
    vector_fp m_coeffs; //! Vector of coefficients
  };


  //! Factory function
  /*!
   * Class facilitates definition of new phases. Derived from Cantera base
   * class 'Factory'.
   */

  class ProfileFactory : public Factory<Profile> {

  public:
    static ProfileFactory* factory() {
      std::unique_lock<std::mutex> lock(profile_mutex);
      if (!s_factory) {
        s_factory = new ProfileFactory;
      }
      return s_factory;
    }

    virtual void deleteFactory() {
      std::unique_lock<std::mutex> lock(profile_mutex);
      delete s_factory;
      s_factory = 0;
    }

  private:
    static ProfileFactory* s_factory;
    static std::mutex profile_mutex;
    ProfileFactory();
  };

  //! Create a Profile object of the specified type using the factory
  /*!
   * @param profileType the type to be created.
   */
  shared_ptr<Profile> newProfile(const std::string& model);

} // namespace CanteraMod

#endif
