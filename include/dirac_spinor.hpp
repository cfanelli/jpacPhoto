// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _SPINOR_
#define _SPINOR_

#include <cmath>
#include <complex>

#include "constants.hpp"
#include "gamma_technology.hpp"
#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// The dirac_spinor object is exactly as it sounds. :)
// Spin-1/2 particles are assumed always particle 2.
//
// That being said, particle 2 phase is already included in definition but
// angle pi for the - z direction is not!
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class dirac_spinor
  {
  private:
    // masses, energies, and momenta
    two_body_state * state;

    //Whether its an anti-particle or not
    const bool ANTI_PARTICLE = false;

    // Energy component
    std::complex<double> omega(int sign, double s);

    // angular component
    double half_angle(int lam, double theta);

  public:
    // Constructor
    dirac_spinor(two_body_state * xstate, bool if_anti = false)
    : state(xstate), ANTI_PARTICLE(if_anti)
    {};

    // Default destructor
    ~dirac_spinor(){};

    // Components
    std::complex<double> component(int i, int lambda, double s, double theta);
  	std::complex<double> adjoint_component(int i, int lambda, double s, double theta);
  };
};

#endif
