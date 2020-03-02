// Base class to work with four vectors.
//
// All code can be easily modified to use TLorentzVector instead
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "two_body_state.hpp"

double two_body_state::get_mass(std::string name)
{
  if (name == particle1)
  {
    return m1;
  }
  else if (name == particle2)
  {
    return m2;
  }
  else
  {
    std::cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

// ---------------------------------------------------------------------------
// return the energy or momentum of a given particle

std::complex<double> two_body_state::energy(std::string name, double s)
{
  if (name == particle1)
  {
    return (s + m1*m1 - m2*m2) / (2. * sqrt(s));
  }
  else if (name == particle2)
  {
    return (s - m1*m1 + m2*m2) / (2. * sqrt(s));
  }
  else
  {
    std::cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

std::complex<double> two_body_state::momentum(std::string name, double s)
{
  if (name == particle1)
  {
    std::complex<double> E1 = energy(name, s);
    return sqrt(E1*E1 - m1*m1);
  }
  else if (name == particle2)
  {
    std::complex<double> E2 = energy(name, s);
    return -sqrt(E2*E2 - m2*m2);
  }
  else
  {
    std::cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

// ---------------------------------------------------------------------------
// The four momenta components in the x-z plane
std::complex<double> two_body_state::component(int i, std::string name, double s, double zs)
{
  switch (i)
  {
    case 0: return energy(name, s);
    case 1: momentum(name, s) * sqrt(1. - zs * zs);
    case 2: return 0.;
    case 3: return momentum(name, s) * zs;
    default: std::cout << "two_body_state: Invalid four vector component! Quiting...\n";
             exit(0);
  }
};
