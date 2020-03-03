// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/pomeron_exchange.hpp"

// ---------------------------------------------------------------------------
// Bottom vertex coupling the target and recoil proton spinors to the vector pomeron
std::complex<double> pomeron_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double zs)
{
  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
    std::complex<double> temp;
    temp = kinematics->recoil.adjoint_component(i, lam_rec, s, zs);
    temp *= gamma_matrices[mu][i][j];
    temp *= kinematics->target.component(j, lam_targ, s, 1.);

    result += temp;
    }
  }

  return result;
};

// ---------------------------------------------------------------------------
// Top vertex coupling the photon, pomeron, and vector meson.
std::complex<double> pomeron_exchange::top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs)
{
  std::complex<double> sum1 = 0., sum2 = 0.;
  for (int nu = 0; nu < 4; nu++)
  {
    std::complex<double> temp1, temp2;

    temp1 = kinematics->initial.component(nu, "beam", s, 1.);
    temp1 *= metric[nu];
    temp1 *= kinematics->eps_vec.conjugate_component(nu, lam_vec, s, zs);
    sum1 += kinematics->eps_gamma.component(mu, lam_gam, s, 1.) * temp1;

    temp2 = kinematics->eps_gamma.component(nu, lam_gam, s, 1.);
    temp2 *= metric[nu];
    temp2 *= kinematics->eps_vec.conjugate_component(nu, lam_vec, s, zs);
    sum2 += kinematics->initial.component(mu, "beam", s, 1.) * temp2;
  }

  return -sum1 + sum2;
};

// ---------------------------------------------------------------------------
// Usual Regge power law behavior, s^alpha(t) with an exponential fall from the forward direction
std::complex<double> pomeron_exchange::regge_factor(double s, double zs)
{
  if (s < kinematics->sth)
  {
    std::cout << "pomeron_exchange: Trying to evaluate below threshold! Quitting... \n";
    exit(0);
  }

  double t = kinematics->t_man(s, zs);
  double t_min = kinematics->t_man(s, 1.);

  std::complex<double> result = exp(b0 * (t - t_min)) / s;
  result *= pow(s - kinematics->sth, trajectory(t));
  result *= xi * norm;

  return result;
};

// ---------------------------------------------------------------------------
// Given a set of helicities for each particle, assemble the helicity amplitude by contracting Lorentz indicies
std::complex<double> pomeron_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp = regge_factor(s, zs);
    temp *= top_vertex(mu, lam_gam, lam_vec, s, zs);
    temp *= metric[mu];
    temp *= bottom_vertex(mu, lam_targ, lam_rec, s, zs);

    result += temp;
  }

  return result;
};