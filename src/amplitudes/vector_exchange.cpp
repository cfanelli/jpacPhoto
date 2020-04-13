// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/vector_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> vector_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    for(int nu = 0; nu < 4; nu++)
    {
      std::complex<double> temp;
      temp = top_vertex(mu, lam_gam, lam_vec, s, zs);
      temp *= metric[mu];
      temp *= vector_propagator(mu, nu, s, zs);
      temp *= metric[nu];
      temp *= bottom_vertex(nu, lam_targ, lam_rec, s, zs);

      result += temp;
    }
  }

  return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> vector_exchange::top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs)
{
  // Contract with LeviCivita
  std::complex<double> result = 0.;
  for (int alpha = 0; alpha < 4; alpha++)
  {
    for (int beta = 0; beta < 4; beta++)
    {
      for (int gamma = 0; gamma < 4; gamma++)
      {
        std::complex<double> temp;
        temp = levi_civita(mu, alpha, beta, gamma);
        temp *= kinematics->initial.component(alpha, "beam", s, 1.);
        temp *= kinematics->eps_gamma.component(beta, lam_gam, s, 1.);
        temp *= kinematics->eps_vec.component(gamma, lam_vec, s, zs);

        result += temp;
      }
    }
  }

  // Multiply by coupling
  return result * gGamma;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> vector_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double zs)
{
  std::complex<double> vector = 0., tensor = 0.;

  // Vector coupling piece
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, -zs);
      temp *= gamma_matrices[mu][i][j];
      temp *= kinematics->target.component(j, lam_targ, s, -1.);

      vector += temp;
    }
  }

  // Tensor coupling piece
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp = 0., temp2 = 0.;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, -zs);

      for (int nu = 0; nu < 4; nu++)
      {
        temp2 += sigma(mu, nu, i, j) * exchange_momenta(nu, s, zs) / (2. * mPro);
      }
      temp *= temp2;

      temp *= kinematics->target.component(j, lam_targ, s, -1.);
    }
  }

  return gV * vector - gT * tensor;
};

// ---------------------------------------------------------------------------
std::complex<double> vector_exchange::vector_propagator(int mu, int nu, double s, double zs)
{
  std::complex<double> result;
  result =  exchange_momenta(mu, s, zs) * exchange_momenta(nu, s, zs) / mEx2;
  if (mu == nu)
  {
    result -= metric[mu];
  }

  return result / (momentum_transfer(s,zs) - mEx2);
};

// ---------------------------------------------------------------------------
// Misc other functions

// ---------------------------------------------------------------------------
// Four-momentum of the exchanged meson.
// Simply the difference of the photon and axial 4-momenta
std::complex<double> vector_exchange::exchange_momenta(int mu, double s, double zs)
{
  std::complex<double> qGamma_mu, qA_mu;
  qGamma_mu = kinematics->initial.component(mu, "beam", s, 1.);
  qA_mu = kinematics->final.component(mu, kinematics->vector_particle, s, zs);

  return (qGamma_mu - qA_mu);
};

// Mandelstam t momentum transfer
double vector_exchange::momentum_transfer(double s, double zs)
{
  double t;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp;
    temp = exchange_momenta(mu, s, zs);
    temp *= metric[mu];
    temp *= exchange_momenta(mu, s, zs);

    t += real(temp);
  }

  return t;
};
