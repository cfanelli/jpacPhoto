// Extention of a vector production amplitude class to incorporate the subsequent 
// decay of j/psi -> l+ l-
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "monte_carlo/jpsi_decay.hpp"

std::complex<double> jpacPhoto::jpsi_decay::helicity_amplitude(std::vector<int> helicities, double s, double t)
{
    int lam_gam  = helicities[0];
    int lam_targ = helicities[1];
    int lam_ee   = helicities[2];
    int lam_rec  = helicities[3];

    // We sum over the intermediate jpsi helicities
    int lam_psi[3] = {-1, 0, 1};
    std::complex<double> amp_sum = 0.;
    for (int i = 0; i < 3; i++)
    {
        amp_sum  = production_amp->helicity_amplitude({lam_gam, lam_targ, lam_psi[i], lam_rec}, s, t);
        amp_sum *= wigner_d_int(1, lam_psi[i], lam_ee, theta_ee);
        amp_sum *= exp(- xi * double(lam_psi[i]) * phi_ee);
    }

    return amp_sum;
};

// make sure caching checks phi and theta dependence as well
void jpacPhoto::jpsi_decay::check_cache(double _s, double _t)
{
  // check if saved version its the one we want
  if (  (abs(cached_s - _s) < 0.00001) && 
        (abs(cached_t - _t) < 0.00001) &&
        (abs(cached_mVec2 - kinematics->mVec2) < 0.00001) &&
        (abs(cached_theta_ee - theta_ee) < 0.00001) &&
        (abs(cached_phi_ee - phi_ee) < 0.00001) 
     )
  {
    return; // do nothing
  }
  else // save a new set
  {
    for (int i = 0; i < 24; i++)
    {
      cached_helicity_amplitude[i] = helicity_amplitude(kinematics->helicities[i], _s, _t);
    }

    // update cache info
    cached_mVec2 = kinematics->mVec2; cached_s = _s; cached_t = _t;
    cached_theta_ee = theta_ee; cached_phi_ee = phi_ee;
  }

  return;
};