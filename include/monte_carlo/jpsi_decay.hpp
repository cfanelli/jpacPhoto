// Extention of a vector production amplitude class to incorporate the subsequent 
// decay of j/psi -> l+ l-
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PSI_DECAY_
#define _PSI_DECAY_

#include "amplitudes/amplitude.hpp"

namespace jpacPhoto
{
    class jpsi_decay : public amplitude
    {
        public: 
        jpsi_decay(amplitude * amp)
        : amplitude(amp->kinematics, amp->identifier), production_amp(amp)
        {};
        
        // else add additional d-func of the jpsi -> l+ l- decay
        std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t);

        inline void set_decay_angles(double _theta, double _phi)
        {
            theta_ee = _theta;
            phi_ee = _phi;
        };

        private:
        // Amplitude describing the gamma p -> jpsi p production
        amplitude * production_amp;

        double theta_ee, phi_ee;
        double cached_theta_ee = 0., cached_phi_ee = 0.;

        void check_cache(double _s, double _t);
    };
};
#endif