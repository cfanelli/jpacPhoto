// Sample event-generator for the gamma p -> R -> Jpsi p -> l+ l- p
// reaction. Events are generated on a flat phase space and weighted by the probabilty distribution
// of a supplied amplitude.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TOY_MC_
#define _TOY_MC_

#include "monte_carlo/jpsi_decay.hpp"
#include "monte_carlo/event.hpp"
#include "monte_carlo/experiment_setup.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1F.h"

namespace jpacPhoto
{
    class jpsi_mc
    {
        public:
        // constructor with only filename
        jpsi_mc(std::string _filename = "mc.root")
        : filename(_filename)
        {
            rnd = new TRandom3(0);
            outfile = new TFile(filename.c_str(), "RECREATE");
            set_up_kin();
        };

        // constructor if you want to specify seed
        jpsi_mc(int seed = 0, std::string _filename = "mc.root")
        : filename(_filename)
        {
            rnd = new TRandom3(seed);
            outfile = new TFile(filename.c_str(), "RECREATE");
            set_up_kin();
        };

        // destructor
        ~jpsi_mc()
        {
            delete outfile;
            delete rnd;
            delete current_event;
            delete default_setup;
        }

        // Setting utilities
        inline void set_amplitude(amplitude * _amp)
        {
            amp = new jpsi_decay(_amp);
            set_up_dyn();
        };

        inline void set_experiment(experiment_setup * _exp)
        {
            exp = _exp;
        };

        // Produce root file with N events
        void generate(int N);

        private:

        // ROOT structures
        std::string filename;
        TFile * outfile = NULL;
        TTree * kin = NULL;
        TTree * dyn = NULL;

        // Set up the branch structure for the output root file
        void set_up_kin();
        void set_up_dyn();

        // Shorthand for generating a random number in a range
        TRandom3 * rnd = NULL;
        inline double random(double min, double max)
        {
            return min + (max - min) * rnd->Rndm();
        }

        // Generate event
        void  generate_event();  // kinematics
        void  generate_weight(); // dynamics

        // pointer with experimental beam profile and acceptances
        event * current_event;
        null_exp * default_setup = new null_exp();
        experiment_setup * exp = default_setup;

        // Pointer to preset amplitude that will be the weighting function for events
        jpsi_decay * amp = NULL;
        bool error_already_triggered = false;

        // Energies
        double beam_energy, W, s, t, beam_pol;
        double DXS, Sigma, weight;

        // Helicities of the photon and 2 x lambda of the target and recoil protons
        int lam_gamma, lam_ptarg, lam_ee, lam_prec;

        // Beam Photon
        double beam_E, beam_P, beam_h;

        // Lepton 4-momenta components
        double ep_px, ep_py, ep_pz, ep_E;
		    double em_px, em_py, em_pz, em_E;

        // Recoil and target proton components
    		double prec_px, prec_py, prec_pz, prec_E;
    		double ptarg_px, ptarg_py, ptarg_pz, ptarg_E;

        // Beam photon 4-momenta
		    double pgamma_px, pgamma_py, pgamma_pz, pgamma_E, pgamma_P;

        // Angles
        double theta_ee, phi_ee;   // of the lepton pair
        double theta_psi, phi_psi; // of the produced jpsi

        double evt_weight;
        bool hit_or_miss_activate;
        float beam_helicity;

        // Histograms
        void fill_histograms();
        void write_histograms();
        TH1F *th1_phi_psi      = new TH1F("th1_phi_psi",      "phi_psi",        20, -M_PI, M_PI);
    		TH1F *th1_costheta_psi = new TH1F("th1_costheta_psi", "costheta_psi",   20, -1., 1.);
    		TH1F *th1_phi_ee       = new TH1F("th1_phi_ee",       "phi_ee",         20, -M_PI, M_PI);
    		TH1F *th1_costheta_ee  = new TH1F("th1_costheta_ee",  "costheta_ee",    20, -1., 1.);
        TH1F *th1_beamE        = new TH1F("th1_beamE",  "beamE",  380, 8.2, 12.);// 380, 8.2, 12.
        TH1F *th1_beamP        = new TH1F("th1_beamP",  "beamP",  380, 8.2, 12.);// 380, 8.2, 12.

    };
};

#endif
