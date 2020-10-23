// Class containing parameters and methods for experiment dependent quantities.
// E.g. beam energy spectrum or event accept/reject
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _EXP_SETUP_
#define _EXP_SETUP_

#include "monte_carlo/event.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLorentzVector.h"


namespace jpacPhoto
{
    class experiment_setup
    {
        public:
        experiment_setup()
        {};

        experiment_setup(std::string name)
        : experiment_name(name)
        {};

        // String identifier
        std::string experiment_name;

        // Acceptance check / cuts for an event
        virtual bool accept_reject(event * _event) = 0;

        // Acceptance check / cuts for an event
        virtual bool hit_or_miss(double weight, double weight_max) = 0;

        // Returns Lab beam energy in GeV
        virtual double beam_energy() = 0;

        // Degree of polarization
        virtual double polarization(double s) = 0;

        // Helicity
        virtual float helicity() = 0;

    };

    // Simplest example
    class null_exp : public experiment_setup
    {
        public:

        null_exp()
        : experiment_setup("none")
        {};

        null_exp(double e, double P)
        : experiment_setup("test"), constant_energy(e), constant_polarization(P)
        {};

        double constant_energy = 15.;
        double constant_polarization = 0.4;

        inline bool accept_reject(event * _event)
        {
            return true; // always accept
        };

        inline bool hit_or_miss(double weight, double weight_max){
          return true;
        }

        inline double polarization(double s)
        {
            return constant_polarization; // constant polarization
        };

        inline double beam_energy()
        {
            return constant_energy; // fixed, mono-energy beam
        };

        inline float helicity()
        {
          return 1.;
        };
    };

    // GlueX class
    class gluex_setup : public experiment_setup
    {
    public:

      // constructor with parameters
      gluex_setup()
      : experiment_setup("gluex")
      {

        //-------------------------------------------------//
        polfile = TFile::Open("./inputfiles/beamP.root");
        polhist = (TH1*)polfile->Get("h");
        //-------------------------------------------------//
        enefile = TFile::Open("./inputfiles/beamE.root");
        enehist = (TH1*)enefile->Get("h");
        //-------------------------------------------------//

        hc = enehist->GetCumulative();
        int n = hc->GetNbinsX();
        double x_c[n], y_x[n];
        max_c = hc->GetMaximum();

        for(int i=0;i<n;i++) {
          x_c[i] = hc->GetBinContent(i);
          y_x[i] = hc->GetBinCenter(i);
        }
        gr = new TGraph(n,x_c,y_x);

      };


      // destructor
      ~gluex_setup()
      {
          delete polfile;
          delete polhist;
          delete enefile;
          delete enehist;
          delete hc;
          delete gr;
      }

      bool accept_reject(event * _event)
      {

          bool pass_selection = true;

          p_ep = _event->p_ep;
          p_em = _event->p_em;
          p_prec = _event->p_prec;


          //SMEARING...

          //SELECTION...

          TLorentzVector MM2;
          TLorentzVector Minv = p_ep+p_em;

          if(p_prec.P()<0.4) pass_selection = false;
          if(p_ep.P()<0.4) pass_selection = false;
          if(p_em.P()<0.4) pass_selection = false;

          if(p_prec.Theta()*TMath::RadToDeg()<2.) pass_selection = false;
          if(p_ep.Theta()*TMath::RadToDeg()<2.) pass_selection = false;
          if(p_em.Theta()*TMath::RadToDeg()<2.) pass_selection = false;

          //this will be continued as analysis on the ROOT tree generated with hit or miss


          return pass_selection; // always accept
      };

      bool hit_or_miss(double weight, double weight_max){

        bool hit_event = false;

        hitmiss_pr = r3.Uniform(0,weight_max);


        if(hitmiss_pr<weight/weight_max) hit_event = true; // pass

        return hit_event;
      }


      inline double polarization(double s)
      {

        Egamma = E_beam(sqrt(s));

        bin_width = polhist->GetBinWidth(1);
        Emin = polhist->GetXaxis()->GetXmin();
        Emax = polhist->GetXaxis()->GetXmax();

        tmpbin = (int) ((Egamma-Emin)/bin_width);

        r_Pol = polhist->GetBinContent(tmpbin);


        if(r_Pol<0. || Egamma<Emin || Egamma>Emax) r_Pol = 0.;



        return r_Pol;


      };

      inline float helicity(){
        hel;
        float tmp = r3.Uniform(-1.,1.);
        if(tmp>=0.) hel =  1.;
        else hel = -1.;

        return hel;


      }


      inline double beam_energy()
      {


        sampledE = -1;
        while(sampledE<8.3 || sampledE>12.){

          tmp_extr = r3.Uniform(0,max_c);
          sampledE = gr->Eval(tmp_extr);

        }

        // assert that the energy is larger than J/Psi threshold

        return sampledE;
      };


    private:

    // ROOT structures
      TFile * polfile = NULL;
      TH1 * polhist = NULL;
      TFile * enefile = NULL;
      TH1 * enehist = NULL;
      TH1 * hc = NULL; //cumulative
      TH1 * hcP = NULL; //cumulative

      TGraph * gr = NULL;
      TGraph * grP = NULL;

      double bin_width, Emin, Emax, Egamma;
      double r_Pol, r_Ene;
      int tmpbin;
      TRandom3 r3;
      float hitmiss_pr;
      int max_c, tmp_extr;
      int max_cP, tmp_extrP;

      float hel;

      double sampledE, sampleP;
      TLorentzVector p_ep, p_em, p_prec;

    };

};
#endif
