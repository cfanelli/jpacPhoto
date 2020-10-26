#include "TLorentzVector.h"
#include "TRandom3.h"

class event {

  public:
  event()
  {};

  //TLorentzVector ep_P4, em_P4, prec_P4, beam_P4, ptar_P4;
  TLorentzVector ep_P4, em_P4, prec_P4, jpsi_P4, ptar_P4, beam_P4, Mand_t, Mand_s;

  double Emin, Emax, tmin, tmax;
  double t_mand, s_mand;
  double ran;

  bool selection(){
    bool res = true;

    // impose cuts for custom analysis in bins of energy and t
    //cout<<"beamP4.E(): "<< beam_P4.E() << endl;
    if(beam_P4.E()<Emin || beam_P4.E()>Emax) res = false;
    if(abs(t_mand)<tmin || abs(t_mand)>tmax) res = false;

    if(beam_P4.E()>11.5) res = false; // fiducial due to polarization

    return res;
  }

  bool detector_acc(double ran){

    bool res = true;

    // analysis cuts
    if(prec_P4.P()< 0.4) res = false;
    if(prec_P4.Theta()< 2.*TMath::DegToRad()) res = false;

    if(ep_P4.P()< 0.4) res = false;
    if(ep_P4.Theta()< 2.*TMath::DegToRad()) res = false;

    if(em_P4.P()< 0.4) res = false;
    if(em_P4.Theta()< 2.*TMath::DegToRad()) res = false;


    // detector efficiency modeled as constant vs energy, see FIG. 26 jpsi_anotes18.pdf

    if(ran<0.25){ //mimicking total efficieny for p, e+, e- (~25% flat vs E)
      res = false;
      //cout<<"eff: "<< ran << endl;
    }

    return res;

  }

  // destructor
  ~event()
  {
      //delete gr;
  }

};
