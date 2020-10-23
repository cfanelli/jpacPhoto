// Sample event-generator for the gamma p -> R -> Jpsi p -> l+ l- p
// reaction. Events are generated on a flat phase space and weighted by the probabilty distribution
// of a supplied amplitude.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "monte_carlo/jpsi_mc.hpp"
#include <chrono>

// ---------------------------------------------------------------------------
// main driver loop to generate N events
void jpacPhoto::jpsi_mc::generate(int N)
{
    std::cout << "\n";
    std::cout << "Generating " << N << " events to " << filename << ". \n";


    if (exp->experiment_name != "none")
    {
        std::cout << std::left << std::setw(20) << "experiment_setup" << std::setw(5) << " = " <<  exp->experiment_name << "\n";
    };

    if (amp != NULL)
    {
        std::cout << std::left << std::setw(20) << "amplitude" << std::setw(5) << " = " << amp->identifier << "\n";
    };

    /*
    beam_energy = exp->beam_energy();
    W = W_cm(beam_energy);
    s = W * W;
    */

    auto start = std::chrono::high_resolution_clock::now();

    /*
    for (int n = 0; n < N; n++)  //
    {


        if(n%1000==0){
          std::cout<< "Evt. "<< n << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
          std::cout << "Time taken by function: " \
          << duration.count() << " [ms]" << std::endl;
          start = stop;

        }
        //moved inside event generation
        beam_energy = exp->beam_energy();
        W = W_cm(beam_energy);
        s = W * W;
        beam_pol = exp->polarization(s);
        //std::cout<<"beam_pol: "<< beam_pol <<std::endl;

        generate_event();
        if (!exp->accept_reject(current_event)) continue;
        generate_weight();
        fill_histograms();

    }
    */

    //-------------------//
    //    HIT or MISS
    //-------------------//

    int count_gen = 0;
    int attempts = 0;
    bool accept_evt = false;
    double weight_max = 1.3; //pre-computed

    while(count_gen<N){

        attempts++;
        accept_evt = false;

        //1) Beam Energy
        beam_energy = exp->beam_energy();
        W = W_cm(beam_energy);
        s = W * W;
        //2) Beam Polarization
        beam_pol = exp->polarization(s);
        beam_helicity = exp->helicity();

        //3) Generate Associated Weight
        generate_weight();
        //4) Hit or MIss
        accept_evt = exp->hit_or_miss(weight, weight_max);

        evt_weight = weight;
        hit_or_miss_activate = true;

        if(!accept_evt) continue; //repeat generation in the loop

        //-------- 1st SURVEY ---------//

        if(count_gen%100==0 && count_gen>0){
          std::cout<< "generated evts: "<< count_gen << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
          std::cout << "Time taken by function: " \
          << duration.count() << " [ms]" << std::endl;
          std::cout<<"evts accepted: "<< count_gen << ", generated: "<< attempts << std::endl;
          start = stop;
        }

        //5) Generate Event
        generate_event();

        //6) Fill Histogram
        fill_histograms();

        count_gen++;

    }

    //-------- 2nd SURVEY ---------//

    std::cout<<"=================================================="<<std::endl;
    std::cout<< std::left << std::setw(20) << "  hit or miss survey: " <<std::endl;
    std::cout<< std::left << std::setw(20) << "  n. of attempts: " << attempts << ", generated: "<< count_gen <<std::endl;
    std::cout<<"=================================================="<<std::endl;

    outfile->cd();
    kin->Write();

    if (dyn != NULL) dyn->Write();

    write_histograms();

    std::cout << "Done! \n\n";

    outfile->Close();
};

// ---------------------------------------------------------------------------
// Set up the tree structure of the output file
void jpacPhoto::jpsi_mc::set_up_kin()
{
    // tree for kinematics quantities
    kin = new TTree("kinematics", "kinematics");


    // Beam
    kin->Branch("beam_E",    &beam_E); //energy
    kin->Branch("beam_P",    &beam_P);  //polarization
    kin->Branch("beam_h",    &beam_h);  //"helicity" //PARA/PERP

    // Lepton pair 4-momentum components
    kin->Branch("ep_px",    &ep_px);
    kin->Branch("ep_py",    &ep_py);
    kin->Branch("ep_pz",    &ep_pz);
    kin->Branch("ep_E",     &ep_E);
    kin->Branch("em_px",    &em_px);
    kin->Branch("em_py",    &em_py);
    kin->Branch("em_pz",    &em_pz);
    kin->Branch("em_E",     &em_E);

    // Recoil proton 4-momenta
    kin->Branch("prec_px",  &prec_px);
    kin->Branch("prec_py",  &prec_py);
    kin->Branch("prec_pz",  &prec_pz);
    kin->Branch("prec_E",   &prec_E);

    // Angles
    kin->Branch("phi_psi",   &phi_psi);
    kin->Branch("theta_psi", &theta_psi);
    kin->Branch("phi_ee",    &phi_ee);
    kin->Branch("theta_ee",  &theta_ee);



    current_event = new event();
};

void jpacPhoto::jpsi_mc::set_up_dyn()
{
    // A second tree for dynamical quantities
    // Only need to have this if we have an amplitude supplied
    dyn = new TTree("dynamics", "dynamics");

    // Invariants
    dyn->Branch("s", &s);
    dyn->Branch("t", &t);

    // Amplitude
    dyn->Branch("DXS",    &DXS);
    dyn->Branch("Sigma",  &Sigma);
    dyn->Branch("weight", &weight);
};

// ---------------------------------------------------------------------------
// Generate an event
void jpacPhoto::jpsi_mc::generate_event()
{

    beam_E = beam_energy;
    beam_P = beam_pol;
    beam_h = beam_helicity;

    // Initial leption momenta in the jpsi decay frame
    TLorentzVector p_psi(0., 0., 0., mJpsi);
    TLorentzVector p_ep(0., 0.,  mJpsi/2., mJpsi/2.);
    TLorentzVector p_em(0., 0., -mJpsi/2., mJpsi/2.);

    // Orient in the right direction before boosting
    p_ep.RotateY(theta_ee);
    p_em.RotateY(theta_ee);
    p_ep.RotateZ(phi_ee);
    p_em.RotateZ(phi_ee);

    // Boost to the psi-proton rest frame
    TVector3 boost_psi(0., 0., sqrt(Kallen(s, mJpsi2, mPro2)) / (s + mJpsi2 - mPro2));
    p_psi.Boost(boost_psi);
    p_ep.Boost(boost_psi);
    p_em.Boost(boost_psi);

    // Rotate by the psi angle
    p_psi.RotateY(theta_psi);
    p_ep.RotateY(theta_psi);
    p_em.RotateY(theta_psi);
    p_psi.RotateZ(phi_psi);
    p_ep.RotateZ(phi_psi);
    p_em.RotateZ(phi_psi);

    TLorentzVector p_prec(-p_psi.X(), -p_psi.Y(), -p_psi.Z(), W - p_psi.E());

    // Boost everything in the lab frame
    TVector3 boost_lab(0., 0., (s - mPro2) / (s + mPro2));
    p_ep.Boost(boost_lab);
    p_em.Boost(boost_lab);
    p_prec.Boost(boost_lab);

    // Save all the values
    ep_px   = p_ep.X();    ep_py   = p_ep.Y();    ep_pz   = p_ep.Z();    ep_E   = p_ep.E();
    em_px   = p_em.X();    em_py   = p_em.Y();    em_pz   = p_em.Z();    em_E   = p_em.E();
    prec_px = p_prec.X();  prec_py = p_prec.Y();  prec_pz = p_prec.Z();  prec_E = p_prec.E();

    // Update the current event struct for accept-reject
    current_event->update(p_ep, p_em, p_prec);

    kin->Fill();
};

// ---------------------------------------------------------------------------
// Generate the weights from the amplitude supplied
void jpacPhoto::jpsi_mc::generate_weight()
{

    //------------- BEGIN moved from generate --------------//
    pgamma_E = beam_energy;
    //std::cout<<"beam_pol2: "<< beam_pol << std::endl;
    pgamma_P = beam_pol;

    // jpsi production angles
    phi_psi     = random(-M_PI, M_PI);
    theta_psi   = acos(random(-1., 1.));

    // jpsi decay angles
    phi_ee      = random(-M_PI, M_PI);
    theta_ee    = acos(random(-1., 1.));
    //------------- END   moved from generate --------------//

    if (dyn == NULL)
    {
        if (error_already_triggered == false)
        {
            std::cout << "mc: No weighting amplitude specified. Skipping weight step... \n";
            error_already_triggered = true;
        }
        return;
    }

    // Momentum transfer squared
    t = amp->kinematics->t_man(s, theta_psi);

    // pass jpsi decay angles to decay model
    amp->set_decay_angles(theta_ee, phi_ee);

    // Calculate unpolarized cross-section
    DXS = amp->differential_xsection(s, t);

    // and beam asymmetry
    Sigma = amp->beam_asymmetry_4pi(s,t);

    // Polarized cross-sections
    weight  = DXS * (1. + exp->helicity() * exp->polarization(s) * Sigma * cos(2. * phi_psi));

    dyn->Fill();
};

// ---------------------------------------------------------------------------
// All histogram related go here

void jpacPhoto::jpsi_mc::fill_histograms()
{

    float weight_h = weight;

    if(hit_or_miss_activate){
      weight_h = 1.;
    }

    th1_phi_psi->Fill(phi_psi, weight_h);
    th1_costheta_psi->Fill(cos(theta_psi), weight_h);
    th1_phi_ee->Fill(phi_ee, weight_h);
    th1_costheta_ee->Fill(cos(theta_ee), weight_h);

    // Beam energy and polarization profile will look different
    // when associated to the J/Psi selection
    th1_beamE->Fill(pgamma_E, weight_h);
    float twidth = th1_beamP->GetBinWidth(1);
    float tdiff = pgamma_E - th1_beamP->GetXaxis()->GetXmin();
    int tmpb = (int )(tdiff/twidth);
    th1_beamP->SetBinContent(tmpb,pgamma_P);

};

void jpacPhoto::jpsi_mc::write_histograms()
{
    th1_phi_psi->Write();
    th1_costheta_psi->Write();
    th1_phi_ee->Write();
    th1_costheta_ee->Write();
    th1_beamE->Write();
    th1_beamP->Write();
};
