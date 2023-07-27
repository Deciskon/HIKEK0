#include "KPIPI_analysis.hh"

#include "Event.hh"
#include "Persistency.hh"
#include "functions.hh"
#include "KaonDecayConstants.hh"
#include "PnnAnalysisTools.hh"
#include "GeometricAcceptance.hh"

#include <iostream>
#include <stdlib.h>
using namespace NA62Analysis;



KPIPI_analysis::KPIPI_analysis(Core::BaseAnalysis *ba): Analyzer(ba, "KPIPI_analysis") {

    RequestTree("Spectrometer",new TRecoSpectrometerEvent);

    AddParam("IsKSBeam", &fIsKSBeam, true);  // if true histos for KS beam mode
    AddParam("IsKLBeam", &fIsKLBeam, false); // if true, histos for KL beam mode
    AddParam("zTarget", &fTargetZ, 113000); // [mm] z position of the target (default 113m from T10 of NA62)

    fTool = PnnAnalysisTools::GetInstance();
}

void KPIPI_analysis::InitOutput() {

}


//TEST

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTANT NOTES $ INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. 'vtx' values (extracted using 'particle->GetEndPos()') are the points where the particle (mostly Kaons in this case) decayed

2. MK0 is Mass of Kaon 0 (Short or Long depends on which mode the analyzer is being operated in)

3. 'gamma' of K is the BOOST   (boost = gamma = Energy/Mass = E/m)

4. Thetak is the angle of the KAON after it has been produced from the proton-on-target collision 

5. 'geom' plots take place after the 1st cut (only particles detected at r>=200mm radially)
   'final' plots take place after the 2nd cut (kaons that decayed before the first straw)

6. Everything is in atomic units
*/



void KPIPI_analysis::InitHist() {

    BookHisto(new TH1I("Ntracks", "Ntracks;No.of Tracks;Count", 20, 0, 20));
    if(fIsKSBeam){

        BookHisto(new TH1I("Prod_KS_zvtx", "True Kaons Decays vs Z;Z_{true vtx} [m];Kaon Decays", 240, 100, 220));
        BookHisto(new TH2F("Prod_KS_xyvtx", "True XY Plane Kaon Decays;X_{true vtx} [mm]; Y_{true vtx} [mm]", 1500 , -1500, 1500,1500 , -1500, 1500));
        BookHisto(new TH2F("Prod_KS_xzvtx", "True ZX Plane Kaon Decays;Z_{true vtx} [m]; X_{true vtx} [mm]", 240, 100, 220, 1500 , -1500, 1500));
        BookHisto(new TH2F("Prod_KS_yzvtx", "True ZY Plane Kaon Decays;Z_{true vtx} [m]; Y_{true vtx} [mm]", 240, 100, 220, 1500 , -1500, 1500));
        BookHisto(new TH1F("Prod_PKS", "True Momentum Distribution;P_{K_{S}} [GeV/c];Count", 800, 0, 400));
        BookHisto(new TH2F("Prod_KS_PKS_vs_zvtx", "True Momentum vs Z;Z_{true vtx} [m]; P_{K_{S}} [GeV/c]",240, 100, 220, 800, 0, 400));
	BookHisto(new TH2F("Prod_tauks_true_vs_Pk_true", "True #tau_{K_{S}} vs True Momentum;P_{K true} [GeV/c];t [#tau_{K_{S} true}]", 800, 0, 400, 200, 0, 50));
    }

    if (fIsKLBeam){

        BookHisto(new TH1I("Prod_KL_zvtx", "True Kaon Decays vs Z;Z_{true vtx} [m];Kaon Decays", 240, 100, 220));
        BookHisto(new TH2F("Prod_KL_xyvtx", "True XY Plane Kaon Decays;X_{true vtx} [mm]; Y_{true vtx} [mm]", 1500 , -1500, 1500,1500 , -1500, 1500));
        BookHisto(new TH2F("Prod_KL_xzvtx", "True ZX Plane Kaon Decays;Z_{true vtx} [m]; X_{true vtx} [mm]", 240, 100, 220, 1500 , -1500, 1500));
        BookHisto(new TH2F("Prod_KL_yzvtx", "True ZY Plane Kaon Decays;Z_{true vtx} [m]; Y_{true vtx} [mm]", 240, 100, 220, 1500 , -1500, 1500));
        BookHisto(new TH1F("Prod_PKL", "True Momentum Distribution;P_{K_{L}} [GeV/c];Count", 800, 0, 400));
        BookHisto(new TH2F("Prod_KL_PKL_vs_zvtx", "True Momentum vs Z;Z_{true vtx} [m]; P_{K_{L}} [GeV/c]",240, 100, 220, 800, 0, 400));
    }

    BookHisto(new TH1I("PDG_all", "True PDG Distribution;PDG;Count", 1000, -500, 500));

    BookHisto(new TH1F("Pdecay_correction_factor", "Kaon Boost vs Decayed Fraction;#gamma_{K true};N_{decayed K}/N_{produced K}", 810,0,810));
    BookHisto(new TH1F("PK_reweighted", "Reweighted Momentum (Momentum vs Inverse of Decayed Fraction);P_{K} [GeV/c];N_{produced K}/N_{decayed K}", 800, 0, 400));

    BookHisto(new TH1F("Reco_Pplus", "Reconstructed #pi^{+} Momentum Distribution;P_{#pi^{+}} [GeV/c];Count", 600, 0, 300));
    BookHisto(new TH1F("Reco_Pminus", "Reconstructed #pi^{-} Momentum Distribution;P_{#pi^{-}} [GeV/c];Count", 600, 0, 300));
    BookHisto(new TH1F("Reco_Pk", "Reconstructed Momentum Distribution;P_{K} [GeV/c];Count", 800, 0, 400));
    BookHisto(new TH1F("Reco_Mk", "Reconstructed Kaon Mass Distribution;M_{K} [GeV/c^{2}];Count", 500, 0, 1));
    BookHisto(new TH2F("Reco_piplus_xy_straw1", "Reconstructed #pi^{+} XY Plane;X_{STRAW1} [mm]; Y_{STRAW1} [mm]", 1500 , -1500, 1500,1500 , -1500, 1500));
    BookHisto(new TH2F("Reco_piminus_xy_straw1", "Reconstructed #pi^{-} XY Plane;X_{STRAW1} [mm]; Y_{STRAW1} [mm]", 1500 , -1500, 1500,1500 , -1500, 1500));
    BookHisto(new TH2F("Reco_piplus_only_xy_straw1", "Reconstructed #pi^{+} ONLY XY Plane;X_{STRAW1} [mm]; Y_{STRAW1} [mm]", 1500 , -1500, 1500,1500 , -1500, 1500));
    BookHisto(new TH2F("Reco_piminus_only_xy_straw1", "Reconstructed #pi^{-} ONLY XY Plane;X_{STRAW1} [mm]; Y_{STRAW1} [mm]", 1500 , -1500, 1500,1500 , -1500, 1500));
    BookHisto(new TH2F("Reco_plus_vs_minus_Rstraw1", "Reconstructed Radial Distance of- #pi^{-} vs #pi^{+};R_{STRAW1 q=-1} [mm]; R_{STRAW1 q=1} [mm]", 1500 , 0, 1500,1500 , 0, 1500));
    BookHisto(new TH2F("Reco_geom_cda_vs_zvtx", "Reconstructed Closest Distance of Approach (CDA) vs Z; Z_{reco vtx} [m];CDA [mm]", 240, 100, 220, 200,0, 200));
    BookHisto(new TH2F("Reco_geom_pk_vs_zvtx", "Reconstructed Momentum vs Z;Z_{reco vtx} [m];P_{K} [GeV/c]", 240, 100, 220, 800,0, 400));
    BookHisto(new TH1F("Reco_geom_dzvtx", "Distribution of Difference between True Z and Reco Z; Z_{reco vtx} - Z_{true vtx} [m];Count", 500, -50, 50));
    BookHisto(new TH1F("Reco_geom_Pplus", "Reconstructed #pi^{+} Momentum Distribution;P_{#pi^{+}} [GeV/c];Count", 600, 0, 300));
    BookHisto(new TH1F("Reco_geom_Pminus", "Reconstructed #pi^{-} Momentum Distribution;P_{#pi^{-}} [GeV/c];Count", 600, 0, 300));
    BookHisto(new TH1F("Reco_geom_Pk", "Reconstructed Momentum Distribution;P_{K} [GeV/c];Count", 800, 0, 400));
    BookHisto(new TH1F("Reco_geom_Mk", "Reconstructed Kaon Mass Distribution;M_{K} [GeV/c^{2}];Count", 500, 0, 1));
    BookHisto(new TH1I("Reco_geom_PDG", "Reconstructed PDG Distribution;PDG;Count", 1000, -500, 500));

    BookHisto(new TH2F("Reco_geom_dPk_vs_Pktrue", "True Kaon Momentum vs Difference in Reco and True Kaon Momentum;P_{K true} [GeV/c];P_{K reco} - P_{K true} [GeV/c]", 800, 0, 400, 1000, -10,10));
    BookHisto(new TH2F("Reco_geom_dThetak_vs_Pktrue", "True Kaon Momentum vs Difference in Reco and True Kaon Angle;P_{K true} [GeV/c];#Theta_{K reco} - #Theta_{K true} [rad]", 800, 0, 400, 1000, -0.001,0.001));

    BookHisto(new TH2F("Reco_geom_dPplus_vs_Pplus_true", "True #pi^{+} Momentum vs Difference in Reco and True #pi^{+} Momentum;P_{#pi^{+} true} [GeV/c];P_{#pi^*{+} reco} - P_{#pi^{+} true} [GeV/c]", 800, 0, 400, 800, -50,50));
    BookHisto(new TH2F("Reco_geom_dPminus_vs_Pminus_true", "True #pi^{-} Momentum vs Difference in Reco and True #pi^{-} Momentum;P_{#pi^{-} true} [GeV/c];P_{#pi^*{-} reco} - P_{#pi^{-} true} [GeV/c]", 800, 0, 400, 1000, -10,10));
    BookHisto(new TH2F("Reco_geom_dTheta_minus_vs_Theta_minus_true", "True #pi^{-} Angle vs Difference in Reco and True #pi^{-} Angle;#Theta_{#pi^{-} true} [rad];#Theta_{#pi^{-} reco} - #Theta_{#pi^{-} true} [rad]", 1000, 0,0.02, 1000, -0.001,0.001));
    BookHisto(new TH2F("Reco_geom_dTheta_minus_vs_Pminus_true", "True #pi^{-} Momentum vs Difference in Reco and True #pi^{-} Angle;P_{#pi^{-} true} [GeV/c];#Theta_{#pi^{-} reco} - #Theta_{#pi^{-} true} [rad]", 800, 0, 400, 1000, -0.001,0.001));
    BookHisto(new TH2F("Reco_geom_dTheta_plus_vs_Theta_plus_true", "True #pi^{+} Angle vs Difference in Reco and True #pi^{+} Angle;#Theta_{#pi^{+} true} [rad];#Theta_{#pi^{+} reco} - #Theta_{#pi^{+} true} [rad]", 1000, 0,0.02, 1000, -0.001,0.001));
    BookHisto(new TH2F("Reco_geom_dTheta_plus_vs_Pplus_true", "True #pi^{+} Momentum vs Difference in Reco and True #pi^{+} Angle;P_{#pi^{+} true} [GeV/c];#Theta_{#pi^{+} reco} - #Theta_{#pi^{+} true} [rad]", 800, 0, 400, 1000, -0.001,0.001));

    BookHisto(new TH1I("Reco_geom_NKineParts", "Number of KineParts Distribution;KineParts;Count", 10, -0.5, 9.5));
    BookHisto(new TH1I("Reco_geom_Npidecays_before_straw4", "Number of Decays Before Each Straw Distribution;Decays Before Straw(x+1);Count", 3, -0.5, 2.5));

    BookHisto(new TH2F("Reco_geom_gammaK_vs_Pk", "Reco Kaon Momentum vs Reco Kaon Boost;P_{K reco} [GeV/c];#gamma_{K reco}", 800, 0, 400, 800,0,400));
    BookHisto(new TH2F("Reco_geom_gammaK_true_vs_Pk_true", "True Kaon Momentum vs True Kaon Boost;P_{K true} [GeV/c];#gamma_{K true}", 800, 0, 400, 800,0,400));
    BookHisto(new TH2F("Reco_geom_dgammaK_vs_Pk_true", "True Kaon Momentum vs Difference in Reco and True Kaon Boost;P_{K true} [GeV/c];#gamma_{K reco} - #gamma_{K true}", 800,0, 400, 200, -20, 20));

    BookHisto(new TH1F("Reco_final_Pplus", "Reconstructed #pi^{+} Momentum Distribution;P_{#pi^{+}} [GeV/c];Count", 600, 0, 300));
    BookHisto(new TH1F("Reco_final_Pminus", "Reconstructed #pi^{-} Momentum Distribution;P_{#pi^{-}} [GeV/c];Count", 600, 0, 300));
    BookHisto(new TH1F("Reco_final_Pk", "Reconstructed Momentum Distribution;P_{K} [GeV/c];Count", 800, 0, 400));
    BookHisto(new TH1F("Reco_final_Mk", "Reconstructed Kaon Mass Distribution;M_{K} [GeV/c^{2}];Count", 500, 0.4, 0.6));
    BookHisto(new TH2F("Reco_final_cda_vs_zvtx", "Closest Distance of Approach (CDA) vs Z;Z_{reco vtx} [m];CDA [mm]", 240, 100, 220, 200,0, 200));
    BookHisto(new TH2F("Reco_final_dPk_vs_Pktrue", "True Kaon Momentum vs Difference in Reco and True Kaon Momentum;P_{K true} [GeV/c];P_{K reco} - P_{K true} [GeV/c]", 800, 0, 400, 1000, -10,10));
    BookHisto(new TH2F("Reco_final_dThetak_vs_Pktrue", "True Kaon Momentum vs Difference in Reco and True Kaon Angle;P_{K true} [GeV/c];#Theta_{K reco} - #Theta_{K true} [rad]", 800, 0, 400, 1000, -0.001,0.001));

    BookHisto(new TH2F("Reco_final_tauks_vs_Pk", "Reconstructed #tau_{K_{S}} vs Reconstructed Momentum;P_{K reco} [GeV/c];t [#tau_{K_{S} reco}]", 800, 0, 400, 200, 0, 50));
    BookHisto(new TH2F("Reco_final_tauks_vs_zvtx", "Reconstructed #tau_{K_{S}} vs Reconstructed Z;Z_{vtx reco} [m];t [#tau_{K_{S} reco}]", 240, 100, 220, 200, 0, 50));

    BookHisto(new TH2F("Reco_final_dtauks_vs_Pk_true", "Difference in Reco and True #tau_{K_{S}} vs True Momentum;P_{K true} [GeV/c];t [#tau_{K_{S} reco} - #tau_{K_{S} true}]", 800, 0, 400, 400,-10,10));
    BookHisto(new TH2F("Reco_final_dtauks_vs_zvtx_true", "Difference in Reco and True #tau_{K_{S}} vs True Z;Z_{vtx true} [m];t [#tau_{K_{S} reco} - #tau_{K_{S} true}]", 240, 100, 220, 400, -10, 10));

    BookHisto(new TH2F("Reco_final_dgammaK_vs_Pk_true", "True Kaon Momentum vs Difference in Reco and True Kaon Boost;P_{K true} [GeV/c];#gamma_{K reco} - #gamma_{K true}", 800,0, 400, 200, -20, 20));

    BookHisto(new TH2F("Reco_final_Mk_vs_Pk", "Reconstructed Momentum vs Mass of Kaons;P_{K reco} [GeV/c]; M_{K reco} [GeV/c^{2}]", 800, 0, 400, 500, 0.4, 0.6));
    BookHisto(new TH2F("Reco_final_dMk_vs_Pk_true", "True Kaon Momentum vs Difference in Reco and True Mass of Kaons;P_{K true} [GeV/c];M_{K reco} - M_{K true} [GeV/c^{2}]", 800, 0, 400, 500, -0.05, 0.05));
    BookHisto(new TH2F("Reco_final_dMk_vs_Thetak_true", "True Kaon Angle vs Difference in Reco and True Mass of Kaons;#Theta_{K true} [rad];M_{K reco} - M_{K true} [GeV/c^{2}]",200, 0,0.001, 500, -0.05, 0.05));

    

    //tau_ks DEBUGGING
    BookHisto(new TH2F("Reco_final_tauks_vs_tauks_true", "Reconstructed #tau_{K_{S}} vs True #tau_{K_{S}};t_true [#tau_{K_{S} true}];t_reco [#tau_{K_{S} reco}]", 200, 0, 50, 200, 0, 50));
    BookHisto(new TH2F("Reco_final_zvtx_vs_zvtx_true", "Reconstructed Z vs True Z;Z_{vtx true} [m];Z_{vtx reco} [m]", 240, 100, 220, 240, 100, 220));
    BookHisto(new TH2F("Reco_geom_posstraw1_plus_vs_minus", "Posstraw1: pi+ vs pi-;posstraw1_plus;posstraw1_minus", 1500, -1500, 1500, 1500, -1500, 1500));
    BookHisto(new TH2F("Reco_geom_fTp_vs_fTm", "fTp vs fTm;fTp;fTm", 200, -1, 1, 200, -1, 1));
    BookHisto(new TH1I("Reco_geom_zvtx", "Kaons Decays vs Reco Z;Z_{reco vtx} [m];Kaon Decays", 240, 100, 220));
}

void KPIPI_analysis::StartOfRunUser() {

}

void KPIPI_analysis::Clear() {

    fTp = TLorentzVector(0.0,0.0,0.0,0.0);
    fTm = TLorentzVector(0.0,0.0,0.0,0.0);
    fTwoTrack = TLorentzVector(0.0,0.0,0.0,0.0);
    fTrueVertex = TLorentzVector(0.0,0.0,0.0,0.0);
    fRecoVertex = TVector3(0.0,0.0,0.0);
    fSpectrometerEvent = NULL;
    fSplus = NULL;
    fSminus = NULL;
    fCDA = 99999.;
    fBeamParticle = NULL;
    fctau_ks = 0;
    fctau_kl = 0;
}

void KPIPI_analysis::StartOfBurstUser() {

}

void KPIPI_analysis::ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType) {

}

void KPIPI_analysis::Process(int iEvent) {

    Clear();
    //fctau_ks = 0.026844; //[m]
    // fctau_kl = 15.34; //[m]
    fctau_ks = 26.844; //[mm]
    fctau_kl = 15340.0; //[mm]

    fSpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");

    FillHisto("Ntracks", fSpectrometerEvent->GetNCandidates());

    Event *evt = GetMCEvent();
    TIter partit_all(evt->GetKineParts());
    while(partit_all.Next()) {  // itterate through ALL particles
        KinePart *particle = static_cast<KinePart *>(*partit_all);
        if(particle->GetPDGcode() == 310) {  // KS
            fTrueVertex = particle->GetEndPos();
            fBeamParticle = particle;
            FillHisto("Prod_KS_zvtx", fTrueVertex.Z()/1e3);
            FillHisto("Prod_KS_xzvtx", particle->GetEndPos().Z()/1e3, particle->GetEndPos().X());
            FillHisto("Prod_KS_yzvtx", particle->GetEndPos().Z()/1e3, particle->GetEndPos().Y());
            FillHisto("Prod_KS_xyvtx", particle->GetEndPos().X(), particle->GetEndPos().Y());
            FillHisto("Prod_PKS", particle->GetInitial4Momentum().P()/1e3);
            FillHisto("Prod_KS_PKS_vs_zvtx", fTrueVertex.Z()/1e3, particle->GetInitial4Momentum().P()/1e3);
        }
        if(particle->GetPDGcode() == 130) {  // KL
            fTrueVertex = particle->GetEndPos();
            fBeamParticle = particle;
            FillHisto("Prod_KL_zvtx", fTrueVertex.Z()/1e3);
            FillHisto("Prod_KL_xzvtx", particle->GetEndPos().Z()/1e3, particle->GetEndPos().X());
            FillHisto("Prod_KL_yzvtx", particle->GetEndPos().Z()/1e3, particle->GetEndPos().Y());
            FillHisto("Prod_KL_xyvtx", particle->GetEndPos().X(), particle->GetEndPos().Y());
            FillHisto("Prod_PKL", particle->GetInitial4Momentum().P()/1e3);
            FillHisto("Prod_KL_PKL_vs_zvtx", fTrueVertex.Z()/1e3, particle->GetInitial4Momentum().P()/1e3);
        }
        FillHisto("PDG_all", particle->GetPDGcode());
        if(particle->GetPDGcode() == 211) {  // pi+
            //fTp = particle->GetFinal4Momentum();
        }
        if(particle->GetPDGcode() == -211) {  // pi-
           //fTm = particle->GetFinal4Momentum();
        }

        //std::cout << "PDG = " << particle->GetPDGcode() << " Endpos = " << particle->GetEndPos().Z() << " P =  " << particle->GetInitial4Momentum().Px() << "  " << particle->GetInitial4Momentum().Py() << "   " <<  particle->GetInitial4Momentum().Pz() << "   " << particle->GetInitial4Momentum().P()  << " M = " << particle->GetFinal4Momentum().M() << std::endl;
    }

    // if(fTrueVertex.Z() < 120000 || fTrueVertex.Z() > 198000) return; //FV between 120m and 198m (target is at 113m and STRAW at 213 m)
    if(fTrueVertex.Z() < fTargetZ || fTrueVertex.Z() > 198000) return; //FV between 120m and 198m (target is at 113m and STRAW at 213 m)
    //if(fTrueVertex.Z() >= 120000) return; //used to compute the fraction of particles lost due to collimation

    double beta_K_true    = (fBeamParticle->GetInitial4Momentum().P()/1e3)/TMath::Sqrt( (fBeamParticle->GetInitial4Momentum().P()/1e3)*(fBeamParticle->GetInitial4Momentum().P()/1e3) + MK0/1000.*MK0/1000.);
    double gamma_K_true   = 1/TMath::Sqrt(1 - beta_K_true*beta_K_true);
    double tau_ks_true    = (fTrueVertex.Z() - fTargetZ)/(fctau_ks*gamma_K_true);

    //Weight used to compute the total amount of KL that correspond to the generated KL decays
    //fdecay : Probability that the kaon that is generated at the primary target would have decayed within the defined fiducial volume
    double fdecay = -1;
    if(fIsKLBeam)
        fdecay = 1 - TMath::Exp(- (198000 - fTargetZ)/(fctau_kl*gamma_K_true));
        //fdecay = 1 - TMath::Exp(- (fTrueVertex.Z() - fTargetZ)/(fctau_kl*gamma_K_true));
    if(fIsKSBeam)
        fdecay = 1 - TMath::Exp(- (198000 - fTargetZ)/(fctau_ks*gamma_K_true));

    //std::cout << "KS - " <<fIsKSBeam << " KL - " << fIsKLBeam << " gamma_K = " << gamma_K_true <<  " zvtx = " << fTrueVertex.Z() << " fdecay = " << fdecay <<  std::endl;
    FillHisto("Pdecay_correction_factor",gamma_K_true, fdecay );
    FillHisto("PK_reweighted", fBeamParticle->GetInitial4Momentum().P()/1e3, 1/fdecay);
    if(fIsKSBeam) {  // KS
	FillHisto("Prod_tauks_true_vs_Pk_true", fBeamParticle->GetInitial4Momentum().P()/1e3, tau_ks_true);
    }

    if(fSpectrometerEvent->GetNCandidates()== 0) return;
    //if(fSpectrometerEvent->GetNCandidates()!= 2) return;
    double qvertex = 0;
    for (int iT(0); iT < fSpectrometerEvent->GetNCandidates();iT++){
        TRecoSpectrometerCandidate *pT = static_cast<TRecoSpectrometerCandidate *>(fSpectrometerEvent->GetCandidate(iT));
        //std::cout << "Q = " << pT->GetCharge() << std::endl;
        qvertex += pT->GetCharge();
        if(pT->GetCharge()==1){
            FillHisto("Reco_Pplus", pT->GetMomentum()/1000.0);

            fTp = GetTrack(pT,MPI/1000.0);
            fSplus  = pT;
        }
        if(pT->GetCharge()==-1){
            FillHisto("Reco_Pminus", pT->GetMomentum()/1000.0);

            fTm = GetTrack(pT,MPI/1000.0);
            fSminus  = pT;
            //std::cout << MPI << "MK0 = " << MK0 << std::endl;
        }
    }


   if(qvertex == 1)
   {
    TVector3 posstraw1_plus_only = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(0));
    /*
    TVector3 posstraw2_plus_only = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(1));
    TVector3 posstraw3_plus_only = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(2));
    TVector3 posstraw4_plus_only = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(3));
    */
    double rstraw1_plus_only  = TMath::Sqrt(TMath::Power(posstraw1_plus_only.X(),2)  + TMath::Power(posstraw1_plus_only.Y(),2) );
    if(rstraw1_plus_only >= 200 ) //GEOMETRIC CUT
    FillHisto("Reco_piplus_only_xy_straw1" , posstraw1_plus_only.X(), posstraw1_plus_only.Y());
   }

   if(qvertex == -1)
   {
    TVector3 posstraw1_minus_only = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(0));
    /*
    TVector3 posstraw2_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(1));
    TVector3 posstraw3_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(2));
    TVector3 posstraw4_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(3));
    */
    double rstraw1_minus_only  = TMath::Sqrt(TMath::Power(posstraw1_minus_only.X(),2)  + TMath::Power(posstraw1_minus_only.Y(),2) );
    if(rstraw1_minus_only >= 200 ) //GEOMETRIC CUT
    FillHisto("Reco_piminus_only_xy_straw1", posstraw1_minus_only.X(), posstraw1_minus_only.Y());
   }


   if(qvertex == 0) 
   {

    TVector3 posstraw1_plus = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(0));
    TVector3 posstraw2_plus = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(1));
    TVector3 posstraw3_plus = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(2));
    TVector3 posstraw4_plus = fTool->GetPositionAtZ(fSplus, GeometricAcceptance::GetInstance()->GetZStraw(3));



    TVector3 posstraw1_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(0));
    TVector3 posstraw2_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(1));
    TVector3 posstraw3_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(2));
    TVector3 posstraw4_minus = fTool->GetPositionAtZ(fSminus, GeometricAcceptance::GetInstance()->GetZStraw(3));

    double rstraw1_plus  = TMath::Sqrt(TMath::Power(posstraw1_plus.X(),2)  + TMath::Power(posstraw1_plus.Y(),2) );
    double rstraw1_minus = TMath::Sqrt(TMath::Power(posstraw1_minus.X(),2) + TMath::Power(posstraw1_minus.Y(),2) );
    fTwoTrack.SetPxPyPzE(fTp.Px()+fTm.Px(), fTp.Py()+fTm.Py(), fTp.Pz()+fTm.Pz(), fTp.E()+fTm.E());

    FillHisto("Reco_Mk", fTwoTrack.M());
    FillHisto("Reco_Pk", fTwoTrack.P());

    //std::cout <<  << std::endl;

    if(rstraw1_minus < 200 || rstraw1_plus < 200) return; //GEOMETRIC CUT
    FillHisto("Reco_piplus_xy_straw1" , posstraw1_plus.X(), posstraw1_plus.Y());
    FillHisto("Reco_piminus_xy_straw1", posstraw1_minus.X(), posstraw1_minus.Y());
    FillHisto("Reco_plus_vs_minus_Rstraw1", rstraw1_minus, rstraw1_plus);
    
    // fTwoTrack = fTp + fTm;

    fCDA = 999.0;
    fRecoVertex = fTool->SingleTrackVertex(fTp.Vect(), fTm.Vect(), posstraw1_plus, posstraw1_minus, fCDA);
    
    
    //std::cout<< "Recovtx = " << fRecoVertex.X() << "  " <<  fRecoVertex.Y() << " " << fRecoVertex.Z()  << "CDA = " <<  fCDA << std::endl;
    //std::cout<< "CDA = " <<  fCDA << std::endl;
    FillHisto("Reco_geom_cda_vs_zvtx", fRecoVertex.Z()/1e3, fCDA);
    FillHisto("Reco_geom_pk_vs_zvtx", fRecoVertex.Z()/1e3, fTwoTrack.P());
    FillHisto("Reco_geom_dzvtx", (fRecoVertex.Z()- fTrueVertex.Z())/1e3);
    FillHisto("Reco_geom_Pk", fTwoTrack.P());
    FillHisto("Reco_geom_Mk", fTwoTrack.M());
    FillHisto("Reco_geom_Pplus", fTp.P());
    FillHisto("Reco_geom_Pminus", fTm.P());

    //Reco zvtx DEBUGGING
    FillHisto("Reco_geom_posstraw1_plus_vs_minus", posstraw1_plus.X(), posstraw1_minus.X());
    FillHisto("Reco_geom_fTp_vs_fTm", fTp.Vect().X(), fTm.Vect().X());
    FillHisto("Reco_geom_zvtx", fRecoVertex.Z()/1e3);

    FillHisto("Reco_geom_dPk_vs_Pktrue",fBeamParticle->GetInitial4Momentum().P()/1e3, fTwoTrack.P() - fBeamParticle->GetInitial4Momentum().P()/1e3);
    FillHisto("Reco_geom_dThetak_vs_Pktrue",fBeamParticle->GetInitial4Momentum().P()/1e3, fTwoTrack.Theta() - fBeamParticle->GetInitial4Momentum().Theta());



    TIter partit_reco_geom(evt->GetKineParts());
    int npidec_bstraw4 = 0; // straw 1,2,3,4 = [212.9, 219.5, 229.8, 236.8] m
    while(partit_reco_geom.Next()) {  // itterate through ALL particles
        KinePart *particle = static_cast<KinePart *>(*partit_reco_geom);
        if(particle->GetPDGcode() != 130 && particle->GetPDGcode() != 310)
            FillHisto("Reco_geom_PDG", particle->GetPDGcode());


        if(particle->GetPDGcode() == 211 ){
            if(particle->GetEndPos().Z() < 237000){
                npidec_bstraw4 += 1;

            } else {

                FillHisto("Reco_geom_dPplus_vs_Pplus_true", particle->GetInitial4Momentum().P()/1e3, fTp.P() - particle->GetInitial4Momentum().P()/1e3 );
                FillHisto("Reco_geom_dTheta_plus_vs_Theta_plus_true", particle->GetInitial4Momentum().Theta(), fTp.Theta() - particle->GetInitial4Momentum().Theta() );
                FillHisto("Reco_geom_dTheta_plus_vs_Pplus_true",particle->GetInitial4Momentum().P()/1e3, fTp.Theta() - particle->GetInitial4Momentum().Theta() );
            }

        }

        if(particle->GetPDGcode() == -211 ){
            if(particle->GetEndPos().Z() < 237000){
                npidec_bstraw4 += 1;

            } else {
                FillHisto("Reco_geom_dPminus_vs_Pminus_true", particle->GetInitial4Momentum().P()/1e3, fTm.P() - particle->GetInitial4Momentum().P()/1e3 );
                FillHisto("Reco_geom_dTheta_minus_vs_Theta_minus_true", particle->GetInitial4Momentum().Theta(), fTm.Theta() - particle->GetInitial4Momentum().Theta() );
                FillHisto("Reco_geom_dTheta_minus_vs_Pminus_true", particle->GetInitial4Momentum().P()/1e3, fTm.Theta() - particle->GetInitial4Momentum().Theta() );
            }

        }



    }

    FillHisto("Reco_geom_Npidecays_before_straw4", npidec_bstraw4);
    FillHisto("Reco_geom_NKineParts", evt->GetNKineParts());

    double beta_K    = fTwoTrack.P()/TMath::Sqrt(fTwoTrack.P()*fTwoTrack.P() + MK0/1000.*MK0/1000.);
    //double beta_K    = pk_sim/TMath::Sqrt(pk_sim*pk_sim + MLAMBDA_GEV*MLAMBDA_GEV);
    double gamma_K   = 1/TMath::Sqrt(1 - beta_K*beta_K);
    double tau_ks    = (fRecoVertex.Z() - fTargetZ)/(fctau_ks*gamma_K);


    FillHisto("Reco_geom_gammaK_vs_Pk", fTwoTrack.P(), gamma_K);
    FillHisto("Reco_geom_gammaK_true_vs_Pk_true", fBeamParticle->GetInitial4Momentum().P()/1e3, gamma_K_true);
    FillHisto("Reco_geom_dgammaK_vs_Pk_true", fBeamParticle->GetInitial4Momentum().P()/1e3, gamma_K - gamma_K_true);



    if(npidec_bstraw4 > 0) return;
    FillHisto("Reco_final_Pk", fTwoTrack.P());
    FillHisto("Reco_final_Mk", fTwoTrack.M());
    FillHisto("Reco_final_Mk_vs_Pk",fTwoTrack.P() ,fTwoTrack.M());
    FillHisto("Reco_final_Pplus", fTp.P());
    FillHisto("Reco_final_Pminus", fTm.P());

    FillHisto("Reco_final_dPk_vs_Pktrue",fBeamParticle->GetInitial4Momentum().P()/1e3, fTwoTrack.P() - fBeamParticle->GetInitial4Momentum().P()/1e3);
    FillHisto("Reco_final_dThetak_vs_Pktrue",fBeamParticle->GetInitial4Momentum().P()/1e3, fTwoTrack.Theta() - fBeamParticle->GetInitial4Momentum().Theta());

    FillHisto("Reco_final_tauks_vs_Pk", fTwoTrack.P(), tau_ks);
    FillHisto("Reco_final_tauks_vs_zvtx", fRecoVertex.Z()/1e3, tau_ks);
    FillHisto("Reco_final_dtauks_vs_zvtx_true", fTrueVertex.Z()/1e3, tau_ks - tau_ks_true);
    FillHisto("Reco_final_dtauks_vs_Pk_true", fBeamParticle->GetInitial4Momentum().P()/1e3, tau_ks - tau_ks_true);
    FillHisto("Reco_final_cda_vs_zvtx", fRecoVertex.Z()/1e3, fCDA);
    FillHisto("Reco_final_dgammaK_vs_Pk_true", fBeamParticle->GetInitial4Momentum().P()/1e3, gamma_K - gamma_K_true);

    FillHisto("Reco_final_dMk_vs_Pk_true",fBeamParticle->GetInitial4Momentum().P()/1e3 ,fTwoTrack.M() - MK0_GEV);
    FillHisto("Reco_final_dMk_vs_Thetak_true",fBeamParticle->GetInitial4Momentum().Theta() ,fTwoTrack.M() - MK0_GEV);


    //tau_ks DEBUGGING
    FillHisto("Reco_final_tauks_vs_tauks_true", tau_ks_true, tau_ks);
    FillHisto("Reco_final_zvtx_vs_zvtx_true", fTrueVertex.Z()/1e3, fRecoVertex.Z()/1e3);


    return;
   }
}

TLorentzVector KPIPI_analysis::GetTrack(TRecoSpectrometerCandidate* tr, Double_t mass){

    TLorentzVector track;
    //double px = tr->GetThreeMomentumBeforeMagnet().X() * 0.001;
    //double py = tr->GetThreeMomentumBeforeMagnet().Y() * 0.001;
    //double pz = tr->GetThreeMomentumBeforeMagnet().Z() * 0.001;
    //
    //double E  = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);
    //std::cout << " Momentum =" << px << "  " << py << "  " << pz << " Compute = " << TMath::Sqrt(px*px+py*py+pz*pz) << "From Reco = " << tr->GetMomentum()/1000.<< "  "<< mass << " E = " << E << std::endl;


    // track.SetPxPyPzE(px,py,pz,E);
    track.SetVectM(tr->GetThreeMomentumBeforeMagnet()*0.001,mass);

    return track;
}


void KPIPI_analysis::PostProcess() {
}

void KPIPI_analysis::EndOfBurstUser() {
}

void KPIPI_analysis::EndOfRunUser() {
}

void KPIPI_analysis::EndOfJobUser() {
    SaveAllPlots();
}

void KPIPI_analysis::DrawPlot() {
}

KPIPI_analysis::~KPIPI_analysis() {
}
