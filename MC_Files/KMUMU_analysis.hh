#ifndef KMUMU_ANALYSIS_HH
#define KMUMU_ANALYSIS_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "PnnAnalysisTools.hh"

#include "TLorentzVector.h"
#include "TVector3.h"


class TH1I;
class TH2F;
class TGraph;
class TTree;

class TRecoSpectrometerEvent;
class TRecoSpectrometerCandidate;

class KMUMU_analysis : public NA62Analysis::Analyzer
{
public:
  explicit KMUMU_analysis(NA62Analysis::Core::BaseAnalysis *ba);
  ~KMUMU_analysis();
  void InitHist();
  void InitOutput();
  void ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType);
  void Process(int iEvent);
  void PostProcess();
  void StartOfBurstUser();
  void EndOfBurstUser();
  void StartOfRunUser();
  void EndOfRunUser();
  void EndOfJobUser();
  void DrawPlot();

  void Clear();
    TLorentzVector GetTrack(TRecoSpectrometerCandidate*, Double_t);
protected:
    TRecoSpectrometerEvent* fSpectrometerEvent;
    TLorentzVector fTp;
    TLorentzVector fTm;
    TLorentzVector fTwoTrack;
    TLorentzVector fTrueVertex;
    TVector3 fRecoVertex;
    double fCDA;
    bool fIsKSBeam;
    bool fIsKLBeam;
    double fTargetZ;
    double fctau_kl;
    double fctau_ks;

    KinePart* fBeamParticle;
    TRecoSpectrometerCandidate* fSplus;
    TRecoSpectrometerCandidate* fSminus;
    PnnAnalysisTools* fTool;
};
#endif
