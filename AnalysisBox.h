#ifndef ACANALYSISBOX
#define ACANALYSISBOX
#include "Analyse.h"
#include "DataForm.h"

using namespace std;
using namespace TMath;

class AnalysisBox
{
 private:
  Analyse* an;
  
 public:
  AnalysisBox(Analyse* analyse);
  Long64_t GetNumAddedEvents() const {return(an->GetNumAddedEvents());}
  string GetCurrentFileName() const {return(an->GetCurrentFileName());}
  
  Double_t Number() const {return(an->Number());}
  UInt_t Run() const {return(an->Run());}
  UInt_t LumiBlock() const {return(an->LumiBlock());}
  UInt_t TimeUnix() const {return(an->TimeUnix());}
  UInt_t TimeMicroSec() const {return(an->TimeMicroSec());}
  Double_t AK5PFRho() const {return(an->AK5PFRho());}
  Double_t AK5PFSigma() const {return(an->AK5PFSigma());}
  
  //RECO-level information
  BeamSpot GetBeamSpot() const {return(an->GetBeamSpot());}
  
  TrackComposedParticle SecVertices(UInt_t n) const {return(an->SecVertices(n));}
  UInt_t NumSecVertices() const {return(an->NumSecVertices());}
  
  Muon Muons(UInt_t n) const {return(an->Muons(n));}
  UInt_t NumMuons() const {return(an->NumMuons());}
  
  Electron Electrons(UInt_t n) const {return(an->Electrons(n));}
  UInt_t NumElectrons() const {return(an->NumElectrons());}
  
  Tau Taus(UInt_t n) const {return(an->Taus(n));}
  UInt_t NumTaus() const {return(an->NumTaus());}
  MuTauTauPair MuTauTauPairs(UInt_t n) const {return(an->MuTauTauPairs(n));}
  UInt_t NumMuTauTauPairs() const {return(an->NumMuTauTauPairs());}
  ElTauTauPair ElTauTauPairs(UInt_t n) const {return(an->ElTauTauPairs(n));}
  UInt_t NumElTauTauPairs() const {return(an->NumElTauTauPairs());}
  
  Photon Photons(UInt_t n) const {return(an->Photons(n));}
  UInt_t NumPhotons() const {return(an->NumPhotons());}
  
  TLorentzVector CaloMET() const {return(an->CaloMET());}
  TLorentzVector CaloMETMuons() const {return(an->CaloMETMuons());}
  TLorentzVector PFMET() const {return(an->PFMET());}
  TLorentzVector PFMETTYPE1() const {return(an->PFMETTYPE1());}
  TLorentzVector TCMET() const {return(an->TCMET());}
  
  Jet AK5CaloJets(UInt_t n) const {return(an->AK5CaloJets(n));}
  UInt_t NumAK5CaloJets() const {return(an->NumAK5CaloJets());}
  
  Jet AK5JPTJets(UInt_t n) const {return(an->AK5JPTJets(n));}
  UInt_t NumAK5JPTJets() const {return(an->NumAK5JPTJets());}
  
  Jet AK5PFJets(UInt_t n) const {return(an->AK5PFJets(n));}
  UInt_t NumAK5PFJets() const {return(an->NumAK5PFJets());}
  
  Track Tracks(UInt_t n) const {return(an->Tracks(n));}
  UInt_t NumTracks() const {return(an->NumTracks());}
  
  Vertex PrimVertex() const {return(an->PrimVertex());}
  UInt_t NumPrimVertices() const {return(an->NumPrimVertices());}
  
  SuperCluster SuperClusters(UInt_t n, TVector3 refpoint) const {return(an->SuperClusters(n, refpoint));}
  SuperCluster SuperClusters(UInt_t n) const {return(an->SuperClusters(n));}
  UInt_t NumSuperClusters() const {return(an->NumSuperClusters());}
  
  Int_t NumPileUpInteractionsMinus() const {return(an->NumPileUpInteractionsMinus());}
  Int_t NumPileUpInteractions() const {return(an->NumPileUpInteractions());}
  Int_t NumPileUpInteractionsPlus() const {return(an->NumPileUpInteractionsPlus());}
  Float_t NumTruePileUpInteractions() const {return(an->NumTruePileUpInteractions());}
  Double_t GetPileUpWeight(Double_t meaninteractions) const {return(an->GetPileUpWeight(meaninteractions));}
  Double_t GetPrimVertexWeight(Double_t meaninteractions) const {return(an->GetPrimVertexWeight(meaninteractions));}
  Double_t GetPrimVertexWeight(vector<Double_t>& datadist) const {return(an->GetPrimVertexWeight(datadist));}
  Int_t NumGoodPrimVertices() const {return(an->NumGoodPrimVertices());}
  
  //generator-level information
  Double_t GenWeight() const {return(an->GenWeight());}
  Double_t GenId1() const {return(an->GenId1());}
  Double_t Genx1() const {return(an->Genx1());}
  Double_t GenId2() const {return(an->GenId2());}
  Double_t Genx2() const {return(an->Genx2());}
  Double_t GenScale() const {return(an->GenScale());}
  
  GenParticle AllGenParticles(UInt_t n) const {return(an->AllGenParticles(n));}
  UInt_t NumAllGenParticles() const {return(an->NumAllGenParticles());}
  
  GenLightParticle GenParticles(UInt_t n) const {return(an->GenParticles(n));}
  UInt_t NumGenParticles() const {return(an->NumGenParticles());}
  
  TLorentzVector GenMETCalo() const {return(an->GenMETCalo());}
  TLorentzVector GenMETTrue() const {return(an->GenMETTrue());}
  
  Long64_t Processed() const {return(an->Processed());}
  
  //trigger information
  bool GetL1Trigger(UInt_t bit) const {return(an->GetL1Trigger(bit));}
  bool GetL1TriggerBits(UInt_t bit) const {return(an->GetL1TriggerBits(bit));}
  bool GetHLTrigger(UInt_t index) const {return(an->GetHLTrigger(index));}
  Int_t GetHLTrigger(vector<string> triggernames) const {return(an->GetHLTrigger(triggernames));}
  Int_t GetHLTriggerIndex(string triggername) {return(an->GetHLTriggerIndex(triggername));}
  string GetHLTriggerName(UInt_t index) {return(an->GetHLTriggerName(index));}
  Int_t GetHLTPrescale(UInt_t triggerindex) {return(an->GetHLTPrescale(triggerindex));}
  Int_t GetNumHLTriggers() {return(an->GetNumHLTriggers());}
  
  TriggerSelection* GetTriggerSelection(string id) {return(an->GetTriggerSelection(id));}
  
  Int_t IsLumiAvailable() {return(an->IsLumiAvailable());}
  Double_t GetInstLumi() {return(an->GetInstLumi());}
  Double_t GetLumi(Int_t format = 0) {return(an->GetLumi(format));}
  Double_t GetLumiBlockLumi() {return(an->GetLumiBlockLumi());}
};

#endif
