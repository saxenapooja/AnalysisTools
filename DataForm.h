#ifndef ACDATA
#define ACDATA

#include <TLorentzVector.h>
#include <TMatrixFSym.h>
#include <iostream>
#include <string>
#include <map>
#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class Analyse;
class Vertex;

using namespace std;
using namespace TMath;

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LVec;


const UInt_t M_trackmaxcount = 1000;
const UInt_t M_superclustermaxcount = 1000;
const UInt_t M_superclustermembermaxcount = 1000;
const UInt_t M_superclusterhitmaxcount = 5000;
const UInt_t M_primvertexmaxcount = 1000;
const UInt_t M_muonmaxcount = 1000;
const UInt_t M_taumaxcount = 1000;
const UInt_t M_mutautaupairmaxcount = 1000;
const UInt_t M_eltautaupairmaxcount = 1000;
const UInt_t M_electronmaxcount = 1000;
const UInt_t M_photonmaxcount = 1000;
const UInt_t M_conversionmaxcount = 1000;
const UInt_t M_jetmaxcount = 1000;
const UInt_t M_secverticesmaxcount = 1000;
const UInt_t M_genallparticlesmaxcount = 10000;
const UInt_t M_genparticlesmaxcount = 1000;
const UInt_t M_genmotherdaughtermaxcount = 50000;
const UInt_t M_btagmax = 6;
const UInt_t M_svmax = 5;

const Double_t MuonMass      = 0.105658;
const Double_t MuonMassQ     = MuonMass*MuonMass;
const Double_t ElectronMass  = 0.000511;
const Double_t ElectronMassQ = ElectronMass*ElectronMass;
const Double_t TauMass       = 1.77699;
const Double_t TauMassQ      = TauMass*TauMass;


class GenLightParticle : public TLorentzVector
{
 private:
  TVector3 vertex;
  UInt_t   info;
  Int_t    status;
  Int_t    pdgid;
 
 public:
  GenLightParticle(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Info);
  GenLightParticle() : pdgid(0) {}
    Int_t Status() const {return(status);}
    Int_t PDGId() const {return(pdgid);}
    TVector3 Vertex() const {return(vertex);}
    enum INFO {FromZG, FromW, FromTop};
    bool MoreInfo(INFO test) const {return(info & 1<<test);}
};

class GenParticle : public TLorentzVector
{
 private:
  const Analyse* myanalyse;
  TVector3 vertex;
  Int_t status;
  Int_t charge_;
  Int_t pdgid;
  LorentzVector p4_;
  UInt_t myindex;
  UInt_t motherfirst;
  UInt_t mothernum;
  UInt_t daughterfirst;
  UInt_t daughternum;
  bool   HasAnyMotherPDGId(vector<UInt_t>& visited, Int_t pdgid, bool antiparticle = true);
  UInt_t GetMotherIndex(UInt_t num) const;
  UInt_t GetDaughterIndex(UInt_t num) const;
  
 public:
  GenParticle(const Analyse* Myanalyse, UInt_t Myindex, Double_t E, Double_t Px, Double_t Py, Double_t Pz, 
	      Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Motherfirst, UInt_t Mothernum, 
	      UInt_t Daughterfirst, UInt_t Daughternum, Int_t charge, LorentzVector p4);
  GenParticle() : myanalyse(0), pdgid(0) {}
    Int_t Status() const {return(status);}
    Int_t PDGId() const {return(pdgid);}
    Int_t Charge() const {return(charge_);}
    LorentzVector P4() const { return(p4_);}
    TVector3 Vertex() const {return(vertex);}
    UInt_t NumMothers() const {return(mothernum);}
    GenParticle GetMother(UInt_t num) const;
    UInt_t NumDaughters() const {return(daughternum);}
    UInt_t GetMyIndex() const {return(myindex);}
    GenParticle GetDaughter(UInt_t num) const;
    bool HasAnyMotherPDGId(Int_t pdgid, bool antiparticle = true);
    bool operator ==(GenParticle& other) const {return(GetMyIndex() == other.GetMyIndex());}
    bool operator !=(GenParticle& other) const {return(GetMyIndex() != other.GetMyIndex());}
};

class EcalHit : public TVector3
{
 private:
  Double_t energy;
 public:
  EcalHit() {}
  EcalHit(Double_t Energy, Double_t X, Double_t Y, Double_t Z);
  Double_t E() const {return(energy);}
};


class Cluster : public TLorentzVector
{
  friend class Analyse;
 private:
  TVector3 position;
  Int_t size;
  vector<EcalHit> hits;
  void AddHit(EcalHit& hit);
 public:
  Cluster() {}
  Cluster(Double_t E, Double_t x, Double_t y, Double_t z, Int_t Size);
  Int_t Size() const {return(size);}
  TVector3 Position() const {return(position);}
  Int_t NumHits() const {return(hits.size());}
  EcalHit Hits(UInt_t n) const {return(hits[n]);};
};

class SuperCluster : public TLorentzVector
{
  friend class Analyse;
 private:
  const Analyse* myanalyse;
  TVector3 position;
  Float_t rawe;
  Float_t phiwidth;
  Float_t etawidth;
  vector<Cluster> clusters;
  vector<Cluster> esclusters;
  void AddCluster(Cluster& newcluster);
  void AddESCluster(Cluster& newcluster);
 public:
  SuperCluster() {}
  SuperCluster(Double_t E, Double_t x, Double_t y, Double_t z, Float_t Rawe, Float_t Phiwidth, Float_t Etawidth);
  Float_t RawEnergy() const {return(rawe);}
  Float_t PhiWidth() const {return(phiwidth);}
  Float_t EtaWidth() const {return(etawidth);}
  TVector3 Position() const {return(position);}
  Cluster Clusters(UInt_t n) const {return(clusters[n]);}
  UInt_t NumClusters() const {return(clusters.size());}
  Cluster ESClusters(UInt_t n) const {return(esclusters[n]);}
  UInt_t NumESClusters() const {return(esclusters.size());}
};

class Track : public TLorentzVector
{
 private:
  TVector3 outerpoint;
  TVector3 closestpoint;
  Double_t chi2;
  Double_t ndof;
  Double_t dxy;
  Double_t dxyerr;
  Double_t dz;
  Double_t dzerr;
  Int_t charge;
  Int_t nhits;
  Int_t nmissinghits;
  Int_t npixelhits;
  Int_t npixellayers;
  Int_t nstriplayers;
  Float_t dedxharmonic2;
 public:
  Track() {}
  Track(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Double_t Ox, Double_t Oy, Double_t Oz, Double_t Cx, Double_t Cy, Double_t Cz, Double_t Chi2, Double_t Ndof, Double_t Dxy, Double_t Dxyerr, Double_t Dz, Double_t Dzerr, Int_t Charge, Int_t Nhits, Int_t Nmissinghits, Int_t Npixelhits, Int_t Npixellayers, Int_t Nstriplayers, Float_t Dedxharmonic2);
  Double_t Chi2() const {return(chi2);}
  Double_t Ndof() const {return(ndof);}
  Double_t Chi2OverNdof() const {return(chi2/ndof);}
  Int_t Charge() const {return(charge);}
  Int_t NHits() const {return(nhits);}
  Int_t NMissingHits() const {return(nmissinghits);}
  Int_t NPixelHits() const {return(npixelhits);}
  Int_t NPixelLayers() const {return(npixellayers);}
  Int_t NStripLayers() const {return(nstriplayers);}
  TVector3 ECalPoint() const {return(outerpoint);}
  TVector3 ClosestPoint() const {return(closestpoint);}
  TVector3 MomentumSpace() const {return(TVector3(Px(), Py(), Pz()));}
  Float_t dEdxHarmonic2() const {return(dedxharmonic2);}
  Double_t Dxy() const {return(dxy);}
  Double_t DxyError() const {return(dxyerr);}
  Double_t DxySig() const {return(dxy/dxyerr);}
  Double_t Dz() const {return(dz);}
  Double_t DzError() const {return(dzerr);}
  Double_t DzSig() const {return(dz/dzerr);}
};


class TrackComposedParticle : public TVector3 
{
  friend class Analyse;
 private:
  TMatrixFSym covmatrix;
  Double_t chi2;
  Double_t ndof;
  vector<Track> tracks;
  void AddTrack(Track& newtrack);
 public:
  TrackComposedParticle() {}
  TrackComposedParticle(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t* Cov);
  Double_t Chi2() const {return(chi2);}
  Double_t Ndof() const {return(ndof);}
  Double_t Chi2OverNdof() const {return(chi2/ndof);}
  Track Tracks(Int_t n) const {return(tracks[n]);}
  Int_t NumTracks() const {return(tracks.size());}
  Double_t XError() const {return(Sqrt(covmatrix(0,0)));}
  Double_t YError() const {return(Sqrt(covmatrix(1,1)));}
  Double_t ZError() const {return(Sqrt(covmatrix(2,2)));}
  TMatrixFSym CovMatrix() const {return(covmatrix);}
  Double_t CovMatrix(Int_t i, Int_t j) const {return(covmatrix(i,j));}
  Double_t VertexSig3D(Vertex& vertex) const;
  Double_t VertexSig2D(Vertex& vertex) const;
};

class TriggerObject
{
 private:
  const Analyse* MA;
  UInt_t trigger;
  const vector<string>* triggernames;
 public:
  TriggerObject() {}
  TriggerObject(const Analyse* ma, const vector<string>* Triggernames, UInt_t Trigger);
  //TriggerObject(const Analyse* ma, const vector<string>* Triggernames, UInt_t Trigger, bool isL1triggerd);
  Int_t Trigger(string name) const; //0: Trigger not available or matched. 1 not prescale fired, -1 not prescaled not fired, 2 prescaled fired, -2 prescalded not fired 
  const vector<string>* TriggerNames() const {return(triggernames);} 
  const Analyse* MyAn() const {return(MA);}
};

class Muon : public TLorentzVector, public TriggerObject
{
  friend class Analyse;
 private:
  Float_t pterror;
  Float_t chi2;
  Float_t ndof;
  Track innertrack;
  Track outertrack;
  Float_t isolationr3track;
  Int_t isolationr3ntrack;
  Float_t isolationr3ecal;
  Float_t isolationr3hcal;
  Float_t pfisolationr4_sumchargedhadronpt;
  Float_t pfisolationr4_sumchargedparticlept;
  Float_t pfisolationr4_sumneutralhadronet;
  Float_t pfisolationr4_sumphotonet;
  Float_t pfisolationr4_sumPUpt;
  Float_t ecalenergy;
  Float_t hcalenergy;
  Int_t charge;
  Int_t numchambers;
  Int_t numchamberswithsegments;
  Int_t numchamberhits;
  UInt_t type;
  UInt_t trackermuonquality;
  Float_t mva_id;
  Float_t mva_iso;
  void SetInnerTrack(Track& Innertrack) {innertrack = Innertrack;}
  void SetOuterTrack(Track& Outertrack) {outertrack = Outertrack;}
 public:
  Muon() {}
  Muon(const Analyse* ma, Double_t Px, Double_t Py, Double_t Pz, Float_t Pterror, Float_t Chi2, Float_t Ndof, Float_t Pfisolationr4_sumchargedhadronpt, Float_t Pfisolationr4_sumchargedparticlept, Float_t Pfisolationr4_sumneutralhadronet,Float_t Pfisolationr4_sumphotonet, Float_t Pfisolationr4_sumPUpt, Float_t Isolationr3track, Int_t Isolationr3ntrack, Float_t Isolationr3ecal, Float_t Isolationr3hcal, Float_t Ecalenergy, Float_t Hcalenergy, Int_t Charge, Int_t Numchambers, Int_t Numchamberswithsegments, Int_t Numchamberhits, UInt_t Type, UInt_t Trigger, UInt_t Trackermuonquality, Float_t MVAId, Float_t MVAIso);
  Float_t PtError() const {return(pterror);}
  Float_t IsoR3Track() const {return(isolationr3track);}
  Float_t IsoR3TrackRel() const {return(isolationr3track/Pt());}
  Int_t IsoR3NTrack() const {return(isolationr3ntrack);}
  Float_t IsoR3ECal() const {return(isolationr3ecal);}
  Float_t IsoR3HCal() const {return(isolationr3hcal);}
  Float_t IsoR3Combined() const {return(IsoR3Track()+IsoR3ECal()+IsoR3HCal());}
  Float_t IsoR3CombinedRel() const {return(IsoR3Combined()/Pt());}
  Float_t PFIsoR4ChargedHad() const {return(pfisolationr4_sumchargedhadronpt);}
  Float_t PFIsoR4NeutralHad() const {return(pfisolationr4_sumneutralhadronet);}
  Float_t PFIsoR4Photon() const {return(pfisolationr4_sumphotonet);}
  Float_t PFIsoR4ChargedParticle() const {return(pfisolationr4_sumchargedparticlept);}
  Float_t PFIsoR4SumPUPt() const {return(pfisolationr4_sumPUpt);}
  Float_t PFIsoR4() const {return(PFIsoR4ChargedParticle()+PFIsoR4NeutralHad()+PFIsoR4Photon());}
  Float_t PFIsoR4DB() const {return(PFIsoR4ChargedParticle()+max(0.0f, PFIsoR4NeutralHad()+PFIsoR4Photon()-0.5f*PFIsoR4SumPUPt()));}
  Float_t PFIsoR4DBCH() const {return(PFIsoR4ChargedHad()+max(0.0f, PFIsoR4NeutralHad()+PFIsoR4Photon()-0.5f*PFIsoR4SumPUPt()));}
  Float_t PFIsoR4Rel() const {return(PFIsoR4()/Pt());}
  Float_t PFIsoR4RelDB() const {return(PFIsoR4DB()/Pt());}
  Float_t PFIsoR4RelDBCH() const {return(PFIsoR4DBCH()/Pt());}
  Float_t MVAId() const { return mva_id; }
  Float_t MVAIso() const { return mva_iso; }
  Float_t ECalEnergy() const {return(ecalenergy);}
  Float_t HCalEnergy() const {return(hcalenergy);}
  Int_t Charge() const {return(charge);}
  Int_t NumChambers() const {return(numchambers);}
  Int_t NumChambersWithSegments() const {return(numchamberswithsegments);}
  Int_t NumChamberHits() const {return(numchamberhits);}
  bool IsGlobal() const {return((type & 1) > 0);}
  bool IsTracker() const {return((type & 1<<1) > 0);}
  bool IsStandAlone() const {return((type & 1<<2) > 0);}
  bool IsCalo() const {return((type & 1<<3) > 0);}
  bool IsPF() const {return((type & 1<<8) > 0);}
  bool HasInnerTrack() const {return((type & 1<<4) > 0);}
  bool HasOuterTrack() const {return((type & 1<<5) > 0);}
  bool InOut() const {return((type & 1<<6) > 0);}
  bool OutIn() const {return((type & 1<<7) > 0);}
  Double_t Dxy() const {return(innertrack.Dxy());}
  Double_t DxyError() const {return(innertrack.DxyError());}
  Double_t Dz() const {return(innertrack.Dz());}
  Double_t DzError() const {return(innertrack.DzError());}
  Float_t Chi2() const {return(chi2);}
  Float_t Ndof() const {return(ndof);}
  Float_t Chi2OverNdof() const {return(chi2/ndof);}
  Track InnerTrack() const {return(innertrack);}
  Track OuterTrack() const {return(outertrack);}
  Int_t NumStations() const;
  //Quality
  bool IsAll() const {return(trackermuonquality & 1 << 0);}
  bool IsAllGlobalMuons() const {return(trackermuonquality & 1 << 1);}
  bool IsAllStandAloneMuons() const {return(trackermuonquality & 1 << 2);}
  bool IsAllTrackerMuons() const {return(trackermuonquality & 1 << 3);}
  bool IsTrackerMuonArbitrated() const {return(trackermuonquality & 1 << 4);}
  bool IsAllArbitrated() const {return(trackermuonquality & 1 << 5);}
  bool IsGlobalMuonPromptTight() const {return(trackermuonquality & 1 << 6);}
  bool IsTMLastStationLoose() const {return(trackermuonquality & 1 << 7);}
  bool IsTMLastStationTight() const {return(trackermuonquality & 1 << 8);}
  bool IsTM2DCompatibilityLoose() const {return(trackermuonquality & 1 << 9);}
  bool IsTM2DCompatibilityTight() const {return(trackermuonquality & 1 << 10);}
  bool IsTMOneStationLoose() const {return(trackermuonquality & 1 << 11);}
  bool IsTMOneStationTight() const {return(trackermuonquality & 1 << 12);}
  bool IsTMLastStationOptimizedLowPtLoose() const {return(trackermuonquality & 1 << 13);}
  bool IsTMLastStationOptimizedLowPtTight() const {return(trackermuonquality & 1 << 14);}
  //only in > CMSSW_3_1_X
  bool IsGMTkChiCompatibility() const {return(trackermuonquality & 1 << 15);}
  bool IsGMStaChiCompatibility() const {return(trackermuonquality & 1 << 16);}
  bool IsGMTkKinkTight() const {return(trackermuonquality & 1 << 17);}
  bool IsTMLastStationAngLoose() const {return(trackermuonquality & 1 << 18);}
  bool IsTMLastStationAngTight() const {return(trackermuonquality & 1 << 19);}
  bool IsTMOneStationAngLoose() const {return(trackermuonquality & 1 << 20);}
  bool IsTMOneStationAngTight() const {return(trackermuonquality & 1 << 21);}
  bool IsTMLastStationOptimizedBarrelLowPtLoose() const {return(trackermuonquality & 1 << 22);}
};

class Electron : public TLorentzVector, public TriggerObject
{
  friend class Analyse;
 private:
  TVector3 outerpoint;
  TVector3 closestpoint;
  Double_t dxy;
  Double_t dxyerr;
  Double_t dz;
  Double_t dzerr;
  Float_t esuperclusterovertrack;
  Float_t eseedclusterovertrack;
  Float_t deltaetasuperclustertrack;
  Float_t deltaphisuperclustertrack;
  Float_t e1x5;
  Float_t e2x5;
  Float_t e5x5;
  Float_t sigmaetaeta;
  Float_t sigmaietaieta;
  Float_t ehcaloverecal;
  Float_t ehcaloverecaldepth1;
  Float_t ehcaloverecaldepth2;
  Float_t isolationr3track;
  Float_t isolationr3ecal;
  Float_t isolationr3hcal;
  Float_t isolationr4track;
  Float_t isolationr4ecal;
  Float_t isolationr4hcal;
  Float_t pfisolationr3_sumchargedhadronpt;
  Float_t pfisolationr3_sumchargedparticlept;
  Float_t pfisolationr3_sumneutralhadronet;
  Float_t pfisolationr3_sumphotonet;
  Float_t pfisolationr3_sumPUpt;
  Float_t pfisolationr4_sumchargedhadronpt;
  Float_t pfisolationr4_sumchargedparticlept;
  Float_t pfisolationr4_sumneutralhadronet;
  Float_t pfisolationr4_sumphotonet;
  Float_t pfisolationr4_sumPUpt;
  Float_t trackchi2;
  Float_t trackndof;
  Int_t nhits;
  Int_t nmissinghits;
  Int_t npixelhits;
  Int_t npixellayers;
  Int_t nstriplayers;
  Float_t convdist;
  Float_t convdcot;
  Float_t convradius;
  UInt_t gapinfo;
  UInt_t chargeinfo;
  Float_t fbrems;
  Int_t numbrems;
  Int_t charge;
  Byte_t info;
  vector<SuperCluster> supercluster;
  Float_t mvaidnontrig;
  Bool_t hasconversion;
  void AddSC(SuperCluster& sc);
  bool WorkingPoint(Int_t Missinghits, Float_t Convdist, Float_t Convdcot, Float_t Sigmaietaieta, Float_t Deltaphisctrack, Float_t Deltaetasctrack, Float_t Ehcaloverecal) const;
 public:
  Electron() {}
  Electron(const Analyse* ma, Double_t Px, Double_t Py, Double_t Pz, Float_t Trackchi2, Float_t Trackndof, Double_t Outerx, Double_t Outery, Double_t Outerz, Double_t Closestx, Double_t Closesty, Double_t Closestz, Double_t Dxy, Double_t Dxyerr, Double_t Dz, Double_t Dzerr, Float_t Esuperclusterovertrack, Float_t Eseedclusterovertrack, Float_t Deltaetasuperclustertrack, Float_t Deltaphisuperclustertrack, Float_t E1x5, Float_t E2x5, Float_t E5x5, Float_t Sigmaetaeta, Float_t Sigmaietaieta, Float_t Ehcaloverecal, Float_t Ehcaloverecaldepth1, Float_t Ehcaloverecaldepth2, Float_t Isolationr3track, Float_t Isolationr3ecal, Float_t Isolationr3hcal, Float_t Isolationr4track, Float_t Isolationr4ecal, Float_t Isolationr4hcal,Float_t Pfisolationr3_sumchargedhadronpt, Float_t Pfisolationr3_sumchargedparticlept, Float_t Pfisolationr3_sumneutralhadronet, Float_t Pfisolationr3_sumphotonet, Float_t Pfisolationr3_sumPUpt, Float_t  Pfisolationr4_sumchargedhadronpt, Float_t Pfisolationr4_sumchargedparticlept, Float_t Pfisolationr4_sumneutralhadronet, Float_t Pfisolationr4_sumphotonet, Float_t Pfisolationr4_sumPUpt, Int_t Nhits, Int_t Nmissinghits, Int_t Npixelhits, Int_t Npixellayers, Int_t Nstriplayers, Float_t Convdist, Float_t Convdcot, Float_t Convradius, UInt_t Gapinfo, UInt_t Chargeinfo, Float_t Fbrems, Int_t Numbrems, Int_t Charge, Byte_t Info, UInt_t Trigger, Float_t MVAIdNontrig, Bool_t HasConversion);
  TVector3 EcalPoint() const {return(outerpoint);}
  TVector3 ClosestPoint() const {return(closestpoint);}
  UInt_t NumSCs() const {return(supercluster.size());}
  SuperCluster SCs(UInt_t n) const {return(supercluster[n]);}
  Float_t ESuperClusterOverTrack() const {return(esuperclusterovertrack);}
  Float_t ESeedClusterOverTrack() const {return(eseedclusterovertrack);}
  Float_t DeltaEtaSuperClusterTrack() const {return(deltaetasuperclustertrack);}
  Float_t DeltaPhiSuperClusterTrack() const {return(deltaphisuperclustertrack);}
  Float_t E1x5() const {return(e1x5);}
  Float_t E2x5() const {return(e2x5);}
  Float_t E5x5() const {return(e5x5);}
  Float_t SigmaEtaEta() const {return(sigmaetaeta);}
  Float_t SigmaIEtaIEta() const {return(sigmaietaieta);}
  Float_t EHcalOverECal() const {return(ehcaloverecal);}
  Float_t EHcalOverECalDepth1() const {return(ehcaloverecaldepth1);}
  Float_t EHcalOverECalDepth2() const {return(ehcaloverecaldepth2);}
  Float_t IsoR3Track() const {return(isolationr3track);}
  Float_t IsoR3TrackRel() const {return(isolationr3track/Pt());}
  Float_t IsoR3ECal() const {return(isolationr3ecal);}
  Float_t IsoR3ECalRel() const {return(IsoR3ECal()/Pt());}
  Float_t IsoR3HCal() const {return(isolationr3hcal);}
  Float_t IsoR3HCalRel() const {return(IsoR3HCal()/Pt());}
  Float_t IsoR3Combined() const {return(IsoR3Track() + IsoR3ECal() + IsoR3HCal());}
  Float_t IsoR3CombinedRel() const {return(IsoR3Combined()/Pt());}
  Float_t IsoR4Track() const {return(isolationr4track);}
  Float_t IsoR4ECal() const {return(isolationr4ecal);}
  Float_t IsoR4HCal() const {return(isolationr4hcal);}
  Float_t Pfisolationr3_sumchargedhadronpt() const {return(pfisolationr3_sumchargedhadronpt);}
  Float_t Pfisolationr3_sumchargedparticlept() const {return(pfisolationr3_sumchargedparticlept);}
  Float_t Pfisolationr3_sumneutralhadronet() const {return(pfisolationr3_sumneutralhadronet);}
  Float_t Pfisolationr3_sumphotonet() const {return(pfisolationr3_sumphotonet);}
  Float_t Pfisolationr3_sumPUpt() const {return(pfisolationr3_sumPUpt);}
  Float_t PFIsoR3() const {return(Pfisolationr3_sumchargedhadronpt()+Pfisolationr3_sumneutralhadronet()+Pfisolationr3_sumphotonet());}
  Float_t PFIsoR3DB() const {return(Pfisolationr3_sumchargedhadronpt()+max(0.0f, Pfisolationr3_sumneutralhadronet()+Pfisolationr3_sumphotonet()-0.5f*Pfisolationr3_sumPUpt()));}
  Float_t PFIsoR3Rel() const {return(PFIsoR3()/Pt());}
  Float_t PFIsoR3RelDB() const {return(PFIsoR3DB()/Pt());}
  Float_t Pfisolationr4_sumchargedhadronpt() const {return(pfisolationr4_sumchargedhadronpt);}
  Float_t Pfisolationr4_sumchargedparticlept() const {return(pfisolationr4_sumchargedparticlept);}
  Float_t Pfisolationr4_sumneutralhadronet() const {return(pfisolationr4_sumneutralhadronet);}
  Float_t Pfisolationr4_sumphotonet() const {return(pfisolationr4_sumphotonet);}
  Float_t Pfisolationr4_sumPUpt() const {return(pfisolationr4_sumPUpt);}
  Float_t PFIsoR4() const {return(Pfisolationr4_sumchargedparticlept()+Pfisolationr4_sumneutralhadronet()+Pfisolationr4_sumphotonet());}
  Float_t PFIsoR4DB() const {return(Pfisolationr4_sumchargedparticlept()+max(0.0f, Pfisolationr4_sumneutralhadronet()+Pfisolationr4_sumphotonet()-0.5f*Pfisolationr4_sumPUpt()));}
  Float_t PFIsoR4DBCH() const {return(Pfisolationr4_sumchargedhadronpt()+max(0.0f, Pfisolationr4_sumneutralhadronet()+Pfisolationr4_sumphotonet()-0.5f*Pfisolationr4_sumPUpt()));}
  Float_t PFIsoR4Rel() const {return(PFIsoR4()/Pt());}
  Float_t PFIsoR4RelDB() const {return(PFIsoR4DB()/Pt());}
  Float_t PFIsoR4RelDBCH() const {return(PFIsoR4DBCH()/Pt());}
  Float_t TrackChi2() const {return(trackchi2);}
  Float_t TrackNdof() const {return(trackndof);}
  Float_t TrackChi2OverNdof() const {return(trackchi2/trackndof);}
  Int_t NHits() const {return(nhits);}
  Int_t NMissingHits() const {return(nmissinghits);}
  Int_t NPixelHits() const {return(npixelhits);}
  Int_t NPixelLayers() const {return(npixellayers);}
  Int_t NStripLayers() const {return(nstriplayers);}
  Float_t ConversionDist() const {return(convdist);}
  Float_t ConversionDCot() const {return(convdcot);}
  Float_t ConversionRadius() const {return(convradius);}
  Float_t FractionBrems() const {return(fbrems);}
  Int_t NumBrems() const {return(numbrems);}
  Int_t Charge() const {return(charge);}
  Byte_t Info() const {return(info);}
  Double_t Dxy() const {return(dxy);}
  Double_t DxyError() const {return(dxyerr);}
  Double_t Dz() const {return(dz);}
  Double_t DzError() const {return(dzerr);}
  Float_t MVAIdNonTrig() const { return(mvaidnontrig); }
  Bool_t HasConversion() const { return(hasconversion); }
  //Working Points: Choose combined 1 or relative isolation 0!
  bool WP95_v1(Int_t combined = 0) const; 
  bool WP90_v1(Int_t combined = 0) const; 
  bool WP85_v1(Int_t combined = 0) const; 
  bool WP80_v1(Int_t combined = 0) const; 
  bool WP70_v1(Int_t combined = 0) const; 
  bool WP60_v1(Int_t combined = 0) const; 
  //Gap?
  bool IsEB() const {return(gapinfo & 1<<0);}
  bool IsEE() const {return(gapinfo & 1<<1);}
  bool IsEBGap() const {return(gapinfo & 1<<2);}
  bool IsEBEtaGap() const {return(gapinfo & 1<<3);}
  bool IsEBPhiGap() const {return(gapinfo & 1<<4);}
  bool IsEEGap() const {return(gapinfo & 1<<5);}
  bool IsEERingGap() const {return(gapinfo & 1<<6);}
  bool IsEEDeeGap() const {return(gapinfo & 1<<7);}
  bool IsEBEEGap() const {return(gapinfo & 1<<8);}
  // Charge?
  bool IsGsfCtfChargeConsistent() const {return(chargeinfo & 1<<0);}
  bool IsGsfCtfScPixChargeConsistent() const {return(chargeinfo & 1<<1);}
  bool IsGsfScPixChargeConsistent() const {return(chargeinfo & 1<<2);}
};

class Vertex : public TVector3
{
 private:
  Double_t chi2;
  Double_t ndof;
  Double_t ptq;
  Int_t ntracks;
  TMatrixFSym covmatrix;
 public:
  Vertex() {}
  Vertex(Double_t X, Double_t Y, Double_t Z, Double_t Chi2, Double_t Ndof, Double_t Ptq, Int_t Ntracks, const Float_t* Cov);
  Double_t Chi2() const {return(chi2);}
  Double_t Ndof() const {return(ndof);}
  Double_t SumPtSquared() const {return(ptq);}
  Double_t Chi2OverNdof() const {return(chi2/ndof);}
  Int_t NTracks() const {return(ntracks);}
  Double_t XError() const {return(Sqrt(covmatrix(0,0)));}
  Double_t YError() const {return(Sqrt(covmatrix(1,1)));}
  Double_t ZError() const {return(Sqrt(covmatrix(2,2)));}
  TMatrixFSym CovMatrix() const {return(covmatrix);}
  Double_t CovMatrix(Int_t i, Int_t j) const {return(covmatrix(i,j));}
  bool IsGood() const {return(ndof > 4 && Z() < 24. && Perp() < 2. && !(chi2 == 0 && ntracks == 0 && ndof == 0));}
};

class Conversion : public TrackComposedParticle
{
 private:
  vector<TVector3> ecalpoints;
  Float_t mvaout;
 public:
  Conversion() {}
  Conversion(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t* Cov, Float_t mymvaout, TVector3 Ecalpoint1, TVector3 Ecalpoint2);
  TVector3 PositionAtEcal(UInt_t n) const {return(ecalpoints[n]);}
  Float_t MVAOut() const {return(mvaout);}
};

class Photon : public TLorentzVector, public TriggerObject
{
  friend class Analyse;
 private:
  Float_t e1x5;
		Float_t e2x5;
		Float_t e3x3;
		Float_t e5x5;
		Float_t sigmaetaeta;
		Float_t sigmaietaieta;
		Float_t ehcaloverecal;
		Float_t ehcaloverecaldepth1;
		Float_t ehcaloverecaldepth2;
		Float_t maxenergyxtal;
		Float_t isolationr3track;
		Float_t isolationr3trackhollow;
		UInt_t isolationr3ntrack;
		UInt_t isolationr3ntrackhollow;
		Float_t isolationr3ecal;
		Float_t isolationr3hcal;
		Float_t isolationr4track;
		Float_t isolationr4trackhollow;
		UInt_t isolationr4ntrack;
		UInt_t isolationr4ntrackhollow;
		Float_t isolationr4ecal;
		Float_t isolationr4hcal;
		UChar_t info;
		UInt_t gapinfo;
		vector<Conversion> conversions;
		void AddConversion(Conversion& newconversion);
		vector<SuperCluster> supercluster;
		void AddSC(SuperCluster& sc);
	public:
		Photon(){}
	 	Photon(const Analyse* ma, Double_t Px, Double_t Py, Double_t Pz, Float_t mye1x5, Float_t mye2x5, Float_t mye3x3, Float_t mye5x5, Float_t mysigmaetaeta, Float_t mysigmaietaieta, Float_t myehcaloverecal, Float_t myehcaloverecaldepth1, Float_t myehcaloverecaldepth2, Float_t mymaxenergyxtal, Float_t myisolationr3track, Float_t myisolationr3trackhollow, UInt_t myisolationr3ntrack, UInt_t myisolationr3ntrackhollow, Float_t myisolationr3ecal, Float_t myisolationr3hcal, Float_t myisolationr4track, Float_t myisolationr4trackhollow, UInt_t myisolationr4ntrack, UInt_t myisolationr4ntrackhollow, Float_t myisolationr4ecal, Float_t myisolationr4hcal, UChar_t myinfo, UInt_t mygapinfo);
		Float_t E1x5() const {return(e1x5);}
		Float_t E2x5() const {return(e2x5);}
		Float_t E3x3() const {return(e3x3);}
		Float_t E5x5() const {return(e5x5);}
		Float_t R9() const {return(e3x3/SCs(0).RawEnergy());}
		Float_t R25() const {return(e5x5/SCs(0).RawEnergy());}
		Float_t SigmaEtaEta() const {return(sigmaetaeta);}
		Float_t SigmaIEtaIEta() const {return(sigmaietaieta);}
		Float_t MaxEnergyXtal() const {return(maxenergyxtal);}
		Float_t EHCalOverECal() const {return(ehcaloverecal);}
		Float_t EHCalOverECalDepth1() const {return(ehcaloverecaldepth1);}
		Float_t EHCalOverECalDepth2() const {return(ehcaloverecaldepth2);}
		Float_t IsoR3Track() const {return(isolationr3track);}
		Float_t IsoR3TrackHollow() const {return(isolationr3trackhollow);}
		UInt_t IsoR3NTrack() const {return(isolationr3ntrack);}
		UInt_t IsoR3NTrackHollow() const {return(isolationr3ntrackhollow);}
		Float_t IsoR3ECal() const {return(isolationr3ecal);}
		Float_t IsoR3HCal() const {return(isolationr3hcal);}
		Float_t IsoR4Track() const {return(isolationr4track);}
		Float_t IsoR4TrackHollow() const {return(isolationr4trackhollow);}
		UInt_t IsoR4NTrack() const {return(isolationr4ntrack);}
		UInt_t IsoR4NTrackHollow() const {return(isolationr4ntrackhollow);}
		Float_t IsoR4ECal() const {return(isolationr4ecal);}
		Float_t IsoR4HCal() const {return(isolationr4hcal);}
		bool isPhoton() const {return((info & 1) > 0);}
		bool HasConversionTracks() const {return((info & 1<<1) > 0);}
		bool HasPixelSeed() const {return((info & 1<<2) > 0);}
		//Gap?
		bool IsEB() const {return((gapinfo & 1<<0) > 0);}
		bool IsEE() const {return((gapinfo & 1<<1) > 0);}
		bool IsEBGap() const {return((gapinfo & 1<<2) > 0);}
		bool IsEBEtaGap() const {return((gapinfo & 1<<3) > 0);}
		bool IsEBPhiGap() const {return((gapinfo & 1<<4) > 0);}
		bool IsEEGap() const {return((gapinfo & 1<<5) > 0);}
		bool IsEERingGap() const {return((gapinfo & 1<<6) > 0);}
		bool IsEEDeeGap() const {return((gapinfo & 1<<7) > 0);}
		bool IsEBEEGap() const {return((gapinfo & 1<<8) > 0);}
		//Conversions
		UInt_t NumConversions() const {return(conversions.size());}
		Conversion Conversions(UInt_t n) const {return(conversions[n]);}
		//Supercluster
		UInt_t NumSCs() const {return(supercluster.size());}
		SuperCluster SCs(UInt_t n) const {return(supercluster[n]);}
		//Trigger
		//Int_t Trigger(string triggername) const;//-1 trigger was not matched, 0 fired trigger not, 1 fired trigger 
};


//Not all Variables are filled for all Jet types.
class Jet : public TLorentzVector
{
	private:
		Float_t hadronicenergy;
		Float_t chargedhadronicenergy;
		Float_t emenergy;
		Float_t chargedemenergy;
		Int_t chargedmulti;
		Int_t neutralmulti;
		Float_t energycorr;
		Float_t btag[M_btagmax];
		Int_t n90;
		Int_t n60;
		Float_t fhpd;
		Float_t restrictedemf;
		Bool_t pujetsimpleloose;
		Bool_t pujetsimplemedium;
		Bool_t pujetsimpletight;
		Float_t pujetsimplemva;
		Bool_t pujetfullloose;
		Bool_t pujetfullmedium;
		Bool_t pujetfulltight;
		Float_t pujetfullmva;
		UInt_t svcount;
	public:
		Jet() {}
		Jet(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t Hadronicenergy, Float_t Chargedhadronicenergy, Float_t Emenergy, Float_t Chargedemenergy, Int_t Chargedmulti, Int_t Neutralmulti, Float_t Energycorr, const Float_t* Btag, Int_t N90, Int_t N60, Float_t Fhpd, Float_t Restrictedemf, Bool_t Pujetsimpleloose, Bool_t Pujetsimplemedium, Bool_t Pujetsimpletight, Float_t Pujetsimplemva, Bool_t Pujetfullloose, Bool_t Pujetfullmedium, Bool_t Pujetfulltight, Float_t Pujetfullmva, UInt_t Svcount);
		Float_t HadEnergyFraction() const {return(hadronicenergy/E());}
		Float_t HadEnergy() const {return(hadronicenergy);}
		Float_t EMEnergyFraction() const {return(emenergy/E());}
		Float_t EMEnergy() const {return(emenergy);}
		Float_t ChargedHadEnergyFraction() const {return(chargedhadronicenergy/E());}
		Float_t ChargedHadEnergy() const {return(chargedhadronicenergy);}
		Float_t ChargedEMEnergyFraction() const {return(chargedemenergy/E());}
		Float_t ChargedEMEnergy() const {return(chargedemenergy);}
		Int_t ChargedMulti() const {return(chargedmulti);}
		Int_t NeutralMulti() const {return(neutralmulti);}
		Float_t BTag(Int_t n) const;
		UInt_t N90() const {return(n90);} 
		UInt_t N60() const {return(n60);}
		Float_t FHPD() const {return(fhpd);}
		Float_t RestrictedEMF() const {return(restrictedemf);}
		Float_t JECRaw() const {return(energycorr);}
		//btag
		Float_t trackCountingHighPurBJetTags() const {return(btag[0]);}
		Float_t trackCountingHighEffBJetTags() const {return(btag[1]);}
		Float_t combinedSecondaryVertexBJetTags() const {return(btag[2]);}
		Float_t combinedSecondaryVertexMVABJetTags() const {return(btag[3]);}
		Float_t simpleSecondaryVertexHighPurBJetTags() const {return(btag[4]);}
		Float_t simpleSecondaryVertexHighEffBJetTags() const {return(btag[5]);}
		// Jet ID
		Bool_t puJetSimpleLoose() const { return pujetsimpleloose; }
		Bool_t puJetSimpleMedium() const { return pujetsimplemedium; }
		Bool_t puJetSimpleTight() const { return pujetsimpletight; }
		Float_t puJetSimpleMVA() const { return pujetsimplemva; }
		Bool_t puJetFullLoose() const { return pujetfullloose; }
		Bool_t puJetFullMedium() const { return pujetfullmedium; }
		Bool_t puJetFullTight() const { return pujetfulltight; }
		Float_t puJetFullMVA() const { return pujetfullmva; }
		// SVs
		UInt_t NumSVs() const { return svcount; }
};

class Tau : public TLorentzVector, public TriggerObject{
  friend class Analyse;
 private:
  Double_t dxy;
  Double_t dz;
  Float_t isolationneutralspt;
  bool    tau_L1trigger_match;  
  UInt_t  isolationneutralsnum;
  Float_t isolationchargedpt;
  UInt_t  isolationchargednum;
  Float_t isolationgammapt;
  UInt_t  isolationgammanum;
  Float_t iso3hits;
  Char_t  charge;
  ULong64_t  dishps;
  Float_t emfraction;
  Float_t newemfraction;
  Float_t hcaltotoverplead;
  Float_t hcal3x3overplead;
  Float_t ecalstripsumeoverplead;
  Float_t bremsrecoveryeoverplead;
  Float_t calocomp;
  Float_t segcomp;
  TLorentzVector jetp4;
  Jet jet;
  Float_t leadpfchargedhadrcandpt_;
  LorentzVector leadpfchargedhadrcandp4_;
  UInt_t   signalPFChargedHadrCands_size_;
  UInt_t   signalPFGammaCands_size_;

  vector<Track> tracks;
  const vector<string>* taudiscriminators;

public:
  Tau() {}
  Tau(const Analyse* ma, UInt_t n);
  Double_t Dxy() const { return(dxy); }
  Double_t Dz() const { return(dz); }
  bool    TauL1trigger_match() const {return(tau_L1trigger_match);}
  Float_t IsoNeutralsPt() const {return(isolationneutralspt);}
  UInt_t  IsoNeutralsNum() const {return(isolationneutralsnum);}
  Float_t IsoChargedPt() const {return(isolationchargedpt);}
  UInt_t  IsoChargedNum() const {return(isolationchargednum);}
  Float_t IsoGammaPt() const {return(isolationgammapt);}
  UInt_t  IsoGammaNum() const {return(isolationgammanum);}
  Float_t Iso3Hits() const { return(iso3hits);}
  Int_t   Charge() const {return(charge);}
  Float_t EMFraction() const {return(emfraction);}
  Float_t NewEMFraction() const {return(newemfraction);}
  Float_t HCalTotOverPLead() const {return(hcaltotoverplead);}
  Float_t HCal3x3OverPLead() const {return(hcal3x3overplead);}
  Float_t ECalStripSumEOverPLead() const {return(ecalstripsumeoverplead);}
  Float_t BremsRecoveryEOverPLead() const {return(bremsrecoveryeoverplead);}
  Float_t CaloComp() const {return(calocomp);}
  Float_t SegComp() const {return(segcomp);}
  Track   LeadingTrack() const;
  UInt_t  NumTracks() const {return(tracks.size());}
  Track   GetTrack(UInt_t num) const {return(tracks[num]);}	
  Float_t leadpfchargedhadrcandpt() const {return(leadpfchargedhadrcandpt_);}
  LorentzVector  leadpfchargedhadrcandp4() const {return(leadpfchargedhadrcandp4_);}
  UInt_t  signalPFChargedHadrCands_size() const {return(signalPFChargedHadrCands_size_);}
  UInt_t signalPFGammaCands_size() const {return(signalPFGammaCands_size_);}


  //TauDiscriminators
  Int_t   TauDiscriminator(string disname) const;
  const TLorentzVector& JetP4() const { return jetp4; }
  const Jet& AK5PFJet() const { return jet; }
};


class diTau : public TLorentzVector, public TriggerObject{
  friend class Analyse;
  friend class Tau;
  
 private:
  UInt_t ditau_Index_;
  UInt_t ditau_leg1_index_;
  UInt_t ditau_leg2_index_;
  TVector3 ditau_Vtx_;
  TVector3 ditau_VtxErr_;
  TVector3 ditau_reFitVtx_;
  TVector3 ditau_reFitVtxErr_;
  TVector3 ditau_leg1_;
  TVector3 ditau_leg1_OPV_;
  TVector3 ditau_leg2_;
  TVector3 ditau_leg2_OPV_;

  Int_t   ditau_reFitVtxRho_;
  Float_t ditau_reFitVtxNdof_;  
  Float_t ditau_leg1_dxy_;
  Float_t ditau_leg1_dz_;
  Float_t ditau_leg1_dxyErr_;
  Float_t ditau_leg1_dzErr_;
  Float_t ditau_leg1_dxy_OPV_;
  Float_t ditau_leg1_dz_OPV_;
  Float_t ditau_leg1_dxyErr_OPV_;
  Float_t ditau_leg1_dzErr_OPV_;
  Float_t ditau_leg2_dxy_;
  Float_t ditau_leg2_dz_;
  Float_t ditau_leg2_dxyErr_;
  Float_t ditau_leg2_dzErr_;
  Float_t ditau_leg2_dxy_OPV_;
  Float_t ditau_leg2_dz_OPV_;
  Float_t ditau_leg2_dxyErr_OPV_;
  Float_t ditau_leg2_dzErr_OPV_;

public:
  diTau() {}
  diTau(const Analyse* ma, UInt_t n);
  Int_t diTau_Leg1Index()   const { return(ditau_leg1_index_); }
  Int_t diTau_Leg2Index()   const { return(ditau_leg2_index_); }
  TVector3 diTau_Vtx()      const { return(ditau_Vtx_);         }
  TVector3 diTau_VtxErr()   const { return(ditau_VtxErr_);      }
  TVector3 diTau_reFitVtx() const {return( ditau_reFitVtx_); }
  TVector3 diTau_reFitVtxErr() const { return( ditau_reFitVtxErr_) ; }
  TVector3 diTau_Leg1_pos() const { return( ditau_leg1_); }
  TVector3 diTau_Leg1_pos_OPV() const { return(ditau_leg1_OPV_); }
  TVector3 diTau_Leg2_pos() const { return( ditau_leg2_); }
  TVector3 diTau_Leg2_pos_OPV() const { return(ditau_leg2_OPV_); }
  Float_t  diTau_Leg1_dxy() const { return(ditau_leg1_dxy_); }
  Float_t  diTau_Leg1_dz () const { return(ditau_leg1_dz_ ); }
  Float_t  diTau_Leg1_dxyErr() const { return(ditau_leg1_dxyErr_); }
  Float_t  diTau_Leg1_dzErr () const { return(ditau_leg1_dzErr_ ); }
  Float_t  diTau_Leg2_dxy() const { return(ditau_leg2_dxy_); }
  Float_t  diTau_Leg2_dz () const { return(ditau_leg2_dz_ ); }
  Float_t  diTau_Leg2_dxyErr() const { return(ditau_leg2_dxyErr_); }
  Float_t  diTau_Leg2_dzErr () const { return(ditau_leg2_dzErr_ ); }
};




class MuTauTauPair : public TLorentzVector{
  friend class Analyse;
private:
  Float_t leg1_px;
  Float_t leg1_py;
  Float_t leg1_pz;
  Float_t leg1_energy;
  Float_t leg2_px;
  Float_t leg2_py;
  Float_t leg2_pz;
  Float_t leg2_energy;
  Float_t mu_px;
  Float_t mu_py;
  Float_t mu_pz;
  Float_t mu_energy;
  Bool_t  svfit_int_valid;
  Float_t svfit_mass_int;
  Float_t svfit_mass_int_err_up;
  Float_t svfit_mass_int_err_down;
public:
  MuTauTauPair() {}
  MuTauTauPair(const Analyse* ma, UInt_t n);
  TLorentzVector Leg1() const { return(TLorentzVector( leg1_px, leg1_py, leg1_pz, leg1_energy));}
  TLorentzVector Leg2() const { return(TLorentzVector( leg2_px, leg2_py, leg2_pz, leg2_energy));}
  TLorentzVector Muon() const { return(TLorentzVector( mu_px, mu_py, mu_pz, mu_energy));}
  Bool_t         IsSVFitIntValid() const { return(svfit_int_valid);}
  Float_t        SVFitMassInt() const { return (svfit_mass_int);}
  Float_t        SVFitMassIntErrUp() const { return (svfit_mass_int_err_up); }
  Float_t        SVFitMassIntErrDown() const { return (svfit_mass_int_err_down); }
};

class ElTauTauPair : public TLorentzVector{
  friend class Analyse;
 private:
  Float_t leg1_px;
  Float_t leg1_py;
  Float_t leg1_pz;
  Float_t leg1_energy;
  Float_t leg2_px;
  Float_t leg2_py;
  Float_t leg2_pz;
  Float_t leg2_energy;
  Float_t el_px;
  Float_t el_py;
  Float_t el_pz;
  Float_t el_energy;
  Bool_t  svfit_int_valid;
  Float_t svfit_mass_int;
  Float_t svfit_mass_int_err_up;
  Float_t svfit_mass_int_err_down;
public:
  ElTauTauPair() {}
  ElTauTauPair(const Analyse* ma, UInt_t n);
  TLorentzVector Leg1() const { return(TLorentzVector( leg1_px, leg1_py, leg1_pz, leg1_energy));}
  TLorentzVector Leg2() const { return(TLorentzVector( leg2_px, leg2_py, leg2_pz, leg2_energy));}
  TLorentzVector Electron() const { return(TLorentzVector( el_px, el_py, el_pz, el_energy));}
  Bool_t         IsSVFitIntValid() const { return(svfit_int_valid);}
  Float_t        SVFitMassInt() const { return (svfit_mass_int);}
  Float_t        SVFitMassIntErrUp() const { return (svfit_mass_int_err_up); }
  Float_t        SVFitMassIntErrDown() const { return (svfit_mass_int_err_down); }
};

class BeamSpot : public TVector3
{
 private:
  TMatrixFSym covmatrix;
  Double_t xwidth;
  Double_t ywidth;
  Double_t zsigma;
  
 public:
  BeamSpot() {}
  BeamSpot(Double_t X, Double_t Y, Double_t Z, const Float_t* Cov, Double_t Xwidth, Double_t Ywidth, Double_t Zsigma);
  Double_t XError() const {return(Sqrt(covmatrix(0,0)));}
  Double_t YError() const {return(Sqrt(covmatrix(1,1)));}
  Double_t ZError() const {return(Sqrt(covmatrix(2,2)));}
  TMatrixFSym CovMatrix() const {return(covmatrix);}
  Double_t CovMatrix(Int_t i, Int_t j) const {return(covmatrix(i,j));}
  Double_t XWidth() const {return(xwidth);}
  Double_t YWidth() const {return(ywidth);}
  Double_t ZSigma() const {return(zsigma);}
};

class Luminosity
{
 private:
  Float_t lumival;
  Float_t lumierr;
  Float_t livefraction;
  Float_t deadfraction;
  UInt_t quality;
  UInt_t numevents;
  UInt_t numeventsorig;
  UInt_t hlttable;
  UInt_t l1techtable;
  UInt_t l1algotable;
  UInt_t counter;
 public:
  Luminosity(Float_t Lumival, Float_t Lumierr, Float_t Livefraction, Float_t Deadfraction, UInt_t Quality, UInt_t Numevents, UInt_t Numeventsprocessed, UInt_t Hlttable, UInt_t L1techtable, UInt_t L1algotable) : lumival(Lumival), lumierr(Lumierr), livefraction(Livefraction), deadfraction(Deadfraction), quality(Quality), numevents(Numevents), numeventsorig(Numeventsprocessed), hlttable(Hlttable), l1techtable(L1techtable), l1algotable(L1algotable), counter(0) {}
    Luminosity() {lumival = -1;}
    bool operator ==(const bool test) {bool val = true; if(lumival == -1){val = false;} return(test == val);}
    void operator ++() {counter++;}
    Luminosity& operator +=(Luminosity& other) {numevents += other.NumEvents(); numeventsorig += other.NumEventsOrig(); return(*this);}
    Float_t LumiValue() const {return(lumival);}
    Float_t LumiError() const {return(lumierr);}
    Float_t LiveFraction() const {return(livefraction);}
    Float_t DeadFraction() const {return(deadfraction);}
    UInt_t Quality() const {return(quality);}
    UInt_t NumEvents() const {return(numevents);}
    UInt_t NumEventsOrig() const {return(numeventsorig);}
    UInt_t NumEventsProcessed() const {return(counter);}
    Float_t ProcessedLumi() const {if(numevents != 0) {return(lumival*counter/numevents);}else{return(lumival);}}
    Float_t ProcessedFraction() const {if(numevents != 0){return(Float_t(counter)/numevents);}else{return(1.);}}
    UInt_t HLTTable() const {return(hlttable);}
    UInt_t L1TechTable() const {return(l1techtable);}
    UInt_t L1AlgoTable() const {return(l1algotable);}
    void Value(Float_t myvalue){lumival = myvalue;}
    void ValueErr(Float_t myvalueerr){lumierr = myvalueerr;}
    void LiveFrac(Float_t mylivefrac){livefraction = mylivefrac;}
    void DeadFrac(Float_t mydeadfrac){deadfraction = mydeadfrac;}
    void Quality(UInt_t myquality){quality = myquality;}
    void EventsFiltered(UInt_t myeventsfiltered){numevents = myeventsfiltered;}
    void EventsOriginal(UInt_t myeventsorig){numeventsorig = myeventsorig;}
};

class RunInfo
{
	private:
		UInt_t runnumber;
		vector<string> hltnames;
		Int_t hltnum;
		Int_t hlttablesnum;
		vector<vector<UInt_t> > hltprescales;
		Int_t l1algonum;
		Int_t l1algotablesnum;
		vector<vector<UInt_t> > l1algoprescales;
		Int_t l1technum;
		Int_t l1techtablesnum;
		vector<vector<UInt_t> > l1techprescales;
		vector<string> hltmunames;
		vector<string> hltelnames;
		vector<string> hltphotonnames;
		vector<string> hlttaunames;
		vector<string> hltjetnames;
		vector<string> taudiscriminators;
	public:
		RunInfo(UInt_t Number, UInt_t Hltcount, string Hltnames, string Hltmunames, string Hltelnames, string Hltphotonnames, string Hlttaunames, string Hltjetnames, string Taudiscriminators, UInt_t Hltprescaletablescount, UInt_t* Hltprescaletables, UInt_t L1algocount, UInt_t L1algoprescaletablescount, UInt_t* L1algoprescaletables, UInt_t L1techcount, UInt_t L1techprescaletablescount, UInt_t* L1techprescaletables);
		RunInfo() : runnumber(0) {}
		UInt_t Run() const {return(runnumber);}
		UInt_t RunNumber() const {return(runnumber);}
		UInt_t NumHLTTables() const {return(hlttablesnum);} 
		UInt_t NumHLT() const {return(hltnum);}
		string HLTAllNames() const;
		string HLTMuonAllNames() const;
		string HLTElectronAllNames() const;
		string HLTTauAllNames() const;
		string HLTPhotonAllNames() const;
		string HLTJetAllNames() const;
		string TauDiscriminatorsAllNames() const;
		UInt_t HLTPrescale(UInt_t trigger, UInt_t table) const {if(trigger < NumHLT() && table < NumHLTTables()){return(hltprescales[trigger][table]);}else{return(1);}}
		string HLTName(UInt_t trigger) const {return(hltnames[trigger]);}
		vector<string> MatchTriggerNames(string name);
		Int_t HLTIndex(string name) const;
		UInt_t NumL1TechTables() const {return(l1techtablesnum);} 
		UInt_t NumL1Tech() const {return(l1technum);}
		UInt_t L1TechPrescale(UInt_t trigger, UInt_t table) const {if(trigger < NumL1Tech() && table < NumL1TechTables()){return(l1techprescales[trigger][table]);}else{return(1);}}
		UInt_t NumL1AlgoTables() const {return(l1algotablesnum);}
		UInt_t NumL1Algo() const {return(l1algonum);}
		UInt_t L1AlgoPrescale(UInt_t trigger, UInt_t table) const {if(trigger < NumL1Algo() && table < NumL1AlgoTables()){return(l1algoprescales[trigger][table]);}else{return(1);}}
		const vector<string>* GetHLTMuonNames() const {return(&hltmunames);}
		const vector<string>* GetHLTElectronNames() const {return(&hltelnames);}
		const vector<string>* GetHLTTauNames() const {return(&hlttaunames);}
		const vector<string>* GetHLTPhotonNames() const {return(&hltphotonnames);}
		const vector<string>* GetHLTJetNames() const {return(&hltjetnames);}
		const vector<string>* GetTauDiscriminators() const {return(&taudiscriminators);}
};

class TriggerLumi
{
	private:
		Int_t index;
		UInt_t prescale;
	public:
		TriggerLumi() : index(-1), prescale(0) {}
		TriggerLumi(Int_t Index, UInt_t Prescale) : index(Index), prescale(Prescale) {}
		Int_t Index() const {return(index);}
		UInt_t Prescale() const {return(prescale);}
};

class TriggerRun
{
	private:
		Float_t lumi;
		map< UInt_t, TriggerLumi> lumiinfo;
	public:
		TriggerRun() : lumi(-1.) {}
		void Lumi(Float_t runlumi) {lumi = runlumi;}
		void SetBlock(UInt_t block, Int_t index, UInt_t prescale);
		TriggerLumi GetBlock(UInt_t block);
		Float_t Lumi() const {return(lumi);}
		map< UInt_t, TriggerLumi>::iterator Begin() {return(lumiinfo.begin());} 
		map< UInt_t, TriggerLumi>::iterator End() {return(lumiinfo.end());} 
};

class TriggerSelection
{
	private:
		Analyse* AN;
		map<UInt_t, TriggerRun> runinfo;
	public:
		TriggerSelection(Analyse* an, vector<string> names, bool useprescaled = false);
		Int_t Result();
		Float_t LumiUsed(Int_t format = 0);
		Float_t LumiBeforeEvent();
		string GetTriggerName(UInt_t run = 0, UInt_t lumiblock = 0);
		void PrintInfo();
};

#endif
