#include "DataForm.h"
#include "Analyse.h"
#include <Math/VectorUtil.h>

GenLightParticle::GenLightParticle(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Info) :
TLorentzVector(Px, Py, Pz, E),
vertex(X, Y, Z),
status(Status),
pdgid(Pdgid),
info(Info)
{
}

GenParticle::GenParticle(const Analyse* Myanalyse, UInt_t Myindex, Double_t E, Double_t Px, Double_t Py, Double_t Pz, 
			 Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Motherfirst, UInt_t Mothernum, 
			 UInt_t Daughterfirst, UInt_t Daughternum, Int_t charge, LorentzVector p4) :
  TLorentzVector(Px, Py, Pz, E),
  myanalyse(Myanalyse),
  vertex(X, Y, Z),
  status(Status),
  pdgid(Pdgid),
  myindex(Myindex),
  motherfirst(Motherfirst),
  mothernum(Mothernum),
  daughterfirst(Daughterfirst),
  daughternum(Daughternum),
  charge_(charge),
  p4_(p4)
{
}

UInt_t GenParticle::GetMotherIndex(UInt_t num) const
{
        return(myanalyse->genallparticles_mothers[motherfirst+num]);
}

UInt_t GenParticle::GetDaughterIndex(UInt_t num) const
{
        return(myanalyse->genallparticles_daughters[daughterfirst+num]);
}

GenParticle GenParticle::GetMother(UInt_t num) const
{
        if(num >= mothernum){cout << "GenParticle::GetMother: mother index out of range. return this." << endl; return(*this);}
        return(myanalyse->AllGenParticles(GetMotherIndex(num)));
}


GenParticle GenParticle::GetDaughter(UInt_t num) const
{
  if(num >= daughternum){cout << "GenParticle::GetDaughter: daughter index out of range. return this." << endl; return(*this);}
  return(myanalyse->AllGenParticles(GetDaughterIndex(num)));
}

bool GenParticle::HasAnyMotherPDGId(Int_t pdgid, bool antiparticle)
{
vector<UInt_t> visited;
return(HasAnyMotherPDGId(visited, pdgid, antiparticle));
}

bool GenParticle::HasAnyMotherPDGId(vector<UInt_t>& visited, Int_t pdgid, bool antiparticle)
{
	for(UInt_t i = 0 ; i < NumMothers() ; i++)
	{
		bool isvisited = false;
		for(UInt_t u = 0 ; u < visited.size() ; u++)
		{
			if(visited[u] == GetMotherIndex(i)){isvisited = true; break;}
		}
		if(isvisited) continue;
		visited.push_back(GetMotherIndex(i));
		if(antiparticle && Abs(pdgid) == GetMother(i).PDGId()) return(true); 
		if(!antiparticle && pdgid == GetMother(i).PDGId()) return(true); 
		if(GetMother(i).HasAnyMotherPDGId(visited, pdgid, antiparticle)) return(true);
	}
	return(false);
}

TrackComposedParticle::TrackComposedParticle(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t* Cov) :
	TVector3(Vx,Vy,Vz),
	covmatrix(3),
	chi2(Chi2),
	ndof(Ndof)
{
	covmatrix(0,0) = Cov[0];
	covmatrix(0,1) = Cov[1];
	covmatrix(0,2) = Cov[2];
	covmatrix(1,0) = Cov[1];
	covmatrix(1,1) = Cov[3];
	covmatrix(1,2) = Cov[4];
	covmatrix(2,0) = Cov[2];
	covmatrix(2,1) = Cov[4];
	covmatrix(2,2) = Cov[5];
}



void TrackComposedParticle::AddTrack(Track& newtrack)
{
tracks.push_back(newtrack);
}

Double_t TrackComposedParticle::VertexSig3D(Vertex& vertex) const
{
TVector3 dist = *this - vertex;
return(dist.Mag2()/Sqrt(((vertex.CovMatrix()+CovMatrix())* dist)*dist));
}
		
Double_t TrackComposedParticle::VertexSig2D(Vertex& vertex) const
{
TVector3 dist = *this - vertex;
dist.SetZ(0);
return(dist.Mag2()/Sqrt(((vertex.CovMatrix()+CovMatrix())* dist)*dist));
}

// TriggerObject::TriggerObject(const Analyse* ma, const vector<string>* Triggernames, UInt_t Trigger, bool isL1triggerd)
// {  
//   ma = MA;
//   Trigger =  trigger;
//   Triggernames = triggernames;
//   isL1triggerd = false;
// }

TriggerObject::TriggerObject(const Analyse* ma, const vector<string>* Triggernames, UInt_t Trigger) : 
MA(ma),
trigger(Trigger),
triggernames(Triggernames)
{

}


Int_t TriggerObject::Trigger(string triggername) const
{
	Int_t index = -1;
	for(UInt_t i = 0 ; i < triggernames->size() ; i++)
	{
		std::string indexedTrigger = (*triggernames)[i];
		std::string::size_type pos = indexedTrigger.find(':');
		if(pos != std::string::npos) indexedTrigger.erase(pos);

		if(triggername == indexedTrigger)
		{
			index = i;
			break;
		}
	}

	if(index == -1) return(0);
	bool result = (trigger & 1<<index) != 0;
	size_t pos = triggername.find(":");
	if(pos != string::npos)
	{
		triggername = triggername.substr(0, pos);
	}
	Int_t glindex = MA->GetHLTriggerIndex(triggername);
	if(MA->GetHLTPrescale(glindex) == 1)
	{
		if(result) return(1);
		else return(-1);
	}
	if(result) return(2);
	else return(-2);
}



Muon::Muon(const Analyse* ma, Double_t Px, Double_t Py, Double_t Pz, Float_t Pterror, Float_t Chi2, Float_t Ndof, Float_t Pfisolationr4_sumchargedhadronpt, Float_t Pfisolationr4_sumchargedparticlept, Float_t Pfisolationr4_sumneutralhadronet,Float_t Pfisolationr4_sumphotonet, Float_t Pfisolationr4_sumPUpt, Float_t Isolationr3track, Int_t Isolationr3ntrack, Float_t Isolationr3ecal, Float_t Isolationr3hcal, Float_t Ecalenergy, Float_t Hcalenergy, Int_t Charge, Int_t Numchambers, Int_t Numchamberswithsegments, Int_t Numchamberhits, UInt_t Type, UInt_t Trigger, UInt_t Trackermuonquality, Float_t MVAId, Float_t MVAIso) :
	TLorentzVector(Px, Py, Pz, sqrt(Px*Px+Py*Py+Pz*Pz+MuonMassQ)),
	TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTMuonNames(), Trigger),
	pterror(Pterror),
	chi2(Chi2),
	ndof(Ndof),
	isolationr3track(Isolationr3track),
	isolationr3ntrack(Isolationr3ntrack),
	isolationr3ecal(Isolationr3ecal),
	isolationr3hcal(Isolationr3hcal),
	pfisolationr4_sumchargedhadronpt(Pfisolationr4_sumchargedhadronpt),
	pfisolationr4_sumchargedparticlept(Pfisolationr4_sumchargedparticlept),
	pfisolationr4_sumneutralhadronet(Pfisolationr4_sumneutralhadronet),
	pfisolationr4_sumphotonet(Pfisolationr4_sumphotonet),
	pfisolationr4_sumPUpt(Pfisolationr4_sumPUpt),
	ecalenergy(Ecalenergy),
	hcalenergy(Hcalenergy),
	charge(Charge),
	numchambers(Numchambers),
	numchamberswithsegments(Numchamberswithsegments),
	numchamberhits(Numchamberhits),
	type(Type),
	trackermuonquality(Trackermuonquality),
	mva_id(MVAId),
	mva_iso(MVAIso)
{
}

Int_t Muon::NumStations() const
{
	Int_t result(0);
	if(trackermuonquality & 1<<24) result++;
	if(trackermuonquality & 1<<25) result++;
	if(trackermuonquality & 1<<26) result++;
	if(trackermuonquality & 1<<27) result++;
	if(trackermuonquality & 1<<28) result++;
	if(trackermuonquality & 1<<29) result++;
	if(trackermuonquality & 1<<30) result++;
	if(trackermuonquality & 1<<31) result++;
	return(result);
}


 


Electron::Electron(const Analyse* ma, Double_t Px, Double_t Py, Double_t Pz, Float_t Trackchi2, Float_t Trackndof, Double_t Outerx, Double_t Outery, Double_t Outerz, Double_t Closestx, Double_t Closesty, Double_t Closestz, Double_t Dxy, Double_t Dxyerr, Double_t Dz, Double_t Dzerr, Float_t Esuperclusterovertrack, Float_t Eseedclusterovertrack, Float_t Deltaetasuperclustertrack, Float_t Deltaphisuperclustertrack, Float_t E1x5, Float_t E2x5, Float_t E5x5, Float_t Sigmaetaeta, Float_t Sigmaietaieta, Float_t Ehcaloverecal, Float_t Ehcaloverecaldepth1, Float_t Ehcaloverecaldepth2, Float_t Isolationr3track, Float_t Isolationr3ecal, Float_t Isolationr3hcal, Float_t Isolationr4track, Float_t Isolationr4ecal, Float_t Isolationr4hcal,Float_t Pfisolationr3_sumchargedhadronpt, Float_t Pfisolationr3_sumchargedparticlept, Float_t Pfisolationr3_sumneutralhadronet, Float_t Pfisolationr3_sumphotonet, Float_t Pfisolationr3_sumPUpt, Float_t  Pfisolationr4_sumchargedhadronpt, Float_t Pfisolationr4_sumchargedparticlept, Float_t Pfisolationr4_sumneutralhadronet, Float_t Pfisolationr4_sumphotonet, Float_t Pfisolationr4_sumPUpt, Int_t Nhits, Int_t Nmissinghits, Int_t Npixelhits, Int_t Npixellayers, Int_t Nstriplayers, Float_t Convdist, Float_t Convdcot, Float_t Convradius, UInt_t Gapinfo, UInt_t Chargeinfo, Float_t Fbrems, Int_t Numbrems, Int_t Charge, Byte_t Info, UInt_t Trigger, Float_t MVAIdNontrig, Bool_t HasConversion) :
	TLorentzVector(Px, Py, Pz, sqrt(Px*Px+Py*Py+Pz*Pz+ElectronMassQ)),
	TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTElectronNames(), Trigger),
	outerpoint(Outerx, Outery, Outerz),
	closestpoint(Closestx, Closesty, Closestz),
	dxy(Dxy),
	dxyerr(Dxyerr),
	dz(Dz),
	dzerr(Dzerr),
	esuperclusterovertrack(Esuperclusterovertrack),
	eseedclusterovertrack(Eseedclusterovertrack),
	deltaetasuperclustertrack(Deltaetasuperclustertrack),
	deltaphisuperclustertrack(Deltaphisuperclustertrack),
	e1x5(E1x5),
	e2x5(E2x5),
	e5x5(E5x5),
	sigmaetaeta(Sigmaetaeta),
	sigmaietaieta(Sigmaietaieta),
	ehcaloverecal(Ehcaloverecal),
	ehcaloverecaldepth1(Ehcaloverecaldepth1),
	ehcaloverecaldepth2(Ehcaloverecaldepth2),
	isolationr3track(Isolationr3track),
	isolationr3ecal(Isolationr3ecal),
	isolationr3hcal(Isolationr3hcal),
	isolationr4track(Isolationr4track),
	isolationr4ecal(Isolationr4ecal),
	isolationr4hcal(Isolationr4hcal),
	pfisolationr3_sumchargedhadronpt(Pfisolationr3_sumchargedhadronpt),
	pfisolationr3_sumchargedparticlept(Pfisolationr3_sumchargedparticlept),
	pfisolationr3_sumneutralhadronet(Pfisolationr3_sumneutralhadronet),
	pfisolationr3_sumphotonet(Pfisolationr3_sumphotonet),
	pfisolationr3_sumPUpt(Pfisolationr3_sumPUpt),
	pfisolationr4_sumchargedhadronpt(Pfisolationr4_sumchargedhadronpt),
	pfisolationr4_sumchargedparticlept(Pfisolationr4_sumchargedparticlept),
	pfisolationr4_sumneutralhadronet(Pfisolationr4_sumneutralhadronet),
	pfisolationr4_sumphotonet(Pfisolationr4_sumphotonet),
	pfisolationr4_sumPUpt(Pfisolationr4_sumPUpt),
	trackchi2(Trackchi2),
	trackndof(Trackndof),
	nhits(Nhits),
	nmissinghits(Nmissinghits),
	npixelhits(Npixelhits),
	npixellayers(Npixellayers),
	nstriplayers(Nstriplayers),
	convdist(Convdist),
	convdcot(Convdcot),
	convradius(Convradius),
	gapinfo(Gapinfo),
	chargeinfo(Chargeinfo),
	fbrems(Fbrems),
	numbrems(Numbrems),
	charge(Charge),
	info(Info),
	mvaidnontrig(MVAIdNontrig),
	hasconversion(HasConversion)
{
}

void Electron::AddSC(SuperCluster& sc) {supercluster.push_back(sc);}

bool Electron::WorkingPoint(Int_t Missinghits, Float_t Convdist, Float_t Convdcot, Float_t Sigmaietaieta, Float_t Deltaphisctrack, Float_t Deltaetasctrack, Float_t Ehcaloverecal) const 
{
	return(NMissingHits() <= Missinghits && !(Abs(ConversionDist()) < Convdist && Abs(ConversionDCot()) < Convdcot) && SigmaIEtaIEta() < Sigmaietaieta && DeltaEtaSuperClusterTrack() < Deltaetasctrack && DeltaPhiSuperClusterTrack() < Deltaphisctrack && EHcalOverECal() < Ehcaloverecal);
}

bool Electron::WP95_v1(Int_t combined) const
{
	if(combined == 1)
	{
		if(IsEE())
		{
			return(WorkingPoint(1, -1., -1., 0.03, 0.7, 0.01, 0.07) && IsoR3CombinedRel() < 0.1);
		}
		else if(IsEB())
		{
			return(WorkingPoint(1, -1., -1., 0.01, 0.8, 0.007, 0.15) && IsoR3CombinedRel() < 0.15);
		}
	}
	else
	{
		if(IsEE())
		{
			return(WorkingPoint(1, -1., -1., 0.03, 0.7, 0.01, 0.07) && IsoR3TrackRel() < 0.08 && IsoR3ECalRel() < 0.06 && IsoR3HCalRel() < 0.05);
		}
		else if(IsEB())
		{
			return(WorkingPoint(1, -1., -1., 0.01, 0.8, 0.007, 0.15) && IsoR3TrackRel() < 0.15 && IsoR3ECalRel() < 2. && IsoR3HCalRel() < 0.12);
		}
	}
	return(false);
}

bool Electron::WP90_v1(Int_t combined) const
{
	if(combined == 1)
	{
		if(IsEE())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.7, 0.009, 0.05) && IsoR3CombinedRel() < 0.07);
		}
		else if(IsEB())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.8, 0.007, 0.12) && IsoR3CombinedRel() < 0.1);
		}
	}
	else
	{
		if(IsEE())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.7, 0.009, 0.05) && IsoR3TrackRel() < 0.05 && IsoR3ECalRel() < 0.06 && IsoR3HCalRel() < 0.03);
		}
		else if(IsEB())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.8, 0.007, 0.12) && IsoR3TrackRel() < 0.12 && IsoR3ECalRel() < 0.09 && IsoR3HCalRel() < 0.1);
		}
	}
	return(false);
}

bool Electron::WP85_v1(Int_t combined) const
{
	if(combined == 1)
	{
		if(IsEE())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.04, 0.007, 0.025) && IsoR3CombinedRel() < 0.06);
		}
		else if(IsEB())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.06, 0.006, 0.04) && IsoR3CombinedRel() < 0.09);
		}
	}
	else
	{
		if(IsEE())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.04, 0.007, 0.025) && IsoR3TrackRel() < 0.05 && IsoR3ECalRel() < 0.05 && IsoR3HCalRel() < 0.025);
		}
		else if(IsEB())
		{
			return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.06, 0.006, 0.04) && IsoR3TrackRel() < 0.09 && IsoR3ECalRel() < 0.08 && IsoR3HCalRel() < 0.1);
		}
	}
	return(false);
}

bool Electron::WP80_v1(Int_t combined) const
{
	if(combined == 1)
	{
		if(IsEE())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.03, 0.007, 0.025) && IsoR3CombinedRel() < 0.06);
		}
		else if(IsEB())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.06, 0.004, 0.04) && IsoR3CombinedRel() < 0.07);
		}
	}
	else
	{
		if(IsEE())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.03, 0.007, 0.025) && IsoR3TrackRel() < 0.04 && IsoR3ECalRel() < 0.05 && IsoR3HCalRel() < 0.025);
		}
		else if(IsEB())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.06, 0.004, 0.04) && IsoR3TrackRel() < 0.09 && IsoR3ECalRel() < 0.07 && IsoR3HCalRel() < 0.1);
		}
	}
	return(false);
}

bool Electron::WP70_v1(Int_t combined) const
{
	if(combined == 1)
	{
		if(IsEE())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3CombinedRel() < 0.03);
		}
		else if(IsEB())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.03, 0.004, 0.025) && IsoR3CombinedRel() < 0.04);
		}
	}
	else
	{
		if(IsEE())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3TrackRel() < 0.025 && IsoR3ECalRel() < 0.025 && IsoR3HCalRel() < 0.02);
		}
		else if(IsEB())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.03, 0.004, 0.025) && IsoR3TrackRel() < 0.05 && IsoR3ECalRel() < 0.06 && IsoR3HCalRel() < 0.03);
		}
	}
	return(false);
}

bool Electron::WP60_v1(Int_t combined) const
{
	if(combined == 1)
	{
		if(IsEE())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3CombinedRel() < 0.02);
		}
		else if(IsEB())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.025, 0.004, 0.025) && IsoR3CombinedRel() < 0.03);
		}
	}
	else
	{
		if(IsEE())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3TrackRel() < 0.025 && IsoR3ECalRel() < 0.02 && IsoR3HCalRel() < 0.02);
		}
		else if(IsEB())
		{
			return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.025, 0.004, 0.025) && IsoR3TrackRel() < 0.04 && IsoR3ECalRel() < 0.04 && IsoR3HCalRel() < 0.03);
		}
	}
	return(false);
}

Conversion::Conversion(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t* Cov, Float_t mymvaout, TVector3 Ecalpoint1, TVector3 Ecalpoint2) : TrackComposedParticle(Vx, Vy, Vz, Chi2, Ndof, Cov),
mvaout(mymvaout)
{
	ecalpoints.push_back(Ecalpoint1);
	ecalpoints.push_back(Ecalpoint2);
}

Photon::Photon(const Analyse* ma, Double_t Px, Double_t Py, Double_t Pz, Float_t mye1x5, Float_t mye2x5, Float_t mye3x3, Float_t mye5x5, Float_t mysigmaetaeta, Float_t mysigmaietaieta, Float_t myehcaloverecal, Float_t myehcaloverecaldepth1, Float_t myehcaloverecaldepth2, Float_t mymaxenergyxtal, Float_t myisolationr3track, Float_t myisolationr3trackhollow, UInt_t myisolationr3ntrack, UInt_t myisolationr3ntrackhollow, Float_t myisolationr3ecal, Float_t myisolationr3hcal, Float_t myisolationr4track, Float_t myisolationr4trackhollow, UInt_t myisolationr4ntrack, UInt_t myisolationr4ntrackhollow, Float_t myisolationr4ecal, Float_t myisolationr4hcal, UChar_t myinfo, UInt_t mygapinfo) :
	TLorentzVector(Px, Py, Pz, sqrt(Px*Px+Py*Py+Pz*Pz)),
	TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTPhotonNames(), 0),
	e1x5(mye1x5),
	e2x5(mye2x5),
	e3x3(mye3x3),
	e5x5(mye5x5),
	sigmaetaeta(mysigmaetaeta),
	sigmaietaieta(mysigmaietaieta),
	ehcaloverecal(myehcaloverecal),
	ehcaloverecaldepth1(myehcaloverecaldepth1),
	ehcaloverecaldepth2(myehcaloverecaldepth2),
	maxenergyxtal(mymaxenergyxtal),
	isolationr3track(myisolationr3track),
	isolationr3trackhollow(myisolationr3trackhollow),
	isolationr3ntrack(myisolationr3ntrack),
	isolationr3ntrackhollow(myisolationr3ntrackhollow),
	isolationr3ecal(myisolationr3ecal),
	isolationr3hcal(myisolationr3hcal),
	isolationr4track(myisolationr4track),
	isolationr4trackhollow(myisolationr4trackhollow),
	isolationr4ntrack(myisolationr4ntrack),
	isolationr4ntrackhollow(myisolationr4ntrackhollow),
	isolationr4ecal(myisolationr4ecal),
	isolationr4hcal(myisolationr4hcal),
	info(myinfo),
	gapinfo(mygapinfo)
{
}

void Photon::AddConversion(Conversion& newconversion) {conversions.push_back(newconversion);}
void Photon::AddSC(SuperCluster& sc) {supercluster.push_back(sc);}

Jet::Jet(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t Hadronicenergy, Float_t Chargedhadronicenergy, Float_t Emenergy, Float_t Chargedemenergy, Int_t Chargedmulti, Int_t Neutralmulti, Float_t Energycorr, const Float_t* Btag, Int_t N90, Int_t N60, Float_t Fhpd, Float_t Restrictedemf, Bool_t Pujetsimpleloose, Bool_t Pujetsimplemedium, Bool_t Pujetsimpletight, Float_t Pujetsimplemva, Bool_t Pujetfullloose, Bool_t Pujetfullmedium, Bool_t Pujetfulltight, Float_t Pujetfullmva, UInt_t Svcount) :
	TLorentzVector(Px, Py, Pz, E),
	hadronicenergy(Hadronicenergy),
	chargedhadronicenergy(Chargedhadronicenergy),
	emenergy(Emenergy),
	chargedemenergy(Chargedemenergy),
	chargedmulti(Chargedmulti),
	neutralmulti(Neutralmulti),
	energycorr(Energycorr),
	n90(N90),
	n60(N60),
	fhpd(Fhpd),
	restrictedemf(Restrictedemf),
	pujetsimpleloose(Pujetsimpleloose),
	pujetsimplemedium(Pujetsimplemedium),
	pujetsimpletight(Pujetsimpletight),
	pujetsimplemva(Pujetsimplemva),
	pujetfullloose(Pujetfullloose),
	pujetfullmedium(Pujetfullmedium),
	pujetfulltight(Pujetfulltight),
	pujetfullmva(Pujetfullmva),
	svcount(Svcount)
{
	if(Btag == 0)
	{
		for(int i = 0 ; i < M_btagmax ; i++)
		{
			btag[i] = -1000000.; 
		}
	}
	else
	{
		for(int i = 0 ; i < M_btagmax ; i++)
		{
			btag[i] = Btag[i]; 
		}
	}
}

Float_t Jet::BTag(Int_t n) const 
{
	if(n< M_btagmax)
	{
		return(btag[n]);
	}
	else
	{
		cerr << "wrong btag number" << endl;
		return(-1);
	}
}


Track::Track(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Double_t Ox, Double_t Oy, Double_t Oz, Double_t Cx, Double_t Cy, Double_t Cz, Double_t Chi2, Double_t Ndof, Double_t Dxy, Double_t Dxyerr, Double_t Dz, Double_t Dzerr, Int_t Charge, Int_t Nhits, Int_t Nmissinghits, Int_t Npixelhits, Int_t Npixellayers, Int_t Nstriplayers, Float_t Dedxharmonic2) :
	TLorentzVector(Px, Py, Pz, E),
	outerpoint(Ox, Oy, Oz),
	closestpoint(Cx, Cy, Cz),
	chi2(Chi2),
	ndof(Ndof),
	dxy(Dxy),
	dxyerr(Dxyerr),
	dz(Dz),
	dzerr(Dzerr),
	charge(Charge),
	nhits(Nhits),
	nmissinghits(Nmissinghits),
	npixelhits(Npixelhits),
	npixellayers(Npixellayers),
	nstriplayers(Nstriplayers),
	dedxharmonic2(Dedxharmonic2)
{
}

EcalHit::EcalHit(Double_t Energy, Double_t X, Double_t Y, Double_t Z) : TVector3(X, Y, Z), energy(Energy)
{
}

Cluster::Cluster(Double_t E, Double_t x, Double_t y, Double_t z, Int_t Size) :
	TLorentzVector(E*x/sqrt(x*x+y*y+z*z), E*y/sqrt(x*x+y*y+z*z), E*z/sqrt(x*x+y*y+z*z), E),
	position(x,y,z),
	size(Size)
{
}

void Cluster::AddHit(EcalHit& hit){hits.push_back(hit);}

SuperCluster::SuperCluster(Double_t E, Double_t x, Double_t y, Double_t z, Float_t Rawe, Float_t Phiwidth, Float_t Etawidth) :
	TLorentzVector(E*x/sqrt(x*x+y*y+z*z), E*y/sqrt(x*x+y*y+z*z), E*z/sqrt(x*x+y*y+z*z), E),
	position(x,y,z),
	rawe(Rawe),
	phiwidth(Phiwidth),
	etawidth(Etawidth)
{
}

void SuperCluster::AddCluster(Cluster& newcluster)
{
	clusters.push_back(newcluster);
}

void SuperCluster::AddESCluster(Cluster& newcluster)
{
	esclusters.push_back(newcluster);
}

Tau::Tau(const Analyse* ma, UInt_t n) :
  TLorentzVector(ma->tau_px[n], ma->tau_py[n], ma->tau_pz[n], ma->tau_e[n] >= 0 ? ma->tau_e[n] : sqrt(ma->tau_px[n]*ma->tau_px[n]+ma->tau_py[n]*ma->tau_py[n]+ma->tau_pz[n]*ma->tau_pz[n]+TauMassQ)),
  TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTTauNames(), ma->tau_trigger[n]),
  dxy(ma->tau_dxy[n]), dz(ma->tau_dz[n]),
  tau_L1trigger_match(ma->tau_L1trigger_match[n]),
  isolationneutralspt(ma->tau_isolationneutralspt[n]),
  isolationneutralsnum(ma->tau_isolationneutralsnum[n]),
  isolationchargedpt(ma->tau_isolationchargedpt[n]),
  isolationchargednum(ma->tau_isolationchargednum[n]),
  isolationgammapt(ma->tau_isolationgammapt[n]),
  isolationgammanum(ma->tau_isolationgammanum[n]),
  iso3hits(ma->tau_bycombinedisolationdeltabetacorrraw3hits[n]),
  charge(ma->tau_charge[n]),
  dishps(ma->tau_dishps[n]),
  emfraction(ma->tau_emfraction[n]),
  newemfraction(ma->tau_newemfraction[n]),
  hcaltotoverplead(ma->tau_hcaltotoverplead[n]),
  hcal3x3overplead(ma->tau_hcal3x3overplead[n]),
  ecalstripsumeoverplead(ma->tau_ecalstripsumeoverplead[n]),
  bremsrecoveryeoverplead(ma->tau_bremsrecoveryeoverplead[n]),
  calocomp(ma->tau_calocomp[n]),
  segcomp(ma->tau_segcomp[n]),
  leadpfchargedhadrcandpt_(ma->tau_leadpfchargedhadrcandpt[n]),
  leadpfchargedhadrcandp4_(ma->tau_leadpfchargedhadrcandp4[n]),
  signalPFChargedHadrCands_size_(ma->tau_signalPFChargedHadrCands_size[n]),
  signalPFGammaCands_size_(ma->tau_signalPFGammaCands_size[n]),

/*  jet(ma->tau_ak5pfjet_e[n], ma->tau_ak5pfjet_px[n], 
      ma->tau_ak5pfjet_py[n], ma->tau_ak5pfjet_pz[n], 
      ma->tau_ak5pfjet_hadronicenergy[n], 
      ma->tau_ak5pfjet_chargedhadronicenergy[n], 
      ma->tau_ak5pfjet_emenergy[n], ma->tau_ak5pfjet_chargedemenergy[n], 
      ma->tau_ak5pfjet_chargedmulti[n], ma->tau_ak5pfjet_neutralmulti[n], -1., -1., -1, 0, -1, -1, -1., -1., false, false, false, 0.0f, false, false, false, 0.0f),*/
  taudiscriminators(ma->runlist.find(ma->Run())->second.GetTauDiscriminators()){
  
  // Set tau four vector
  jetp4.SetPxPyPzE(ma->tau_ak5pfjet_px[n], ma->tau_ak5pfjet_py[n], ma->tau_ak5pfjet_pz[n], ma->tau_ak5pfjet_e[n]);

  unsigned int i;
  for(i = 0; i < ma->ak5pfjet_count; ++i)
  {
    // Note the pt/eta/phi between jets and taus are not exactly the same... maybe because
    // jet energy corrections have been applied to jets, and not to taus?
    // TODO: That's what we have the tau_ak5pfjet_e/px/py/pz for, no!? This should really match, actually...
    Jet j = ma->AK5PFJets(i);
    if(ROOT::Math::VectorUtil::DeltaR(j, *this) < 0.3)
    {
      jet = j;
      break;
    }
  }  // for(i = 0; i < ma->ak5pfjet_count; ++i)

  if(i == ma->ak5pfjet_count)
    {
      jet = Jet(0.0, 0.0, 0.0, 0.0, -1.0f, -1.0f, -1.0f, -1.0f, -1, -1, -1.0f, NULL, -1, -1, -1.0f, -1.0f, false, false, false, 0.0f, false, false, false, 0.0f, 0);
    }
  
  UInt_t begin = ma->tau_chargedbegin[n];
  UInt_t end = ma->tau_charged_count; 
  if(n < ma->tau_count - 1){
    end = ma->tau_chargedbegin[n+1];
  }
  for(UInt_t i = begin ; i < end ; i++){
    tracks.push_back( Track(Sqrt(ma->tau_charged_px[i]*ma->tau_charged_px[i]+
				 ma->tau_charged_py[i]*ma->tau_charged_py[i]+
				 ma->tau_charged_pz[i]*ma->tau_charged_pz[i]), 
			    ma->tau_charged_px[i], 
			    ma->tau_charged_py[i], 
			    ma->tau_charged_pz[i], 
			    ma->tau_charged_outerx[i], 
			    ma->tau_charged_outery[i], 
			    ma->tau_charged_outerz[i], 
			    ma->tau_charged_closestpointx[i], 
			    ma->tau_charged_closestpointy[i], 
			    ma->tau_charged_closestpointz[i], 
			    ma->tau_charged_chi2[i], 
			    ma->tau_charged_ndof[i], 
			    ma->tau_charged_dxy[i], 
			    ma->tau_charged_dxyerr[i], 
			    ma->tau_charged_dz[i], 
			    ma->tau_charged_dzerr[i], 
			    ma->tau_charged_charge[i], 
			    ma->tau_charged_nhits[i], 
			    ma->tau_charged_nmissinghits[i], 
			    ma->tau_charged_npixelhits[i], 
			    ma->tau_charged_npixellayers[i], 
			    ma->tau_charged_nstriplayers[i], 
			    ma->tau_charged_dedxharmonic2[i]));
  }  
}


// diTau
diTau::diTau(const Analyse* ma, UInt_t n): 
  ditau_Index_(ma->ditau_Index),
  ditau_leg1_index_(ma->ditau_leg1_index[n]),
  ditau_leg2_index_(ma->ditau_leg2_index[n]),
  ditau_Vtx_(TVector3(ma->ditau_VtxX, ma->ditau_VtxY, ma->ditau_VtxZ)),
  ditau_VtxErr_(TVector3(ma->ditau_VtxXErr, ma->ditau_VtxYErr, ma->ditau_VtxZErr)),
  ditau_reFitVtx_(TVector3( ma->ditau_reFitVtxX[n], ma->ditau_reFitVtxY[n], ma->ditau_reFitVtxZ[n])),
  ditau_reFitVtxErr_(TVector3( ma->ditau_reFitVtxXErr[n], ma->ditau_reFitVtxYErr[n], ma->ditau_reFitVtxZErr[n])),
  ditau_leg1_OPV_(TVector3( ma->ditau_leg1_X_OPV[n], ma->ditau_leg1_Y_OPV[n], ma->ditau_leg1_Z_OPV[n])),
  ditau_leg1_(TVector3( ma->ditau_leg1_X[n], ma->ditau_leg1_Y[n], ma->ditau_leg1_Z[n])),
  ditau_leg2_(TVector3( ma->ditau_leg2_X[n], ma->ditau_leg2_Y[n], ma->ditau_leg2_Z[n])),
  ditau_leg2_OPV_(TVector3( ma->ditau_leg2_X_OPV[n], ma->ditau_leg2_Y_OPV[n], ma->ditau_leg2_Z_OPV[n])),
  ditau_reFitVtxRho_(ma->ditau_reFitVtxRho[n]),
  ditau_reFitVtxNdof_(ma->ditau_reFitVtxNdof[n]),
  ditau_leg1_dxy_(ma->ditau_leg1_dxy[n]),
  ditau_leg1_dz_(ma->ditau_leg1_dz[n]),
  ditau_leg1_dxyErr_(ma->ditau_leg1_dxyErr[n]),
  ditau_leg1_dzErr_(ma->ditau_leg1_dzErr[n]),
  ditau_leg1_dxy_OPV_(ma->ditau_leg1_dxy_OPV[n]),
  ditau_leg1_dz_OPV_(ma->ditau_leg1_dz_OPV[n]),
  ditau_leg1_dxyErr_OPV_(ma->ditau_leg1_dxyErr_OPV[n]),
  ditau_leg1_dzErr_OPV_(ma->ditau_leg1_dzErr_OPV[n]),
  ditau_leg2_dxy_(ma->ditau_leg2_dxy[n]),
  ditau_leg2_dz_(ma->ditau_leg2_dz[n]),
  ditau_leg2_dxyErr_(ma->ditau_leg2_dxyErr[n]),
  ditau_leg2_dzErr_(ma->ditau_leg2_dzErr[n]),
  ditau_leg2_dxy_OPV_(ma->ditau_leg2_dxy_OPV[n]),
  ditau_leg2_dz_OPV_(ma->ditau_leg2_dz_OPV[n]),
  ditau_leg2_dxyErr_OPV_(ma->ditau_leg2_dxyErr_OPV[n]),
  ditau_leg2_dzErr_OPV_(ma->ditau_leg2_dzErr_OPV[n])
{
}


MuTauTauPair::MuTauTauPair( const Analyse* ma, UInt_t n) :
  TLorentzVector( ma->mutautaupair_leg1_px[n]+ma->mutautaupair_leg2_px[n], 
		  ma->mutautaupair_leg1_py[n]+ma->mutautaupair_leg2_py[n], 
		  ma->mutautaupair_leg1_pz[n]+ma->mutautaupair_leg2_pz[n],
		  ma->mutautaupair_leg1_energy[n]+ma->mutautaupair_leg2_energy[n]),
  leg1_px(ma->mutautaupair_leg1_px[n]),
  leg1_py(ma->mutautaupair_leg1_py[n]),
  leg1_pz(ma->mutautaupair_leg1_pz[n]),
  leg1_energy(ma->mutautaupair_leg1_energy[n]),
  leg2_px(ma->mutautaupair_leg2_px[n]),
  leg2_py(ma->mutautaupair_leg2_py[n]),
  leg2_pz(ma->mutautaupair_leg2_pz[n]),
  leg2_energy(ma->mutautaupair_leg2_energy[n]),
  mu_px(ma->mutautaupair_mu_px[n]),
  mu_py(ma->mutautaupair_mu_py[n]),
  mu_pz(ma->mutautaupair_mu_pz[n]),
  mu_energy(ma->mutautaupair_mu_energy[n]),
  svfit_int_valid(ma->mutautaupair_svfit_int_valid[n]),
  svfit_mass_int(ma->mutautaupair_svfit_mass_int[n]),
  svfit_mass_int_err_up(ma->mutautaupair_svfit_mass_int_err_up[n]),
  svfit_mass_int_err_down(ma->mutautaupair_svfit_mass_int_err_down[n])
{

}


Int_t Tau::TauDiscriminator(string disname) const{
  Int_t pos = -1;
  for(UInt_t i = 0 ; i < taudiscriminators->size() ; i++){
    if(disname == taudiscriminators->at(i)){
      pos = i;
      break;
    }
  }
  if(pos == -1) return(-1);
  if((dishps & ((ULong64_t)1<<(ULong64_t)pos)) != 0) return(1);
  return(0);
}


ElTauTauPair::ElTauTauPair( const Analyse* ma, UInt_t n) :
  TLorentzVector( ma->eltautaupair_leg1_px[n]+ma->eltautaupair_leg2_px[n], 
		  ma->eltautaupair_leg1_py[n]+ma->eltautaupair_leg2_py[n], 
		  ma->eltautaupair_leg1_pz[n]+ma->eltautaupair_leg2_pz[n],
		  ma->eltautaupair_leg1_energy[n]+ma->eltautaupair_leg2_energy[n]),
  leg1_px(ma->eltautaupair_leg1_px[n]),
  leg1_py(ma->eltautaupair_leg1_py[n]),
  leg1_pz(ma->eltautaupair_leg1_pz[n]),
  leg1_energy(ma->eltautaupair_leg1_energy[n]),
  leg2_px(ma->eltautaupair_leg2_px[n]),
  leg2_py(ma->eltautaupair_leg2_py[n]),
  leg2_pz(ma->eltautaupair_leg2_pz[n]),
  leg2_energy(ma->eltautaupair_leg2_energy[n]),
  el_px(ma->eltautaupair_el_px[n]),
  el_py(ma->eltautaupair_el_py[n]),
  el_pz(ma->eltautaupair_el_pz[n]),
  el_energy(ma->eltautaupair_el_energy[n]),

  svfit_int_valid(ma->eltautaupair_svfit_int_valid[n]),
  svfit_mass_int(ma->eltautaupair_svfit_mass_int[n]),
  svfit_mass_int_err_up(ma->eltautaupair_svfit_mass_int_err_up[n]),
  svfit_mass_int_err_down(ma->eltautaupair_svfit_mass_int_err_down[n])
{

}

Vertex::Vertex(Double_t X, Double_t Y, Double_t Z, Double_t Chi2, Double_t Ndof, Double_t Ptq, Int_t Ntracks, const Float_t* Cov):
  TVector3(X, Y, Z),
  chi2(Chi2),
  ndof(Ndof),
  ptq(Ptq),
  ntracks(Ntracks),
  covmatrix(3)
{
  covmatrix(0,0) = Cov[0];
  covmatrix(0,1) = Cov[1];
  covmatrix(0,2) = Cov[2];
  covmatrix(1,0) = Cov[1];
  covmatrix(1,1) = Cov[3];
  covmatrix(1,2) = Cov[4];
  covmatrix(2,0) = Cov[2];
  covmatrix(2,1) = Cov[4];
  covmatrix(2,2) = Cov[5];
}

BeamSpot::BeamSpot(Double_t X, Double_t Y, Double_t Z, const Float_t* Cov, Double_t Xwidth, Double_t Ywidth, Double_t Zsigma) :
  TVector3(X, Y, Z),
  covmatrix(3),
  xwidth(Xwidth),
  ywidth(Ywidth),
  zsigma(Zsigma)
{
	covmatrix(0,0) = Cov[0];
	covmatrix(0,1) = Cov[1];
	covmatrix(0,2) = Cov[2];
	covmatrix(1,0) = Cov[1];
	covmatrix(1,1) = Cov[3];
	covmatrix(1,2) = Cov[4];
	covmatrix(2,0) = Cov[2];
	covmatrix(2,1) = Cov[4];
	covmatrix(2,2) = Cov[5];
}

void splitstring(string input, vector<string>& output)
{
	UInt_t posstart = 0;
	for(UInt_t i = 0 ; i < input.size() ; i++)
	{
		if(input[i] == ' ')
		{
			output.push_back(input.substr(posstart, i-posstart));
			posstart=i+1;
		}
	}
}

string combinestring(const vector<string>& input)
{
	string output;
	for(UInt_t i = 0 ; i < input.size() ; i++)
	{
		output += input[i] + string(" ");
	}
	return(output);
}

//RunInfo
RunInfo::RunInfo(UInt_t Number, UInt_t Hltcount, string Hltnames, string Hltmunames, string Hltelnames, string Hltphotonnames, string Hlttaunames, string Hltjetnames, string Taudiscriminators, UInt_t Hltprescaletablescount, UInt_t* Hltprescaletables, UInt_t L1algocount, UInt_t L1algoprescaletablescount, UInt_t* L1algoprescaletables, UInt_t L1techcount, UInt_t L1techprescaletablescount, UInt_t* L1techprescaletables) : 
runnumber(Number),
hltnum(Hltcount),
hlttablesnum(Hltprescaletablescount/Hltcount),
l1algonum(L1algocount),
l1algotablesnum(L1algoprescaletablescount/L1algocount),
l1technum(L1techcount),
l1techtablesnum(L1techprescaletablescount/L1techcount)
{
	splitstring(Hltnames, hltnames);
	splitstring(Hltmunames, hltmunames);
	splitstring(Hltelnames, hltelnames);
	splitstring(Hlttaunames, hlttaunames);
	splitstring(Hltphotonnames, hltphotonnames);
	splitstring(Hltjetnames, hltjetnames);
	splitstring(Taudiscriminators, taudiscriminators);
	hltprescales.resize(NumHLT());
	for(UInt_t i = 0 ; i < NumHLT() ; i++)
	{
		for(UInt_t j = 0 ; j < NumHLTTables() ; j++)
		{
			hltprescales[i].push_back(Hltprescaletables[i+Hltcount*j]);
		}
	}

	l1algoprescales.resize(NumL1Algo());
	for(UInt_t i = 0 ; i < NumL1Algo() ; i++)
	{
		for(UInt_t j = 0 ; j < NumL1AlgoTables() ; j++)
		{
			l1algoprescales[i].push_back(L1algoprescaletables[i+L1algocount*j]);
		}
	}

	l1techprescales.resize(NumL1Tech());
	for(UInt_t i = 0 ; i < NumL1Tech() ; i++)
	{
		for(UInt_t j = 0 ; j < NumL1TechTables() ; j++)
		{
			l1techprescales[i].push_back(L1techprescaletables[i+L1techcount*j]);
		}
	}

}

vector<string> RunInfo::MatchTriggerNames(string name)
{
	boost::cmatch what;
	vector<string> result;

	boost::regex trigregex(name.c_str());
	for(int i = 0 ; i < NumHLT() ; i++)
	{
		if(boost::regex_match(HLTName(i).c_str(), what, trigregex))
		{
			result.push_back(HLTName(i));
		}
	}
	return(result);
}

Int_t RunInfo::HLTIndex(string name) const
{
	for(Int_t i = 0 ; i < Int_t(hltnames.size()) ; i++)
	{
		if(hltnames[i] == name) return(i);
	}
	return(-1);
}

string RunInfo::HLTAllNames() const
{
	return(combinestring(hltnames));
}

string RunInfo::HLTMuonAllNames() const
{
	return(combinestring(hltmunames));
}

string RunInfo::HLTElectronAllNames() const
{
	return(combinestring(hltelnames));
}

string RunInfo::HLTTauAllNames() const
{
	return(combinestring(hlttaunames));
}

string RunInfo::HLTPhotonAllNames() const
{
	return(combinestring(hltphotonnames));
}

string RunInfo::HLTJetAllNames() const
{
	return(combinestring(hltjetnames));
}

string RunInfo::TauDiscriminatorsAllNames() const
{
  return(combinestring(taudiscriminators));
}

void TriggerRun::SetBlock(UInt_t block, Int_t index, UInt_t prescale)
{
	if(lumiinfo.find(block) == lumiinfo.end())
	{
		lumiinfo[block] = TriggerLumi(index, prescale);
	}
}

TriggerLumi TriggerRun::GetBlock(UInt_t block)
{
	//return(lumiinfo[block]);
	map<UInt_t, TriggerLumi>::const_iterator blockinfo = lumiinfo.find(block);
	if(blockinfo == lumiinfo.end())
	{
		return(TriggerLumi(-1, 0));
	}
	return(blockinfo->second);
}

TriggerSelection::TriggerSelection(Analyse* an, vector<string> names, bool useprescaled) : AN(an)
{
	for(map<UInt_t, RunInfo>::iterator a = AN->runlist.begin() ; a != AN->runlist.end(); ++a)
	{
		UInt_t runnumber = a->first;
		Float_t runlumi = 0.;
		vector<UInt_t> indices;
		for(UInt_t i = 0 ; i < names.size() ; i++)
		{
			vector<string> goodnames = a->second.MatchTriggerNames(names[i]);
			if(goodnames.size() == 1)
			{
				indices.push_back(a->second.HLTIndex(goodnames[0]));
			}
			else if(goodnames.size() > 1)
			{
				cout << "WARNING: TriggerSelection::TriggerSelection: " << names[i] << " is not a unique selection!" << endl;
			}
		}

		map<UInt_t, map<UInt_t, Luminosity> >::iterator blocklist = AN->lumilist.find(runnumber);
		for(map<UInt_t, Luminosity>::iterator b = blocklist->second.begin() ; b != blocklist->second.end() ; ++b)
		{
			Int_t minindex = -1;
			UInt_t minprescale = 10000000;

			for(UInt_t i = 0 ; i < indices.size() ; i++)
			{
				UInt_t prescale = a->second.HLTPrescale(indices[i], b->second.HLTTable());
				if(prescale <= 0)
				{
					continue;
				}
				if(prescale < minprescale)
				{
					minindex = indices[i];
					minprescale = prescale;
				}
				if(prescale == 1)
				{
					break;
				}
			}
			if(useprescaled == false && minprescale != 1)
			{
				runinfo[runnumber].SetBlock(b->first, -1, minprescale);
				continue;
			}
			runlumi += b->second.LumiValue() * minprescale;
			runinfo[runnumber].SetBlock(b->first, minindex, minprescale);
		}
		runinfo[runnumber].Lumi(runlumi);
	}
}

Int_t TriggerSelection::Result()
{
	const TriggerLumi& triggerlumi = runinfo[AN->Run()].GetBlock(AN->LumiBlock());
	if(triggerlumi.Index() == -1){return(0);}
	if(AN->GetHLTrigger(triggerlumi.Index()))
	{
		return(triggerlumi.Prescale());
	}
	else
	{
		return(triggerlumi.Prescale() * -1);
	}
}

Float_t TriggerSelection::LumiUsed(Int_t format)
{
	Double_t lumi = 0.;
	Double_t alllumi = 0.;
	Double_t zerolumi = 0.;
	UInt_t nolumiinfo = 0;
	Double_t events = 0.;
	Double_t eventsprocessed = 0.;
	for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = AN->lumilist.begin() ; a != AN->lumilist.end() ; ++a)
	{
		Double_t runlumi = 0.;
		Double_t runalllumi = 0.;
		Double_t runzerolumi = 0.;
		Double_t runevents = 0.;
		Double_t runeventsprocessed = 0.;
		UInt_t runnolumiinfo = 0;
		UInt_t numblocks = 0;	
		for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++)
		{
			const TriggerLumi& triggerlumi = runinfo[a->first].GetBlock(b->first);
			if(b->second.LumiValue() == -1 || triggerlumi.Index() == -1){continue;}
			if(AN->IsInRange(a->first, b->first))
			{
				numblocks++;
			//	cout << a->first << " " << b->first << " " << minLumi << " " << maxLumi << endl;
				if(b->second == true)
				{
					if(b->second.NumEvents() > 0)
					{
						runlumi += b->second.ProcessedLumi()*triggerlumi.Prescale();
						runalllumi += b->second.LumiValue()*triggerlumi.Prescale();
						runevents += b->second.NumEvents();
						runeventsprocessed += b->second.NumEventsProcessed();
						if(b->second.ProcessedFraction() > 1)
						{
							cerr << "WARNING GetLumi: You ran on " << b->second.ProcessedFraction()*100 <<"\% of available events in Run " << a->first << ", Lumiblock " << b->first <<". :-O" << endl;  
						}
						if(format >= 3) cout << "    Block: " << b->first << ", fraction: " << b->second.ProcessedFraction() << ", lumi: " << b->second.ProcessedLumi() << " pb^-1" << endl;
					}
					else
					{
						runzerolumi += b->second.LumiValue()*triggerlumi.Prescale(); 
						runalllumi += b->second.LumiValue()*triggerlumi.Prescale();
					}
				}
				else if(b->second == false)
				{
					runnolumiinfo += b->second.NumEventsProcessed();
				}
			}
		}
		if(a->first <= AN->maxRun && a->first >= AN->minRun)
		{
			if(runevents != 0) zerolumi *= runeventsprocessed/runevents;
			if(format >= 2)  cout << "Run: " << a->first << ", blocks: " << numblocks << ", fraction: " << runeventsprocessed/runevents << ", total lumi: " << runalllumi << " pb^-1, processed lumi: " << runlumi+zerolumi << " zero event lumi: " << zerolumi << " pb^-1, no lumi info for: " << runnolumiinfo << " event(s)."<< endl;
			lumi += runlumi;
			zerolumi += runzerolumi;
			alllumi += runalllumi;
			events += runevents;
			eventsprocessed += runeventsprocessed;
			nolumiinfo += runnolumiinfo; 
		}
	}
	if(format >= 1)
	{
		cout << "From run " << AN->minRun << ", block " << AN->minLumi << " to run " << AN->maxRun << ", block " << AN->maxLumi << "." << endl;
		cout << "Total lumi: " << lumi+zerolumi << " pb^-1, zero event lumi: " << zerolumi << " pb^-1, (lumi in Range: " << alllumi << " pb^-1), no lumi-info for " << nolumiinfo << " event(s), events: " << eventsprocessed << "(" << events << ")"  << endl;
	}
return(lumi+zerolumi);

}

Float_t TriggerSelection::LumiBeforeEvent()
{
	Float_t lumi;
	for(map<UInt_t, TriggerRun>::iterator a = runinfo.begin() ; a != runinfo.end(); ++a)
	{
		if(a->first < AN->Run()) 
		{
			lumi+=a->second.Lumi();
		}
		else
		{
			break;
		}
	}

	map<UInt_t, map<UInt_t, Luminosity> >::iterator blocklist = AN->lumilist.find(AN->Run());
	for(map<UInt_t, Luminosity>::iterator b = blocklist->second.begin() ; b != blocklist->second.end() ; ++b)
	{
		if(b->first <= AN->LumiBlock())
		{
			const TriggerLumi& triggerlumi = runinfo[AN->Run()].GetBlock(b->first);
			if(triggerlumi.Index() != -1)
			{
				lumi+=b->second.LumiValue()*triggerlumi.Prescale();
			}
		}
		else
		{
			break;
		}
	}
	return(lumi);
}


string TriggerSelection::GetTriggerName(UInt_t run, UInt_t lumiblock)
{
	if(run == 0 && lumiblock == 0)
	{
		run = AN->Run();
		lumiblock = AN->LumiBlock();
	}

	Int_t index = runinfo[run].GetBlock(lumiblock).Index();
	if(index != -1)
	{
		return(AN->runlist[run].HLTName(index));
	}
	cerr << "TriggerSelection::GetTriggerName: Invalid run and/or lumiblock." << endl;
	return("ERROR: TriggerSelection::GetTriggerName: Invalid run and/or lumiblock.");
}

void TriggerSelection::PrintInfo()
{
	Float_t lumi = 0.;
	Float_t lumiuse = 0.;
	Float_t curprescale;
	string curname;
	bool beg = true;
	UInt_t begrun, begblock;
	Float_t beglumi;
	UInt_t run, block;
	UInt_t prevrun, prevblock;

	for(map<UInt_t, TriggerRun>::iterator runit = runinfo.begin() ; runit != runinfo.end() ; ++runit)
	{
		for(map< UInt_t, TriggerLumi >::iterator lumiit = runit->second.Begin() ; lumiit != runit->second.End() ; ++lumiit)
		{
			run = runit->first;
			block = lumiit->first;
			Int_t index = lumiit->second.Index();
			Float_t prescale = 0;
			string name("NA");
			if(index != -1)
			{
				name = AN->runlist[run].HLTName(index);
				prescale = lumiit->second.Prescale();
				lumiuse += AN->lumilist[run][block].LumiValue();
			}
			if(beg)
			{
				begrun = run;
				begblock = block;
				beglumi = lumi;
				curname = name;
				curprescale = prescale;
				beg = false;
			}
			lumi += AN->lumilist[run][block].LumiValue();
			if(curname != name || curprescale != prescale)
			{
				cout << begrun << ":" << begblock << " - " << prevrun << ":" << prevblock << ", " << curname << "(" << curprescale << ")" << endl;
				begrun = run;
				begblock = block;
				beglumi = lumi;
				curname = name;
				curprescale = prescale;
			}

			prevrun = run;
			prevblock = block;
		}
	}

	cout << begrun << ":" << begblock << " - " << prevrun << ":" << prevblock << ", " << curname << "(" << curprescale << ")" << endl;
}




