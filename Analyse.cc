#include "Analyse.h"
#include "DataForm.h"
#include <TH1D.h>

Analyse::Analyse() : 
  currentfile(0),
  currentmin(0),
  currentmax(0),
  currentloadtype(-1),
  loadbeamspot(0),
  loadsecvertices(0),
  loadmuons(0),
  loadelectrons(0),
  loadphotons(0),
  loadtaus(0),
  loadtautrackiso(0),
  loadmet(0),
  loadak5calojets(0),
  loadak5jptjets(0),
  loadak5pfjets(0),
  loadtracks(0),
  loadsuperclusters(0),
  loadprimvertices(0),
  loadtrigger(0),
  loadgeninfo(0),
  loadgenparticles(0),
  loadallgenparticles(0),
  printinfo(1),
  processed(0),
  duplicatecheck(false),
  skimtree(0),
  skimfilename("None"),
  skimfile(0),
  lumicalculation(false),
  minRun(1000000000),
  minLumi(100000000),
  maxRun(0),
  maxLumi(0),
  jsonfilter(false),
  usepileupinfo(false),
  pileUpDistMinus(200),
  pileUpDist(200),
  pileUpDistPlus(200),
  useprimvertexinfo(false),
  primVertexDist(200),
  track_count(0),
  primvertex_count(0),
  supercluster_count(0),
  supercluster_basiccluster_count(0),
  supercluster_basiccluster_hit_count(0),
  supercluster_escluster_count(0),
  supercluster_escluster_hit_count(0),
  muon_count(0),
  ak5calojet_count(0),
  ak5jptjet_count(0),
  ak5pfjet_count(0),
  electron_count(0),
  photon_count(0),
  conversion_count(0),
  tau_count(0),
  tau_charged_count(0),
  secvertices_count(0),
  genparticles_count(0),
  genallparticles_count(0),
  genallparticlesmother_count(0),
  genallparticlesdaughter_count(0)
{
	for(int i = 0; i < M_taumaxcount; ++i)
		tau_e[i] = -1.0f;
}

Analyse::~Analyse()
{
	if(skimfilename != string("None"))
	{
		TFile* file = skimtree->GetCurrentFile();
		file->Write();
		file->Close();
		string dirname;
		string outfilename;
		string lumifilename;

		UInt_t slashpos = skimfilename.find_last_of("/");
		if(slashpos != skimfilename.size())
		  {
		    dirname = skimfilename.substr(0, slashpos);
		    outfilename = skimfilename.substr(slashpos+1);
		  }
		else
		{
		  outfilename = skimfilename;
		}
		
		lumifilename = string("LUMI_") + outfilename;
		cout << "Writing LUMI-file: " << lumifilename << endl;
		
		WriteLumiFile(dirname + string("/") + lumifilename);
	}
}

void Analyse::SetLoad()
{
  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("errors", true);
  tree->SetBranchStatus("event_*", true);
  tree->SetBranchStatus("numpileupinteractions*", true);
  tree->SetBranchStatus("numtruepileupinteractions", true);
  tree->SetBranchStatus("ak5pfjet_rho", true);
  tree->SetBranchStatus("ak5pfjet_sigma", true);
  
  if(loadbeamspot == 1)         tree->SetBranchStatus("beamspot_*", true);
  else 		                tree->SetBranchStatus("beamspot_*", false);
  
  if(loadsecvertices == 1)	tree->SetBranchStatus("secvertices_*", true);
  else 		                tree->SetBranchStatus("secvertices_*", false);
  
  if(loadmuons == 1) 		tree->SetBranchStatus("muon_*", true);
  else		                tree->SetBranchStatus("muon_*", false);
  
  if(loadelectrons == 1) 	tree->SetBranchStatus("electron_*", true);
  else 		                tree->SetBranchStatus("electron_*", false);
  
  if(loadphotons == 1)		tree->SetBranchStatus("photon_*", true);
  else 		                tree->SetBranchStatus("photon_*", false);

  if(loadtaus == 1){
    tree->SetBranchStatus("tau_*", true);
    tree->SetBranchStatus("ditau*", true);
    tree->SetBranchStatus("mutautaupair_*", false);
    tree->SetBranchStatus("eltautaupair_*", false);
  } else {
    tree->SetBranchStatus("tau_*", false);
    tree->SetBranchStatus("mutautaupair_*", false);
    tree->SetBranchStatus("eltautaupair_*", false);
  }
  
  if(loadmet == 1)
    {
      tree->SetBranchStatus("calomet_*"   , true);
      tree->SetBranchStatus("tcmet_*"     , true);
      tree->SetBranchStatus("pfmet_*"     , true);
      tree->SetBranchStatus("pfmettype1_*", true);
      tree->SetBranchStatus("pfmetmva_*"  , true);
    }
  else
    {
      tree->SetBranchStatus("calomet_*"   , false);
      tree->SetBranchStatus("tcmet_*"     , false);
      tree->SetBranchStatus("pfmet_*"     , false);
      tree->SetBranchStatus("pfmettype1_*", false);
      tree->SetBranchStatus("pfmetmvaD_*" , false);
    }
  
  if(loadak5calojets == 1) 		tree->SetBranchStatus("ak5calojet_*", true);
  else        		                tree->SetBranchStatus("ak5calojet_*", false);
  
  if(loadak5jptjets == 1) 		tree->SetBranchStatus("ak5jptjet_*", true);
  else 		                        tree->SetBranchStatus("ak5jptjet_*", false);
  
  if(loadak5pfjets == 1)  		tree->SetBranchStatus("ak5pfjet_*", true);
  else
    {
      tree->SetBranchStatus("ak5pfjet_*", false);
      tree->SetBranchStatus("ak5pfjet_rho", true);
      tree->SetBranchStatus("ak5pfjet_sigma", true);
    }
  
  if(loadtracks == 1) 		tree->SetBranchStatus("track_*", true);
  else     		        tree->SetBranchStatus("track_*", false);
  
  if(loadsuperclusters > 0)
    {
      tree->SetBranchStatus("supercluster*", true);
      if(loadsuperclusters != 3 && loadsuperclusters != 4)
	{
	  tree->SetBranchStatus("supercluster_basiccluster_hit*", false);
	  tree->SetBranchStatus("supercluster_escluster_hit*", false);
	}
      if(loadsuperclusters != 2 && loadsuperclusters != 4)
	{
	  tree->SetBranchStatus("supercluster_basiccluster*", false);
	  tree->SetBranchStatus("supercluster_escluster*", false);
	}
    }
  else
    {
      tree->SetBranchStatus("supercluster*", false);
    }
  
  if(loadprimvertices == 1)     tree->SetBranchStatus("primvertex_*", true);
  else 		tree->SetBranchStatus("primvertex_*", false);
  
  if(loadtrigger == 1) 		tree->SetBranchStatus("trigger_*", true);
  else 		tree->SetBranchStatus("trigger_*", false);
  
  if(loadgeninfo == 1)
    {
      tree->SetBranchStatus("genweight", true);
      tree->SetBranchStatus("genid1", true);
      tree->SetBranchStatus("genx1", true);
      tree->SetBranchStatus("genid2", true);
      tree->SetBranchStatus("genx2", true);
      tree->SetBranchStatus("genScale", true);
      tree->SetBranchStatus("genmet*", true);
    }
  else
    {
      tree->SetBranchStatus("genweight", false);
      tree->SetBranchStatus("genid1", false);
      tree->SetBranchStatus("genx1", false);
      tree->SetBranchStatus("genid2", false);
      tree->SetBranchStatus("genx2", false);
      tree->SetBranchStatus("genScale", false);
      tree->SetBranchStatus("genmet*", false);
    }
  
  if(loadallgenparticles == 1) 		tree->SetBranchStatus("genall*", true);
  else 		tree->SetBranchStatus("genall*", false);
  
  if(loadgenparticles == 1) 		tree->SetBranchStatus("genparticle*", true);
  else 		tree->SetBranchStatus("genparticle*", false);
}

bool sortac1b(string a, string b)
{
  vector<string> strs;
  boost::split(strs, a, boost::is_any_of("_"));
  int na = atoi(strs[1].c_str());
  boost::split(strs, b, boost::is_any_of("_"));
  int nb = atoi(strs[1].c_str());
  return(na < nb);
}

bool sortroot(string a, string b)
{
	vector<string> strs;
	boost::split(strs, a, boost::is_any_of("_."));
	int na = atoi(strs[strs.size()-2].c_str());
	boost::split(strs, b, boost::is_any_of("_."));
	int nb = atoi(strs[strs.size()-2].c_str());
	return(na < nb);
}

Int_t Analyse::AddDir(string path)
{
	DIR* dir = opendir(path.c_str());
	if(dir == 0)
	{cerr << "ERROR AddDir: could not open directory!" << endl; exit(0);}
	dirent* currfile;
	vector<string> filenames;
	int orig = 0;
	while((currfile = readdir(dir)) != 0)
	{
	  string currname(currfile->d_name);
	  if(currname.find(".root") != string::npos)
	    {
	      if(currname.find("LUMI") != string::npos)
		{
		  AddLumiFile(path+string("/")+currname);
		  cout << "Added LUMI-file: " << currname << endl;
		  continue;
		}
	      if(currname.find("AC1B") != string::npos){orig++;}
	      filenames.push_back(currname);
	    }
	}
	
	if(filenames.size() == orig)
	  {
	    sort(filenames.begin(), filenames.end(), sortac1b);
	  }
	else
	  {
	    sort(filenames.begin(), filenames.end(), sortroot);
	    
	  }
	
	for(int i = 0 ; i < filenames.size() ; i++)
	  {
	    cout << i+1 << ". Add: " << filenames[i] << endl;
	    AddFile(path+string("/")+filenames[i]);
	  }
	
	closedir(dir);
}

Int_t Analyse::AddFile(string filename)
{
  TTree* test = 0;
  TFile* file = new TFile(filename.c_str(), "READ");
  if(file->IsZombie())
    {
      cerr << "ERROR AddFile: " << filename << " is not a valid file." << endl;
    }
  file->GetObject("makeroottree/AC1B", test);
  
  if(test != 0)
    {
      if(treename.size() == 0)
	{
	  treename = string("makeroottree/AC1B");
	}
      filenames.push_back(filename);
      
      if(filelimits.size() == 0)
	{
	  filelimits.push_back(test->GetEntries());
	}
      else
	{
	  filelimits.push_back(test->GetEntries() + filelimits[filelimits.size()-1]);
	}
      file->Close();
      return(1);
    }
  
  file->GetObject("AC1B", test);
  if(test != 0)
    {
      if(treename.size() == 0)
	{
	  treename = string("AC1B");
	}
      filenames.push_back(filename);
      if(filelimits.size() == 0)
	{
	  filelimits.push_back(test->GetEntries());
	}
      else
	{
	  filelimits.push_back(test->GetEntries() + filelimits[filelimits.size()-1]);
	}
      file->Close();
      return(1);
    }
  
  cerr <<"ERROR Addfile: " << filename << " no AC1B tree found." << endl;
  return(-1);
  
}

void Analyse::GetEvent(Long64_t num, Int_t loadtype)
{
  bool changed = false;
  
  if(!(num >= currentmin && num < currentmax))
    {
      if(currentfile != 0)  delete currentfile;
      
      for(UInt_t i = 0 ; i < filelimits.size() ; i++)
	{
	  if(num < filelimits[i])
	    {
	      currentfile = new TFile(filenames[i].c_str(), "READ");
	      currentfile->GetObject(treename.c_str(), tree);
	      currentmax = filelimits[i];
	      if(i == 0)  currentmin = 0;
	      else        currentmin = filelimits[i-1];
	      Load();
	      break;
	    }
	}
      changed = true;
    }//if(!(num >= currentmin && num < currentmax))
  
  
  if(changed || currentloadtype != loadtype)
    {
      if(loadtype == 0)  SetLoad();
      else if(loadtype == 1)
	{
	  tree->SetBranchStatus("*", false);
	  tree->SetBranchStatus("numpileupinteractions*", true);
	}
      else if(loadtype == 2)
	{
	  tree->SetBranchStatus("*", false);
	  tree->SetBranchStatus("primvertex_*", true);
	}
      else  cerr << "ERROR GetEvent: This is an unknown loadtype." << endl;
      currentloadtype = loadtype;
    }
  tree->GetEntry(num - currentmin);
}

void Analyse::Load(){
  tree->SetBranchAddress("errors", &errors);
  tree->SetBranchAddress("event_nr", &event_nr);
  tree->SetBranchAddress("event_luminosityblock", &event_luminosityblock);
  tree->SetBranchAddress("event_run", &event_run);
  tree->SetBranchAddress("event_timeunix", &event_timeunix);
  tree->SetBranchAddress("event_timemicrosec", &event_timemicrosec);
  tree->SetBranchAddress("trigger_level1bits", trigger_level1bits);
  tree->SetBranchAddress("trigger_level1", trigger_level1);
  tree->SetBranchAddress("trigger_HLT", trigger_HLT);
  
  tree->SetBranchAddress("beamspot_x", &beamspot_x);
  tree->SetBranchAddress("beamspot_y", &beamspot_y);
  tree->SetBranchAddress("beamspot_z", &beamspot_z);
  tree->SetBranchAddress("beamspot_xwidth", &beamspot_xwidth);
  tree->SetBranchAddress("beamspot_ywidth", &beamspot_ywidth);
  tree->SetBranchAddress("beamspot_zsigma", &beamspot_zsigma);
  tree->SetBranchAddress("beamspot_cov", &beamspot_cov);
  
  tree->SetBranchAddress("track_count", &track_count);
  tree->SetBranchAddress("track_px", track_px);
  tree->SetBranchAddress("track_py", track_py);
  tree->SetBranchAddress("track_pz", track_pz);
  tree->SetBranchAddress("track_outerx", track_outerx);
  tree->SetBranchAddress("track_outery", track_outery);
  tree->SetBranchAddress("track_outerz", track_outerz);
  tree->SetBranchAddress("track_closestpointx", track_closestpointx);
  tree->SetBranchAddress("track_closestpointy", track_closestpointy);
  tree->SetBranchAddress("track_closestpointz", track_closestpointz);
  tree->SetBranchAddress("track_chi2", track_chi2);
  tree->SetBranchAddress("track_ndof", track_ndof);
  tree->SetBranchAddress("track_dxy", track_dxy);
  tree->SetBranchAddress("track_dxyerr", track_dxyerr);
  tree->SetBranchAddress("track_dz", track_dz);
  tree->SetBranchAddress("track_dzerr", track_dzerr);
  tree->SetBranchAddress("track_dedxharmonic2", track_dedxharmonic2);
  tree->SetBranchAddress("track_charge", track_charge);
  tree->SetBranchAddress("track_nhits", track_nhits);
  tree->SetBranchAddress("track_npixelhits", track_npixelhits);
  tree->SetBranchAddress("track_nmissinghits", track_nmissinghits);
  tree->SetBranchAddress("track_npixellayers", track_npixellayers);
  tree->SetBranchAddress("track_nstriplayers", track_nstriplayers);
  
  tree->SetBranchAddress("primvertex_count", &primvertex_count);
  tree->SetBranchAddress("primvertex_x", &primvertex_x);
  tree->SetBranchAddress("primvertex_y", &primvertex_y);
  tree->SetBranchAddress("primvertex_z", &primvertex_z);
  tree->SetBranchAddress("primvertex_chi2", &primvertex_chi2);
  tree->SetBranchAddress("primvertex_ndof", &primvertex_ndof);
  tree->SetBranchAddress("primvertex_ptq", &primvertex_ptq);
  tree->SetBranchAddress("primvertex_ntracks", &primvertex_ntracks);
  tree->SetBranchAddress("primvertex_cov", primvertex_cov);
  
  tree->SetBranchAddress("supercluster_count", &supercluster_count);
  tree->SetBranchAddress("supercluster_e", supercluster_e);
  tree->SetBranchAddress("supercluster_x", supercluster_x);
  tree->SetBranchAddress("supercluster_y", supercluster_y);
  tree->SetBranchAddress("supercluster_z", supercluster_z);
  tree->SetBranchAddress("supercluster_rawe", supercluster_rawe);
  tree->SetBranchAddress("supercluster_phiwidth", supercluster_phiwidth);
  tree->SetBranchAddress("supercluster_etawidth", supercluster_etawidth);
  tree->SetBranchAddress("supercluster_nbasiccluster", supercluster_nbasiccluster);
  tree->SetBranchAddress("supercluster_basicclusterbegin", supercluster_basicclusterbegin);
  tree->SetBranchAddress("supercluster_esclusterbegin", supercluster_esclusterbegin);
  
  tree->SetBranchAddress("supercluster_basiccluster_count", &supercluster_basiccluster_count);
  tree->SetBranchAddress("supercluster_basiccluster_e", supercluster_basiccluster_e);
  tree->SetBranchAddress("supercluster_basiccluster_x", supercluster_basiccluster_x);
  tree->SetBranchAddress("supercluster_basiccluster_y", supercluster_basiccluster_y);
  tree->SetBranchAddress("supercluster_basiccluster_z", supercluster_basiccluster_z);
  tree->SetBranchAddress("supercluster_basiccluster_nhit", supercluster_basiccluster_nhit);
  tree->SetBranchAddress("supercluster_basiccluster_hitbegin", supercluster_basiccluster_hitbegin);
  
  tree->SetBranchAddress("supercluster_basiccluster_hit_count", &supercluster_basiccluster_hit_count);
  tree->SetBranchAddress("supercluster_basiccluster_hit_e", supercluster_basiccluster_hit_e);
  tree->SetBranchAddress("supercluster_basiccluster_hit_x", supercluster_basiccluster_hit_x);
  tree->SetBranchAddress("supercluster_basiccluster_hit_y", supercluster_basiccluster_hit_y);
  tree->SetBranchAddress("supercluster_basiccluster_hit_z", supercluster_basiccluster_hit_z);
  
  tree->SetBranchAddress("supercluster_escluster_count", &supercluster_escluster_count);
  tree->SetBranchAddress("supercluster_escluster_e", supercluster_escluster_e);
  tree->SetBranchAddress("supercluster_escluster_x", supercluster_escluster_x);
  tree->SetBranchAddress("supercluster_escluster_y", supercluster_escluster_y);
  tree->SetBranchAddress("supercluster_escluster_z", supercluster_escluster_z);
  tree->SetBranchAddress("supercluster_escluster_nhit", supercluster_escluster_nhit);
  tree->SetBranchAddress("supercluster_escluster_hitbegin", supercluster_escluster_hitbegin);
  
  tree->SetBranchAddress("supercluster_escluster_hit_count", &supercluster_escluster_hit_count);
  tree->SetBranchAddress("supercluster_escluster_hit_e", supercluster_escluster_hit_e);
  tree->SetBranchAddress("supercluster_escluster_hit_x", supercluster_escluster_hit_x);
  tree->SetBranchAddress("supercluster_escluster_hit_y", supercluster_escluster_hit_y);
  tree->SetBranchAddress("supercluster_escluster_hit_z", supercluster_escluster_hit_z);
  
  tree->SetBranchAddress("muon_count", &muon_count);
  tree->SetBranchAddress("muon_px", muon_px);
  tree->SetBranchAddress("muon_py", muon_py);
  tree->SetBranchAddress("muon_pz", muon_pz);
  tree->SetBranchAddress("muon_pterror", muon_pterror);
  tree->SetBranchAddress("muon_chi2", muon_chi2);
  tree->SetBranchAddress("muon_ndof", muon_ndof);
  tree->SetBranchAddress("muon_nchambers", muon_nchambers);
  tree->SetBranchAddress("muon_innertrack_px", muon_innertrack_px);
  tree->SetBranchAddress("muon_innertrack_py", muon_innertrack_py);
  tree->SetBranchAddress("muon_innertrack_pz", muon_innertrack_pz);
  tree->SetBranchAddress("muon_innertrack_outerx", muon_innertrack_outerx);
  tree->SetBranchAddress("muon_innertrack_outery", muon_innertrack_outery);
  tree->SetBranchAddress("muon_innertrack_outerz", muon_innertrack_outerz);
  tree->SetBranchAddress("muon_innertrack_closestpointx", muon_innertrack_closestpointx);
  tree->SetBranchAddress("muon_innertrack_closestpointy", muon_innertrack_closestpointy);
  tree->SetBranchAddress("muon_innertrack_closestpointz", muon_innertrack_closestpointz);
  tree->SetBranchAddress("muon_innertrack_chi2", muon_innertrack_chi2);
  tree->SetBranchAddress("muon_innertrack_ndof", muon_innertrack_ndof);
  tree->SetBranchAddress("muon_innertrack_dxy", muon_innertrack_dxy);
  tree->SetBranchAddress("muon_innertrack_dxyerr", muon_innertrack_dxyerr);
  tree->SetBranchAddress("muon_innertrack_dz", muon_innertrack_dz);
  tree->SetBranchAddress("muon_innertrack_dzerr", muon_innertrack_dzerr);
  tree->SetBranchAddress("muon_innertrack_dedxharmonic2", muon_innertrack_dedxharmonic2);
  tree->SetBranchAddress("muon_innertrack_charge", muon_innertrack_charge);
  tree->SetBranchAddress("muon_innertrack_nhits", muon_innertrack_nhits);
  tree->SetBranchAddress("muon_innertrack_npixelhits", muon_innertrack_npixelhits);
  tree->SetBranchAddress("muon_innertrack_nmissinghits", muon_innertrack_nmissinghits);
  tree->SetBranchAddress("muon_innertrack_npixellayers", muon_innertrack_npixellayers);
  tree->SetBranchAddress("muon_innertrack_nstriplayers", muon_innertrack_nstriplayers);
  tree->SetBranchAddress("muon_outertrack_px", muon_outertrack_px);
  tree->SetBranchAddress("muon_outertrack_py", muon_outertrack_py);
  tree->SetBranchAddress("muon_outertrack_pz", muon_outertrack_pz);
  tree->SetBranchAddress("muon_outertrack_hits", muon_outertrack_hits);
  tree->SetBranchAddress("muon_outertrack_missinghits", muon_outertrack_missinghits);
  tree->SetBranchAddress("muon_outertrack_chi2", muon_outertrack_chi2);
  tree->SetBranchAddress("muon_outertrack_ndof", muon_outertrack_ndof);
  tree->SetBranchAddress("muon_isolationr3track", muon_isolationr3track);
  tree->SetBranchAddress("muon_isolationr3ntrack", muon_isolationr3ntrack);
  tree->SetBranchAddress("muon_isolationr3ecal", muon_isolationr3ecal);
  tree->SetBranchAddress("muon_isolationr3hcal", muon_isolationr3hcal);
  tree->SetBranchAddress("muon_pfisolationr3_sumchargedhadronpt", muon_pfisolationr3_sumchargedhadronpt);
  tree->SetBranchAddress("muon_pfisolationr3_sumchargedparticlept", muon_pfisolationr3_sumchargedparticlept);
  tree->SetBranchAddress("muon_pfisolationr3_sumneutralhadronet", muon_pfisolationr3_sumneutralhadronet);
  tree->SetBranchAddress("muon_pfisolationr3_sumphotonet", muon_pfisolationr3_sumphotonet);
  tree->SetBranchAddress("muon_pfisolationr3_sumPUpt", muon_pfisolationr3_sumPUpt);
  tree->SetBranchAddress("muon_pfisolationr4_sumchargedhadronpt", muon_pfisolationr4_sumchargedhadronpt);
  tree->SetBranchAddress("muon_pfisolationr4_sumchargedparticlept", muon_pfisolationr4_sumchargedparticlept);
  tree->SetBranchAddress("muon_pfisolationr4_sumneutralhadronet", muon_pfisolationr4_sumneutralhadronet);
  tree->SetBranchAddress("muon_pfisolationr4_sumphotonet", muon_pfisolationr4_sumphotonet);
  tree->SetBranchAddress("muon_pfisolationr4_sumPUpt", muon_pfisolationr4_sumPUpt);
  tree->SetBranchAddress("muon_ecalenergy", muon_ecalenergy);
  tree->SetBranchAddress("muon_hcalenergy", muon_hcalenergy);
  tree->SetBranchAddress("muon_charge", muon_charge);
  tree->SetBranchAddress("muon_numchambers", muon_numchambers);
  tree->SetBranchAddress("muon_numchamberswithsegments", muon_numchamberswithsegments);
  tree->SetBranchAddress("muon_type", muon_type);
  tree->SetBranchAddress("muon_trigger", muon_trigger);
  tree->SetBranchAddress("muon_trackermuonquality", muon_trackermuonquality);
  tree->SetBranchAddress("muon_mva_id", muon_mva_id);
  tree->SetBranchAddress("muon_mva_iso", muon_mva_iso);
  
  tree->SetBranchAddress("ak5calojet_count", &ak5calojet_count);
  tree->SetBranchAddress("ak5calojet_e", ak5calojet_e);
  tree->SetBranchAddress("ak5calojet_px", ak5calojet_px);
  tree->SetBranchAddress("ak5calojet_py", ak5calojet_py);
  tree->SetBranchAddress("ak5calojet_pz", ak5calojet_pz);
  tree->SetBranchAddress("ak5calojet_hadronicenergy", ak5calojet_hadronicenergy);
  tree->SetBranchAddress("ak5calojet_emenergy", ak5calojet_emenergy);
  tree->SetBranchAddress("ak5calojet_energycorr", ak5calojet_energycorr);
  tree->SetBranchAddress("ak5calojet_fhpd", ak5calojet_fhpd);
  tree->SetBranchAddress("ak5calojet_restrictedemf", ak5calojet_restrictedemf);
  tree->SetBranchAddress("ak5calojet_btag", ak5calojet_btag);
  tree->SetBranchAddress("ak5calojet_n90", ak5calojet_n90);
  tree->SetBranchAddress("ak5calojet_n60", ak5calojet_n60);

  tree->SetBranchAddress("ak5jptjet_count", &ak5jptjet_count);
  tree->SetBranchAddress("ak5jptjet_e", ak5jptjet_e);
  tree->SetBranchAddress("ak5jptjet_px", ak5jptjet_px);
  tree->SetBranchAddress("ak5jptjet_py", ak5jptjet_py);
  tree->SetBranchAddress("ak5jptjet_pz", ak5jptjet_pz);
  tree->SetBranchAddress("ak5jptjet_hadronicenergy", ak5jptjet_hadronicenergy);
  tree->SetBranchAddress("ak5jptjet_chargedhadronicenergy", ak5jptjet_chargedhadronicenergy);
  tree->SetBranchAddress("ak5jptjet_emenergy", ak5jptjet_emenergy);
  tree->SetBranchAddress("ak5jptjet_chargedemenergy", ak5jptjet_chargedemenergy);
  tree->SetBranchAddress("ak5jptjet_chargedmulti", ak5jptjet_chargedmulti);
  tree->SetBranchAddress("ak5jptjet_energycorr", ak5jptjet_energycorr);
  tree->SetBranchAddress("ak5jptjet_fhpd", ak5jptjet_fhpd);
  tree->SetBranchAddress("ak5jptjet_restrictedemf", ak5jptjet_restrictedemf);
  tree->SetBranchAddress("ak5jptjet_btag", ak5jptjet_btag);
  tree->SetBranchAddress("ak5jptjet_n90", ak5jptjet_n90);
  
  tree->SetBranchAddress("ak5pfjet_count", &ak5pfjet_count);
  tree->SetBranchAddress("ak5pfjet_e", ak5pfjet_e);
  tree->SetBranchAddress("ak5pfjet_px", ak5pfjet_px);
  tree->SetBranchAddress("ak5pfjet_py", ak5pfjet_py);
  tree->SetBranchAddress("ak5pfjet_pz", ak5pfjet_pz);
  tree->SetBranchAddress("ak5pfjet_hadronicenergy", ak5pfjet_hadronicenergy);
  tree->SetBranchAddress("ak5pfjet_chargedhadronicenergy", ak5pfjet_chargedhadronicenergy);
  tree->SetBranchAddress("ak5pfjet_emenergy", ak5pfjet_emenergy);
  tree->SetBranchAddress("ak5pfjet_chargedemenergy", ak5pfjet_chargedemenergy);
  tree->SetBranchAddress("ak5pfjet_chargedmulti", ak5pfjet_chargedmulti);
  tree->SetBranchAddress("ak5pfjet_neutralmulti", ak5pfjet_neutralmulti);
  tree->SetBranchAddress("ak5pfjet_energycorr", ak5pfjet_energycorr);
  tree->SetBranchAddress("ak5pfjet_btag", ak5pfjet_btag);
  tree->SetBranchAddress("ak5pfjet_trigger", ak5pfjet_trigger);
  tree->SetBranchAddress("ak5pfjet_sv_count", ak5pfjet_sv_count);
  tree->SetBranchAddress("ak5pfjet_sv_x", ak5pfjet_sv_x);
  tree->SetBranchAddress("ak5pfjet_sv_y", ak5pfjet_sv_y);
  tree->SetBranchAddress("ak5pfjet_sv_z", ak5pfjet_sv_z);
  tree->SetBranchAddress("ak5pfjet_sv_dx", ak5pfjet_sv_dx);
  tree->SetBranchAddress("ak5pfjet_sv_dy", ak5pfjet_sv_dy);
  tree->SetBranchAddress("ak5pfjet_sv_dz", ak5pfjet_sv_dz);
  tree->SetBranchAddress("ak5pfjet_pu_jet_simple_loose", ak5pfjet_pu_jet_simple_loose);
  tree->SetBranchAddress("ak5pfjet_pu_jet_simple_medium", ak5pfjet_pu_jet_simple_medium);
  tree->SetBranchAddress("ak5pfjet_pu_jet_simple_tight", ak5pfjet_pu_jet_simple_tight);
  tree->SetBranchAddress("ak5pfjet_pu_jet_simple_mva", ak5pfjet_pu_jet_simple_mva);
  tree->SetBranchAddress("ak5pfjet_pu_jet_full_loose", ak5pfjet_pu_jet_full_loose);
  tree->SetBranchAddress("ak5pfjet_pu_jet_full_medium", ak5pfjet_pu_jet_full_medium);
  tree->SetBranchAddress("ak5pfjet_pu_jet_full_tight", ak5pfjet_pu_jet_full_tight);
  tree->SetBranchAddress("ak5pfjet_pu_jet_full_mva", ak5pfjet_pu_jet_full_mva);
  
  tree->SetBranchAddress("electron_count", &electron_count);
  tree->SetBranchAddress("electron_px", electron_px);
  tree->SetBranchAddress("electron_py", electron_py);
  tree->SetBranchAddress("electron_pz", electron_pz);
  tree->SetBranchAddress("electron_trackchi2", electron_trackchi2);
  tree->SetBranchAddress("electron_trackndof", electron_trackndof);
  tree->SetBranchAddress("electron_outerx", electron_outerx);
  tree->SetBranchAddress("electron_outery", electron_outery);
  tree->SetBranchAddress("electron_outerz", electron_outerz);
  tree->SetBranchAddress("electron_closestpointx", electron_closestpointx);
  tree->SetBranchAddress("electron_closestpointy", electron_closestpointy);
  tree->SetBranchAddress("electron_closestpointz", electron_closestpointz);
  tree->SetBranchAddress("electron_esuperclusterovertrack", electron_esuperclusterovertrack);
  tree->SetBranchAddress("electron_eseedclusterovertrack", electron_eseedclusterovertrack);
  tree->SetBranchAddress("electron_deltaetasuperclustertrack", electron_deltaetasuperclustertrack);
  tree->SetBranchAddress("electron_deltaphisuperclustertrack", electron_deltaphisuperclustertrack);
  tree->SetBranchAddress("electron_e1x5", electron_e1x5);
  tree->SetBranchAddress("electron_e2x5", electron_e2x5);
  tree->SetBranchAddress("electron_e5x5", electron_e5x5);
  tree->SetBranchAddress("electron_sigmaetaeta", electron_sigmaetaeta);
  tree->SetBranchAddress("electron_sigmaietaieta", electron_sigmaietaieta);
  tree->SetBranchAddress("electron_ehcaloverecal", electron_ehcaloverecal);
  tree->SetBranchAddress("electron_ehcaloverecaldepth1", electron_ehcaloverecaldepth1);
  tree->SetBranchAddress("electron_ehcaloverecaldepth2", electron_ehcaloverecaldepth2);
  tree->SetBranchAddress("electron_isolationr3track", electron_isolationr3track);
  tree->SetBranchAddress("electron_isolationr3ecal", electron_isolationr3ecal);
  tree->SetBranchAddress("electron_isolationr3hcal", electron_isolationr3hcal);
  tree->SetBranchAddress("electron_isolationr4track", electron_isolationr4track);
  tree->SetBranchAddress("electron_isolationr4ecal", electron_isolationr4ecal);
  tree->SetBranchAddress("electron_isolationr4hcal", electron_isolationr4hcal);
  tree->SetBranchAddress("electron_pfisolationr3_sumchargedhadronpt", electron_pfisolationr3_sumchargedhadronpt);
  tree->SetBranchAddress("electron_pfisolationr3_sumchargedparticlept", electron_pfisolationr3_sumchargedparticlept);
  tree->SetBranchAddress("electron_pfisolationr3_sumneutralhadronet", electron_pfisolationr3_sumneutralhadronet);
  tree->SetBranchAddress("electron_pfisolationr3_sumphotonet", electron_pfisolationr3_sumphotonet);
  tree->SetBranchAddress("electron_pfisolationr3_sumPUpt", electron_pfisolationr3_sumPUpt);
  tree->SetBranchAddress("electron_pfisolationr4_sumchargedhadronpt", electron_pfisolationr4_sumchargedhadronpt);
  tree->SetBranchAddress("electron_pfisolationr4_sumchargedparticlept", electron_pfisolationr4_sumchargedparticlept);
  tree->SetBranchAddress("electron_pfisolationr4_sumneutralhadronet", electron_pfisolationr4_sumneutralhadronet);
  tree->SetBranchAddress("electron_pfisolationr4_sumphotonet", electron_pfisolationr4_sumphotonet);
  tree->SetBranchAddress("electron_pfisolationr4_sumPUpt", electron_pfisolationr4_sumPUpt);
  tree->SetBranchAddress("electron_nhits", electron_nhits);
  tree->SetBranchAddress("electron_npixelhits", electron_npixelhits);
  tree->SetBranchAddress("electron_nmissinghits", electron_nmissinghits);
  tree->SetBranchAddress("electron_npixellayers", electron_npixellayers);
  tree->SetBranchAddress("electron_nstriplayers", electron_nstriplayers);
  tree->SetBranchAddress("electron_dxy", electron_dxy);
  tree->SetBranchAddress("electron_dxyerr", electron_dxyerr);
  tree->SetBranchAddress("electron_dz", electron_dz);
  tree->SetBranchAddress("electron_dzerr", electron_dzerr);
  tree->SetBranchAddress("electron_convdist", electron_convdist);
  tree->SetBranchAddress("electron_convdcot", electron_convdcot);
  tree->SetBranchAddress("electron_convradius", electron_convradius);
  tree->SetBranchAddress("electron_gapinfo", electron_gapinfo);
  tree->SetBranchAddress("electron_chargeinfo", electron_chargeinfo);
  tree->SetBranchAddress("electron_fbrems", electron_fbrems);
  tree->SetBranchAddress("electron_numbrems", electron_numbrems);
  tree->SetBranchAddress("electron_charge", electron_charge);
  tree->SetBranchAddress("electron_superclusterindex", electron_superclusterindex);
  tree->SetBranchAddress("electron_info", electron_info);
  tree->SetBranchAddress("electron_trigger", electron_trigger);
  tree->SetBranchAddress("electron_mva_id_trig", electron_mva_id_trig);
  tree->SetBranchAddress("electron_mva_id_nontrig", electron_mva_id_nontrig);
  tree->SetBranchAddress("electron_mva_iso", electron_mva_iso);
  tree->SetBranchAddress("electron_has_conversion", electron_has_conversion);
  
  tree->SetBranchAddress("photon_count", &photon_count);
  tree->SetBranchAddress("photon_px", photon_px);
  tree->SetBranchAddress("photon_py", photon_py);
  tree->SetBranchAddress("photon_pz", photon_pz);
  tree->SetBranchAddress("photon_e1x5", photon_e1x5);
  tree->SetBranchAddress("photon_e2x5", photon_e2x5);
  tree->SetBranchAddress("photon_e3x3", photon_e3x3);
  tree->SetBranchAddress("photon_e5x5", photon_e5x5);
  tree->SetBranchAddress("photon_sigmaetaeta", photon_sigmaetaeta);
  tree->SetBranchAddress("photon_sigmaietaieta", photon_sigmaietaieta);
  tree->SetBranchAddress("photon_ehcaloverecal", photon_ehcaloverecal);
  tree->SetBranchAddress("photon_ehcaloverecaldepth1", photon_ehcaloverecaldepth1);
  tree->SetBranchAddress("photon_ehcaloverecaldepth2", photon_ehcaloverecaldepth2);
  tree->SetBranchAddress("photon_maxenergyxtal", photon_maxenergyxtal);
  tree->SetBranchAddress("photon_isolationr3track", photon_isolationr3track);
  tree->SetBranchAddress("photon_isolationr3trackhollow", photon_isolationr3trackhollow);
  tree->SetBranchAddress("photon_isolationr3ntrack", photon_isolationr3ntrack);
  tree->SetBranchAddress("photon_isolationr3ntrackhollow", photon_isolationr3ntrackhollow);
  tree->SetBranchAddress("photon_isolationr3ecal", photon_isolationr3ecal);
  tree->SetBranchAddress("photon_isolationr3hcal", photon_isolationr3hcal);
  tree->SetBranchAddress("photon_isolationr4track", photon_isolationr4track);
  tree->SetBranchAddress("photon_isolationr4trackhollow", photon_isolationr4trackhollow);
  tree->SetBranchAddress("photon_isolationr4ntrack", photon_isolationr4ntrack);
  tree->SetBranchAddress("photon_isolationr4ntrackhollow", photon_isolationr4ntrackhollow);
  tree->SetBranchAddress("photon_isolationr4ecal", photon_isolationr4ecal);
  tree->SetBranchAddress("photon_isolationr4hcal", photon_isolationr4hcal);
  tree->SetBranchAddress("photon_superclusterindex", photon_superclusterindex);
  tree->SetBranchAddress("photon_info", photon_info);
  tree->SetBranchAddress("photon_gapinfo", photon_gapinfo);
  tree->SetBranchAddress("photon_conversionbegin", photon_conversionbegin);
  
  tree->SetBranchAddress("conversion_count", &conversion_count);
  tree->SetBranchAddress("conversion_info", conversion_info);
  tree->SetBranchAddress("conversion_vx", conversion_vx);
  tree->SetBranchAddress("conversion_vy", conversion_vy);
  tree->SetBranchAddress("conversion_vz", conversion_vz);
  tree->SetBranchAddress("conversion_chi2", conversion_chi2);
  tree->SetBranchAddress("conversion_ndof", conversion_ndof);
  tree->SetBranchAddress("conversion_cov", conversion_cov);
  tree->SetBranchAddress("conversion_mvaout", conversion_mvaout);
  tree->SetBranchAddress("conversion_trackecalpointx", conversion_trackecalpointx);
  tree->SetBranchAddress("conversion_trackecalpointy", conversion_trackecalpointy);
  tree->SetBranchAddress("conversion_trackecalpointz", conversion_trackecalpointz);
  tree->SetBranchAddress("conversion_trackpx", conversion_trackpx);
  tree->SetBranchAddress("conversion_trackpy", conversion_trackpy);
  tree->SetBranchAddress("conversion_trackpz", conversion_trackpz);
  tree->SetBranchAddress("conversion_trackclosestpointx", conversion_trackclosestpointx);
  tree->SetBranchAddress("conversion_trackclosestpointy", conversion_trackclosestpointy);
  tree->SetBranchAddress("conversion_trackclosestpointz", conversion_trackclosestpointz);
  tree->SetBranchAddress("conversion_trackchi2", conversion_trackchi2);
  tree->SetBranchAddress("conversion_trackndof", conversion_trackndof);
  tree->SetBranchAddress("conversion_trackdxy", conversion_trackdxy);
  tree->SetBranchAddress("conversion_trackdxyerr", conversion_trackdxyerr);
  tree->SetBranchAddress("conversion_trackdz", conversion_trackdz);
  tree->SetBranchAddress("conversion_trackdzerr", conversion_trackdzerr);
  tree->SetBranchAddress("conversion_trackcharge", conversion_trackcharge);
  tree->SetBranchAddress("conversion_tracknhits", conversion_tracknhits);
  tree->SetBranchAddress("conversion_tracknmissinghits", conversion_tracknmissinghits);
  tree->SetBranchAddress("conversion_tracknpixelhits", conversion_tracknpixelhits);
  tree->SetBranchAddress("conversion_tracknpixellayers", conversion_tracknpixellayers);
  tree->SetBranchAddress("conversion_tracknstriplayers", conversion_tracknstriplayers);
  // taus
  tree->SetBranchAddress("tau_count", &tau_count);
  tree->SetBranchAddress("tau_e", tau_e);
  tree->SetBranchAddress("tau_px", tau_px);
  tree->SetBranchAddress("tau_py", tau_py);
  tree->SetBranchAddress("tau_pz", tau_pz);
  tree->SetBranchAddress("tau_isolationneutralspt", tau_isolationneutralspt);
  tree->SetBranchAddress("tau_isolationneutralsnum", tau_isolationneutralsnum);
  tree->SetBranchAddress("tau_isolationchargedpt", tau_isolationchargedpt);
  tree->SetBranchAddress("tau_isolationchargednum", tau_isolationchargednum);
  tree->SetBranchAddress("tau_isolationgammapt", tau_isolationgammapt);
  tree->SetBranchAddress("tau_isolationgammanum", tau_isolationgammanum);
  tree->SetBranchAddress("tau_leadpfchargedhadrcandecalenergy", tau_leadpfchargedhadrcandecalenergy);
  tree->SetBranchAddress("tau_leadpfchargedhadrcandhcalenergy", tau_leadpfchargedhadrcandhcalenergy);
  tree->SetBranchAddress("tau_leadpfchargedhadrcandp", tau_leadpfchargedhadrcandp);
  tree->SetBranchAddress("tau_leadpfchargedhadrcandp4", &tau_leadpfchargedhadrcandp4);
  tree->SetBranchAddress("tau_leadpfchargedhadrcandpt", tau_leadpfchargedhadrcandpt);
  tree->SetBranchAddress("tau_dxy", tau_dxy);
  tree->SetBranchAddress("tau_dz", tau_dz);
  tree->SetBranchAddress("tau_ip3d", tau_ip3d);
  tree->SetBranchAddress("tau_ip3dSig", tau_ip3dSig);
  tree->SetBranchAddress("tau_nprongs", tau_nprongs);
  tree->SetBranchAddress("tau_charge", tau_charge);
  tree->SetBranchAddress("tau_newemfraction", tau_newemfraction);
  tree->SetBranchAddress("tau_hcaltotoverplead", tau_hcaltotoverplead);
  tree->SetBranchAddress("tau_hcal3x3overplead", tau_hcal3x3overplead);
  tree->SetBranchAddress("tau_ecalstripsumeoverplead", tau_ecalstripsumeoverplead);
  tree->SetBranchAddress("tau_bremsrecoveryeoverplead", tau_bremsrecoveryeoverplead);
  tree->SetBranchAddress("tau_calocomp", tau_calocomp);
  tree->SetBranchAddress("tau_segcomp", tau_segcomp);
  //  tree->SetBranchAddress("tau_disbyisolationmvaraw", tau_disbyisolationmvaraw);
  //  tree->SetBranchAddress("tau_disbyisolationmva2raw", tau_disbyisolationmva2raw);
  //  tree->SetBranchAddress("tau_againstelectronmva2raw", tau_againstelectronmva2raw);
  //  tree->SetBranchAddress("tau_againstelectronmva2category", tau_againstelectronmva2category);
  //  tree->SetBranchAddress("tau_againstelectronmva3raw", tau_againstelectronmva3raw);
  //  tree->SetBranchAddress("tau_againstelectronmva3category", tau_againstelectronmva3category);
  tree->SetBranchAddress("tau_bycombinedisolationdeltabetacorrraw3hits", tau_bycombinedisolationdeltabetacorrraw3hits);
  tree->SetBranchAddress("tau_dishps", tau_dishps);
  tree->SetBranchAddress("tau_trigger", tau_trigger);
  tree->SetBranchAddress("tau_L1trigger_match", tau_L1trigger_match);
  tree->SetBranchAddress("tau_signalPFChargedHadrCands_size", tau_signalPFChargedHadrCands_size);
  tree->SetBranchAddress("tau_signalPFGammaCands_size",tau_signalPFGammaCands_size);

  // ditau
  tree->SetBranchAddress("ditau_Index", &ditau_Index);
  tree->SetBranchAddress("ditau_leg1_index", ditau_leg1_index);
  tree->SetBranchAddress("ditau_leg2_index", ditau_leg2_index);
  tree->SetBranchAddress("ditau_VtxX", &ditau_VtxX);
  tree->SetBranchAddress("ditau_VtxY", &ditau_VtxY);
  tree->SetBranchAddress("ditau_VtxZ", &ditau_VtxZ);
  tree->SetBranchAddress("ditau_VtxXErr", &ditau_VtxXErr);
  tree->SetBranchAddress("ditau_VtxYErr", &ditau_VtxYErr);
  tree->SetBranchAddress("ditau_VtxZErr", &ditau_VtxZErr);
  tree->SetBranchAddress("ditau_reFitVtxX", ditau_reFitVtxX);
  tree->SetBranchAddress("ditau_reFitVtxY", ditau_reFitVtxY);
  tree->SetBranchAddress("ditau_reFitVtxZ", ditau_reFitVtxZ);
  tree->SetBranchAddress("ditau_reFitVtxXErr", ditau_reFitVtxXErr);
  tree->SetBranchAddress("ditau_reFitVtxYErr", ditau_reFitVtxYErr);
  tree->SetBranchAddress("ditau_reFitVtxZErr", ditau_reFitVtxZErr);
  tree->SetBranchAddress("ditau_reFitVtxRho", ditau_reFitVtxRho);
  tree->SetBranchAddress("ditau_leg1_dxy", ditau_leg1_dxy);
  tree->SetBranchAddress("ditau_leg1_dz", ditau_leg1_dz);
  tree->SetBranchAddress("ditau_leg1_dxyErr", ditau_leg1_dxyErr);
  tree->SetBranchAddress("ditau_leg1_dzErr", ditau_leg1_dzErr);
  tree->SetBranchAddress("ditau_leg1_dxy_OPV", ditau_leg1_dxy_OPV);
  tree->SetBranchAddress("ditau_leg1_dz_OPV", ditau_leg1_dz_OPV);
  tree->SetBranchAddress("ditau_leg1_dxyErr_OPV", ditau_leg1_dxyErr_OPV);
  tree->SetBranchAddress("ditau_leg1_dzErr_OPV", ditau_leg1_dzErr_OPV);
  tree->SetBranchAddress("ditau_leg1_X", ditau_leg1_X);
  tree->SetBranchAddress("ditau_leg1_Y", ditau_leg1_Y);
  tree->SetBranchAddress("ditau_leg1_Z", ditau_leg1_Z);
  tree->SetBranchAddress("ditau_leg1_X_OPV", ditau_leg1_X_OPV);
  tree->SetBranchAddress("ditau_leg1_Y_OPV", ditau_leg1_Y_OPV);
  tree->SetBranchAddress("ditau_leg1_Z_OPV", ditau_leg1_Z_OPV);
  tree->SetBranchAddress("ditau_leg2_dxy", ditau_leg2_dxy);
  tree->SetBranchAddress("ditau_leg2_dz", ditau_leg2_dz);
  tree->SetBranchAddress("ditau_leg2_dxyErr", ditau_leg2_dxyErr);
  tree->SetBranchAddress("ditau_leg2_dzErr", ditau_leg2_dzErr);
  tree->SetBranchAddress("ditau_leg2_dxy_OPV", ditau_leg2_dxy_OPV);
  tree->SetBranchAddress("ditau_leg2_dz_OPV", ditau_leg2_dz_OPV);
  tree->SetBranchAddress("ditau_leg2_dxyErr_OPV", ditau_leg2_dxyErr_OPV);
  tree->SetBranchAddress("ditau_leg2_dzErr_OPV", ditau_leg2_dzErr_OPV);
  tree->SetBranchAddress("ditau_leg2_X", ditau_leg2_X);
  tree->SetBranchAddress("ditau_leg2_Y", ditau_leg2_Y);
  tree->SetBranchAddress("ditau_leg2_Z", ditau_leg2_Z);
  tree->SetBranchAddress("ditau_leg2_X_OPV", ditau_leg2_X_OPV);
  tree->SetBranchAddress("ditau_leg2_Y_OPV", ditau_leg2_Y_OPV);
  tree->SetBranchAddress("ditau_leg2_Z_OPV", ditau_leg2_Z_OPV);

  // tau pairs
  tree->SetBranchAddress("mutautaupair_count", &mutautaupair_count);
  tree->SetBranchAddress("mutautaupair_leg1_px", mutautaupair_leg1_px);
  tree->SetBranchAddress("mutautaupair_leg1_py", mutautaupair_leg1_py);
  tree->SetBranchAddress("mutautaupair_leg1_pz", mutautaupair_leg1_pz);
  tree->SetBranchAddress("mutautaupair_leg1_energy", mutautaupair_leg1_energy);
  tree->SetBranchAddress("mutautaupair_leg2_px", mutautaupair_leg2_px);
  tree->SetBranchAddress("mutautaupair_leg2_py", mutautaupair_leg2_py);
  tree->SetBranchAddress("mutautaupair_leg2_pz", mutautaupair_leg2_pz);
  tree->SetBranchAddress("mutautaupair_leg2_energy", mutautaupair_leg2_energy);
  tree->SetBranchAddress("mutautaupair_mu_px", mutautaupair_mu_px);
  tree->SetBranchAddress("mutautaupair_mu_py", mutautaupair_mu_py);
  tree->SetBranchAddress("mutautaupair_mu_pz", mutautaupair_mu_pz);
  tree->SetBranchAddress("mutautaupair_mu_energy", mutautaupair_mu_energy);
  tree->SetBranchAddress("mutautaupair_svfit_int_valid", mutautaupair_svfit_int_valid);
  tree->SetBranchAddress("mutautaupair_svfit_mass_int", mutautaupair_svfit_mass_int);
  tree->SetBranchAddress("mutautaupair_svfit_mass_int_err_up", mutautaupair_svfit_mass_int_err_up);
  tree->SetBranchAddress("mutautaupair_svfit_mass_int_err_down", mutautaupair_svfit_mass_int_err_down);
  tree->SetBranchAddress("mutautaupair_dca2d", mutautaupair_dca2d);
  tree->SetBranchAddress("mutautaupair_dca2d_err", mutautaupair_dca2d_err);
  tree->SetBranchAddress("mutautaupair_dca3d", mutautaupair_dca3d);
  tree->SetBranchAddress("mutautaupair_dca3d_err", mutautaupair_dca3d_err);

  tree->SetBranchAddress("eltautaupair_count", &eltautaupair_count);
  tree->SetBranchAddress("eltautaupair_leg1_px", eltautaupair_leg1_px);
  tree->SetBranchAddress("eltautaupair_leg1_py", eltautaupair_leg1_py);
  tree->SetBranchAddress("eltautaupair_leg1_pz", eltautaupair_leg1_pz);
  tree->SetBranchAddress("eltautaupair_leg1_energy", eltautaupair_leg1_energy);
  tree->SetBranchAddress("eltautaupair_leg2_px", eltautaupair_leg2_px);
  tree->SetBranchAddress("eltautaupair_leg2_py", eltautaupair_leg2_py);
  tree->SetBranchAddress("eltautaupair_leg2_pz", eltautaupair_leg2_pz);
  tree->SetBranchAddress("eltautaupair_leg2_energy", eltautaupair_leg2_energy);
  tree->SetBranchAddress("eltautaupair_el_px", eltautaupair_el_px);
  tree->SetBranchAddress("eltautaupair_el_py", eltautaupair_el_py);
  tree->SetBranchAddress("eltautaupair_el_pz", eltautaupair_el_pz);
  tree->SetBranchAddress("eltautaupair_el_energy", eltautaupair_el_energy);
  tree->SetBranchAddress("eltautaupair_svfit_int_valid", eltautaupair_svfit_int_valid);
  tree->SetBranchAddress("eltautaupair_svfit_mass_int", eltautaupair_svfit_mass_int);
  tree->SetBranchAddress("eltautaupair_svfit_mass_int_err_up", eltautaupair_svfit_mass_int_err_up);
  tree->SetBranchAddress("eltautaupair_svfit_mass_int_err_down", eltautaupair_svfit_mass_int_err_down);
  tree->SetBranchAddress("eltautaupair_dca2d", eltautaupair_dca2d);
  tree->SetBranchAddress("eltautaupair_dca2d_err", eltautaupair_dca2d_err);
  tree->SetBranchAddress("eltautaupair_dca3d", eltautaupair_dca3d);
  tree->SetBranchAddress("eltautaupair_dca3d_err", eltautaupair_dca3d_err);

	
  tree->SetBranchAddress("tau_ak5pfjet_e", tau_ak5pfjet_e);
  tree->SetBranchAddress("tau_ak5pfjet_px", tau_ak5pfjet_px);
  tree->SetBranchAddress("tau_ak5pfjet_py", tau_ak5pfjet_py);
  tree->SetBranchAddress("tau_ak5pfjet_pz", tau_ak5pfjet_pz);

  tree->SetBranchAddress("tau_chargedbegin", tau_chargedbegin);
  tree->SetBranchAddress("tau_charged_count", &tau_charged_count);
  tree->SetBranchAddress("tau_charged_px", tau_charged_px);
  tree->SetBranchAddress("tau_charged_py", tau_charged_py);
  tree->SetBranchAddress("tau_charged_pz", tau_charged_pz);
  tree->SetBranchAddress("tau_charged_outerx", tau_charged_outerx);
  tree->SetBranchAddress("tau_charged_outery", tau_charged_outery);
  tree->SetBranchAddress("tau_charged_outerz", tau_charged_outerz);
  tree->SetBranchAddress("tau_charged_closestpointx", tau_charged_closestpointx);
  tree->SetBranchAddress("tau_charged_closestpointy", tau_charged_closestpointy);
  tree->SetBranchAddress("tau_charged_closestpointz", tau_charged_closestpointz);
  tree->SetBranchAddress("tau_charged_chi2", tau_charged_chi2);
  tree->SetBranchAddress("tau_charged_ndof", tau_charged_ndof);
  tree->SetBranchAddress("tau_charged_dxy", tau_charged_dxy);
  tree->SetBranchAddress("tau_charged_dxyerr", tau_charged_dxyerr);
  tree->SetBranchAddress("tau_charged_dz", tau_charged_dz);
  tree->SetBranchAddress("tau_charged_dzerr", tau_charged_dzerr);
  tree->SetBranchAddress("tau_charged_dedxharmonic2", tau_charged_dedxharmonic2);
  tree->SetBranchAddress("tau_charged_charge", tau_charged_charge);
  tree->SetBranchAddress("tau_charged_nhits", tau_charged_nhits);
  tree->SetBranchAddress("tau_charged_nmissinghits", tau_charged_nmissinghits);
  tree->SetBranchAddress("tau_charged_npixelhits", tau_charged_npixelhits);
  tree->SetBranchAddress("tau_charged_npixellayers", tau_charged_npixellayers);
  tree->SetBranchAddress("tau_charged_nstriplayers", tau_charged_nstriplayers);
  
	tree->SetBranchAddress("ak5pfjet_rho", &ak5pfjet_rho);
	tree->SetBranchAddress("ak5pfjet_sigma", &ak5pfjet_sigma);

	tree->SetBranchAddress("calomet_ex", &calomet_ex);
	tree->SetBranchAddress("calomet_ey", &calomet_ey);
	tree->SetBranchAddress("calomet_exmuons", &calomet_exmuons);
	tree->SetBranchAddress("calomet_eymuons", &calomet_eymuons);

	tree->SetBranchAddress("tcmet_ex", &tcmet_ex);
	tree->SetBranchAddress("tcmet_ey", &tcmet_ey);

	tree->SetBranchAddress("pfmet_ex", &pfmet_ex);
	tree->SetBranchAddress("pfmet_ey", &pfmet_ey);
	tree->SetBranchAddress("pfmettype1_ex", &pfmettype1_ex);
	tree->SetBranchAddress("pfmettype1_ey", &pfmettype1_ey);
	tree->SetBranchAddress("pfmetmva_ex", &pfmetmva_ex);
	tree->SetBranchAddress("pfmetmva_ey", &pfmetmva_ey);

	tree->SetBranchAddress("secvertices_count", &secvertices_count);
	tree->SetBranchAddress("secvertices_vx", secvertices_vx);
	tree->SetBranchAddress("secvertices_vy", secvertices_vy);
	tree->SetBranchAddress("secvertices_vz", secvertices_vz);
	tree->SetBranchAddress("secvertices_chi2", secvertices_chi2);
	tree->SetBranchAddress("secvertices_ndof", secvertices_ndof);
	tree->SetBranchAddress("secvertices_cov", secvertices_cov);
	tree->SetBranchAddress("secvertices_track_px", secvertices_track_px);
	tree->SetBranchAddress("secvertices_track_py", secvertices_track_py);
	tree->SetBranchAddress("secvertices_track_pz", secvertices_track_pz);
	tree->SetBranchAddress("secvertices_track_closestpointx", secvertices_track_closestpointx);
	tree->SetBranchAddress("secvertices_track_closestpointy", secvertices_track_closestpointy);
	tree->SetBranchAddress("secvertices_track_closestpointz", secvertices_track_closestpointz);
	tree->SetBranchAddress("secvertices_track_chi2", secvertices_track_chi2);
	tree->SetBranchAddress("secvertices_track_ndof", secvertices_track_ndof);
	tree->SetBranchAddress("secvertices_track_dxy", secvertices_track_dxy);
	tree->SetBranchAddress("secvertices_track_dxyerr", secvertices_track_dxyerr);
	tree->SetBranchAddress("secvertices_track_dz", secvertices_track_dz);
	tree->SetBranchAddress("secvertices_track_dzerr", secvertices_track_dzerr);
	tree->SetBranchAddress("secvertices_track_dedxharmonic2", secvertices_track_dedxharmonic2);
	tree->SetBranchAddress("secvertices_track_charge", secvertices_track_charge);
	tree->SetBranchAddress("secvertices_track_nhits", secvertices_track_nhits);
	tree->SetBranchAddress("secvertices_track_nmissinghits", secvertices_track_nmissinghits);
	tree->SetBranchAddress("secvertices_track_npixelhits", secvertices_track_npixelhits);
	tree->SetBranchAddress("secvertices_track_npixellayers", secvertices_track_npixellayers);
	tree->SetBranchAddress("secvertices_track_nstriplayers", secvertices_track_nstriplayers);

	tree->SetBranchAddress("genweight", &genweight);
	tree->SetBranchAddress("genx1", &genx1);
	tree->SetBranchAddress("genid1", &genid1);
	tree->SetBranchAddress("genx2", &genx2);
	tree->SetBranchAddress("genid2", &genid2);
	tree->SetBranchAddress("genScale", &genScale);

	tree->SetBranchAddress("numpileupinteractionsminus", &numpileupinteractionsminus);
	tree->SetBranchAddress("numpileupinteractions", &numpileupinteractions);
	tree->SetBranchAddress("numpileupinteractionsplus", &numpileupinteractionsplus);
	tree->SetBranchAddress("numtruepileupinteractions", &numtruepileupinteractions);

	tree->SetBranchAddress("genparticles_count", &genparticles_count);
	tree->SetBranchAddress("genparticles_e", genparticles_e);
	tree->SetBranchAddress("genparticles_px", genparticles_px);
	tree->SetBranchAddress("genparticles_py", genparticles_py);
	tree->SetBranchAddress("genparticles_pz", genparticles_pz);
	tree->SetBranchAddress("genparticles_vx", genparticles_vx);
	tree->SetBranchAddress("genparticles_vy", genparticles_vy);
	tree->SetBranchAddress("genparticles_vz", genparticles_vz);
	tree->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
	tree->SetBranchAddress("genparticles_status", genparticles_status);
	tree->SetBranchAddress("genparticles_info", genparticles_info);

	tree->SetBranchAddress("genallparticles_count", &genallparticles_count);
	tree->SetBranchAddress("genallparticles_e", genallparticles_e);
	tree->SetBranchAddress("genallparticles_px", genallparticles_px);
	tree->SetBranchAddress("genallparticles_py", genallparticles_py);
	tree->SetBranchAddress("genallparticles_pz", genallparticles_pz);
	tree->SetBranchAddress("genallparticles_vx", genallparticles_vx);
	tree->SetBranchAddress("genallparticles_vy", genallparticles_vy);
	tree->SetBranchAddress("genallparticles_vz", genallparticles_vz);
	tree->SetBranchAddress("genallparticles_pdgid", genallparticles_pdgid);
	tree->SetBranchAddress("genallparticles_charge", genallparticles_charge);
	tree->SetBranchAddress("genallparticles_status", genallparticles_status);
	tree->SetBranchAddress("genallparticles_motherbeg", genallparticles_motherbeg);
	tree->SetBranchAddress("genallparticles_daughterbeg", genallparticles_daughterbeg);

	tree->SetBranchAddress("genallparticlesmother_count", &genallparticlesmother_count);
	tree->SetBranchAddress("genallparticles_mothers", genallparticles_mothers);

	tree->SetBranchAddress("genallparticlesdaughter_count", &genallparticlesdaughter_count);
	tree->SetBranchAddress("genallparticles_daughters", genallparticles_daughters);

	tree->SetBranchAddress("genmetcalo_ex", &genmetcalo_ex);
	tree->SetBranchAddress("genmetcalo_ey", &genmetcalo_ey);
	tree->SetBranchAddress("genmettrue_ex", &genmettrue_ex);
	tree->SetBranchAddress("genmettrue_ey", &genmettrue_ey);
}

void Analyse::UsePileUpInfo()
{
  usepileupinfo = true;
}

void Analyse::UsePrimVertexInfo()
{
	useprimvertexinfo = true;
}

void Analyse::LoadBeamSpot(bool select)
{
if(select){loadbeamspot = 1;}else{loadbeamspot = 0;}
}

void Analyse::LoadSecVertices(bool select)
{
if(select){loadsecvertices = 1;}else{loadsecvertices = 0;}
}

void Analyse::LoadMuons(bool select)
{
if(select){loadmuons = 1;}else{loadmuons = 0;}
}

void Analyse::LoadElectrons(bool select)
{
if(select){loadelectrons = 1;}else{loadelectrons = 0;}
}

void Analyse::LoadPhotons(bool select)
{
if(select){loadphotons = 1;}else{loadphotons = 0;}
}

void Analyse::LoadTaus(bool select)
{
  if(select){
    loadtaus = 1;
  }
  else {
    loadtaus = 0;
  }
}

void Analyse::LoadMET(bool select)
{
if(select){loadmet = 1;}else{loadmet = 0;}
}

void Analyse::LoadAK5CaloJets(bool select)
{
if(select){loadak5calojets = 1;}else{loadak5calojets = 0;}
}

void Analyse::LoadAK5JPTJets(bool select)
{
if(select){loadak5jptjets = 1;}else{loadak5jptjets = 0;}
}

void Analyse::LoadAK5PFJets(bool select)
{
if(select){loadak5pfjets = 1;}else{loadak5pfjets = 0;}
}

void Analyse::LoadTracks(bool select)
{
if(select){loadtracks = 1;}else{loadtracks = 0;}
}

void Analyse::LoadPrimVertices(bool select)
{
if(select){loadprimvertices = 1;}else{loadprimvertices = 0;}
}

void Analyse::LoadTrigger(bool select)
{
  if(select) loadtrigger = 1;
  else loadtrigger = 0;
}

void Analyse::LoadGenInfo(bool select)
{
if(select){loadgeninfo = 1;}else{loadgeninfo = 0;}
}

void Analyse::LoadAllGenParticles(bool select)
{
if(select){loadallgenparticles = 1;}else{loadallgenparticles = 0;}
}

void Analyse::LoadGenParticles(bool select)
{
if(select){loadgenparticles = 1;}else{loadgenparticles = 0;}
}

void Analyse::LoadSuperClusters(bool select, bool usebasiccluster, bool usebasicclusterhit)
{
	if(!select){loadsuperclusters = 1; return;}
	loadsuperclusters = 1;
	if(usebasicclusterhit)
	{
		loadsuperclusters = 3;
	}
	if(usebasiccluster)
	{
		loadsuperclusters = 2;
	}
	if(usebasicclusterhit && usebasiccluster)
	{
		loadsuperclusters = 4;
	}
}


bool Analyse::GetL1Trigger(UInt_t bit) const
{
	if((trigger_level1[bit/8] & 1<<(bit % 8)) > 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

bool Analyse::GetL1TriggerBits(UInt_t bit) const
{
	if((trigger_level1bits[bit/8] & 1<<(bit % 8)) > 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

bool Analyse::GetHLTrigger(UInt_t index) const
{
	if((trigger_HLT[index/8] & 1<<(index % 8)) > 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

BeamSpot Analyse::GetBeamSpot() const
{
	return(BeamSpot(beamspot_x, beamspot_y, beamspot_z, beamspot_cov, beamspot_xwidth, beamspot_ywidth, beamspot_zsigma));
}

TrackComposedParticle Analyse::SecVertices(UInt_t n) const
{
  TrackComposedParticle newsecvertices = 
    TrackComposedParticle(secvertices_vx[n], secvertices_vy[n], secvertices_vz[n], secvertices_chi2[n], secvertices_ndof[n], secvertices_cov[n]);
  for(Int_t i = 0 ; i < 2 ; i++)
    {
      Double_t E = sqrt(secvertices_track_px[n][i]*secvertices_track_px[n][i] + 
			secvertices_track_py[n][i]*secvertices_track_py[n][i] + 
			secvertices_track_pz[n][i]*secvertices_track_pz[n][i]);
      Track newtrack(E, secvertices_track_px[n][i], secvertices_track_py[n][i], secvertices_track_pz[n][i], 0., 0., 0., 
		     secvertices_track_closestpointx[n][i], secvertices_track_closestpointy[n][i], secvertices_track_closestpointz[n][i], 
		     secvertices_track_chi2[n][i], secvertices_track_ndof[n][i], secvertices_track_dxy[n][i], secvertices_track_dxyerr[n][i], 
		     secvertices_track_dz[n][i], secvertices_track_dzerr[n][i], secvertices_track_charge[n][i], secvertices_track_nhits[n][i], 
		     secvertices_track_nmissinghits[n][i], secvertices_track_npixelhits[n][i], secvertices_track_npixellayers[n][i], 
		     secvertices_track_nstriplayers[n][i], secvertices_track_dedxharmonic2[n][i]);
      newsecvertices.AddTrack(newtrack);
    }
  return(newsecvertices);
}

Muon Analyse::Muons(UInt_t n) const
{
  //vector<string> a = runlist.find(Run())->second.GetHLTMuonNames();
  Muon newmuon(this, muon_px[n], muon_py[n], muon_pz[n], muon_pterror[n], muon_chi2[n], muon_ndof[n],muon_pfisolationr4_sumchargedhadronpt[n],muon_pfisolationr4_sumchargedparticlept[n],muon_pfisolationr4_sumneutralhadronet[n],muon_pfisolationr4_sumphotonet[n],muon_pfisolationr4_sumPUpt[n], muon_isolationr3track[n], muon_isolationr3ntrack[n], muon_isolationr3ecal[n], muon_isolationr3hcal[n], muon_ecalenergy[n], muon_hcalenergy[n], muon_charge[n], muon_numchambers[n], muon_numchamberswithsegments[n], muon_nchambers[n], muon_type[n], muon_trigger[n], muon_trackermuonquality[n], muon_mva_id[n], muon_mva_iso[n]);
  
  Track Outertrack(sqrt(pow(muon_outertrack_px[n],2) + pow(muon_outertrack_py[n],2) + pow(muon_outertrack_pz[n],2)), muon_outertrack_px[n], muon_outertrack_py[n], muon_outertrack_pz[n], 10000., 10000., 10000., 1000., 1000., 1000., muon_outertrack_chi2[n], muon_outertrack_ndof[n], 0., 0., 0., 0., muon_charge[n], muon_outertrack_hits[n], muon_outertrack_missinghits[n], 0, 0, 0, 0.);
  newmuon.SetOuterTrack(Outertrack);
  
  Track Innertrack(sqrt(muon_innertrack_px[n]*muon_innertrack_px[n]+muon_innertrack_py[n]*muon_innertrack_py[n]+muon_innertrack_pz[n]*muon_innertrack_pz[n]), muon_innertrack_px[n], muon_innertrack_py[n], muon_innertrack_pz[n], muon_innertrack_outerx[n], muon_innertrack_outery[n], muon_innertrack_outerz[n], muon_innertrack_closestpointx[n], muon_innertrack_closestpointy[n], muon_innertrack_closestpointz[n], muon_innertrack_chi2[n], muon_innertrack_ndof[n], muon_innertrack_dxy[n], muon_innertrack_dxyerr[n], muon_innertrack_dz[n], muon_innertrack_dzerr[n], muon_innertrack_charge[n], muon_innertrack_nhits[n], muon_innertrack_nmissinghits[n], muon_innertrack_npixelhits[n], muon_innertrack_npixellayers[n], muon_innertrack_nstriplayers[n], muon_innertrack_dedxharmonic2[n]);
  newmuon.SetInnerTrack(Innertrack);
  return(newmuon);
}


Electron Analyse::Electrons(UInt_t n) const
{
  Electron newelectron(this, electron_px[n], electron_py[n], electron_pz[n], electron_trackchi2[n], electron_trackndof[n], electron_outerx[n], electron_outery[n], electron_outerz[n], electron_closestpointx[n], electron_closestpointy[n], electron_closestpointz[n], electron_dxy[n], electron_dxyerr[n], electron_dz[n], electron_dzerr[n], electron_esuperclusterovertrack[n], electron_eseedclusterovertrack[n], electron_deltaetasuperclustertrack[n], electron_deltaphisuperclustertrack[n], electron_e1x5[n], electron_e2x5[n], electron_e5x5[n], electron_sigmaetaeta[n], electron_sigmaietaieta[n], electron_ehcaloverecal[n], electron_ehcaloverecaldepth1[n], electron_ehcaloverecaldepth2[n], electron_isolationr3track[n], electron_isolationr3ecal[n], electron_isolationr3hcal[n], electron_isolationr4track[n], electron_isolationr4ecal[n], electron_isolationr4hcal[n], electron_pfisolationr3_sumchargedhadronpt[n] ,electron_pfisolationr3_sumchargedparticlept[n] ,electron_pfisolationr3_sumneutralhadronet[n] ,electron_pfisolationr3_sumphotonet[n] ,electron_pfisolationr3_sumPUpt[n] , electron_pfisolationr4_sumchargedhadronpt[n] ,electron_pfisolationr4_sumchargedparticlept[n] ,electron_pfisolationr4_sumneutralhadronet[n] ,electron_pfisolationr4_sumphotonet[n] ,electron_pfisolationr4_sumPUpt[n] , electron_nhits[n], electron_nmissinghits[n], electron_npixelhits[n], electron_npixellayers[n], electron_nstriplayers[n], electron_convdist[n], electron_convdcot[n], electron_convradius[n], electron_gapinfo[n], electron_chargeinfo[n], electron_fbrems[n], electron_numbrems[n], electron_charge[n], electron_info[n], electron_trigger[n], electron_mva_id_nontrig[n], electron_has_conversion[n]);
  
  if(loadsuperclusters > 0 && electron_superclusterindex[n] != -1)
    {
      SuperCluster newsc = SuperClusters(electron_superclusterindex[n]);
      newelectron.AddSC(newsc);
    }
  
  return(newelectron);
}

Photon Analyse::Photons(UInt_t n) const
{
  Photon newphoton(this, photon_px[n], photon_py[n], photon_pz[n], photon_e1x5[n], photon_e2x5[n], photon_e3x3[n], photon_e5x5[n], photon_sigmaetaeta[n], photon_sigmaietaieta[n], photon_ehcaloverecal[n], photon_ehcaloverecaldepth1[n], photon_ehcaloverecaldepth2[n], photon_maxenergyxtal[n], photon_isolationr3track[n], photon_isolationr3trackhollow[n], photon_isolationr3ntrack[n], photon_isolationr3ntrackhollow[n], photon_isolationr3ecal[n], photon_isolationr3hcal[n], photon_isolationr4track[n], photon_isolationr4trackhollow[n], photon_isolationr4ntrack[n], photon_isolationr4ntrackhollow[n], photon_isolationr4ecal[n], photon_isolationr4hcal[n], photon_info[n], photon_gapinfo[n]); 
  
  UInt_t begin = photon_conversionbegin[n];
  UInt_t end = conversion_count;
  if(n < photon_count - 1)
    {
      end = photon_conversionbegin[n+1];
    }
  
  for(UInt_t i = begin ; i < end ; i++)
    {
      Conversion newconversion(conversion_vx[i], conversion_vy[i], conversion_vz[i], conversion_chi2[i], conversion_ndof[i], conversion_cov[i], conversion_mvaout[i], TVector3(conversion_trackecalpointx[i][0], conversion_trackecalpointy[i][0], conversion_trackecalpointz[i][0]), TVector3(conversion_trackecalpointx[i][1], conversion_trackecalpointy[i][1], conversion_trackecalpointz[i][1]));
      if(conversion_trackndof[i][0] != -1)
	{
	  Double_t E = sqrt(conversion_trackpx[i][0]*conversion_trackpx[i][0]+conversion_trackpy[i][0]*conversion_trackpy[i][0]+conversion_trackpz[i][0]*conversion_trackpz[i][0]);
	  Track newtrack(E, conversion_trackpx[i][0], conversion_trackpy[i][0], conversion_trackpz[i][0], conversion_trackecalpointx[i][0], conversion_trackecalpointy[i][0], conversion_trackecalpointz[i][0], conversion_trackclosestpointx[i][0], conversion_trackclosestpointy[i][0], conversion_trackclosestpointz[i][0], conversion_trackchi2[i][0], conversion_trackndof[i][0], conversion_trackdxy[i][0], conversion_trackdxyerr[i][0], conversion_trackdz[i][0], conversion_trackdzerr[i][0], conversion_trackcharge[i][0], conversion_tracknhits[i][0], conversion_tracknmissinghits[i][0], conversion_tracknpixelhits[i][0], conversion_tracknpixellayers[i][0], conversion_tracknstriplayers[i][0], 0.);
	  newconversion.AddTrack(newtrack);
	}
      if(conversion_trackndof[i][1] != -1)
	{
	  Double_t E = sqrt(conversion_trackpx[i][1]*conversion_trackpx[i][1]+conversion_trackpy[i][1]*conversion_trackpy[i][1]+conversion_trackpz[i][1]*conversion_trackpz[i][1]);
	  Track newtrack(E, conversion_trackpx[i][1], conversion_trackpy[i][1], conversion_trackpz[i][1], conversion_trackecalpointx[i][1], conversion_trackecalpointy[i][1], conversion_trackecalpointz[i][1], conversion_trackclosestpointx[i][1], conversion_trackclosestpointy[i][1], conversion_trackclosestpointz[i][1], conversion_trackchi2[i][1], conversion_trackndof[i][1], conversion_trackdxy[i][1], conversion_trackdxyerr[i][1], conversion_trackdz[i][1], conversion_trackdzerr[i][1], conversion_trackcharge[i][1], conversion_tracknhits[i][1], conversion_tracknmissinghits[i][1], conversion_tracknpixelhits[i][1], conversion_tracknpixellayers[i][1], conversion_tracknstriplayers[i][1], 0.);
	  newconversion.AddTrack(newtrack);
	}
      newphoton.AddConversion(newconversion);
      
    } 
  if(loadsuperclusters > 0 && photon_superclusterindex[n] != -1)
    {
      SuperCluster newsc = SuperClusters(photon_superclusterindex[n]);
      newphoton.AddSC(newsc);
    }
  
  return(newphoton);
}

Tau Analyse::Taus(UInt_t n) const{
  return(Tau(this, n));
}

diTau Analyse::diTaus(UInt_t n) const{
  return(diTau(this, n));
}


MuTauTauPair Analyse::MuTauTauPairs(UInt_t n) const{
  return(MuTauTauPair(this, n));
}

ElTauTauPair Analyse::ElTauTauPairs(UInt_t n) const{
  return(ElTauTauPair(this, n));
}

TLorentzVector Analyse::CaloMET() const
{
	return(TLorentzVector(calomet_ex, calomet_ey, 0., sqrt(pow(calomet_ex, 2) + pow(calomet_ey, 2))));
}

TLorentzVector Analyse::CaloMETMuons() const
{
	return(TLorentzVector(calomet_exmuons, calomet_eymuons, 0., sqrt(pow(calomet_exmuons, 2) + pow(calomet_eymuons, 2))));
}

TLorentzVector Analyse::PFMET() const
{
  return(TLorentzVector(pfmet_ex, pfmet_ey, 0., sqrt(pow(pfmet_ex, 2) + pow(pfmet_ey, 2))));
}

TLorentzVector Analyse::PFMETTYPE1() const
{
  return(TLorentzVector(pfmettype1_ex, pfmettype1_ey, 0., sqrt(pow(pfmettype1_ex, 2) + pow(pfmettype1_ey, 2))));
}

TLorentzVector Analyse::PFMETMVA() const
{
  return(TLorentzVector(pfmetmva_ex, pfmetmva_ey, 0., sqrt(pow(pfmetmva_ex, 2) + pow(pfmetmva_ey, 2))) );
}

TLorentzVector Analyse::TCMET() const
{
  return(TLorentzVector(tcmet_ex, tcmet_ey, 0., sqrt(pow(tcmet_ex, 2) + pow(tcmet_ey, 2))));
}

Jet Analyse::AK5CaloJets(UInt_t n) const
{
  return(Jet(ak5calojet_e[n], ak5calojet_px[n], ak5calojet_py[n], ak5calojet_pz[n], ak5calojet_hadronicenergy[n], -1., ak5calojet_emenergy[n], -1., -1, -1, ak5calojet_energycorr[n], ak5calojet_btag[n], ak5calojet_n90[n], ak5calojet_n60[n], ak5calojet_fhpd[n], ak5calojet_restrictedemf[n], false, false, false, 0.0f, false, false, false, 0.0f, 0));
}

Jet Analyse::AK5JPTJets(UInt_t n) const
{
  return(Jet(ak5jptjet_e[n], ak5jptjet_px[n], ak5jptjet_py[n], ak5jptjet_pz[n], ak5jptjet_hadronicenergy[n], ak5jptjet_chargedhadronicenergy[n], ak5jptjet_emenergy[n], ak5jptjet_chargedemenergy[n], ak5jptjet_chargedmulti[n], -1, ak5jptjet_energycorr[n], ak5jptjet_btag[n], ak5jptjet_n90[n], -1, ak5jptjet_fhpd[n], ak5jptjet_restrictedemf[n], false, false, false, 0.0f, false, false, false, 0.0f, 0));
}

Jet Analyse::AK5PFJets(UInt_t n) const
{
  return(Jet(ak5pfjet_e[n], ak5pfjet_px[n], ak5pfjet_py[n], ak5pfjet_pz[n], ak5pfjet_hadronicenergy[n], ak5pfjet_chargedhadronicenergy[n], ak5pfjet_emenergy[n], ak5pfjet_chargedemenergy[n], ak5pfjet_chargedmulti[n], ak5pfjet_neutralmulti[n], ak5pfjet_energycorr[n], ak5pfjet_btag[n], -1, -1, -1., -1., ak5pfjet_pu_jet_simple_loose[n], ak5pfjet_pu_jet_simple_medium[n], ak5pfjet_pu_jet_simple_tight[n], ak5pfjet_pu_jet_simple_mva[n], ak5pfjet_pu_jet_full_loose[n], ak5pfjet_pu_jet_full_medium[n], ak5pfjet_pu_jet_full_tight[n], ak5pfjet_pu_jet_full_mva[n], ak5pfjet_sv_count[n]));
}

Track Analyse::Tracks(UInt_t n) const
{
  return(Track(sqrt(track_px[n]*track_px[n]+track_py[n]*track_py[n]+track_pz[n]*track_pz[n]), track_px[n], track_py[n], track_pz[n], track_outerx[n], track_outery[n], track_outerz[n], track_closestpointx[n], track_closestpointy[n], track_closestpointz[n], track_chi2[n], track_ndof[n], track_dxy[n], track_dxyerr[n], track_dz[n], track_dzerr[n], track_charge[n], track_nhits[n], track_nmissinghits[n], track_npixelhits[n], track_npixellayers[n], track_nstriplayers[n], track_dedxharmonic2[n]));
}

Vertex Analyse::PrimVertex() const
{
  return( Vertex(primvertex_x, primvertex_y, primvertex_z, primvertex_chi2, primvertex_ndof, primvertex_ptq, primvertex_ntracks, primvertex_cov));
}

SuperCluster Analyse::SuperClusters(UInt_t n, TVector3 refpoint) const
{
	SuperCluster newsupercluster(supercluster_e[n], supercluster_x[n] - refpoint.X(), supercluster_y[n] - refpoint.Y(), supercluster_z[n] - refpoint.Z(), supercluster_rawe[n], supercluster_phiwidth[n], supercluster_etawidth[n]);

if(loadsuperclusters == 2 || loadsuperclusters == 4)
{
	UInt_t endcluster = supercluster_basiccluster_count;
	if(n != supercluster_count) endcluster = supercluster_basicclusterbegin[n+1];
	for(UInt_t i = supercluster_basicclusterbegin[n] ; i < endcluster ; i++)
	{
		Cluster newcluster(supercluster_basiccluster_e[i], supercluster_basiccluster_x[i] - refpoint.X(), supercluster_basiccluster_y[i] - refpoint.Y(), supercluster_basiccluster_z[i] - refpoint.Z(), supercluster_basiccluster_nhit[i]);
		if(loadsuperclusters == 3 || loadsuperclusters == 4)
		{
			UInt_t endhits = supercluster_basiccluster_hit_count;
			if(i != supercluster_basiccluster_count) endhits = supercluster_basiccluster_hitbegin[i+1];
			for(UInt_t u = supercluster_basiccluster_hitbegin[i] ; u < endhits ; u++)
			{
				EcalHit newhit(supercluster_basiccluster_hit_e[u],supercluster_basiccluster_hit_x[u],supercluster_basiccluster_hit_y[u],supercluster_basiccluster_hit_z[u]);
				newcluster.AddHit(newhit);
			}

		}
		newsupercluster.AddCluster(newcluster);
	}
	if(false)
	{
		endcluster = supercluster_escluster_count;
		if(n != supercluster_count) endcluster = supercluster_esclusterbegin[n+1];
		for(UInt_t i = supercluster_esclusterbegin[n] ; i < endcluster ; i++)
		{
			Cluster newcluster(supercluster_escluster_e[i], supercluster_escluster_x[i] - refpoint.X(), supercluster_escluster_y[i] - refpoint.Y(), supercluster_escluster_z[i] - refpoint.Z(), supercluster_escluster_nhit[i]);
			if(loadsuperclusters == 3 || loadsuperclusters == 4)
			{
				UInt_t endhits = supercluster_escluster_hit_count;
				if(i != supercluster_escluster_count) endhits = supercluster_escluster_hitbegin[i+1];
				for(UInt_t u = supercluster_escluster_hitbegin[i] ; u < endhits ; u++)
				{
					EcalHit newhit(supercluster_escluster_hit_e[u],supercluster_escluster_hit_x[u],supercluster_escluster_hit_y[u],supercluster_escluster_hit_z[u]);
					newcluster.AddHit(newhit);
				}

			}
			newsupercluster.AddESCluster(newcluster);
		}
	}
}
	return(newsupercluster);
}

SuperCluster Analyse::SuperClusters(UInt_t n) const
{
  if(NumPrimVertices() > 0)
    {
      return(SuperClusters(n, PrimVertex()));
    }
  else
    {
      return(SuperClusters(n, TVector3(0., 0., 0.)));
    }
}

GenLightParticle Analyse::GenParticles(UInt_t n) const
{        
  return(GenLightParticle(genparticles_e[n], genparticles_px[n], genparticles_py[n], genparticles_pz[n], genparticles_vx[n], genparticles_vy[n], genparticles_vz[n], genparticles_status[n], genparticles_pdgid[n], genparticles_info[n]));
}

GenParticle Analyse::AllGenParticles(UInt_t n) const
{        
  UInt_t motherend(n == genallparticles_count-1 ? genallparticlesmother_count : genallparticles_motherbeg[n+1]);
  UInt_t daughterend(n == genallparticles_count-1 ? genallparticlesdaughter_count : genallparticles_daughterbeg[n+1]);
  LorentzVector TV_p4(genallparticles_px[n], genallparticles_py[n], genallparticles_pz[n], genallparticles_e[n]);
  GenParticle newpart(this, n, genallparticles_e[n], genallparticles_px[n], genallparticles_py[n], genallparticles_pz[n], 
		      genallparticles_vx[n], genallparticles_vy[n], genallparticles_vz[n], 
		      genallparticles_status[n], genallparticles_pdgid[n], genallparticles_motherbeg[n], 
		      motherend - genallparticles_motherbeg[n], 
		      genallparticles_daughterbeg[n], 
		      daughterend - genallparticles_daughterbeg[n], genallparticles_charge[n], TV_p4);
//   GenParticle newpart(this, n, genallparticles_e[n], genallparticles_px[n], genallparticles_py[n], genallparticles_pz[n], 
// 		      genallparticles_vx[n], genallparticles_vy[n], genallparticles_vz[n], 
// 		      genallparticles_status[n], genallparticles_pdgid[n], genallparticles_motherbeg[n], 
// 		      motherend - genallparticles_motherbeg[n], 
// 		      genallparticles_daughterbeg[n], 
// 		      genallparticles_daughters[n] - genallparticles_daughterbeg[n], genallparticles_charge[n], TV_p4);

  //genallparticles_daughters
  return(newpart);
}

TLorentzVector Analyse::GenParticleP4(UInt_t n) const
{
  return(TLorentzVector(genallparticles_px[n], genallparticles_py[n], genallparticles_pz[n], genallparticles_e[n]) );
}


TLorentzVector Analyse::GenMETTrue() const
{
  return(TLorentzVector(genmettrue_ex, genmettrue_ey, 0., sqrt(pow(genmettrue_ex, 2) + pow(genmettrue_ey, 2))));
}


//main loop
Long64_t Analyse::Loop(Long64_t start, Long64_t end)
{
  //  cout<<"Inside the Analyse::Loop "<< endl;
  if(end == -1 || end > GetNumAddedEvents())    end = GetNumAddedEvents();
  //cout<<"GetNumAddedEvents() : "<< GetNumAddedEvents() << endl;  

  //Pileup Reweighting
  if(usepileupinfo && useprimvertexinfo) cerr << "Don't use UsePileUpInfo() and UsePrimVertexInfo() at the same time." << endl;
  GetEvent(0, 1);
  
  if(NumPileUpInteractions() != -1) {
    cout<<"Passed the NumPileUpInteractions criteria"   <<endl;

    if(usepileupinfo) {
      cerr << "Reading PileUp-Information" << endl;
      for(unsigned i = 0 ; i < pileUpDist.size() ; i++) 	{
	pileUpDistMinus[i] = 0.;
	pileUpDist[i] = 0.;
	pileUpDistPlus[i] = 0.;
      }//for(unsigned i = 0 ; i < pileUpDist.size() ; i++)
      
      Double_t totnum = 0.;
      for(Long64_t i = start ; i < end ; i++)  {
	cout<<"Loop("<<i << ")"<<endl;
	GetEvent(i, 1);

	if(unsigned(NumPileUpInteractionsMinus()) < pileUpDistMinus.size()  && 
	   unsigned(NumPileUpInteractions())      < pileUpDist.size()       && 
	   unsigned(NumPileUpInteractionsPlus())  < pileUpDistPlus.size()) 	  {
	  
	  pileUpDistMinus[NumPileUpInteractionsMinus()]++;
	  pileUpDist[NumPileUpInteractions()]++;
	  pileUpDistPlus[NumPileUpInteractionsPlus()]++;
	  totnum++;
	}
      }// for(Long64_t i = start ; i < end ; i++)

      for(unsigned i = 0 ; i < pileUpDist.size() ; i++) {
	pileUpDistMinus[i] /= totnum;
	pileUpDist[i] /= totnum;
	pileUpDistPlus[i] /= totnum;
      }
    }//if(usepileupinfo)

    cout<<"Moving to Vertex Reweighting"<< endl;
    //Vertex Reweighting
    if(useprimvertexinfo)     {
      cerr << "Reading PrimVertex-Information" << endl;
      for(unsigned i = 0 ; i < primVertexDist.size() ; i++)  primVertexDist[i] = 0.;
      
      Double_t totnum = 0.;
      for(Long64_t i = start ; i < end ; i++) 	{
	GetEvent(i, 2);
	int numgoodprimvertices = NumGoodPrimVertices();
	totnum++;
	primVertexDist[numgoodprimvertices]++;
      }
      for(unsigned i = 0 ; i < primVertexDist.size() ; i++)  primVertexDist[i] /= totnum;
    }
  }//if(NumPileUpInteractions() != -1)
  else
    {
      useprimvertexinfo = false;
      usepileupinfo = false;
    }
  
  processed = 0;
  cerr << GetNumAddedEvents() << " Events are added." << endl;
  
  //Looping
  cerr << "Events " << start << " - " << end << " will be processed (" << end - start << " Events)." << endl;
  BeginLoop();
  
  TStopwatch watch;
  Int_t evres = 0; 
  double timeremain;
  watch.Start();




  //Main Loop over events
  for(Long64_t i = start ; i < end ; i++) {
    if(doDebug)    cout<<"   "<< endl;
    if(doDebug)    cout<<"@@@@@@@ Processing event ::"<< i << endl;
    GetEvent(i);

    if(errors != 0) {
      cerr << "ERROR " << errors << "in Event: Nr: " << UInt_t(Number()) << ", Run: " << Run() << ", LumiBlock: " << LumiBlock() << ". It's rejected!" << endl;
      continue;
    }
    
    if(jsonfilter) {
      if(jsonlist[Run()][LumiBlock()]) evres = AnalyzeEvent();
      else evres = 0;
    }
    else evres = AnalyzeEvent();
    processed++;
    
    if(duplicatecheck)	eventlist[Run()][LumiBlock()][Number()]++;
    

    if(lumicalculation && evres > 0) {
      ++lumilist[Run()][LumiBlock()];
      //   cout<<"inside the lumicalculation with event "<<evres<< endl;//" and lumiList"<< lumilist[Run()][LumiBlock()] <<endl;

      //Setting the Range
      if(Run() < minRun){minRun = Run(); minLumi = LumiBlock();}
      if(Run() == minRun && LumiBlock() < minLumi) {minLumi = LumiBlock();} 
      if(Run() > maxRun){maxRun = Run(); maxLumi = LumiBlock();}
      if(Run() == maxRun && LumiBlock() > maxLumi) {maxLumi = LumiBlock();} 
    } // if(lumicalculation && evres > 0)
    
    if(printinfo > 0 && processed%printinfo == 0) {
      //watch.Stop();
      timeremain = watch.RealTime()/processed*(end-start-processed);
      watch.Continue();
      cerr << "Processed: " << processed << "/" << end - start << ", Event: " << UInt_t(Number()) << ", Run: " << Run() << ", LumiBlock: " << LumiBlock() << ", Time: " << int(timeremain)/3600 << ":" << int(timeremain)%3600/60 << ":" << int(timeremain)%60 << endl;
    } // if(printinfo > 0 && processed%printinfo == 0)
    
    if(evres == -1) break;
  } //for(Long64_t i = start ; i < end ; i++)
  EndLoop();
  //  cout<<"Existing the Loop() with processed_event"<< processed << endl;
  return(processed);
}//Loop(Long64_t start, Long64_t end)

Int_t Analyse::PrepareSkimming(string filename)
{
	skimfilename = filename;
	if(skimfile != 0) {cerr << "YOU CAN USE PREPARESKIMMING() ONLY ONCE." << endl; return(-1);}
	
	skimfile = new TFile(skimfilename.c_str(), "recreate");
	skimfile->cd();
	skimtree = new TTree("AC1B", "AC1B", 1);
	skimtree->Branch("errors", &errors, "errors/i");
	skimtree->Branch("event_nr", &event_nr, "event_nr/i");
	skimtree->Branch("event_run", &event_run, "event_run/i");
	skimtree->Branch("event_timeunix", &event_timeunix, "event_timeunix/i");
	skimtree->Branch("event_timemicrosec", &event_timemicrosec, "event_timemicrosec/i");
	skimtree->Branch("event_luminosityblock", &event_luminosityblock, "event_luminosityblock/i");
	skimtree->Branch("trigger_level1bits", &trigger_level1bits, "trigger_level1bits[8]/b");
	skimtree->Branch("trigger_level1", &trigger_level1, "trigger_level1[128]/b");
	skimtree->Branch("trigger_HLT", &trigger_HLT, "trigger_HLT[128]/b");

	skimtree->Branch("beamspot_x", &beamspot_x, "beamspot_x/F");
	skimtree->Branch("beamspot_y", &beamspot_y, "beamspot_y/F");
	skimtree->Branch("beamspot_z", &beamspot_z, "beamspot_z/F");
	skimtree->Branch("beamspot_xwidth", &beamspot_xwidth, "beamspot_xwidth/F");
	skimtree->Branch("beamspot_ywidth", &beamspot_ywidth, "beamspot_ywidth/F");
	skimtree->Branch("beamspot_zsigma", &beamspot_zsigma, "beamspot_zsigma/F");
	skimtree->Branch("beamspot_cov", &beamspot_cov, "beamspot_cov[6]/F");

	skimtree->Branch("track_count", &track_count, "track_count/i"); 
	skimtree->Branch("track_px", track_px, "track_px[track_count]/F");
	skimtree->Branch("track_py", track_py, "track_py[track_count]/F");
	skimtree->Branch("track_pz", track_pz, "track_pz[track_count]/F");
	skimtree->Branch("track_outerx", track_outerx, "track_outerx[track_count]/F");
	skimtree->Branch("track_outery", track_outery, "track_outery[track_count]/F");
	skimtree->Branch("track_outerz", track_outerz, "track_outerz[track_count]/F");
	skimtree->Branch("track_closestpointx", track_closestpointx, "track_closestpointx[track_count]/F");
	skimtree->Branch("track_closestpointy", track_closestpointy, "track_closestpointy[track_count]/F");
	skimtree->Branch("track_closestpointz", track_closestpointz, "track_closestpointz[track_count]/F");
	skimtree->Branch("track_chi2", track_chi2, "track_chi2[track_count]/F");
	skimtree->Branch("track_ndof", track_ndof, "track_ndof[track_count]/F");
	skimtree->Branch("track_dxy", track_dxy, "track_dxy[track_count]/F");
	skimtree->Branch("track_dxyerr", track_dxyerr, "track_dxyerr[track_count]/F");
	skimtree->Branch("track_dz", track_dz, "track_dz[track_count]/F");
	skimtree->Branch("track_dzerr", track_dzerr, "track_dzerr[track_count]/F");
	skimtree->Branch("track_dedxharmonic2", track_dedxharmonic2, "track_dedxharmonic2[track_count]/F");
	skimtree->Branch("track_charge", track_charge, "track_charge[track_count]/I");
	skimtree->Branch("track_nhits", track_nhits, "track_nhits[track_count]/b");
	skimtree->Branch("track_npixelhits", track_npixelhits, "track_npixelhits[track_count]/b");
	skimtree->Branch("track_nmissinghits", track_nmissinghits, "track_nmissinghits[track_count]/b");
	skimtree->Branch("track_npixellayers", track_npixellayers, "track_npixellayers[track_count]/b");
	skimtree->Branch("track_nstriplayers", track_nstriplayers, "track_nstriplayers[track_count]/b");

	skimtree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i"); 
	skimtree->Branch("primvertex_x", &primvertex_x, "primvertex_x/F");
	skimtree->Branch("primvertex_y", &primvertex_y, "primvertex_y/F");
	skimtree->Branch("primvertex_z", &primvertex_z, "primvertex_z/F");
	skimtree->Branch("primvertex_chi2", &primvertex_chi2, "primvertex_chi2/F");
	skimtree->Branch("primvertex_ndof", &primvertex_ndof, "primvertex_ndof/F");
	skimtree->Branch("primvertex_ptq", &primvertex_ptq, "primvertex_ptq/F");
	skimtree->Branch("primvertex_ntracks", &primvertex_ntracks, "primvertex_ntracks/I");
	skimtree->Branch("primvertex_cov", primvertex_cov, "primvertex_cov[6]/F");

	skimtree->Branch("supercluster_count", &supercluster_count, "supercluster_count/i");
	skimtree->Branch("supercluster_e", supercluster_e, "supercluster_e[supercluster_count]/F");
	skimtree->Branch("supercluster_x", supercluster_x, "supercluster_x[supercluster_count]/F");
	skimtree->Branch("supercluster_y", supercluster_y, "supercluster_y[supercluster_count]/F");
	skimtree->Branch("supercluster_z", supercluster_z, "supercluster_z[supercluster_count]/F");
	skimtree->Branch("supercluster_rawe", supercluster_rawe, "supercluster_rawe[supercluster_count]/F");
	skimtree->Branch("supercluster_phiwidth", supercluster_phiwidth, "supercluster_phiwidth[supercluster_count]/F");
	skimtree->Branch("supercluster_etawidth", supercluster_etawidth, "supercluster_etawidth[supercluster_count]/F");
	skimtree->Branch("supercluster_nbasiccluster", supercluster_nbasiccluster, "supercluster_nbasiccluster[supercluster_count]/I");
	skimtree->Branch("supercluster_basicclusterbegin", supercluster_basicclusterbegin, "supercluster_basicclusterbegin[supercluster_count]/I");
	skimtree->Branch("supercluster_esclusterbegin", supercluster_esclusterbegin, "supercluster_esclusterbegin[supercluster_count]/I");

	skimtree->Branch("supercluster_basiccluster_count", &supercluster_basiccluster_count, "supercluster_basiccluster_count/i");
	skimtree->Branch("supercluster_basiccluster_e", supercluster_basiccluster_e, "supercluster_basiccluster_e[supercluster_basiccluster_count]/F");
	skimtree->Branch("supercluster_basiccluster_x", supercluster_basiccluster_x, "supercluster_basiccluster_x[supercluster_basiccluster_count]/F");
	skimtree->Branch("supercluster_basiccluster_y", supercluster_basiccluster_y, "supercluster_basiccluster_y[supercluster_basiccluster_count]/F");
	skimtree->Branch("supercluster_basiccluster_z", supercluster_basiccluster_z, "supercluster_basiccluster_z[supercluster_basiccluster_count]/F");
	skimtree->Branch("supercluster_basiccluster_nhit", supercluster_basiccluster_nhit, "supercluster_basiccluster_nhit[supercluster_basiccluster_count]/I");
	skimtree->Branch("supercluster_basiccluster_hitbegin", supercluster_basiccluster_hitbegin, "supercluster_basiccluster_hitbegin[supercluster_basiccluster_count]/I");

	skimtree->Branch("supercluster_basiccluster_hit_count", &supercluster_basiccluster_hit_count, "supercluster_basiccluster_hit_count/i");
	skimtree->Branch("supercluster_basiccluster_hit_e", supercluster_basiccluster_hit_e, "supercluster_basiccluster_hit_e[supercluster_basiccluster_hit_count]/F");
	skimtree->Branch("supercluster_basiccluster_hit_x", supercluster_basiccluster_hit_x, "supercluster_basiccluster_hit_x[supercluster_basiccluster_hit_count]/F");
	skimtree->Branch("supercluster_basiccluster_hit_y", supercluster_basiccluster_hit_y, "supercluster_basiccluster_hit_y[supercluster_basiccluster_hit_count]/F");
	skimtree->Branch("supercluster_basiccluster_hit_z", supercluster_basiccluster_hit_z, "supercluster_basiccluster_hit_z[supercluster_basiccluster_hit_count]/F");

	skimtree->Branch("supercluster_escluster_count", &supercluster_escluster_count, "supercluster_escluster_count/i");
	skimtree->Branch("supercluster_escluster_e", supercluster_escluster_e, "supercluster_escluster_e[supercluster_escluster_count]/F");
	skimtree->Branch("supercluster_escluster_x", supercluster_escluster_x, "supercluster_escluster_x[supercluster_escluster_count]/F");
	skimtree->Branch("supercluster_escluster_y", supercluster_escluster_y, "supercluster_escluster_y[supercluster_escluster_count]/F");
	skimtree->Branch("supercluster_escluster_z", supercluster_escluster_z, "supercluster_escluster_z[supercluster_escluster_count]/F");
	skimtree->Branch("supercluster_escluster_nhit", supercluster_escluster_nhit, "supercluster_escluster_nhit[supercluster_escluster_count]/I");
	skimtree->Branch("supercluster_escluster_hitbegin", supercluster_escluster_hitbegin, "supercluster_escluster_hitbegin[supercluster_escluster_count]/I");

	skimtree->Branch("supercluster_escluster_hit_count", &supercluster_escluster_hit_count, "supercluster_escluster_hit_count/i");
	skimtree->Branch("supercluster_escluster_hit_e", supercluster_escluster_hit_e, "supercluster_escluster_hit_e[supercluster_escluster_hit_count]/F");
	skimtree->Branch("supercluster_escluster_hit_x", supercluster_escluster_hit_x, "supercluster_escluster_hit_x[supercluster_escluster_hit_count]/F");
	skimtree->Branch("supercluster_escluster_hit_y", supercluster_escluster_hit_y, "supercluster_escluster_hit_y[supercluster_escluster_hit_count]/F");
	skimtree->Branch("supercluster_escluster_hit_z", supercluster_escluster_hit_z, "supercluster_escluster_hit_z[supercluster_escluster_hit_count]/F");

	skimtree->Branch("muon_count", &muon_count, "muon_count/i");
	skimtree->Branch("muon_px", muon_px, "muon_px[muon_count]/F");
	skimtree->Branch("muon_py", muon_py, "muon_py[muon_count]/F");
	skimtree->Branch("muon_pz", muon_pz, "muon_pz[muon_count]/F");
	skimtree->Branch("muon_pterror", muon_pterror, "muon_pterror[muon_count]/F");
	skimtree->Branch("muon_chi2", muon_chi2, "muon_chi2[muon_count]/F");
	skimtree->Branch("muon_ndof", muon_ndof, "muon_ndof[muon_count]/F");
	skimtree->Branch("muon_nchambers", muon_nchambers, "muon_nchambers[muon_count]/i");
	skimtree->Branch("muon_innertrack_px", muon_innertrack_px, "muon_innertrack_px[muon_count]/F");
	skimtree->Branch("muon_innertrack_py", muon_innertrack_py, "muon_innertrack_py[muon_count]/F");
	skimtree->Branch("muon_innertrack_pz", muon_innertrack_pz, "muon_innertrack_pz[muon_count]/F");
	skimtree->Branch("muon_innertrack_outerx", muon_innertrack_outerx, "muon_innertrack_outerx[muon_count]/F");
	skimtree->Branch("muon_innertrack_outery", muon_innertrack_outery, "muon_innertrack_outery[muon_count]/F");
	skimtree->Branch("muon_innertrack_outerz", muon_innertrack_outerz, "muon_innertrack_outerz[muon_count]/F");
	skimtree->Branch("muon_innertrack_closestpointx", muon_innertrack_closestpointx, "muon_innertrack_closestpointx[muon_count]/F");
	skimtree->Branch("muon_innertrack_closestpointy", muon_innertrack_closestpointy, "muon_innertrack_closestpointy[muon_count]/F");
	skimtree->Branch("muon_innertrack_closestpointz", muon_innertrack_closestpointz, "muon_innertrack_closestpointz[muon_count]/F");
	skimtree->Branch("muon_innertrack_chi2", muon_innertrack_chi2, "muon_innertrack_chi2[muon_count]/F");
	skimtree->Branch("muon_innertrack_ndof", muon_innertrack_ndof, "muon_innertrack_ndof[muon_count]/F");
	skimtree->Branch("muon_innertrack_dxy", muon_innertrack_dxy, "muon_innertrack_dxy[muon_count]/F");
	skimtree->Branch("muon_innertrack_dxyerr", muon_innertrack_dxyerr, "muon_innertrack_dxyerr[muon_count]/F");
	skimtree->Branch("muon_innertrack_dz", muon_innertrack_dz, "muon_innertrack_dz[muon_count]/F");
	skimtree->Branch("muon_innertrack_dzerr", muon_innertrack_dzerr, "muon_innertrack_dzerr[muon_count]/F");
	skimtree->Branch("muon_innertrack_dedxharmonic2", muon_innertrack_dedxharmonic2, "muon_innertrack_dedxharmonic2[muon_count]/F");
	skimtree->Branch("muon_innertrack_charge", muon_innertrack_charge, "muon_innertrack_charge[muon_count]/I");
	skimtree->Branch("muon_innertrack_nhits", muon_innertrack_nhits, "muon_innertrack_nhits[muon_count]/b");
	skimtree->Branch("muon_innertrack_npixelhits", muon_innertrack_npixelhits, "muon_innertrack_npixelhits[muon_count]/b");
	skimtree->Branch("muon_innertrack_nmissinghits", muon_innertrack_nmissinghits, "muon_innertrack_nmissinghits[muon_count]/b");
	skimtree->Branch("muon_innertrack_npixellayers", muon_innertrack_npixellayers, "muon_innertrack_npixellayers[muon_count]/b");
	skimtree->Branch("muon_innertrack_nstriplayers", muon_innertrack_nstriplayers, "muon_innertrack_nstriplayers[muon_count]/b");
	skimtree->Branch("muon_outertrack_px", muon_outertrack_px, "muon_outertrack_px[muon_count]/F");
	skimtree->Branch("muon_outertrack_py", muon_outertrack_py, "muon_outertrack_py[muon_count]/F");
	skimtree->Branch("muon_outertrack_pz", muon_outertrack_pz, "muon_outertrack_pz[muon_count]/F");
	skimtree->Branch("muon_outertrack_hits", muon_outertrack_hits, "muon_outertrack_hits[muon_count]/b");
	skimtree->Branch("muon_outertrack_missinghits", muon_outertrack_missinghits, "muon_outertrack_missinghits[muon_count]/b");
	skimtree->Branch("muon_outertrack_chi2", muon_outertrack_chi2, "muon_outertrack_chi2[muon_count]/F");
	skimtree->Branch("muon_outertrack_ndof", muon_outertrack_ndof, "muon_outertrack_ndof[muon_count]/F");
	skimtree->Branch("muon_isolationr3track", muon_isolationr3track, "muon_isolationr3track[muon_count]/F");
	skimtree->Branch("muon_isolationr3ntrack", muon_isolationr3ntrack, "muon_isolationr3ntrack[muon_count]/I");
	skimtree->Branch("muon_isolationr3ecal", muon_isolationr3ecal, "muon_isolationr3ecal[muon_count]/F");
	skimtree->Branch("muon_isolationr3hcal", muon_isolationr3hcal, "muon_isolationr3hcal[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr3_sumchargedhadronpt", muon_pfisolationr3_sumchargedhadronpt, "muon_pfisolationr3_sumchargedhadronpt[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr3_sumchargedparticlept", muon_pfisolationr3_sumchargedparticlept, "muon_pfisolationr3_sumchargedparticlept[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr3_sumneutralhadronet", muon_pfisolationr3_sumneutralhadronet, "muon_pfisolationr3_sumneutralhadronet[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr3_sumphotonet", muon_pfisolationr3_sumphotonet, "muon_pfisolationr3_sumphotonet[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr3_sumPUpt", muon_pfisolationr3_sumPUpt, "muon_pfisolationr3_sumPUpt[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr4_sumchargedhadronpt", muon_pfisolationr4_sumchargedhadronpt, "muon_pfisolationr4_sumchargedhadronpt[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr4_sumchargedparticlept", muon_pfisolationr4_sumchargedparticlept, "muon_pfisolationr4_sumchargedparticlept[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr4_sumneutralhadronet", muon_pfisolationr4_sumneutralhadronet, "muon_pfisolationr4_sumneutralhadronet[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr4_sumphotonet", muon_pfisolationr4_sumphotonet, "muon_pfisolationr4_sumphotonet[muon_count]/F");
  	skimtree->Branch("muon_pfisolationr4_sumPUpt", muon_pfisolationr4_sumPUpt, "muon_pfisolationr4_sumPUpt[muon_count]/F");
	skimtree->Branch("muon_ecalenergy", muon_ecalenergy, "muon_ecalenergy[muon_count]/F");
	skimtree->Branch("muon_hcalenergy", muon_hcalenergy, "muon_hcalenergy[muon_count]/F");
	skimtree->Branch("muon_charge", muon_charge, "muon_charge[muon_count]/I");
        skimtree->Branch("muon_numchambers", muon_numchambers, "muon_numchambers[muon_count]/I");
        skimtree->Branch("muon_numchamberswithsegments", muon_numchamberswithsegments, "muon_numchamberswithsegments[muon_count]/I");
	skimtree->Branch("muon_type", muon_type, "muon_type[muon_count]/i");
	skimtree->Branch("muon_trigger", muon_trigger, "muon_trigger[muon_count]/i");
	skimtree->Branch("muon_trackermuonquality", muon_trackermuonquality, "muon_trackermuonquality[muon_count]/i");
	skimtree->Branch("muon_mva_id", muon_mva_id, "muon_mva_id[muon_count]/F");
	skimtree->Branch("muon_mva_iso", muon_mva_iso, "muon_mva_iso[muon_count]/F");

	skimtree->Branch("ak5calojet_count", &ak5calojet_count, "ak5calojet_count/i");
	skimtree->Branch("ak5calojet_e", ak5calojet_e, "ak5calojet_e[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_px", ak5calojet_px, "ak5calojet_px[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_py", ak5calojet_py, "ak5calojet_py[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_pz", ak5calojet_pz, "ak5calojet_pz[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_hadronicenergy", ak5calojet_hadronicenergy, "ak5calojet_hadronicenergy[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_emenergy", ak5calojet_emenergy, "ak5calojet_emenergy[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_energycorr", ak5calojet_energycorr, "ak5calojet_energycorr[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_fhpd", ak5calojet_fhpd, "ak5calojet_fhpd[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_restrictedemf", ak5calojet_restrictedemf, "ak5calojet_restrictedemf[ak5calojet_count]/F");
	skimtree->Branch("ak5calojet_btag", ak5calojet_btag, "ak5calojet_btag[ak5calojet_count][6]/F");
	skimtree->Branch("ak5calojet_n90", ak5calojet_n90, "ak5calojet_n90[ak5calojet_count]/i");
	skimtree->Branch("ak5calojet_n60", ak5calojet_n60, "ak5calojet_n60[ak5calojet_count]/i");

	skimtree->Branch("ak5jptjet_count", &ak5jptjet_count, "ak5jptjet_count/i");
	skimtree->Branch("ak5jptjet_e", ak5jptjet_e, "ak5jptjet_e[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_px", ak5jptjet_px, "ak5jptjet_px[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_py", ak5jptjet_py, "ak5jptjet_py[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_pz", ak5jptjet_pz, "ak5jptjet_pz[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_hadronicenergy", ak5jptjet_hadronicenergy, "ak5jptjet_hadronicenergy[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_chargedhadronicenergy", ak5jptjet_chargedhadronicenergy, "ak5jptjet_chargedhadronicenergy[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_emenergy", ak5jptjet_emenergy, "ak5jptjet_emenergy[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_chargedemenergy", ak5jptjet_chargedemenergy, "ak5jptjet_chargedemenergy[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_chargedmulti", ak5jptjet_chargedmulti, "ak5jptjet_chargedmulti[ak5jptjet_count]/i");	
	skimtree->Branch("ak5jptjet_energycorr", ak5jptjet_energycorr, "ak5jptjet_energycorr[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_fhpd", ak5jptjet_fhpd, "ak5jptjet_fhpd[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_restrictedemf", ak5jptjet_restrictedemf, "ak5jptjet_restrictedemf[ak5jptjet_count]/F");
	skimtree->Branch("ak5jptjet_btag", ak5jptjet_btag, "ak5jptjet_btag[ak5jptjet_count][6]/F");
	skimtree->Branch("ak5jptjet_n90", ak5jptjet_n90, "ak5jptjet_n90[ak5jptjet_count]/i");

	skimtree->Branch("ak5pfjet_count", &ak5pfjet_count, "ak5pfjet_count/i");
	skimtree->Branch("ak5pfjet_e", ak5pfjet_e, "ak5pfjet_e[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_px", ak5pfjet_px, "ak5pfjet_px[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_py", ak5pfjet_py, "ak5pfjet_py[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_pz", ak5pfjet_pz, "ak5pfjet_pz[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_hadronicenergy", ak5pfjet_hadronicenergy, "ak5pfjet_hadronicenergy[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_chargedhadronicenergy", ak5pfjet_chargedhadronicenergy, "ak5pfjet_chargedhadronicenergy[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_emenergy", ak5pfjet_emenergy, "ak5pfjet_emenergy[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_chargedemenergy", ak5pfjet_chargedemenergy, "ak5pfjet_chargedemenergy[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_chargedmulti", ak5pfjet_chargedmulti, "ak5pfjet_chargedmulti[ak5pfjet_count]/i");	
	skimtree->Branch("ak5pfjet_neutralmulti", ak5pfjet_neutralmulti, "ak5pfjet_neutralmulti[ak5pfjet_count]/i");	
	skimtree->Branch("ak5pfjet_energycorr", ak5pfjet_energycorr, "ak5pfjet_energycorr[ak5pfjet_count]/F");
	skimtree->Branch("ak5pfjet_btag", ak5pfjet_btag, "ak5pfjet_btag[ak5pfjet_count][6]/F");
        skimtree->Branch("ak5pfjet_trigger", ak5pfjet_trigger, "ak5pfjet_trigger[ak5pfjet_count]/i");
        skimtree->Branch("ak5pfjet_sv_count", ak5pfjet_sv_count, "ak5pfjet_sv_count[ak5pfjet_count]/i");
        skimtree->Branch("ak5pfjet_sv_x", ak5pfjet_sv_x, "ak5pfjet_sv_x[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_sv_y", ak5pfjet_sv_y, "ak5pfjet_sv_y[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_sv_z", ak5pfjet_sv_z, "ak5pfjet_sv_z[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_sv_dx", ak5pfjet_sv_dx, "ak5pfjet_sv_dx[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_sv_dy", ak5pfjet_sv_dy, "ak5pfjet_sv_dy[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_sv_dz", ak5pfjet_sv_dz, "ak5pfjet_sv_dz[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_pu_jet_simple_loose", ak5pfjet_pu_jet_simple_loose, "ak5pfjet_pu_jet_simple_loose[ak5pfjet_count]/O");
	skimtree->Branch("ak5pfjet_pu_jet_simple_medium", ak5pfjet_pu_jet_simple_medium, "ak5pfjet_pu_jet_simple_medium[ak5pfjet_count]/O");
	skimtree->Branch("ak5pfjet_pu_jet_simple_tight", ak5pfjet_pu_jet_simple_tight, "ak5pfjet_pu_jet_simple_tight[ak5pfjet_count]/O");
	skimtree->Branch("ak5pfjet_pu_jet_simple_mva", ak5pfjet_pu_jet_simple_mva, "ak5pfjet_pu_jet_simple_mva[ak5pfjet_count]/F");
        skimtree->Branch("ak5pfjet_pu_jet_full_loose", ak5pfjet_pu_jet_full_loose, "ak5pfjet_pu_jet_full_loose[ak5pfjet_count]/O");
	skimtree->Branch("ak5pfjet_pu_jet_full_medium", ak5pfjet_pu_jet_full_medium, "ak5pfjet_pu_jet_full_medium[ak5pfjet_count]/O");
	skimtree->Branch("ak5pfjet_pu_jet_full_tight", ak5pfjet_pu_jet_full_tight, "ak5pfjet_pu_jet_full_tight[ak5pfjet_count]/O");
	skimtree->Branch("ak5pfjet_pu_jet_full_mva", ak5pfjet_pu_jet_full_mva, "ak5pfjet_pu_jet_full_mva[ak5pfjet_count]/F");

	skimtree->Branch("electron_count", &electron_count, "electron_count/i");
	skimtree->Branch("electron_px", electron_px, "electron_px[electron_count]/F");
	skimtree->Branch("electron_py", electron_py, "electron_py[electron_count]/F");
	skimtree->Branch("electron_pz", electron_pz, "electron_pz[electron_count]/F");
	skimtree->Branch("electron_trackchi2", electron_trackchi2, "electron_trackchi2[electron_count]/F");
	skimtree->Branch("electron_trackndof", electron_trackndof, "electron_trackndof[electron_count]/F");
	skimtree->Branch("electron_outerx", electron_outerx, "electron_outerx[electron_count]/F");
	skimtree->Branch("electron_outery", electron_outery, "electron_outery[electron_count]/F");
	skimtree->Branch("electron_outerz", electron_outerz, "electron_outerz[electron_count]/F");
	skimtree->Branch("electron_closestpointx", electron_closestpointx, "electron_closestpointx[electron_count]/F");
	skimtree->Branch("electron_closestpointy", electron_closestpointy, "electron_closestpointy[electron_count]/F");
	skimtree->Branch("electron_closestpointz", electron_closestpointz, "electron_closestpointz[electron_count]/F");
	skimtree->Branch("electron_esuperclusterovertrack", electron_esuperclusterovertrack, "electron_esuperclusterovertrack[electron_count]/F");
	skimtree->Branch("electron_eseedclusterovertrack", electron_eseedclusterovertrack, "electron_eseedclusterovertrack[electron_count]/F");
	skimtree->Branch("electron_deltaetasuperclustertrack", electron_deltaetasuperclustertrack, "electron_deltaetasuperclustertrack[electron_count]/F");
	skimtree->Branch("electron_deltaphisuperclustertrack", electron_deltaphisuperclustertrack, "electron_deltaphisuperclustertrack[electron_count]/F");
	skimtree->Branch("electron_e1x5", electron_e1x5, "electron_e1x5[electron_count]/F");
	skimtree->Branch("electron_e2x5", electron_e2x5, "electron_e2x5[electron_count]/F");
	skimtree->Branch("electron_e5x5", electron_e5x5, "electron_e5x5[electron_count]/F");
	skimtree->Branch("electron_sigmaetaeta", electron_sigmaetaeta, "electron_sigmaetaeta[electron_count]/F");
	skimtree->Branch("electron_sigmaietaieta", electron_sigmaietaieta, "electron_sigmaietaieta[electron_count]/F");
	skimtree->Branch("electron_ehcaloverecal", electron_ehcaloverecal, "electron_ehcaloverecal[electron_count]/F");
	skimtree->Branch("electron_ehcaloverecaldepth1", electron_ehcaloverecaldepth1, "electron_ehcaloverecaldepth1[electron_count]/F");
	skimtree->Branch("electron_ehcaloverecaldepth2", electron_ehcaloverecaldepth2, "electron_ehcaloverecaldepth2[electron_count]/F");
	skimtree->Branch("electron_isolationr3track", electron_isolationr3track, "electron_isolationr3track[electron_count]/F");
	skimtree->Branch("electron_isolationr3ecal", electron_isolationr3ecal, "electron_isolationr3ecal[electron_count]/F");
	skimtree->Branch("electron_isolationr3hcal", electron_isolationr3hcal, "electron_isolationr3hcal[electron_count]/F");
	skimtree->Branch("electron_isolationr4track", electron_isolationr4track, "electron_isolationr4track[electron_count]/F");
	skimtree->Branch("electron_isolationr4ecal", electron_isolationr4ecal, "electron_isolationr4ecal[electron_count]/F");
	skimtree->Branch("electron_isolationr4hcal", electron_isolationr4hcal, "electron_isolationr4hcal[electron_count]/F");
	skimtree->Branch("electron_pfisolationr3_sumchargedhadronpt", electron_pfisolationr3_sumchargedhadronpt, "electron_pfisolationr3_sumchargedhadronpt[electron_count]/F");
    skimtree->Branch("electron_pfisolationr3_sumchargedparticlept", electron_pfisolationr3_sumchargedparticlept, "electron_pfisolationr3_sumchargedparticlept[electron_count]/F");
    skimtree->Branch("electron_pfisolationr3_sumneutralhadronet", electron_pfisolationr3_sumneutralhadronet, "electron_pfisolationr3_sumneutralhadronet[electron_count]/F");
    skimtree->Branch("electron_pfisolationr3_sumphotonet", electron_pfisolationr3_sumphotonet, "electron_pfisolationr3_sumphotonet[electron_count]/F");
    skimtree->Branch("electron_pfisolationr3_sumPUpt", electron_pfisolationr3_sumPUpt, "electron_pfisolationr3_sumPUpt[electron_count]/F");
    skimtree->Branch("electron_pfisolationr4_sumchargedhadronpt", electron_pfisolationr4_sumchargedhadronpt, "electron_pfisolationr4_sumchargedhadronpt[electron_count]/F");
    skimtree->Branch("electron_pfisolationr4_sumchargedparticlept", electron_pfisolationr4_sumchargedparticlept, "electron_pfisolationr4_sumchargedparticlept[electron_count]/F");
    skimtree->Branch("electron_pfisolationr4_sumneutralhadronet", electron_pfisolationr4_sumneutralhadronet, "electron_pfisolationr4_sumneutralhadronet[electron_count]/F");
    skimtree->Branch("electron_pfisolationr4_sumphotonet", electron_pfisolationr4_sumphotonet, "electron_pfisolationr4_sumphotonet[electron_count]/F");
    skimtree->Branch("electron_pfisolationr4_sumPUpt", electron_pfisolationr4_sumPUpt, "electron_pfisolationr4_sumPUpt[electron_count]/F");
	skimtree->Branch("electron_nhits", electron_nhits, "electron_nhits[electron_count]/b");
	skimtree->Branch("electron_npixelhits", electron_npixelhits, "electron_npixelhits[electron_count]/b"); 
	skimtree->Branch("electron_nmissinghits", electron_nmissinghits, "electron_nmissinghits[electron_count]/b");
	skimtree->Branch("electron_npixellayers", electron_npixellayers, "electron_npixellayers[electron_count]/b");
	skimtree->Branch("electron_nstriplayers", electron_nstriplayers, "electron_nstriplayers[electron_count]/b");
	skimtree->Branch("electron_dxy", electron_dxy, "electron_dxy[electron_count]/F");
	skimtree->Branch("electron_dxyerr", electron_dxyerr, "electron_dxyerr[electron_count]/F");
	skimtree->Branch("electron_dz", electron_dz, "electron_dz[electron_count]/F");
	skimtree->Branch("electron_dzerr", electron_dzerr, "electron_dzerr[electron_count]/F"); 
	skimtree->Branch("electron_convdist", electron_convdist, "electron_convdist[electron_count]/F");
	skimtree->Branch("electron_convdcot", electron_convdcot, "electron_convdcot[electron_count]/F");
	skimtree->Branch("electron_convradius", electron_convradius, "electron_convradius[electron_count]/F");
	skimtree->Branch("electron_gapinfo", electron_gapinfo, "electron_gapinfo[electron_count]/i");
	skimtree->Branch("electron_chargeinfo", electron_chargeinfo, "electron_chargeinfo[electron_count]/i");
	skimtree->Branch("electron_fbrems", electron_fbrems, "electron_fbrems[electron_count]/F");
	skimtree->Branch("electron_numbrems", electron_numbrems, "electron_numbrems[electron_count]/I");
	skimtree->Branch("electron_charge", electron_charge, "electron_charge[electron_count]/I");
	skimtree->Branch("electron_superclusterindex", electron_superclusterindex, "electron_superclusterindex[electron_count]/I");
	skimtree->Branch("electron_info", electron_info, "electron_info[electron_count]/b");
	skimtree->Branch("electron_trigger", electron_trigger, "electron_trigger[electron_count]/i");
	skimtree->Branch("electron_mva_id_trig", electron_mva_id_trig, "electron_mva_id_trig[electron_count]/F");
	skimtree->Branch("electron_mva_id_nontrig", electron_mva_id_nontrig, "electron_mva_id_nontrig[electron_count]/F");
	skimtree->Branch("electron_mva_iso", electron_mva_iso, "electron_mva_iso[electron_count]/F");
	skimtree->Branch("electron_has_conversion", electron_has_conversion, "electron_has_conversion[electron_count]/O");

	skimtree->Branch("photon_count", &photon_count, "photon_count/i");
	skimtree->Branch("photon_px", photon_px, "photon_px[photon_count]/F");
	skimtree->Branch("photon_py", photon_py, "photon_py[photon_count]/F");
	skimtree->Branch("photon_pz", photon_pz, "photon_pz[photon_count]/F");
	skimtree->Branch("photon_e1x5", photon_e1x5, "photon_e1x5[photon_count]/F");
	skimtree->Branch("photon_e2x5", photon_e2x5, "photon_e2x5[photon_count]/F");
	skimtree->Branch("photon_e3x3", photon_e3x3, "photon_e3x3[photon_count]/F");
	skimtree->Branch("photon_e5x5", photon_e5x5, "photon_e5x5[photon_count]/F");
	skimtree->Branch("photon_sigmaetaeta", photon_sigmaetaeta, "photon_sigmaetaeta[photon_count]/F");
	skimtree->Branch("photon_sigmaietaieta", photon_sigmaietaieta, "photon_sigmaietaieta[photon_count]/F");
	skimtree->Branch("photon_ehcaloverecal", photon_ehcaloverecal, "photon_ehcaloverecal[photon_count]/F");
	skimtree->Branch("photon_ehcaloverecaldepth1", photon_ehcaloverecaldepth1, "photon_ehcaloverecaldepth1[photon_count]/F");
	skimtree->Branch("photon_ehcaloverecaldepth2", photon_ehcaloverecaldepth2, "photon_ehcaloverecaldepth2[photon_count]/F");
	skimtree->Branch("photon_maxenergyxtal", photon_maxenergyxtal, "photon_maxenergyxtal[photon_count]/F");
	skimtree->Branch("photon_isolationr3track", photon_isolationr3track, "photon_isolationr3track[photon_count]/F");
	skimtree->Branch("photon_isolationr3trackhollow", photon_isolationr3trackhollow, "photon_isolationr3trackhollow[photon_count]/F");
	skimtree->Branch("photon_isolationr3ntrack", photon_isolationr3ntrack, "photon_isolationr3ntrack[photon_count]/i");
	skimtree->Branch("photon_isolationr3ntrackhollow", photon_isolationr3ntrackhollow, "photon_isolationr3ntrackhollow[photon_count]/i");
	skimtree->Branch("photon_isolationr3ecal", photon_isolationr3ecal, "photon_isolationr3ecal[photon_count]/F");
	skimtree->Branch("photon_isolationr3hcal", photon_isolationr3hcal, "photon_isolationr3hcal[photon_count]/F");
	skimtree->Branch("photon_isolationr4track", photon_isolationr4track, "photon_isolationr4track[photon_count]/F");
	skimtree->Branch("photon_isolationr4trackhollow", photon_isolationr4trackhollow, "photon_isolationr4trackhollow[photon_count]/F");
	skimtree->Branch("photon_isolationr4ntrack", photon_isolationr4ntrack, "photon_isolationr4ntrack[photon_count]/i");
	skimtree->Branch("photon_isolationr4ntrackhollow", photon_isolationr4ntrackhollow, "photon_isolationr4ntrackhollow[photon_count]/i");
	skimtree->Branch("photon_isolationr4ecal", photon_isolationr4ecal, "photon_isolationr4ecal[photon_count]/F");
	skimtree->Branch("photon_isolationr4hcal", photon_isolationr4hcal, "photon_isolationr4hcal[photon_count]/F");
	skimtree->Branch("photon_superclusterindex", photon_superclusterindex, "photon_superclusterindex[photon_count]/I");
	skimtree->Branch("photon_info", photon_info, "photon_info[photon_count]/b");
	skimtree->Branch("photon_gapinfo", photon_gapinfo, "photon_gapinfo[photon_count]/i");
	skimtree->Branch("photon_conversionbegin", photon_conversionbegin, "photon_conversionbegin[photon_count]/i");

	skimtree->Branch("conversion_count", &conversion_count, "conversion_count/i");
	skimtree->Branch("conversion_info", conversion_info, "conversion_info[conversion_count]/b");
	skimtree->Branch("conversion_vx", conversion_vx, "conversion_vx[conversion_count]/F");
	skimtree->Branch("conversion_vy", conversion_vy, "conversion_vy[conversion_count]/F");
	skimtree->Branch("conversion_vz", conversion_vz, "conversion_vz[conversion_count]/F");
	skimtree->Branch("conversion_chi2", conversion_chi2, "conversion_chi2[conversion_count]/F");
	skimtree->Branch("conversion_ndof", conversion_ndof, "conversion_ndof[conversion_count]/F");
	skimtree->Branch("conversion_cov", conversion_cov, "conversion_cov[conversion_count][6]/F");
	skimtree->Branch("conversion_mvaout", conversion_mvaout, "conversion_mvaout[conversion_count]/F");
	skimtree->Branch("conversion_trackecalpointx", conversion_trackecalpointx, "conversion_trackecalpointx[conversion_count][2]/F");
	skimtree->Branch("conversion_trackecalpointy", conversion_trackecalpointy, "conversion_trackecalpointy[conversion_count][2]/F");
	skimtree->Branch("conversion_trackecalpointz", conversion_trackecalpointz, "conversion_trackecalpointz[conversion_count][2]/F");
	skimtree->Branch("conversion_trackpx", conversion_trackpx, "conversion_trackpx[conversion_count][2]/F");
	skimtree->Branch("conversion_trackpy", conversion_trackpy, "conversion_trackpy[conversion_count][2]/F");
	skimtree->Branch("conversion_trackpz", conversion_trackpz, "conversion_trackpz[conversion_count][2]/F");
	skimtree->Branch("conversion_trackclosestpointx", conversion_trackclosestpointx, "conversion_trackclosestpointx[conversion_count][2]/F");
	skimtree->Branch("conversion_trackclosestpointy", conversion_trackclosestpointy, "conversion_trackclosestpointy[conversion_count][2]/F");
	skimtree->Branch("conversion_trackclosestpointz", conversion_trackclosestpointz, "conversion_trackclosestpointz[conversion_count][2]/F");
	skimtree->Branch("conversion_trackchi2", conversion_trackchi2, "conversion_trackchi2[conversion_count][2]/F");
	skimtree->Branch("conversion_trackndof", conversion_trackndof, "conversion_trackndof[conversion_count][2]/F");
	skimtree->Branch("conversion_trackdxy", conversion_trackdxy, "conversion_trackdxy[conversion_count][2]/F");
	skimtree->Branch("conversion_trackdxyerr", conversion_trackdxyerr, "conversion_trackdxyerr[conversion_count][2]/F");
	skimtree->Branch("conversion_trackdz", conversion_trackdz, "conversion_trackdz[conversion_count][2]/F");
	skimtree->Branch("conversion_trackdzerr", conversion_trackdzerr, "conversion_trackdzerr[conversion_count][2]/F");
	skimtree->Branch("conversion_trackcharge", conversion_trackcharge, "conversion_trackcharge[conversion_count][2]/I");
	skimtree->Branch("conversion_tracknhits", conversion_tracknhits, "conversion_tracknhits[conversion_count][2]/b");
	skimtree->Branch("conversion_tracknmissinghits", conversion_tracknmissinghits, "conversion_tracknmissinghits[conversion_count][2]/b");
	skimtree->Branch("conversion_tracknpixelhits", conversion_tracknpixelhits, "conversion_tracknpixelhits[conversion_count][2]/b");
	skimtree->Branch("conversion_tracknpixellayers", conversion_tracknpixellayers, "conversion_tracknpixellayers[conversion_count][2]/b");
	skimtree->Branch("conversion_tracknstriplayers", conversion_tracknstriplayers, "conversion_tracknstriplayers[conversion_count][2]/b");
	//tau
  skimtree->Branch("tau_count", &tau_count, "tau_count/i");
  skimtree->Branch("tau_e", tau_e, "tau_e[tau_count]/F");
  skimtree->Branch("tau_px", tau_px, "tau_px[tau_count]/F");
  skimtree->Branch("tau_py", tau_py, "tau_py[tau_count]/F");
  skimtree->Branch("tau_pz", tau_pz, "tau_pz[tau_count]/F");
  skimtree->Branch("tau_isolationneutralspt", tau_isolationneutralspt, "tau_isolationneutralspt[tau_count]/F");
  skimtree->Branch("tau_isolationneutralsnum", tau_isolationneutralsnum, "tau_isolationneutralsnum[tau_count]/i");
  skimtree->Branch("tau_isolationchargedpt", tau_isolationchargedpt, "tau_isolationchargedpt[tau_count]/F");
  skimtree->Branch("tau_isolationchargednum", tau_isolationchargednum, "tau_isolationchargednum[tau_count]/i");
  skimtree->Branch("tau_isolationgammapt", tau_isolationgammapt, "tau_isolationgammapt[tau_count]/F");
  skimtree->Branch("tau_isolationgammanum", tau_isolationgammanum, "tau_isolationgammanum[tau_count]/i");
  skimtree->Branch("tau_leadpfchargedhadrcandecalenergy", tau_leadpfchargedhadrcandecalenergy, "tau_leadpfchargedhadrcandecalenergy[tau_count]/F");
  skimtree->Branch("tau_leadpfchargedhadrcandhcalenergy", tau_leadpfchargedhadrcandhcalenergy, "tau_leadpfchargedhadrcandhcalenergy[tau_count]/F");
  skimtree->Branch("tau_leadpfchargedhadrcandp", tau_leadpfchargedhadrcandp, "tau_leadpfchargedhadrcandp[tau_count]/F");
  skimtree->Branch("tau_leadpfchargedhadrcandpt", tau_leadpfchargedhadrcandpt, "tau_leadpfchargedhadrcandpt[tau_count]/F");
  //  skimtree->Branch("tau_leadpfchargedhadrcandp4", tau_leadpfchargedhadrcandp4, "tau_leadpfchargedhadrcandp4[tau_count]/F");
  skimtree->Branch("tau_dxy", tau_dxy, "tau_dxy[tau_count]/F");
  skimtree->Branch("tau_dz", tau_dz, "tau_dz[tau_count]/F");
  skimtree->Branch("tau_ip3d", tau_ip3d, "tau_ip3d[tau_count]/F");
  skimtree->Branch("tau_ip3dSig", tau_ip3dSig, "tau_ip3dSig[tau_count]/F");
  skimtree->Branch("tau_nprongs", tau_nprongs, "tau_nprongs[tau_count]/i");
  skimtree->Branch("tau_charge", tau_charge, "tau_charge[tau_count]/I");
  skimtree->Branch("tau_emfraction", tau_emfraction, "tau_emfraction[tau_count]/F");
  skimtree->Branch("tau_newemfraction", tau_newemfraction, "tau_newemfraction[tau_count]/F");
  skimtree->Branch("tau_hcaltotoverplead", tau_hcaltotoverplead, "tau_hcaltotoverplead[tau_count]/F");
  skimtree->Branch("tau_hcal3x3overplead", tau_hcal3x3overplead, "tau_hcal3x3overplead[tau_count]/F");
  skimtree->Branch("tau_ecalstripsumeoverplead", tau_ecalstripsumeoverplead, "tau_ecalstripsumeoverplead[tau_count]/F");
  skimtree->Branch("tau_bremsrecoveryeoverplead", tau_bremsrecoveryeoverplead, "tau_bremsrecoveryeoverplead[tau_count]/F");
  skimtree->Branch("tau_calocomp", tau_calocomp, "tau_calocomp[tau_count]/F");
  skimtree->Branch("tau_segcomp", tau_segcomp, "tau_segcomp[tau_count]/F");
  //  skimtree->Branch("tau_disbyisolationmvaraw", tau_disbyisolationmvaraw, "tau_disbyisolationmvaraw[tau_count]/F");
  //  skimtree->Branch("tau_disbyisolationmva2raw", tau_disbyisolationmva2raw,"tau_disbyisolationmva2raw[tau_count]/F");
  //  skimtree->Branch("tau_againstelectronmva2raw", tau_againstelectronmva2raw, "tau_againstelectronmva2raw[tau_count]/F");
  //  skimtree->Branch("tau_againstelectronmva2category", tau_againstelectronmva2category, "tau_againstelectronmva2category[tau_count]/F");
  //  skimtree->Branch("tau_againstelectronmva3raw", tau_againstelectronmva3raw, "tau_againstelectronmva3raw[tau_count]/F");
  //  skimtree->Branch("tau_againstelectronmva3category", tau_againstelectronmva3category, "tau_againstelectronmva3category[tau_count]/F");
  skimtree->Branch("tau_bycombinedisolationdeltabetacorrraw3hits", tau_bycombinedisolationdeltabetacorrraw3hits, "tau_bycombinedisolationdeltabetacorrraw3hits[tau_count]/F");
  skimtree->Branch("tau_dishps", tau_dishps, "tau_dishps[tau_count]/l");
  skimtree->Branch("tau_trigger", tau_trigger, "tau_trigger[tau_count]/i");
  skimtree->Branch("tau_L1trigger_match", tau_L1trigger_match, "tau_L1trigger_match[tau_count]/O");

  //ditau
  skimtree->Branch("ditau_Index", &ditau_Index, "ditau_Index/i");

  // ditaupair
  skimtree->Branch("mutautaupair_count", &mutautaupair_count, "mutautaupair_count/i");
  skimtree->Branch("mutautaupair_leg1_px", mutautaupair_leg1_px, "mutautaupair_leg1_px[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg1_py", mutautaupair_leg1_py, "mutautaupair_leg1_py[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg1_pz", mutautaupair_leg1_pz, "mutautaupair_leg1_pz[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg1_energy", mutautaupair_leg1_energy, "mutautaupair_leg1_energy[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg2_px", mutautaupair_leg2_px, "mutautaupair_leg2_px[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg2_py", mutautaupair_leg2_py, "mutautaupair_leg2_py[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg2_pz", mutautaupair_leg2_pz, "mutautaupair_leg2_pz[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_leg2_energy", mutautaupair_leg2_energy, "mutautaupair_leg2_energy[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_mu_px", mutautaupair_mu_px, "mutautaupair_mu_px[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_mu_py", mutautaupair_mu_py, "mutautaupair_mu_py[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_mu_pz", mutautaupair_mu_pz, "mutautaupair_mu_pz[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_mu_energy", mutautaupair_mu_energy, "mutautaupair_mu_energy[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_svfit_int_valid", mutautaupair_svfit_int_valid, "mutautaupair_svfit_int_valid[mutautaupair_count]/O");
  skimtree->Branch("mutautaupair_svfit_mass_int", mutautaupair_svfit_mass_int, "mutautaupair_svfit_mass_int[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_svfit_mass_int_err_up", mutautaupair_svfit_mass_int_err_up, "mutautaupair_svfit_mass_int_err_up[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_svfit_mass_int_err_down", mutautaupair_svfit_mass_int_err_down, "mutautaupair_svfit_mass_int_err_down[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_dca2d", mutautaupair_dca2d, "mutautaupair_dca2d[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_dca2d_err", mutautaupair_dca2d_err, "mutautaupair_dca2d_err[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_dca3d", mutautaupair_dca3d, "mutautaupair_dca3d[mutautaupair_count]/F");
  skimtree->Branch("mutautaupair_dca3d_err", mutautaupair_dca3d_err, "mutautaupair_dca3d_err[mutautaupair_count]/F");

  skimtree->Branch("eltautaupair_count", &eltautaupair_count, "eltautaupair_count/i");
  skimtree->Branch("eltautaupair_leg1_px", eltautaupair_leg1_px, "eltautaupair_leg1_px[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg1_py", eltautaupair_leg1_py, "eltautaupair_leg1_py[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg1_pz", eltautaupair_leg1_pz, "eltautaupair_leg1_pz[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg1_energy", eltautaupair_leg1_energy, "eltautaupair_leg1_energy[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg2_px", eltautaupair_leg2_px, "eltautaupair_leg2_px[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg2_py", eltautaupair_leg2_py, "eltautaupair_leg2_py[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg2_pz", eltautaupair_leg2_pz, "eltautaupair_leg2_pz[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_leg2_energy", eltautaupair_leg2_energy, "eltautaupair_leg2_energy[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_el_px", eltautaupair_el_px, "eltautaupair_el_px[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_el_py", eltautaupair_el_py, "eltautaupair_el_py[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_el_pz", eltautaupair_el_pz, "eltautaupair_el_pz[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_el_energy", eltautaupair_el_energy, "eltautaupair_el_energy[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_svfit_int_valid", eltautaupair_svfit_int_valid, "eltautaupair_svfit_int_valid[eltautaupair_count]/O");
  skimtree->Branch("eltautaupair_svfit_mass_int", eltautaupair_svfit_mass_int, "eltautaupair_svfit_mass_int[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_svfit_mass_int_err_up", eltautaupair_svfit_mass_int_err_up, "eltautaupair_svfit_mass_int_err_up[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_svfit_mass_int_err_down", eltautaupair_svfit_mass_int_err_down, "eltautaupair_svfit_mass_int_err_down[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_dca2d", eltautaupair_dca2d, "eltautaupair_dca2d[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_dca2d_err", eltautaupair_dca2d_err, "eltautaupair_dca2d_err[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_dca3d", eltautaupair_dca3d, "eltautaupair_dca3d[eltautaupair_count]/F");
  skimtree->Branch("eltautaupair_dca3d_err", eltautaupair_dca3d_err, "eltautaupair_dca3d_err[eltautaupair_count]/F");

  skimtree->Branch("tau_ak5pfjet_e", tau_ak5pfjet_e, "tau_ak5pfjet_e[tau_count]/F");
  skimtree->Branch("tau_ak5pfjet_px", tau_ak5pfjet_px, "tau_ak5pfjet_px[tau_count]/F");
  skimtree->Branch("tau_ak5pfjet_py", tau_ak5pfjet_py, "tau_ak5pfjet_py[tau_count]/F");
  skimtree->Branch("tau_ak5pfjet_pz", tau_ak5pfjet_pz, "tau_ak5pfjet_pz[tau_count]/F");

  skimtree->Branch("tau_chargedbegin", tau_chargedbegin, "tau_chargedbegin[tau_count]/i");
  skimtree->Branch("tau_charged_count", &tau_charged_count, "tau_charged_count/i");
  skimtree->Branch("tau_charged_px", tau_charged_px, "tau_charged_px[tau_charged_count]/F");
  skimtree->Branch("tau_charged_py", tau_charged_py, "tau_charged_py[tau_charged_count]/F");
  skimtree->Branch("tau_charged_pz", tau_charged_pz, "tau_charged_pz[tau_charged_count]/F");
  skimtree->Branch("tau_charged_outerx", tau_charged_outerx, "tau_charged_outerx[tau_charged_count]/F");
  skimtree->Branch("tau_charged_outery", tau_charged_outery, "tau_charged_outery[tau_charged_count]/F");
  skimtree->Branch("tau_charged_outerz", tau_charged_outerz, "tau_charged_outerz[tau_charged_count]/F");
  skimtree->Branch("tau_charged_closestpointx", tau_charged_closestpointx, "tau_charged_closestpointx[tau_charged_count]/F");
  skimtree->Branch("tau_charged_closestpointy", tau_charged_closestpointy, "tau_charged_closestpointy[tau_charged_count]/F");
  skimtree->Branch("tau_charged_closestpointz", tau_charged_closestpointz, "tau_charged_closestpointz[tau_charged_count]/F");
  skimtree->Branch("tau_charged_chi2", tau_charged_chi2, "tau_charged_chi2[tau_charged_count]/F");
  skimtree->Branch("tau_charged_ndof", tau_charged_ndof, "tau_charged_ndof[tau_charged_count]/F");
  skimtree->Branch("tau_charged_dxy", tau_charged_dxy, "tau_charged_dxy[tau_charged_count]/F");
  skimtree->Branch("tau_charged_dxyerr", tau_charged_dxyerr, "tau_charged_dxyerr[tau_charged_count]/F");
  skimtree->Branch("tau_charged_dz", tau_charged_dz, "tau_charged_dz[tau_charged_count]/F");
  skimtree->Branch("tau_charged_dzerr", tau_charged_dzerr, "tau_charged_dzerr[tau_charged_count]/F");
  skimtree->Branch("tau_charged_dedxharmonic2", tau_charged_dedxharmonic2, "tau_charged_dedxharmonic2[tau_charged_count]/F");
  skimtree->Branch("tau_charged_charge", tau_charged_charge, "tau_charged_charge[tau_charged_count]/I");
  skimtree->Branch("tau_charged_nhits", tau_charged_nhits, "tau_charged_nhits[tau_charged_count]/b");
  skimtree->Branch("tau_charged_nmissinghits", tau_charged_nmissinghits, "tau_charged_nmissinghits[tau_charged_count]/b");
  skimtree->Branch("tau_charged_npixelhits", tau_charged_npixelhits, "tau_charged_npixelhits[tau_charged_count]/b");
  skimtree->Branch("tau_charged_npixellayers", tau_charged_npixellayers, "tau_charged_npixellayers[tau_charged_count]/b");
  skimtree->Branch("tau_charged_nstriplayers", tau_charged_nstriplayers, "tau_charged_nstriplayers[tau_charged_count]/b");


 
        skimtree->Branch("ak5pfjet_rho", &ak5pfjet_rho, "ak5pfjet_rho/F");
        skimtree->Branch("ak5pfjet_sigma", &ak5pfjet_sigma, "ak5pfjet_sigma/F");

	skimtree->Branch("calomet_ex", &calomet_ex, "calomet_ex/F");
	skimtree->Branch("calomet_ey", &calomet_ey, "calomet_ey/F");
	skimtree->Branch("calomet_exmuons", &calomet_exmuons, "calomet_exmuons/F");
	skimtree->Branch("calomet_eymuons", &calomet_eymuons, "calomet_eymuons/F");

	skimtree->Branch("tcmet_ex", &tcmet_ex, "tcmet_ex/F");
	skimtree->Branch("tcmet_ey", &tcmet_ey, "tcmet_ey/F");

	skimtree->Branch("pfmet_ex", &pfmet_ex, "pfmet_ex/F");
	skimtree->Branch("pfmet_ey", &pfmet_ey, "pfmet_ey/F");
	skimtree->Branch("pfmettype1_ex", &pfmettype1_ex, "pfmettype1_ex/F");
	skimtree->Branch("pfmettype1_ey", &pfmettype1_ey, "pfmettype1_ey/F");
	skimtree->Branch("pfmetmva_ex", &pfmetmva_ex, "pfmetmva_ex/F");
	skimtree->Branch("pfmetmva_ey", &pfmetmva_ey, "pfmetmva_ey/F");

	skimtree->Branch("secvertices_count", &secvertices_count, "secvertices_count/i");
	skimtree->Branch("secvertices_vx", secvertices_vx, "secvertices_vx[secvertices_count]/F");
	skimtree->Branch("secvertices_vy", secvertices_vy, "secvertices_vy[secvertices_count]/F");
	skimtree->Branch("secvertices_vz", secvertices_vz, "secvertices_vz[secvertices_count]/F");
	skimtree->Branch("secvertices_chi2", secvertices_chi2, "secvertices_chi2[secvertices_count]/F");
	skimtree->Branch("secvertices_ndof", secvertices_ndof, "secvertices_ndof[secvertices_count]/F");
	skimtree->Branch("secvertices_cov", secvertices_cov, "secvertices_cov[secvertices_count][6]/F");
	skimtree->Branch("secvertices_track_px", secvertices_track_px, "secvertices_track_px[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_py", secvertices_track_py, "secvertices_track_py[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_pz", secvertices_track_pz, "secvertices_track_pz[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_closestpointx", secvertices_track_closestpointx, "secvertices_track_closestpointx[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_closestpointy", secvertices_track_closestpointy, "secvertices_track_closestpointy[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_closestpointz", secvertices_track_closestpointz, "secvertices_track_closestpointz[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_chi2", secvertices_track_chi2, "secvertices_track_chi2[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_ndof", secvertices_track_ndof, "secvertices_track_ndof[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_dxy", secvertices_track_dxy, "secvertices_track_dxy[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_dxyerr", secvertices_track_dxyerr, "secvertices_track_dxyerr[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_dz", secvertices_track_dz, "secvertices_track_dz[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_dzerr", secvertices_track_dzerr, "secvertices_track_dzerr[secvertices_count][2]/F");
	skimtree->Branch("secvertices_track_dedxharmonic2", secvertices_track_dedxharmonic2, "secvertices_track_dedxharmonic2[track_count][2]/F");
	skimtree->Branch("secvertices_track_charge", secvertices_track_charge, "secvertices_track_charge[secvertices_count][2]/I");
	skimtree->Branch("secvertices_track_nhits", secvertices_track_nhits, "secvertices_track_nhits[secvertices_count][2]/b");
	skimtree->Branch("secvertices_track_nmissinghits", secvertices_track_nmissinghits, "secvertices_track_nmissinghits[secvertices_count][2]/b");
	skimtree->Branch("secvertices_track_npixelhits", secvertices_track_npixelhits, "secvertices_track_npixelhits[secvertices_count][2]/b");
	skimtree->Branch("secvertices_track_npixellayers", secvertices_track_npixellayers, "secvertices_track_npixellayers[secvertices_count][2]/b");
	skimtree->Branch("secvertices_track_nstriplayers", secvertices_track_nstriplayers, "secvertices_track_nstriplayers[secvertices_count][2]/b");

	skimtree->Branch("genweight", &genweight, "genweight/F");
	skimtree->Branch("genid1", &genid1, "genid1/F");
	skimtree->Branch("genx1", &genx1, "genx1/F");
	skimtree->Branch("genid2", &genid2, "genid2/F");
	skimtree->Branch("genx2", &genx2, "genx2/F");
	skimtree->Branch("genScale", &genScale, "genScale/F");

	skimtree->Branch("numpileupinteractionsminus", &numpileupinteractionsminus, "numpileupinteractionsminus/I");
	skimtree->Branch("numpileupinteractions", &numpileupinteractions, "numpileupinteractions/I");
	skimtree->Branch("numpileupinteractionsplus", &numpileupinteractionsplus, "numpileupinteractionsplus/I");
	skimtree->Branch("numtruepileupinteractions", &numtruepileupinteractions, "numtruepileupinteractions/F");

	skimtree->Branch("genparticles_count", &genparticles_count, "genparticles_count/i");
	skimtree->Branch("genparticles_e", genparticles_e, "genparticles_e[genparticles_count]/F");
	skimtree->Branch("genparticles_px", genparticles_px, "genparticles_px[genparticles_count]/F");
	skimtree->Branch("genparticles_py", genparticles_py, "genparticles_py[genparticles_count]/F");
	skimtree->Branch("genparticles_pz", genparticles_pz, "genparticles_pz[genparticles_count]/F");
	skimtree->Branch("genparticles_vx", genparticles_vx, "genparticles_vx[genparticles_count]/F");
	skimtree->Branch("genparticles_vy", genparticles_vy, "genparticles_vy[genparticles_count]/F");
	skimtree->Branch("genparticles_vz", genparticles_vz, "genparticles_vz[genparticles_count]/F");
	skimtree->Branch("genparticles_pdgid", genparticles_pdgid, "genparticles_pdgid[genparticles_count]/I");
	skimtree->Branch("genparticles_status", genparticles_status, "genparticles_status[genparticles_count]/I");
	skimtree->Branch("genparticles_info", genparticles_info, "genparticles_info[genparticles_count]/i");

	skimtree->Branch("genallparticles_count", &genallparticles_count, "genallparticles_count/i");
	skimtree->Branch("genallparticles_e", genallparticles_e, "genallparticles_e[genallparticles_count]/F");
	skimtree->Branch("genallparticles_px", genallparticles_px, "genallparticles_px[genallparticles_count]/F");
	skimtree->Branch("genallparticles_py", genallparticles_py, "genallparticles_py[genallparticles_count]/F");
	skimtree->Branch("genallparticles_pz", genallparticles_pz, "genallparticles_pz[genallparticles_count]/F");
	skimtree->Branch("genallparticles_vx", genallparticles_vx, "genallparticles_vx[genallparticles_count]/F");
	skimtree->Branch("genallparticles_vy", genallparticles_vy, "genallparticles_vy[genallparticles_count]/F");
	skimtree->Branch("genallparticles_vz", genallparticles_vz, "genallparticles_vz[genallparticles_count]/F");
	skimtree->Branch("genallparticles_pdgid", genallparticles_pdgid, "genallparticles_pdgid[genallparticles_count]/I");
	skimtree->Branch("genallparticles_charge", genallparticles_charge, "genallparticles_charge[genparticles_count]/I");
	skimtree->Branch("genallparticles_status", genallparticles_status, "genallparticles_status[genallparticles_count]/I");
	skimtree->Branch("genallparticles_motherbeg", genallparticles_motherbeg, "genallparticles_motherbeg[genallparticles_count]/i");
	skimtree->Branch("genallparticles_daughterbeg", genallparticles_daughterbeg, "genallparticles_daughterbeg[genallparticles_count]/i");

	skimtree->Branch("genallparticlesmother_count", &genallparticlesmother_count, "genallparticlesmother_count/i");
	skimtree->Branch("genallparticles_mothers", genallparticles_mothers, "genallparticles_mothers[genallparticlesmother_count]/i");

	skimtree->Branch("genallparticlesdaughter_count", &genallparticlesdaughter_count, "genallparticlesdaughter_count/i");
	skimtree->Branch("genallparticles_daughters", genallparticles_daughters, "genallparticles_daughters[genallparticlesdaughter_count]/i");

	skimtree->Branch("genmetcalo_ex", &genmetcalo_ex, "genmetcalo_ex/F");
	skimtree->Branch("genmetcalo_ey", &genmetcalo_ey, "genmetcalo_ey/F");
	skimtree->Branch("genmettrue_ex", &genmettrue_ex, "genmettrue_ex/F");
	skimtree->Branch("genmettrue_ey", &genmettrue_ey, "genmettrue_ey/F");

	return(1);
}

Int_t Analyse::SkimEvent()
{
  if(skimtree == 0)
    {
      cerr << "use PrepareSkimming() before SkimEvent()" << endl;
      return(-1);
    }
  selected[Run()][LumiBlock()]++;
  skimtree->Fill();
  return(1);
}

Int_t Analyse::SetLumi(UInt_t run, UInt_t block, Float_t lumival)
{
  if(lumilist.find(run) != lumilist.end() && lumilist.find(run)->second.find(block) != lumilist.find(run)->second.end())
    {
      lumilist[run][block].Value(lumival);
      return(1);
    }
  return(0);
}

void Analyse::WriteLumiFile(string filename)
{
	TFile* file = new TFile(filename.c_str(), "recreate");

	UInt_t run_number;
	UInt_t run_hltcount;
	Char_t run_hltnames[20000];
	Char_t run_hltmunames[10000];
	Char_t run_hltelnames[10000];
	Char_t run_hlttaunames[10000];
	Char_t run_hltphotonnames[10000];
	Char_t run_hltjetnames[10000];
	Char_t run_taudiscriminators[10000];
	UInt_t run_hltprescaletablescount;
	UInt_t run_hltprescaletables[10000];
	UInt_t run_l1algocount;
	UInt_t run_l1algoprescaletablescount;
	UInt_t run_l1algoprescaletables[10000];
	UInt_t run_l1techcount;
	UInt_t run_l1techprescaletablescount;
	UInt_t run_l1techprescaletables[10000];
	
	TTree* runtree = new TTree("AC1Brun", "AC1Brun", 1);
	runtree->Branch("run_number", &run_number, "run_number/i");
	runtree->Branch("run_hltcount", &run_hltcount, "run_hltcount/i");
	runtree->Branch("run_hltnames", run_hltnames, "run_hltnames[run_hltcount]/C");
	runtree->Branch("run_hltmunames", run_hltmunames, "run_hltmunames/C");
	runtree->Branch("run_hltelnames", run_hltelnames, "run_hltelnames/C");
	runtree->Branch("run_hlttaunames", run_hlttaunames, "run_hlttaunames/C");
	runtree->Branch("run_hltphotonnames", run_hltphotonnames, "run_hltphotonnames/C");
	runtree->Branch("run_hltjetnames", run_hltjetnames, "run_hltjetnames/C");
	runtree->Branch("run_taudiscriminators", run_taudiscriminators, "run_taudiscriminators/C");
	runtree->Branch("run_hltprescaletablescount", &run_hltprescaletablescount, "run_hltprescaletablescount/i");
	runtree->Branch("run_hltprescaletables", run_hltprescaletables, "run_hltprescaletables[run_hltprescaletablescount]/i");
	runtree->Branch("run_l1algoprescaletablescount", &run_l1algoprescaletablescount, "run_l1algoprescaletablescount/i");
	runtree->Branch("run_l1algocount", &run_l1algocount, "run_l1algocount/i");
	runtree->Branch("run_l1algoprescaletables", run_l1algoprescaletables, "run_l1algoprescaletables[run_l1algoprescaletablescount]/i");
	runtree->Branch("run_l1techprescaletablescount", &run_l1techprescaletablescount, "run_l1techprescaletablescount/i");
	runtree->Branch("run_l1techcount", &run_l1techcount, "run_l1techcount/i");
	runtree->Branch("run_l1techprescaletables", run_l1techprescaletables, "run_l1techprescaletables[run_l1techprescaletablescount]/i");

	for(map<UInt_t, RunInfo>::iterator b = runlist.begin() ; b != runlist.end() ; ++b)
	{
		run_number = b->second.Run();
		run_hltcount = b->second.NumHLT();
		strcpy(run_hltnames, b->second.HLTAllNames().c_str());
		strcpy(run_hltmunames, b->second.HLTMuonAllNames().c_str());
		strcpy(run_hltelnames, b->second.HLTElectronAllNames().c_str());
		strcpy(run_hlttaunames, b->second.HLTTauAllNames().c_str());
		strcpy(run_hltphotonnames, b->second.HLTPhotonAllNames().c_str());
		strcpy(run_hltjetnames, b->second.HLTJetAllNames().c_str());
		strcpy(run_taudiscriminators, b->second.TauDiscriminatorsAllNames().c_str());
		run_hltprescaletablescount = b->second.NumHLT()*b->second.NumHLTTables();
		for(UInt_t i = 0 ; i < b->second.NumHLTTables() ; i++)
		{
			for(UInt_t j = 0 ; j < b->second.NumHLT() ; j++)
			{
				run_hltprescaletables[j+b->second.NumHLT()*i] = b->second.HLTPrescale(j,i);
			}
		}

		run_l1algocount = b->second.NumL1Algo();
		run_l1algoprescaletablescount = b->second.NumL1Algo()*b->second.NumL1AlgoTables();
		for(UInt_t i = 0 ; i < b->second.NumL1AlgoTables() ; i++)
		{
			for(UInt_t j = 0 ; j < b->second.NumL1Algo() ; j++)
			{
				run_l1algoprescaletables[j+b->second.NumL1Algo()*i] = b->second.L1AlgoPrescale(j,i);
			}
		}

		run_l1techcount = b->second.NumL1Tech();
		run_l1techprescaletablescount = b->second.NumL1Tech()*b->second.NumL1TechTables();
		for(UInt_t i = 0 ; i < b->second.NumL1TechTables() ; i++)
		{
			for(UInt_t j = 0 ; j < b->second.NumL1Tech() ; j++)
			{
				run_l1techprescaletables[j+b->second.NumL1Tech()*i] = b->second.L1TechPrescale(j,i);
			}
		}
		runtree->Fill();
	}
	runtree->Write();

	UInt_t lumi_run;
	UInt_t lumi_block;
	Float_t lumi_value;
	Float_t lumi_valueerr;
	Float_t lumi_livefrac;
	Float_t lumi_deadfrac;
	UInt_t lumi_quality;
	UInt_t lumi_eventsprocessed;
	UInt_t lumi_eventsfiltered;
	UInt_t lumi_hltprescaletable;
	UInt_t lumi_l1algoprescaletable;
	UInt_t lumi_l1techprescaletable;

	TTree* lumitree = new TTree("AC1Blumi", "AC1Blumi", 1);
	lumitree->Branch("lumi_run", &lumi_run, "lumi_run/i");
	lumitree->Branch("lumi_block", &lumi_block, "lumi_block/i");
	lumitree->Branch("lumi_value", &lumi_value, "lumi_value/F");
	lumitree->Branch("lumi_valueerr", &lumi_valueerr, "lumi_valueerr/F");
	lumitree->Branch("lumi_livefrac", &lumi_livefrac, "lumi_livefrac/F");
	lumitree->Branch("lumi_deadfrac", &lumi_deadfrac, "lumi_deadfrac/F");
	lumitree->Branch("lumi_quality", &lumi_quality, "lumi_quality/i");
	lumitree->Branch("lumi_eventsprocessed", &lumi_eventsprocessed, "lumi_eventsprocessed/i");
	lumitree->Branch("lumi_eventsfiltered", &lumi_eventsfiltered, "lumi_eventsfiltered/i");
	lumitree->Branch("lumi_hltprescaletable", &lumi_hltprescaletable, "lumi_hltprescaletable/i");
	lumitree->Branch("lumi_l1algoprescaletable", &lumi_l1algoprescaletable, "lumi_l1algoprescaletable/i");
	lumitree->Branch("lumi_l1techprescaletable", &lumi_l1techprescaletable, "lumi_l1techprescaletable/i");

	for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a)
	{
		for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; ++b)
		{
			lumi_run = a->first;
			lumi_block = b->first;
			lumi_value = b->second.LumiValue();
			lumi_valueerr = b->second.LumiError();
			lumi_livefrac = b->second.LiveFraction();
			lumi_deadfrac = b->second.DeadFraction();
			lumi_quality = b->second.Quality();
			lumi_eventsprocessed = b->second.NumEventsOrig();

			if(skimtree != 0) lumi_eventsfiltered = selected[a->first][b->first];
			else      	  lumi_eventsfiltered = b->second.NumEvents();
			
			lumi_hltprescaletable = b->second.HLTTable();
			lumi_l1algoprescaletable = b->second.L1AlgoTable();
			lumi_l1techprescaletable = b->second.L1TechTable();
			if(jsonfilter && jsonlist[lumi_run][lumi_block]) lumitree->Fill();
			if(!jsonfilter) lumitree->Fill();
		}
	}
	lumitree->Write();
	
	file->Write();
	file->Close();
}

void Analyse::AddLumiFile(string filename, string dir)
{
  TFile* lumifile = new TFile(filename.c_str());
  if(lumifile->IsZombie()) cerr << "ERROR AddLumiFile: " << filename << " is not a valid file." << endl;

  TTree* lumitree = 0;
  lumifile->GetObject((dir + string("AC1Blumi")).c_str(), lumitree);
  if(lumitree == 0)      cerr << "ERROR AddLumiFile: " << filename << " " << dir + string("AC1Blumi") << " is no TTree." << endl;

  TTree* runtree = 0;
  lumifile->GetObject((dir + string("AC1Brun")).c_str(), runtree);
  if(runtree == 0)      cerr << "ERROR AddLumiFile: " << filename << " " << dir + string("AC1Brun") << " is no TTree." << endl;

  UInt_t lumi_run;
  UInt_t lumi_block;
  Float_t lumi_value;
  Float_t lumi_valueerr;
  Float_t lumi_livefrac;
  Float_t lumi_deadfrac;
  UInt_t lumi_quality;
  UInt_t lumi_eventsprocessed;
  UInt_t lumi_eventsfiltered;
  UInt_t lumi_hltprescaletable;
  UInt_t lumi_l1algoprescaletable;
  UInt_t lumi_l1techprescaletable;
  
  lumitree->SetBranchAddress("lumi_run", &lumi_run);
  lumitree->SetBranchAddress("lumi_block", &lumi_block);
  lumitree->SetBranchAddress("lumi_value", &lumi_value);
  lumitree->SetBranchAddress("lumi_valueerr", &lumi_valueerr);
  lumitree->SetBranchAddress("lumi_livefrac", &lumi_livefrac);
  lumitree->SetBranchAddress("lumi_deadfrac", &lumi_deadfrac);
  lumitree->SetBranchAddress("lumi_quality", &lumi_quality);
  lumitree->SetBranchAddress("lumi_eventsprocessed", &lumi_eventsprocessed);
  lumitree->SetBranchAddress("lumi_eventsfiltered", &lumi_eventsfiltered);
  lumitree->SetBranchAddress("lumi_hltprescaletable", &lumi_hltprescaletable);
  lumitree->SetBranchAddress("lumi_l1algoprescaletable", &lumi_l1algoprescaletable);
  lumitree->SetBranchAddress("lumi_l1techprescaletable", &lumi_l1techprescaletable);
  
  for(UInt_t i = 0 ; i < lumitree->GetEntries() ; i++) {
    lumitree->GetEntry(i);
    if(lumilist.find(lumi_run)                          ==  lumilist.end() || 
       lumilist.find(lumi_run)->second.find(lumi_block) ==  lumilist.find(lumi_run)->second.end())
      {
	lumilist[lumi_run][lumi_block] = Luminosity(lumi_value, lumi_valueerr, lumi_livefrac, lumi_deadfrac, lumi_quality, lumi_eventsfiltered, lumi_eventsprocessed, lumi_l1algoprescaletable, lumi_l1algoprescaletable, lumi_l1techprescaletable);
      }
    else
      {
	Luminosity newlumi(lumi_value, lumi_valueerr, lumi_livefrac, lumi_deadfrac, lumi_quality, lumi_eventsfiltered, lumi_eventsprocessed, lumi_l1algoprescaletable, lumi_l1algoprescaletable, lumi_l1techprescaletable);
	lumilist[lumi_run][lumi_block] += newlumi;
	
	cout << "Warning: Information of run: " << lumi_run << ", block: " << lumi_block << " has already been filled. For data this means there are duplicated files for MC it might be ok. Check for duplicated events anyway!" << endl;
      }
  }
  
  UInt_t run_number;
  UInt_t run_hltcount;
  Char_t run_hltnames[20000];
  Char_t run_hltmunames[10000];
  Char_t run_hltelnames[10000];
  Char_t run_hlttaunames[10000];
  Char_t run_hltphotonnames[10000];
  Char_t run_hltjetnames[10000];
  Char_t run_taudiscriminators[10000];
  UInt_t run_hltprescaletablescount;
  UInt_t run_hltprescaletables[10000];
  UInt_t run_l1algocount;
  UInt_t run_l1algoprescaletablescount;
  UInt_t run_l1algoprescaletables[10000];
  UInt_t run_l1techcount;
  UInt_t run_l1techprescaletablescount;
  UInt_t run_l1techprescaletables[10000];
  
  runtree->SetBranchAddress("run_number", &run_number);
  runtree->SetBranchAddress("run_hltcount", &run_hltcount);
  runtree->SetBranchAddress("run_hltnames", run_hltnames);
  runtree->SetBranchAddress("run_hltmunames", run_hltmunames);
  runtree->SetBranchAddress("run_hltelnames", run_hltelnames);
  runtree->SetBranchAddress("run_hlttaunames", run_hlttaunames);
  runtree->SetBranchAddress("run_hltphotonnames", run_hltphotonnames);
  runtree->SetBranchAddress("run_hltjetnames", run_hltjetnames);
  runtree->SetBranchAddress("run_taudiscriminators", run_taudiscriminators);
  runtree->SetBranchAddress("run_hltprescaletablescount", &run_hltprescaletablescount);
  runtree->SetBranchAddress("run_hltprescaletables", run_hltprescaletables);
  runtree->SetBranchAddress("run_l1algocount", &run_l1algocount);
  runtree->SetBranchAddress("run_l1algoprescaletablescount", &run_l1algoprescaletablescount);
  runtree->SetBranchAddress("run_l1algoprescaletables", run_l1algoprescaletables);
  runtree->SetBranchAddress("run_l1techcount", &run_l1techcount);
  runtree->SetBranchAddress("run_l1techprescaletablescount", &run_l1techprescaletablescount);
  runtree->SetBranchAddress("run_l1techprescaletables", run_l1techprescaletables);
  
  for(UInt_t i = 0 ; i < runtree->GetEntries() ; i++)
    {
      runtree->GetEntry(i);
      if(runlist.find(run_number) == runlist.end())
	{
	  runlist[run_number] = RunInfo(run_number, run_hltcount, run_hltnames, run_hltmunames, run_hltelnames, run_hlttaunames, run_hltphotonnames, run_hltjetnames, run_taudiscriminators, run_hltprescaletablescount, run_hltprescaletables, run_l1algocount, run_l1algoprescaletablescount, run_l1algoprescaletables, run_l1techcount, run_l1techprescaletablescount, run_l1techprescaletables);
	}
    }
  
  lumicalculation = true;
  lumifile->Close();
}

Int_t Analyse::IsLumiAvailable() const
{
  if(lumilist.find(Run()) != lumilist.end() && 
     lumilist.find(Run())->second.find(LumiBlock()) != lumilist.find(Run())->second.end())
    {
      if(lumilist.at(Run()).at(LumiBlock()).LumiValue() != -1)
	{
	  return(2);
	}
      return(1);
    }
  return(0);
}

Int_t Analyse::GetNumHLTriggers() const
{
	if(runlist.find(Run()) != runlist.end())
	{
		return(runlist.at(Run()).NumHLT());
	}
	else
	{
		return(-1);
	}
}

Int_t Analyse::GetHLTriggerIndex(string triggername) const
{
	if(runlist.find(Run()) != runlist.end())
	{
		return(runlist.at(Run()).HLTIndex(triggername));
	}
	else
	{
		cerr << "ERROR GetHLTriggerIndex: Trigger " << triggername << " not found." << endl; 
		return(-1);
	}
}

TriggerSelection* Analyse::AddTriggerSelection(string id, vector<string> triggernames, bool useprescaled)
{
	TriggerSelection* triggerselection = new TriggerSelection(this, triggernames, useprescaled);
	triggerselections[id] = triggerselection;
	return(triggerselection);
}

TriggerSelection* Analyse::GetTriggerSelection(string id)
{
  return(triggerselections[id]);
}

Int_t Analyse::GetHLTrigger(vector<string> triggernames) const
{
	Int_t result = -1;
	for(UInt_t i = 0 ; i < triggernames.size() ; i++)
	{
		Int_t index = GetHLTriggerIndex(triggernames[i]);

		if(index != -1 && GetHLTPrescale(index) == 1)
		{
			result = 0;
			if(GetHLTrigger(index)) return(1);
		}
	}
	return(result);
}

string Analyse::GetHLTriggerName(UInt_t index) const
{
	if(runlist.find(Run()) != runlist.end())
	{
		return(runlist.at(Run()).HLTName(index));
	}
	else
	{
		cerr << "ERROR GetHLTriggerName: Could not find name of trigger with index " << index << endl;
		return(string("NotFound"));
	}
}

Int_t Analyse::GetHLTPrescale(UInt_t triggerindex) const
{
	if(IsLumiAvailable())
	{
		Int_t triggertable = lumilist.at(Run()).at(LumiBlock()).HLTTable();
		return(runlist.find(Run())->second.HLTPrescale(triggerindex, triggertable));
	}
	else
	{
		return(-1);
	}
}

void Analyse::PrintPrescaleInfo(string triggername)
{
	for(map< UInt_t , RunInfo>::iterator a = runlist.begin() ; a != runlist.end() ; ++a)
	{
		Int_t index = a->second.HLTIndex(triggername);
		if(index < 0) {cout << a->first << ": " << triggername << " not available." << endl;}
		cout << a->first << ":" << index;
		for(UInt_t i = 0 ; i < a->second.NumHLTTables() ; i++)
		{
			cout << " " << a->second.HLTPrescale(index, i);
		}
		cout << lumilist.at(a->first).size() <<  endl;
		for(map<UInt_t, Luminosity>::iterator b = lumilist.at(a->first).begin() ; b != lumilist.at(a->first).end() ; b++)
		{
			Int_t triggertable = b->second.HLTTable();
			Int_t prescale = a->second.HLTPrescale(index, triggertable);
			if(prescale != 1)
			{
				cout << "(" << a->first << ":" << prescale << "),";
			}
		}
		cout << endl;
	}

}


void Analyse::PrintPrescaleInfoB(string triggername)
{
	triggername += string("_v.*");
	boost::cmatch what;
	boost::regex triggernameregex = boost::regex(triggername.c_str());
	Int_t precur = -2;
	UInt_t runcur = 0;
	UInt_t lumicur = 0;
	UInt_t runlast = 0;
	UInt_t lumilast = 0;
	Double_t lumi = 0.;
	string triggernamecur;
	Int_t triggerindexcur;
	for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a)
	{
		string triggernamenow;
		RunInfo& runinfo = runlist.find(a->first)->second;
		for(UInt_t i = 0 ; i < runinfo.NumHLT() ; i++)
		{
			if(boost::regex_match(runinfo.HLTName(i).c_str(), what, triggernameregex))
			{
				triggernamenow = runinfo.HLTName(i);
				break;
			}
		}
		Int_t triggerindex = runinfo.HLTIndex(triggernamenow);
		for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++)
		{
			lumi+=b->second.LumiValue();
			Int_t triggertable = b->second.HLTTable();
			Int_t prenow = -1;
			if(triggerindex != -1)
			{
				prenow = runinfo.HLTPrescale(triggerindex, triggertable);
			}
	
			
			if(precur == -2)
			{
				triggernamecur = triggernamenow;
				triggerindexcur = triggerindex;
				precur = prenow;
				runcur = a->first;
				lumicur = b->first;
			}			

			if(precur != prenow || triggernamecur != triggernamenow || triggerindexcur != triggerindex)
			{
				cout << triggerindexcur << " " << triggernamecur << ", " << runcur << ":" << lumicur << " - " << runlast << ":" << lumilast << ", " << precur << ", Lumi: " << lumi << "/pb"<< endl; 
				lumi = 0.;
				triggerindexcur = triggerindex;
				triggernamecur = triggernamenow;
				precur = prenow;
				runcur = a->first;
				lumicur = b->first;
			}
			runlast = a->first;
			lumilast = b->first; 
		}
	}
	cout << triggerindexcur << " " << triggernamecur << ", " << runcur << ":" << lumicur << " - " << runlast << ":" << lumilast << ", " << precur << ", Lumi: " << lumi << "/pb"<< endl; 
}

Double_t Analyse::GetInstLumi()
{
return(lumilist[Run()][LumiBlock()].LumiValue());
}

Double_t Analyse::GetLumi(Int_t format)
{
	Double_t lumi = 0.;
	Double_t alllumi = 0.;
	Double_t zerolumi = 0.;
	UInt_t nolumiinfo = 0;
	Double_t events = 0.;
	Double_t eventsprocessed = 0.;
	for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a)
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
			if(b->second.LumiValue() == -1){continue;}
			if(IsInRange(a->first, b->first))
			{
				numblocks++;
			//	cout << a->first << " " << b->first << " " << minLumi << " " << maxLumi << endl;
				if(b->second == true)
				{
					if(b->second.NumEvents() > 0)
					{
						runlumi += b->second.ProcessedLumi();
						runalllumi += b->second.LumiValue();
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
						runzerolumi += b->second.LumiValue(); 
						runalllumi += b->second.LumiValue();
					}
				}
				else if(b->second == false)
				{
					runnolumiinfo += b->second.NumEventsProcessed();
				}
			}
		}
		if(a->first <= maxRun && a->first >= minRun)
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
		cout << "From run " << minRun << ", block " << minLumi << " to run " << maxRun << ", block " << maxLumi << "." << endl;
		cout << "Total lumi: " << lumi + zerolumi << " pb^-1, zero event lumi: " << zerolumi << " pb^-1, (lumi in Range: " << alllumi << " pb^-1), no lumi-info for " << nolumiinfo << " event(s), events: " << eventsprocessed << "(" << events << ")"  << endl;
	}
return(lumi+zerolumi);
}

bool Analyse::IsInRange(UInt_t theRun, UInt_t theLumiBlock)
{
	if(theRun < maxRun && theRun > minRun) 
	{
		return(true);
	}
	else if(maxRun != minRun && theRun == maxRun && theLumiBlock <= maxLumi)
	{
		return(true);
	}
	else if(maxRun != minRun && theRun == minRun && theLumiBlock >= minLumi)
	{
		return(true);
	}
	else if(maxRun == minRun && theRun == minRun && theLumiBlock <= maxLumi && theLumiBlock >= minLumi)
	{
		return(true);
	}
	else 
	{
		return(false);
	}
}

Double_t Analyse::GetLumiBlockLumi() 
{
}

void Analyse::PrintLumiOfRuns()
{
	Double_t totlumi = 0.;
	Double_t filteredevents = 0;
	Double_t allevents = 0;
	for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a)
	{
		Double_t runlumi = 0;
		Double_t runfilteredevents = 0;
		Double_t runallevents = 0;
		for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++)
		{
			if(b->second.LumiValue() != -1.)
			{
				runlumi+=b->second.LumiValue();
				totlumi+=b->second.LumiValue();
				filteredevents+=b->second.NumEvents();
				runfilteredevents+=b->second.NumEvents();
				allevents+=b->second.NumEventsOrig();
				runallevents+=b->second.NumEventsOrig();
			}
			else
			{
				cout << "Missing lumi info for section: " << b->first << endl;
			}
		}
		cout << "Run: " << a->first << " -> " << runlumi << ", Events: " << (UInt_t)runfilteredevents << "/" << (UInt_t)runallevents <<"="<< runfilteredevents/runallevents << "  sum: " << totlumi << ", Events: " << (UInt_t)filteredevents << "/" << (UInt_t)allevents <<"="<< filteredevents/allevents<< endl;
	}
}

void Analyse::PrintLumiOfLumiSectionsInRun(UInt_t runnumber)
{
	Double_t runlumi = 0.;
	Double_t runfilteredevents = 0.;
	Double_t runallevents = 0.;
	map<UInt_t, map<UInt_t, Luminosity> >::iterator blocklist = lumilist.find(runnumber);
	for(map<UInt_t, Luminosity>::iterator b = blocklist->second.begin() ; b != blocklist->second.end() ; ++b)
	{
		if(b->second.LumiValue() != -1.)
		{
			runlumi+=b->second.LumiValue();
			runfilteredevents+=b->second.NumEvents();
			runallevents+=b->second.NumEventsOrig();
			cout << "Section: " << b->first << ", Lumi: " << b->second.LumiValue() << "("<<runlumi<<"), Events: "<< b->second.NumEvents() << "(" << runfilteredevents << ")" << endl; 
		}
		else
		{
			cout << "Missing lumi info for section: " << b->first << endl;
		}
	}

}

bool Analyse::LoadJSON(string filename)
{
	jsonfilter = true;
	char buf[100];
	string strjson;
	fstream json(filename.c_str());
	while(true)
	{
		json.read(buf, 100);
		int read = json.gcount();
		strjson += string(buf, read);
		if(read < 100) break;
	}

	UInt_t run;
	UInt_t lumibeg;
	UInt_t lumiend;
	size_t posquota = 0;
	size_t posquotb = 0;
	size_t runbraon = 0;
	size_t runbraoff = 0;
	size_t lumibraon = 0;
	size_t lumibraoff = 0;
	size_t sep;
	bool result = true;
	while(true)
	{
		posquota = strjson.find("\"", runbraoff+1);
		posquotb = strjson.find("\"", posquota+1);
		if(posquota == string::npos) break;
		run = atoi(strjson.substr(posquota+1, posquotb-1).c_str());
		runbraon = strjson.find("[", posquotb+1);
		lumibraoff = runbraon;
		while(true)
		{
			lumibraon = strjson.find("[", lumibraoff+1);
			lumibraoff = strjson.find("]", lumibraoff+1);
			if(lumibraoff < lumibraon) break;
			sep = strjson.find(",", lumibraon+1);
			lumibeg = atoi(strjson.substr(lumibraon+1, sep-1).c_str());
			lumiend = atoi(strjson.substr(sep+1, lumibraoff-1).c_str());

			for(UInt_t u = lumibeg ; u <= lumiend ; u++)
			{
				if(lumilist.find(run) != lumilist.end() && lumilist.find(run)->second.find(u) != lumilist.find(u)->second.end())
				{
					jsonlist[run][u] = true;
				}
				else
				{
					cout << "JSONFilter: Run: " << run << ", LS: " << u << " is not part of your input files." << endl;
					result = false;
				}
			}  

		}
		runbraoff = lumibraoff;
	}
	return(result);
}

Double_t Analyse::GetPileUpWeight(Double_t meaninteractions) const
{
	if(!usepileupinfo) return(1.);
	return(Poisson(NumPileUpInteractions(), meaninteractions)/pileUpDist[NumPileUpInteractions()]);
}

Double_t Analyse::GetPrimVertexWeight(Double_t meaninteractions) const
{
	if(!useprimvertexinfo) return(1.);
	Int_t numgoodprimvertices = NumGoodPrimVertices();
	return(Poisson(numgoodprimvertices, meaninteractions)/primVertexDist[numgoodprimvertices]);
}

Double_t Analyse::GetPrimVertexWeight(vector<Double_t>& datadist) const
{
	if(!useprimvertexinfo) return(1.);
	Int_t numgoodprimvertices = NumGoodPrimVertices();
	Double_t dataval = 0.;
	if(datadist.size() > numgoodprimvertices)
	{
		dataval = datadist[numgoodprimvertices];
	}
	return(dataval/primVertexDist[numgoodprimvertices]);
}

Int_t Analyse::NumGoodPrimVertices() const
{
	return NumPrimVertices();
	/*Int_t numgoodprimvertices = 0;
	for(UInt_t i = 0 ; i < NumPrimVertices() ; i++)
	{
		if(PrimVertices(i).IsGood()) numgoodprimvertices++;
	}
	return(numgoodprimvertices);*/
}

unsigned int Analyse::GetOrigEvents() const
{
  assert(currentfile);
  TH1D* d = dynamic_cast<TH1D*>(currentfile->Get("nEvents"));
  if(!d) d = dynamic_cast<TH1D*>(currentfile->Get("makeroottree/nEvents"));
  if(!d) throw std::runtime_error("Failed to obtain nEvents histogram");
  
  return static_cast<unsigned int>(d->GetBinContent(1));
}
