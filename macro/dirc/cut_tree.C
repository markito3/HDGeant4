#define glx__sim
#include "DrcHit.h"
#include "DrcEvent.h"
#include "glxtools.C"
#include "TSystemDirectory.h"
#include "TSystemFile.h"


// void cut_tree(TString infile="/volatile/halld/home/gxproj7/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass04/hists",TString out="ver08_pass04_1c.root"){
void cut_tree(TString infile="/volatile/halld/home/jrsteven/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass05/merged/hd_root_*.root",TString out="ver08_pass05_1c.root"){
  if(!glx_initc(infile,0,"")) return;

  if(!infile.EndsWith("root") && !infile.Contains("*")){
    TSystemDirectory dir(infile, infile);
    TList *files = dir.GetListOfFiles();
    if(files){
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
	fname = file->GetName();
	if(file->IsDirectory() && fname.BeginsWith("06")) {
	  std::cout<<"adding "<<infile+"/"+fname+"/*.root"<<std::endl;	  
	  glx_ch->Add(infile+"/"+fname+"/*.root");	  
	}
      }
    }
    std::cout<<"Entries "<<glx_ch->GetEntries()<<std::endl; 
  }
  
  TFile *newFile = new TFile(out,"recreate");
  TTree *newTree = new TTree("dirc", "dirc tree");
  TClonesArray *newEvent = new TClonesArray("DrcEvent");
  newTree->Branch("DrcEvent",&newEvent,256000,2);

  for(int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);

    TClonesArray& cevt = *newEvent;
    cevt.Clear();      

    for(int t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,10000);      
      double mom = glx_event->GetMomentum().Mag();

      // custom selection
      const int nbins=20;
      // if(mom<3 || mom>4.5 ) continue;
      int bar = glx_event->GetId();
      int bin = (100+glx_event->GetPosition().X())/200.*nbins;
      if((bar<0 || bar>23)) continue;
      if((bin<9 || bin>11)) continue; 
      if(fabs(glx_event->GetPdg())==2212) continue;
      

      if(glx_entries<5) continue;
      if(fabs(glx_event->GetPdg())==211){
	if(fabs(glx_event->GetInvMass()-0.77)>0.1) continue;
	if(glx_event->GetChiSq()>5) continue;
      }else if(fabs(glx_event->GetPdg())==321){
	if(fabs(glx_event->GetInvMass()-1.02)>0.05) continue;
	if(glx_event->GetChiSq()>20) continue;
      }else if(fabs(glx_event->GetPdg())==2212){
	if(glx_event->GetChiSq()>20) continue;
      }else continue;

      new (cevt[cevt.GetEntriesFast()]) DrcEvent(*glx_event);
    }    
    if(cevt.GetEntriesFast()>0) newTree->Fill();    
  }
  newTree->AutoSave();

}
