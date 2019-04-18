#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void cut_tree(TString infile="vol/ver08/ver08/tree/all_tree.root",TString out="beam_tree_23.root"){
  if(!glx_initc(infile,0,"")) return;
  
  TFile *newfile = new TFile(out,"recreate");
  TTree *newtree = glx_ch->CloneTree(0);
  int count=0;
  for(int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);

    bool good = false;
    for(int t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,10000);    
      if(glx_event->GetParent()>0) continue;
      if(glx_event->GetType()!=2) continue; //1-LED 2-beam 0-rest

      if (fabs(glx_event->GetPdg())==321){
	count++;	      
	good=true;
      }
            
      // for(auto hit : glx_event->GetHits()){
      // 	int ch = hit.GetChannel();
      // 	int pmt = hit.GetPmtId();
      // 	int pix = hit.GetPixelId();
      // 	double time = hit.GetLeadTime();
      // }

    }
    
    if(good || e%25==0){
      newtree->Fill();
      glx_event->Clear();
      //if(count>300000) break;
    }
    
  }
  newtree->Print();
  newtree->AutoSave();

}
