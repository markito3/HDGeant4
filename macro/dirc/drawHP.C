#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawHP(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawHP")) return;
    
  for(int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for(int t=0; t<glx_events->GetEntriesFast(); t++){      
      glx_nextEventc(e,t,100);
      if(e > glx_ch->GetEntries()-2) cout<<"particle is "<<glx_names[glx_findPdgId(glx_event->GetPdg())]<<endl;
      if(glx_event->GetParent()>0) continue;
      for(auto hit : glx_event->GetHits()){
  	int ch = hit.GetChannel();
    	int pmt = hit.GetPmtId();
    	int pix = hit.GetPixelId();
    	double time = hit.GetLeadTime();	

	if(time>70) continue;
  	if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8,1/(double)glx_ch->GetEntries());
  	//if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8);
      }      
    }
  }
  
  glx_drawDigi("",0,0.0001);
  glx_canvasAdd(glx_cdigi);
  glx_cdigi->SetName(TString("cdigi_")+infile);
  glx_canvasSave(1,0);
}
