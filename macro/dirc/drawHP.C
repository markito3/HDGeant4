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
      if(glx_event->GetType()!=2) continue; //1-LED 2-beam 0-rest
      if(glx_event->GetMomentum().Mag()<4) continue;
      
      for(auto hit : glx_event->GetHits()){
  	int ch = hit.GetChannel();
    	int pmt = hit.GetPmtId();
    	int pix = hit.GetPixelId();
    	double time = hit.GetLeadTime();	

  	if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8);
      }

      // if(glx_event->GetHitSize()>5){
      // 	glx_drawDigi("",0,1,0);
      // 	glx_canvasAdd(glx_cdigi);
      // 	glx_cdigi->SetName(Form("evt_cdigi_%d",e)+infile);    
      // 	glx_canvasSave(0);
      // 	glx_resetDigi();
      // }
    }
  }
  
  glx_drawDigi("",0);
  glx_canvasAdd(glx_cdigi);
  glx_cdigi->SetName(TString("cdigi_")+infile);
  glx_canvasSave(1,0);
}
