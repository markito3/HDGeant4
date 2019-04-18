#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawMult_timeline(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawMult_timeline")) return;

  TH1F *hMult[glx_nch],*hMultFit[glx_nch];    
  for(int i=0; i<glx_nch; i++){
    hMult[i] = new TH1F(Form("ch_%d",i),Form("ch %d",i),200,0,200);
    hMultFit[i] = new TH1F(Form("fch_%d",i),Form("ch %d",i),200,0,200);
  }
  glx_canvasAdd("hMultFit",800,500);    
  
  int step=0,mult[glx_nch]={0};
  int chs[]={2240,2247};
  for (int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (int t=0; t<glx_events->GetEntriesFast(); t++){      
      glx_nextEventc(e,t,1000);      
      if(glx_event->GetParent()>0) continue;
      for(auto hit : glx_event->GetHits()){
  	int ch = hit.GetChannel();
	mult[ch]++;
      }
      if(e%2000==0){
	for(auto i : chs) hMultFit[i]->Fill(mult[i]);
	memset(mult, 0, sizeof(mult));
      }
      if(e%200000==0){
	for(auto i : chs){
	  double nph = glx_fit(hMultFit[i],500,50,500,1).X();
	  hMultFit[i]->Draw();
	  glx_waitPrimitive("hMultFit");
	  hMult[i]->Fill(step,nph);
	}
	step++;
	memset(mult, 0, sizeof(mult));
      }

    }
  }
  
  glx_canvasAdd("mult",800,500);
  hMult[2240]->Draw("h");
  hMult[2247]->Draw("h same");
  glx_canvasSave(1,0);
}
