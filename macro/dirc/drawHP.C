#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawHP(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawHP")) return;

  
  TH1F *h = new TH1F("h","h",7000,0,7000);

  int c0[7000]={0};
  int c1[7000]={0};
  int count=0;
  DrcHit hit;
  for (Int_t e=0; e<glx_ch->GetEntries() && e<5000; e++){
    glx_ch->GetEntry(e);
    for (Int_t t=0; t<glx_events->GetEntriesFast(); t++){      
      glx_nextEventc(e,t,100);
      if(e > glx_ch->GetEntries()-2) cout<<"particle is "<<glx_names[glx_findPdgId(glx_event->GetPdg())]<<endl;
      if(glx_event->GetParent()>0) continue;
      for(Int_t h=0; h<glx_event->GetHitSize(); h++){
    	hit = glx_event->GetHit(h);
	Int_t ch = hit.GetChannel();
    	Int_t pmt = hit.GetPmtId();
    	Int_t pix = hit.GetPixelId();
    	TVector3 gpos = hit.GetPosition();
    	Double_t time = hit.GetLeadTime();
	
	if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8,1/(double)5000);
	//if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8);
	
	if(hit.GetTotTime()==0) c0[ch]++;
	if(hit.GetTotTime()==1) c1[ch]++;	
      }      
    }
  }
  
  // for(int i=0; i<7000; i++) if(c1[i]>0) h->Fill(i,c0[i]/(double)c1[i]);
  // glx_canvasAdd("ratio",800,500);
  // h->Draw("hist");
  
  glx_drawDigi("m,p,v\n",0);
  glx_canvasAdd(glx_cdigi);
  glx_cdigi->SetName(TString("cdigi_")+infile);
  // glx_canvasSave(1,0);
}
