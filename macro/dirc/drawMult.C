#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

//void drawMult(TString infile="data/out100k.root"){
//void drawMult(TString infile="data/sim_etaprime2300.root"){
void drawMult(TString infile="data/sim_hprime2600.root"){
  if(!glx_initc(infile,1,"data/drawHP")) return;

  //  gStyle->SetOptStat(0);
  TH2F *hPoint = new TH2F("hPoint",";x [cm]; y [cm]",200,-120,120,200,-120,120);
  TH1F *hMult = new TH1F("hMult",";detected photons [#]; [#]",500,0,500);
  
  const auto xmax(100);
  TH1F *hMultX[xmax];
  for(auto i=0; i<xmax; i++){
    hMultX[i] = new TH1F(Form("hMultX_%d",i),Form("hMultX_%d;x [cm]; stat [#]",i),500,0,500);
  }
  
  TVector3 hpos,gpos;
  DrcHit hit;
  for (auto e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (auto t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,100);
      // if(glx_event->GetParent()>0) continue;
      hpos = glx_event->GetPosition();
      double x(hpos.X()), y(hpos.Y());      
      hPoint->Fill(x, y);
      
      int nhits=glx_event->GetHitSize();
      for(auto h=0; h<nhits; h++){
    	hit = glx_event->GetHit(h);
    	Int_t pmt = hit.GetPmtId();
    	Int_t pix = hit.GetPixelId();
    	gpos = hit.GetPosition();
    	Double_t time = hit.GetLeadTime();
    	if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, 7-pix/8);
    	//if(pmt>=108) glx_hdigi[pmt-108]->Fill(pix%8, 7-pix/8);
      }

      hMult->Fill(nhits);
      if(fabs(y-12)<4){
	int xid=50+x/2.;
	if(xid>=0 && xid<100) hMultX[xid]->Fill(nhits);
      }
    }
  }

  TGaxis::SetMaxDigits(2);
  glx_drawDigi();
  glx_canvasAdd(glx_cdigi);

  glx_canvasSave(1,0);
  
  TGaxis::SetMaxDigits(4);
  glx_canvasAdd("hPoint",500,500);
  hPoint->Draw("colz");

  glx_canvasAdd("hMult",800,500);
  hMult->Draw();
  
  // glx_canvasAdd("hMultX");
  // for(auto i=0; i<xmax; i++){
  //   hMultX[i]->Draw();
  //   glx_waitPrimitive("hMultX");
  // }
  glx_canvasSave(1,0);
  
}

