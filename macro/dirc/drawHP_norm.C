#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawHP_norm(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawHP_norm")) return;
  
  int ch;
  double w,weights[7000]={0};
  ifstream in("weights_magnet_0.txt");
  if (in && !in.bad()){
    while (!in.eof()) {
      in>>ch>>w;
      weights[ch]=w;  
    }
  }
  in.close();
  
  for (int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (int t=0; t<glx_events->GetEntriesFast(); t++){      
      glx_nextEventc(e,t,1000);
      if(glx_event->GetParent()>0) continue;
      for(auto hit : glx_event->GetHits()){
  	int ch = hit.GetChannel();
    	int pmt = hit.GetPmtId();
    	int pix = hit.GetPixelId();
    	double time = hit.GetLeadTime();
	
  	if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8,1/(double)glx_ch->GetEntries());
      }      
    }
  }
  
  TString str="";
  for(int i=0; i<108*64; i++){
    int pmt=i/64;
    int pix=i%64;
    double rel = glx_hdigi[pmt]->GetBinContent(pix%8+1,pix/8+1)/weights[i];
    if(rel>0.00001) glx_hdigi[pmt]->SetBinContent(pix%8+1, pix/8+1,rel);
    //str += Form("%d %f\n", i, glx_hdigi[pmt]->GetBinContent(pix%8+1,pix/8+1));
  }
  //glx_writeString("weights_magnet_0.txt",str);
  
  glx_drawDigi("",0,1.2,0.8);
  glx_canvasAdd(glx_cdigi);
  glx_cdigi->SetName(TString("cdigi_")+infile);
  //glx_canvasSave(1,0);
}
