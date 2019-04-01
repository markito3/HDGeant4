#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawTime(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawTime")) return;

  int ch;
  double w,mean[glx_nch]={0};
  ifstream in("mean.txt");
  if (in && !in.bad()){
    while (!in.eof()) {
      in>>ch>>w;
      mean[ch]=w;  
    }
  }
  in.close();


  TH1F *htime[glx_nch];
  for(int i=0; i<glx_nch; i++){
    htime[i] = new TH1F(Form("h_%d",i),";time [ns];entries [#]",500,110,160);
    //htime[i] = new TH1F(Form("h_%d",i),";time [ns];entries [#]",200,124,126);
  }
  TH1F *hSigma = new TH1F("hSigma",";#sigma [ns];entries [#]",200,0,3);
  TH1F *hMean = new TH1F("hMean",";mean [ns];entries [#]",500,110,160);
  //TH1F *hMean = new TH1F("hMean",";mean [ns];entries [#]",200,124,126);
  
  for (int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (int t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,1000);
      if(glx_event->GetParent()>0) continue;
      double etime = glx_event->GetTime();
      
      for(auto hit : glx_event->GetHits()){	
	int ch = hit.GetChannel();
    	double time = hit.GetLeadTime()-etime; //-mean[ch]+125;
	double tot = hit.GetTotTime();
	//if(tot<40 || tot>42) continue;
	if((ch/64)%18>6) time -=10;
	if((ch/64)%18>12) time -=10; 
	htime[ch]->Fill(time);
      }      
    }
  }

  TString str="";
  TCanvas *cc = new TCanvas("cc","cc",800,500);
  for(int i=0; i<glx_nch; i++){
    auto res = glx_fit(htime[i],3,20,10,1,"0");
    hMean->Fill(res.X());
    if(res.Y()>0.01) hSigma->Fill(res.Y());
    int pix=i%64;
    glx_hdigi[i/64]->Fill(pix/8,pix%8,res.X());
    // htime[i]->Draw();
    // glx_waitPrimitive(cc);
    str += Form("%d %f\n", i, res.X());
  }
  //glx_writeString("mean.txt",str);

  
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  
  TString nid=infile+"_";
  glx_canvasAdd(nid+"hmean",800,400);
  // hMean->Fit("gaus");
  glx_fit(hMean);
  hMean->Draw();
  glx_canvasAdd(nid+"hsigma",800,400);
  hSigma->Draw();
  
  glx_drawDigi("",0,-2,-2);
  glx_canvasAdd(glx_cdigi);
  glx_cdigi->SetName(nid+"cdigi_mean_diff_");
  glx_canvasSave(0);
}
