#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawMult_data(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawMult_data_fiber")) return;

  int thr(0),dac(0), pmt(108),hv(0);
  TObjArray *sarr = infile.Tokenize("_");
  if(sarr->GetEntries()==5){
    TString s1 = ((TObjString *) sarr->At(2))->GetName();
    thr = s1.Atof();
    TString s2 = ((TObjString *) sarr->At(3))->GetName();
    dac = s2.Atof();
  }
  if(sarr->GetEntries()==6){
    TString s1 = ((TObjString *) sarr->At(2))->GetName();
    thr = s1.Atof();
    TString s2 = ((TObjString *) sarr->At(3))->GetName();
    dac = s2.Atof();
    TString s3 = ((TObjString *) sarr->At(4))->GetName();
    hv = s3.Atoi();
  }
  std::cout<<"dac "<<dac << " "<<thr<<std::endl;
  

  glx_createMap();
  TH1F *hMultA = new TH1F("hMultA","hMultA",2000,0,2000);
  TH1F *hMult[glx_npmt];
  for(int i=0; i<glx_npmt; i++){
    hMult[i] = new TH1F(Form("pmt_%d",i),Form("pmt_%d",i),2000,0,10000);
  }
  
  TH1F *hMult1 = new TH1F("hMult1",";multiplicity per 10 triggers; entries [#]",70,0,70);
  TH1F *hMult2 = new TH1F("hMult2","",70,0,70);
  TH1F *hMult3 = new TH1F("hMult3","",70,0,70);
  
  int evt=0,mult1=0,mult2=0,mult3=0,mult[glx_npmt]={0};
  for(int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for(int t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,1000);      
      if(glx_event->GetParent()>0) continue;      
      hMultA->Fill(glx_event->GetHitSize());
      if(glx_event->GetType()!=15) continue; //LED
      evt++;
      for(auto hit : glx_event->GetHits()){
      	int ch = hit.GetChannel();
      	int pmt = hit.GetPmtId();
      	int pix = hit.GetPixelId();
      	double time = hit.GetLeadTime();
	double tot = hit.GetTotTime();
	if(time>50 && time<100)// && tot>20 && tot<70)
	  mult[pmt]++;
	if(pmt==76 || pmt==58) mult1++;
	if(pmt==82 || pmt==64) mult2++;
	if(pmt==88 || pmt==70) mult3++;

	// if(pmt==73 || pmt==55) mult1++;
	// if(pmt==61 || pmt==79) mult2++;
	// if(pmt==85 || pmt==67) mult3++;

	
      	//if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, pix/8);
      }
      if(evt%10==0){
	for(int i=0; i<glx_npmt; i++) if(mult[i]>0) hMult[i]->Fill(mult[i]);
	memset(mult, 0, sizeof(mult));
	hMult1->Fill(mult1);
	hMult2->Fill(mult2);
	hMult3->Fill(mult3);
	mult1=0;
	mult2=0;
	mult3=0;
      }
    }
  }


  glx_canvasAdd("hMult_fiber123_right_led",800,500);
  gStyle->SetOptStat(0);
  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  glx_normalize(hMult1,hMult3);
  glx_normalize(hMult2,hMult1);
  
  hMult1->SetLineColor(kRed+1);
  hMult1->Draw();
  hMult2->SetLineColor(kBlue+1);
  hMult2->Draw("same");
  hMult3->SetLineColor(kGreen+2);
  hMult3->Draw("same");
  
  leg->AddEntry(hMult1,"fiber 1 right, ssp 5","lp");
  leg->AddEntry(hMult2,"fiber 2 right, ssp 5","lp");
  leg->AddEntry(hMult3,"fiber 3 right, ssp 5","lp");
  leg->Draw();

  
  glx_canvasAdd("hMultA",800,500);
  hMultA->Draw();
  double nph = glx_fit(hMultA,40,50,30).X();

  TString rpath = infile;
  rpath.ReplaceAll(".root",".res.root");
  TFile fc(rpath,"recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("thr",&thr,"thr/I");
  tc->Branch("dac",&dac,"dac/I");
  tc->Branch("nph",&nph,"nph/D");
  tc->Branch("nph",&nph,"nph/D");
  tc->Branch("pmt",&pmt,"pmt/I");
  tc->Branch("hv",&hv,"hv/I");
  tc->Fill();
  
  for(int i=0; i<glx_npmt; i++){
    //nph = hMult[i]->GetMean();
    nph = glx_fit(hMult[i],500,50,500,1).X();
    pmt=i;
    tc->Fill();
    // hMult[i]->Draw();
    // glx_waitPrimitive("hMultA");
  }
  
  tc->Write();
  fc.Write();
  fc.Close();
  glx_canvasSave(1,0);    
}
