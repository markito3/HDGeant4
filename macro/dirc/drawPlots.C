#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawPlots(TString infile="KalongBar77/KalongBar77.root"){
//void drawPlots(TString infile="data/etaprime2300_all.root"){
//void drawPlots(TString infile="data/sim_hprime2600_100k_tr_nodc.root"){

  TFile *f = new TFile("drawPlots_out.root","recreate");
  
  if(!glx_initc(infile,1,"data/drawHP")) return;

  //  gStyle->SetOptStat(0);
  TH2F *hPoint = new TH2F("hPoint",";x [cm]; y [cm]",200,-120,120,200,-120,120);
  TH1F *hMult = new TH1F("hMult",";detected photons [#]; [#]",500,0,500);
  TH2F* hMult2 = new TH2F("hMult2","; x [cm]; y[cm]", 24, -120.,120.,24,-120.,120.);
  TH1F* hLam = new TH1F("hLam","; wavelength [nm]", 200, 0., 1000.);

  TH1F* hPdg = new TH1F("hPdg","; pdg [#]",5000,0.,5000.);
  TH2F* hPTh1 = new TH2F("hPTh1","; pion momentum [GeV/c]; polar angle [degrees]",100,0.,10.,100,0.,15.);
  TH2F* hPTh2 = new TH2F("hPTh2","; kaon momentum [GeV/c]; polar angle [degrees]",100,0.,10.,100,0.,15.);
  TH2F* hPTh3 = new TH2F("hPTh3","; proton momentum [GeV/c]; polar angle [degrees]",100,0.,10.,100,0.,15.);
  
  const auto nmax(40);
  double minx=-95, maxx=95;
  TH1F *hMultX[nmax];
  for(auto i=0; i<nmax; i++){
    hMultX[i] = new TH1F(Form("hMultX_%d",i),Form("hMultX_%d;N photons [cm]; stat [#]",i),100,0,300);
  }
  
  double miny=-100, maxy=100;
  TH1F *hMultXY[nmax][nmax];
  for(auto i=0; i<nmax; i++){
    for(auto j=0; j<nmax; j++){
      hMultXY[i][j] = new TH1F(Form("hMultXY_%d_%d",i,j),Form("hMultXY_%d_%d;x [cm]; stat [#]",i,j),100,0,200);
    }
  }
  
  TVector3 hpos,gpos;
  TVector3 mom;
  DrcHit hit;
  Double_t theta = 0.;
  Int_t pdg = 0.;
  for (auto e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (auto t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,100);
      if(glx_event->GetParent()>0) continue;
      mom = glx_event->GetMomentum();
      theta = TMath::ACos(mom.Z()/mom.Mag())/3.1415*180.;
      pdg = glx_event->GetPdg();
      // choose only kaons with p> 3
      //      if(pdg != 321 && mom.Mag() < 3.)continue;
      //      if(pdg != 211) continue;
      
      hPdg->Fill(pdg);
      if(pdg == 211){
	hPTh1->Fill(mom.Mag(), theta);
      }else if (pdg == 321){
	hPTh2->Fill(mom.Mag(), theta);
      }else if (pdg == 2212){
	hPTh3->Fill(mom.Mag(), theta);
      }
      hpos = glx_event->GetPosition();
      double x(hpos.X()), y(hpos.Y());
      //      cout<<"x = "<<x<<", y = "<<y<<endl;
      hPoint->Fill(x, y);
      
      int nhits=glx_event->GetHitSize();
      for(auto h=0; h<nhits; h++){
	// cut out two long rows of PMTs:
	//	if(pmt > 17 && pmt < 90){
	  hit = glx_event->GetHit(h);
	  Int_t pmt = hit.GetPmtId();
	  Int_t pix = hit.GetPixelId();
	  gpos = hit.GetPosition();
	  Double_t time = hit.GetLeadTime();
	  Double_t E = hit.GetEnergy();
	  hLam->Fill(1239.84/(E*1E9));
	  if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, 7-pix/8);
	  //if(pmt>=108) glx_hdigi[pmt-108]->Fill(pix%8, 7-pix/8);
	  //}
      }

      hMult->Fill(nhits);
      if(fabs(fabs(y)-12)<4){
	int xid = std::round(nmax*(x-minx)/(maxx - minx));
	if(xid>=0 && xid<nmax) hMultX[xid]->Fill(nhits);
      }
      int xid = std::round(nmax*(x-minx)/(maxx - minx));
      int yid = std::round(nmax*(y-miny)/(maxy - miny));
      if(xid >= 0 && yid >=0 && xid<nmax && yid < nmax){
	hMultXY[xid][yid]->Fill(nhits);
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

  glx_canvasAdd("hMult",800,400);
  hMult->Draw();
  
  glx_canvasAdd("cMultXY");
  
  TGraph *gMult = new TGraph();
  for(auto i=0; i<nmax; i++){
    //    double nph = glx_fit(hMultX[i],40,50,30).X();
    double nph = hMultX[i]->GetBinCenter(hMultX[i]->GetMaximumBin());
    hMultX[i]->Draw();
    double xpos = minx + 0.5*(maxx - minx)/nmax + i*(maxx - minx)/nmax;
    gMult->SetPoint(i,xpos,nph);
    cout<<"xpos = "<<xpos<<", nph - "<<nph<<endl;
    for(auto j=0; j<nmax; j++){
      if(hMultXY[i][j]->GetEntries() < 1)continue;
      //      double nph2 = hMultXY[i][j]->GetBinCenter(hMultXY[i][j]->GetMaximumBin());
           double nph2 = glx_fit(hMultXY[i][j],40,50,30).X();
      double ypos = miny + 0.5*(maxy - miny)/nmax + j*(maxy - miny)/nmax;
      //hMultXY[i][j]->Draw();
      //      glx_waitPrimitive("cMultX");
      //  hMultXY[i][j]->Write();
      hMult2->SetBinContent(int((xpos+120.)/10.)+1, int((ypos+120.)/10.)+1, nph2);
    }
   }

  // f->Close();
      glx_canvasAdd("x1");
      hMultX[3]->Draw();
      glx_canvasAdd("x2");
      hMultX[30]->Draw();

  
  glx_canvasAdd("cMultX",800,400);
  gMult->GetXaxis()->SetRangeUser(-110,110);
  gMult->GetYaxis()->SetRangeUser(0,150);
  gMult->GetXaxis()->SetTitle("x [cm]");
  gMult->GetYaxis()->SetTitle("detected photons [#]");
  gMult->SetMarkerStyle(20);
  gMult->SetMarkerSize(0.8);
  gMult->Draw("APL");
  
  glx_canvasSave(1,0);

  glx_canvasAdd("cMultXY", 800,800);
   hMult2->Draw("hcolz");

  glx_canvasAdd("cLam", 800,800);
  hLam->Draw();

  glx_canvasAdd("pdg");
  hPdg->Draw();
  
  glx_canvasAdd("ptheta_pion",800,800);
  hPTh1->Draw("hcolz");

  glx_canvasAdd("ptheta_kaon",800,800);
  hPTh2->Draw("hcolz");

  glx_canvasAdd("ptheta_proton",800,800);
  hPTh3->Draw("hcolz");
  
}



