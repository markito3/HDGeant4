#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

//void drawMult(TString infile="KalongBar0/kaons1.root"){
//void drawMult(TString infile="data/sim_etaprime2300.root"){
void drawPlots(TString infile="data/sim_hprime2600_100k_tr_nodc.root"){
  if(!glx_initc(infile,1,"data/drawHP")) return;

  //  gStyle->SetOptStat(0);
  TH2F *hPoint = new TH2F("hPoint",";x [cm]; y [cm]",200,-120,120,200,-120,120);
  TH1F *hMult = new TH1F("hMult",";detected photons [#]; [#]",500,0,500);
  TH2F* hMult2 = new TH2F("hMult2","; x [cm]; y[cm]", 24, -120.,120.,24,-120.,120.);
  
  const auto nmax(40);
  double minx=-100, maxx=100;
  TH1F *hMultX[nmax];
  for(auto i=0; i<nmax; i++){
    hMultX[i] = new TH1F(Form("hMultX_%d",i),Form("hMultX_%d;N photons [cm]; stat [#]",i),300,0,300);
  }
  
  double miny=-100, maxy=100;
  TH1F *hMultXY[nmax][nmax];
  for(auto i=0; i<nmax; i++){
    for(auto j=0; j<nmax; j++){
      hMultXY[i][j] = new TH1F(Form("hMultXY_%d_%d",i,j),Form("hMultXY_%d_%d;x [cm]; stat [#]",i,j),100,0,200);
    }
  }
  
  TVector3 hpos,gpos;
  DrcHit hit;
  for (auto e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (auto t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,100);
      if(glx_event->GetParent()>0) continue;
      hpos = glx_event->GetPosition();
      double x(hpos.X()), y(hpos.Y());
      //      cout<<"x = "<<x<<", y = "<<y<<endl;
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
      if(fabs(fabs(y)-12)<4){
	int xid = int(nmax*(x-minx)/(maxx - minx));
	if(xid>=0 && xid<nmax) hMultX[xid]->Fill(nhits);
      }
      int xid = std::round(nmax*(x-minx)/(maxx - minx));
      int yid = std::round(nmax*(y-miny)/(maxy - miny));
      //	cout<<"xid = "<<xid<<", yid = "<<yid<<", nhits = "<<nhits<<endl;
      if(xid >= 0 && yid >=0 && xid<nmax && yid < nmax){
	hMultXY[xid][yid]->Fill(nhits);
	if(fabs(yid) == 5 && x < 50. && x > -50.){
	  cout<<"x = "<<x<<", xid = "<<xid<<", y = "<<y<<", yid = "<<yid<<endl;
	}
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
  
  glx_canvasAdd("hMultX");

  //glx_canvasAdd("cMultXY");
  
  TGraph *gMult = new TGraph();
  for(auto i=0; i<nmax; i++){
    double nph = glx_fit(hMultX[i],40,50,30).X();
    // hMultX->cd();
    hMultX[i]->Draw();

    double xpos = minx + 0.5*(maxx - minx)/nmax + i*(maxx - minx)/nmax;
    //    glx_waitPrimitive("cMultX");
    gMult->SetPoint(i,xpos,nph);
   
    for(auto j=0; j<nmax; j++){
      double nph2 = glx_fit(hMultXY[i][j],40,50,30).X();
      double ypos = miny + 0.5*(maxy - miny)/nmax + j*(maxy - miny)/nmax;
      //cMultXY->cd();
      //hMultXY[i][j]->Draw();
      if(int((ypos+120.)/10.)+1 == 12 || int((ypos+120.)/10.)+1 == 13){
	cout<<"xpos = "<<xpos<<"ypos = "<<ypos<<endl;
	cout<<int((ypos+120.)/10.)+1<<" "<<int((xpos+120.)/10.)+1<<endl;
      }
      hMult2->SetBinContent(int((xpos+120.)/10.)+1, int((ypos+120.)/10.)+1, nph2);
    }
   }
 
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
  hMult2->Draw("hcol");
  
}


