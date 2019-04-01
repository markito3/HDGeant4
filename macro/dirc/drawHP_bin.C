#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawHP_bin(TString infile="drc.root"){
  if(!glx_initc(infile,1,"data/drawHP_bin")) return;
  TH2F *hPos = new TH2F("hPos",";x [cm];y [cm]",200,-100,100,100,-100,0);
  TH1F *hTime = new TH1F("hTime",";propagation time [ns];entries [#]",300,0,300);
  
  const int xbins=40;
  const int ybins=24;
  int count[ybins][xbins]={0};  
  long int *occup[ybins][xbins];
  for(int y=0; y<ybins; y++){
    for(int x=0; x<xbins; x++){
      occup[y][x] = new long int[glx_nch];
      for(int c=0; c<glx_nch; c++){
      	occup[y][x][c]=0;
      }
    }
  }
  
  TVector3 pos,mom;
  for(int e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for(int t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,10000);
      if(glx_event->GetParent()>0) continue;
      //if(glx_event->GetType()!=2) continue; //1-LED 2-beam 0-rest

      pos = glx_event->GetPosition();
      mom = glx_event->GetMomentum();
      int y = glx_event->GetId(); // bar

      if(mom.Mag()<4.0) continue;
      int x = (100+pos.X())/200.*xbins;
      if(x<0 || x>xbins) continue;
      if(y<0 || y>ybins) continue;

      hPos->Fill(pos.X(),pos.Y());

      bool good=0;
      for(auto hit : glx_event->GetHits()){
  	int ch = hit.GetChannel();
    	int pmt = hit.GetPmtId();
    	int pix = hit.GetPixelId();
    	double time = hit.GetLeadTime();	
	if(pmt<=10 || (pmt>=90 && pmt<=96)) continue; // dummy pmts
  	if(pmt<108){
	  hTime->Fill(time);
	  occup[y][x][ch]++;
	  //glx_hdigi[pmt]->Fill(pix%8, pix/8);
	  good=1;
	}
      }
      if(good) count[y][x]++;

      // if(glx_event->GetHitSize()>5){
      // 	glx_drawDigi("",0);
      // 	glx_canvasAdd(glx_cdigi);
      // 	glx_cdigi->SetName(Form("evt_cdigi_%d",e)+infile);    
      // 	glx_canvasSave(1,0);
      // 	glx_resetDigi();
      // }
    }
  }

  std::cout<<"xbins "<<xbins << " ybins "<<ybins<<std::endl;
  
  for(int y=0; y<ybins; y++){
    for(int x=0; x<xbins; x++){      
      if(count[y][x]>0){
      	for(int i=0; i<glx_nch; i++){
      	  int pmt=i/64;
      	  int pix=i%64;
	  double w = occup[y][x][i]/(double)count[y][x];
      	  glx_hdigi[pmt]->SetBinContent(pix%8+1, pix/8+1,w);
      	}
      }
      glx_drawDigi("",0);
      glx_canvasAdd(glx_cdigi);
      glx_cdigi->SetName(Form("cdigi_bin_%d_%d_beam",y,x));    
      glx_canvasSave(1,0);
      glx_resetDigi();
    }
  }
  
  gStyle->SetOptStat(0);

  TString nid="_beam";
  glx_canvasAdd("time"+nid,800,400);
  hTime->Draw();

  glx_canvasAdd("pos"+nid,800,400);
  hPos->Draw("colz");
  
  glx_drawDigi("",0);
  glx_canvasAdd(glx_cdigi);
  glx_cdigi->SetName(TString("cdigi"+nid));
  glx_canvasSave(0);
}
