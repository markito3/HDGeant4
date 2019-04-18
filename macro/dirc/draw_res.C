#include "glxtools.C"

void draw_res(TString in="res_12.4_1.root"){

  glx_savepath = "data/draw_res1";
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(99);

  double thr,dac,nph;
  
  TChain ch("reco"); ch.Add(in);  
  ch.SetBranchAddress("thr",&thr);
  ch.SetBranchAddress("dac",&dac);
  ch.SetBranchAddress("nph",&nph);

  TGraph *gMult = new TGraph();
  //g->SetLineColor(c[id%5]);
  //g->SetMarkerColor(c[id%5+5]);
  gMult->SetMarkerStyle(20);
  gMult->SetMarkerSize(0.7);
  //g->SetLineStyle(id>4 ? 2:1);
  gMult->SetTitle(";DAC [#];multiplicity [#]");
  
  
  for(int i=0; i<ch.GetEntries(); i++){
    ch.GetEvent(i);
    gMult->SetPoint(i,dac,nph);
    std::cout<<"i "<<i << " "<<dac<<std::endl;
    
  }

  gMult->Sort();

  glx_canvasAdd("gMult",800,500);
  gMult->Draw("apl");
  gMult->GetXaxis()->SetRangeUser(200,1100);
  
  glx_canvasSave(0);
}
