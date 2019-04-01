#include "glxtools.C"

void draw_res_hv(TString in="run_12.4_400_hvscan.root"){

  glx_savepath = "data/draw_res_hv_1";
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(99);

  double nph;
  int thr,dac,pmt,hv;
  
  TChain ch("reco"); ch.Add(in);  
  ch.SetBranchAddress("thr",&thr);
  ch.SetBranchAddress("dac",&dac);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("pmt",&pmt);
  ch.SetBranchAddress("hv",&hv);



  TGraph *gMult[6], *gMultP[glx_npmt+1];
  int hvi[]={850,900,950,1000,1050,1100};
  int colors[]={1,kGreen+1,kRed+2,kRed,4,kMagenta+2,6,7,8,9,10};
  int hvid[2000]={0};
  for(int i=0; i<6; i++){
    gMult[i]= new TGraph();
    gMult[i]->SetMarkerStyle(20);
    gMult[i]->SetMarkerSize(0.7);
    gMult[i]->SetLineColor(colors[i]);
    gMult[i]->SetMarkerColor(colors[i]);
    gMult[i]->SetTitle(";PMT [#];multiplicity/multiplicity_{1000V}");
    hvid[hvi[i]]=i;
  }

  for(int i=0; i<glx_npmt+1; i++){
    gMultP[i]= new TGraph();
    gMultP[i]->SetMarkerStyle(20);
    gMultP[i]->SetMarkerSize(0.7);
    // gMultP[i]->SetLineColor(colors[i]);
    // gMultP[i]->SetMarkerColor(colors[i]);
    gMultP[i]->SetTitle(";HV [V];multiplicity [#]");
  }

  double mult[108]={0};
  for(int i=0; i<ch.GetEntries(); i++){
    ch.GetEvent(i);
    if(hv==1000) mult[pmt]=nph;
  }

  int id[6]={0};
  int ind[glx_npmt+1]={0};
  for(int i=0; i<ch.GetEntries(); i++){
    ch.GetEvent(i);
    gMultP[pmt]->SetPoint(ind[pmt]++,hv,nph);
    if(mult[pmt]>0) gMult[hvid[hv]]->SetPoint(id[hvid[hv]]++,pmt,nph/(double)mult[pmt]);    
  }

  glx_canvasAdd("gMultP",800,500);
  TLegend *le = new TLegend(0.15,0.7,0.8,0.9);
  le->SetFillColor(0);
  le->SetFillStyle(0);
  le->SetBorderSize(0);
  le->SetFillStyle(0);
  le->SetNColumns(6);

  for(int i=0; i<glx_npmt; i++){
    gMultP[i]->Sort();
    gMultP[i]->Draw(i==0 ? "apl" : "pl same");
    gMultP[i]->GetYaxis()->SetRangeUser(0,5);
    gMultP[i]->SetLineColor(kGray);
    gMultP[i]->SetMarkerColor(kGray+1);
    if(i==17) gMultP[i]->SetLineColor(kBlue+1);
    if(i==53) gMultP[i]->SetLineColor(kBlue+1);
    if(i==89) gMultP[i]->SetLineColor(kBlue+1);
    if(i==16) gMultP[i]->SetLineColor(kRed+1);
    if(i==31) gMultP[i]->SetLineColor(kRed+1);
    if(i==84) gMultP[i]->SetLineColor(kRed+1);

    if(i==17) gMultP[i]->SetMarkerColor(kBlue+2);
    if(i==53) gMultP[i]->SetMarkerColor(kBlue+2);
    if(i==89) gMultP[i]->SetMarkerColor(kBlue+2);
    if(i==16) gMultP[i]->SetMarkerColor(kRed+2);
    if(i==31) gMultP[i]->SetMarkerColor(kRed+2);
    if(i==84) gMultP[i]->SetMarkerColor(kRed+2);
    
    if(i==17) gMultP[i]->SetLineWidth(2);
    if(i==53) gMultP[i]->SetLineWidth(2);
    if(i==89) gMultP[i]->SetLineWidth(2);
    if(i==16) gMultP[i]->SetLineWidth(2);
    if(i==31) gMultP[i]->SetLineWidth(2);
    if(i==84) gMultP[i]->SetLineWidth(2);

    // le->AddEntry(gMultP[i],Form("PMT %d ",i),"lp");
    // le->Draw();    
  }

  gMultP[17]->Draw("pl same");
  gMultP[53]->Draw("pl same");
  gMultP[89]->Draw("pl same");
  gMultP[16]->Draw("pl same");
  gMultP[31]->Draw("pl same");
  gMultP[84]->Draw("pl same");
  gPad->SetLogy();

  
  TLegend *leg = new TLegend(0.15,0.7,0.7,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  
  glx_canvasAdd("gMult",800,500);
  for(int i=0; i<6; i++){
    gMult[i]->Sort();
    gMult[i]->Draw(i==0 ? "apl" : "pl same");
    gMult[i]->GetXaxis()->SetRangeUser(200,1100);
    gMult[i]->GetYaxis()->SetRangeUser(-0.5,2);
    leg->AddEntry(gMult[i],Form("%d V",hvi[i]),"lp");
    leg->Draw();    
  }  
  
  glx_canvasSave(0);
}
