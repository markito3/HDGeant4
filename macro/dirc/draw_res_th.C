#include "glxtools.C"

void draw_res_th(TString in="res_thscan_03_02.root"){

  TString nid="mult_dac_letimecut_per1k_";
  glx_savepath = "data/draw_res_th_03_02";
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

  const int max=7;//14;
  TGraph *gMult[max], *gMultP[glx_npmt+1];
  //int hvi[]={230,235,240,250,300,350,400,450,500,550,600,650,750,1000};
  int hvi[]={25,50,100,200,300,400,500};
  int colors[]={1,kGreen+1,kRed+2,kRed,kBlue,kMagenta+2,kBlue+3,kOrange,kGray+2,kOrange+2,kGreen+3,kGray+4,kRed+4,kGreen+5,kYellow+2};
  int hvid[2000]={0};
  for(int i=0; i<max; i++){
    gMult[i]= new TGraph();
    gMult[i]->SetMarkerStyle(20);
    gMult[i]->SetMarkerSize(0.7);
    gMult[i]->SetLineColor(colors[i]);
    gMult[i]->SetMarkerColor(colors[i]);
    gMult[i]->SetTitle(";PMT [#];multiplicity/multiplicity_{dac=100}");
    hvid[hvi[i]]=i;
  }

  for(int i=0; i<glx_npmt+1; i++){
    gMultP[i]= new TGraph();
    gMultP[i]->SetName(Form("pmt_%d",i));
    gMultP[i]->SetMarkerStyle(20);
    gMultP[i]->SetMarkerSize(0.7);
    //gMultP[i]->SetTitle(Form("PMT# %d;DAC threshold; multiplicity_{scaled to 4000 @ dac=100} [#]",i));
    gMultP[i]->SetTitle(Form("PMT# %d;DAC threshold; multiplicity_{per 1k} [#]",i));
  }

  TGaxis::SetMaxDigits(3);
    
  double mult[108]={0};
  for(int i=0; i<ch.GetEntries(); i++){
    ch.GetEvent(i);
    if(dac==100) mult[pmt]=nph;
  }

  int id[max]={0};
  int ind[glx_npmt+1]={0};
  for(int i=0; i<ch.GetEntries(); i++){
    ch.GetEvent(i);
    gMultP[pmt]->SetPoint(ind[pmt]++,dac,nph);
    if(pmt==59) std::cout<<"nph "<<nph<<" "<< mult[pmt]<<" "<<dac<<std::endl;
    
    if(mult[pmt]>0) gMult[hvid[dac]]->SetPoint(id[hvid[dac]]++,pmt,nph/(double)mult[pmt]);
  }

  for(int i=0; i<glx_npmt; i++){
    glx_canvasAdd(nid+Form("all_%d",i),800,500);
    gMultP[i]->Sort();
    gMultP[i]->Draw("apl");
    glx_canvasSave(0);
  }
  
  if(0){
    gStyle->SetOptTitle(0);
    glx_canvasAdd(nid+"p_scaled",1200,800);
    TLegend *le = new TLegend(0.15,0.7,0.8,0.9);
    le->SetFillColor(0);
    le->SetFillStyle(0);
    le->SetBorderSize(0);
    le->SetFillStyle(0);
    le->SetNColumns(6);

    for(int i=0; i<glx_npmt; i++){
      double shift,x,y;
      for(int p=0; p<gMultP[i]->GetN(); p++){
	gMultP[i]->GetPoint(p,x,y);
	if(x==100) shift=4000-y;
      }  
      for(int p=0; p<gMultP[i]->GetN(); p++){
	gMultP[i]->GetPoint(p,x,y);
	gMultP[i]->SetPoint(p,x,y+shift);
      }

      gMultP[i]->Draw(i==0 ? "apl" : "pl same");
      gMultP[i]->GetYaxis()->SetRangeUser(0,7000);
      gMultP[i]->GetXaxis()->SetRangeUser(20,600);    
    }
    //gPad->SetLogy();    
  }

  int lis[]={57,38,62,56,39,75,63,74};
  int ll=0;

  TLegend *la = new TLegend(0.25,0.7,0.8,0.9);
  la->SetFillColor(0);
  la->SetFillStyle(0);
  la->SetBorderSize(0);
  la->SetFillStyle(0);
  la->SetNColumns(6);

  for(auto i: lis){
    gMultP[i]->SetLineColor(colors[1+ll]);
    gMultP[i]->SetMarkerColor(colors[1+ll++]);
    gMultP[i]->GetYaxis()->SetRangeUser(0,0.5);
    gMultP[i]->Draw(i==13? "apl": "pl same");
    la->AddEntry(gMultP[i],Form("%d",i),"lp");
    la->Draw();
  }

  //gPad->SetLogy();
  
  glx_canvasAdd(nid+"div",800,500);
  
  TLegend *l = new TLegend(0.25,0.7,0.8,0.9);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetNColumns(6);

  int list[]={13,17,31,66,84,101};
  int ii=0;
  for(auto i: list){
    // gMultP[i]->SetLineColor(colors[ii]);
    // gMultP[i]->SetMarkerColor(colors[ii++]);
    gMultP[i]->GetYaxis()->SetRangeUser(0,0.5);
    gMultP[i]->Draw(i==13? "apl": "pl same");
    l->AddEntry(gMultP[i],Form("%d",i),"lp");
    l->Draw();
  }
  
  TLegend *leg = new TLegend(0.15,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(4);

  glx_canvasAdd(nid+"examples",800,500);
  for(int i=0; i<max; i++){
    gMult[i]->Sort();
    gMult[i]->Draw(i==0 ? "apl" : "pl same");
    gMult[i]->GetXaxis()->SetRangeUser(200,800);
    gMult[i]->GetYaxis()->SetRangeUser(0,2.9);
    leg->AddEntry(gMult[i],Form("%d",hvi[i]),"lp");
    leg->Draw();
  }

  glx_canvasSave(0);
  
}
