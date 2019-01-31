#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"
#include <TVirtualFitter.h>
#include <TKey.h>

TLine *gLine = new TLine(0,0,0,1000);

void recoPdf(TString path="data/pdf_test.root", TString data="data/test.root",  Double_t sigma=1000,Bool_t debug=0, Double_t r1=0, Double_t r2=0, Int_t nforpdf=0){
// void recoPdf(TString path="momscan/pdf_p0.5.root", TString data="momscan/data_p0.5.root",  Double_t sigma=200,Bool_t debug=false, Double_t r1=0, Double_t r2=0, Int_t nforpdf=0){
TCanvas *cc;
  if(debug) cc = new TCanvas("cc","cc",800,400);
  
  TH1F *hpdff[glx_npixtot],*hpdfs[glx_npixtot], *hl[5],*hnph[5],*hll[5][5], *hnphof, *hnphos;
  TGraph *gpdff[glx_npixtot],*gpdfs[glx_npixtot];
  TF1 *gnphof, *gnphos;
  for(Int_t i=0; i<5; i++){
    hl[i] = new TH1F(Form("hl_%d",i),"pdf; LE time [ns]; entries [#]", 1000,0,50);
    hnph[i] = new TH1F(Form("hnph_%d",i),"; detected photons [#]; entries [#]", 160,0,160);
    hll[i][0] = new TH1F(Form("hll_%d_mom0",i),"hll; ln L(K) - ln L(#pi); entries [#]",200,-300,300);
    hll[i][1] = new TH1F(Form("hll_%d_mom1",i),"hll; ln L(K) - ln L(#pi); entries [#]",200,-200,200);
    hll[i][2] = new TH1F(Form("hll_%d_mom2",i),"hll; ln L(K) - ln L(#pi); entries [#]",200,-200,200);
    hll[i][3] = new TH1F(Form("hll_%d_mom3",i),"hll; ln L(K) - ln L(#pi); entries [#]",200,-100,100);
    hll[i][4] = new TH1F(Form("hll_%d_mom4",i),"hll; ln L(K) - ln L(#pi); entries [#]",200,-100,100);
  }
  TH1F *hnphf =  new TH1F("hnphf","hnphf",200,0,200);
  TH1F *hnphs =  new TH1F("hnphs","hnphs",200,0,200);

  TFile f(path);

  TRandom rand;

  std::cout<<"total number of pixels = "<<glx_npixtot<<std::endl;
  TF1* fitf = new TF1("fitf","[0]*TMath::Poisson(x,[1])");
  TF1* fits = new TF1("fits","[0]*TMath::Poisson(x,[1])");
  glx_canvasAdd("NphoN",500,500);
  hnphof = (TH1F*)f.Get(Form("npho_pion"));
  hnphos = (TH1F*)f.Get(Form("npho_kaon"));
  hnphof->SetDirectory(0);
  hnphos->SetDirectory(0);
  hnphof->SetLineColor(4);
  hnphof->SetLineWidth(2);
  hnphof->Scale(1./hnphof->GetMaximum());
  hnphos->SetLineColor(2);
  hnphos->SetLineWidth(2);
  hnphos->Scale(1./hnphos->GetMaximum());
  fitf->SetParameter(0, 20);
  fits->SetParameter(0, 5);
  fitf->SetParameter(1, 50);
  fits->SetParameter(1, 2);
  hnphos->GetXaxis()->SetRangeUser(0.,100.);
  hnphos->Draw("h");
  hnphof->Draw("h same");
 
  //  cout<<"nbins = "<<hnphof->GetXaxis()->GetNbins()<<", x max = "<<hnphof->GetXaxis()->GetXmax()<<", x min = "<<hnphof->GetXaxis()->GetXmin()<<", cont = "<<hdphof->GetBinContent( () )<<endl;
  
    for(Int_t i=0; i<glx_npixtot; i++){
    hpdff[i] = (TH1F*)f.Get(Form("h_%d_pion",i));
    hpdfs[i] = (TH1F*)f.Get(Form("h_%d_kaon",i));
    if(hpdff[i]){
      if(sigma > 0.){
	hpdff[i]->Rebin((Int_t)(sigma/50.+0.1));
      }
      hpdff[i]->SetLineColor(2);
      gpdff[i] = new TGraph(hpdff[i]);
    }
    if(hpdfs[i]){
      if(sigma > 0.){
	hpdfs[i]->Rebin((Int_t)(sigma/50.+0.1));
      }
      hpdfs[i]->SetLineColor(4);
      gpdfs[i] = new TGraph(hpdfs[i]);
    }
  }
 
  TFile f1(data);
  if(!glx_initc(data,1,"data/dataForPdf")) return;
  Double_t theta(0);
  TVirtualFitter *fitter;
  Double_t nphf, nphs,time, time0;//,timeres(-1);
  DrcHit fHit;
  Int_t totalf(0),totals(0),mcp,pix,ch;
  std::cout<<"entries = "<<glx_ch->GetEntries()<<std::endl;
  Int_t firstCount = 0;
  Int_t intmom = 0;
  Double_t Nf, Ns;
  Int_t nparticle[5]={0};
  
  for (Int_t ievent=0; ievent<glx_ch->GetEntries(); ievent++){
    glx_ch->GetEntry(ievent);
    if(nparticle[2] > 10000 && nparticle[3] > 10000)break;
    for(Int_t t=0; t < glx_events->GetEntriesFast(); t++){
      glx_nextEventc(ievent, t, 10);
      if(glx_event->GetParent()>0)continue;
      Double_t aminf,amins, sum(0),sumf(0),sums(0);
      Int_t nGoodHits(0), nHits =glx_event->GetHitSize();
      if (fabs(glx_event->GetPdg()) != 321 && fabs(glx_event->GetPdg()) != 211) continue;
      Int_t pid = glx_findPdgId(glx_event->GetPdg());
      if(nparticle[pid]++ > 10000)continue;

      time0 = glx_event->GetTime();
      //if(debug)std::cout<<"===================== event === "<< ievent <<", pid - "<<pid<<std::endl;
      intmom = glx_event->GetMomentum().Mag() + 0.1;
      
      Int_t mult[glx_npixtot];
      memset(mult, 0, sizeof(mult));
      for(Int_t i=0; i<nHits; i++){
	fHit = glx_event->GetHit(i);
	mcp = fHit.GetPmtId();
	pix = fHit.GetPixelId();
	ch = glx_getChNum(mcp, pix);
	
	//	if(ch != 2304) continue; /// !!!!!!!!
	time = fHit.GetLeadTime() - time0;
      	time+= rand.Gaus(0,sigma/1000.); // 1 ns, 600 ps timing resolution
	//		std::cout<<"=== event === "<< ievent <<", pid - "<<pid<<", ch = "<<ch<<", pix = "<<pix<<", mcp = "<<mcp<<", time = "<<time<<std::endl;	
	nGoodHits++;
	
       	if(!gpdff[ch] || !gpdfs[ch]) continue;
	aminf = gpdff[ch]->Eval(time);
	amins = gpdfs[ch]->Eval(time);
	//	std::cout<<"aminf = "<<aminf<<", amins = "<<amins<<std::endl;
      	if(amins < 1e-20 || aminf < 1e-20) continue;
        if(debug){
	  firstCount++;
	  TString x=(aminf>amins)? " <====== KAON" : "";
          std::cout<<Form("f %1.6f s %1.6f mcp %d pix %d   pid %d",aminf,amins,mcp,pix,pid)<<"  "<<x <<std::endl;
	  cc->cd();
	  std::cout<<"hist 1 "<<hpdff[ch]<<", hist2 "<<hpdfs[ch]<<std::endl;
	  if(hpdff[ch] && hpdfs[ch]) std::cout<<"exist"<<std::endl;
	  if(!hpdff[ch] || !hpdfs[ch]) std::cout<<"NOOOOO"<<std::endl;
       	  if(firstCount == 1) glx_normalize(hpdff[ch],hpdfs[ch]); //!!!!
	  hpdff[ch]->SetLineColor(2);
	  hpdfs[ch]->SetLineColor(4);
	  hpdff[ch]->Draw();
	  hpdff[ch]->SetTitle(Form("mcp=%d  pix=%d",mcp,pix));
	  hpdff[ch]->GetXaxis()->SetTitle("LE time [ns]");
	  hpdff[ch]->GetYaxis()->SetTitle("PDF value");
	  hpdfs[ch]->Draw("same");
	  gpdff[ch]->Draw("PL same");
	  gpdfs[ch]->Draw("PL same");
	  cc->Update();
	  gLine->SetLineWidth(2);
	  gLine->SetX1(time);
	  gLine->SetX2(time);
	  // gLine->SetY1(cc->GetUymin());
	  //gLine->SetY2(cc->GetUymax());
	  hpdff[ch]->GetXaxis()->SetRangeUser(20.,60.);
	  gLine->Draw();
	  cc->Update();
	  cc->WaitPrimitive();
	}
	
	Double_t noise = 1e-6; //1e-7; // nHits //1e-5
	
	sumf+=TMath::Log((aminf+noise));
	sums+=TMath::Log((amins+noise));

	hl[pid]->Fill(time);
	
      }
      //            cout<<"n good hits = "<<nGoodHits<<endl; 
      if(nGoodHits<1) continue; // !!!!!!! was 5
      hnph[pid]->Fill(nGoodHits);
      Nf = hnphof->GetBinContent( hnphof->FindFixBin(nHits) );
      Ns = hnphos->GetBinContent( hnphos->FindFixBin(nHits) );
      //      cout<<"nHits = "<<nHits<<", Nf = "<<Nf<<", Ns = "<<Ns<<", pdg = "<<pid<<endl;
      sum =  sumf - sums + ( TMath::Log(Nf) - TMath::Log(Ns) );
      hll[pid][intmom]->Fill(sum);
      //std::cout<<"sumf "<<sumf<<", sums = "<<sums<<". sum = "<<sum<<endl;//", nHits = "<<nHits<<", Nf = "<<Nf<<", Ns = "<<Ns<<", logf = "<<TMath::Log(Nf)<<", logs = "<<TMath::Log(Ns)<<std::endl;      
 
    }
  }
  
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);

  glx_canvasAdd("ll",800,400);

  TF1 *ff;
  Double_t sep(0),esep(0),m1,m2,s1,s2,dm1,dm2,ds1,ds2;
  if(hll[3][intmom]->GetEntries()>10 && hll[2][intmom]->GetEntries()>10){
    hll[3][intmom]->Fit("gaus","Sq","");
    ff = hll[3][intmom]->GetFunction("gaus");
    ff->SetLineColor(1);
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
    dm1=ff->GetParError(1);
    ds1=ff->GetParError(2);

    hll[2][intmom]->Fit("gaus","Sq");
    ff = hll[2][intmom]->GetFunction("gaus");
    ff->SetLineColor(1);
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
    dm2=ff->GetParError(1);
    ds2=ff->GetParError(2);

    sep = (fabs(m1-m2))/(0.5*(s1+s2));

    Double_t e1,e2,e3,e4;
    e1=2/(s1+s2)*dm1;
    e2=2/(s1+s2)*dm2;
    e3=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds1;
    e4=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds2;
    esep=sqrt(e1*e1+e2*e2+e3*e3+e4*e4);
  }
  
  hll[3][intmom]->SetTitle(Form("separation = %1.2f",sep));
  hll[3][intmom]->SetLineWidth(2);
  hll[3][intmom]->SetLineColor(2);
  hll[3][intmom]->Draw();
  hll[2][intmom]->SetLineWidth(2);
  hll[2][intmom]->SetLineColor(4);
  hll[2][intmom]->Draw("same");

  TLegend *leg = new TLegend(0.65,0.65,0.83,0.87);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hll[2][intmom],"pions ","lp");
  leg->AddEntry(hll[3][intmom],"kaons","lp");
  leg->Draw();

  glx_canvasAdd("hnph",800,500);

  // glx_normalize(hnph[3],hnph[2]);
  nphf = hnph[2]->GetMaximum();
  nphs = hnph[3]->GetMaximum();
  hnph[3]->SetLineColor(2);
  hnph[3]->Draw();
  hnph[2]->SetLineColor(4);
  hnph[2]->Draw("same");

  TLegend *leg1 = new TLegend(0.65,0.65,0.83,0.87);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hnph[2],"pions ","lp");
  leg1->AddEntry(hnph[3],"kaons","lp");
  leg1->Draw();
  
  std::cout<<dm1<<" "<<dm2<<" "<<ds1 <<" "<<ds2<<std::endl;
  std::cout<<path<<" separation "<< sep <<" +/- "<<esep <<std::endl;
  std::cout<<"entries:  pi "<<hll[2][intmom]->GetEntries()<<" p "<<hll[3][intmom]->GetEntries() <<std::endl;

  TString output=data;
  TString add(Form("_%1.1f",r2));
  output.ReplaceAll(".root", add+"_res.root");
  TFile fc(output,"recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("esep",&esep,"esep/D");
  tc->Branch("sigma",&sigma,"sigma/D");
  tc->Branch("nphf",&nphf,"nphf/D");
  tc->Branch("nphs",&nphs,"nphs/D");
  tc->Branch("r1",&r1,"r1/D");
  tc->Branch("r2",&r2,"r2/D");
  tc->Fill();
  tc->Write();
  fc.Write();
  fc.Close();
    
}
