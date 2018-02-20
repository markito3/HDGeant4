#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "../../../../sim-recon/master/src/plugins/Analysis/lut_dirc/DrcLutNode.h"
#include "glxtools.C"

//void reco_lut(TString infile="hit_k_pi.root",TString inlut="lut_50M_cs_avr.root"){
void reco_lut(TString infile="hit_k_pi.root",TString inlut="lut_50M_avr.root"){
  if(!glx_initc(infile,1,"data/reco_lut")) return;

  TFile *fLut = new TFile(inlut);
  TTree *tLut=(TTree *) fLut->Get("lut_dirc") ;
  TClonesArray *cLut = new TClonesArray("DrcLutNode");
  tLut->SetBranchAddress("LUT_1",&cLut); 
  tLut->GetEntry(0);
  const int nodes = 30000;
  DrcLutNode *lutNode[nodes];
  for(int i=0; i<nodes; i++) lutNode[i] = (DrcLutNode*) cLut->At(i);
  TGaxis::SetMaxDigits(4);
  
  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  TVector3 fnZ1 = TVector3( 0,0,1);
  double radiatorL = 4*1225; 
  double barend = 2938; //3036.05;
  double minChangle = 0.6;
  double maxChangle = 0.9;
  double sum1,sum2,noise = 0.3;
  double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
  double pathid,evtime,luttheta,tangle,lenz;
  TVector3 momInBar,dir,dird;

  TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
  TSpectrum *spect = new TSpectrum(10);   
  TH1F *hAngle[5], *hLnDiff[5];
  TF1  *fAngle[5];
  double mAngle[5];
  TH1F *hDiff = new TH1F("hDiff",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
  TH1F *hTime = new TH1F("hTime",";propagation time [ns];entries [#]",   1000,0,100);
  TH1F *hCalc = new TH1F("hCalc",";calculated time [ns];entries [#]",   1000,0,100);
  TH1F *hNph = new TH1F("hNph",";detected photons [#];entries [#]", 150,0,150);
  TH1F *hNphC = new TH1F("hNphC",";detected photons [#];entries [#]", 150,0,150);
  
  double momentum=4;
  for(int i=0; i<5; i++){
    hAngle[i]= new TH1F(Form("hAngle_%d",i),  "chrenkov angle;#theta_{C} [rad];entries/N_{max} [#]", 250,0.6,1);
    mAngle[i] = acos(sqrt(momentum * momentum + glx_mass[i]*glx_mass[i])/momentum/1.4738);  //1.4738
    fAngle[i] = new TF1(Form("fAngle_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
    fAngle[i]->SetParameter(0,1);        // const
    fAngle[i]->SetParameter(1,mAngle[i]);// mean
    fAngle[i]->SetParameter(2,0.007);    // sigma
    hAngle[i]->SetMarkerStyle(20);
    hAngle[i]->SetMarkerSize(0.8);
 
    hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(p) - ln L(#pi);entries [#]",200,-30,30);
  }
  hAngle[2]->SetLineColor(4);
  hAngle[3]->SetLineColor(2);
  hAngle[2]->SetMarkerColor(kBlue+1);
  hAngle[3]->SetMarkerColor(kRed+1);
  fAngle[2]->SetLineColor(4);
  fAngle[3]->SetLineColor(2);

  hLnDiff[2]->SetLineColor(4);
  hLnDiff[3]->SetLineColor(2);
  
  DrcHit hit;
  for (int e = 0; e < glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    
    for (int t = 0; t < glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,10);
      momInBar = glx_event->GetMomentum();
      int pdgId = glx_findPdgId(glx_event->GetPdg());
      if(glx_event->GetParent()>0) continue;
									    
      if(glx_event->GetPosition().Y()<0) lenz = fabs(barend+glx_event->GetPosition().X()*10);
      else lenz =fabs(glx_event->GetPosition().X()*10-barend);

      int ncount=0;
      sum1=0;
      sum2=0;
      int nph=0;
      hNphC->Fill(glx_event->GetHitSize());
      
      for(int h = 0; h < glx_event->GetHitSize(); h++){
    	hit = glx_event->GetHit(h);
    	int pmt = hit.GetPmtId();
    	int pix = hit.GetPixelId();
    	double hitTime = hit.GetLeadTime()-glx_event->GetTime();
	TVector3 gpos = hit.GetPosition();

	//	if(hitTime<40) continue;
	
	bool reflected = hitTime>40;
	if(reflected) lenz = 2*radiatorL - lenz;
	int sensorId = 100*pmt + pix;
	bool isGood(false);
	
	for(int i = 0; i < lutNode[sensorId]->Entries(); i++){
	  dird   = lutNode[sensorId]->GetEntry(i);
	  evtime = lutNode[sensorId]->GetTime(i);
	  pathid = lutNode[sensorId]->GetPathId(i);
	  bool samepath(false);
	  if(fabs(pathid-hit.GetPathId())<0.0001) samepath=true;	  	  
	  if(!samepath) continue;
	  
	  for(int u = 0; u < 4; u++){
	    if(u == 0) dir = dird;
	    if(u == 1) dir.SetXYZ( dird.X(), -dird.Y(),  dird.Z());
	    if(u == 2) dir.SetXYZ( dird.X(),  dird.Y(), -dird.Z());
	    if(u == 3) dir.SetXYZ( dird.X(),-dird.Y(), -dird.Z());
	    if(reflected) dir.SetXYZ( -dir.X(), dir.Y(), dir.Z());
	    if(dir.Angle(fnY1) < criticalAngle || dir.Angle(fnZ1) < criticalAngle) continue;

	    TVector3 vt = dir;
	    vt.RotateY(TMath::PiOver2());
	    luttheta = vt.Theta();	  
	    if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
	    
	    tangle = momInBar.Angle(dir);//-0.002; //correction
	    double bartime = fabs(lenz/cos(luttheta)/203.76); //198 //203.767 for 1.47125
	    double totalTime = bartime+evtime;

	    if(fabs(totalTime-hitTime)>10) continue;
	    if(fabs(tangle-0.82)>0.3) continue;
	      
	    isGood=true;
	    ncount++;
	    //	    if(ncount>1) continue;
	   
	    hAngle[pdgId]->Fill(tangle);
	    hDiff->Fill(totalTime-hitTime);
	    hTime->Fill(hitTime);
	    hCalc->Fill(totalTime);

	    sum1 += TMath::Log(fAngle[2]->Eval(tangle)+0.0001);
	    sum2 += TMath::Log(fAngle[3]->Eval(tangle)+0.0001);	    
	    
	  }
	}
	if(!isGood) {
	  nph++;
	  if(pmt<108) glx_hdigi[pmt]->Fill(pix%8, 7-pix/8);
	}
      }
      hNph->Fill(nph);
      
      // std::cout<<"ncount "<<ncount<<std::endl;            
      double sum = sum1-sum2;
      hLnDiff[pdgId]->Fill(sum);
      
      // TCanvas *c = new TCanvas("c","c",800,500);
      // hAngle[3]->Draw();
      // c->Update();
      // c->WaitPrimitive();
      // hAngle[3]->Reset();
    }
  }

  glx_drawDigi("m,p,v\n",0);
  glx_canvasAdd(glx_cdigi);

  glx_canvasAdd("hAngle",800,500);

  hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
  hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
  
  double cherenkovreco[5],spr[5];

  for(auto i=0; i<5; i++){
    if(hAngle[i]->GetEntries()<100) continue;
    
    int nfound = spect->Search(hAngle[i],1,"",0.9);
    if(nfound>0) cherenkovreco[i] = spect->GetPositionX()[0];
    else cherenkovreco[i] =  hAngle[i]->GetXaxis()->GetBinCenter(hAngle[i]->GetMaximumBin());
    if(cherenkovreco[i]>0.85) cherenkovreco[i]=0.82;
    
    if(i==2)  fit->SetLineColor(kBlue);
    if(i==3)  fit->SetLineColor(kRed);
    fit->SetParameters(100,cherenkovreco[i],0.010,10);
    fit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
    fit->SetParLimits(0,0.1,1E6);
    fit->SetParLimits(1,cherenkovreco[i]-0.04,cherenkovreco[i]+0.04);
    fit->SetParLimits(2,0.005,0.030); // width
    hAngle[i]->Fit("fgaus","I","",cherenkovreco[i]-0.08,cherenkovreco[i]+0.08);
    hAngle[i]->Fit("fgaus","M","",cherenkovreco[i]-0.08,cherenkovreco[i]+0.08);

    cherenkovreco[i] = fit->GetParameter(1);
    spr[i] = fit->GetParameter(2);
    
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  TF1 *ff;
  double sep, m1,m2,s1,s2; 
  if(hLnDiff[3]->GetEntries()>10){
    hLnDiff[3]->Fit("gaus","S");
    ff = hLnDiff[3]->GetFunction("gaus");
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
  }
  if(hLnDiff[2]->GetEntries()>10){
    hLnDiff[2]->Fit("gaus","S");
    ff = hLnDiff[2]->GetFunction("gaus");
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
  }
  sep = (fabs(m2-m1))/(0.5*(s1+s2));
  std::cout<<"separation "<< sep <<std::endl;

  hAngle[2]->GetXaxis()->SetRangeUser(0.7,0.9);
  hAngle[2]->GetYaxis()->SetRangeUser(0,1.2);
  hAngle[2]->Draw();
  hAngle[3]->Draw("same");
  fAngle[3]->Draw("same");
  fAngle[2]->Draw("same");

  TLine *line = new TLine(0,0,0,1000);
  line->SetX1(mAngle[3]);
  line->SetX2(mAngle[3]);
  line->SetY1(0);
  line->SetY2(1.2);
  line->SetLineColor(kRed);
  line->Draw();

  TLine *line2 = new TLine(0,0,0,1000);
  line2->SetX1(mAngle[2]);
  line2->SetX2(mAngle[2]);
  line2->SetY1(0);
  line2->SetY2(1.2);
  line2->SetLineColor(kBlue);
  line2->Draw();

  TLegend *leg = new TLegend(0.15,0.5,0.45,0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hAngle[2],Form("#theta_{c}^{#pi} = %2.4f rad",cherenkovreco[2]),"");
  leg->AddEntry(hAngle[3],Form("#theta_{c}^{K} = %2.4f rad",cherenkovreco[3]),"");
  leg->AddEntry(hAngle[2],Form("#sigma_{c}^{#pi} = %2.2f mrad",spr[2]*1000),"");
  leg->AddEntry(hAngle[3],Form("#sigma_{c}^{K} = %2.2f mrad",spr[3]*1000),"");
  leg->Draw();
  
  // fAngle[2]->Draw("same");
  // fAngle[3]->Draw("same");
  
  glx_canvasAdd("hTime",800,400);
  
  hTime->Draw();
  hCalc->SetLineColor(2);
  hCalc->Draw("same");
  TLegend *leg1 = new TLegend(0.5,0.6,0.85,0.80);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hTime,"measured in geant","lp");
  leg1->AddEntry(hCalc,"calculated","lp");
  leg1->Draw();


  glx_canvasAdd("hDiff",800,400);
  hDiff->Draw();


  
  glx_canvasAdd("hLnDiff",800,400);
  hLnDiff[2]->Draw();
  hLnDiff[3]->Draw("same");

  glx_canvasAdd("hNph",800,400);
  hNph->Draw();
  hNphC->SetLineColor(kRed);
  hNphC->Draw("same");

  //glx_canvasSave(1,0);
  
}
