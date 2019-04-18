#define glx__sim
#include "../../../halld_recon/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "../../../halld_recon/src/plugins/Analysis/lut_dirc/DrcLutNode.h"
#include "glxtools.C"

void reco_lut(TString infile="vol/tree_060772.root",TString inlut="lut/lut_12/lut_all_avr.root",int xbar=-1, int ybar=-1, double moms=4){

  if(!glx_initc(infile,1,"data/reco_lut_tt_sim")) return;
  const int nodes = glx_maxch;
  const int luts = 24;
  
  TFile *fLut = new TFile(inlut);
  TTree *tLut=(TTree *) fLut->Get("lut_dirc") ;
  TClonesArray* cLut[luts];
  for(int l=0; l<luts; l++){
    cLut[l] = new TClonesArray("DrcLutNode");
    tLut->SetBranchAddress(Form("LUT_%d",l),&cLut[l]); 
  }  
  tLut->GetEntry(0);

  DrcLutNode *lutNode[luts][nodes];
  for(int l=0; l<luts; l++){
    for(int i=0; i<nodes; i++) lutNode[l][i] = (DrcLutNode*) cLut[l]->At(i);
  }
  TGaxis::SetMaxDigits(4);
  
  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  TVector3 fnZ1 = TVector3( 0,0,1);
  double radiatorL = 489.712; //4*122.5;
  double barend = -294.022; // 4*1225-1960; -294.022

  double minChangle = 0.6;
  double maxChangle = 0.9;
  double sum1,sum2,noise = 0.5;
  double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
  double evtime,luttheta,tangle,lenz;
  int64_t pathid;
  TVector3 posInBar,momInBar,dir,dird;
  double cherenkovreco[5],spr[5];

  TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
  TSpectrum *spect = new TSpectrum(10);   
  TH1F *hAngle[5], *hLnDiff[5], *hNph[5];
  TF1  *fAngle[5];
  double mAngle[5];
  TH1F *hDiff = new TH1F("hDiff",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
  TH1F *hDiffT = new TH1F("hDiffT",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
  TH1F *hDiffD = new TH1F("hDiffD",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
  TH1F *hDiffR = new TH1F("hDiffR",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
  TH1F *hTime = new TH1F("hTime",";propagation time [ns];entries [#]",   1000,0,200);
  TH1F *hCalc = new TH1F("hCalc",";calculated time [ns];entries [#]",   1000,0,200);
  TH1F *hNphC = new TH1F("hNphC",";detected photons [#];entries [#]", 150,0,150);
  hDiff->SetMinimum(0);

  
  TGaxis::SetMaxDigits(3);
  
  double sigma[]={0.01,0.01,0.01,0.010,0.01,0.01};
  for(int i=0; i<5; i++){
    double momentum=4;
    hAngle[i]= new TH1F(Form("hAngle_%d",i),  "cherenkov angle;#theta_{C} [rad];entries/N_{max} [#]", 250,0.6,1);
    hNph[i] = new TH1F(Form("hNph_%d",i),";detected photons [#];entries [#]", 150,0,150);
    mAngle[i] = acos(sqrt(momentum * momentum + glx_mass[i]*glx_mass[i])/momentum/1.473);  //1.4738
    fAngle[i] = new TF1(Form("fAngle_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
    fAngle[i]->SetParameter(0,1);        // const
    fAngle[i]->SetParameter(1,mAngle[i]);// mean
    fAngle[i]->SetParameter(2,sigma[i]);    // sigma
    hAngle[i]->SetMarkerStyle(20);
    hAngle[i]->SetMarkerSize(0.8);
    hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",120,-60,60);
  }
  
  hAngle[2]->SetLineColor(4);
  hAngle[3]->SetLineColor(2);
  hAngle[2]->SetMarkerColor(kBlue+1);
  hAngle[3]->SetMarkerColor(kRed+1);
  fAngle[2]->SetLineColor(4);
  fAngle[3]->SetLineColor(2);

  hLnDiff[2]->SetLineColor(4);
  hLnDiff[3]->SetLineColor(2);
  int evtcount=0,count[5]={0};
  TCanvas *cc = new TCanvas("cc","cc",800,500);
  TLine *gLine = new TLine();
  // cuts
  double cut_cangle=0.02;
  double cut_tdiffd=4;//3;
  double cut_tdiffr=4;//3.5;

  const int nbins=20;
  DrcHit hit;
  for (int e = 0; e < glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    
    for (int t = 0; t < glx_events->GetEntriesFast(); t++){      
      glx_nextEventc(e,t,1000);
      posInBar = glx_event->GetPosition();
      momInBar = glx_event->GetMomentum();
      double momentum = momInBar.Mag();
      int pdgId = glx_findPdgId(glx_event->GetPdg());
      int bar = glx_event->GetId();
      //if(count[pdgId]>1000) continue;	
      
      // selection
      if(glx_event->GetType()!=2) continue; //1-LED 2-beam 0-rest
      if(momInBar.Mag()<3.9 || momInBar.Mag()>4.15 ) continue;
      //if(momInBar.Mag()<2.0 || momInBar.Mag()>2.5 ) continue;

      int bin = (100+posInBar.X())/200.*nbins;
      // if(bar<0 || bar>=luts || (bar!=ybar && ybar!=-1)) continue;
      // if(bin<0 || bin>nbins || (bin!=xbar && xbar!=-1)) continue; 

      if(bar<0 || bar>=luts || (bar<3 || bar>3)) continue;
      if(bin<0 || bin>nbins || (bin<10 || bin>10)) continue; 

      // if(bar<0 || bar>=luts) continue;
      // if(bin<0 || bin>nbins) continue; 
      
      if(glx_event->GetParent()>0) continue;
      // if(hLnDiff[pdgId]->GetEntries()>200) continue;

      for(int p=0; p<5; p++){
	mAngle[p] = acos(sqrt(momentum * momentum + glx_mass[p]*glx_mass[p])/momentum/1.473);  //1.4738
	fAngle[p]->SetParameter(1,mAngle[p]);// mean
      }
      
      sum1=0;
      sum2=0;
      int nph=0;
      int nphc=0;
      //      hNphC->Fill(glx_event->GetHitSize());

      bool goodevt=0;
      for(int h = 0; h < glx_event->GetHitSize(); h++){
    	hit = glx_event->GetHit(h);
	int ch = hit.GetChannel();		
	int pmt = hit.GetPmtId();
    	int pix = hit.GetPixelId();
	
    	double hitTime = hit.GetLeadTime()-glx_event->GetTime();
	
	if(pmt<=10 || (pmt>=90 && pmt<=96)) continue; // dummy pmts
	if(ch>glx_nch) continue;
	//if(hitTime>40) continue;
	nphc++;
	
	bool reflected = hitTime>40;	
	lenz = fabs(barend-posInBar.X());
	double rlenz = 2*radiatorL - lenz;
	double dlenz = lenz;	
	if(reflected) lenz = 2*radiatorL - lenz;

	bool isGood(false);

	double p1,p2;
	
	for(int i = 0; i < lutNode[bar][ch]->Entries(); i++){
	  dird   = lutNode[bar][ch]->GetEntry(i);
	  evtime = lutNode[bar][ch]->GetTime(i);
	  pathid = lutNode[bar][ch]->GetPathId(i);
	  bool samepath(false);
	  if(fabs(pathid-hit.GetPathId())<0.0001) samepath=true;
	  p1=hit.GetPathId();
	  if(samepath) p2=pathid;	  
	  //if(!samepath) continue;
	  
	  for(int r=0; r<2; r++){
	    if(!reflected && r==1) continue;
		    
	    if(r) lenz = rlenz;	      
	    else lenz = dlenz;
	    
	    for(int u = 0; u < 4; u++){
	      if(u == 0) dir = dird;
	      if(u == 1) dir.SetXYZ( dird.X(),-dird.Y(),  dird.Z());
	      if(u == 2) dir.SetXYZ( dird.X(), dird.Y(), -dird.Z());
	      if(u == 3) dir.SetXYZ( dird.X(),-dird.Y(), -dird.Z());
	      if(r) dir.SetXYZ( -dir.X(), dir.Y(), dir.Z());	   
	      if(dir.Angle(fnY1) < criticalAngle || dir.Angle(fnZ1) < criticalAngle) continue;

	      luttheta = dir.Angle(TVector3(-1,0,0));
	      if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
	      tangle = momInBar.Angle(dir);//-0.004; //correction
	    
	      //double bartime = lenz/cos(luttheta)/20.4; //198 //203.767 for 1.47125
	      //double bartime = lenz/cos(luttheta)/19.8; //198 //203.767 for 1.47125
	      double bartime = lenz/cos(luttheta)/19.6; //203.767 for 1.47125
	      double totalTime = bartime+evtime;
	      // hTime->Fill(hitTime);
	      // hCalc->Fill(totalTime);
	      
	      if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))<cut_cangle){
		hDiff->Fill(totalTime-hitTime);
		//if(samepath)
		  {
		  hDiffT->Fill(totalTime-hitTime);
		  if(r) hDiffR->Fill(totalTime-hitTime);
		  else hDiffD->Fill(totalTime-hitTime);
		}
	      }
	      
	      if(!r && fabs(totalTime-hitTime)>cut_tdiffd) continue;
	      if(r && fabs(totalTime-hitTime) >cut_tdiffr) continue;

	      hAngle[pdgId]->Fill(tangle);
	      if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))>0.05) continue;
	      isGood=true;	  

	      hTime->Fill(hitTime);
	      hCalc->Fill(totalTime);

	      //	      std::cout<<pdgId<<" TMath::Log(fAngle[2]->Eval(tangle)+0.001) "<<TMath::Log(fAngle[2]->Eval(tangle)+noise)<<"    "<< TMath::Log(fAngle[3]->Eval(tangle)+noise)<< " "<< tangle<<std::endl;


	      
	      sum1 += TMath::Log(fAngle[2]->Eval(tangle)+noise);
	      sum2 += TMath::Log(fAngle[3]->Eval(tangle)+noise);


	      if(0){
		TString x=(sum1>sum2)? " <====== PION" : "";
		std::cout<<Form("%1.6f  %1.6f | %1.6f  %1.6f        pid %d",TMath::Log(fAngle[2]->Eval(tangle)+noise),TMath::Log(fAngle[3]->Eval(tangle)+noise), sum1, sum2,pdgId)<<"  " <<std::endl;
	
		cc->cd();	
		fAngle[2]->Draw("");
		fAngle[3]->Draw("same");
	
		cc->Update();
		gLine->SetLineWidth(2);
		gLine->SetX1(tangle);
		gLine->SetX2(tangle);
		gLine->SetY1(cc->GetUymin());
		gLine->SetY2(cc->GetUymax());
		gLine->Draw();
		cc->Update();
		cc->WaitPrimitive();
	      }

	    }
	  }
	}
	
	if(isGood){
	  nph++;
	  if(pmt<108) {
	    glx_hdigi[pmt]->Fill(pix%8, pix/8);
	    goodevt=1;
	  }
	}
      }

      if(goodevt) evtcount++;
      
      if(nph<5) continue;
      hNph[pdgId]->Fill(nph);
      hNphC->Fill(nphc);
      
      double sum = sum1-sum2;
      hLnDiff[pdgId]->Fill(sum);
      
      count[pdgId]++;

      if(0 && pdgId==3){
	//	if(!cc)
	TString x=(sum1>sum2)? " <====== Pion" : "";
	// std::cout<<Form("f %1.6f s %1.6f mcp %d pix %d   pid %d",aminf,amins,mcp,pix  ,prt_particle)<<"  "<<x <<std::endl;

	std::cout<<"PID "<< glx_event->GetPdg() <<" sum1 "<<sum1<<" sum2 "<<sum2<<" sum "<<sum<<" "<<x<<std::endl;
	
	cc->cd();

	if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
	if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
	
	hAngle[2]->Draw("hist");
	hAngle[3]->Draw("hist same");
	fAngle[2]->Draw("same");
	fAngle[3]->Draw("same");

	// hAngle[2]->GetYaxis()->SetRangeUser(0,20);
	// hAngle[3]->GetYaxis()->SetRangeUser(0,20);
	
	cc->Update();
	TLine *line = new TLine(0,0,0,1000);
	line->SetX1(mAngle[3]);
	line->SetX2(mAngle[3]);
	line->SetY1(cc->GetUymin());
	line->SetY2(cc->GetUymax());
	line->SetLineColor(kRed);
	line->Draw();

	TLine *line2 = new TLine(0,0,0,1000);
	line2->SetX1(mAngle[2]);
	line2->SetX2(mAngle[2]);
	line2->SetY1(cc->GetUymin());
	line2->SetY2(cc->GetUymax());
	line2->SetLineColor(kBlue);
	line2->Draw();
	
	cc->Update();
	cc->WaitPrimitive();
      }

      // hAngle[2]->Reset();
      // hAngle[3]->Reset();

    }
  }

  if(evtcount>0){
    for(int i=0; i<glx_nch; i++){
      int pmt=i/64;
      int pix=i%64;
      double rel = glx_hdigi[pmt]->GetBinContent(pix%8+1,pix/8+1)/(double)evtcount;
      glx_hdigi[pmt]->SetBinContent(pix%8+1, pix/8+1,rel);
    }
  }
  
  //TString nid=Form("_%2.2f_%2.2f",theta,phi);
  TString nid=Form("_%d_%d",xbar,ybar);
  glx_drawDigi("m,p,v\n",0);
  glx_cdigi->SetName("hp"+nid);
  glx_canvasAdd(glx_cdigi);

  glx_canvasAdd("hAngle"+nid,800,400);

  if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
  if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
  
  for(int i=0; i<5; i++){
    if(hAngle[i]->GetEntries()<20) continue;
    
    int nfound = spect->Search(hAngle[i],1,"goff",0.9);
    if(nfound>0) cherenkovreco[i] = spect->GetPositionX()[0];
    else cherenkovreco[i] =  hAngle[i]->GetXaxis()->GetBinCenter(hAngle[i]->GetMaximumBin());
    if(cherenkovreco[i]>0.85) cherenkovreco[i]=0.82;
    
    if(i==2)  fit->SetLineColor(kBlue);
    if(i==3)  fit->SetLineColor(kRed);
    fit->SetParameters(100,cherenkovreco[i],0.010,10);
    fit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
    fit->SetParLimits(0,0.1,1E6);
    fit->SetParLimits(1,cherenkovreco[i]-2*cut_cangle,cherenkovreco[i]+2*cut_cangle);
    fit->SetParLimits(2,0.005,0.030); // width
    hAngle[i]->Fit("fgaus","I","",cherenkovreco[i]-2*cut_cangle,cherenkovreco[i]+2*cut_cangle);
    hAngle[i]->Fit("fgaus","M","",cherenkovreco[i]-2*cut_cangle,cherenkovreco[i]+2*cut_cangle);

    cherenkovreco[i] = fit->GetParameter(1);
    spr[i] = fit->GetParameter(2);    
  }
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  TF1 *ff;
  double sep=0,esep=0, m1=0,m2=0,s1=0,s2=0; 
  if(hLnDiff[3]->GetEntries()>200){
    hLnDiff[3]->Fit("gaus","S");
    ff = hLnDiff[3]->GetFunction("gaus");
    ff->SetLineColor(1);
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
  }
  if(hLnDiff[2]->GetEntries()>200){
    hLnDiff[2]->Fit("gaus","S");
    ff = hLnDiff[2]->GetFunction("gaus");
    ff->SetLineColor(1);
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
  }
  if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));

  hAngle[2]->GetXaxis()->SetRangeUser(0.7,0.9);
  hAngle[2]->GetYaxis()->SetRangeUser(0,1.2);
  hAngle[2]->Draw();
  hAngle[3]->Draw("same");
  // fAngle[3]->Draw("same");
  // fAngle[2]->Draw("same");


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

  TLine *line3 = new TLine(0,0,0,1000);
  line3->SetLineStyle(2);
  line3->SetX1(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
  line3->SetX2(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
  line3->SetY1(0);
  line3->SetY2(1.2);
  line3->SetLineColor(1);
  line3->Draw();

  TLine *line4 = new TLine(0,0,0,1000);
  line4->SetLineStyle(2);
  line4->SetX1(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
  line4->SetX2(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
  line4->SetY1(0);
  line4->SetY2(1.2);
  line4->SetLineColor(1);
  line4->Draw();


  TLegend *leg = new TLegend(0.1,0.5,0.4,0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hAngle[2],Form("#theta_{c}^{#pi} = %2.4f rad",cherenkovreco[2]),"");
  leg->AddEntry(hAngle[3],Form("#theta_{c}^{K} = %2.4f rad",cherenkovreco[3]),"");
  leg->AddEntry(hAngle[2],Form("#sigma_{c}^{#pi} = %2.1f mrad",spr[2]*1000),"");
  leg->AddEntry(hAngle[3],Form("#sigma_{c}^{K} = %2.1f mrad",spr[3]*1000),"");
  leg->Draw();

  TLegend *lnpa = new TLegend(0.7,0.67,0.9,0.85);
  lnpa->SetFillColor(0);
  lnpa->SetFillStyle(0);
  lnpa->SetBorderSize(0);
  lnpa->SetFillStyle(0);
  lnpa->AddEntry(hAngle[2],"pions","lp");
  lnpa->AddEntry(hAngle[3],"kaons","lp");
  lnpa->Draw();
  
  // fAngle[2]->Draw("same");
  // fAngle[3]->Draw("same");
  
  glx_canvasAdd("hTime"+nid,800,400);
  
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

  glx_canvasAdd("hDiff"+nid,800,400);
  hDiff->SetLineColor(kBlack);
  hDiff->Draw();
  
  // hDiffT->SetLineColor(kRed+1);
  // hDiffT->Draw("same");
  hDiffD->SetLineColor(kGreen+2);
  hDiffD->Draw("same");
  hDiffR->SetLineColor(kBlue+1);
  hDiffR->Draw("same");

  double maxTD= hDiffD->GetXaxis()->GetBinCenter(hDiffD->GetMaximumBin());
  double maxTR= hDiffR->GetXaxis()->GetBinCenter(hDiffR->GetMaximumBin());
  double maxTT= hTime->GetXaxis()->GetBinCenter(hTime->GetMaximumBin());
  
  line = new TLine(0,0,0,1000);
  line->SetLineStyle(2);
  line->SetX1(-cut_tdiffd);
  line->SetX2(-cut_tdiffd);
  line->SetY1(0);
  line->SetY2(hDiff->GetMaximum()+0.05*hDiff->GetMaximum());
  line->SetLineColor(1);
  //line->Draw();

  line2 = new TLine(0,0,0,1000);
  line2->SetLineStyle(2);
  line2->SetX1(cut_tdiffd);
  line2->SetX2(cut_tdiffd);
  line2->SetY1(0);
  line2->SetY2(hDiff->GetMaximum()+0.05*hDiff->GetMaximum());
  line2->SetLineColor(1);
  //line2->Draw();

  TLegend *leg2 = new TLegend(0.6,0.57,0.9,0.85);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hDiff,"all","lp");
  // leg2->AddEntry(hDiffT,"MC path in EV","lp");
  // leg2->AddEntry(hDiffD,"MC path in EV for direct photons","lp");
  // leg2->AddEntry(hDiffR,"MC path in EV for reflected photons","lp");
  leg2->AddEntry(hDiffD,"direct photons","lp");
  leg2->AddEntry(hDiffR,"reflected photons","lp");
  
  leg2->Draw();
  
  glx_canvasAdd("hLnDiff"+nid,800,400);
  hLnDiff[2]->SetTitle(Form("sep = %2.2f s.d.",sep));
  hLnDiff[2]->Draw();
  hLnDiff[3]->Draw("same");

  TLegend *lnpl = new TLegend(0.7,0.67,0.9,0.85);
  lnpl->SetFillColor(0);
  lnpl->SetFillStyle(0);
  lnpl->SetBorderSize(0);
  lnpl->SetFillStyle(0);
  lnpl->AddEntry(hLnDiff[2],"pions","lp");
  lnpl->AddEntry(hLnDiff[3],"kaons","lp");
  lnpl->Draw();

  glx_canvasAdd("hNph"+nid,800,400);
    
  double nph = 0;      
  if(hNph[2]->GetEntries()>50){
    // nph = glx_fit(hNph[2],40,100,40).X();
    // auto rfit = hNph[2]->GetFunction("glx_gaust");
    // if(rfit) rfit->SetLineColor(kBlue);
    hNph[2]->SetLineColor(kBlue);
    hNph[2]->Draw();
    //glx_fit(hNph[3],40,100,40).X();
    //hNph[3]->GetFunction("glx_gaust")->SetLineColor(kRed);
    hNph[3]->SetLineColor(kRed);
    hNph[3]->Draw("same");
  }
  // hNphC->SetLineColor(kBlack);
  // hNphC->Draw("same");


  TLegend *lnph = new TLegend(0.6,0.65,0.9,0.85);
  lnph->SetFillColor(0);
  lnph->SetFillStyle(0);
  lnph->SetBorderSize(0);
  lnph->SetFillStyle(0);
  // lnph->AddEntry(hNphC,"simulated","lp");
  lnph->AddEntry(hNph[2],"pions","lp");
  lnph->AddEntry(hNph[3],"kaons","lp");
  lnph->Draw();
  
  std::cout<<"separation = "<< sep << "  nph = "<<nph <<std::endl;
  std::cout<<"maxTD "<<maxTD<<"  maxTR "<<maxTR<<std::endl;
  
  //TFile fc(infile+"_res"+nid+".root","recreate");
  // TFile fc("data/reco_lut_res/res"+nid+".root","recreate");
  // TTree *tc = new TTree("reco","reco");
  // // tc->Branch("theta",&theta,"theta/D");
  // // tc->Branch("phi",&phi,"prt_phi/D");
  // tc->Branch("sep",&sep,"sep/D");
  // tc->Branch("esep",&esep,"esep/D");
  // tc->Branch("moms",&moms,"prt_mom/D");
  // tc->Branch("xbar",&xbar,"xbar/D");
  // tc->Branch("ybar",&ybar,"ybar/D");
  // tc->Branch("nph",&nph,"nph/D");
  // tc->Branch("spr",&spr[3],"spr/D");
  // tc->Branch("maxTD",&maxTD,"maxTD/D");
  // tc->Branch("maxTR",&maxTR,"maxTR/D");
  // tc->Branch("maxTT",&maxTT,"maxTT/D");
  // tc->Fill();
  // tc->Write();
  // fc.Write();
  // fc.Close();

  glx_canvasSave(0);
  
}
