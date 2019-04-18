// gdisplay - tool to plot different quantities from the *M.root file
// original author: Roman Dzhygadlo - GSI Darmstadt 

#include "DrcHit.h"
#include "DrcEvent.h"
#define glx__sim
#include "gdisplay.h"
#include "glxtools.C"

MyMainFrame *gMain;
TGHProgressBar *pbar;
MSelector *fSelector;
TGraph *gGrDiff[glx_maxch];
TGraph *gWalk[glx_maxch];
TCanvas *cTime;

const Int_t maxMult = 30;
Int_t mult[glx_maxch]={0},gComboId(0), gTrigger(0), gMode(0), gWorkers(4), gEntries(0); //3for2015
Double_t gTimeCutMin(-10000),gTimeCutMax(10000),gTofMin(0),gTofMax(0);
Double_t gMultCutMin(0),gMultCutMax(0),gTimeCuts[glx_npmt][glx_npix][2], gTotMean[glx_npmt][glx_npix];
TString ginFile(""), gPath(""), gInfo(""),gsTimeCuts("0"), gsTotMean("0");

TH1F *hTotM[maxMult], *hLeM[maxMult];
TH1F *hTot,*hLe,*hLes,*hMult,*hCh,*hTof,*hMultEvtNum1,*hMultEvtNum2;
TH1F *hPTime[glx_npmt][glx_npix],*hPiTime[glx_npmt][glx_npix],*hSTime[glx_npmt][glx_npix],*hPTot[glx_npmt][glx_npix],*hPMult[glx_npmt][glx_npix],*hEMult[glx_npmt+1];
TH2F *hLeTot[glx_npmt][glx_npix],*hShape[glx_npmt][glx_npix];

void MSelector::Init(TTree *tree){
  fChain = tree; 
  glx_events = new TClonesArray("DrcEvent");
  fChain->SetBranchAddress("DrcEvent", &glx_events);
}

void PrintStressProgress(Long64_t total, Long64_t processed, Float_t, Long64_t){
  pbar->SetPosition(100*processed/(Float_t)total);
}

void init(){
  if(!glx_initc(ginFile,1,"")) return;
  glx_createMap();
  TString insim =  ginFile;
  insim.ReplaceAll("C.root","S.root");
  
  Long_t *id(0), *size(0), *flags(0), *modtime(0);
  if(gMode>=100 && !gSystem->GetPathInfo(insim,id,size,flags,modtime)){
    std::cout<<"Add Sim file: "<<insim <<std::endl;
    glx_ch->Add(insim);
    glx_entries = glx_ch->GetEntries();
  }
  if(gEntries>0) glx_entries=gEntries;
  
  TProof *proof;
  if(gWorkers>1){
    proof= TProof::Open(Form("workers=%d",gWorkers));
    TString dir = gSystem->pwd();
    gProof->Exec("gSystem->Load(\""+ dir+"/DrcHit_cc.so\")");
    gProof->Exec("gSystem->Load(\""+ dir+"/DrcEvent_cc.so\")");
    proof->Load("gdisplay.C+");
    glx_ch->SetProof();
    proof->SetPrintProgress(&PrintStressProgress);
  }

  fSelector = new MSelector();
  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit(1111);
}

void MSelector::SlaveBegin(TTree *){
  
  TString option = GetOption();
  std::istringstream source(option.Data());
  Int_t bins1 =100, bins2 = 100;
  Double_t min1=-5, max1=5, min2=-5, max2=5;
  source>>gMode>>gTrigger>>bins1>>min1>>max1>>bins2>>min2>>max2>>gTimeCutMin>>gTimeCutMax>>gMultCutMin>>gMultCutMax>>gsTimeCuts>>gsTotMean>>gTofMin>>gTofMax;

  TObjArray *sarr = gsTimeCuts.Tokenize(";");
  if(sarr->GetEntries()>1){
    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
	TString cut = ((TObjString *) sarr->At(2*(m*glx_npix+p)))->GetName();
	gTimeCuts[m][p][0] = cut.Atof();
	cut = ((TObjString *) sarr->At(2*(m*glx_npix+p)+1))->GetName();
	gTimeCuts[m][p][1] = cut.Atof();
      }
    }
  }
  sarr = gsTotMean.Tokenize(";");
  if(sarr->GetEntries()>1){
    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
	TString totmean = ((TObjString *) sarr->At(m*glx_npix+p))->GetName();
	gTotMean[m][p] = totmean.Atof();
	std::cout<<"gTotMean[m][p]  "<<gTotMean[m][p] <<std::endl;	
      }
    }
  }
  
  for(Int_t i=0; i<maxMult; i++){
    hTotM[i]=new TH1F(Form("hTot%d",i),Form("hTot%d",i),500,min2,max2);
    hLeM[i]=new TH1F(Form("hLe%d",i),Form("hLe%d",i),10000,-300,300);
    fOutput->Add(hTotM[i]);
    fOutput->Add(hLeM[i]);
  }

  hMultEvtNum1=new TH1F("hMultEvtNum1","",1000001,-0.5,1000000.5);
  hMultEvtNum2=new TH1F("hMultEvtNum2","",1000001,-0.5,1000000.5);
  glx_axisTime800x500(hMultEvtNum1,"event number  [#]");
  hMultEvtNum1->GetYaxis()->SetTitle("multiplicity [#]");
  glx_axisTime800x500(hMultEvtNum2,"event number  [#]");
  hMultEvtNum2->GetYaxis()->SetTitle("multiplicity [#]");
  fOutput->Add(hMultEvtNum1);
  fOutput->Add(hMultEvtNum2);

  for(Int_t m=0; m<glx_npmt; m++){
    glx_hdigi[m] = new TH2F( Form("pmt_%d", m),Form("pmt%d", m),8,0.,8.,8,0.,8.);
    glx_hdigi[m]->SetStats(0);
    glx_hdigi[m]->SetTitle(0);
    glx_hdigi[m]->GetXaxis()->SetNdivisions(10);
    glx_hdigi[m]->GetYaxis()->SetNdivisions(10);
    glx_hdigi[m]->GetXaxis()->SetLabelOffset(100);
    glx_hdigi[m]->GetYaxis()->SetLabelOffset(100);
    glx_hdigi[m]->GetXaxis()->SetTickLength(1);
    glx_hdigi[m]->GetYaxis()->SetTickLength(1);
    glx_hdigi[m]->GetXaxis()->SetAxisColor(15);
    glx_hdigi[m]->GetYaxis()->SetAxisColor(15);
    fOutput->Add(glx_hdigi[m]);

    hEMult[m]   = new TH1F(Form("emult_m%d",m),Form("pmt %d",m),  500,0,500);
    glx_axisTime800x500(hEMult[m],"multiplicity per event [#]");
    fOutput->Add(hEMult[m]);

    for(Int_t p=0; p<glx_npix; p++){
      TString name=Form("pmt %d, pixel %d, ch %d",m, p,m*64+p);
      hPTime[m][p]  = new TH1F(Form("le_pmt%dpix%d",m,p),name,bins1,min1,max1);
      hPiTime[m][p] = new TH1F(Form("lepi_pmt%dpix%d",m,p),name,bins1,min1,max1);
      hSTime[m][p]  = new TH1F(Form("les_pmt%dpix%d",m,p),name,bins1,min1,max1);
      hPTot[m][p]   = new TH1F(Form("tot_pmt%dpix%d",m,p),name,bins2,min2,max2);
      hPMult[m][p]  = new TH1F(Form("mult_pmt%dpix%d",m,p),name,50,0,50);
      
      if(gMode==1){
	hShape[m][p] = new TH2F(Form("hShape_pmt%dpix%d",m,p), Form("hShape_%d_%d;LE [ns];offset [mV]",m,p) , bins1,min1,max1,100,-5.1,15);
	hLeTot[m][p] = new TH2F(Form("hLeTot_pmt%dpix%d" ,m,p), Form("pmt %d, pixel %d;LE [ns];TOT [ns]",m, p), 200,min1,max1, bins2,min2,max2);
	glx_axisTime800x500(hShape[m][p],"time [ns]");
	hShape[m][p]->GetYaxis()->SetRangeUser(-4,10);
	hShape[m][p]->GetYaxis()->SetTitle("offset to the threshold [mV]");
      
	fOutput->Add(hLeTot[m][p]);
	fOutput->Add(hShape[m][p]);
      }

      glx_axisTime800x500(hPTime[m][p]);
      glx_axisTime800x500(hPiTime[m][p]);
      glx_axisTime800x500(hPTot[m][p],"TOT time [ns]");
      glx_axisTime800x500(hPMult[m][p],"multiplicity [#]");

      hSTime[m][p]->SetLineColor(2);
      hPTime[m][p]->SetLineColor((gMode==100)? 1: kRed+1);
      hPiTime[m][p]->SetLineColor(kBlue+1);
      
      
      fOutput->Add(hPTime[m][p]);
      fOutput->Add(hPiTime[m][p]);
      fOutput->Add(hSTime[m][p]);
      fOutput->Add(hPTot[m][p]);
      fOutput->Add(hPMult[m][p]);
    }
  }

  hTot=new TH1F("hTotA","",100,min2,max2);
  hLe=new TH1F("hLeA","",600,-300,300);
  hLes=new TH1F("hLeAs","",100,0,100);
  hMult=new TH1F("hMultA","",50,0,50);
  hCh=new TH1F("hChA","",glx_maxch,0,glx_maxch);
  //hTof=new TH1F("hTof","",2000,25,40);
  hTof=new TH1F("hTof","",2000,25,70);
 
  glx_axisTime800x500(hTot,"TOT time [ns]");
  glx_axisTime800x500(hLe,"LE time [ns]");
  glx_axisTime800x500(hMult,"multiplicity [#]");
  hCh->SetLineColor(1);
  
  gStyle->SetOptStat(1001111);
  
  fOutput->Add(hLe);
  fOutput->Add(hLes);
  fOutput->Add(hTot);
  fOutput->Add(hMult);
  fOutput->Add(hEMult[glx_npmt]);
  fOutput->Add(hCh);
  fOutput->Add(hTof);
}

Bool_t MSelector::Process(Long64_t entry){
  GetEntry(entry);
  
  for (Int_t t=0; t<glx_events->GetEntriesFast(); t++){
    glx_nextEventc(entry,t,100);
    Int_t particleId=glx_event->GetPdg();
    if(glx_event->GetType()!=2) continue; //beam
    //if(glx_event->GetType()!=1) continue; //LED

    Double_t timeres(0);
    TString fileid("");
    bool bsim(false);
    Double_t offset=0;
    Double_t evtTime=glx_event->GetTime();
    
    Double_t le,tot, triggerLe(0),toftime(0), tof(0), tof1(0),tof2(0);
    DrcHit hit;
    Int_t pmt,pix,col,row,ch,chMultiplicity(0);
    Int_t thitCount1(0), thitCount2(0), hitCount1(0), hitCount2(0);
    memset(mult, 0, sizeof(mult));

    Int_t nhits = glx_event->GetHitSize();
    
    if(gTrigger>-1){
      for(Int_t h=0; h<nhits; h++){
	hit = glx_event->GetHit(h);
	ch  = hit.GetChannel();
	mult[ch]++;
	if(ch == gTrigger && gTrigger>0) triggerLe = hit.GetLeadTime();
      }
    }
    
    for(Int_t h=0; h<nhits; h++){
      hit = glx_event->GetHit(h);
      ch  = hit.GetChannel();
      pmt = ch/64;
      pix = ch%64;
      row = pix/8;
      col = pix%8;
      
      if(!bsim){
	hCh->Fill(ch);
	chMultiplicity = mult[ch];
	hMult->Fill(chMultiplicity);
      }
    
      if(gMultCutMin!=gMultCutMax && (thitCount2<gMultCutMin || thitCount2>gMultCutMax)) continue;
      if(gTofMin != gTofMax && (tof<gTofMin || tof>gTofMax)) continue;
      
      le = hit.GetLeadTime()-evtTime;      
      // if((ch/64)%18>6) le -=10;
      // if((ch/64)%18>12) le -=10; 

      tot = hit.GetTotTime();   
      //if(tot<40 || tot>42) continue;
      
      Double_t timeDiff = le-triggerLe;
    
      if(gsTimeCuts!="0" && (timeDiff<gTimeCuts[pmt][pix][0] || timeDiff>gTimeCuts[pmt][pix][1])) continue;
      else if(gTimeCutMin!=gTimeCutMax &&  (timeDiff<gTimeCutMin || timeDiff>gTimeCutMax)) continue;
      hitCount2++;
    
      if(gsTotMean!="0"){
	timeDiff += 0.3*(tot - gTotMean[pmt][pix]);
      }
      // if(timeDiff>30 || timeDiff<10 ) continue;
      if(triggerLe!=-1 || gTrigger==0) {
	//Double_t offset = 284.89;
	if(bsim){
	  if(particleId==321) hSTime[pmt][pix]->Fill(timeDiff);
	  hLes->Fill(le);
	  continue;
	}
	//timeDiff-=offset;
	if(gMode==1){
	  hLeTot[pmt][pix]->SetTitle(Form("ch %d",ch));
	  hLeTot[pmt][pix]->Fill(timeDiff,tot);
	  hShape[pmt][pix]->Fill(timeDiff,offset);
	  hShape[pmt][pix]->Fill(timeDiff + tot,offset);
	}	
	glx_hdigi[pmt]->Fill(col, row);
	if(particleId==321 || particleId==0) hPTime[pmt][pix]->Fill(timeDiff);
	else hPiTime[pmt][pix]->Fill(timeDiff);
	
      }
      hPTot[pmt][pix]->Fill(tot);
      hPTot[pmt][pix]->SetTitle(Form("pmt %d, pixel %d, channel %d",pmt, pix, ch));
      
      hPMult[pmt][pix]->Fill(chMultiplicity);
      hPMult[pmt][pix]->SetTitle(Form("pmt %d, pixel %d, channel %d",pmt, pix, ch));
      
      if(chMultiplicity<maxMult){
	hTotM[chMultiplicity]->Fill(tot);
	hLeM[chMultiplicity]->Fill(le);
      }
      hTot->Fill(tot);
      hLe->Fill(le);
    }

    if(!bsim){
      if(hitCount1+hitCount2 !=0) hEMult[0]->Fill(hitCount1+hitCount2);
      if(hitCount1!=0) hEMult[1]->Fill(hitCount1);
      if(hitCount2!=0) hEMult[2]->Fill(hitCount2);
      hMultEvtNum1->Fill(entry,hitCount1);
      hMultEvtNum2->Fill(entry,hitCount2);
    }

  }
  return kTRUE;
}

void calculateTimeCut(){
  TVector3 res;
  gsTimeCuts="";
  for (Int_t m=0; m <glx_npmt; m++) {
    for(Int_t p=0; p<glx_npix; p++){
      res = glx_fit(hPTime[m][p]);
      gTimeCuts[m][p][0]=res.X() - 3*res.Y();
      gTimeCuts[m][p][1]=res.X() + 3*res.Y();
      gsTimeCuts+=Form("%f;%f;",gTimeCuts[m][p][0],gTimeCuts[m][p][1]);
    }
  }
}

void getTimeOffset(){
  std::cout<<"Creating calibration"<<std::endl;
  TH1D* h;
  TH2F* hh;
  for (Int_t m=0; m <glx_npmt; m++) {
    for(Int_t p=0; p<glx_npix; p++){
      hh =(TH2F*) hLeTot[m][p]->Clone("hh");

      Double_t mean = hPTime[m][p]->GetXaxis()->GetBinCenter(hPTime[m][p]->GetMaximumBin());
      
      Int_t threshold =  hPTime[m][p]->GetMaximum()*0.3;      
      Int_t firstbin = hPTime[m][p]->FindFirstBinAbove(threshold);
      Double_t front = hPTime[m][p]->GetXaxis()->GetBinCenter(firstbin);
      
      hh->RebinY(2);

      TCutG *cutg = new TCutG("onepeakcut",5);
      cutg->SetVarX("y");
      cutg->SetVarY("x");
      Double_t range=1.85;
      cutg->SetPoint(0,front-range,0.5);
      cutg->SetPoint(1,front-range,8.5);
      cutg->SetPoint(2,front+range,8.5);
      cutg->SetPoint(3,front+range,0.5);
      cutg->SetPoint(4,front-range,0.5);
  
      // TF1 *gaust = new TF1("gaust","gaus(0)",85,105);
      // gaust->SetParameter(1,90);
      // gaust->SetParameter(2,0.3);
      // hh->FitSlicesX(gaust,0,-1,10,"");
      // TGraph * gg = new TGraph((TH1D*)gDirectory->Get("hh_1")); 
      Int_t ch = map_mpc[m][p];

      gGrDiff[ch] = new TGraph();      
      gWalk[ch] = new TGraph();
      Double_t pvx(mean);
      TGraph* gsmooth = new TGraph();
      TGraph* gsmoothw = new TGraph();
      for (int i=0;i<hh->GetNbinsY();i++){
	Double_t x = hh->GetYaxis()->GetBinCenter(i);
	Double_t vx(0);
 
	h = hh->ProjectionX(Form("bin%d",i),i,i,"[onepeakcut]");
	Int_t fb = h->FindFirstBinAbove(h->GetMaximum()*0.6);
	vx = h->GetXaxis()->GetBinCenter(fb);

        if(vx==0 || fabs(vx-front)>1.5) vx = front;
	
	if(fabs(pvx-vx)>0.05) vx += 0.5*(pvx-vx);
	
	gsmooth->SetPoint(i,x,vx);
	gsmoothw->SetPoint(i,x,vx-front);
	
	pvx=vx;

	// cTime->cd();
	// h->Draw();	
	// cTime->Update();
	// cTime->WaitPrimitive();
      }
      
      // gWalk[ch] = glx_smooth(gsmoothw,10);
      // gGrDiff[ch] = glx_smooth(gsmooth,10);
      
      gGrDiff[ch]->SetName(Form("gCalib_ch%d",ch));
      gGrDiff[ch]->GetXaxis()->SetTitle("fine bin [#]");
      gGrDiff[ch]->GetYaxis()->SetTitle("fine time [ns]");
    
    }
  }
}

void MyMainFrame::DoExportOffsets(){
  if(gMode==1) {
    getTimeOffset();

    TString filedir=ginFile;
    filedir.Remove(filedir.Last('/'));
    TFile efile(filedir+"/calib_walk.root","RECREATE");
    Int_t c;
    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
	c = map_mpc[m][p];
	gWalk[c]->SetName(Form("walk_%d",c));
	gWalk[c]->Write();
      }
    }
    efile.Write();
    efile.Close();
    std::cout<<"Exporting .. Done"<<std::endl;
  }else{
    std::cout<<"For exporting use -a1 flag"<<std::endl;
  }
}

TLine *gLine1 = new TLine(0,0,0,1000);
Bool_t glock= false;
void exec3event(Int_t event, Int_t gx, Int_t gy, TObject *selected){
  if(gComboId==0 || gComboId==2 || gComboId==5 || gComboId==4 || gComboId==10 || gComboId==11 || gComboId==7){
    TCanvas *c = (TCanvas *) gTQSender;
    TPad *pad = (TPad *) c->GetSelectedPad();
    if (!pad) return;
    Float_t x = pad->AbsPixeltoX(gx);
    Float_t y = pad->AbsPixeltoY(gy);
    x = pad->PadtoX(x);
    y = pad->PadtoY(y);
    if(event ==1 && glock) glock = false;
    else if(event ==1) glock = true;
    if(glock) return;
    if (selected->InheritsFrom(TH2::Class())){
      TH2F *hDigi = (TH2F *) selected;
      Int_t binx = hDigi->GetXaxis()->FindBin(x);
      Int_t biny = hDigi->GetYaxis()->FindBin(y);
      TString spmt = selected->GetName();
      spmt = spmt(4,spmt.Sizeof());
      Int_t m = spmt.Atoi();      
      Int_t p = 8*(biny-1)+binx-1;
      Int_t ch = m*64+p; //map_mpc[m][p];
      //printf("Canvas %s: event=%d, x=%d, y=%d, p=%d, selected=%d\n", spmt.Data(), event, binx, biny, p,spmt.Atoi());
      cTime->cd();
      if(gComboId==0) {
	TH1F * hh[] = {hPTime[m][p],(gMode==100)? hSTime[m][p]: hPiTime[m][p]};
	glx_normalize(hh,2);
	glx_fit(hh[0],1,20,10);	
	hh[0]->Draw();
	if(hh[0]->GetEntries()>10) hh[1]->Draw("same");
      }
      if(gComboId==2) hPTot[m][p]->Draw();   
      if(gComboId==5) hPMult[m][p]->Draw();      
      if(gComboId==4) hLeTot[m][p]->Draw("colz");
      if(gComboId==10) hShape[m][p]->Draw("colz");
      if(gComboId==11){
	hLeTot[m][p]->Draw("colz");
	if(gGrDiff[ch]){
	  Double_t* xx = gGrDiff[ch]->GetX();
	  Double_t* yy = gGrDiff[ch]->GetY();

	  TGraph* gr = new TGraph(gGrDiff[ch]->GetN(),yy,xx);
	  gr->SetMarkerStyle(7);
	  gr->SetMarkerColor(2); 
	
	  gr->Draw("PL same");
	  // glx_smooth(gr,10)->Draw("PL same");
	}else{
	  std::cout<<"Press the btn \"Export offsets\" "<<std::endl;
	  
	}
      }
      if(gMain->fCheckBtn2->GetState() == kButtonDown){
	gMain->fEdit3->SetText(Form("%2.2f %2.2f", gTimeCuts[m][p][0], gTimeCuts[m][p][1]));
      }
      if(gComboId==7){
	gLine1->SetX1(ch+0.5);
	gLine1->SetX2(ch+0.5);
	gLine1->SetY1(cTime->GetUymin());
	gLine1->SetY2(cTime->GetUymax());
	gLine1->SetLineColor(kRed);
	gLine1->Draw();

	hCh->SetBit(TH1::kNoStats);
	TPaveStats *ps = (TPaveStats*)hCh->GetListOfFunctions()->FindObject("stats");
	TList *list = ps->GetListOfLines();
	list->RemoveAll();

        ps->AddText(Form("entries %d",(Int_t)hCh->GetEntries()));
	ps->AddText(Form("global ch = %d",ch));
	ps->AddText(Form("pmt = %d",ch/64));
	ps->AddText(Form("pixel = %d",ch%64));
	ps->AddText(Form("ssp slot = %d",map_ssp_slot[ch]));
	ps->AddText(Form("ssp fiber = %d",map_ssp_fiber[ch]));
      }
      cTime->Update();
      gMain->fNumber2->SetIntNumber(ch);
    }
  }
}

void exec4event(Int_t event, Int_t gx, Int_t gy, TObject *selected){
  if(gComboId!=7) return;
  TCanvas *c = (TCanvas *) gTQSender;
  TPad *pad = (TPad *) c->GetSelectedPad();
  if (!pad) return;
  Float_t x = pad->AbsPixeltoX(gx);
  Float_t y = pad->AbsPixeltoY(gy);
  x = pad->PadtoX(x);
  y = pad->PadtoY(y);
  if (selected->InheritsFrom(TAxis::Class())) return;
  
  TH1F *hChannel = (TH1F *) pad->FindObject("hChA");
    
  Int_t ch = hChannel->GetXaxis()->FindBin(x);
  cTime->cd();
  gLine1->SetX1(ch+0.5);
  gLine1->SetX2(ch+0.5);
  gLine1->SetY1(cTime->GetUymin());
  gLine1->SetY2(cTime->GetUymax());
  gLine1->SetLineColor(kRed);
  gLine1->Draw();
  hCh->SetBit(TH1::kNoStats);
  TPaveStats *ps = (TPaveStats*)hCh->GetListOfFunctions()->FindObject("stats");
  TList *list = ps->GetListOfLines();
  list->RemoveAll();
  
  ps->AddText(Form("entries %d",(Int_t)hCh->GetEntries()));
  ps->AddText(Form("global ch = %d",ch));
  ps->AddText(Form("pmt = %d",ch/64));
  ps->AddText(Form("pixel = %d",ch%64));
  ps->AddText(Form("ssp slot = %d",map_ssp_slot[ch]));
  ps->AddText(Form("ssp fiber = %d",map_ssp_fiber[ch]));

  cTime->Update();
}

void MyMainFrame::InterestChanged(){
  Int_t ch = fNumber2->GetIntNumber();
  std::cout<<"ch of interest "<<ch <<std::endl;
  Int_t pmt = map_pmt[ch];
  Int_t pix = map_pix[ch];

  if(gComboId==0) {
    TH1F * hh[] = {hPiTime[pmt][pix],hSTime[pmt][pix]}; 
    glx_normalize(hh,2);
    hh[0]->Draw();
    if(gMode==100) hSTime[pmt][pix]->Draw("same");
    else hPiTime[pmt][pix]->Draw("same");
    glx_fit(hh[0],10,100);
    if(hh[0]->GetEntries()>10) hh[1]->Draw("same");
  }
  if(gComboId==2) hPTot[pmt][pix]->Draw();   
  if(gComboId==5) hPMult[pmt][pix]->Draw();      
  if(gComboId==4) hLeTot[pmt][pix]->Draw("colz");
  if(gComboId==10) hShape[pmt][pix]->Draw("colz");
  if(gComboId==11){
    Int_t ch = map_mpc[pmt][pix];
    hLeTot[pmt][pix]->Draw("colz");
    Double_t* xx = gGrDiff[ch]->GetX();
    Double_t* yy = gGrDiff[ch]->GetY();

    TGraph* gr = new TGraph(gGrDiff[ch]->GetN(),yy,xx);
    gr->SetMarkerStyle(7);
    gr->SetMarkerColor(2);
    gr->Draw("P same");
  }
  cTime->Update();
}

void MSelector::Terminate(){


  for (Int_t m=0; m <glx_npmt; m++) {
    glx_hdigi[m] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("pmt_%d",m), fOutput));
    hEMult[m] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("emult_m%d",m), fOutput)); 
    for(Int_t p=0; p<glx_npix; p++){
      hPTime[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("le_pmt%dpix%d",m,p), fOutput));
      hPiTime[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("lepi_pmt%dpix%d",m,p), fOutput)); 
      hSTime[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("les_pmt%dpix%d",m,p), fOutput)); 
      hPTot[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("tot_pmt%dpix%d",m,p), fOutput)); 
      hPMult[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("mult_pmt%dpix%d",m,p), fOutput)); 
      if(gMode==1){
	hLeTot[m][p] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hLeTot_pmt%dpix%d",m,p), fOutput)); 
	hShape[m][p] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hShape_pmt%dpix%d",m,p), fOutput));
      }
    }
  }
  for(Int_t i=0; i<maxMult; i++){
    hLeM[i] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hLe%d",i), fOutput)); 
    hTotM[i] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTot%d",i), fOutput)); 
  }

  hLe = dynamic_cast<TH1F *>(TProof::GetOutput("hLeA", fOutput));   
  hLes = dynamic_cast<TH1F *>(TProof::GetOutput("hLeAs", fOutput));   
  hTot = dynamic_cast<TH1F *>(TProof::GetOutput("hTotA", fOutput)); 
  hMult = dynamic_cast<TH1F *>(TProof::GetOutput("hMultA", fOutput)); 
  hEMult[glx_npmt] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("emult_m%d",glx_npmt), fOutput)); 
  hMultEvtNum1 = dynamic_cast<TH1F *>(TProof::GetOutput("hMultEvtNum1", fOutput));
  hMultEvtNum2 = dynamic_cast<TH1F *>(TProof::GetOutput("hMultEvtNum2", fOutput));

  hCh = dynamic_cast<TH1F *>(TProof::GetOutput("hChA", fOutput));
  hTof = dynamic_cast<TH1F *>(TProof::GetOutput("hTof", fOutput)); 
}

void MyMainFrame::DoDraw(){
  fCheckBtn2->SetText("3 sigma time cut                         ");
  fCheckBtn2->SetTextColor(Pixel_t(0x000000),kFALSE);   
  fCheckBtn3->SetText("Walk correction                         ");
  fCheckBtn3->SetTextColor(Pixel_t(0x000000),kFALSE);  
  if(gTrigger>-1) gTrigger = fNumber->GetIntNumber();
  

  for(Int_t m=0; m<glx_npmt; m++){
    if(glx_hdigi[m]) glx_hdigi[m]->Reset();
    if(hEMult[m]) hEMult[m]->Reset();
    for(Int_t p=0; p<glx_npix; p++){
      if(hPTime[m][p]) hPTime[m][p]->Reset();
      if(hPiTime[m][p]) hPiTime[m][p]->Reset();
      if(hSTime[m][p]) hSTime[m][p]->Reset();
      if(hPTot[m][p]) hPTot[m][p]->Reset();
      if(hPMult[m][p]) hPMult[m][p]->Reset();
    }
  }
  if(hLe) hLe->Reset();
  if(hLes) hLes->Reset();
  if(hTot) hTot->Reset();
  if(hMult) hMult->Reset();
  if(hCh) hCh->Reset();
  if(hTof) hTof->Reset(); 

  fHProg3->Reset();
  //glx_entries = 10000;

  TString option = Form("%d %d %s %s %s %s %s %s %s",gMode,gTrigger,fEdit1->GetText(),fEdit2->GetText(),fEdit3->GetText(),fEdit4->GetText(),gsTimeCuts.Data(), gsTotMean.Data(),fEdit5->GetText());

  gROOT->SetBatch(1);
  glx_ch->Process(fSelector,option,glx_entries);
  gROOT->SetBatch(0);
 
  if(gTrigger<0) gTrigger=0;

  Int_t tmax, max=0;
  for(Int_t p=0; p<glx_npmt;p++){
    tmax = glx_hdigi[p]->GetMaximum();
    if(max<tmax) max = tmax;
  }
  max = (max>0)? max : 1;
  fHslider1->SetRange(0,max);
  fHslider1->SetPosition(max);

  fCheckBtn1->SetState(kButtonUp);

  glx_drawDigi("m,p,v\n",glx_geometry);

  updatePlot(gComboId);

  if(gMode==10) {
    gComboId=7;
    DoExport();
  }

  if(gROOT->IsBatch() && gMode==100) {
    DoExport();
  }

  if(gMode==101){
    double xmax1 = hLe->GetXaxis()->GetBinCenter(hLe->GetMaximumBin());
    double xmax2 = hLes->GetXaxis()->GetBinCenter(hLes->GetMaximumBin());
    
    TString fileid = ginFile;
    TString dirid = ginFile;
    fileid.Remove(0,fileid.Last('_')+1);
    fileid.Remove(fileid.Last('C'));
    dirid.Remove(dirid.Last('/'));
    TString line = Form("if(\"%s\" == fileid) offset = %f;\n", fileid.Data(), xmax1-xmax2);
    std::cout<<"line  "<<line <<std::endl;

    std::ofstream out;
    out.open(dirid+"/offsets.txt", std::ios::app);
    out << line;
    out.close();
    
    TFile efile(dirid+"/off_"+fileid+ ".root","RECREATE");
    TGraph *gr = new TGraph();
    gr->SetPoint(0,xmax1-xmax2,  xmax1-xmax2);
    gr->SetName(fileid);
    gr->Write();
    efile.Write();
    efile.Close();

    gApplication->Terminate(0);
  }

}

TH2F* glx_hdigi_temp_updateplot[glx_npmt];
TString MyMainFrame::updatePlot(Int_t id, TCanvas *cT){
  if(!cT) cT = cTime;
  TString histname="";
  gComboId = id;
  cT->cd();
  cT->SetLogy(0);
  Int_t max = 0;
  TLegend *leg;

  if(fBackToHp){
    for(Int_t m=0; m<glx_npmt; m++){
      glx_hdigi[m] = glx_hdigi_temp_updateplot[m];
    }
    glx_drawDigi("m,p,v\n",glx_geometry);
    fBackToHp=false;
  }


  switch(id) {
  case 0: // LE 
    break;
  case 1: // LE All
    cT->SetLogy(1);
    hLe->Draw();
    hLes->Draw("same");
    if(gMode!=100){
      for(Int_t i=0; i<maxMult; i++){
	Int_t colorid = i+1;
	if(colorid ==5) colorid=9;
	if(colorid ==3) colorid=8;
	if(hTotM[i]->GetEntries()>0){
	  hLeM[i]->SetLineColor(colorid);
	  if(i==0) hLeM[i]->SetLineWidth(2);
	  hLeM[i]->Draw("same");
	}
      }
    }
    histname=hLe->GetName();
    break;
  case 2: // ToT
    break;
  case 3: // ToT All
    hTot->Draw();
    for(Int_t i=0; i<maxMult; i++){
      Int_t colorid = i+1;
      if(colorid ==5) colorid=9;
      if(colorid ==3) colorid=8;
      if(hTotM[i]->GetEntries()>0){
	hTotM[i]->SetLineColor(colorid);
	if(i==0) hTotM[i]->SetLineWidth(2);
	hTotM[i]->Draw("same");
      }
    } 
    histname=hTot->GetName();
    break;
  case 4: // Le vs. ToT
    break;
  case 5: // Mult
    break;
  case 6: // Mult All
    cT->SetLogy(1);
    hMult->Draw();
    histname=hMult->GetName();
    break;
  case 7: // Channels
    hCh->Draw();
    histname=hCh->GetName();
    break;
  case 12: // TOF
    hTof->Draw();
    histname=hTof->GetName();
    break;
  case 8: // Event Mult
    if(hEMult[0]->GetMaximum()>max) max = hEMult[0]->GetMaximum();
    if(hEMult[1]->GetMaximum()>max) max = hEMult[1]->GetMaximum();
    if(hEMult[2]->GetMaximum()>max) max = hEMult[2]->GetMaximum();
    
    hEMult[0]->SetMaximum(max+0.05*max);
    hEMult[0]->Draw();
    hEMult[1]->SetLineColor(2);
    hEMult[1]->Draw("same");
    hEMult[2]->SetLineColor(4);
    hEMult[2]->Draw("same");
    leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hEMult[0],"Total","l");
    leg->AddEntry(hEMult[2],"Dirc PMTs","l");
    leg->AddEntry(hEMult[1],"Rest","l");
    leg->Draw();
    histname=hEMult[0]->GetName();
    break;
  case 9: // Mult vs. Event number
    if(hMultEvtNum1->GetMaximum()>max) max = hMultEvtNum1->GetMaximum();
    if(hMultEvtNum2->GetMaximum()>max) max = hMultEvtNum2->GetMaximum();
    
    hMultEvtNum1->SetMaximum(max+0.05*max);
    hMultEvtNum1->SetMinimum(0);
    hMultEvtNum1->SetLineColor(2);
    hMultEvtNum1->Draw("hist");
    hMultEvtNum2->SetLineColor(4);
    hMultEvtNum2->Draw("same hist");
    leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hMultEvtNum2,"Dirc PMTs","l");
    leg->AddEntry(hMultEvtNum1,"Rest","l");
    leg->Draw();
    histname=hMultEvtNum1->GetName();
    break;
  case 13:
    TVector3 res;
    for(Int_t m=0; m<glx_npmt; m++){
      glx_hdigi_temp_updateplot[m] = (TH2F*)glx_hdigi[m]->Clone();
      if(glx_hdigi[m]) glx_hdigi[m]->Reset("M");
    }
    TH1F *hSigma = new TH1F("hSigma",";#sigma [ns];entries [#]",1000,0,10);

    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
	Int_t col = p/8;
	Int_t row = p%8;
	Double_t sigma = glx_fit(hPTime[m][p],10,20,10,1,"0").Y();
	hSigma->Fill(sigma);
	glx_hdigi[m]->Fill(row,col,sigma);
      }
    }
    glx_drawDigi("m,p,v\n",glx_geometry,2,0);
    
    fBackToHp=true;

    cT->cd();
    //hSigma->Fit("gaus");
    hSigma->Draw();
    
    break;    
  }
  cT->Modified();
  cT->Update();
  return histname;
}

void MyMainFrame::DoMore(){
  if(GetState(fHm)) {
    HideFrame(fHm);
    fBtnMore->SetText("&More");
  }
  else{
    ShowFrame(fHm);
    fBtnMore->SetText("&Less");
  }
}
TCanvas *cExport;
void MyMainFrame::DoSavePng(){
  gROOT->SetBatch(1);
  cExport = new TCanvas("cExport","cExport",0,0,800,400);
  cExport->SetName("current");
  cExport->SetCanvasSize(800,400);
  
  Int_t saveFlag = 1;
  TString histname="", filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  if(glx_savepath == "") glx_createSubDir(filedir+"/auto");
  
  TObject *obj; 
  TIter next(cTime->GetListOfPrimitives());
  obj=next(); obj=next(); obj->Draw();
  while ((obj=next())) obj->Draw("same");
  
  glx_canvasAdd(cExport);
  glx_canvasSave(1,0);
  gROOT->SetBatch(0);
  std::cout<<"Save current png .. done"<<std::endl;
  
}

void MyMainFrame::DoExport(){

  gROOT->SetBatch(1);
  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  cExport->SetCanvasSize(800,400);
  Int_t saveFlag = 1;
  TString histname="", filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  if(glx_savepath == "") glx_createSubDir(filedir+"/auto");

  //glx_canvasAdd("hp",800,400);
  glx_canvasAdd(glx_cdigi); 
  //glx_cdigi->DrawClonePad();
  glx_canvasSave(1,0);
  
  std::cout<<"Exporting into  "<<glx_savepath <<std::endl;
  glx_writeString(glx_savepath+"/digi.csv", glx_drawDigi("m,p,v\n",glx_geometry));
  
  pbar->Reset();
  Float_t total = (glx_npmt-1)*(glx_npix-1);
  if(gComboId==0 || gComboId==2 || gComboId==5 || gComboId==4 || gComboId==10 || gComboId==11){
    for(Int_t m=0; m<glx_npmt; m++){
      for(Int_t p=0; p<glx_npix; p++){
	cExport->cd();
	if(gComboId==0) {
	  TH1F * hh[] = {hPTime[m][p],(gMode==100)? hSTime[m][p]: hPiTime[m][p]};
	  glx_normalize(hh,2);
	  glx_fit(hh[0],10,1);
	  hh[0]->Draw();
	  if(hh[1]->GetEntries()>10) hh[1]->Draw("same");
	  histname=hh[0]->GetName();
	}
	if(gComboId==2){
	  hPTot[m][p]->Draw(); 
	  histname=hPTot[m][p]->GetName();
	}  
	if(gComboId==5){
	  hPMult[m][p]->Draw();  
	  histname=hPMult[m][p]->GetName();    
	}
	if(gComboId==4){
	  hLeTot[m][p]->Draw("colz");
	  histname=hLeTot[m][p]->GetName();
	}
	if(gComboId==10){
	  hShape[m][p]->Draw("colz");
	  histname=hShape[m][p]->GetName();
	}
	
	cExport->SetName(histname);
	glx_canvasAdd(cExport);
	glx_canvasSave(1,0);
	
	pbar->SetPosition(100*(m*p)/total);
	gSystem->ProcessEvents();
      }
    }
  }else{
    histname = updatePlot(gComboId,cExport);
    cExport->SetName(histname);
    glx_canvasAdd(cExport);
    glx_canvasSave(1,0);
  }

  // if( gMode==100){
  //   for(Int_t m=0; m<glx_npmt; m++){
  //     for(Int_t p=0; p<glx_npix; p++){
  // 	cExport->cd();
  // 	TH1F * hh[] = {hPTime[m][p],hPiTime[m][p]}; 
  // 	glx_normalize(hh,2);
  // 	hh[0]->Draw();
  // 	glx_fit(hh[0],1,1);
  // 	hh[0]->Draw("same");
  // 	if(hh[1]->GetEntries()>10) hh[1]->Draw("same");

  // 	cExport->SetName(hh[0]->GetName());
  // 	canvasAdd(cExport);
  // 	canvasSave(1,0);
  //     }
  //   }
  // }

  gROOT->SetBatch(0);
  std::cout<<"Exporting .. Done"<<std::endl;
  
}

TLine *gLine = new TLine(0,0,3000,0);
void MyMainFrame::DoSlider(Int_t pos){
  if(fCheckBtn1->GetState() != kButtonDown)
    glx_drawDigi("m,p,v\n",glx_geometry,pos);

  if(gComboId==7){
    cTime->cd();  
    gLine->SetY1(pos);
    gLine->SetY2(pos);
    gLine->SetLineColor(kRed);
    gLine->Draw();
    cTime->Update();
  }
  
}

void MyMainFrame::DoCheckBtnClecked1(){
  Int_t state = -1;
  if(fCheckBtn1->GetState() != kButtonDown){
    state = fHslider1->GetPosition();
    gMain->fHslider1->SetEnabled(kTRUE);
  }else{
    gMain->fHslider1->SetEnabled(kFALSE);
  }
  glx_drawDigi("m,p,v\n",glx_geometry,state);
}

void MyMainFrame::DoCheckBtnClecked2(){
  Int_t state = -1;
  fCheckBtn2->SetTextColor(Pixel_t(0xff0000),kFALSE);      
  fCheckBtn2->SetText("3 sigma time cut (Press Evaluate)");
  if(fCheckBtn2->GetState() == kButtonDown){
    calculateTimeCut();    
    gMain->fEdit3->SetEnabled(kFALSE);
  }
  if(fCheckBtn2->GetState() == kButtonUp){  
    gsTimeCuts="0";
    gMain->fEdit3->SetText("0 0");
    gMain->fEdit3->SetEnabled(kTRUE);
    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
	gTimeCuts[m][p][0]=0;
	gTimeCuts[m][p][1]=0;
      }
    }
  }
}

void MyMainFrame::DoCheckBtnClecked3(){
  Int_t state = -1;
  fCheckBtn3->SetTextColor(Pixel_t(0xff0000),kFALSE);       
  fCheckBtn3->SetText("Walk correction (Press Evaluate)");
  if(fCheckBtn3->GetState() == kButtonDown){
    // get mean TOT for each pixel
    gsTotMean="";
    TVector3 res;
    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
	res = glx_fit(hPTot[m][p]);
	gsTotMean+=Form("%f;",res.X());
      }
    }  
  }
  if(fCheckBtn3->GetState() == kButtonUp){      
    gsTotMean="0";
  }
}

TH2F* glx_hdigi_temp[glx_npmt];
void MyMainFrame::DoCheckBtnClecked4(){
  if(fCheckBtn4->GetState() == kButtonDown){
    for(Int_t m=0; m<glx_npmt; m++){
      glx_hdigi_temp[m] = (TH2F*)glx_hdigi[m]->Clone();
      if(glx_hdigi[m]) glx_hdigi[m]->Reset();
    }
    for (Int_t m=0; m <glx_npmt; m++) {
      for(Int_t p=0; p<glx_npix; p++){
    	Int_t col = p/8;
    	Int_t row = p%8;
    	Double_t mean = glx_fit(hPiTime[m][p],10,20,10,1,"0").X();
    	//if(mean>90) mean = 90; 
    	glx_hdigi[m]->Fill(row,col,mean);
      }
    }
    glx_drawDigi("m,p,v\n",glx_geometry,-2,-2);
  }
  if(fCheckBtn4->GetState() == kButtonUp){
    for(Int_t m=0; m<glx_npmt; m++){
      glx_hdigi[m] = glx_hdigi_temp[m];
    }
    glx_drawDigi("m,p,v\n",glx_geometry);
  }
}

TH2F* glx_hdigi_history[glx_npmt];
void MyMainFrame::DoHistory(){
    for(Int_t m=0; m<glx_npmt; m++){
      glx_hdigi_temp[m] = (TH2F*)glx_hdigi[m]->Clone();
      if(glx_hdigi_history[m]){
	glx_hdigi[m] = (TH2F*)glx_hdigi_history[m]->Clone();   
      }
      glx_hdigi_history[m] = (TH2F*)glx_hdigi_temp[m]->Clone();
    }
    glx_drawDigi("m,p,v\n",glx_geometry);
}

void MyMainFrame::DoExit(){
  gApplication->Terminate(0);
}

void MyMainFrame::SetStatusText(const char *txt, Int_t pi){
  fStatusBar->SetText(txt,pi);
}

void MyMainFrame::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected){
  const char *text0, *text1, *text3;
  char text2[50];
  text0 = selected->GetTitle();
  SetStatusText(text0,0);
  text1 = selected->GetName();
  SetStatusText(text1,1);
  if (event == kKeyPress)
    sprintf(text2, "%c", (char) px);
  else
    sprintf(text2, "%d,%d", px, py);
  SetStatusText(text2,2);
  text3 = selected->GetObjectInfo(px,py);
  SetStatusText(text3,3);
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h) : TGMainFrame(p, w, h){

  // Create the embedded canvas
  fEcan = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid0 = fEcan->GetCanvasWindowId();
  glx_cdigi = new TCanvas("glx_cdigi",10,10,wid0);
  glx_cdigi->SetMargin(0,0,0,0);
  fEcan->AdoptCanvas(glx_cdigi);

  fTime = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid1 = fTime->GetCanvasWindowId();
  cTime = new TCanvas("cTime",10,10,wid1);
  fTime->AdoptCanvas(cTime);

  AddFrame(fEcan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  AddFrame(fTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  // =============== More frame ===================
  fHm = new TGVerticalFrame(this, 200, 120, kFixedSize);

  TGHorizontalFrame *fHm0 = new TGHorizontalFrame(fHm, 400, 40);
  TGLabel *fL1 = new TGLabel(fHm0, "LE hist: ");
  fEdit1 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit1->SetToolTipText("bins min max");
  fEdit1->Resize(80, fEdit1->GetDefaultHeight());
  fHm0->AddFrame(fL1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));

  TGLabel *fL2 = new TGLabel(fHm0, "TOT hist: ");
  fEdit2 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit2->SetToolTipText("bins min max");
  fEdit2->Resize(80, fEdit2->GetDefaultHeight());
  fHm0->AddFrame(fL2, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit2, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));

  TGLabel *fL4 = new TGLabel(fHm0, "Time cut: ");
  fEdit3 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit3->SetToolTipText("min max");
  fEdit3->Resize(80, fEdit3->GetDefaultHeight());
  fHm0->AddFrame(fL4, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit3, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));


  TGLabel *fL5 = new TGLabel(fHm0, "Multiplicity cut: ");
  fEdit4 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit4->SetToolTipText("min max");
  fEdit4->Resize(80, fEdit4->GetDefaultHeight());
  fHm0->AddFrame(fL5, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit4, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
 
  TGLabel *fL6 = new TGLabel(fHm0, "TOF cut: ");
  fEdit5 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit5->SetToolTipText("min max");
  fEdit5->Resize(80, fEdit5->GetDefaultHeight());
  fHm0->AddFrame(fL6, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit5, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnHistory = new TGTextButton(fHm0, "M");
  fBtnHistory->Connect("Clicked()", "MyMainFrame", this, "DoHistory()");
  fHm0->AddFrame(fBtnHistory, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));
  
  fHm->AddFrame(fHm0, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX,5, 5, 5, 5));

  TGHorizontalFrame *fHm1 = new TGHorizontalFrame(fHm, 400, 40);

  TGLabel *fL3 = new TGLabel(fHm1, "Max in HP: ");
  fHm1->AddFrame(fL3, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));
  fHslider1 = new TGHSlider(fHm1, 183, kSlider1 | kScaleBoth, 0);
  fHslider1->Connect("PositionChanged(Int_t)", "MyMainFrame", this, "DoSlider(Int_t)");
  fHslider1->SetRange(0,50);
  fHm1->AddFrame(fHslider1, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fCheckBtn1  = new TGCheckButton(fHm1, new TGHotString("Local max"),        -1);
  fCheckBtn1->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked1()");
  fHm1->AddFrame(fCheckBtn1, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fCheckBtn2  = new TGCheckButton(fHm1, new TGHotString("3 sigma time cut                         "),        -1);
  fCheckBtn2->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked2()");
  fHm1->AddFrame(fCheckBtn2, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fCheckBtn4  = new TGCheckButton(fHm1, new TGHotString("show time map"),-1);
  fCheckBtn4->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked4()");
  fHm1->AddFrame(fCheckBtn4, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnExport = new TGTextButton(fHm1, "Export for the &blog");
  fBtnExport->Connect("Clicked()", "MyMainFrame", this, "DoExport()");
  fHm1->AddFrame(fBtnExport, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fHm->AddFrame(fHm1, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX,5, 5, 5, 5));

  TGHorizontalFrame *fHm2 = new TGHorizontalFrame(fHm, 400, 40);

  fCheckBtn3  = new TGCheckButton(fHm2, new TGHotString("Walk correction                         "),        -1);
  fCheckBtn3->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked3()");
  fHm2->AddFrame(fCheckBtn3, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnExportOffsets = new TGTextButton(fHm2, "Export &offsets");
  fBtnExportOffsets->Connect("Clicked()", "MyMainFrame", this, "DoExportOffsets()");
  fHm2->AddFrame(fBtnExportOffsets, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));
  
  fHm2->AddFrame(new TGLabel(fHm2, "Ch of interest: "), new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 3, 4));
  fNumber2 = new TGNumberEntry(fHm2, 0, 9,999, TGNumberFormat::kNESInteger,
			       TGNumberFormat::kNEANonNegative,
			       TGNumberFormat::kNELLimitMinMax,
			       0, glx_maxch);

  fNumber2->Connect("ValueSet(Long_t)", "MyMainFrame", this, "InterestChanged()");
  (fNumber2->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "InterestChanged()");
  fHm2->AddFrame(fNumber2, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnSavetPng = new TGTextButton(fHm2, "Save &current");
  fBtnSavetPng->Connect("Clicked()", "MyMainFrame", this, "DoSavePng()");
  fHm2->AddFrame(fBtnSavetPng, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));
  
  fHm->AddFrame(fHm2, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX,5, 5, 5, 5));


  AddFrame(fHm,   new TGLayoutHints(kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2));
  // =============== More frame ===================


  // status bar
  //Int_t parts[] = {45, 15, 10, 30};
  //fStatusBar = new TGStatusBar(this, 50, 10, kVerticalFrame);
  //fStatusBar->SetParts(parts, 4);
  //fStatusBar->Draw3DCorner(kFALSE);
  //AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
   
  // Create a horizontal frame containing two buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 400, 40);

  fLabel = new TGLabel(hframe, "Refference channel: ");
  hframe->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 3, 4));

  fNumber = new TGNumberEntry(hframe, 0, 9,999, TGNumberFormat::kNESInteger,
			      TGNumberFormat::kNEANonNegative, 
			      TGNumberFormat::kNELLimitMinMax,
			      0, glx_maxch);
  fNumber->SetNumber(0);
  hframe->AddFrame(fNumber, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 3, 4));
   
  TGTextButton *draw = new TGTextButton(hframe, "&Evaluate");
  draw->Connect("Clicked()", "MyMainFrame", this, "DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   
  fHProg3 = new TGHProgressBar(hframe, TGProgressBar::kStandard, 100);
  fHProg3->SetFillType(TGProgressBar::kBlockFill);
  hframe->AddFrame(fHProg3, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  pbar = fHProg3;


  fLabel1 = new TGLabel(hframe, "   Show: ");
  hframe->AddFrame(fLabel1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));


  fComboMode = new TGComboBox(hframe);
  hframe->AddFrame(fComboMode, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));
  fComboMode->AddEntry("Leading Edge", 0);
  fComboMode->AddEntry("Leading Edge All", 1);
  fComboMode->AddEntry("ToT", 2);
  fComboMode->AddEntry("ToT All", 3);
  fComboMode->AddEntry("TOF",12);
  if(gMode==1){
    fComboMode->AddEntry("Le vs. ToT", 4);
    fComboMode->AddEntry("Signal shape",10);
    fComboMode->AddEntry("Offset graph",11);
  }
  fComboMode->AddEntry("Multiplicity",5);
  fComboMode->AddEntry("Multiplicity All",6);
  fComboMode->AddEntry("Event Mult",8);
  fComboMode->AddEntry("Mult vs. Evt#",9);
  fComboMode->AddEntry("Time resolution",13);
  fComboMode->AddEntry("Channels",7);

  fComboMode->Resize(300, 20);
  fComboMode->Connect("Selected(Int_t)", "MyMainFrame", this, "updatePlot(Int_t)");
  

  TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
  exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  fBtnMore = new TGTextButton(hframe, "&More ");
  fBtnMore->Connect("Pressed()", "MyMainFrame", this, "DoMore()");
  hframe->AddFrame(fBtnMore, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  AddFrame(hframe, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2));

  fBackToHp=false;
  
  // Set a name to the main frame   
  SetWindowName("gdisplay " + ginFile);
  MapSubwindows();

  // Initialize the layout algorithm via Resize()
  Resize(GetDefaultSize());
  fComboMode->Resize(300, 20);
  HideFrame(fHm);

  // fEdit1->SetText("300 -100 200");
  // fEdit2->SetText("100 0 100");  

  fEdit1->SetText("1000 110 150");
  fEdit2->SetText("100 0 100");  

  
  if(ginFile.Contains("beam")) fEdit1->SetText("400 0 50");
  if(ginFile.Contains("aug2017")) fEdit1->SetText("400 -200 -150");
  if(ginFile.Contains("th_")) fEdit1->SetText("400 24 50");
  if(ginFile.Contains("pilas")) fEdit1->SetText("400 28 34");
  if(ginFile.Contains("hits.root")) fEdit1->SetText("400 0 50");
  
  fEdit3->SetText("0 0");
  fEdit4->SetText("0 0");
  fEdit5->SetText("0 0");
  fNumber->SetIntNumber(gTrigger);
  fCheckBtn2->SetState(kButtonUp);

  init();
  DoDraw();
  fComboMode->Select(7);

  glx_cdigi->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		 "exec3event(Int_t,Int_t,Int_t,TObject*)");


  cTime->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		 "exec4event(Int_t,Int_t,Int_t,TObject*)");
  
  MapWindow();  // Map main frame
}

MyMainFrame::~MyMainFrame(){
  Cleanup();
  delete fEcan;
  delete fTime;
}

void gdisplay(TString inFile= "pilasM.root", Int_t trigger=0, Int_t mode=0, TString path="", TString info="0",Int_t entries=50000, Int_t workers = 2){
  ginFile = inFile;
  gEntries=entries;
  gTrigger= trigger;
  gMode=mode;
  gPath=path;
  gInfo=info;
  glx_addInfo("Program = gdisplay");
  glx_addInfo("In file = " + ginFile);
  glx_addInfo("In path = " + gPath);
  gWorkers = workers;
  
  gMain = new MyMainFrame(gClient->GetRoot(), 800, 800);
}
