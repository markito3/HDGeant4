#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void producePDFs(TString infile="data/test.root",TString outfile="data/pdf_test.root"){
  if(!glx_initc(infile,1,"data/drawHP")) return;

  Int_t nPho[5]={0};

  // histograms for storing time pdfs for different particle species (5)
  TH1F *htime[5][glx_npixtot];
  TH1F* hnpho[5];
  for(Int_t j=0; j<5; j++){
    hnpho[j] = new TH1F(Form("npho_")+glx_names[j],"N photons",200, 0.,200.);
    for(Int_t i=0; i<glx_npixtot; i++){
      htime[j][i] = new TH1F(Form("h_%d_",i)+glx_names[j],"pdf; hit time [ns]; entries [#]", 1000, 0., 50.);
      //      cout<<Form("time_%d for ",i) + glx_names[j]<<endl;
    }
  }

  Double_t time;
  Int_t ch, pid, pmt, pix;
  DrcHit hit;
  Int_t nparticle[5]={0};
  for (Int_t e=0; e<glx_ch->GetEntries(); e++){
    glx_ch->GetEntry(e);
    for (Int_t t=0; t<glx_events->GetEntriesFast(); t++){
      glx_nextEventc(e,t,10);
      if(glx_event->GetParent()>0) continue;
      pid = glx_findPdgId(glx_event->GetPdg());
      //      cout <<"pid = "<<pid <<", hits - "<<glx_event->GetHitSize()<<endl;
      if(nparticle[pid]++ < 40000)continue;
      hnpho[pid]->Fill(glx_event->GetHitSize());
      for(Int_t h=0; h<glx_event->GetHitSize(); h++){
    	hit = glx_event->GetHit(h);
    	pmt = hit.GetPmtId();
        pix = hit.GetPixelId();
    	time = hit.GetLeadTime();
	ch = glx_getChNum(pmt, pix);
	nPho[pid]++;
	htime[pid][ch]->Fill(time);
      }
    }
  }

  TFile efile(outfile, "RECREATE");
  for(Int_t i=0; i<5; i++){
    cout<<"Npho in pix = "<<(Double_t)nPho[i]<<endl;;
    for(Int_t j=0; j<glx_npixtot; j++){
      if(htime[i][j]->GetEntries() > 0){
	htime[i][j]->Scale(1./(Double_t)nPho[i]);
	htime[i][j]->Write();
      }
    }
  }

  for(Int_t ipdg=0; ipdg < 5; ipdg++){
    cout<<"pdg = "<<ipdg<<" has "<<nPho[ipdg]<<" photons, events = "<<nparticle[ipdg] <<endl;
    hnpho[ipdg]->Write();
  }
  
  efile.Write();
  efile.Close();
}
