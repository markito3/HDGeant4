#include "glxtools.C"
#include "../../../halld_recon/src/plugins/Analysis/lut_dirc/DrcLutNode.h"

//lut/lut_11/lut_all_avr.root
void drawLut(TString infile = "lut/lut_04/3/lut_3_avr.root"){
  gInterpreter->GenerateDictionary("vector<TVector3>","TVector3.h"); 

  if(infile=="") return;  
  TFile* f = new TFile(infile);

  TTree *t=(TTree *) f->Get("lut_dirc") ;
  TClonesArray* fLut[48];
  for(Int_t l=0; l<48; l++){
    fLut[l] = new TClonesArray("DrcLutNode");
    t->SetBranchAddress(Form("LUT_%d",l),&fLut[l]); 
  }
  glx_initDigi();
  glx_setRootPalette(1);    
  glx_savepath="data/drawLut";
  
  int nsum;
  TVector3 dir,dir2,osum,sum;
  Double_t angle,minangle,time,sumt;
  Long64_t pathid;
  DrcLutNode *node;
  t->GetEntry(0);
  
  for(size_t l=3; l<4; l++){
    for (int inode=0; inode<fLut[l]->GetEntriesFast(); inode++){
      if(inode%1000==0) cout<<"Node # "<<inode << "  L "<<l<<endl;
      node= (DrcLutNode*) fLut[l]->At(inode);
      int pmt=inode/64;
      int pix=inode%64;
      int size = node->Entries();      
      if(size<1) continue;
      std::cout<<l<<" size "<<size<<std::endl;
      for(int i=0; i<size; i++){
	dir = node->GetEntry(i);
	time = node->GetTime(i);
	pathid = node->GetPathId(i);
	if(time>10) continue;
	
	glx_hdigi[pmt]->Fill(pix%8, pix/8);
      } 
    }
  }
  glx_drawDigi("",0);
  glx_canvasAdd(glx_cdigi);
}
