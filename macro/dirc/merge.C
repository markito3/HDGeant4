void merge(TString in="vol/bggen/bggen11/tree/060000/tree_060000_*.root", TString out="bggen_tree_all.root"){

  TChain ch("dirc");

  if(in.Contains("*")){
    TString tname = in;
    TString tname1 = in;
    TString tname2 = in;
    TString end = tname1.Remove(0,in.Last('*')+1);
    TString start = tname2.Remove(in.Last('*'));
    
    TString dirname("./");
    if(in.Last('/')!=-1) dirname= tname.Remove(in.Last('/')) + "/";

    const char *ext=".root";
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
	fname = file->GetName();
	if (!file->IsDirectory() && fname.EndsWith(ext)) {
	  TString path = dirname+fname;
	  TString substr = path.SubString(start);
	  if( substr.Length()>0 && path.EndsWith(end)){
	    std::cout<<"adding "<<path<<std::endl;
	    ch.Add(path);	    
	  }
	}
      }
    }
  }
  
  std::cout<<"merging into "<<out<<std::endl;
  ch.Merge(out);
}
