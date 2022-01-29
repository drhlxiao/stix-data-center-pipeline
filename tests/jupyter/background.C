void background() {
   Int_t i;
   const Int_t nbins = 1024;
   Double_t xmin     = 0;
   Double_t xmax     = nbins;
   Double_t source[nbins];
   Double_t back_data[nbins];
   gROOT->ForceStyle();
   TH1F *d    = new TH1F("d","",nbins,xmin,xmax);
   TString file = "expspec.root";
   TFile *f     = new TFile(file.Data());
   TH1F *back = (TH1F*) f->Get("hexp");
   back->SetTitle("Estimation of background with decreasing window");
   back->GetXaxis()->SetRange(1,nbins);
   back->Draw("L");
   TSpectrum *s = new TSpectrum();
   TF1 fbkg("fbkg","0.5e3*exp(-x/400)",0,nbins);
   for (i = 0; i < nbins; i++) {
	   back->SetBinContent(i+1, back->GetBinContent(i+1)+fbkg.Eval(i)+gRandom->Uniform(10));
	   source[i]=back->GetBinContent(i + 1);
	   back_data[i]=back->GetBinContent(i + 1);
   }
   // Estimate the background
   s->Background(source,nbins,6,TSpectrum::kBackDecreasingWindow,
                 TSpectrum::kBackOrder2,kFALSE,
                 TSpectrum::kBackSmoothing3,kFALSE);
   // Draw the estimated background
   cout<<TSpectrum::kBackDecreasingWindow<<" "<<TSpectrum::kBackOrder2<<" "<<
                 TSpectrum::kBackSmoothing3<<endl;
   for (i = 0; i < nbins; i++) d->SetBinContent(i + 1,source[i]);
   d->SetLineColor(kRed);
   d->Draw("SAME L");

   TFile fo("exspc_out.root","recreate");
   fo.cd();
   d->Write("hbkg");
   TH1F *hsubbkg=(TH1F*)back->Clone();
   hsubbkg->Reset();

   for (i = 0; i < nbins; i++) {
	   hsubbkg->SetBinContent(i + 1,back_data[i]-source[i]);
   }
   hsubbkg->Write("hsubbkg");
   TCanvas c2;
   hsubbkg->Draw();
   back->Write();
   fo.Close();
   


}
