#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
Int_t npeaks = 2;

Double_t XMIN=350;
Double_t XMAX=700;
Int_t    SPECTRUM_NBINS=1024;

Double_t fpeaks(Double_t *x, Double_t *par) 
{
	Double_t result = par[0] ;
	for (Int_t p = 0; p < npeaks; p++) 
	{
		Double_t norm = par[3 * p + 1]; // "height" or "area"
		Double_t mean = par[3 * p + 2];
		Double_t sigma = par[3 * p + 3];
		result += norm * TMath::Gaus(x[0], mean, sigma);
	}
	return result;
}
void findMaxBinX(TH1F *h, Double_t min, Double_t max,  Double_t &maxX, Double_t &maxY)
{
	if(!h){
		maxX=-1;
		maxY=-1;
		return;
	}

	Double_t y;
	maxY=0;

	for(int i=1;i<=SPECTRUM_NBINS;i++)
	{
		Double_t x=h->GetXaxis()->GetBinCenter(i);
		if(x<min || x>max)continue;
		y=h->GetBinContent(i);
		if(y>maxY)
		{
			maxX=x;
			maxY=y;
		}
	}
	return;
}
void peaks(TFile *_file0,  TFile *_file1,  Int_t detector, Int_t channel, Double_t *par)
{
	TH1F *h = (TH1F *)_file0->Get(Form("hcal_%d_%d",detector,channel));
	if(!h)
	{
		return;
	}
	Int_t num_bins=h->GetNbinsX();
	// Generate n peaks at random
	par[0] = 0;
	Double_t maxX, maxY;
	findMaxBinX(h, XMIN, XMAX,  maxX, maxY);
	XMIN=maxX-20;

	Double_t maxX2, maxY2;

	findMaxBinX(h, maxX+40, XMAX,  maxX2, maxY2);

	cout<<"MaxX:"<<maxX<<" "<<maxY<<endl;
	par[1]=maxY;
	par[2]=maxX;
	par[3]=5.5;

	par[4]=0.3*maxY;
	par[5]=maxX+100;
	par[6]=10;
	TF1 *f = new TF1("f", fpeaks, XMIN, XMAX, 1 + 3 * npeaks);
	f->SetParNames("offset", "norm0","mean0","sigma0","norm1","mean1","sigma1");
	f->SetNpx(1000);
	f->SetParameters(par);
	TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1000, 900);
	h->Draw();
	h->Fit(f, "R", "", XMIN, XMAX);
	if(_file1)
	{
		_file1->cd();
		c1->Write(Form("c_%d_%d",detector, channel));
		h->Write(Form("h2_%d_%d",detector, channel));
	}
	f->GetParameters(&par[0]);
}
void find_peaks(TString filename="spectra.root", TString filename_out="test.root")
{
	TFile *fin=new TFile(filename);
	TFile *fout=new TFile(filename_out,"recreate");
	Double_t par[10]={0};
	cout<<"detector, pixel, 35 keV peak, sigma, 80 keV, sigma, Emin, Emax, factor"<<endl;
	cout<<"Output:"<<filename_out<<endl;
	TH2F hcal("hcal","calibration factors", 32,0,32,12,0,12);
	for(int i=0;i<32;i++)
		for(int j=0;j<12;j++)
		{
			peaks(fin, fout,i, j,par);
			Double_t factor= (par[5]-par[2])/(81-35.);
			if(par[5]>par[2])
			{
				cout<<i<<","<<j<<","<<par[2]<<","<<par[3]<<","<<par[5]<<","<<par[6]
					<< ","<<350/factor<<","<<1000./factor<<","<<factor<<endl;
				hcal.SetBinContent(i+1,j+1,factor);
			}
		}
	fout->cd();
	hcal.Write();

}
