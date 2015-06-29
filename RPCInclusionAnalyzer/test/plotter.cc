//
//


void plotter(TString Input="csctf_mc.root")
{

    TFile *infile = TFile::Open(Input, "READ");

    vector<TString> All1Dhists;
    vector<TString> All2Dhists;
    vector<TString> AllEffhists;

    //Add 1D histograms to the list to be plotted
    All1Dhists.push_back("dphi_csc1_csc2_all");
    All1Dhists.push_back("dphi_csc1_rpc2");
    All1Dhists.push_back("dphi_csc1_csc2");
    All1Dhists.push_back("dphi_rpc2_csc3");
    All1Dhists.push_back("dphi_rpc2_csc2");
    //All1Dhists.push_back("");



    //Add 2D histograms to the appropriate list
    All2Dhists.push_back("dphi_csc1_rpc2_invpt");


    AllEffhists.push_back("pt_turnon_threehit");
    AllEffhists.push_back("pt_turnon_rpc");

    //Set the plot style
    //gStyle->SetOptStat(0);
    gStyle->SetOptStat("emruoi");


    //Loop over the 1D histograms
    for (size_t ihist = 0; ihist < All1Dhists.size(); ihist++) {

	//Read in the histogram
        TH1F* hist_2 = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/"+All1Dhists[ihist]);

	//Check if histogram actually exists
        if (hist_2 == NULL) {
            cout << "hist is NULL!" << endl;
            continue;
        }

	//Set up the canvas
        TCanvas *c = new TCanvas("c", "c", 700, 550);
        //TCanvas *c = new TCanvas("c", "c", 900, 550);
        TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
        t1->Draw();
        t1->cd();
        //t1->SetBottomMargin(0.3);
        t1->SetRightMargin(0.03);
        t1->SetLogy(1);
        //c->Divide(1,2);

	//Draw the histogram on the canvas
        hist_2->Draw("hist");

	//Make plot look nice
        //hx_->GetXaxis()->SetTitle("Vertices");
        hist_2->GetYaxis()->SetTitle("Events");
        hist_2->SetFillColor(kYellow);
        if(All1Dhists[ihist].Contains("MEp") || All1Dhists[ihist].Contains("_p_")) hist_2->SetFillColor(kGreen-7);
        if(All1Dhists[ihist].Contains("MEm") || All1Dhists[ihist].Contains("_m_")) hist_2->SetFillColor(kOrange-3);


        c->Modified();
        c->Update();
        TPaveStats *stats =  (TPaveStats*) hist_2->GetListOfFunctions()->FindObject("stats");
        stats->SetFillStyle(0);
        stats->SetName("");
        stats->SetX1NDC(.75);
        stats->SetY1NDC(.60);
        stats->SetX2NDC(.95);
        stats->SetY2NDC(.93);
        stats->SetTextColor(2);

        c->Update();

	//Save plots
        c->SaveAs(All1Dhists[ihist]+".png");
        c->SaveAs(All1Dhists[ihist]+".pdf");

	//delete pointers
        delete c;
        delete hist_2;
    }



    //Loop over the 2D histograms
    for(size_t ihist=0; ihist<All2Dhists.size(); ihist++) {

        TH2F* hist_temp = (TH2F*) infile->Get("RPCInclusionAnalyzer/csctf/"+All2Dhists[ihist]);
        if(hist_temp==NULL) {
            cout << "hist is NULL!" << endl;
            continue;
        }

        TCanvas *c = new TCanvas("c", "c", 700, 550);
        TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.0);
        t1->Draw();
        t1->cd();
        //t1->SetBottomMargin(0.3);
        t1->SetRightMargin(0.12);
        //c->Divide(1,2);
        hist_temp->Draw("col Z");
        //hx_->GetXaxis()->SetTitle("Vertices");
        //hx_->GetYaxis()->SetTitle("<E_{Y}^{miss}>");


        c->Modified();
        c->Update();
        TPaveStats *stats =  (TPaveStats*) hist_temp->GetListOfFunctions()->FindObject("stats");
	stats->SetFillStyle(0);
        stats->SetName("");
        //stats->SetFillColor();
        stats->SetX1NDC(.70);
        stats->SetY1NDC(.67);
        stats->SetX2NDC(.85);
        stats->SetY2NDC(.93);
        stats->SetTextColor(2);




        c->SaveAs(All2Dhists[ihist]+".png");
        c->SaveAs(All2Dhists[ihist]+".pdf");

        delete c;
        delete hist_temp;

    }


    //Loop over the 1D Eff histograms
    for (size_t ihist = 0; ihist < AllEffhists.size(); ihist++) {

	gStyle->SetOptStat(0);
	//gROOT->ProcessLine(".L TurnOnFit.C");

	//Read in the histogram
        TH1F* hist_denominator = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/"+AllEffhists[ihist]+"_all");
	TH1F* hist_numerator = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/"+AllEffhists[ihist]);

	//Check if histogram actually exists
        if (hist_denominator == NULL || hist_numerator == NULL) {
            cout << "hist is NULL!" << endl;
            continue;
        }

	//Set up the canvas
        TCanvas *c = new TCanvas("c", "c", 700, 550);
        //TCanvas *c = new TCanvas("c", "c", 900, 550);
        TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
        t1->Draw();
        t1->cd();
        //t1->SetBottomMargin(0.3);
        t1->SetRightMargin(0.03);
        t1->SetLogx(1);
        //c->Divide(1,2);

	//Draw the histogram on the canvas
	hist_numerator->Divide(hist_denominator);
        hist_numerator->Draw("EP");

	hist_numerator->SetMaximum(1.1);
	hist_numerator->SetMinimum(0);
	//Make plot look nice
        //hx_->GetXaxis()->SetTitle("Vertices");
        hist_numerator->GetYaxis()->SetTitle("Efficiency");
        hist_numerator->SetFillColor(kYellow);


  	TF1 *fit = new TF1("fit",turnon_func,10.,18.,3);
  	fit->SetParameters(13.0,1.0,0.75);
  	hist_numerator->Fit("fit");


        c->Modified();
        c->Update();
        TPaveStats *stats =  (TPaveStats*) hist_numerator->GetListOfFunctions()->FindObject("stats");
	stats->SetFillStyle(0);
        stats->SetName("");
        //stats->SetFillColor();
        stats->SetX1NDC(.70);
        stats->SetY1NDC(.27);
        stats->SetX2NDC(.85);
        stats->SetY2NDC(.53);
        stats->SetTextColor(2);



	//Save plots
        c->SaveAs("Eff_"+AllEffhists[ihist]+".png");
        c->SaveAs("Eff_"+AllEffhists[ihist]+".pdf");

	//delete pointers
        delete c;
        delete hist_denominator;
	delete hist_numerator;
    }




    infile->Close();

}



void TurnOnFit(TH1 *h1, TH1 *h2, TH1 *h) {
  // Calculate the efficiency bin by bin as the ratio of two histograms.
  // Errors on the efficiency are calculated properly.

  int nbins = h->GetNbinsX();

  // Check histogram compatibility
  if (nbins != h1->GetNbinsX() || nbins != h2->GetNbinsX()) {
    std::cout << "Attempt to divide histograms with different number of bins" << std::endl;
    return;
  }

  // Issue a Warning if histogram limits are different
  if (h->GetXaxis()->GetXmin() != h1->GetXaxis()->GetXmin() ||
      h->GetXaxis()->GetXmax() != h1->GetXaxis()->GetXmax() ){
    std::cout << "Attempt to divide histograms with different axis limits" << std::endl;
  }

  if (h->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin() ||
      h->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax() ){
    std::cout << "Attempt to divide histograms with different axis limits" << std::endl;
  }

  //if (h->GetDimension() < 2) nbinsy = -1;

  // Loop on bins (including underflows/overflows)
  int bin;
  double b1,b2,effi,erro;
  for (bin=0;bin<=nbins+1;bin++) {

    b1  = h1->GetBinContent(bin);
    b2  = h2->GetBinContent(bin);
    //std::cout << "b1: " << b1 << endl;
    //std::cout << "b2: " << b2 << endl;

    if ( b2 ) {
      effi = b1/b2;
      erro = TMath::Sqrt((1-b1/b2)/b2*(b1/b2));
      //std::cout << "effi: " << effi << endl;
      //std::cout << "erro: " << erro << endl;
    } else {
      effi = 0;
      erro = 0;
    }
    h->SetBinContent(bin,effi);
    h->SetBinError(bin,erro);
  }
}

Double_t turnon_func(Double_t *x, Double_t *par)
{
  double halfpoint = par[0];
  double slope = par[1];
  double plateau = par[2];

  //double offset = par[3];
  //double plateau = 1.0;
  double offset = 0;

  double pt = TMath::Max(x[0],0.000001);

   double arg = 0;
   //cout << pt <<", "<< halfpoint <<", " << slope <<endl;
   arg = (pt - halfpoint)/(TMath::Sqrt(pt)*slope);
   double fitval = offset+0.5*plateau*(1+TMath::Erf(arg));
   return fitval;
}

Double_t turnon_func2(Double_t *x, Double_t *par)
{
  return turnon_func(x,par) + turnon_func(x,&par[3]);
}


