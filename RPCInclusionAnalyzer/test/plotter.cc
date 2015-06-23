//
//


void plotter(TString Input="csctf_mc.root")
{

    TFile *infile = TFile::Open(Input, "READ");

    vector<TString> All1Dhists;
    vector<TString> All2Dhists;

    //Add 1D histograms to the list to be plotted
    All1Dhists.push_back("dphi_rpc2_csc2");
    //Add 2D histograms to the appropriate list
    All2Dhists.push_back("dphi_csc1_rpc2_invpt");


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


    //Efficiency Histograms
    TCanvas *c = new TCanvas("c", "c", 700, 550);
    TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
    t1->Draw();
    t1->cd();
    t1->SetRightMargin(0.03);
    t1->SetLogx();
    
    TH1F* hist_1 = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/pt_turnon_threehit_all");
    TH1F* hist_2 = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/pt_turnon_threehit");
    heff = new TGraphAsymmErrors(hist_2, hist_1);

    heff->Draw();

    //Make plot look nice
    heff->GetYaxis()->SetTitle("Efficiency");
    heff->SetFillColor(kYellow);
    
    c->Modified();
    c->Update();
    /*TPaveStats *stats =  (TPaveStats*) heff->GetListOfFunctions()->FindObject("stats");
    stats->SetFillStyle(0);
    stats->SetName("");
    stats->SetX1NDC(.75);
    stats->SetY1NDC(.60);
    stats->SetX2NDC(.95);
    stats->SetY2NDC(.93);
    stats->SetTextColor(2);*/
    
    c->Update();
    
    //Save plots
    c->SaveAs("efficiency_threehit.png");
    c->SaveAs("efficiency_threehit.pdf");
    
    //delete pointers
    delete c;
    delete hist_1;
    delete hist_2;
    delete heff;

    TCanvas *c = new TCanvas("c", "c", 700, 550);
    TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
    t1->Draw();
    t1->cd();
    t1->SetRightMargin(0.03);
    t1->SetLogx();

    TH1F* hist_1 = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/pt_turnon_rpc_all");
    TH1F* hist_2 = (TH1F*) infile->Get("RPCInclusionAnalyzer/csctf/pt_turnon_rpc");
    heff = new TGraphAsymmErrors(hist_2, hist_1);

    heff->Draw();

    //Make plot look nice
    heff->GetYaxis()->SetTitle("Efficiency");
    heff->SetFillColor(kBlue);
    
    c->Modified();
    c->Update();    
    /*TPaveStats *stats =  (TPaveStats*) heff->GetListOfFunctions()->FindObject("stats");
    stats->SetFillStyle(0);
    stats->SetName("");
    stats->SetX1NDC(.75);
    stats->SetY1NDC(.60);
    stats->SetX2NDC(.95);
    stats->SetY2NDC(.93);
    stats->SetTextColor(2);*/
    c->Update();
    
    //Save plots
    c->SaveAs("efficiency_rpc.png");
    c->SaveAs("efficiency_rpc.pdf");
    
    //delete pointers
    delete c;
    delete hist_1;
    delete hist_2;
    delete heff;


    infile->Close();

}
