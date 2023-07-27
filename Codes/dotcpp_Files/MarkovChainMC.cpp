//ROOFit Parameter Estimation - {Cs,Cint,Phi0} using Hamiltonian Markov Chain Monte-Carlos

NOT FINISHED WHATSOEVER

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TVector3.h>
#include <TF1.h>
#include <RooMCStudy.h>



std::array<double,2> roof_CiPh2(int bins, double startt, double endd, int events, double DiluCon, double CSv, double Cintv, double Phi0v, double CLfrac, int loops)
{
    using namespace RooFit;
    int sorint = 0; //0 to fix CS, 1 to fix Cint, otherwie to fix neither
    
    
    double CL = 1.0;
    double tauKS = 0.08954, tauKL = 51.16; //in nano-seconds
    double MK = 497.614 /*in MeV/c^2*/,  dMK = 2.1890e-11 /*in MeV/c^2*/;

    double mom = 75000.0; //in MeV/c
    double c_light = 299.792; //in 10^6 m/s
    double tauK = 2.0/((1.0/tauKS)+(1.0/tauKL)), gamma = mom/MK;

    double lambdaL = gamma*c_light*tauKL, lambdaS = gamma*c_light*tauKS, lambda = gamma*c_light*tauK;

    ostringstream decayform, decayfit, testfit, testfit2;
    decayform<<CL<<"*exp(-x*"<<tauKS/tauKL<<") + CS*exp(-x*"<<tauKS/tauKS<<") + 2*exp(-x*"<<tauKS/tauK<<")*"<<DiluCon<<"*Cint*cos(Phi0-("<<dMK<<"*x/"<<tauKS<<"e-9))";                            
    //decayform<<CL<<"*exp(-x*"<<tauKS/tauKL<<")";
    //decayform<<"CS*exp(-x*"<<tauKS/tauKS<<")";
    //decayform<<CL<<"*exp(-x*"<<tauKS/tauKL<<") + CS*exp(-x*"<<tauKS/tauKS<<")";
    //decayform<<"2*exp(-x*"<<tauKS/tauK<<")*"<<DiluCon<<"*Cint*cos(Phi0-("<<dMK<<"*x/"<<tauKS<<"e-9))";
    //decayform<<CL<<"*exp(-x*"<<tauKS/tauKL<<") + 2*exp(-x*"<<tauKS/tauK<<")*"<<DiluCon<<"*Cint*cos(Phi0-("<<dMK<<"*x/"<<tauKS<<"e-9))";
    //decayform<<"CS*exp(-x*"<<tauKS/tauKS<<") + 2*exp(-x*"<<tauKS/tauK<<")*"<<DiluCon<<"*Cint*cos(Phi0-("<<dMK<<"*x/"<<tauKS<<"e-9))";

    
    //////////////////////////////////////////CLendd Evaluation////////////////////////////////////////////////////////
    ostringstream cldecayform;
    cldecayform<<CL<<"*exp(-x*"<<tauKS/tauKL<<")";
    TF1 *cldecayf = new TF1("cldecayf",cldecayform.str().c_str(),startt,endd);
    double CLendd = cldecayf->EvalPar(&endd);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    int i=0, j=0;

    RooRealVar x("x", "x", startt, endd);
    RooRealVar CS("CS", "CS", CSv, 0.00, 1.0);
    RooRealVar Cint("Cint", "Cint", Cintv, 0.00, 1.0);
    RooRealVar Phi0("Phi0", "Phi0", gRandom->Gaus(Phi0v,0.1), -1.58, 1.58);
    //Phi0.setConstant(true);
    RooGenericPdf decaypdf("decaypdf", "decaypdf", decayform.str().c_str(), RooArgSet(x, CS, Cint, Phi0));
    
    RooGaussian Phicon("Phicon","Phicon", Phi0, RooConst(Phi0v), RooConst(0.1/*From Precision Measurement Paper*/)) ;

    decayfit<<CL*CLfrac<<"*exp(-x*"<<tauKS/tauKL<<") + CS*exp(-x*"<<tauKS/tauKS<<") + 2*exp(-x*"<<tauKS/tauK<<")*"<<DiluCon<<"*Cint*cos(Phi0-("<<dMK<<"*x/"<<tauKS<<"e-9))";                                                       
    RooGenericPdf woaht("woaht", "woaht", decayfit.str().c_str(), RooArgSet(x, CS, Cint, Phi0));

    
    RooDataSet *data = decaypdf.generate(x, events);
    TH1 *roohist = data->createHistogram("roohist", x, Binning(bins,startt,endd));
    for(j=1;j<=bins;j++)
    {
        roohist->SetBinContent(j, (roohist->GetBinContent(j)-((CL-CLfrac)*roohist->GetBinContent(bins)*(CL/CLendd)*exp(-(j*(endd/bins))*tauKS/tauKL))));
        if(roohist->GetBinContent(j) < 0.0)
        {
            roohist->SetBinContent(j, 0.0);
        }
    }
    RooDataHist rdh("rdh", "rdh", x, roohist);

    
    TGraph csgr;
    TGraph cintgr;
    TGraph onlyint;
  
    for(i=0;i<(2*loops);i++)
    {
        sorint = (i%2);
        if(i==0)
        {
            CS.setVal(0.21);
            Cint.setVal(0.06);
        }
        if(sorint==1)
        {
            CS.setConstant(true);
            Cint.setConstant(false);
            //std::unique_ptr<RooFitResult> fitResult{decaypdf.fitTo(rdh, ExternalConstraints(Phicon), Save(), PrintLevel(-1))};
            std::unique_ptr<RooFitResult> fitResult{woaht.fitTo(rdh, ExternalConstraints(Phicon), PrintLevel(-1), PrintEvalErrors(-1), Warnings(false))};
            //std::unique_ptr<RooFitResult> fitResult{woaht.fitTo(*data, ExternalConstraints(Phicon), Save(), PrintLevel(-1))};
            cintgr.SetPoint((i/2), (i/2), (Cint.getValV()-Cintv)/Cintv);
            onlyint.SetPoint((i/2), (i/2), Cint.getValV());
        }
        if(sorint==0)
        {
            CS.setConstant(false);
            Cint.setConstant(true);
            //std::unique_ptr<RooFitResult> fitResult{decaypdf.fitTo(rdh, ExternalConstraints(Phicon), Save(), PrintLevel(-1))};
            std::unique_ptr<RooFitResult> fitResult{woaht.fitTo(rdh, ExternalConstraints(Phicon), PrintLevel(-1), PrintEvalErrors(-1), Warnings(false))};
            //std::unique_ptr<RooFitResult> fitResult{woaht.fitTo(*data, ExternalConstraints(Phicon), Save(), PrintLevel(-1))};
            csgr.SetPoint((i/2), (i/2), (CS.getValV()-CSv)/CSv);
        }
    }  
    
    /*
    TCanvas c("c","Canvase",1920,1080);
    
    cout<<"\n\n";
    csgr.SetTitle("Difference Between Fitted CS and True CS vs Recursion;Loop Number;(Fit(CS)-0.43)/0.43");
    csgr.Draw();
    gPad->Print("Recursive_CS.png");
    c.Clear();
    cintgr.SetTitle("Difference Between Fitted Cint and True Cint vs Recursion;Loop Number;(Fit(Cint)-0.12)/0.12");
    cintgr.Draw();
    gPad->Print("Recursive_Cint.png");
    c.Clear();
    onlyint.SetTitle("Fitted Cint vs Recursion;Loop Number;Fit(Cint)");
    onlyint.Draw();
    gPad->Print("Recursive_Cint_ActualVal.png");
    c.Clear();
    
    
        
    
    testfit<<CL*CLfrac<<"*exp(-x*"<<tauKS/tauKL<<") + "<<CSv<<"*exp(-x*"<<tauKS/tauKS<<") + 2*exp(-x*"<<tauKS/tauK<<")*"<<DiluCon<<"*"<<Cintv<<"*cos("<<Phi0v<<"-("<<dMK<<"*x/"<<tauKS<<"e-9))";                            
    RooGenericPdf testw("testw", "testw", testfit.str().c_str(), RooArgSet(x));

    testfit2<<CL*CLfrac<<"*exp(-x*"<<tauKS/tauKL<<") + "<<CSv<<"*exp(-x*"<<tauKS/tauKS<<")";                            
    RooGenericPdf testw2("testw", "testw", testfit2.str().c_str(), RooArgSet(x));
 
    
    RooPlot *xframe = x.frame(Title("K0 Decay PDF"));
    rdh.plotOn(xframe);
    //data->plotOn(xframe);
    woaht.plotOn(xframe, LineColor(kBlue));
    testw.plotOn(xframe, LineColor(kRed), LineStyle(kDashed));
    testw2.plotOn(xframe, LineColor(kGreen));
    xframe->Draw();
    gPad->Print("ROOHIST.png");
    */
    
    
    
    std::array<double,2> carray = {CS.getValV(),Cint.getValV()};
    return carray;
}









void RecOp()
{
	int bins = 60, noh = 10, events = (300000), loops = 50, i;
	double DiluCon = 0.285714286;//0.285714286;   //Dilucon = 2.0/7.0;  N(K0)/N(K0bar) = 1.8
	double CS = 0.43;//0.43
	double Cint = 0.12;//0.12;
	double Phi0 = 0.201357921;//0.201357921;

	cout<<"\nNumber of Histograms (noh) = ";
	cin>>noh;
	cout<<"\n\n";

	TH1D rechistcs("rechistcs","Recursive Optimization: C_{S};C_{S};Count",40,-0.0,1.0);
	TH1D rechistcint("rechistcint","Recursive Optimization: C_{int};C_{int};Count",40,-0.0,1.0);

	double startt = 0, endd = 6;
	double CLfrac = 1.0;
	std::array<double,2> carray;
	for(i=0;i<noh;i++)
	{
	carray = roof_CiPh2(bins, startt, endd, events, DiluCon, CS, Cint, Phi0, CLfrac, loops);
	rechistcs.Fill(carray[0]);
	rechistcint.Fill(carray[1]);
	cout<<endl<<(double(i)/double(noh))*100<<"% Done";
	}

	cout<<"\n\n";
	gStyle->SetOptStat(1111111);
	rechistcs.Draw();
	gPad->Print("RecursiveHisto_CS.png");
	rechistcint.Draw();
	gPad->Print("RecursiveHisto_Cint.png");
}

