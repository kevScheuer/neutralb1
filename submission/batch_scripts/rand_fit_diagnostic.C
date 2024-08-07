/*Evaluate the success of the AmpTools randomized fit procedure with various metrics

Plot the sum of squares of error (with respect to the best fit) divided by the MINUIT
uncertainty squared for every production and fit parameter, as a function of the
difference in likelihood between each fit and the best fit. Plotted on log-log scale

NOTE: Once ROOT is updated so that TScatter is available, this script will become much
simpler. Currently has a lot of hacky ways to get around not being able to natively
make a scatter plot
*/
#include <cmath>
#include <dirent.h>
#include <string>

#include "glueXstyle.C"
#include "IUAmpTools/FitResults.h"
#include "TString.h"

// forward declarations
std::vector<std::string> GetFileList(std::string path, std::string reaction);

void rand_fit_diagnostic(std::string path, std::string reaction = "omegapi")
{
    gluex_style();
    gStyle->SetOptStat(0);
    auto file_list = GetFileList(path, reaction);

    // load the best fit result first to make comparisons to
    FitResults best_fit(path + "/best.fit");
    if (!best_fit.valid())
    {
        std::cout << "Invalid fit results in file: " << path + "best.fit" << "\n";
        exit(1);
    }
    // get number of fit parameters by assuming any errors=0 means the param is fixed
    int num_free_params = 0;
    for (std::string name : best_fit.parNameList())
    {
        if (best_fit.parError(name) != 0)
            num_free_params += 1;
    }

    double best_likelihood = best_fit.likelihood(); // min likelihood to compare to
    double convergence_rate = 0;                    // tracks percentage fits that converged successfully

    // params for making the plot later
    float p[2];
    float l_max = 0;
    float chi2_max = 0;

    TNtuple *tp = new TNtuple("tp", "tp", "delta_l:chi2");

    for (std::string file : file_list)
    {
        FitResults rand_fit(file);
        if (!rand_fit.valid())
        {
            std::cout << "Invalid fit results in file: " << file << "\n";
            exit(1);
        }

        // skip fits that failed to converge
        if (rand_fit.lastMinuitCommandStatus() == 4)
            continue;

        convergence_rate += 1;
        float chi2 = 0; // tracks sum_{params} ((best_fit - rand_fit)/(best_error))^2

        // best_likelihood is always smallest number, so subtract it from rand_fit to
        // get a nice positive # for plotting
        float delta_l = rand_fit.likelihood() - best_likelihood;
        p[0] = delta_l;

        if (delta_l > l_max)
            l_max = delta_l;

        for (unsigned int i = 0; i < rand_fit.parNameList().size(); i++)
        {
            std::string rand_par_name = rand_fit.parNameList()[i];
            std::string best_par_name = best_fit.parNameList()[i];

            if (rand_par_name != best_par_name)
            {
                std::cout << "best and randomized fit have different param set. Exiting"
                          << "\n";
                exit(1);
            }

            double rand_par_val = rand_fit.parValue(rand_par_name);
            double best_par_val = best_fit.parValue(best_par_name);
            double best_par_err = best_fit.parError(best_par_name);

            if (best_par_err == 0) // avoid discontinuity from fixed params
                continue;

            chi2 += std::pow((best_par_val - rand_par_val) / best_par_err, 2);
        }

        float reduced_chi2 = chi2 / num_free_params;
        p[1] = reduced_chi2;
        tp->Fill(p);        

        if (reduced_chi2 > chi2_max)
            chi2_max = reduced_chi2;
    }

    // convergence rate has count of successful fits, so divide by total # of rand fits
    convergence_rate = convergence_rate / file_list.size();

    // avoid plotting issues if all rand fits at best value
    if (l_max == 0.0)
        l_max = 0.001;
    if (chi2_max == 0.0)
        chi2_max = 0.001;

    auto c1 = new TCanvas();
    c1->SetLogx(1);
    c1->SetLogy(1);
    TH2F *frame = new TH2F("frame", "frame", 100, 0, l_max, 100, 0, chi2_max);
    frame->GetXaxis()->SetTitle("\\Delta (-2ln L)");
    frame->GetYaxis()->SetTitle("\\chi^{2}");
    frame->Draw();

    tp->SetMarkerSize(1.5);
    tp->SetMarkerStyle(20);
    tp->Draw("chi2:delta_l", "", "same");

    TLatex t(.6, .96, TString::Format("%.2f%% Convergence", convergence_rate * 100).Data());
    t.SetNDC(kTRUE);
    t.Draw();

    gPad->Modified();
    gPad->Update();
    c1->SaveAs((path + "rand_fit_diagnostic.pdf").c_str());
}

// Returns unordered list of all randomized .fit files with {path}/rand/. Assumes
// they're of the form "{reaction}_#.fit"
std::vector<std::string> GetFileList(std::string path, std::string reaction)
{
    std::vector<std::string> file_list;
    DIR *dir;
    struct dirent *ent;

    dir = opendir((path + "rand/").c_str());
    if (dir == NULL)
    {
        cout << "Directory " << path << " could not be opened" << "\n";
        exit(1);
    }
    while ((ent = readdir(dir)))
    {
        std::string file = ent->d_name;

        if (file.find(reaction) != std::string::npos &&
            file.find(".fit") != std::string::npos)
        {
            file_list.push_back(path + "rand/" + file);
        }
    }

    return file_list;
}