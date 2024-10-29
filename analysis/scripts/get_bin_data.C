/*Aggregate binned data from mass independent fits into one csv

Currently when mass independent fits are submitted, they source pre-cut ROOT tress for
the fitting. In order to plot the fit results of the entire mass range and compare the
fit to data, we need to gather all the pre-cut ROOT trees into one csv. We of course
want the total events in each bin, but also the mass range, mean, and rms of the bin.
*/
#include <dirent.h>
#include <fstream>
#include <string>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

// forward declaration
void GetDirList(std::vector<std::string> &dir_list, std::string parent_dir, std::string target_file);
float GetLastNumOfString(std::string str);
bool NaturalComp(std::string s1, std::string s2);

/*
    Get data from every mass bin and store results in csv file
    Args:
        parent_dir: searches all subdirectories of this directory
        target_file: file name containing ROOT tree of interest
        csv_name (optional): output csv name. Defaults to
*/
void get_bin_data(std::string parent_dir, std::string target_file, std::string csv_name)
{
    // create csv file
    ofstream csv;
    csv.open(csv_name);

    // write header of bins to csv
    csv << "mass_low_edge,mass_high_edge,t_low_edge,t_high_edge,"
        << "bin_contents,bin_error,mass_mean,mass_rms,"
        << "t_mean,t_rms\n";

    std::vector<std::string> dir_list;
    GetDirList(dir_list, parent_dir, target_file);

    if (dir_list.size() == 0)
    {
        cout << "File not found in any directories" << "\n";
        exit(1);
    }

    // sort directories list by last integer in path. Ensures each row is ordered by
    // mass bin (typically the last bin in directory path)
    std::sort(dir_list.begin(), dir_list.end(), NaturalComp);

    // since these directories have all the cuts applied already, we just need the data
    for (auto dir : dir_list)
    {
        std::cout << "Processing " << (dir + target_file).c_str() << "\n";
        std::unique_ptr<TFile> f(TFile::Open((dir + target_file).c_str()));
        TTree *tree = f->Get<TTree>("kin");
        if (!tree)
        {
            cout << "Tree DNE" << "\n";
            exit(1);
        }
        TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1800);
        tree->Draw("M4Pi>>h_mass(100)", "Weight");
        TH1F *h_mass = (TH1F *)gPad->GetPrimitive("h_mass");

        tree->Draw("t>>h_t(100)", "Weight");
        TH1F *h_t = (TH1F *)gPad->GetPrimitive("h_t");

        // Get mass bin edges
        // NOTE: this assumes mass binning is neatly rounded to 2nd decimal (10s of MeV)
        int nbins = h_mass->GetNbinsX();
        double low_mass = roundf(h_mass->GetXaxis()->GetBinLowEdge(1) * 100) / 100;
        double high_mass = roundf(h_mass->GetXaxis()->GetBinUpEdge(nbins - 1) * 100) / 100;
        double bin_contents_error;
        double bin_contents = h_mass->IntegralAndError(0, nbins+1, bin_contents_error);
        double mean_mass = h_mass->GetMean();
        double rms_mass = h_mass->GetRMS();

        // Get t bin edges
        // NOTE: similar idea, but rounded to 1st decimal
        nbins = h_t->GetNbinsX();
        double low_t = roundf(h_t->GetXaxis()->GetBinLowEdge(1) * 10) / 10;
        double high_t = roundf(h_t->GetXaxis()->GetBinUpEdge(nbins - 1) * 10) / 10;
        double mean_t = h_t->GetMean();
        double rms_t = h_t->GetRMS();

        csv << low_mass << "," << high_mass << "," << low_t << "," << high_t << ","
            << bin_contents << "," << bin_contents_error << "," 
            << mean_mass << "," << rms_mass << ","
            << mean_t << "," << rms_t << "\n";

        delete c1;
    }

    return;
}

void GetDirList(std::vector<std::string> &dir_list, std::string parent_dir, std::string target_file)
{
    DIR *dir;
    struct dirent *ent;

    dir = opendir(parent_dir.c_str());

    if (dir == NULL)
    {
        cout << "Directory " << parent_dir << " could not be opened"
             << "\n";
        exit(1);
    }

    // read directory (directories) and add ROOT info to csv
    while ((ent = readdir(dir)))
    {
        if (ent->d_name[0] == '.') // skip reading current or parent dir again
        {
            continue;
        }
        // recursively search the subdirs
        if (ent->d_type == DT_DIR)
        {
            std::string path = parent_dir + ent->d_name + '/';
            GetDirList(dir_list, path, target_file);
        }

        std::string file = ent->d_name;

        if (file != target_file)
        {
            continue;
        }
        dir_list.push_back(parent_dir);
    }
}

// Returns last integer in string. If none is found, returns maximum int
// allowed. This ensures in "NaturalComp" that files not indexed by an
// integer at the end of the csv
float GetLastNumOfString(std::string str)
{
    float last_num = 2147483647;
    size_t begin = 0, end = 0;

    std::string numbers = ".0123456789";

    begin = str.find_first_of(numbers);

    while (begin != std::string::npos)
    {
        // avoid cases where is just single "." in substring
        if (str.at(begin) == '.' &&
            numbers.find(str.at(begin + 1)) == std::string::npos)
        {
            str = str.substr(begin + 1, str.length() - (begin + 1));
            begin = str.find_first_of(numbers);
            continue;
        }

        end = str.find_first_not_of(numbers, begin);
        std::string num = str.substr(begin, end - begin);
        last_num = std::atof(num.c_str());
        str = str.substr(end, str.length() - end);
        begin = str.find_first_of(numbers);
    }

    return last_num;
}

bool NaturalComp(std::string s1, std::string s2)
{
    return (GetLastNumOfString(s1) < GetLastNumOfString(s2));
}