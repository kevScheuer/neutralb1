/**
 * @file check_systematic.cc
 * @author Kevin Scheuer
 * @brief Ensure systematic variation does not exceed 10% change in statistics per bin
 *
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

// Names of cut branches that may appear in the systematic files.
static const std::vector<std::string> CUT_BRANCH_NAMES = {
    "unusedE", "unusedTracks", "chi2", "PzCMrecoilPi",
    "shQualityP2a", "shQualityP2b", "shQualityP3a", "shQualityP3b"
};

/**
 * @brief In-memory representation of events loaded from a set of ROOT files.
 *
 * Each element of the parallel vectors corresponds to a single event.
 * For background events the @p weight is the sideband weight (stored negative
 * so that @p count_events can simply sum all weights); for signal events it
 * is +1.0.
 */
struct EventCollection
{
    std::vector<Double_t> M4Pi;
    std::vector<Double_t> t;
    std::vector<Double_t> weight;                              ///< +1.0 for signal, -w for background
    std::map<std::string, std::vector<Double_t>> cut_vars;     ///< keyed by branch name
};

// forward declarations
EventCollection load_events(
    const std::vector<TString> &signal_files,
    const std::vector<TString> &bkgd_files,
    bool load_cut_vars = false);

double count_events(
    const EventCollection &events,
    float t_min, float t_max,
    float mass_min, float mass_max,
    const std::map<TString, std::pair<float, float>> &cuts = {});

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <systematic_variation>\n";
        return 1;
    }
    std::string systematic = argv[1];

    // -----------------------------------------------------------------------
    // Build file lists
    // -----------------------------------------------------------------------
    TString nominal_data_dir    = "/w/halld-scshelf2101/kscheuer/neutralb1/data/FSRoot/GlueX/";
    TString systematic_data_dir = "/w/halld-scshelf2101/kscheuer/neutralb1/data/systematics/GlueX/";

    std::vector<TString> nominal_data_files = {
        nominal_data_dir + "PARA_0_allPeriods_data_signal.root",
        nominal_data_dir + "PARA_135_allPeriods_data_signal.root",
        nominal_data_dir + "PERP_45_allPeriods_data_signal.root",
        nominal_data_dir + "PERP_90_allPeriods_data_signal.root"
    };

    std::vector<TString> nominal_bkgd_files = {
        nominal_data_dir + "PARA_0_allPeriods_data_background.root",
        nominal_data_dir + "PARA_135_allPeriods_data_background.root",
        nominal_data_dir + "PERP_45_allPeriods_data_background.root",
        nominal_data_dir + "PERP_90_allPeriods_data_background.root"
    };

    std::vector<TString> systematic_data_files = {
        systematic_data_dir + "PARA_0_allPeriods_data_signal.root",
        systematic_data_dir + "PARA_135_allPeriods_data_signal.root",
        systematic_data_dir + "PERP_45_allPeriods_data_signal.root",
        systematic_data_dir + "PERP_90_allPeriods_data_signal.root"
    };

    std::vector<TString> systematic_bkgd_files = {
        systematic_data_dir + "PARA_0_allPeriods_data_background.root",
        systematic_data_dir + "PARA_135_allPeriods_data_background.root",
        systematic_data_dir + "PERP_45_allPeriods_data_background.root",
        systematic_data_dir + "PERP_90_allPeriods_data_background.root"
    };

    // -----------------------------------------------------------------------
    // Build nominal cuts map
    // The map is built as branch name -> pair of (low cut, high cut).
    // -----------------------------------------------------------------------
    std::map<TString, std::pair<float, float>> nominal_cuts = {
        {"unusedE",      {-0.1f,   0.1f}},  // nominal is unusedE < 0.1; lower bound keeps events with unusedE = 0
        {"unusedTracks", {-0.1f,   1.0f}},  // nominal is unusedTracks < 1, see above
        {"chi2",         {-0.1f,   5.0f}},  // nominal is chi2 < 5, see above
        {"PzCMrecoilPi", {-0.1f, 100.0f}},  // nominal is PzCMrecoilPi > -0.1; upper bound is arbitrary
        {"shQualityP2a", { 0.5f,   1.1f}},  // nominal is shQualityP2a > 0.5, see above
        {"shQualityP2b", { 0.5f,   1.1f}},  // nominal is shQualityP2b > 0.5, see above
        {"shQualityP3a", { 0.5f,   1.1f}},  // nominal is shQualityP2c > 0.5, see above
        {"shQualityP3b", { 0.5f,   1.1f}}   // nominal is shQualityP2d > 0.5, see above
    };

    // Parse the systematic argument: "variable:min:max"
    size_t first_colon  = systematic.find(':');
    size_t second_colon = systematic.find(':', first_colon + 1);

    if (first_colon == std::string::npos || second_colon == std::string::npos)
    {
        std::cerr << "Exiting: Skipping malformed optional cut: " << systematic << "\n";
        std::cerr << "  Expected format: variable:min:max\n";
        return 1;
    }

    std::string var_name = systematic.substr(0, first_colon);
    std::string min_val  = systematic.substr(first_colon + 1, second_colon - first_colon - 1);
    std::string max_val  = systematic.substr(second_colon + 1);

    // shQuality is a special case, as it needs to be applied across all 4 photons
    if (var_name == "shQuality")
    {
        for (const auto &suffix : {"P2a", "P2b", "P3a", "P3b"})
        {
            TString cut_name = TString::Format("shQuality%s", suffix);
            if (nominal_cuts.find(cut_name) == nominal_cuts.end())
            {
                std::cerr << "Exiting: Systematic variation branch not found in nominal cuts: " << cut_name << "\n";
                return 1;
            }
            nominal_cuts[cut_name] = {std::stof(min_val), std::stof(max_val)};

            std::cout << "Replaced nominal cut for " << cut_name
                      << " with new cut: " << nominal_cuts[cut_name].first
                      << " < " << cut_name << " < " << nominal_cuts[cut_name].second << "\n";
        }
    }
    else
    {
        if (nominal_cuts.find(var_name) == nominal_cuts.end())
        {
            std::cerr << "Exiting: Systematic variation branch not found in nominal cuts: " << var_name << "\n";
            return 1;
        }
        nominal_cuts[var_name] = {std::stof(min_val), std::stof(max_val)};

        std::cout << "Replaced nominal cut for " << var_name
                << " with new cut: " << nominal_cuts[var_name].first
                << " < " << var_name << " < " << nominal_cuts[var_name].second << "\n";
    }

    // Load all events into memory ONCE.
    std::cout << "Loading nominal events into memory...\n";
    EventCollection nominal_ec = load_events(nominal_data_files, nominal_bkgd_files, false);

    std::cout << "Loading systematic events into memory...\n";
    EventCollection systematic_ec = load_events(systematic_data_files, systematic_bkgd_files, true);

    
    // Bin definitions
    std::vector<float> mass_bin_edges;
    for (float edge = 1.0f; edge <= 1.8f + 1e-5f; edge += 0.02f)
        mass_bin_edges.push_back(edge);

    std::vector<float> t_bin_edges = {0.10f, 0.16f, 0.23f, 0.35f, 1.00f}; // GeV^2

    // -----------------------------------------------------------------------
    // Bin-by-bin comparison — purely in-memory, no file I/O inside the loop
    // -----------------------------------------------------------------------
    for (size_t t_idx = 0; t_idx < t_bin_edges.size() - 1; ++t_idx)
    {
        float t_min = t_bin_edges[t_idx];
        float t_max = t_bin_edges[t_idx + 1];

        std::cout << "Checking systematic variation for t bin: " << t_min << " < t < " << t_max << "\n";

        for (size_t m_idx = 0; m_idx < mass_bin_edges.size() - 1; ++m_idx)
        {
            float mass_min_bin = mass_bin_edges[m_idx];
            float mass_max_bin = mass_bin_edges[m_idx + 1];

            double n_nominal    = count_events(nominal_ec,    t_min, t_max, mass_min_bin, mass_max_bin);
            double n_systematic = count_events(systematic_ec, t_min, t_max, mass_min_bin, mass_max_bin, nominal_cuts);

            if (n_nominal == 0)
            {
                std::cerr << "Warning: No nominal events in bin t: [" << t_min << ", " << t_max
                          << "], mass: [" << mass_min_bin << ", " << mass_max_bin << "]. Skipping.\n";
                continue;
            }
            if (n_systematic == 0)
            {
                std::cerr << "Warning: No systematic events in bin t: [" << t_min << ", " << t_max
                          << "], mass: [" << mass_min_bin << ", " << mass_max_bin << "]. Skipping.\n";
                continue;
            }

            double relative_change = std::abs(n_systematic - n_nominal) / n_nominal;
            if (relative_change > 0.10)
            {
                std::cerr << "Error: Systematic variation exceeds 10% change in bin t: ["
                          << t_min << ", " << t_max << "], mass: ["
                          << mass_min_bin << ", " << mass_max_bin
                          << "]. Relative change: " << relative_change * 100.0 << "%\n";
                return 1;
            }
        }
    }

    std::cout << "Systematic variation check passed: No bins have more than 10% change in statistics compared to nominal.\n";
    return 0;
}

/**
 * @brief Load all events from signal and background files into an EventCollection.
 *
 * Signal events are stored with weight = +1.0. Background events are stored with
 * weight = -w (negated) so that @p count_events can simply accumulate weights
 * without special-casing.
 *
 * Only @p M4Pi and @p t are loaded unconditionally. When @p load_cut_vars is true,
 * the branches listed in @p CUT_BRANCH_NAMES are also read — this is needed for
 * the systematic collection where cuts must be re-applied in-memory.
 *
 * @param signal_files    Paths to signal ROOT files.
 * @param bkgd_files      Paths to background ROOT files.
 * @param load_cut_vars   If true, also read cut branches.
 */
EventCollection load_events(
    const std::vector<TString> &signal_files,
    const std::vector<TString> &bkgd_files,
    bool load_cut_vars)
{
    EventCollection ec;

    // Helper lambda: reads one TTree into ec.
    // fixed_weight == +1 for signal trees, -1 as sentinel for background trees
    // (background weight is then read from the "weight" branch and negated).
    auto read_tree = [&](TFile *file, const TString &fname, float fixed_weight)
    {
        TTree *tree = dynamic_cast<TTree *>(file->Get("ntFSGlueX_100_112"));
        if (!tree)
        {
            std::cerr << "Error: Could not find TTree in file: " << fname << "\n";
            return;
        }

        Double_t M4Pi_val = 0.0, t_val = 0.0, w_val = 1.0;
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("M4Pi", 1);  tree->SetBranchAddress("M4Pi", &M4Pi_val);
        tree->SetBranchStatus("t",    1);  tree->SetBranchAddress("t",    &t_val);

        // The "weight" branch only exists in background files.
        const bool is_background = (fixed_weight < 0);
        if (is_background && tree->GetBranch("weight"))
        {
            tree->SetBranchStatus("weight", 1);
            tree->SetBranchAddress("weight", &w_val);
        }

        // Cut branches are only needed for the systematic collection.
        std::map<std::string, Double_t> cut_vals;
        if (load_cut_vars)
        {
            for (const auto &name : CUT_BRANCH_NAMES)
            {
                if (tree->GetBranch(name.c_str()))
                {
                    tree->SetBranchStatus(name.c_str(), 1);
                    cut_vals[name] = 0.0;
                    tree->SetBranchAddress(name.c_str(), &cut_vals[name]);
                }
            }
        }

        const Long64_t n_entries = tree->GetEntries();
        ec.M4Pi.reserve(ec.M4Pi.size()     + static_cast<size_t>(n_entries));
        ec.t.reserve(ec.t.size()           + static_cast<size_t>(n_entries));
        ec.weight.reserve(ec.weight.size() + static_cast<size_t>(n_entries));

        for (Long64_t i = 0; i < n_entries; ++i)
        {
            tree->GetEntry(i);
            ec.M4Pi.push_back(M4Pi_val);
            ec.t.push_back(t_val);
            ec.weight.push_back(is_background ? -w_val : 1.0);

            if (load_cut_vars)
            {
                for (const auto &kv : cut_vals)
                    ec.cut_vars[kv.first].push_back(kv.second);
            }
        }
    };

    for (const auto &fname : signal_files)
    {
        TFile *file = TFile::Open(fname);
        if (!file || file->IsZombie()) { std::cerr << "Error: Could not open file: " << fname << "\n"; continue; }
        read_tree(file, fname, /*fixed_weight=*/+1.0f);
        file->Close();
        delete file;
    }

    for (const auto &fname : bkgd_files)
    {
        TFile *file = TFile::Open(fname);
        if (!file || file->IsZombie()) { std::cerr << "Error: Could not open file: " << fname << "\n"; continue; }
        read_tree(file, fname, /*fixed_weight=*/-1.0f); // -1 signals "use negated weight branch"
        file->Close();
        delete file;
    }

    return ec;
}

// ---------------------------------------------------------------------------

/**
 * @brief Sum (background-subtracted) event weights within a (t, mass) bin.
 *
 * Iterates over the in-memory @p events collection applying t, mass, and
 * optional @p cuts selections. Background events already carry a negative
 * weight from @p load_events, so a plain weighted sum gives the
 * background-subtracted integral.
 *
 * @param events    In-memory event collection from @p load_events.
 * @param t_min     Lower edge of the t bin.
 * @param t_max     Upper edge of the t bin.
 * @param mass_min  Lower edge of the mass bin.
 * @param mass_max  Upper edge of the mass bin.
 * @param cuts      Optional map of additional branch cuts (branch -> [lo, hi]).
 * @return          Background-subtracted integral over the bin.
 */
double count_events(
    const EventCollection &events,
    float t_min, float t_max,
    float mass_min, float mass_max,
    const std::map<TString, std::pair<float, float>> &cuts)
{
    double total = 0.0;
    const size_t n = events.M4Pi.size();

    for (size_t i = 0; i < n; ++i)
    {
        if (events.t[i]    <= t_min    || events.t[i]    >= t_max)    continue;
        if (events.M4Pi[i] <= mass_min || events.M4Pi[i] >= mass_max) continue;

        if (!cuts.empty())
        {
            bool pass = true;
            for (const auto &cut : cuts)
            {
                auto it = events.cut_vars.find(cut.first.Data());
                if (it == events.cut_vars.end()) continue; // branch was not loaded
                const Double_t val = it->second[i];
                if (val <= cut.second.first || val >= cut.second.second) { pass = false; break; }
            }
            if (!pass) continue;
        }

        total += events.weight[i];
    }

    return total;
}
