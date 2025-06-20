#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>
#include <cassert>

using namespace std;

//------------------------------
// Data Structures
//------------------------------
struct Instance { string inst_name, type_name; double x, y; };
struct CellInfo { double width, height; int bits; double power = 0.0; double qpin_delay = 0.0; };
struct PlacementRow { double startX, startY, siteWidth, siteHeight; int totalNumOfSites; };
struct Point { double x, y; int cluster_id = -1; };

// A temporary structure to hold information about a cell before it's legally placed
struct CandidateCell {
    string original_name; // For pass-through cells
    string final_type;
    double ideal_x, ideal_y;
    vector<int> original_indices; // List of initialFFs indices that form this cell
};

// Placement structures adapted from the reference code
struct stPlacementMap {
    vector<vector<uint8_t>> byteMap;
    int bmRows = 0, bmCols = 0;
    double startX = 0, startY = 0, siteWidth = 0, siteHeight = 0, endX = 0, endY = 0;
};

//------------------------------
// Global Variables
//------------------------------
map<string, CellInfo> cellLibrary;
map<int, string> bits_minAreaFFType_map;
map<string, int> clk_net_map;
vector<Instance> initialFFs;
vector<Instance> initialGates;
vector<PlacementRow> placementRows;
vector<stPlacementMap> placementMaps;
int max_bit = 0;
double ALPHA = 1.0, BETA = 1.0, GAMMA = 1.0;

//------------------------------
// Function Declarations
//------------------------------
void InputParsing(const char* fname);
vector<int> k_means(vector<Point>& points, int k);
void pmCreatePlacementMap();
// **FIXED**: Changed function signature to accept CandidateCell
bool pmFFSearchPlacementMap(CandidateCell& candidate);


//------------------------------
// Main Program
//------------------------------
int main(int argc, char* argv[]) {
    if (argc < 3) { return 1; }

    InputParsing(argv[1]);
    if (!bits_minAreaFFType_map.empty()) {
        max_bit = bits_minAreaFFType_map.rbegin()->first;
    }

    vector<CandidateCell> cells_to_legalize;
    vector<bool> processed(initialFFs.size(), false);

    map<int, vector<int>> clock_groups;
    for (size_t i = 0; i < initialFFs.size(); ++i) {
        clk_net_map.count(initialFFs[i].inst_name) ? clock_groups[clk_net_map[initialFFs[i].inst_name]].push_back(i) : clock_groups[-1].push_back(i);
    }
    
    for (auto const& [clk_id, indices] : clock_groups) {
        if (indices.size() <= 1) continue;

        vector<Point> points_for_kmeans;
        int total_bits_in_group = 0;
        for (int idx : indices) {
            points_for_kmeans.push_back({initialFFs[idx].x, initialFFs[idx].y});
            total_bits_in_group += cellLibrary.at(initialFFs[idx].type_name).bits;
        }
        
        int k = (max_bit > 0) ? ceil((double)total_bits_in_group / max_bit) : indices.size();
        k = max(1, k);

        map<int, vector<int>> clusters;
        vector<int> labels = k_means(points_for_kmeans, k);
        for (size_t i = 0; i < labels.size(); ++i) clusters[labels[i]].push_back(indices[i]);
        
        for (auto const& [cluster_id, cluster_indices] : clusters) {
            if (cluster_indices.size() <= 1) continue;

            int total_bits = 0; double sum_x = 0, sum_y = 0;
            double original_cost_power = 0, original_cost_area = 0, original_cost_timing = 0;

            for (int idx : cluster_indices) {
                const auto& old_info = cellLibrary.at(initialFFs[idx].type_name);
                total_bits += old_info.bits; sum_x += initialFFs[idx].x; sum_y += initialFFs[idx].y;
                original_cost_power += old_info.power;
                original_cost_area += old_info.width * old_info.height;
                original_cost_timing += old_info.qpin_delay;
            }
            
            auto it = bits_minAreaFFType_map.lower_bound(total_bits);
            bool banking_is_beneficial = false;
            if (it != bits_minAreaFFType_map.end() && it->first <= max_bit) {
                const auto& new_ff_type = it->second;
                const auto& new_ff_info = cellLibrary.at(new_ff_type);

                double new_cost_power = new_ff_info.power;
                double new_cost_area = new_ff_info.width * new_ff_info.height;
                double new_cost_timing = new_ff_info.qpin_delay;
                double cost_change = (ALPHA * (new_cost_timing - original_cost_timing)) + (BETA * (new_cost_power - original_cost_power)) + (GAMMA * (new_cost_area - original_cost_area));
                
                if (cost_change < 0) banking_is_beneficial = true;
            }

            if(banking_is_beneficial) {
                cells_to_legalize.push_back({"", it->second, sum_x / cluster_indices.size(), sum_y / cluster_indices.size(), cluster_indices});
                for (int idx : cluster_indices) processed[idx] = true;
            }
        }
    }

    for (size_t i = 0; i < initialFFs.size(); ++i) {
        if (!processed[i]) {
            const auto& old_inst = initialFFs[i];
            cells_to_legalize.push_back({old_inst.inst_name, old_inst.type_name, old_inst.x, old_inst.y, {(int)i}});
        }
    }

    vector<Instance> final_instances;
    vector<string> pin_mappings;
    
    pmCreatePlacementMap();
    int new_inst_counter = initialFFs.size() + initialGates.size() + 1;

    for (auto& candidate : cells_to_legalize) {
        // **FIXED**: Call to pmFFSearchPlacementMap now uses a CandidateCell
        if (pmFFSearchPlacementMap(candidate)) {
            string final_name = "C" + to_string(new_inst_counter++);
            
            bool is_pass_through = (candidate.original_indices.size() == 1 && initialFFs[candidate.original_indices[0]].type_name == candidate.final_type);
            
            if (!is_pass_through) { // Banked cell
                int pin_offset = 0;
                for (int old_idx : candidate.original_indices) {
                    const auto& old_inst = initialFFs[old_idx];
                    int old_bits = cellLibrary.at(old_inst.type_name).bits;
                    int new_total_bits = cellLibrary.at(candidate.final_type).bits;
                    for (int k = 0; k < old_bits; ++k) {
                        string old_suf = (old_bits > 1) ? to_string(k) : "";
                        string new_suf = (new_total_bits > 1) ? to_string(pin_offset + k) : "";
                        pin_mappings.push_back(old_inst.inst_name + "/D" + old_suf + " map " + final_name + "/D" + new_suf);
                        pin_mappings.push_back(old_inst.inst_name + "/Q" + old_suf + " map " + final_name + "/Q" + new_suf);
                    }
                    pin_mappings.push_back(old_inst.inst_name + "/CLK map " + final_name + "/CLK");
                    pin_offset += old_bits;
                }
            } else { // Pass-through cell
                const auto& old_inst = initialFFs[candidate.original_indices[0]];
                int bits = cellLibrary.at(old_inst.type_name).bits;
                for (int k = 0; k < bits; ++k) {
                    string pin_suffix = (bits > 1) ? to_string(k) : "";
                    pin_mappings.push_back(old_inst.inst_name + "/D" + pin_suffix + " map " + final_name + "/D" + pin_suffix);
                    pin_mappings.push_back(old_inst.inst_name + "/Q" + pin_suffix + " map " + final_name + "/Q" + pin_suffix);
                }
                pin_mappings.push_back(old_inst.inst_name + "/CLK map " + final_name + "/CLK");
            }
            final_instances.push_back({final_name, candidate.final_type, candidate.ideal_x, candidate.ideal_y});
        }
    }
    
    ofstream fout(argv[2]);
    fout << "CellInst " << final_instances.size() << endl;
    for (const auto& inst : final_instances) fout << "Inst " << inst.inst_name << " " << inst.type_name << " " << inst.x << " " << inst.y << endl;
    for (const auto& mapping : pin_mappings) fout << mapping << endl;
    fout.close();

    return 0;
}


// --- Function Implementations ---

void InputParsing(const char* fname) {
    ifstream file(fname); string line, key; int clk_group_idx = 0; map<string, double> ff_areas;
    while (getline(file, line)) {
        istringstream iss(line); iss >> key;
        if (key == "Alpha") { iss >> ALPHA; } else if (key == "Beta") { iss >> BETA; } else if (key == "Gamma") { iss >> GAMMA; }
        else if (key == "FlipFlop" || key == "Gate") {
            CellInfo info; string name; int p_count;
            if (key == "FlipFlop") {
                iss >> info.bits >> name >> info.width >> info.height >> p_count;
                double area = info.width * info.height;
                if (!bits_minAreaFFType_map.count(info.bits) || area < ff_areas[bits_minAreaFFType_map[info.bits]]) {
                    bits_minAreaFFType_map[info.bits] = name; ff_areas[name] = area;
                }
            } else { iss >> name >> info.width >> info.height >> p_count; info.bits = 0; }
            cellLibrary[name] = info;
            for (int i = 0; i < p_count; ++i) getline(file, line);
        } else if (key == "GatePower") {
            string name; double power; iss >> name >> power; if(cellLibrary.count(name)) cellLibrary.at(name).power = power;
        } else if (key == "QpinDelay") {
            string name; double delay; iss >> name >> delay; if(cellLibrary.count(name)) cellLibrary.at(name).qpin_delay = delay;
        } else if (key == "Inst") {
            Instance inst; iss >> inst.inst_name >> inst.type_name >> inst.x >> inst.y;
            if (cellLibrary.count(inst.type_name)) {
                if (cellLibrary.at(inst.type_name).bits > 0) initialFFs.push_back(inst); else initialGates.push_back(inst);
            }
        } else if (key == "Net") {
            string netName; int numPins; iss >> netName >> numPins; bool is_clk_net = false; vector<string> connected_ffs;
            for (int i = 0; i < numPins; ++i) {
                getline(file, line);
                string pin_line_str(line); size_t first_char = pin_line_str.find_first_not_of(" \t\r\n");
                if (string::npos != first_char) pin_line_str = pin_line_str.substr(first_char);
                istringstream pinIss(pin_line_str); string pin_type, connect_to; pinIss >> pin_type >> connect_to;
                size_t slash_pos = connect_to.find('/');
                if (slash_pos != string::npos && connect_to.substr(slash_pos + 1) == "CLK") {
                    is_clk_net = true; connected_ffs.push_back(connect_to.substr(0, slash_pos));
                }
            }
            if (is_clk_net) { for (const auto& ff_name : connected_ffs) clk_net_map[ff_name] = clk_group_idx; clk_group_idx++; }
        } else if (key == "PlacementRows") {
            PlacementRow pr; iss >> pr.startX >> pr.startY >> pr.siteWidth >> pr.siteHeight >> pr.totalNumOfSites;
            placementRows.push_back(pr);
        }
    }
}

vector<int> k_means(vector<Point>& points, int k) {
    if (k <= 0 || points.empty()) return vector<int>(points.size(), -1);
    vector<Point> centroids(k); vector<int> initial_indices(points.size());
    iota(initial_indices.begin(), initial_indices.end(), 0);
    random_device rd; mt19937 g(rd());
    shuffle(initial_indices.begin(), initial_indices.end(), g);
    for (int i = 0; i < k; ++i) centroids[i] = {points[initial_indices[i % points.size()]].x, points[initial_indices[i % points.size()]].y};
    bool changed = true;
    for (int iter = 0; iter < 100 && changed; ++iter) {
        changed = false;
        for (auto& p : points) {
            double min_dist = numeric_limits<double>::max(); int best_cluster = p.cluster_id;
            for (int i = 0; i < k; ++i) {
                double dist = pow(p.x - centroids[i].x, 2) + pow(p.y - centroids[i].y, 2);
                if (dist < min_dist) { min_dist = dist; best_cluster = i; }
            }
            if (p.cluster_id != best_cluster) { p.cluster_id = best_cluster; changed = true; }
        }
        if (!changed) break;
        vector<Point> sums(k, {0, 0, -1}); vector<int> counts(k, 0);
        for (const auto& p : points) if (p.cluster_id != -1) { sums[p.cluster_id].x += p.x; sums[p.cluster_id].y += p.y; counts[p.cluster_id]++; }
        for (int i = 0; i < k; ++i) if (counts[i] > 0) centroids[i] = {sums[i].x / counts[i], sums[i].y / counts[i]};
    }
    vector<int> labels; for (const auto& p : points) labels.push_back(p.cluster_id);
    return labels;
}

// --- Placement Functions (Adapted from Reference Code) ---
void pmUpdatePlacementMap(int startR, int rowCnt);
void pmGateSetPlacementMap(const Instance& gate);
int SetNewInstPoint(CandidateCell& candidate);

void pmCreatePlacementMap(void) {
    if (placementRows.empty()) return;
    sort(placementRows.begin(), placementRows.end(), [](const PlacementRow& a, const PlacementRow& b){
        if (a.startY != b.startY) return a.startY < b.startY;
        return a.startX < b.startX;
    });
    int totalRows = placementRows.size(); int startR = 0;
    for (int r = 1; r < totalRows; r++) {
        if ((abs(placementRows[r-1].startX - placementRows[r].startX) > 1e-6)
          || (placementRows[r-1].totalNumOfSites != placementRows[r].totalNumOfSites)
          || (abs((placementRows[r-1].startY + placementRows[r-1].siteHeight) - placementRows[r].startY) > 1e-6 )) {
            pmUpdatePlacementMap(startR, r - startR); startR = r;
        }
    }
    pmUpdatePlacementMap(startR, totalRows - startR);
    // **FIXED**: Call to pmGateSetPlacementMap now exists
    for (const auto& instance : initialGates) pmGateSetPlacementMap(instance);
}

void pmUpdatePlacementMap(int startR, int rowCnt) {
    stPlacementMap pm;
    pm.bmCols = placementRows[startR].totalNumOfSites; pm.bmRows = rowCnt;
    pm.startX = placementRows[startR].startX; pm.startY = placementRows[startR].startY;
    pm.siteWidth = placementRows[startR].siteWidth; pm.siteHeight = placementRows[startR].siteHeight;
    pm.endX = pm.startX + pm.siteWidth * pm.bmCols; pm.endY = pm.startY + pm.siteHeight * pm.bmRows;
    pm.byteMap.resize(pm.bmRows, vector<uint8_t>(pm.bmCols, 0));
    placementMaps.push_back(pm);
}

// **FIXED**: Added the missing function definition
void pmGateSetPlacementMap(const Instance& gate){
    if (cellLibrary.find(gate.type_name) == cellLibrary.end()) return;
    const auto& gateValue = cellLibrary.at(gate.type_name);
    double gateEndX = gate.x + gateValue.width; double gateEndY = gate.y + gateValue.height;
    for(size_t idx=0; idx<placementMaps.size(); idx++){
        stPlacementMap &pm = placementMaps[idx];
        if ((gateEndX < pm.startX || pm.endX < gate.x) || (gateEndY < pm.startY || pm.endY < gate.y)) continue;
        int c0 = (gate.x > pm.startX) ? floor((gate.x - pm.startX) / pm.siteWidth) : 0;
        int r0 = (gate.y > pm.startY) ? floor((gate.y - pm.startY) / pm.siteHeight) : 0;
        int c1 = ceil((gateEndX - pm.startX) / pm.siteWidth) - 1;
        if (c1 >= pm.bmCols) c1 = pm.bmCols - 1;
        int r1 = ceil((gateEndY - pm.startY) / pm.siteHeight) - 1;
        if (r1 >= pm.bmRows) r1 = pm.bmRows - 1;
        for (int r=r0; r<=r1; r++) for (int c=c0; c<=c1; c++) if(r < pm.bmRows && c < pm.bmCols) pm.byteMap[r][c] = 'G';
        return;
    }
}

bool pmIsEmpty(int idx, int c0, int r0, int sitesC, int sitesR){
    stPlacementMap &pm = placementMaps[idx];
    if ((c0 + sitesC) > pm.bmCols || (r0 + sitesR) > pm.bmRows || c0 < 0 || r0 < 0) return false;
    for (int r=r0; r < r0 + sitesR; r++) for (int c=c0; c < c0 + sitesC; c++) if (pm.byteMap[r][c]) return false;
    return true;
}

int SetNewInstPoint(CandidateCell& candidate){
    vector<pair<double, double>> possible_points;
    for(const auto& pm : placementMaps){
        double new_x, new_y;
        if (candidate.ideal_x <= pm.startX) new_x = pm.startX;
        else if (candidate.ideal_x >= pm.endX) new_x = pm.endX - pm.siteWidth;
        else new_x = pm.startX + round((candidate.ideal_x - pm.startX) / pm.siteWidth) * pm.siteWidth;
        if (candidate.ideal_y <= pm.startY) new_y = pm.startY;
        else if (candidate.ideal_y >= pm.endY) new_y = pm.endY - pm.siteHeight;
        else new_y = pm.startY + round((candidate.ideal_y - pm.startY) / pm.siteHeight) * pm.siteHeight;
        possible_points.push_back({new_x, new_y});
    }
    int min_idx = 0; double min_dist_sq = -1;
    for(size_t i = 0; i < possible_points.size(); ++i) {
        double dist_sq = pow(possible_points[i].first - candidate.ideal_x, 2) + pow(possible_points[i].second - candidate.ideal_y, 2);
        if(min_dist_sq < 0 || dist_sq < min_dist_sq) { min_dist_sq = dist_sq; min_idx = i; }
    }
    candidate.ideal_x = possible_points[min_idx].first;
    candidate.ideal_y = possible_points[min_idx].second;
    return min_idx;
}

// **FIXED**: Signature changed to take CandidateCell&
bool pmFFSearchPlacementMap(CandidateCell& candidate){
    int idx = SetNewInstPoint(candidate);
    stPlacementMap &pm = placementMaps[idx];
    if (cellLibrary.find(candidate.final_type) == cellLibrary.end()) return false;
    const auto& ffValue = cellLibrary.at(candidate.final_type);
    int sitesC = ceil(ffValue.width / pm.siteWidth);
    int sitesR = ceil(ffValue.height / pm.siteHeight);
    int c0 = floor((candidate.ideal_x - pm.startX) / pm.siteWidth);
    int r0 = floor((candidate.ideal_y - pm.startY) / pm.siteHeight);
    int ct = c0, rt = r0;
    
    bool bEmpty = pmIsEmpty(idx, ct, rt, sitesC, sitesR);
    if (!bEmpty) {
        int dmax = max({c0, pm.bmCols-c0-1, r0, pm.bmRows-r0-1});
        for (int d=1; !bEmpty && d<=dmax ; d++) {
            for(int side = 0; side < 4 && !bEmpty; ++side) {
                for(int step = 0; step < d * 2; ++step) {
                    if (pmIsEmpty(idx, ct, rt, sitesC, sitesR)) { bEmpty = true; break; }
                    switch(side) {
                        case 0: ct++; break; case 1: rt--; break;
                        case 2: ct--; break; case 3: rt++; break;
                    }
                }
            }
        }
    }
    if (bEmpty) {
        candidate.ideal_x = pm.startX + pm.siteWidth * ct;
        candidate.ideal_y = pm.startY + pm.siteHeight * rt;
        uint8_t bits_char = (uint8_t)ffValue.bits + '0';
        for (int r=rt; r < rt + sitesR && r < pm.bmRows; r++) for (int c=ct; c < ct + sitesC && c < pm.bmCols; c++) pm.byteMap[r][c] = bits_char;
        return true;
    }
    return false;
}

