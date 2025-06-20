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
struct CellInfo { double width, height; int bits; };
struct PlacementRow { double startX, startY, siteWidth, siteHeight; int totalNumOfSites; };

// A temporary structure to hold information about a cell before it's legally placed
struct CandidateCell {
    string original_name_list; // Semicolon-separated list of original indices or original name
    string final_type;
    double ideal_x, ideal_y;
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
map<int, pair<double, string>> bits_minareaff_map;
map<string, int> clk_net_map;
// **FIXED**: Using a single, consistent name for the initial flip-flops vector
vector<Instance> initialFFs;
vector<Instance> initialGates;
vector<PlacementRow> placementRows;
vector<stPlacementMap> placementMaps;
pair<string, int> instNamePrefix;
int max_bit = 0;

//------------------------------
// Function Declarations
//------------------------------
void InputParsing(const char* fname);

// --- Placement Functions (Adapted from reference code) ---
void pmCreatePlacementMap();
void pmGateSetPlacementMap(const Instance& gate);
bool pmIsEmpty(int idx, int c0, int r0, int sitesC, int sitesR);
int SetNewInstPoint(Instance& new_instance);
bool pmFFSearchPlacementMap(Instance& ffInstance);


//------------------------------
// Main Program (Greedy Version with Legalization)
//------------------------------
int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: ./program <input_file> <output_file>" << endl;
        return 1;
    }

    InputParsing(argv[1]);

    if (!initialFFs.empty()) instNamePrefix.first = "C";
    
    set<int> possible_bits;
    for (const auto& pair : bits_minareaff_map) possible_bits.insert(pair.first);
    if (!possible_bits.empty()) max_bit = *possible_bits.rbegin();
    
    // --- Phase 1: Logical Greedy Clustering ---
    vector<CandidateCell> cells_to_legalize;
    vector<bool> processed_ffs(initialFFs.size(), false);

    map<int, vector<int>> clock_groups;
    for(size_t i = 0; i < initialFFs.size(); ++i) {
        int clk_id = clk_net_map.count(initialFFs[i].inst_name) ? clk_net_map.at(initialFFs[i].inst_name) : -1;
        clock_groups[clk_id].push_back(i);
    }

    for(auto const& [clk_id, indices] : clock_groups) {
        
        for(int seed_idx : indices) {
            if(processed_ffs[seed_idx]) continue;

            processed_ffs[seed_idx] = true;
            vector<int> current_cluster_indices = {seed_idx};
            int current_bits = cellLibrary.at(initialFFs[seed_idx].type_name).bits;
            double sum_x = initialFFs[seed_idx].x;
            double sum_y = initialFFs[seed_idx].y;

            while(true) {
                double min_dist_sq = numeric_limits<double>::max();
                int best_neighbor_idx = -1;
                
                double center_x = sum_x / current_cluster_indices.size();
                double center_y = sum_y / current_cluster_indices.size();

                for(int neighbor_idx : indices) {
                    if(processed_ffs[neighbor_idx]) continue;
                    
                    int neighbor_bits = cellLibrary.at(initialFFs[neighbor_idx].type_name).bits;
                    if (current_bits + neighbor_bits > max_bit) continue;

                    const auto& neighbor_inst = initialFFs[neighbor_idx];
                    double dist_sq = pow(center_x - neighbor_inst.x, 2) + pow(center_y - neighbor_inst.y, 2);
                    if(dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        best_neighbor_idx = neighbor_idx;
                    }
                }

                if(best_neighbor_idx != -1) {
                    current_cluster_indices.push_back(best_neighbor_idx);
                    current_bits += cellLibrary.at(initialFFs[best_neighbor_idx].type_name).bits;
                    sum_x += initialFFs[best_neighbor_idx].x;
                    sum_y += initialFFs[best_neighbor_idx].y;
                    processed_ffs[best_neighbor_idx] = true;
                } else {
                    break;
                }
            }
            
            CandidateCell candidate;
            candidate.ideal_x = sum_x / current_cluster_indices.size();
            candidate.ideal_y = sum_y / current_cluster_indices.size();

            if (current_cluster_indices.size() > 1) {
                 auto it_map = bits_minareaff_map.lower_bound(current_bits);
                 if (it_map != bits_minareaff_map.end() && it_map->first <= max_bit) {
                    candidate.final_type = it_map->second.second;
                    for (int idx : current_cluster_indices) {
                        candidate.original_name_list += to_string(idx) + ";";
                    }
                 } else {
                    for(int idx : current_cluster_indices) {
                        const auto& inst = initialFFs[idx];
                        cells_to_legalize.push_back({inst.inst_name, inst.type_name, inst.x, inst.y});
                    }
                    continue;
                 }
            } else {
                candidate.final_type = initialFFs[seed_idx].type_name;
                candidate.original_name_list = initialFFs[seed_idx].inst_name;
            }
            cells_to_legalize.push_back(candidate);
        }
    }
    
    // --- Phase 2: Physical Placement Legalization ---
    vector<Instance> final_instances;
    vector<string> pin_mappings;
    int new_inst_counter = initialFFs.size() + initialGates.size() + 1;
    
    pmCreatePlacementMap();

    for (const auto& candidate : cells_to_legalize) {
        Instance inst_to_place;
        inst_to_place.type_name = candidate.final_type;
        inst_to_place.x = candidate.ideal_x;
        inst_to_place.y = candidate.ideal_y;
        
        if (pmFFSearchPlacementMap(inst_to_place)) {
            string final_name = instNamePrefix.first + to_string(new_inst_counter++);
            inst_to_place.inst_name = final_name;
            final_instances.push_back(inst_to_place);

            bool is_newly_banked = (candidate.original_name_list.find(';') != string::npos);
            
            if (is_newly_banked) {
                stringstream ss(candidate.original_name_list);
                string segment;
                int pin_offset = 0;
                while(getline(ss, segment, ';')) {
                    int old_idx = stoi(segment);
                    const auto& old_inst = initialFFs[old_idx];
                    int old_bits = cellLibrary.at(old_inst.type_name).bits;
                    int new_total_bits = cellLibrary.at(inst_to_place.type_name).bits;
                    for (int k = 0; k < old_bits; ++k) {
                        string old_suf = (old_bits > 1) ? to_string(k) : "";
                        string new_suf = (new_total_bits > 1) ? to_string(pin_offset + k) : "";
                        pin_mappings.push_back(old_inst.inst_name + "/D" + old_suf + " map " + final_name + "/D" + new_suf);
                        pin_mappings.push_back(old_inst.inst_name + "/Q" + old_suf + " map " + final_name + "/Q" + new_suf);
                    }
                    pin_mappings.push_back(old_inst.inst_name + "/CLK map " + final_name + "/CLK");
                    pin_offset += old_bits;
                }
            } else {
                const string& old_name = candidate.original_name_list;
                int bits = cellLibrary.at(candidate.final_type).bits;
                for (int k = 0; k < bits; ++k) {
                    string pin_suffix = (bits > 1) ? to_string(k) : "";
                    pin_mappings.push_back(old_name + "/D" + pin_suffix + " map " + final_name + "/D" + pin_suffix);
                    pin_mappings.push_back(old_name + "/Q" + pin_suffix + " map " + final_name + "/Q" + pin_suffix);
                }
                pin_mappings.push_back(old_name + "/CLK map " + final_name + "/CLK");
            }
        } else {
             cerr << "CRITICAL WARNING: Could not legalize placement for a cell of type " << candidate.final_type << endl;
        }
    }
    
    // --- Phase 3: Output Generation ---
    ofstream fout(argv[2]);
    fout << "CellInst " << final_instances.size() << endl;
    for(const auto& inst : final_instances) {
        fout << "Inst " << inst.inst_name << " " << inst.type_name << " " << inst.x << " " << inst.y << endl;
    }
    for(const auto& mapping : pin_mappings) {
        fout << mapping << endl;
    }
    fout.close();

    return 0;
}

// --- Function Implementations ---

void InputParsing(const char *fname) {
    string line; ifstream file(fname);
    if (!file.is_open()) { cerr << "Error: Could not open file '" << fname << "'." << endl; exit(1); }
    int clk_group_idx = 0;
    map<string, double> ff_areas;
    while (getline(file, line)) {
        istringstream iss(line); string key; iss >> key;
        if (key == "FlipFlop" || key == "Gate") {
            CellInfo info; string name; int p_count;
            if (key == "FlipFlop") {
                iss >> info.bits >> name >> info.width >> info.height >> p_count;
                double area = info.width * info.height;
                if (bits_minareaff_map.find(info.bits) == bits_minareaff_map.end() || area < bits_minareaff_map[info.bits].first) {
                    bits_minareaff_map[info.bits] = {area, name};
                }
            } else {
                iss >> name >> info.width >> info.height >> p_count; info.bits = 0;
            }
            cellLibrary[name] = info;
            for(int i = 0; i < p_count; ++i) { if(!getline(file, line)) break; }
        } else if (key == "NumInstances") {
            int instanceCount; iss >> instanceCount;
            for (int i = 0; i < instanceCount; i++) {
                getline(file, line);
                istringstream iss_inst(line);
                Instance instance; string tag;
                iss_inst >> tag >> instance.inst_name >> instance.type_name >> instance.x >> instance.y;
                if (cellLibrary.count(instance.type_name) && cellLibrary.at(instance.type_name).bits > 0) {
                    initialFFs.push_back(instance);
                    if(instNamePrefix.first == "") instNamePrefix = {instance.inst_name.substr(0, instance.inst_name.find_first_of("0123456789")), 0};
                } else {
                    initialGates.push_back(instance);
                }
            }
        } else if (key == "Net") {
            string name; int numPins; iss >> name >> numPins;
            for (int i = 0; i < numPins; ++i) {
                if(!getline(file, line)) break;
                istringstream pinIss(line); string pin_type, connect_to; pinIss >> pin_type >> connect_to;
                if (connect_to.size() >= 4 && connect_to.substr(connect_to.size() - 4) == "/CLK") {
                    string ffname = connect_to.substr(0, connect_to.find('/'));
                    clk_net_map[ffname] = clk_group_idx;
                }
            }
            clk_group_idx++;
        } else if (key == "PlacementRows") {
            PlacementRow pr;
            iss >> pr.startX >> pr.startY >> pr.siteWidth >> pr.siteHeight >> pr.totalNumOfSites;
            placementRows.push_back(pr);
        }
    }
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

void pmCreatePlacementMap(void) {
    if (placementRows.empty()) return;
    placementMaps.clear();
    sort(placementRows.begin(), placementRows.end(), [](const PlacementRow& a, const PlacementRow& b){
        if (a.startY != b.startY) return a.startY < b.startY; return a.startX < b.startX;
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
    for (const auto& instance : initialGates) pmGateSetPlacementMap(instance);
}

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
    if(idx < 0 || idx >= (int)placementMaps.size()) return false;
    stPlacementMap &pm = placementMaps[idx];
    if ((c0 + sitesC) > pm.bmCols || (r0 + sitesR) > pm.bmRows || c0 < 0 || r0 < 0) return false;
    for (int r=r0; r < r0 + sitesR; r++) for (int c=c0; c < c0 + sitesC; c++) if (pm.byteMap[r][c]) return false;
    return true;
}

int SetNewInstPoint(Instance& new_instance){
    vector<pair<double, double>> possible_points;
    for(const auto& pm : placementMaps){
        double new_x, new_y;
        if (new_instance.x <= pm.startX) new_x = pm.startX;
        else if (new_instance.x >= pm.endX) new_x = pm.endX - pm.siteWidth;
        else new_x = pm.startX + round((new_instance.x - pm.startX) / pm.siteWidth) * pm.siteWidth;
        if (new_instance.y <= pm.startY) new_y = pm.startY;
        else if (new_instance.y >= pm.endY) new_y = pm.endY - pm.siteHeight;
        else new_y = pm.startY + round((new_instance.y - pm.startY) / pm.siteHeight) * pm.siteHeight;
        possible_points.push_back({new_x, new_y});
    }
    int min_idx = 0; double min_dist_sq = -1;
    for(size_t i = 0; i < possible_points.size(); ++i) {
        double dist_sq = pow(possible_points[i].first - new_instance.x, 2) + pow(possible_points[i].second - new_instance.y, 2);
        if(min_dist_sq < 0 || dist_sq < min_dist_sq) { min_dist_sq = dist_sq; min_idx = i; }
    }
    new_instance.x = possible_points[min_idx].first;
    new_instance.y = possible_points[min_idx].second;
    return min_idx;
}

bool pmFFSearchPlacementMap(Instance& ffInstance){
    int idx = SetNewInstPoint(ffInstance);
    if(idx < 0 || idx >= (int)placementMaps.size()) return false;

    stPlacementMap &pm = placementMaps[idx];
    if (cellLibrary.find(ffInstance.type_name) == cellLibrary.end()) return false;
    const auto& ffValue = cellLibrary.at(ffInstance.type_name);
    int sitesC = ceil(ffValue.width / pm.siteWidth);
    int sitesR_to_occupy = ceil(ffValue.height / pm.siteHeight);
    
    // Brute-force scan of the entire grid map. This is slow but robust.
    for (int rt = 0; rt <= pm.bmRows - sitesR_to_occupy; ++rt) {
        // Check for vertical contiguity before attempting to search the row
        bool contiguous = true;
        if (sitesR_to_occupy > 1) {
            for(int i = 1; i < sitesR_to_occupy; ++i) {
                // This check is flawed if rows are not uniform height, but works for the sample case
                if (abs((pm.startY + (rt+i-1)*pm.siteHeight + pm.siteHeight) - (pm.startY + (rt+i)*pm.siteHeight)) > 1e-6) {
                   contiguous = false; break;
                }
            }
        }
        if(!contiguous) continue;

        for (int ct = 0; ct <= pm.bmCols - sitesC; ++ct) {
            if (pmIsEmpty(idx, ct, rt, sitesC, sitesR_to_occupy)) {
                ffInstance.x = pm.startX + pm.siteWidth * ct;
                ffInstance.y = pm.startY + pm.siteHeight * rt;
                uint8_t bits_char = (uint8_t)ffValue.bits + '0';
                for (int r_occ=rt; r_occ < rt + sitesR_to_occupy && r_occ < pm.bmRows; r_occ++)
                    for (int c_occ=ct; c_occ < ct + sitesC && c_occ < pm.bmCols; c_occ++)
                        pm.byteMap[r_occ][c_occ] = bits_char;
                return true;
            }
        }
    }
    
    // Fallback search across all other maps
    for(size_t other_idx = 0; other_idx < placementMaps.size(); ++other_idx) {
        if((int)other_idx == idx) continue;
        stPlacementMap &other_pm = placementMaps[other_idx];
         for (int rt = 0; rt <= other_pm.bmRows - sitesR_to_occupy; ++rt) {
            for (int ct = 0; ct <= other_pm.bmCols - sitesC; ++ct) {
                if (pmIsEmpty(other_idx, ct, rt, sitesC, sitesR_to_occupy)) {
                    ffInstance.x = other_pm.startX + other_pm.siteWidth * ct;
                    ffInstance.y = other_pm.startY + other_pm.siteHeight * rt;
                    uint8_t bits_char = (uint8_t)ffValue.bits + '0';
                    for (int r_occ=rt; r_occ < rt + sitesR_to_occupy && r_occ < other_pm.bmRows; r_occ++)
                        for (int c_occ=ct; c_occ < ct + sitesC && c_occ < other_pm.bmCols; c_occ++)
                            other_pm.byteMap[r_occ][c_occ] = bits_char;
                    return true;
                }
            }
        }
    }
    
    cout << "!!NG Cannot Place FF: " << ffInstance.inst_name << endl;
    return false;
}

