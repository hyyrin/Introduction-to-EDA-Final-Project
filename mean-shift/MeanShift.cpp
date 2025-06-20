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
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cmath>
#include "utils.h"

using namespace std;

#define _CJDBG              1   // separate instances into ffInstances and gateInstances
#define _DBG_PlacementMAP   0

//------------------------------
// Data Structures
//------------------------------
struct diesize { double x_left, y_bottom, x_right, y_up; };
struct Pin { string name; double x, y; };
struct FlipFlop { int bits; string name; double width, height; int pinCount; vector<Pin> pins; };
struct Gate { string name; double width, height; int pinCount; vector<Pin> pins; };
struct Instance { string inst_name, type_name; double x, y; };
struct Net { string name; int numPins; vector<string> pins; };
struct PlacementRow { double startX, startY, siteWidth, siteHeight; int totalNumOfSites; };
struct Qpindelay { string flipflopname; double value; };
struct TimingSlack { string instanceCellName, pinName; double slack; };
struct GatePower { string libCellName; double powerConsumption; };
struct stFFValue { int bits; double ffWidth, ffHeight; int pinCount; int sitesC, sitesR; };
struct stGateValue { double gateWidth, gateHeight; int pinCount; };
struct stPlacementMap { vector<vector<uint8_t>> byteMap; int bmRows, bmCols; double startX, startY, siteWidth, siteHeight, endX, endY; };

//------------------------------
// Global Variables
//------------------------------
map<string, int> type_bits_map;
map<int, pair<double, string>> bits_minareaff_map;
map<string, stFFValue> kmFFLib;
map<string, stGateValue> kmGateLib;
set<int> possible_bits;
int max_bit = 0;
string ffTag = "";
unordered_map<string, double> weights;
vector<Instance> ffInstances;
vector<Instance> gateInstances;
map<string, int> clk_net_map;
vector<PlacementRow> placementRows;
vector<stPlacementMap> placementMaps;
pair<string, int> instName;


//------------------------------
// Function Declarations
//------------------------------
void InputParsing(char* fname);
pair<string, int> splitString(const string& str);
void pmCreatePlacementMap(void);
bool pmFFSearchPlacementMap(Instance& ffInstance);

bool CompareInstancesInCluster(int idx_a, int idx_b) {
    const auto& inst_a = ffInstances[idx_a];
    const auto& inst_b = ffInstances[idx_b];
    int clk_a = clk_net_map.count(inst_a.inst_name) ? clk_net_map[inst_a.inst_name] : -1;
    int clk_b = clk_net_map.count(inst_b.inst_name) ? clk_net_map[inst_b.inst_name] : -1;
    if (clk_a != clk_b) return clk_a < clk_b;
    int bits_a = type_bits_map.count(inst_a.type_name) ? type_bits_map[inst_a.type_name] : 0;
    int bits_b = type_bits_map.count(inst_b.type_name) ? type_bits_map[inst_b.type_name] : 0;
    return bits_a > bits_b;
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    // --- 1. Parsing ---
    cout << "Parsing input file: " << argv[1] << endl;
    InputParsing(argv[1]);
    if (!possible_bits.empty()) {
        max_bit = *max_element(possible_bits.begin(), possible_bits.end());
    }
    cout << "Max flip-flop bit size: " << max_bit << endl;

    // --- 2. Clustering ---
    cout << "\nPreparing points for clustering..." << endl;
    vector<MeanShift::Point> points;
    points.reserve(ffInstances.size());
    for(const auto& instance : ffInstances) points.push_back({instance.x, instance.y});
    cout << "Starting Mean Shift clustering..." << endl;
    MeanShift ms;
    double kernel_bandwidth = 1000.0;
    vector<Cluster> clusters = ms.cluster(points, kernel_bandwidth);
    cout << "Found " << clusters.size() << " clusters." << endl;

    // --- 3. Banking and Transformation ---
    cout << "\nProcessing clusters for banking..." << endl;
    vector<Instance> final_instances;
    vector<string> output_mapping_lines;
    vector<bool> processed_ffs(ffInstances.size(), false);
    int new_idx_counter = ffInstances.size() + gateInstances.size() + 1;

    for (auto& cluster : clusters) {
        sort(cluster.original_reg_idx.begin(), cluster.original_reg_idx.end(), CompareInstancesInCluster);

        for (size_t i = 0; i < cluster.original_reg_idx.size(); ) {
            int start_idx = i;
            int current_clk_net = clk_net_map.count(ffInstances[cluster.original_reg_idx[i]].inst_name) ? clk_net_map[ffInstances[cluster.original_reg_idx[i]].inst_name] : -1;
            int total_bits = 0;
            vector<double> x_coords, y_coords;
            vector<int> bankable_group;

            size_t j = i;
            while(j < cluster.original_reg_idx.size()) {
                const auto& current_inst = ffInstances[cluster.original_reg_idx[j]];
                int inst_clk = clk_net_map.count(current_inst.inst_name) ? clk_net_map[current_inst.inst_name] : -1;
                if(inst_clk != current_clk_net) break;
                
                int bits = type_bits_map[current_inst.type_name];
                if (total_bits + bits > max_bit) break;
                
                total_bits += bits;
                bankable_group.push_back(cluster.original_reg_idx[j]);
                x_coords.push_back(current_inst.x);
                y_coords.push_back(current_inst.y);
                j++;
            }
            
            // *** FIX: Only bank if the group has more than one FF ***
            // This prevents re-instantiating single flip-flops and causing pin name errors.
            if (bankable_group.size() > 1 && total_bits > 0) {
                Instance new_banked_inst;
                string new_inst_name = instName.first + to_string(new_idx_counter++);
                new_banked_inst.inst_name = new_inst_name;
                
                auto it = bits_minareaff_map.lower_bound(total_bits);
                new_banked_inst.type_name = (it != bits_minareaff_map.end()) ? it->second.second : bits_minareaff_map[max_bit].second;
                
                new_banked_inst.x = accumulate(x_coords.begin(), x_coords.end(), 0.0) / x_coords.size();
                new_banked_inst.y = accumulate(y_coords.begin(), y_coords.end(), 0.0) / y_coords.size();
                final_instances.push_back(new_banked_inst);
                
                for(int original_idx : bankable_group) processed_ffs[original_idx] = true;
                
                int new_bits = type_bits_map[new_banked_inst.type_name];
                int pin_offset = 0;
                for (int original_idx : bankable_group) {
                    const auto& old_inst = ffInstances[original_idx];
                    int old_bits = type_bits_map[old_inst.type_name];
                    
                    for (int k = 0; k < old_bits; ++k) {
                        string old_pin_suffix = (old_bits > 1) ? to_string(k) : "";
                        string new_pin_suffix = (new_bits > 1) ? to_string(pin_offset + k) : "";
                        output_mapping_lines.push_back(old_inst.inst_name + "/D" + old_pin_suffix + " map " + new_inst_name + "/D" + new_pin_suffix);
                        output_mapping_lines.push_back(old_inst.inst_name + "/Q" + old_pin_suffix + " map " + new_inst_name + "/Q" + new_pin_suffix);
                    }
                    output_mapping_lines.push_back(old_inst.inst_name + "/CLK map " + new_inst_name + "/CLK");
                    pin_offset += old_bits;
                }
            }
            i = j;
        }
    }

    // --- 4. Add Unchanged Flip-Flops and Generate Self-Mappings ---
    cout << "Handling unprocessed flip-flops..." << endl;
    for(size_t i = 0; i < ffInstances.size(); ++i) {
        if (!processed_ffs[i]) {
            const auto& inst = ffInstances[i];
            final_instances.push_back(inst);
            
            int bits = type_bits_map[inst.type_name];
            for(int k = 0; k < bits; ++k) {
                string pin_suffix = (bits > 1) ? to_string(k) : "";
                output_mapping_lines.push_back(inst.inst_name + "/D" + pin_suffix + " map " + inst.inst_name + "/D" + pin_suffix);
                output_mapping_lines.push_back(inst.inst_name + "/Q" + pin_suffix + " map " + inst.inst_name + "/Q" + pin_suffix);
            }
            output_mapping_lines.push_back(inst.inst_name + "/CLK map " + inst.inst_name + "/CLK");
        }
    }

    // --- 5. Placement Legalization ---
    cout << "\nCreating placement map and legalizing positions..." << endl;
    pmCreatePlacementMap();
    for (auto& instance : final_instances) {
        // Only legalize flip-flops, as gates are already considered fixed.
        if (splitString(instance.type_name).first == ffTag) {
            pmFFSearchPlacementMap(instance);
        }
    }
    cout << "Placement finished." << endl;

    // --- 6. Output Generation ---
    cout << "Writing output file: " << argv[2] << endl;
    ofstream fout(argv[2]);
    if (!fout) {
        cerr << "Error: Could not open output file '" << argv[2] << "'." << endl;
        return 1;
    }
    fout << "CellInst " << final_instances.size() << endl;
    for (const auto& instance : final_instances) {
        fout << "Inst " << instance.inst_name << " " << instance.type_name << " " << fixed << instance.x << " " << instance.y << endl;
    }
    for (const auto& line : output_mapping_lines) {
        fout << line << endl;
    }
    fout.close();
    cout << "Done." << endl;

    return 0;
}


//------------------------------
// Utility and Parsing Functions (mostly unchanged from your original)
//------------------------------

void InputParsing(char *fname) {
    string line;
    int clk_group_idx = 0;
    ifstream file(fname);
    if (!file.is_open()) { cerr << "Error: Could not open file '" << fname << "'." << endl; exit(1); }
    while (getline(file, line)) {
        istringstream iss(line);
        string key;
        iss >> key;
        if (key == "Alpha" || key == "Beta" || key == "Gamma") {
            double value; iss >> value; weights[key] = value;
        } else if (key == "FlipFlop") {
            FlipFlop ff;
            iss >> ff.bits >> ff.name >> ff.width >> ff.height >> ff.pinCount;
            double area = ff.width * ff.height;
            if (bits_minareaff_map.find(ff.bits) == bits_minareaff_map.end() || area < bits_minareaff_map[ff.bits].first) {
                bits_minareaff_map[ff.bits] = {area, ff.name};
            }
            possible_bits.insert(ff.bits);
            type_bits_map[ff.name] = ff.bits;
            for(int i = 0; i < ff.pinCount; ++i) getline(file, line);
            stFFValue values = {ff.bits, ff.width, ff.height, ff.pinCount, 0, 0};
            kmFFLib[ff.name] = values;
            if (ffTag == "") ffTag = splitString(ff.name).first;
        } else if (key == "Gate") {
            Gate gate;
            iss >> gate.name >> gate.width >> gate.height >> gate.pinCount;
            for (int i = 0; i < gate.pinCount; ++i) getline(file, line);
            stGateValue values = {gate.width, gate.height, gate.pinCount};
            kmGateLib[gate.name] = values;
        } else if (key == "NumInstances") {
            int instanceCount; iss >> instanceCount;
            for (int i = 0; i < instanceCount; i++) {
                getline(file, line);
                istringstream iss_inst(line);
                Instance instance; string tag;
                iss_inst >> tag >> instance.inst_name >> instance.type_name >> instance.x >> instance.y;
                if (splitString(instance.type_name).first == ffTag) {
                    ffInstances.push_back(instance);
                    if(instName.first == "") instName = splitString(instance.inst_name);
                } else {
                    gateInstances.push_back(instance);
                }
            }
        } else if (key == "Net") {
            Net net;
            iss >> net.name >> net.numPins;
            for (int i = 0; i < net.numPins; ++i) {
                getline(file, line);
                istringstream pinIss(line);
                string pin_type, connect_to;
                pinIss >> pin_type >> connect_to;
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
    file.close();
}

pair<string, int> splitString(const string& str) {
    size_t pos = str.find_last_not_of("0123456789");
    if (pos == string::npos || pos == str.length() - 1) return {str, -1};
    return {str.substr(0, pos + 1), stoi(str.substr(pos + 1))};
}

// --- Placement Functions (unchanged) ---
void pmGateSetPlacementMap(Instance gate){
    if (kmGateLib.find(gate.type_name) == kmGateLib.end()) return;
    stGateValue gateValue = kmGateLib[gate.type_name];
    double gateEndX = gate.x + gateValue.gateWidth;
    double gateEndY = gate.y + gateValue.gateHeight;
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

int SetNewInstPoint(Instance& new_instance){
    vector<vector<double>> possible_points;
    for(auto& pm : placementMaps){
        vector<double> possible_point;
        double new_x, new_y;
        if (new_instance.x <= pm.startX) new_x = pm.startX;
        else if (new_instance.x >= pm.endX) new_x = pm.endX - pm.siteWidth;
        else new_x = pm.startX + round((new_instance.x - pm.startX) / pm.siteWidth) * pm.siteWidth;

        if (new_instance.y <= pm.startY) new_y = pm.startY;
        else if (new_instance.y >= pm.endY) new_y = pm.endY - pm.siteHeight;
        else new_y = pm.startY + round((new_instance.y - pm.startY) / pm.siteHeight) * pm.siteHeight;
        
        possible_points.push_back({new_x, new_y});
    }
    auto min_it = min_element(possible_points.begin(), possible_points.end(), [&](const vector<double>& a, const vector<double>& b) {
        return (pow(a[0] - new_instance.x, 2) + pow(a[1] - new_instance.y, 2)) < (pow(b[0] - new_instance.x, 2) + pow(b[1] - new_instance.y, 2));
    });
    new_instance.x = (*min_it)[0];
    new_instance.y = (*min_it)[1];
    return distance(possible_points.begin(), min_it);
}

bool pmFFSearchPlacementMap(Instance& ffInstance){
    int idx = SetNewInstPoint(ffInstance);
    stPlacementMap &pm = placementMaps[idx];
    if (kmFFLib.find(ffInstance.type_name) == kmFFLib.end()) return false;
    stFFValue ffValue = kmFFLib[ffInstance.type_name];
    uint8_t bits_char = (uint8_t)ffValue.bits + '0';
    int sitesC = ceil(ffValue.ffWidth / pm.siteWidth);
    int sitesR = ceil(ffValue.ffHeight / pm.siteHeight);
    int c0 = floor((ffInstance.x - pm.startX) / pm.siteWidth);
    int r0 = floor((ffInstance.y - pm.startY) / pm.siteHeight);
    int ct = c0, rt = r0;
    bool bEmpty = pmIsEmpty(idx, ct, rt, sitesC, sitesR);
    if (!bEmpty) {
        int dmax = max({c0, pm.bmCols-c0-1, r0, pm.bmRows-r0-1});
        for (int d=1; !bEmpty && d<=dmax ; d++) {
            for(int side = 0; side < 4 && !bEmpty; ++side) {
                for(int step = 0; step < 2 * d; ++step) {
                    if (pmIsEmpty(idx, ct, rt, sitesC, sitesR)) { bEmpty = true; break; }
                    switch(side) {
                        case 0: ct++; break; case 1: rt++; break;
                        case 2: ct--; break; case 3: rt--; break;
                    }
                }
            }
        }
    }
    if (bEmpty) {
        ffInstance.x = pm.startX + pm.siteWidth * ct;
        ffInstance.y = pm.startY + pm.siteHeight * rt;
        for (int r=rt; r < rt + sitesR && r < pm.bmRows; r++) for (int c=ct; c < ct + sitesC && c < pm.bmCols; c++) pm.byteMap[r][c] = bits_char;
        return true;
    } else {
        cout << "!!NG Cannot Place FF: " << ffInstance.inst_name << endl;
        return false;
    }
}

void pmUpdatePlacementMap(int startR, int rowCnt) {
    stPlacementMap pm;
    pm.bmCols = placementRows[startR].totalNumOfSites;
    pm.bmRows = rowCnt;
    pm.startX = placementRows[startR].startX;
    pm.startY = placementRows[startR].startY;
    pm.siteWidth = placementRows[startR].siteWidth;
    pm.siteHeight = placementRows[startR].siteHeight;
    pm.endX = pm.startX + pm.siteWidth * pm.bmCols;
    pm.endY = pm.startY + pm.siteHeight * pm.bmRows;
    pm.byteMap.resize(pm.bmRows, vector<uint8_t>(pm.bmCols, 0));
    placementMaps.push_back(pm);
}

void pmCreatePlacementMap(void) {
    if (placementRows.empty()) return;
    sort(placementRows.begin(), placementRows.end(), [](const PlacementRow& a, const PlacementRow& b){
        if (a.startY != b.startY) return a.startY < b.startY;
        return a.startX < b.startX;
    });
    int totalRows = placementRows.size();
    int startR = 0;
    for (int r = 1; r < totalRows; r++) {
        if ((abs(placementRows[r-1].startX - placementRows[r].startX) > 1e-6)
          || (placementRows[r-1].totalNumOfSites != placementRows[r].totalNumOfSites)
          || (abs((placementRows[r-1].startY + placementRows[r-1].siteHeight) - placementRows[r].startY) > 1e-6 )) {
            pmUpdatePlacementMap(startR, r - startR);
            startR = r;
        }
    }
    pmUpdatePlacementMap(startR, totalRows - startR);
    for (const auto& instance : gateInstances) pmGateSetPlacementMap(instance);
}

