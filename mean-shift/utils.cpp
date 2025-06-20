#include "utils.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>
#include <queue>
#include <stdexcept>
#include <chrono>
#include <omp.h> // For OpenMP parallelization

// Default kernel function if none is provided.
double gaussian_kernel(double dist_sq, double bandwidth) {
    if (bandwidth == 0) return 0;
    return exp(-0.5 * dist_sq / (bandwidth * bandwidth));
}


// --- K-D Tree for Fast Spatial Queries (Internal helper class) ---
namespace { 
using PqElement = std::pair<double, int>;
class KdTree {
    struct KdNode {
        int point_index;
        std::unique_ptr<KdNode> left = nullptr;
        std::unique_ptr<KdNode> right = nullptr;
    };
public:
    KdTree(const std::vector<MeanShift::Point>& points) : points_ref(points) {
        if (points.empty()) return;
        std::vector<int> indices(points.size());
        std::iota(indices.begin(), indices.end(), 0);
        root = build_recursive(indices.begin(), indices.end(), 0);
    }
    std::vector<int> find_knn(const MeanShift::Point& query, int k) const {
        std::priority_queue<PqElement> pq;
        find_knn_recursive(root.get(), query, k, 0, pq);
        std::vector<int> result;
        result.reserve(k);
        while (!pq.empty()) {
            result.push_back(pq.top().second);
            pq.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }
    std::vector<int> radius_search(const MeanShift::Point& query, double radius) const {
        std::vector<int> result;
        double radius_sq = radius * radius;
        radius_search_recursive(root.get(), query, radius_sq, 0, result);
        return result;
    }
private:
    const std::vector<MeanShift::Point>& points_ref;
    std::unique_ptr<KdNode> root = nullptr;
    double dist_sq(const MeanShift::Point& p1, const MeanShift::Point& p2) const {
        double d = 0.0;
        for (size_t i = 0; i < p1.size(); ++i) d += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        return d;
    }
    template<typename It>
    std::unique_ptr<KdNode> build_recursive(It start, It end, int depth) {
        if (start >= end) return nullptr;
        const int k = points_ref[0].size();
        int axis = depth % k;
        size_t n = std::distance(start, end);
        It mid = start + n / 2;
        std::nth_element(start, mid, end, [&](int a, int b) { return points_ref[a][axis] < points_ref[b][axis]; });
        auto node = std::make_unique<KdNode>();
        node->point_index = *mid;
        node->left = build_recursive(start, mid, depth + 1);
        node->right = build_recursive(mid + 1, end, depth + 1);
        return node;
    }
    void find_knn_recursive(const KdNode* node, const MeanShift::Point& query, int k, int depth, std::priority_queue<PqElement>& pq) const {
        if (!node) return;
        const int dim = query.size();
        int axis = depth % dim;
        double d2 = dist_sq(query, points_ref[node->point_index]);
        if (pq.size() < k) {
            pq.push({d2, node->point_index});
        } else if (d2 < pq.top().first) {
            pq.pop();
            pq.push({d2, node->point_index});
        }
        double axis_dist = query[axis] - points_ref[node->point_index][axis];
        auto first_branch = (axis_dist < 0) ? node->left.get() : node->right.get();
        auto second_branch = (axis_dist < 0) ? node->right.get() : node->left.get();
        find_knn_recursive(first_branch, query, k, depth + 1, pq);
        if (pq.size() < k || (axis_dist * axis_dist) < pq.top().first) {
            find_knn_recursive(second_branch, query, k, depth + 1, pq);
        }
    }
    void radius_search_recursive(const KdNode* node, const MeanShift::Point& query, double radius_sq, int depth, std::vector<int>& result) const {
        if (!node) return;
        const int dim = query.size();
        int axis = depth % dim;
        double d2 = dist_sq(query, points_ref[node->point_index]);
        if (d2 <= radius_sq) result.push_back(node->point_index);
        double axis_dist = query[axis] - points_ref[node->point_index][axis];
        auto first_branch = (axis_dist < 0) ? node->left.get() : node->right.get();
        auto second_branch = (axis_dist < 0) ? node->right.get() : node->left.get();
        radius_search_recursive(first_branch, query, radius_sq, depth, result);
        if ((axis_dist * axis_dist) <= radius_sq) {
            radius_search_recursive(second_branch, query, radius_sq, depth, result);
        }
    }
};
}

// --- MeanShift Class Implementation ---

// *** FIX: Added definitions for the constructors here. ***
MeanShift::MeanShift() {
    set_kernel(nullptr); // Call the private method
}

MeanShift::MeanShift(double (*_kernel_func)(double,double)) {
    set_kernel(_kernel_func);
}

void MeanShift::set_kernel(double (*_kernel_func)(double,double)) {
    kernel_func = (_kernel_func == nullptr) ? gaussian_kernel : _kernel_func;
}

std::vector<Cluster> MeanShift::cluster(const std::vector<Point>& points, double kernel_bandwidth) {
    int default_m_neighbors = 15;
    return cluster_optimized(points, kernel_bandwidth, default_m_neighbors);
}
void MeanShift::legalization(const std::vector<Point>& original_placement, const std::vector<Point>& new_placement) {
    std::cout << "Legalization function called. (Not implemented)" << std::endl;
}
std::vector<Cluster> MeanShift::cluster_optimized(const std::vector<Point>& points, double kernel_bandwidth, int m_neighbors) {
    if (points.empty()) return {};
    double search_radius = kernel_bandwidth * 4.0;
    std::vector<double> var_h = calculate_variable_bandwidth(points, m_neighbors);
    std::vector<Point> shifted_points = meanshift(points, var_h, search_radius);
    return group_clusters(points, shifted_points);
}
std::vector<double> MeanShift::calculate_variable_bandwidth(const std::vector<Point>& points, int m) {
    KdTree tree(points);
    std::vector<double> var_h(points.size());
    #pragma omp parallel for
    for (int i = 0; i < points.size(); ++i) {
        std::vector<int> neighbors = tree.find_knn(points[i], m + 1);
        if (neighbors.size() > m) {
            const Point& neighbor_pt = points[neighbors[m]];
            double dist_sq = 0.0;
            for(size_t d = 0; d < points[i].size(); ++d) dist_sq += (points[i][d] - neighbor_pt[d]) * (points[i][d] - neighbor_pt[d]);
            var_h[i] = std::sqrt(dist_sq);
        } else if (!neighbors.empty()){
            const Point& last_neighbor = points[neighbors.back()];
            double dist_sq = 0.0;
            for(size_t d = 0; d < points[i].size(); ++d) dist_sq += (points[i][d] - last_neighbor[d]) * (points[i][d] - last_neighbor[d]);
            var_h[i] = std::sqrt(dist_sq);
        } else {
            var_h[i] = 1.0;
        }
    }
    return var_h;
}
std::vector<MeanShift::Point> MeanShift::meanshift(const std::vector<Point>& points, const std::vector<double>& var_h, double search_radius) {
    KdTree tree(points);
    std::vector<Point> shifted_points = points;
    const double CONVERGENCE_EPSILON_SQ = 1e-4 * 1e-4;
    const int MAX_ITERATIONS = 100;
    const int DIMS = points[0].size();
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        double max_shift_dist_sq = 0;
        std::vector<Point> next_shifted_points = shifted_points;
        #pragma omp parallel for reduction(max:max_shift_dist_sq)
        for (int i = 0; i < points.size(); ++i) {
            Point current_pos = shifted_points[i];
            double total_weight = 0;
            Point sum_points(DIMS, 0.0);
            std::vector<int> neighbor_indices = tree.radius_search(current_pos, search_radius);
            if (neighbor_indices.empty()) continue;
            for (int neighbor_idx : neighbor_indices) {
                const Point& original_neighbor_pt = points[neighbor_idx];
                double dist_sq = 0.0;
                for(int d = 0; d < DIMS; ++d) dist_sq += (current_pos[d] - original_neighbor_pt[d]) * (current_pos[d] - original_neighbor_pt[d]);
                double weight = this->kernel_func(dist_sq, var_h[neighbor_idx]);
                for(int d = 0; d < DIMS; ++d) sum_points[d] += original_neighbor_pt[d] * weight;
                total_weight += weight;
            }
            if (total_weight > 1e-9) {
                Point new_pos(DIMS);
                for(int d=0; d < DIMS; ++d) new_pos[d] = sum_points[d] / total_weight;
                double shift_dist_sq = 0.0;
                for(int d=0; d < DIMS; ++d) shift_dist_sq += (new_pos[d] - current_pos[d]) * (new_pos[d] - current_pos[d]);
                next_shifted_points[i] = new_pos;
                if (shift_dist_sq > max_shift_dist_sq) max_shift_dist_sq = shift_dist_sq;
            }
        }
        shifted_points = next_shifted_points;
        DBGPRN("Iteration %d: Max shift distance = %f\n", iter, std::sqrt(max_shift_dist_sq));
        if (max_shift_dist_sq < CONVERGENCE_EPSILON_SQ) {
            std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
            break;
        }
    }
    return shifted_points;
}
std::vector<Cluster> MeanShift::group_clusters(const std::vector<Point>& original_points, const std::vector<Point>& shifted_points) {
    std::vector<Cluster> clusters;
    const double CLUSTER_EPSILON_SQ = 0.5 * 0.5;
    if (original_points.empty()) return clusters;
    const int DIMS = original_points[0].size();
    for (int i = 0; i < shifted_points.size(); ++i) {
        int target_cluster_idx = -1;
        for (size_t j = 0; j < clusters.size(); ++j) {
            double dist_sq = 0.0;
            for(int d=0; d < DIMS; ++d) dist_sq += (shifted_points[i][d] - clusters[j].mode[d]) * (shifted_points[i][d] - clusters[j].mode[d]);
            if (dist_sq <= CLUSTER_EPSILON_SQ) {
                target_cluster_idx = j;
                break;
            }
        }
        if (target_cluster_idx == -1) {
            Cluster new_cluster;
            new_cluster.mode = shifted_points[i];
            new_cluster.original_points.push_back(original_points[i]);
            new_cluster.shifted_points.push_back(shifted_points[i]);
            new_cluster.original_reg_idx.push_back(i);
            new_cluster.shifted_reg_idx.push_back(i);
            clusters.push_back(new_cluster);
        } else {
            clusters[target_cluster_idx].original_points.push_back(original_points[i]);
            clusters[target_cluster_idx].shifted_points.push_back(shifted_points[i]);
            clusters[target_cluster_idx].original_reg_idx.push_back(i);
            clusters[target_cluster_idx].shifted_reg_idx.push_back(i);
        }
    }
    return clusters;
}
