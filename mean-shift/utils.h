#ifndef MEAN_SHIFT_H
#define MEAN_SHIFT_H

#include <vector>
#include <cstdio> // For printf

// Using a commented-out pragma as requested, with standard include guards.
// #pragma once 

// Debugging macro as specified.
#define _DEBUG  0
#if _DEBUG
#define DBGPRN(fmt,args...)          printf(fmt, ##args)
#else
#define DBGPRN(fmt,args...)          ((void)0)
#endif

// Represents a cluster of flip-flops (registers) after the clustering process.
struct Cluster {
    std::vector<double> mode;
    std::vector<std::vector<double>> original_points;
    std::vector<std::vector<double>> shifted_points;
    std::vector<int> original_reg_idx;
    std::vector<int> shifted_reg_idx;
};

// A high-performance, parallelized implementation of the Mean Shift clustering algorithm.
class MeanShift {
public:
    using Point = std::vector<double>;

    // *** FIX: Constructor is now only DECLARED here. ***
    // The definition has been moved to the .cpp file.
    MeanShift();
    
    // This constructor was unused, but for good practice, it should also be defined in the .cpp file.
    MeanShift(double (*_kernel_func)(double,double));

    // --- Main Public API ---
    std::vector<Cluster> cluster(const std::vector<Point>& points, double kernel_bandwidth);
    void legalization(const std::vector<Point>& original_placement, const std::vector<Point>& new_placement);

private:
    double (*kernel_func)(double dist_sq, double bandwidth);
    void set_kernel(double (*_kernel_func)(double, double));
    std::vector<Cluster> cluster_optimized(const std::vector<Point>& points, double kernel_bandwidth, int m_neighbors);
    std::vector<Point> meanshift(const std::vector<Point>& points, const std::vector<double>& var_h, double search_radius);
    std::vector<double> calculate_variable_bandwidth(const std::vector<Point>& points, int m);
    std::vector<Cluster> group_clusters(const std::vector<Point>& original_points, const std::vector<Point>& shifted_points);
};

#endif // MEAN_SHIFT_H
