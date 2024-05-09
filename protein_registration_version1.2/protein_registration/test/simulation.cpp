//
// Created by Xin Sui on 1/29/18.
//

#include <array>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iomanip>
#include <ios>
#include "../src/dSSD_sort.h"
#include "../src/SSD_sort.h"
#include "../src/hungarian_algorithm.h"
#include "../src/DelaunayWrapper.h"
#include "../src/ThreadPool.h"

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

double align_this(const int number_of_points, const double noise_level,
                      const Eigen::Matrix<double,3,3,Eigen::ColMajor>& rand_rotation, const int dSSD_cutoff, const int SSD_cutoff){

    static Eigen::Matrix<int, 12, 4, Eigen::RowMajor> Permu;
    Permu << 0,1,2,3, 0,2,3,1, 0,3,1,2, 1,0,3,2, 1,2,0,3, 1,3,2,0,
            2,0,1,3, 2,1,3,0, 2,3,0,1, 3,0,2,1, 3,1,0,2, 3,2,1,0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> uniform(0.,1.);
    std::uniform_int_distribution<int> rint(0,number_of_points - 1);

    Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand2;
    Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand1;
    rand2 = rand2.Random(3,number_of_points);
    rand1.resize(3,number_of_points);
    // Apply scaling to rand2
    // Apply rotation on rand2 to get rand1
    rand1 = rand_rotation * rand2;
    // Apply transposition
    Eigen::VectorXi indices = Eigen::VectorXi::LinSpaced(number_of_points, 0, number_of_points);
    std::shuffle(indices.data(), indices.data() + rand1.cols(), gen);
    // Add noise to our target
    Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> noise;
    noise = noise.Random(3,number_of_points) * noise_level;
    rand1 += noise;


    Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean_a = rand1.rowwise().mean();
    Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean_b = rand2.rowwise().mean();
    rand1.colwise() -= mean_a;
    rand2.colwise() -= mean_b; // You don't copy. Should be a way there.
    // TODO: write unary expression for extracting, and binary expression for before svd.
    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > svd(rand1 * rand2.transpose(),
            Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix<double, 3, 3> RotMat = svd.matrixU() * svd.matrixV().transpose();

    rand1 = rand1 * indices.asPermutation();





    DelaunayWrapper Pro1(rand1);
    DelaunayWrapper Pro2(rand2);

    std::vector<dSSD_vec> dSSD_result = dSSD_sort(Pro1.Delaunay_data.Edge_matrix, Pro2.Delaunay_data.Edge_matrix, dSSD_cutoff);
    std::vector<Ini_Trans> SSD_result = SSD_sort(dSSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix,
                                                 Pro1.Delaunay_data.Vertex_reference, Pro2.Delaunay_data.Vertex_reference, SSD_cutoff);

    //return Hungarian_algorithm(SSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix, 20, 0);
    Alignment_result this_result = Hungarian_algorithm(SSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix, 20, 0);


    return (RotMat - this_result.Rotation_matrix).norm();

}

double standard_deviation(const std::vector<double>& data_pool, const double mean){
    double var = 0.f;
    for(auto iter = data_pool.begin(); iter != data_pool.end(); ++iter){
        var += (*iter - mean) * (*iter - mean);
    }
    return std::sqrt(var/data_pool.size());
}


void align(const int iters, const double noise_level, const int number_of_points, int dSSD_cutoff, int SSD_cutoff){
    
    auto start = std::chrono::steady_clock::now();

    std::vector<double> pool_recovery_score(iters);
    std::vector<int> pool_converged(iters);
    std::vector<double> pool_fitness_score(iters);
    double iters_f = iters;

    ThreadPool::ParallelFor(0, iters, [&] (int j) {
        Eigen::Matrix<double,3,3,Eigen::ColMajor> rand_rotation = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();
        pool_recovery_score[j] = align_this(number_of_points, noise_level, rand_rotation, dSSD_cutoff, SSD_cutoff);
    });

    /*
     *     ThreadPool::ParallelFor(0, iters, [&] (int j) {
        Eigen::Matrix<double,3,3,Eigen::ColMajor> rand_rotation = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();
        pool_recovery_score[j] = align_this(number_of_points, noise_level, rand_rotation, dSSD_cutoff, SSD_cutoff);
    });

     */

/*
    ThreadPool::ParallelFor(0, iters, [&] (int j) {
        Eigen::Matrix<double,3,3,Eigen::ColMajor> rand_rotation = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();
        Alignment_result this_result = align_this(number_of_points, noise_level, rand_rotation, dSSD_cutoff, SSD_cutoff);
        pool_recovery_score[j] = (this_result.Rotation_matrix - rand_rotation).squaredNorm();
        pool_converged[j] = this_result.has_converged;
        pool_fitness_score[j] = this_result.fitness_score;
    });
*/
    auto end = std::chrono::steady_clock::now();
    double avg_time = std::chrono::duration_cast<std::chrono::seconds> (end-start).count() / iters_f;
    //int number_of_convergence = std::accumulate(pool_converged.begin(), pool_converged.end(), 0);
    double avg_rotation_recovery = std::accumulate(pool_recovery_score.begin(), pool_recovery_score.end(), 0.f)/ iters_f;
    double rotation_recovery_std = standard_deviation(pool_recovery_score, avg_rotation_recovery);
    //double avg_fitness_score = std::accumulate(pool_fitness_score.begin(), pool_fitness_score.end(), 0.f)/ iters_f;
    //double fitness_score_std = standard_deviation(pool_fitness_score, avg_fitness_score);


    std::cout << number_of_points << ' ' << noise_level << ' ' << avg_rotation_recovery << ' '
              << rotation_recovery_std << ' ' << avg_time << std::endl;


    //std::cout << std::setw(16) << number_of_points << std::setw(15) << noise_level << std::setw(15) << number_of_convergence << std::setw(25) << avg_rotation_recovery << std::setw(22)
    //          << rotation_recovery_std << std::setw(21) << avg_fitness_score << std::setw(22) << fitness_score_std << std::setw(13) << avg_time << std::endl;
}


int main(){

    const int iters = 0200;
    std::array<int, 4> number_of_points = {10, 50, 100, 200};
    std::array<double, 6> noise_level = {0, .05, .1, .15, .2, .25};
    int dSSD_cutoff, SSD_cutoff;

    std::cout << "number of points    noise level    # converged    avg rotation recovery    standard deviation       avg. time" << std::endl;
    for(int cutoff_strat = 0; cutoff_strat != 5; ++cutoff_strat) {
        switch (cutoff_strat){
            case 0:
                dSSD_cutoff = 100;
                SSD_cutoff = 50;
                break;
            case 1:
                dSSD_cutoff = 200;
                SSD_cutoff = 100;
                break;
            case 2:
                dSSD_cutoff = 400;
                SSD_cutoff = 200;
                break;
            case 3:
                dSSD_cutoff = 800;
                SSD_cutoff = 400;
                break;
            case 4:
                dSSD_cutoff = 1200;
                SSD_cutoff = 1000;
                break;
            default:
                return 1;
        }
        std::cout << "Cutoff Parameters: dRMSD cutoff=" << dSSD_cutoff << ", RMSD cutoff=" << SSD_cutoff << std::endl;
        std::cout << "=====================================================================================================================================================" << std::endl;
        for (auto this_noise_level = noise_level.begin(); this_noise_level != noise_level.end(); ++this_noise_level) {
            for (auto this_num_o_points = number_of_points.begin();
                 this_num_o_points != number_of_points.end(); ++this_num_o_points) {
                align(iters, *this_noise_level, *this_num_o_points, dSSD_cutoff, SSD_cutoff);
            }
            std::cout << std::endl;
        }
        std::cout << "=========================================================================================================================================" << std::endl;
    }
    return 0;
}