//
// Created by Xin Sui on 2/8/18.
//

#include "dSSD_sort.h"
#include "SSD_sort.h"
#include "hungarian_algorithm.h"
#include "DelaunayWrapper.h"
#include "ThreadPool.h"
#include "file_IO.h"
#include <iostream>
#include "configuration.h"

double one_vs_one(const boost::filesystem::path& target, const boost::filesystem::path& source, const unsigned length, const unsigned dSSD_cutoff,
                const unsigned SSD_cutoff)
{
    const DelaunayWrapper source_protein(source, length);
    const DelaunayWrapper target_protein(target, length);

    const std::vector<dSSD_vec> dSSD_result = dSSD_sort(target_protein.Delaunay_data.Edge_matrix, source_protein.Delaunay_data.Edge_matrix, dSSD_cutoff);
    const std::vector<Ini_Trans> SSD_result = SSD_sort(dSSD_result, target_protein.Delaunay_data.Vertex_matrix, source_protein.Delaunay_data.Vertex_matrix,
                                                 target_protein.Delaunay_data.Vertex_reference, source_protein.Delaunay_data.Vertex_reference, SSD_cutoff);

    const Alignment_result ret = Hungarian_algorithm(SSD_result, target_protein.Delaunay_data.Vertex_matrix, source_protein.Delaunay_data.Vertex_matrix, 20, false);

    return ret.fitness_score;
}

Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> one_vs_all(const boost::filesystem::path& target, const boost::filesystem::path& directory, const unsigned length, const unsigned dSSD_cutoff,
                               const unsigned SSD_cutoff)
{
//    const DelaunayWrapper target_protein(target, length);
    const std::vector<boost::filesystem::path> sources = scan_directory(directory);
    int number_of_proteins = sources.size();
    Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> fitness_scores;
    fitness_scores.resize(1, number_of_proteins);

    ThreadPool::ParallelFor(0, number_of_proteins, [&] (int j) {
        fitness_scores[j] = one_vs_one(target, sources[j], length, dSSD_cutoff, SSD_cutoff);
    });

//    for (int j = 0; j != number_of_proteins; ++j){
//        fitness_scores[j] = one_vs_one(target, sources[j], length, dSSD_cutoff, SSD_cutoff);
//    }

    return fitness_scores;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> all_vs_all(const boost::filesystem::path& directory,
                                                                                  const unsigned length, const unsigned dSSD_cutoff,
                                                                                  const unsigned SSD_cutoff)
{
    const std::vector<boost::filesystem::path> proteins = scan_directory(directory);
    int number_of_proteins = proteins.size();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> fitness_scores(number_of_proteins, number_of_proteins);
    fitness_scores.setZero();

    for(int i = 0; i != number_of_proteins; ++i){
        if( number_of_proteins > 10 && i % (number_of_proteins / 10) == 0)
            std::cout << i << " of " << number_of_proteins << " computed. " << std::endl;
        ThreadPool::ParallelFor(i, number_of_proteins, [&] (int j) {
            fitness_scores(i, j) = one_vs_one(proteins[i], proteins[j], length, dSSD_cutoff, SSD_cutoff);
        });
//        for(int j = i; j!= number_of_proteins; ++j)
//            fitness_scores(i,j) = one_vs_one(proteins[i], proteins[j], length, dSSD_cutoff, SSD_cutoff);
    }
    return fitness_scores;
}

int main(){
//    boost::filesystem::path data1 = "1bf4A-10.pdb";
//    boost::filesystem::path data2 = "1mxbA-10.pdb";
//    boost::filesystem::path dataa = "1bf4A-30.pdb";
//    boost::filesystem::path datab = "1mxbA-30.pdb";
//    boost::filesystem::path directory = "/Users/xinsui/CLionProjects/protein_registration/dataset/2.test_dataset/1.10";
//    boost::filesystem::path output = "result.txt";

    auto start = std::chrono::steady_clock::now();

//    std::cout << one_vs_one(directory / data1, directory / data2, 10, 100, 50) << std::endl;
//    std::cout << one_vs_one(directory / data2, directory / data1, 10, 100, 50) << std::endl;
//    const Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> result = one_vs_all(directory / data2, directory, 10, 100, 50, .5);

    ConfigData Cfg = configuration("config.txt");

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> result = all_vs_all(Cfg.data_directory, Cfg.N, Cfg.dSSD_cutoff, Cfg.SSD_cutoff);

    auto end = std::chrono::steady_clock::now();
    double time = std::chrono::duration_cast<std::chrono::seconds> (end-start).count();

    std::cout << "Computation completed in " << time << " seconds." << std::endl;

    Eigen::IOFormat NumpyFormat(Eigen::StreamPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
    boost::filesystem::ofstream os;
    os.open(Cfg.data_directory / Cfg.out_file_name, std::ios::out);
    os << result.format(NumpyFormat);
    os.close();


    return 0;
}
