//
// Created by Xin Sui on 2/8/18.
//
/*
#include <iostream>
#include "../src/dSSD_sort.h"
#include "../src/SSD_sort.h"
#include "../src/Icp_algorithm.h"
#include "../src/DelaunayWrapper.h"
#include "../src/ThreadPool.h"

int dSSD_cutoff = 100;
int SSD_cutoff = 50;
double max_corres_dist = 0.5;

int main(){


    DelaunayWrapper Pro1("/Users/xinsui/CLionProjects/protein_registration/dataset/2.test_dataset/1.10/1a22A-10.pdb", 10, "pdb");
    DelaunayWrapper Pro2("/Users/xinsui/CLionProjects/protein_registration/dataset/2.test_dataset/1.10/1axiA-10.pdb", 10, "pdb");

    std::vector<dSSD_vec> dSSD_result = dSSD_sort(Pro1.Delaunay_data.Edge_matrix, Pro2.Delaunay_data.Edge_matrix, dSSD_cutoff);
    std::vector<Ini_Trans> SSD_result = SSD_sort(dSSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix,
                                                 Pro1.Delaunay_data.Vertex_reference, Pro2.Delaunay_data.Vertex_reference, SSD_cutoff);

    ICP_result ret = Icp_algorithm(SSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix, max_corres_dist, 20, 0);

    std::cout << ret.Transformation_matrix << std::endl;

}
 */

#include <string>
#include <iostream>
#include "file_IO.h"


std::vector<boost::filesystem::path> write_protein_names(const boost::filesystem::path& directory)
{
    long pdb_asa_file_count = 0;
    std::vector<boost::filesystem::path> proteins;
    for (auto& p : boost::filesystem::directory_iterator(directory)) {
        if (p.path().extension() == ".pdb" || p.path().extension() == ".asa") {
            proteins.emplace_back(p.path());
            pdb_asa_file_count++;
        }
    }

    std::cout << pdb_asa_file_count << " pdb/asa files found." << std::endl;

    // Write protein name and order to txt file
    long j = 0;
    boost::filesystem::ofstream log;
    boost::filesystem::path txt_file_name = "protein_names.txt";
    log.open(directory / txt_file_name, std::ios::out);

    for(auto& p : proteins){
        log << j++ << ' ' << p.filename() << std::endl;
    }
    log.close();

    return proteins;
}

std::vector<boost::filesystem::path> read_protein_names(const boost::filesystem::path& directory)
{
    std::vector<boost::filesystem::path> proteins;

    // Write protein name and order to txt file
    long j = 0;
    boost::filesystem::ifstream log;
    boost::filesystem::path txt_file_name = directory / "protein_names.txt";

    log.open(txt_file_name, std::ios::in);
    long number;
    boost::filesystem::path file_name;
    while(log >> number >> file_name)
        proteins.emplace_back(directory / file_name);

    std::cout << proteins.size() << " pdb/asa files found." << std::endl;
    return proteins;
}

std::vector<boost::filesystem::path> scan_directory(const boost::filesystem::path& directory)
{

    const boost::filesystem::path txt_file = directory / "protein_names.txt";
    std::vector<boost::filesystem::path> proteins;
    if (!boost::filesystem::exists(txt_file))
        proteins = write_protein_names(directory);
    else
        proteins = read_protein_names(directory);
    return proteins;
}
