//
// Created by Xin Sui on 1/29/18.
//

#ifndef PROTEIN_REGISTRATION_DELAUNAYWRAPPER_H
#define PROTEIN_REGISTRATION_DELAUNAYWRAPPER_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>


struct Delaunay_info {
    typedef double FT;
    typedef Eigen::Matrix<FT, 6, Eigen::Dynamic, Eigen::ColMajor> Edge_matrix_type;
    typedef Eigen::Matrix<FT, 3, Eigen::Dynamic, Eigen::ColMajor> Vertex_matrix_type;
    typedef Eigen::Matrix<int, 4, Eigen::Dynamic, Eigen::ColMajor> Vertex_reference_type;

    Edge_matrix_type Edge_matrix;
    Vertex_matrix_type Vertex_matrix;
    Vertex_reference_type Vertex_reference;
};


struct DelaunayWrapper{
private:
    void read_asa (const boost::filesystem::path&, unsigned); //"surface/1a0ia-surf.asa"
    void read_pdb (const boost::filesystem::path&, unsigned);
    void Delaunay_Transformation();
public:
    DelaunayWrapper() {};
    DelaunayWrapper(const boost::filesystem::path&, unsigned);
    template <typename Derived> explicit DelaunayWrapper(const Eigen::MatrixBase<Derived>& mat) { Delaunay_data.Vertex_matrix = mat; Delaunay_Transformation();} ;
    Delaunay_info Delaunay_data;
};



#endif //PROTEIN_REGISTRATION_DELAUNAYWRAPPER_H
