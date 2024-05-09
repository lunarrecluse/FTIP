//
// Created by Xin Sui on 2/20/18.
//

#ifndef PROTEIN_REGISTRATION_HUNGARIAN_ALGORITHM_H
#define PROTEIN_REGISTRATION_HUNGARIAN_ALGORITHM_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include "SSD_sort.h"

struct Alignment_result{
    bool has_converged;
    double fitness_score;
    Eigen::MatrixXd Rotation_matrix;
};

Alignment_result Hungarian_algorithm(const std::vector<Ini_Trans>&, const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& ,
                                     const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& , int , bool);



#endif //PROTEIN_REGISTRATION_HUNGARIAN_ALGORITHM_H
