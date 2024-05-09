//
// Created by Xin Sui on 2/20/18.
//

#include "hungarian_algorithm.h"
#include <dlib/optimization/max_cost_assignment.h>
#include <dlib/matrix/matrix_mat_abstract.h>

Alignment_result Hungarian_algorithm(const std::vector<Ini_Trans>& trans, const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& protein_target,
                                    const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& protein_source, int max_iter, bool DEBUG=0)
{

    const Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean_target = protein_target.rowwise().mean();
    const Eigen::MatrixXd target_demean = protein_target.colwise() - mean_target;
    const Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean_source = protein_source.rowwise().mean();

    Alignment_result result;
    result.fitness_score = std::numeric_limits<double>::max();

    for (std::vector<Ini_Trans>::const_iterator it = trans.begin(); it != trans.end(); ++it) // Iterate over Init
    {

        Eigen::MatrixXd source_demean = protein_source.colwise() - mean_source;
        source_demean = it->RotMat * source_demean;
        Eigen::Matrix3d accumulated_rotation_matrix = it->RotMat;

        int iters = 0;
        double residual_new = std::numeric_limits<double>::max();
        double residual = 0.;

        while( residual_new > 1e-6 && abs(residual - residual_new) > 1e-6 && iters < max_iter) {
            residual = residual_new;

            Eigen::MatrixXd distances = apply_kernel(source_demean, target_demean, my_kernel);  // job(tgt),col and person(src),row
            distances = distances * (- 128 * 128);
            std::vector<long> assignment = dlib::max_cost_assignment(dlib::mat(distances.cast<long>().eval()));


            Eigen::MatrixXd source_temp = source_demean;
            for (long j = 0; j != assignment.size(); ++j) {
                source_temp.col(assignment[j]) = source_demean.col(j);
            }

            Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > svd(
                    target_demean * source_temp.transpose(),
                    Eigen::ComputeFullU | Eigen::ComputeFullV);
            const Eigen::DiagonalMatrix<double, 3> diag(1, 1, svd.singularValues()[0] > 0 ? 1 : -1);
            const Eigen::Matrix<double, 3, 3, Eigen::ColMajor> RotMat = svd.matrixU() * diag * svd.matrixV().transpose();

            accumulated_rotation_matrix = RotMat * accumulated_rotation_matrix;
            source_demean = RotMat * source_temp;

            residual_new = (target_demean - source_demean).squaredNorm();
            ++iters;
        }
        if (residual_new < result.fitness_score * result.fitness_score)
        {
            result.fitness_score = std::sqrt(residual_new);
            result.has_converged = (iters < max_iter || abs(residual - residual_new) <= 1e-6); // Need to consider convergence at last iter
            result.Rotation_matrix = accumulated_rotation_matrix;
        }

    }

    return result;
}

