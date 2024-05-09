//
// Created by Xin Sui on 2/9/18.
//

//
#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/icp_nl.h>
#include <pcl/registration/transformation_estimation_svd_scale.h>
#include <pcl/common/pca.h>
#include <pcl/common/transforms.h>
#include <pcl/registration/correspondence_rejection_one_to_one.h>
#include <Eigen/Core>
#include <Eigen/Dense>
/*
typedef pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ, double> ICP;


double Icp_algorithm(const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& pro1,
                         const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& pro2,
                         double max_corres_dist = .00001, int max_iter = 20, bool DEBUG=0)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_1 (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2 (new pcl::PointCloud<pcl::PointXYZ>);
    Eigen::Matrix<double,4,4,Eigen::ColMajor> ret;

    // Fill in the CloudIn data
    cloud_1->width    = pro1.cols();
    cloud_1->height   = 1; // Our data is unorganized.
    cloud_1->is_dense = true; // TRUE when no points is NaN or invalid
    cloud_1->points.resize (cloud_1->width * cloud_1->height); // We can use this technique as well, build a matrix of 1x1 and resize it later XD
    cloud_2->width    = pro2.cols();
    cloud_2->height   = 1;
    cloud_2->is_dense = true; // TRUE when no points is NaN or invalid
    cloud_2->points.resize (cloud_2->width * cloud_2->height);
    for (unsigned i = 0; i != cloud_1->points.size (); ++i)
    {
        cloud_1->points[i].x = pro1(0,i);
        cloud_1->points[i].y = pro1(1,i);
        cloud_1->points[i].z = pro1(2,i);
        cloud_2->points[i].x = pro2(0,i);
        cloud_2->points[i].y = pro2(1,i);
        cloud_2->points[i].z = pro2(2,i);
    } // TODO: Copy overhead. But Let's deal with it now.

    pcl::registration::CorrespondenceEstimation<pcl::PointXYZ, pcl::PointXYZ, double> est;
    est.setInputSource(cloud_2);
    est.setInputTarget(cloud_1);
    pcl::Correspondences all_correspondences;
// Determine all reciprocal correspondences
    est.determineCorrespondences (all_correspondences);

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_1t (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2t (new pcl::PointCloud<pcl::PointXYZ>);

    pcl::PCA<pcl::PointXYZ> pca;
    pca.setInputCloud(cloud_1);
    Eigen::Vector3f ev_A = pca.getEigenValues();
    float scale_1 = 1 / std::sqrt(ev_A[0]);
    Eigen::DiagonalMatrix<float, 4> transform_1;
    transform_1.diagonal() << scale_1, scale_1, scale_1, 1;
    pcl::transformPointCloud (*cloud_1, *cloud_1t, transform_1);

    pca.setInputCloud(cloud_2);
    Eigen::Vector3f ev_B = pca.getEigenValues();
    float scale_2 = 1 / std::sqrt(ev_B[0]);
    Eigen::DiagonalMatrix<float, 4> transform_2;
    transform_2.diagonal() << scale_2, scale_2, scale_2, 1;
    pcl::transformPointCloud (*cloud_2, *cloud_2t, transform_2);



    ICP icp;
    icp.setInputSource(cloud_2t);
    icp.setInputTarget(cloud_1t);

    pcl::registration::CorrespondenceRejectorOneToOne::Ptr rej (new pcl::registration::CorrespondenceRejectorOneToOne);
    icp.addCorrespondenceRejector(rej);

    pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZ>::Ptr rej_samp (new pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZ>);
    icp.addCorrespondenceRejector (rej_samp);

    icp.setMaximumIterations(888);
    icp.setEuclideanFitnessEpsilon(1e-6);
    icp.setTransformationEpsilon(1e-6);
 //   icp.setMaxCorrespondenceDistance(5);

//    pcl::registration::TransformationEstimationSVDScale<pcl::PointXYZ, pcl::PointXYZ, double>::Ptr trans_svd (new pcl::registration::TransformationEstimationSVDScale<pcl::PointXYZ, pcl::PointXYZ, double>);
//    icp.setTransformationEstimation (trans_svd);
    pcl::PointCloud<pcl::PointXYZ> Final;

    Eigen::Matrix<double,4,4,Eigen::ColMajor> trans_init;
    icp.align(Final); // TODO: Catch ERROR and process it

    Eigen::IOFormat NumpyFormat(Eigen::StreamPrecision, 0, ", ", ",\n", "[", "]", "[", "]");

    std::cout << "Converged?" << std::endl << icp.hasConverged() << std::endl;
    std::cout << "getEuclideanFitnessEpsilon" << std::endl << icp.getEuclideanFitnessEpsilon() << std::endl;


    std::cout << "target" << std::endl << cloud_1t->getMatrixXfMap().format(NumpyFormat) << std::endl << std::endl;

    std::cout << "transformed" << std::endl << Final.getMatrixXfMap().format(NumpyFormat) << std::endl << std::endl;

    std::cout << icp.getFinalTransformation() << std::endl;

    return icp.getFitnessScore();

}
*/

typedef pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ, double> ICP;


pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud_demean_scaling(const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& protein, double* scale){
    const Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean = protein.rowwise().mean();
    Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor> protein_demean = protein.colwise() - mean;
    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > svd(protein_demean * protein_demean.transpose(),
                                                                        Eigen::ComputeFullU | Eigen::ComputeFullV);
    *scale = std::sqrt(svd.singularValues()[0]);
    protein_demean /= *scale;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    cloud->width    = protein_demean.cols();
    cloud->height   = 1; // Our data is unorganized.
    cloud->is_dense = true; // TRUE when no points is NaN or invalid
    cloud->points.resize (cloud->width * cloud->height); // We can use this technique as well, build a matrix of 1x1 and resize it later XD
    for (unsigned i = 0; i != cloud->points.size(); ++i)
    {
        cloud->points[i].x = protein_demean(0,i);
        cloud->points[i].y = protein_demean(1,i);
        cloud->points[i].z = protein_demean(2,i);
    } // TODO: Copy overhead. But Let's deal with it now.

    return cloud;
}


double Icp_algorithm(const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& pro1,
                         const Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>& pro2,
                         double max_corres_dist=99.9, int max_iter=10, bool DEBUG=0)
{
    Eigen::Matrix<double,4,4,Eigen::ColMajor> ret;

    double scale1;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_1 = point_cloud_demean_scaling(pro1, &scale1);
    std::cout << "Pro1A" << std::endl << cloud_1->getMatrixXfMap() << std::endl << std::endl;

    double scale2;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2 = point_cloud_demean_scaling(pro2, &scale2);
    std::cout << "Pro2A" << std::endl << cloud_2->getMatrixXfMap() << std::endl << std::endl;

    // ICP settings
    ICP icp;
    icp.setInputSource(cloud_2);
    icp.setInputTarget(cloud_1);
    icp.setMaximumIterations(max_iter);
    icp.setMaxCorrespondenceDistance(max_corres_dist);
//        pcl::registration::TransformationEstimationSVDScale<pcl::PointXYZ, pcl::PointXYZ, double>::Ptr trans_svd (new pcl::registration::TransformationEstimationSVDScale<pcl::PointXYZ, pcl::PointXYZ, double>);
//        icp.setTransformationEstimation (trans_svd);
//        pcl::registration::CorrespondenceRejectorOneToOne::Ptr rej (new pcl::registration::CorrespondenceRejectorOneToOne);
//        icp.addCorrespondenceRejector(rej);
    icp.setUseReciprocalCorrespondences(true);
    pcl::PointCloud<pcl::PointXYZ> Final;

    Eigen::Matrix<double,4,4,Eigen::ColMajor> trans_init;
    trans_init.setIdentity();
    icp.align(Final, trans_init); // TODO: Catch ERROR and process it

    if(DEBUG) {
        std::cout << "Initial Transformation" << std::endl << trans_init << std::endl << std::endl;
        std::cout << "Final Transformation" << std::endl << icp.getFinalTransformation() << std::endl << std::endl;
    }

    return icp.getFitnessScore();
}

int main(){
    int N = 10;
    Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor> pro1, pro2;
    pro1 = pro1.Random(3, N);
    pro2 = pro1 * 3;
    std::cout << "Pro1" << std::endl << pro1 << std::endl << std::endl;
    std::cout << "Pro2" << std::endl << pro2 << std::endl << std::endl;
    std::cout << Icp_algorithm(pro2, pro1) << std::endl;
    std::cout << Icp_algorithm(pro1, pro2) << std::endl;
    return 0;
}
