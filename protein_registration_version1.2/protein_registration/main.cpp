
#include <valarray>
#include <random>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <valarray> // TODO: Not updated for long. Change it?
#include <chrono>
#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "src/dSSD_sort.h"
#include "src/SSD_sort.h"
#include "src/Icp_algorithm.h"
#include "src/DelaunayWrapper.h"
#include "src/ThreadPool.h"

const int N = 200;

int main()
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.,1.);
    std::uniform_int_distribution<int> rint(0,N-1);

    Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand1;
    Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand2;
    rand2 = rand2.Random(3,N);
    rand1 = rand1.Random(3,N);


    // Apply rotation
    Eigen::Matrix<double,3,3,Eigen::ColMajor> rand_rotation = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();
    Eigen::Matrix<double,3,1,Eigen::ColMajor> rand_translation(dist(gen), dist(gen), dist(gen));
    Eigen::Matrix<double,4,4,Eigen::ColMajor> expect, ret;
    expect << rand_rotation, rand_translation, 0., 0., 0., 1.;
 //   rand1 = rand_rotation * rand2;
 //   rand1.colwise() += rand_translation;

    // Apply transposition
 //   Eigen::Transpositions<Eigen::Dynamic> transposition;
 //   transposition.resize(N);
 //   for(int k=0;k!=N;++k)
  //      transposition[k] = rint(gen);
 //   rand2 = rand2 * transposition;
    auto start = std::chrono::steady_clock::now();

    for(int j=0;j!=3;++j) {
        DelaunayWrapper Pro1(rand1);
        DelaunayWrapper Pro2(rand2);

        std::vector<dSSD_vec> dSSD_result = dSSD_sort(Pro1.Delaunay_data.Edge_matrix, Pro2.Delaunay_data.Edge_matrix, 200);
        std::vector<Ini_Trans> SSD_result = SSD_sort(dSSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix,
                                                     Pro2.Delaunay_data.Vertex_reference, Pro2.Delaunay_data.Vertex_reference, 100);

        ICP_result result = Icp_algorithm(SSD_result, Pro1.Delaunay_data.Vertex_matrix, Pro2.Delaunay_data.Vertex_matrix, 0.5, 20, 0);

    }

    auto end = std::chrono::steady_clock::now();

    std::cout << "Total time in seconds is " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count() <<std::endl;


    return 0;

}













/*
void reg_algorithm(protein& pro1, protein& pro2, int dSSD_cutoff, int SSD_cutoff,
                  double max_corres_dist, int max_iter, bool DEBUG=0) // TODO: by the current def. they cannot be const
{

    pro1.Delaunay_Transformation();
    pro2.Delaunay_Transformation();

    pro1.build_table(DEBUG);
    pro2.build_table(DEBUG);

    // find dSSD_cutoff min dSSD val for all pairs.
    std::vector<tet_pair> dSSD_tet_tet_permu;
    dSSD_sort(pro1, pro2, dSSD_tet_tet_permu, dSSD_cutoff, DEBUG);


    // find SSD_cutoff min SSD val for all pairs and compute TraVec etc.
    std::vector<Ini_Trans> trans;// Iterate over dSSDs to get Init
    SSD_sort(trans, dSSD_tet_tet_permu, SSD_cutoff, DEBUG);

    return Icp_algorithm(trans, pro1, pro2, max_corres_dist, max_iter, DEBUG);

}
*/

/*
int parallel(int M, int N, int dSSD_cutoff=200, int SSD_cutoff=50, double max_corres_dist=0.5, int max_iter=10, bool DEBUG=0)
{

    auto start = std::chrono::steady_clock::now();
    std::mutex critical;
    std::vector<double> error;
    error.resize(M);
    ThreadPool::ParallelFor(0, M, [&] (int i) {

        srand(time(0));
        Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand1;
        Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand2;
        rand2 = rand2.Random(3,N);
        rand1 = rand1.Random(3,N);

        protein pro1(rand1);
        protein pro2(rand2);

        pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ, double> icp_ret;
        icp_ret = reg_algorithm(pro1, pro2, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter, DEBUG);
        double this_error = icp_ret.getFitnessScore();
        error[i] = this_error;
    });
    auto end = std::chrono::steady_clock::now();
    std::cout << "Model parameters: dSSD_cutoff=" << dSSD_cutoff << " SSD_cutoff=" << SSD_cutoff << std::endl;
    std::cout << "Max correspondence distance=" << max_corres_dist << " Max iterations=" << max_iter << std::endl;
    std::cout << "Total squared error is  " << std::accumulate(error.begin(), error.end(), 0.) << " on " << M << " samples of size " << N << std::endl;
    std::cout << "Total time in seconds is " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count() <<std::endl;
    return 0;
}
*/
/*
int rand_transform(int M, int N, int dSSD_cutoff=200, int SSD_cutoff=50, double max_corres_dist=0.5, int max_iter=10, bool DEBUG=0)
{ // M proteins pairs of size N

    std::cout << "performing random transformation on " << M << " samples of size " << N << std::endl;
    std::clock_t begin = clock();

    double error = 0.;
    double max_error =0.;

    Eigen::Matrix<double,4,4,Eigen::ColMajor> max_error_trans, max_error_ret;
    protein max_error_pro1, max_error_pro2;
    for(int j=0;j!=M;++j)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.,1.);
        std::uniform_int_distribution<int> rint(0,N-1);

        Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand1;
        Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand2;
        rand2 = rand2.Random(3,N);

        // Apply rotation
        Eigen::Matrix<double,3,3,Eigen::ColMajor> rand_rotation = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();
        Eigen::Matrix<double,3,1,Eigen::ColMajor> rand_translation(dist(gen), dist(gen), dist(gen));
        Eigen::Matrix<double,4,4,Eigen::ColMajor> expect, ret;
        expect << rand_rotation, rand_translation, 0., 0., 0., 1.;
        rand1 = rand_rotation * rand2;
        rand1.colwise() += rand_translation;

        // Apply transposition
        Eigen::Transpositions<Eigen::Dynamic> transposition;
        transposition.resize(N);
        for(int k=0;k!=N;++k)
            transposition[k] = rint(gen);
        rand2 = rand2 * transposition;


        protein pro1(rand1);
        protein pro2(rand2);

        if(DEBUG)
        {
            std::cout << "The rotation matrix R is (A = RB + t)" << std::endl;
            std::cout << rand_rotation << std::endl << std::endl;
            std::cout << "And the translation vector t is" << std::endl;
            std::cout << rand_translation << std::endl << std::endl;
            std::cout << "And the protein A/B are" << std::endl << pro1 << std::endl << std::endl << pro2 << std::endl << std::endl;
            std::cout << "Expected" << std::endl << expect << std::endl << std::endl; // 0: run

        }
        pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ, double> icp_ret;
        icp_ret = reg_algorithm(pro1, pro2, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter, DEBUG);
        ret = icp_ret.getFinalTransformation().cast<double>();
        double this_error = (ret - expect).norm();
        error += this_error;
        if(this_error > max_error)
        {
            max_error = this_error;
            max_error_pro1 = pro1;
            max_error_pro2 = pro2;
            max_error_trans = expect;
            max_error_ret = ret;
        }
    }
    std::cout << "Model parameters: dSSD_cutoff=" << dSSD_cutoff << " SSD_cutoff=" << SSD_cutoff << std::endl;
    std::cout << "Max correspondence distance=" << max_corres_dist << " Max iterations=" << max_iter << std::endl;
    std::cout << "Total Error on transformation matrix is " << error << std::endl;
    std::cout << "The max error " << max_error << " happened on protein pair" << std::endl << max_error_pro1 << std::endl << std::endl
              << max_error_pro2 << std::endl << std::endl;
    std::cout << "The expected transformation matrix is" << std::endl << max_error_trans << std::endl << std::endl;
    std::cout << "And returned transformation matrix is" << std::endl << max_error_ret << std::endl << std::endl;
    std::clock_t end = clock(); // NOt accurate when computing in parallel
    std::cout << "Calculation finished in " << (double)(end - begin)/ CLOCKS_PER_SEC << " seconds. " << std::endl;
    return 0;
}
*/
/*
int rand_pair(int M, int N, int dSSD_cutoff=200, int SSD_cutoff=50, double max_corres_dist=0.5, int max_iter=10, bool DEBUG=0)
{ // M proteins pairs of size N

    std::cout << "performing random alignment on " << M << " samples of size " << N << std::endl;

    std::clock_t begin = clock();

    double error = 0.;

    for(int j=0;j!=M;++j)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.,1.);

        Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand1;
        Eigen::Matrix<double,3,Eigen::Dynamic,Eigen::ColMajor> rand2;
        rand2 = rand2.Random(3,N);
        rand1 = rand1.Random(3,N);

        protein pro1(rand1);
        protein pro2(rand2);

        pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ, double> icp_ret;
        icp_ret = reg_algorithm(pro1, pro2, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter, DEBUG);
        double this_error = icp_ret.getFitnessScore();
         error += this_error;
    }

    std::cout << "Model parameters: dSSD_cutoff=" << dSSD_cutoff << " SSD_cutoff=" << SSD_cutoff << std::endl;
    std::cout << "Max correspondence distance=" << max_corres_dist << " Max iterations=" << max_iter << std::endl;
    std::cout << "Total Error on Protein dismatch? is  " << error << " on " << M << " samples of size " << N << std::endl;
    std::clock_t end = clock(); // NOt accurate when computing in parallel
    std::cout << "Calculation finished in " << (double)(end - begin)/ CLOCKS_PER_SEC << " seconds. " << std::endl;
    return 0;
}
*/
/*
int protein_pair(std::string protein1, std::string protein2, int dSSD_cutoff=200, int SSD_cutoff=50,
                 double max_corres_dist=0.5, int max_iter=10, bool DEBUG=0)
{
    std::cout << "performing alignment on " << protein1 << " and " << protein2 << " datasets " << std::endl;

    std::clock_t begin = clock();

    protein pro1(protein1);
    protein pro2(protein2);
    //TODO: Is tuple a good idea? Hard to compute distance here.
    ICP ret;
    ret = reg_algorithm(pro1, pro2, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter, DEBUG);

    std::cout << "Model parameters: dSSD_cutoff=" << dSSD_cutoff << " SSD_cutoff=" << SSD_cutoff << std::endl;
    std::cout << "Max correspondence distance=" << max_corres_dist << " Max iterations=" << max_iter << std::endl;
    std::cout << "Protein A:" << std::endl << pro1 << std::endl << std::endl;
    std::cout << "Protein B:" << std::endl << pro2 << std::endl << std::endl;
    std::cout << "Final transformation matrix:"  << std::endl << ret.getFinalTransformation() << std::endl << std::endl;
    std::cout << "Has the algorithm converged? " << ret.hasConverged() << " Final fitness score: " << ret.getFitnessScore() << std::endl;
    std::clock_t end = clock(); // NOt accurate when computing in parallel
    std::cout << "Calculation finished in " << (double)(end - begin)/ CLOCKS_PER_SEC << " seconds. " << std::endl;

    return 0;

}
*/

/*
int main()
{
    std::string mode, file_1, file_2;
    int M;
    int N; int dSSD_cutoff=200;
    int SSD_cutoff=50;
    double max_corres_dist=0.5;
    int max_iter=10;


    std::ifstream infile;
    infile.open("/Users/xinsui/CLionProjects/protein_registration/config.txt", std::ios::in);

    infile >> mode;
    if(mode == "file") {
        infile >> file_1 >> file_2 >> dSSD_cutoff >> SSD_cutoff >> max_corres_dist >> max_iter;
        protein_pair(file_1, file_2, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter);
    } else if(mode == "transform") {
        infile >> M >> N >> dSSD_cutoff >> SSD_cutoff >> max_corres_dist >> max_iter;
        rand_transform(M, N, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter);
    } else if(mode == "random") {
        infile >> M >> N >> dSSD_cutoff >> SSD_cutoff >> max_corres_dist >> max_iter;
        rand_pair(M, N, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter);
    } else if(mode == "parallel"){
        infile >> M >> N >> dSSD_cutoff >> SSD_cutoff >> max_corres_dist >> max_iter;
        parallel(M, N, dSSD_cutoff, SSD_cutoff, max_corres_dist, max_iter);
    } else {
        std::cout << " No such mode! " << std::endl;
        return 1;
    }
    infile.close();

    return 0;
}
*/