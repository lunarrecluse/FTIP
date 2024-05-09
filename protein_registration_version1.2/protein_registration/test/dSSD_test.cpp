//
// Created by Xin Sui on 1/29/18.
//

#include "../src/dSSD_sort.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

void dSSD_test(const int permu, const int length){
    static Eigen::Matrix<int, 12, 4, Eigen::RowMajor> Permu;
    Permu << 0,1,2,3, 0,2,3,1, 0,3,1,2, 1,0,3,2, 1,2,0,3, 1,3,2,0,
            2,0,1,3, 2,1,3,0, 2,3,0,1, 3,0,2,1, 3,1,0,2, 3,2,1,0;

    Eigen::Matrix<int, 6, 12> edge_matrix;
    edge_matrix << 0 ,1 ,2 ,0 ,3 ,4 ,1 ,3 ,5 ,2 ,4 ,5,
                   1 ,2 ,0 ,4 ,0 ,3 ,3 ,5 ,1 ,5 ,2 ,4,
                   2 ,0 ,1 ,3 ,4 ,0 ,5 ,1 ,3 ,4 ,5 ,2,
                   3 ,5 ,4 ,2 ,1 ,5 ,0 ,4 ,2 ,1 ,0 ,3,
                   4 ,3 ,5 ,1 ,5 ,2 ,2 ,0 ,4 ,0 ,3 ,1,
                   5 ,4 ,3 ,5 ,2 ,1 ,4 ,2 ,0 ,3 ,1 ,0;

    Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::ColMajor> rand1;
    Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::ColMajor> rand2;
    rand2 = rand2.Random(6, length);
    rand1.resize(6, length);
    for(int j = 0; j != 6; ++j){
        rand1.row(j) = rand2.row( edge_matrix(j, permu) );
    }
    std::cout << "Protein A is " << std::endl << rand1 << std::endl << std::endl;
    std::cout << "which is derived from permutation" << Permu.row(permu) << "of Protein B"  << std::endl << rand2 << std::endl << std::endl;
    std::vector<dSSD_vec> ret = dSSD_sort(rand1, rand2, 2 * length);
    std::cout << "After calling dSSD_sort, the top " << 2 * length << "elements are" << std::endl;
    std::cout << "dSSD" << '\t' << "index_a" << '\t' << "index_b" << '\t' << "permu" << std::endl;
    for(auto iter = ret.begin(); iter != ret.end(); ++iter){
        std::cout << iter->dSSD << '\t' << iter->index_a << '\t' << iter->index_b << '\t' << iter->permu << std::endl;
    }
}


int main(){
    for (int i=0; i!=12; ++i){
        dSSD_test(i, 10);
    }
    return 0;
}