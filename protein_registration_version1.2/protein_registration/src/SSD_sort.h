//
// Created by Xin Sui on 1/19/18.
//

#ifndef PROTEIN_REGISTRATION_SSD_SORT_H
#define PROTEIN_REGISTRATION_SSD_SORT_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include "dSSD_sort.h"

template <class T> void push_back_w_threshold(T& container, const typename T::value_type& ele, typename T::size_type cutoff)
{
    static double threshold = std::numeric_limits<double>::infinity();

    if(container.empty())
    {
        container.reserve(5 * cutoff + 1);
        threshold = std::numeric_limits<double>::infinity();
    }

    if(ele.value() < threshold)
        container.push_back(ele);

    if(container.size() == 5*cutoff)
    {
        std::sort(container.begin(), container.end());
        typename T::iterator it;
        for(it=container.end();it!=container.begin()+cutoff;--it)
        {
            container.pop_back();
        }
        threshold = container[cutoff].value();
    }

} // TODO::When cutoff is really big, and the total elements can come in, there is some problem. One more reason to write a class


class Ini_Trans{ // Not a good idea. Should define a base cell class with permu op. Together with ops used in Eigen::
public:
    typedef Eigen::Matrix<double, 3, 3, Eigen::ColMajor> RotMatType;
    typedef Eigen::Matrix<double, 3, 1, Eigen::ColMajor> TraVecType;
    Ini_Trans() {};
    Ini_Trans(const double this_SSD, const RotMatType& this_RotMat, const TraVecType& this_TraVec): SSD(this_SSD),
                                                                                                   RotMat(this_RotMat),
                                                                                                   TraVec(this_TraVec) {};
    double SSD;
    RotMatType RotMat;
    TraVecType TraVec;
    double value() const {return SSD;}
}; // Bad name. The name is SSD_Rot_Tra_pair.

inline bool operator< (const Ini_Trans& a, const Ini_Trans& b)
{
    return a.SSD < b.SSD;
}


template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
std::vector<Ini_Trans> SSD_sort(const std::vector<dSSD_vec>& dSSD_container, const Eigen::MatrixBase<Derived1>& vertex_matrix_a,
                                const Eigen::MatrixBase<Derived2>& vertex_matrix_b, const Eigen::MatrixBase<Derived3>& vertex_reference_a,
                                const Eigen::MatrixBase<Derived4>& vertex_reference_b, unsigned short SSD_cutoff){

    static Eigen::Matrix<int, 12, 4, Eigen::RowMajor> Permu;
    Permu << 0,1,2,3, 0,2,3,1, 0,3,1,2, 1,0,3,2, 1,2,0,3, 1,3,2,0,
            2,0,1,3, 2,1,3,0, 2,3,0,1, 3,0,2,1, 3,1,0,2, 3,2,1,0;

    if(SSD_cutoff > dSSD_container.size())
        SSD_cutoff = dSSD_container.size();

    std::vector<Ini_Trans> trans;
    trans.reserve(SSD_cutoff);

    for(std::vector<dSSD_vec>::const_iterator iter = dSSD_container.begin(); iter != dSSD_container.end(); ++iter){
        Eigen::Matrix<double, 3, 4, Eigen::ColMajor> tet_a;
        Eigen::Matrix<double, 3, 4, Eigen::ColMajor> tet_b;
        for(int i = 0; i != 4; ++i){
            tet_a.col(i) = vertex_matrix_a.col(vertex_reference_a(i,iter->index_a));
            tet_b.col(i) = vertex_matrix_b.col(vertex_reference_b(Permu(iter->permu, i), iter->index_b));
        }

        Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean_a = tet_a.rowwise().mean();
        Eigen::Matrix<double, 3, 1, Eigen::ColMajor> mean_b = tet_b.rowwise().mean();
        tet_a.colwise() -= mean_a;
        tet_b.colwise() -= mean_b; // You don't copy. Should be a way there.
        // TODO: write unary expression for extracting, and binary expression for before svd.
        Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > svd(tet_a * tet_b.transpose(),
                                                                           Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix<double, 3, 3> RotMat = svd.matrixU() * svd.matrixV().transpose();
        Eigen::Matrix<double, 3, 1> TraVec = mean_a - RotMat * mean_b;
        double SSD = (tet_a - RotMat * tet_b).squaredNorm();

        push_back_w_threshold(trans, Ini_Trans(SSD, RotMat, TraVec), SSD_cutoff);
    }
    
    if(trans.begin() + SSD_cutoff < trans.end()) {
        std::partial_sort(trans.begin(), trans.begin() + SSD_cutoff, trans.end());
        trans.resize(SSD_cutoff);
    }

    return trans;

}


#endif //PROTEIN_REGISTRATION_SSD_SORT_H
