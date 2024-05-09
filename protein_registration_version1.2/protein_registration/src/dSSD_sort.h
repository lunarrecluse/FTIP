//
// Created by Xin Sui on 1/19/18.
//

#ifndef PROTEIN_REGISTRATION_DSSD_SORT_H
#define PROTEIN_REGISTRATION_DSSD_SORT_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include <numeric>

struct dSSD_vec{
    dSSD_vec(){};
    dSSD_vec(double dSSDval, unsigned index_aval, unsigned index_bval, unsigned per): dSSD(dSSDval), index_a(index_aval), index_b(index_bval), permu(per) {};
    double dSSD;
    unsigned short index_a;
    unsigned short index_b;
    unsigned short permu;
};

inline bool operator< (const dSSD_vec& lhs, const dSSD_vec& rhs) {return lhs.dSSD < rhs.dSSD;}


inline double my_kernel(const Eigen::MatrixXd::ConstColXpr &x, const Eigen::MatrixXd::ConstColXpr &y) {
    return (x-y).squaredNorm();
}

template<typename Kernel>
Eigen::MatrixXd apply_kernel(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y, Kernel kernel) {
    return Eigen::MatrixXd::NullaryExpr(x.cols(), y.cols(),
                                        [&x,&y,&kernel](int i,int j) { return kernel(x.col(i), y.col(j)); });
}

template <typename Derived>
std::vector<unsigned> paritial_sort_indexes(const Eigen::ArrayBase<Derived>& v, int cutoff) {

    // initialize original index locations
    std::vector<unsigned> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    std::partial_sort(idx.begin(), idx.begin() + cutoff, idx.end(), [&v](unsigned i1, unsigned i2) { return v(i1) < v(i2); });
    idx.resize(cutoff);

    return idx;
}





template <typename Derived, typename OtherDerived>
std::vector<dSSD_vec> dSSD_sort(const Eigen::MatrixBase<Derived>& edge_matrix_a,
                                const Eigen::MatrixBase<OtherDerived>& edge_matrix_b, unsigned short dSSD_cutoff) {


    Eigen::Matrix<int, 12, 6, Eigen::RowMajor> transposition;
    transposition << 0, 1, 2, 3, 4, 5,
            2, 0, 2, 5, 5, 5,
            2, 2, 2, 4, 5, 5,
            0, 4, 3, 3, 4, 5,
            3, 3, 4, 3, 5, 5,
            4, 3, 4, 5, 4, 5,
            1, 3, 5, 3, 5, 5,
            3, 5, 5, 4, 4, 5,
            5, 1, 3, 3, 4, 5,
            2, 5, 4, 5, 4, 5,
            4, 2, 5, 4, 4, 5,
            5, 4, 2, 3, 4, 5;

    std::vector<dSSD_vec> dSSD_container;

    if(dSSD_cutoff > edge_matrix_a.cols() * edge_matrix_b.cols())
        dSSD_cutoff = edge_matrix_a.cols() * edge_matrix_b.cols();

    dSSD_container.reserve(2*dSSD_cutoff);
    for(int per = 0; per != 12; ++per) {
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ret( apply_kernel(edge_matrix_a, Eigen::Transpositions<6>(transposition.row(per)) * edge_matrix_b, my_kernel) );
        for (auto index: paritial_sort_indexes(ret.array(), dSSD_cutoff)) {
            unsigned row = index / ret.cols();
            unsigned col = index - row * ret.cols();
            dSSD_container.emplace_back(dSSD_vec(ret.array()(index), row, col, per));
        }
        std::partial_sort(dSSD_container.begin(), dSSD_container.begin() + dSSD_cutoff, dSSD_container.end());
        dSSD_container.resize(dSSD_cutoff);
    }
    return dSSD_container;
}

#endif //PROTEIN_REGISTRATION_DSSD_SORT_H
