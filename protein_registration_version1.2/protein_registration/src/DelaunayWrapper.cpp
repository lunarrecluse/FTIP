//
// Created by Xin Sui on 1/29/18.
//

#include "DelaunayWrapper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <tuple>
#include <ios>
#include <string>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, Traits> Vertex_base;
typedef CGAL::Delaunay_triangulation_cell_base_3<Traits> Cell_base;
typedef CGAL::Triangulation_data_structure_3<Vertex_base, Cell_base> TDS;
typedef CGAL::Fast_location LP;
typedef CGAL::Delaunay_triangulation_3<Traits, TDS, LP> Delaunay;
typedef Delaunay::Finite_cells_iterator Finite_cells_iterator;
typedef Delaunay::Finite_edges_iterator Finite_edges_iterator;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Point Point;
// typedef Delaunay::Vector Vector; // Point - Point = Diff



void DelaunayWrapper::read_asa (const boost::filesystem::path& file_name, unsigned length) //"surface/1a0ia-surf.asa"
{

    Delaunay_data.Vertex_matrix.resize(3,length);
    boost::filesystem::ifstream in;
    in.open(file_name, std::ios::in);

    std::string line;
    unsigned j = 0;
    while(getline(in, line)){
        Delaunay_data.Vertex_matrix(0,j) = std::stod(line.substr(30,8));
        Delaunay_data.Vertex_matrix(1,j) = std::stod(line.substr(38,8));
        Delaunay_data.Vertex_matrix(2,j) = std::stod(line.substr(46,8));
        j++;
    }
    if(j!=length)
        std::cout << file_name << " has " << j << " elemenets!" << std::endl;
    in.close();
}

void DelaunayWrapper::read_pdb (const boost::filesystem::path& file_name, unsigned length)
{

    Delaunay_data.Vertex_matrix.resize(3,length);
    boost::filesystem::ifstream in;
    in.open(file_name, std::ios::in);

    std::string line;
    unsigned j = 0;
    while(getline(in, line)){
        Delaunay_data.Vertex_matrix(0,j) = std::stod(line.substr(30,8));
        Delaunay_data.Vertex_matrix(1,j) = std::stod(line.substr(38,8));
        Delaunay_data.Vertex_matrix(2,j) = std::stod(line.substr(46,8));
        j++;
    }
    if(j!=length)
        std::cout << file_name << " has " << j << " elemenets!" << std::endl;
    in.close();
}

DelaunayWrapper::DelaunayWrapper(const boost::filesystem::path& file_name, unsigned length){
    //TODO: find a way to remove the length variable, that is, support data with different length.
    if(!boost::filesystem::exists(file_name)){
        std::cout << "File " << file_name << " does not exist." << std::endl;
        //TODO: handle this error
    } else if (file_name.extension() == ".asa"){
        read_asa(file_name, length);
    } else if (file_name.extension() == ".pdb"){
        read_pdb(file_name, length);
    } else {
        std::cout << file_name << "is not a asa/pdb file." << std::endl;
        //TODO: handle this error
    }
    Delaunay_Transformation();
}

void DelaunayWrapper::Delaunay_Transformation()
{
    long length = Delaunay_data.Vertex_matrix.cols();
    assert(length > 3);
    std::vector<std::pair<Point,unsigned> > P;
    P.reserve(length);
    for(int i=0; i != length; ++i){
        P.emplace_back( std::make_pair( Point(Delaunay_data.Vertex_matrix(0,i), Delaunay_data.Vertex_matrix(1,i),
                                              Delaunay_data.Vertex_matrix(2,i)), i) );
    }
    Delaunay Triangulation;
    Triangulation.insert(P.begin(), P.end()); // TODO: try to not copy here

    unsigned n_cells = Triangulation.number_of_finite_cells();
    Delaunay_data.Edge_matrix.resize(6, n_cells);
    Delaunay_data.Vertex_reference.resize(4, n_cells);

    unsigned j = 0;
    for (Finite_cells_iterator cell_iterator = Triangulation.finite_cells_begin();
    cell_iterator != Triangulation.finite_cells_end(); ++cell_iterator) {

        Vertex_handle v0 = cell_iterator->vertex(0);
        Vertex_handle v1 = cell_iterator->vertex(1);
        Vertex_handle v2 = cell_iterator->vertex(2);
        Vertex_handle v3 = cell_iterator->vertex(3);

        Delaunay_data.Edge_matrix.col(j) << std::sqrt(CGAL::squared_distance(v0->point(), v1->point())),
                              std::sqrt(CGAL::squared_distance(v0->point(), v2->point())),
                              std::sqrt(CGAL::squared_distance(v0->point(), v3->point())),
                              std::sqrt(CGAL::squared_distance(v1->point(), v2->point())),
                              std::sqrt(CGAL::squared_distance(v1->point(), v3->point())),
                              std::sqrt(CGAL::squared_distance(v2->point(), v3->point()));

        Delaunay_data.Vertex_reference.col(j) << v0->info(), v1->info(), v2->info(), v3->info();
        ++j;
    }
}

