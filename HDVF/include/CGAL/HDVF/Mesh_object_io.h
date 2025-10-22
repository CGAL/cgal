// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_HDVF_MESH_OBJECT_H
#define CGAL_HDVF_MESH_OBJECT_H

#include <CGAL/license/HDVF.h>

#include <CGAL/centroid.h>
#include <CGAL/bounding_box.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <array>

namespace CGAL {
namespace Homological_discrete_vector_field {

/** \brief Type of cells of Mesh_object_io.
 *
 * *Sorted* vector of the vertex indices.
 */
typedef std::vector<size_t> Io_cell_type ;
/** \brief Type of pre-chains in Mesh_object_io (list of cells without coefficients). */
typedef std::vector<Io_cell_type> Io_chain_type ;


// ----- VTK format -----

// For triangular meshes
// Associated to any dimension the type number of associated VTK cells
static std::vector<size_t> VTK_types_IO = {1, 3, 5, 10} ;


// ----- Generic -----

inline bool check_sanity_line(const std::string &line, const std::string &file)
{
    // Check that line is sanitized. If not, throw.
    for ( size_t i = 0; i < line.size(); ++ i ) {
        if ( ! ( std::isspace(line[i]) || std::isdigit(line[i]) || (line[i] == '-') || (line[i] == '.') || (line[i] == 'e')) ) {
            std::cerr << "Error:\n  Cannot parse file " << file << std::endl ;
            std::cerr << "--- " << line[i] << std::endl ;
            return false;
        }
    }
    return true ;
}

inline bool get_next_uncommented_line(std::ifstream &infile, std::string &result) {
    while(getline(infile,result)) {
        if(result.length() > 1 && result[0] != '#') {
            return true;
        }
    }
    return false;
}

/** \brief Load nodes from a .nodes file.
 *
 * Load vertices coordinates from a .nodes file.
 *
 * \param[in] nodes_file Name of the input file.
 * \param[in] nodes Pointer to a vector of points into which nodes are outputed.
 * \param[in] adapt If `fill` is false, nodes must have the same dimension as the traits Point, if true, nodes dimension can be lower (and missing coordinates are filled with zeros) or higher (and coordinates are truncated to the traits dimension).
 **/

template <typename Traits>
inline size_t read_nodes(const std::string &node_file, std::vector<typename Traits::Point> *nodes, bool adapt = false)
{
    typedef typename Traits::Point Point;
    std::ifstream in_file (node_file) ;
    if ( ! in_file . good () ) {
        std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << node_file << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    // First line is the number of nodes
    size_t nnodes, nnodes_tmp, nodes_dim, padding, d ;
    bool fill;
    if ( ! in_file.eof())
    {
        std::string line;
        getline( in_file, line );
        check_sanity_line(line, node_file) ;
        // First number is the number of nodes, then dimension of nodes
        std::istringstream is (line);
        is >> nnodes ;
        is >> nodes_dim ;
        if (!adapt) {
            if (nodes_dim != Traits::Dimension::value){
                std::cerr << "read_nodes error: dimension of nodes incompatible with Traits::Dimension" << std::endl;
                throw("read_nodes error: dimension of nodes incompatible with Traits::Dimension");
            }
        }
        else {
            if (nodes_dim < Traits::Dimension::value){
                fill = true;
                padding = Traits::Dimension::value - nodes_dim;
            }
            else
                fill = false;
        }
        d = (nodes_dim <= Traits::Dimension::value)?nodes_dim:Traits::Dimension::value;
    }
    nnodes_tmp =  nnodes ;
    while ( !(in_file.eof()) && (nnodes_tmp>0))
    {
        size_t trash ;
        double x ;
        std::vector<double> node ;
        --nnodes_tmp ;
        std::string line;
        getline( in_file, line );
        check_sanity_line(line, node_file) ;
        std::istringstream is (line);
        is >> trash ;
        for (int i = 0; i<d; ++i)
        {
            is >> x ;
            node.push_back(x) ;
        }
        
        if constexpr (Traits::Dimension::value == 2){
            nodes->push_back(Point(typename Traits::FT(node[0]), typename Traits::FT(node[1]))) ;
        }else if constexpr (Traits::Dimension::value == 3){
            if ((padding <= 0) || (!adapt))
                nodes->push_back(Point(typename Traits::FT(node[0]), typename Traits::FT(node[1]), typename Traits::FT(node[2]))) ;
            else
                nodes->push_back(Point(typename Traits::FT(node[0]), typename Traits::FT(node[1]), 0.)) ;
        }else{
            std::vector<typename Traits::FT> res(Traits::Dimension::value) ;
            for (int i=0; i<d; ++i)
                res[i] = typename Traits::FT(node[i]);
            if (adapt)
            {
                for (int i=d; i<Traits::Dimension::value; ++i)
                    res[i] = 0. ; ;
            }
            nodes->push_back(Point(res.begin(), res.end()));
        }
    }
    in_file.close() ;
    return nnodes;
}

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Mesh_object_io` is an intermediate IO class, used to load triangular/tetraedral meshes and produce simplicial complexes.

 \tparam Traits a geometric traits class model of the `HDVFTraits` concept.
 */

// Generic Mesh_object_io class - for 3D triangular meshes
template <typename Traits>
class Mesh_object_io
{
    public:
    typedef typename Traits::Point Point ;
    typedef typename Traits::Bbox Bbox ;
private:
    // Write vtk file
    template <typename CoefficientRing>
    void write_vtk(const std::string &filename, const std::vector<Point> &nodes, const std::vector<Io_chain_type> &chains, const std::vector<CoefficientRing> *labels=NULL, const std::string scalar_type="none")
    {
        bool with_scalars = (labels != NULL) ;
        // Load ...
        std::ofstream out ( filename, std::ios::out | std::ios::trunc);

        if ( !out . good () ) {
            std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        // Header
        out << "# vtk DataFile Version 2.0" << std::endl ;
        out << "generators" << std::endl ;
        out << "ASCII" << std::endl ;
        out << "DATASET  UNSTRUCTURED_GRID" << std::endl ;

        // Points
        size_t nnodes = nodes.size() ;
        out << "POINTS " << nnodes << " double" << std::endl ;
        for (Point n : nodes) {
            out << Traits::to_point3(n);
            out << std::endl ;
        }

        // Cells
        // Number of cells : for each chain, each number of vertices for each cell
        // Size : size of a cell defined by d vertices : d+1
        size_t ncells_tot = 0, size_cells_tot = 0 ;
        for (size_t i = 0; i<chains.size(); ++i)
        {
            // for each cells in the ith chain
            for (Io_cell_type c : chains.at(i))
            {
                ncells_tot += 1 ;
                size_cells_tot += (c.size()+1) ;
            }

        }
        out << "CELLS " << ncells_tot << " " << size_cells_tot << std::endl ;
        // Output cells
        std::vector<size_t> types ;
        std::vector<CoefficientRing> scalars ;
        for (size_t i = 0; i<chains.size(); ++i)
        {
            const Io_chain_type cc = chains.at(i) ;

            for (Io_cell_type c : cc)
            {
                const int d = c.size() ;
                const size_t cell_type = VTK_types_IO.at(d-1) ;

                out << d << " " ;
                for (size_t j : c)
                    out << j << " " ;
                out << std::endl ;
                types.push_back(cell_type) ;
                if (with_scalars)
                    scalars.push_back(labels->at(i)) ;
            }
        }

        // CELL_TYPES
        out << "CELL_TYPES " << ncells_tot << std::endl ;
        for (size_t t : types)
            out << t << " " ;
        out << std::endl ;

        if (with_scalars)
        {
            // CELL_TYPES
            out << "CELL_DATA " << ncells_tot << std::endl ;
            out << "SCALARS CriticalCellsId " << scalar_type << " 1" << std::endl ;
            out << "LOOKUP_TABLE default" << std::endl ;
            for (CoefficientRing s : scalars)
                out << s << " " ;
            out << std::endl ;
        }
        out.close() ;
    }

public:
    // The variable `dim` is used to encode both the dimension of the object loaded and wether it encodes a complex (with cells of various dimensions) or a mesh (a collection of triangles)
    // - if dim > 0 : Mesh_object_io encodes a mesh and all cells have dimension d
    // - if dim < 0 : Mesh_object_io encodes a complex (possibly incomplete) of dimension d
    int dim = 0 ;
    size_t nvertices, ncells, nedges ;
    std::vector<Point> nodes ; // Coordinates of vertices (optional)
    std::vector<Io_cell_type> cells ;

    /* \brief Default constructor.
     *
     * Create an empty Mesh_object_io.
     */
    Mesh_object_io(int d = 0) : dim(d), nvertices(0), ncells(0), nedges(0) {}

    /** \brief Constructor from a vector of Point (vertex coordinates) and a vector of simplices.
     *
     * Simplices are described by the list of vertex indices.
     *
     * \param[in] d The dimension `d` can be positive or negative:
     * - if positive: the set of simplicial cells loaded is a "mesh" and all cells have the same dimension
     * - if negative: the set of simplicial cells loaded have various dimensions and `d` must be the maximum of these dimensions.
     * \param[in] vnodes Vector of vertex coordinates.
     * \param[in] vcells Vector of cells (described by a sorted vector of indices)
     * \param[in] sort_data If `true` the vectors of vertex indices are sorted, if `false` they are assumed to be sorted (faster).
     */
    Mesh_object_io(int d, const std::vector<Point> &vnodes, const std::vector<Io_cell_type> &vcells, bool sort_data = false) : dim(d), nvertices(vnodes.size()), ncells(vcells.size()), nedges(0), nodes(vnodes), cells(vcells) {
        check_dimension() ;
        if (sort_data) {
            // Sort the vector encoding each cell
            for (int i=0; i<ncells; ++i) {
                std::sort(cells.at(i).begin(), cells.at(i).end());
            }
        }
    }


    Mesh_object_io(const Point &BBmin, const Point &BBmax)
    : dim(3), nvertices(8), ncells(12), nedges(0)
{
    std::vector points = {BBmin, BBmax} ;
    Bbox ic = bounding_box(points.begin(), points.end()) ;
    nodes.resize(8) ;
    nodes[0] = ic[0] ;
    nodes[1] = ic[1] ;
    nodes[2] = ic[2] ;
    nodes[3] = ic[3] ;
    nodes[4] = ic[5] ;
    nodes[5] = ic[6] ;
    nodes[6] = ic[7] ;
    nodes[7] = ic[4] ;

    cells.resize(12) ;
    cells[0] = Io_cell_type({0, 1, 4}) ;
    cells[1] = Io_cell_type({1, 4, 5}) ;
    cells[2] = Io_cell_type({1, 2, 6}) ;
    cells[3] = Io_cell_type({1, 5, 6}) ;
    cells[4] = Io_cell_type({0, 1, 3}) ;
    cells[5] = Io_cell_type({1, 2, 3}) ;
    cells[6] = Io_cell_type({2, 3, 6}) ;
    cells[7] = Io_cell_type({3, 6, 7}) ;
    cells[8] = Io_cell_type({0, 3, 4}) ;
    cells[9] = Io_cell_type({3, 4, 7}) ;
    cells[10] = Io_cell_type({4, 5, 6}) ;
    cells[11] = Io_cell_type({4, 6, 7}) ;
}

// TODO : check that de dimensions of nodes are consistent...


    std::vector<Point> get_nodes () const
    {
        return nodes ;
    }

    // Mesh operations
    // Add a Mesh_object_io to the current Mesh_object_io
    void push_back(const Mesh_object_io &mesh)
    {
        size_t off = nvertices ; // The index of all the cells of mesh has to be incremented by off
        nvertices += mesh.nvertices ;
        ncells += mesh.ncells ;
        // Append all the vertices
        for (size_t i=0; i<mesh.nodes.size(); ++i)
            nodes.push_back(mesh.nodes.at(i)) ;
        // Append all the cells (and increment their indices by off)
        for (size_t i=0; i<mesh.cells.size(); ++i)
        {
            Io_cell_type tmp ;
            for (size_t c : mesh.cells.at(i))
                tmp.push_back(c+off) ;
            cells.push_back(tmp) ;
        }
    }

    void add_node(const Point &v) {nodes.push_back(v); ++nvertices ;}

    void clear_cells() { cells.clear() ; ncells = 0 ; }

    void clear_nodes() { nodes.clear() ; nvertices = 0 ; }

    void clear() { clear_nodes() ; clear_cells() ;}

    void add_cell(Io_cell_type &c, bool sort_indices = false) {
        if (sort_indices)
            std::sort(c.begin(), c.end());
        cells.push_back(c); ++ncells ;
    }

    size_t cells_of_dim (int q) const
    {
        size_t n = 0 ;
        for (Io_cell_type c : cells)
        {
            if (c.size() == (q+1))
                ++n ;
        }
        return n ;
    }
    // OFF
    bool read_off(const std::string &filename)
    {
        dim = 0 ;
        // 0 - open input file
        std::ifstream infile(filename);
        if(!infile.is_open()) {
            // failed to open the file
            std::cerr << "Warning: file " << filename << " does not exist" << std::endl ;
            return false;
        }
        // 1 - header
        std::string header;
        if (!get_next_uncommented_line(infile, header)) {
            return false;
        }
        // todo : check for header == "off"

        // 2 - number of vertices, number of faces, number of edges (can be ignored)
        std::string info;
        if (!get_next_uncommented_line(infile, info)) {
            return false;
        }
        std::istringstream info_stream;
        info_stream.str(info);

        info_stream >> nvertices >> ncells >> nedges;

        nodes.resize(nvertices) ;
        for(auto i=0; i < nvertices; ++i) {
            if (!get_next_uncommented_line(infile,info)) {
                return false;
            }
            std::istringstream info_stream(info);
            std::array<double,3> p;
            info_stream >> p[0] >> p[1] >> p[2] ;
            nodes[i] = Point(p[0], p[1], p[2]) ;
        }

        // 4 - the actual faces
        cells.resize(ncells);
        for(auto i=0; i < ncells; ++i) {
            if (!get_next_uncommented_line(infile,info)) {
                return false;
            }
            std::istringstream info_stream(info);
            unsigned long n;
            unsigned long index;
            info_stream >> n;
            Io_cell_type c ;
            // Read vertex indices
            for (auto j = 0; j < n; ++j) {
                info_stream >> index;
                c.push_back(index) ;
            }
            // Sort the vector
            std::sort(c.begin(), c.end());
            // Insert the cell
            cells[i] = c ;
            dim = (c.size()-1>dim)?c.size()-1:dim ;
        }

        infile.close();
        return true;
    }

    bool write_off(const std::string &filename)
    {
        // 0 - open input file
        std::ofstream outfile(filename);
        if(!outfile.is_open()) {
            std::cerr << "Warning: cannot open file " << filename << std::endl ;
            // failed to open the file
            return false;
        }
        // 1 - header
        outfile << "OFF" << std::endl ;
        outfile << nvertices << " " << cells_of_dim(2) << " " << nedges << std::endl ;
        // 2 - nodes
        for (size_t i=0; i<nvertices; ++i)
        {
            outfile << nodes.at(i)[0] << " " << nodes.at(i)[1] << " " << nodes.at(i)[2] << std::endl ;
        }
        // 3 - cells (export only triangles)
        for (size_t i=0; i<ncells; ++i)
        {
            const size_t ni = cells.at(i).size() ;
            if (ni == 3)
            {
                outfile << ni << " " ;
                for (size_t k : cells.at(i))
                    outfile << k << " " ;
                outfile << std::endl ;
            }
        }
        outfile.close() ;
        return true;
    }

    // VTK
    void write_to_vtk(const std::string &filename)
    {
        std::vector<Io_chain_type> chains{cells} ;
        write_vtk<int>(filename, nodes, chains) ;
    }

    // SIMP
    bool write_simp(const std::string &filename)
    {
        // 0 - open input file
        std::ofstream outfile(filename);
        if(!outfile.is_open()) {
            std::cerr << "Warning: cannot open file " << filename << std::endl ;
            // failed to open the file
            return false;
        }
        // 1 - write cells
        for (size_t i=0; i<ncells; ++i)
        {
            const Io_cell_type cell = cells.at(i) ;
            for (size_t c : cell)
                outfile << c << " " ;
            outfile << std::endl ;
        }
        outfile.close() ;
        return true;
    }

    bool read_simp(const std::string &filename)
    {
        int d = 0 ;
        // 0 - open input file
        std::ifstream infile(filename);
        if(!infile.is_open()) {
            // failed to open the file
            std::cerr << "Warning: file " << filename << " does not exist" << std::endl ;
            return false;
        }
        std::size_t line_number = 0;
        while ( !(infile.eof()) )
        {
            ++line_number;
            std::string line;
            getline( infile, line );
            // Check that line is sanitized. If not, throw.
            for ( size_t i = 0; i < line.size(); ++ i ) {
                if ( ! ( std::isspace(line[i]) || std::isdigit(line[i]) ) ) {
                    std::cerr << "Fatal Error:\n  Cannot parse line #" << line_number << " of " << filename << "\n";
                    std::cerr << " --> " << line << "\n";
                    throw std::runtime_error("File Parsing Error: Invalid file");
                }
            }
            Io_cell_type cell ;
            std::istringstream is( line );
            size_t v;
            // Read vertex indices
            while ( is >> v )
                cell.push_back(v);
            // Sort the vector of indices
            std::sort(cell.begin(), cell.end());
            // Add this simplex to cells
            if (!(cell.empty()))
            {
                ++ncells ;
                cells.push_back(cell) ;
                const int dcell = cell.size()-1 ;
                if (dcell > d)
                    d = dcell ;
            }
        }
        dim = d ;
        infile.close() ;
        return true ;
    }

    bool read_nodes_file(const std::string &filename) ;

    void print_infos () const
    {
        std::cout << "Mesh_object_io infos - dim : "<< dim << ", nodes : " << nodes.size() << ", cells : " << cells.size() << std::endl ;
        for (int q = 0; q <= dim; ++q)
            std::cout << "cells of dim " << q << " : " << cells_of_dim(q) << std::endl ;
    }

    // Mesh computations
    Point centroid()
    {
        return CGAL::centroid(nodes.begin(), nodes.end()) ;
    }

    double radius(const Point &bary)
    {
        double r = 0 ;
        for (Point v : nodes)
        {
            r = (std::max) (r, sqrt(squared_distance(v,bary))) ;
        }
        return r ;
    }

    Bbox bbox(double ratio=1.)
    {
        return bounding_box(nodes.begin(), nodes.end()) .bbox();
        // @todo deal with ratio
    }
private:
    void check_dimension()
    {
        if (dim > 0) // exact dimension
        {
            for (Io_cell_type cell : cells)
            {
                if (cell.size() != (dim+1))
                    throw "Mesh has a cell of inconsistent dimension" ;
            }
        }
        else // max dim
        {
            for (Io_cell_type cell : cells)
            {
                if ((cell.size() > (-dim+1)))
                    throw "Mesh has a cell of inconsistent dimension" ;
            }
        }
    }
} ;




} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */


#endif // CGAL_HDVF_MESH_OBJECT_H
