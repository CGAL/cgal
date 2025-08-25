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

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>

namespace CGAL {
namespace HDVF {

/** \brief Type of cells of Mesh_object_io.
 *
 * *Sorted* vector of the vertex indices.
 */
typedef std::vector<size_t> Io_cell_type ;
/** \brief Type of pre-chains in Mesh_object_io (list of cells without coefficients). */
typedef std::vector<Io_cell_type> Io_chain_type ;

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Io_node_type` implements a simple data type used to import vertices coordinates (nodes) in various dimensions.
 Hence, coordinates are loaded as vectors of `double`.

 The class provides standard affine geometry functions on such points.
 */

// Type of points (for vertices coordinates in R^d)
struct Io_node_type {
private:
    std::vector<double> _coords;

public:
    Io_node_type(size_t d = 0, double x = 0.) : _coords(d,x) {}
    Io_node_type(std::vector<double> v) : _coords(v) {}
    Io_node_type(const Io_node_type& v) : _coords(v._coords) {}

    size_t size() const { return _coords.size(); }
    double at(size_t i) const { return _coords.at(i) ;}
    double& operator[](size_t i) { return _coords.at(i) ;}
    Io_node_type& operator=(const Io_node_type& v) { _coords = v._coords ; return *this ; }
    void push_back(double x) { _coords.push_back(x) ; }
    std::vector<double> get_coords() const { return _coords; }

    friend Io_node_type & operator+(const Io_node_type &v1, const Io_node_type &v2)
    {
        assert(v1.size() == v2.size()) ;
        Io_node_type &tmp = *new(Io_node_type)(v1.size(),0) ;
        for (size_t i =0; i<v1.size(); ++i)
            tmp[i] = (v1.at(i) + v2.at(i)) ;
        return tmp ;
    }

    friend Io_node_type & operator-(const Io_node_type &v1, const Io_node_type &v2)
    {
        assert(v1.size() == v2.size()) ;
        Io_node_type &tmp = *new(Io_node_type)(v1.size(),0) ;
        for (size_t i =0; i<v1.size(); ++i)
            tmp[i] = (v1.at(i) - v2.at(i)) ;
        return tmp ;
    }

    friend Io_node_type & operator+= (Io_node_type &v1, const Io_node_type &v2)
    {
        assert(v1.size() == v2.size()) ;
        for (size_t i =0; i<v1.size(); ++i)
            v1[i] += v2.at(i) ;
        return v1 ;
    }

    friend Io_node_type operator/(const Io_node_type &v1, double d)
    {
        Io_node_type &tmp = *new(Io_node_type)(v1.size(),0) ;
        for (size_t i =0; i<v1.size(); ++i)
            tmp[i] = v1.at(i)/d ;
        return tmp ;
    }

    friend Io_node_type & operator/=(Io_node_type &v1, double d)
    {
        for (size_t i =0; i<v1.size(); ++i)
            v1[i] /= d ;
        return v1 ;
    }

    friend Io_node_type & operator*=(Io_node_type &v, double d)
    {
        for (size_t i =0; i<v.size(); ++i)
            v[i] *= d ;
        return v ;
    }

    friend Io_node_type & operator*=(Io_node_type &v, Io_node_type &lambda)
    {
        for (size_t i =0; i<v.size(); ++i)
            v[i] *= lambda.at(i) ;
        return v ;
    }

    friend double dist(const Io_node_type &v1, const Io_node_type &v2)
    {
        assert(v1.size() == v2.size()) ;
        double x=0., tmp ;
        for (size_t i =0; i<v1.size(); ++i)
        {
            tmp = v2.at(i) - v1.at(i) ;
            x += tmp*tmp ;
        }

        return sqrt(x) ;
    }

    friend Io_node_type max(const Io_node_type &v1, const Io_node_type &v2)
    {
        assert(v1.size() == v2.size()) ;
        Io_node_type tmp(v1) ;
        for (size_t i=0; i<v1.size(); ++i)
            tmp[i] = std::max(tmp.at(i), v2.at(i)) ;
        return tmp ;
    }

    friend Io_node_type min(const Io_node_type &v1, const Io_node_type &v2)
    {
        assert(v1.size() == v2.size()) ;
        Io_node_type tmp(v1) ;
        for (size_t i=0; i<v1.size(); ++i)
            tmp[i] = std::min(tmp.at(i), v2.at(i)) ;
        return tmp ;
    }

    friend void normalize(Io_node_type &v)
    {
        const Io_node_type zero = Io_node_type(v.size(),0) ;
        double n = dist(zero, v) ;
        for (size_t i =0; i<v.size(); ++i)
            v[i] /= n ;
    }
};




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

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Mesh_object_io` is an intermediate IO class, used to load triangular/tetraedral meshes and produce simplicial complexes.

 */

// Generic Mesh_object_io class - for 3D triangular meshes
class Mesh_object_io
{
private:
    // Write vtk file
    template <typename CoefficientRing>
    void write_vtk(const std::string &filename, const std::vector<Io_node_type> &nodes, const std::vector<Io_chain_type> &chains, const std::vector<CoefficientRing> *labels=NULL, const std::string scalar_type="none")
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
        for (Io_node_type n : nodes)
            out << n.at(0) << " " << n.at(1) << " " << n.at(2) << std::endl ;

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
    std::vector<Io_node_type> nodes ; // Coordinates of vertices (optional)
    std::vector<Io_cell_type> cells ;

    /* \brief Default constructor.
     *
     * Create an empty Mesh_object_io.
     */
    Mesh_object_io() : dim(0), nvertices(0), ncells(0), nedges(0) {}

    /** \brief Constructor from a vector of Io_node_type (vertices coordinates) and a vector of simplices.
     *
     * Simplices are described by the list of vertices indices.
     *
     * \param[in] d The dimension `d` can be positive or negative:
     * - if positive: the set of simplicial cells loaded is a "mesh" and all cells have the same dimension
     * - if negative: the set of simplicial cells loaded have various dimensions and `d` must be the maximum of these dimensions.
     * \param[in] vnodes Vector of vertices coordinates.
     * \param[in] vcells Vector of cells (described by a sorted vector of indices)
     * \param[in] sort_data If `true` the vectors of vertex indices are sorted, if `false` they are assumed to be sorted (faster).
     *
     */
    Mesh_object_io(int d, const std::vector<Io_node_type> &vnodes, const std::vector<Io_cell_type> &vcells, bool sort_data = false) : dim(d), nvertices(vnodes.size()), ncells(vcells.size()), nedges(0), nodes(vnodes), cells(vcells) {
        check_dimension() ;
        if (sort_data) {
            // Sort the vector encoding each cell
            for (int i=0; i<ncells; ++i) {
                std::sort(cells.at(i).begin(), cells.at(i).end());
            }
        }
    }

    /* \brief Copy constructor. */
    Mesh_object_io(const Mesh_object_io &m) : dim(m.dim), nvertices(m.nvertices), ncells(m.ncells), nedges(m.nedges), nodes(m.nodes), cells(m.cells) {}

    std::vector<std::vector<double> > get_nodes ()
    {
        std::vector<std::vector<double> > res ;
        for (Io_node_type v : nodes)
            res.push_back(v.get_coords()) ;
        return res ;
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

    void add_node(const Io_node_type &v) {nodes.push_back(v); ++nvertices ;}

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
            Io_node_type p(3) ;
            info_stream >> p[0] >> p[1] >> p[2] ;
            nodes[i] = p ;
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
            // Read vertices indices
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
            outfile << nodes.at(i).at(0) << " " << nodes.at(i).at(1) << " " << nodes.at(i).at(2) << std::endl ;
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
            // Read vertices indices
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
    Io_node_type barycenter()
    {
        // Init the barycenter
        Io_node_type bary(3,0) ;
        // Compute the barycenter
        for (Io_node_type v : nodes)
        {
            bary += v ;
        }
        bary /= nodes.size() ;
        return bary;
    }

    double radius(const Io_node_type &bary)
    {
        double r = 0 ;
        for (Io_node_type v : nodes)
        {
            r = std::max (r, dist(v,bary)) ;
        }
        return r ;
    }
    std::pair<Io_node_type, Io_node_type> BB(double ratio=1.)
    {
        Io_node_type minBB(nodes.at(0)), maxBB(nodes.at(0)) ;
        for (size_t i=1; i<nodes.size(); ++i)
        {
            minBB = min(minBB, nodes.at(i)) ;
            maxBB = max(maxBB, nodes.at(i)) ;
        }
        Io_node_type c = (minBB+maxBB)/2. ;
        Io_node_type rad = (maxBB-minBB)/2. ;
        rad *= ratio ;
        return std::pair<Io_node_type, Io_node_type>(c-rad, c+rad) ;
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

inline Mesh_object_io mesh_BB(const Io_node_type &BBmin, const Io_node_type &BBmax)
{
    Mesh_object_io m ;
    Io_node_type delta = BBmax-BBmin ;
    m.nvertices = 8 ;
    m.ncells = 12 ;
    m.nodes.resize(8) ;
    m.nodes[0] = Io_node_type({0, 0, 0}) ;
    m.nodes[1] = Io_node_type({1, 0, 0}) ;
    m.nodes[2] = Io_node_type({1, 1, 0}) ;
    m.nodes[3] = Io_node_type({0, 1, 0}) ;
    m.nodes[4] = Io_node_type({0, 0, 1}) ;
    m.nodes[5] = Io_node_type({1, 0, 1}) ;
    m.nodes[6] = Io_node_type({1, 1, 1}) ;
    m.nodes[7] = Io_node_type({0, 1, 1}) ;
    for (size_t i=0; i<8; ++i)
    {
        m.nodes[i] *= delta ;
        m.nodes[i] += BBmin ;
    }

    m.cells.resize(12) ;
    m.cells[0] = Io_cell_type({0, 1, 4}) ;
    m.cells[1] = Io_cell_type({1, 4, 5}) ;
    m.cells[2] = Io_cell_type({1, 2, 6}) ;
    m.cells[3] = Io_cell_type({1, 5, 6}) ;
    m.cells[4] = Io_cell_type({0, 1, 3}) ;
    m.cells[5] = Io_cell_type({1, 2, 3}) ;
    m.cells[6] = Io_cell_type({2, 3, 6}) ;
    m.cells[7] = Io_cell_type({3, 6, 7}) ;
    m.cells[8] = Io_cell_type({0, 3, 4}) ;
    m.cells[9] = Io_cell_type({3, 4, 7}) ;
    m.cells[10] = Io_cell_type({4, 5, 6}) ;
    m.cells[11] = Io_cell_type({4, 6, 7}) ;
    return m;
}

// TODO : check that de dimensions of nodes are consistent...


} /* end namespace HDVF */
} /* end namespace CGAL */


#endif // CGAL_HDVF_MESH_OBJECT_H
