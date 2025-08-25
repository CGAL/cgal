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

#ifndef CGAL_HDVF_TET_OBJECT_H
#define CGAL_HDVF_TET_OBJECT_H

#include <CGAL/license/HDVF.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <CGAL/HDVF/Mesh_object_io.h>

namespace CGAL {
namespace HDVF {


// Tetgen related

inline size_t read_nodes(const std::string &node_file, bool load_nodes, std::vector<Io_node_type> *nodes)
{
    std::ifstream in_file (node_file) ;
    if ( ! in_file . good () ) {
        std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << node_file << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    // First line is the number of nodes
    size_t nnodes, nnodes_tmp ;
    if ( ! in_file.eof())
    {
        std::string line;
        getline( in_file, line );
        check_sanity_line(line, node_file) ;
        // First number is the number of nodes
        std::istringstream is (line);
        is >> nnodes ;
    }
    nnodes_tmp =  nnodes ;
    while ( !(in_file.eof()) && (nnodes_tmp>0))
    {
        size_t trash ;
        double x ;
        Io_node_type node ;
        --nnodes_tmp ;
        std::string line;
        getline( in_file, line );
        check_sanity_line(line, node_file) ;
        std::istringstream is (line);
        is >> trash ;
        for (int i = 0; i<3; ++i)
        {
            is >> x ;
            node.push_back(x) ;
        }
        if (load_nodes)
            nodes->push_back(node) ;
    }
    in_file.close() ;
    return nnodes;
}

inline bool Mesh_object_io::read_nodes_file(const std::string &filename)
{
    nvertices = read_nodes(filename, true, &nodes) ;
}

// Tetgen

/* Class used to load tetgen outputs (for Alexander duality). */
class Tet_object_io : public Mesh_object_io
{
public:
    Tet_object_io(const std::string & prefix) : Mesh_object_io(), _prefix(prefix)
    {
        dim = -3 ;
        add_nodes() ;
        create_nodes() ;
        // If the user wished to keep tetgen indices, insert edges and faces below. Otherwise, they will be created by the Abstract_simplicial_chain_complex constructor
//        add_edges() ;
//        add_faces() ;
        add_tets() ;
        ncells = cells.size() ;
    }

    void add_nodes()
    {
        const std::string file_node = fnodes_from_prefix(_prefix) ;
        std::cout << "file_node : " << file_node << std::endl ;
        read_nodes_file(file_node) ;
    }

    void create_nodes()
    {
        for (size_t i=0; i<nvertices; ++i)
        {
            Io_cell_type cell({i}) ;
            cells.push_back(cell) ;
        }
        std::cout << "--- " << nvertices << "vert" << std::endl ;
    }

    void add_edges()
    {
        const std::string file_edge = fedges_from_prefix(_prefix) ;
        std::ifstream input_file ( file_edge );
        size_t f_nedges ;
        // Open the input file
        if ( ! input_file . good () )
        {
            std::cerr << "File " << file_edge << " not found.\n";
            throw std::runtime_error("File Parsing Error: File ! found");
        }
        // First line is the number of edges
        std::size_t line_number = 0;
        if (! input_file.eof())
        {
            ++line_number ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_edge) ;
            // First number is the number of nodes
            std::istringstream is (line);
            is >> f_nedges ;
        }
        while ( !(input_file.eof()) && (f_nedges>0))
        {
            size_t trash, i, j ;
            ++ line_number;
            --f_nedges ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_edge) ;
            std::istringstream is (line);
            is >> trash ;
            is >> i ;
            is >> j ;
            Io_cell_type cell({i,j}) ;
            add_cell(cell) ;
        }
        input_file.close() ;
        std::cout << "--- " << f_nedges << " edges" << std::endl ;
    }

    void add_faces()
    {
        const std::string file_face = ffaces_from_prefix(_prefix) ;
        std::ifstream input_file ( file_face );
        size_t f_nfaces ;
        // Open the input file
        if ( ! input_file . good () )
        {
            std::cerr << "File " << file_face << " not found.\n";
            throw std::runtime_error("File Parsing Error: File ! found");
        }
        // First line is the number of edges
        std::size_t line_number = 0;
        if (! input_file.eof())
        {
            ++line_number ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_face) ;
            // First number is the number of nodes
            std::istringstream is (line);
            is >> f_nfaces ;
        }
        while ( !(input_file.eof()) && (f_nfaces>0))
        {
            size_t trash, i, j, k ;
            ++ line_number;
            --f_nfaces ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_face) ;
            std::istringstream is (line);
            is >> trash ;
            is >> i ;
            is >> j ;
            is >> k ;
            Io_cell_type cell({i, j, k}) ;
            add_cell(cell) ;
        }
        input_file.close() ;
        std::cout << "--- " << f_nfaces << " faces" << std::endl ;
    }

    void add_tets()
    {
        const std::string file_ele = ftets_from_prefix(_prefix) ;
        std::ifstream input_file ( file_ele );
        size_t f_ntets, tmp ;
        // Open the input file
        if ( ! input_file . good () )
        {
            std::cerr << "File " << file_ele << " not found.\n";
            throw std::runtime_error("File Parsing Error: File ! found");
        }
        // First line is the number of edges
        std::size_t line_number = 0;
        if (! input_file.eof())
        {
            ++line_number ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_ele) ;
            // First number is the number of nodes
            std::istringstream is (line);
            is >> f_ntets ;
            
        }
        tmp = f_ntets ;
        while ( !(input_file.eof()) && (tmp>0))
        {
            size_t trash, i, j, k, l ;
            ++ line_number;
            --tmp ;
            std::string line;
            getline( input_file, line );
            check_sanity_line(line, file_ele) ;
            std::istringstream is (line);
            is >> trash ;
            is >> i ;
            is >> j ;
            is >> k ;
            is >> l ;
            Io_cell_type tmp_cell({i, j, k, l});
            add_cell(tmp_cell, true) ;
        }
        input_file.close() ;
        std::cout << "--- " << f_ntets << " tets" << std::endl ;
    }
private:
    std::string _prefix ;
    std::string fnodes_from_prefix(const std::string &prefix) {return prefix+".1.node"; }
    std::string fedges_from_prefix(const std::string &prefix) {return prefix+".1.edge"; }
    std::string ffaces_from_prefix(const std::string &prefix) {return prefix+".1.face"; }
    std::string ftets_from_prefix(const std::string &prefix) {return prefix+".1.ele"; }
} ;

} /* end namespace HDVF */
} /* end namespace CGAL */


#endif // CGAL_HDVF_TET_OBJECT_H
