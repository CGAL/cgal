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

#ifndef CGAL_HDVF_CUB_OBJECT_H
#define CGAL_HDVF_CUB_OBJECT_H

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>

namespace CGAL {
namespace HDVF {

using namespace std ;

// ------ For cubical complexes
/** \brief Type of cells coordinates in Cub_object_io (Khalimsky or voxel coordinates) */
typedef std::vector<size_t> IOCubCellType ;
/** \brief Type of pre-chains in Cub_object_io (list of cells without coefficients). */
typedef std::vector<IOCubCellType> IOCubChainType ;

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Cub_object_io` is an intermediate IO class, used to load binary volumes and produce cubical complexes.

 */

class Cub_object_io
{
public:
    int dim = 0 ; // Dimension of the complex
    vector<size_t> ncubs ; // Number of cubs in each dimension
    vector<size_t> N ; // Size of BB along each dimension
    std::vector<IOCubCellType> cubs ;
    bool khalimsky ;

    /* \brief Default constructor.
     *
     * Create an empty Cub_object_io of dimension 3.
     */
    Cub_object_io() : dim(3), ncubs(vector<size_t>(3)), N(vector<size_t>(3)), khalimsky(false) {}

    /**
     * \brief Constructor from a vector of cells.
     *
     * Cells coordinates are given in Khalimsky coordinates if the boolean `khal` is `true`, and as integer indices of voxels otherwise.
     *
     */
    Cub_object_io(int d, const std::vector<IOCubCellType> &vcubs, bool khal = false) : dim(d), ncubs(vector<size_t>(d)), N(vector<size_t>(d)), cubs(vcubs), khalimsky(khal)
    { check_dimension() ;}

    /* \brief Copy constructor. */
    Cub_object_io(const Cub_object_io &m) : dim(m.dim), ncubs(m.ncubs), N(m.N), cubs(m.cubs), khalimsky(m.khalimsky) {}

    // Mesh operations
    void clear_cubs() { cubs.clear() ; for (size_t i=0; i<dim; ++i) ncubs[i] = 0 ; }
    void add_cub(const IOCubCellType &c) {cubs.push_back(c); ++ncubs[cub_dim(c)] ;}
    void frame() // Enlarge the bouding box to add 1 voxel around
    {
        // Enlarge the BB along each dimension
        for (size_t i=0; i<N.size(); ++i)
        {
            if (khalimsky)
                N.at(i) += 4 ;
            else
                N.at(i) += 2 ;
        }
        // Shift cells coordinates
        for (size_t i = 0; i<cubs.size(); ++i)
        {
            for (int k = 0; k<dim; ++k)
            {
                if (khalimsky)
                    cubs.at(i).at(k) += 2 ;
                else
                    ++cubs.at(i).at(k) ;
            }
        }
    }

    // PGM
    bool read_pgm(const std::string &filename, bool khal = false)
    {
        khalimsky = khal ;
        size_t row = 0, col = 0, numrows = 0, numcols = 0;
        ifstream infile(filename);
        if(!infile.is_open()) {
            cerr << "Warning: cannot open file " << filename << endl ;
            // failed to open the file
            return false;
        }

        stringstream ss;
        string inputLine = "";

        // First line : version
        getline(infile,inputLine);
        if(inputLine.compare("P2") != 0) cerr << "Version error" << endl;
        else cout << "Version : " << inputLine << endl;


        // Second line: dimensions
        getline(infile,inputLine);
        vector<size_t> sizes ;

        size_t tmp ;
        stringstream sseizes (inputLine);
        while (sseizes >> tmp)
            sizes.push_back(tmp) ;

        cout << "dimensions : " ;
        for (size_t i=0; i<sizes.size(); ++i)
            cout << sizes.at(i) << " " ;
        cout << endl ;
        dim = sizes.size() ;
        N = vector<size_t>(dim) ;
        for (size_t i=0; i<dim; ++i)
            N.at(i) = sizes.at(i) ;

        // Throw away next data (255)

        getline(infile,inputLine);

        // Remainder: data
        // Read by decreasing dimension
        // In order to get : read row, then column, then depth...

        // Continue with a stringstream
        ss << infile.rdbuf();

        size_t NN = 1 ;
        for (int i=0; i<dim; ++i)
            NN = NN * N.at(i) ;

        // Following lines : data
        for(size_t i = 0; i < NN; ++i)
        {
            ss >> tmp ;
            if (tmp) // pixel on
            {
                cubs.push_back(index_to_coords(i, khal)) ;
                ++ncubs[dim] ;
            }
        }

        if (khal)
        {
            // Change N to Khalimsky maxima
            for(int i=0; i<dim; ++i)
                N.at(i) = 2*N.at(i)+1 ;
        }

        infile.close();
    }

    bool write_pgm(const std::string &filename) ;
    // CUB
    bool read_cub(const std::string &filename, bool khalimsky = false)
    {
        // 0 - open input file
        std::ifstream infile(filename);
        if(!infile.is_open()) {
            // failed to open the file
            cerr << "Warning: file " << filename << " does not exist" << endl ;
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
            std::istringstream is( line );

            if (line_number == 1) // Header 1
                is >> dim ;
            else if (line_number == 2) // Header 2
            {
                for (int i=0; i<dim; ++i)
                    is >> N[i] ;
            }
            else // Cub line
            {
                IOCubCellType cub ;

                size_t v;
                while ( is >> v )
                    cub.push_back(v);
                if (cub.size() != dim)
                    cout << "Discard line " << line_number << " cub of invalid dimension" << endl ;
                else
                {
                    // Add this cub to cubs
                    ++ncubs[cub_dim(cub)] ;

                    if (khalimsky)
                        cubs.push_back(cub) ;
                    else
                    {
                        size_t ind(khal_to_index(cub)) ;
                        cubs.push_back(index_to_coords(ind, false)) ;
                    }
                }
            }
        }
        infile.close() ;
        return true ;
    }

    bool write_cub(const std::string &filename) ;

    void print_infos (size_t level = 0) const
    {
        cout << "Cub_object_io infos - dim : " << dim << ", cubs : " << cubs.size() << endl ;
        for (int q=0; q < dim; ++q)
            cout << "\tSize along dim " << q << " : " << N[q] << endl ;
        for (int q = 0; q <= dim; ++q)
            cout << "Cubs of dim " << q << " : " << ncubs[q] << endl ;
        if (level == 1) // Print coordinates of cubes
        {
            for (size_t i=0; i<cubs.size(); ++i)
            {
                for (size_t j : cubs.at(i))
                    cout << j << " " ;
                cout << endl ;
            }
        }
    }

private:
    void check_dimension()
    {
        if (khalimsky)
        {
            for (IOCubCellType c : cubs)
                ++ncubs[cub_dim(c)] ;
        }
        else
        {
            ncubs[dim] = cubs.size() ;
        }
    }

    inline int cub_dim (IOCubCellType c)
    {
        int q = 0 ;
        for (size_t i : c)
        {
            if (i%2 == 1)
                ++q ;
        }
        return q ;
    }

    inline IOCubCellType index_to_coords(size_t i, bool khalimsky = false)
    {
        IOCubCellType coords ;
        // Convert index in binary object to size_t coordinates
        for (int q = 0; q < dim; ++q)
        {
            coords.push_back(i % N[q]) ;
            i = i / (N[q]) ;
        }
        if (khalimsky)
        {
            // Convert int coordinates to Khalimsky
            for (int q = 0; q < dim; ++q)
            {
                coords.at(q) = 2*coords.at(q)+1 ;
            }
        }
        return coords ;
    }

    inline size_t khal_to_index (const IOCubCellType& coords)
    {
        size_t cell_index(0);
        for (int i = 0; i < dim; ++i) {
            cell_index += coords[i] * N[i];
        }
    }
} ;

} /* end namespace HDVF */
} /* end namespace CGAL */


#endif // CGAL_HDVF_CUB_OBJECT_H
