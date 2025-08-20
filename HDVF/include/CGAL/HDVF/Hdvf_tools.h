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

#ifndef CGAL_HDVF_HDVF_TOOLS_H
#define CGAL_HDVF_HDVF_TOOLS_H

#include <CGAL/license/HDVF.h>

#include <iostream>

#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

// IO tools

namespace CGAL {
namespace HDVF {

inline std::ostream & operator<< (std::ostream &out, const Cell_pair& p)
{
    out << "(" << p.sigma << "/" << p.tau << " - dim " << p.dim << ")" ;
    return out ;
}

/**
 * \brief Runs an interaction loop to iterated M, W or MW operations and export the results to vtk.
 *
 * The loop runs until the key `Q` is pressed. Otherwise, the loop asks for an operation (M, W or MW) and a cell (index and dimension).
 * Then all possible paired cells are listed and the user can chose one of them (or none).
 */

template <typename ComplexType>
void interaction_loop(Hdvf<ComplexType> &hdvf,
                      ComplexType &complex,
                      const std::function<void(Hdvf<ComplexType> &hdvf, ComplexType &complex)> &output_vtk)
{
    bool over = false ;
    output_vtk(hdvf, complex) ;
    while (!over)
    {
        std::string instr ;
        int ipair ;
        std::cout << "Next instruction : M, W, MW or Q (to quit)" << std::endl ;
        std::cin >> instr ;
        if (instr == std::string("M"))
        {
            int sigma, q ;
            std::cout << "Provide sigma / q (separate with space):" << std::endl ;
            std::cin >> sigma >> q ;
            bool found = false ;
            std::vector<Cell_pair> possM(hdvf.find_pairs_M(q, found, sigma)) ;
            std::cout << "-> possible pairings for M:" << std::endl ;
            for (int i = 0; i<possM.size(); ++i)
                std::cout << i << " : " << possM.at(i) << std::endl ;
            if (possM.size() > 0)
            {
                std::cout << "chose indice of pairing (out the range = quit): " ;
                std::cin >> ipair ;
                if ((ipair >= 0) && (ipair < possM.size()))
                {
                    Cell_pair p(possM.at(ipair)) ;
                    std::cout << "M(" << p << ")" << std::endl ;
                    hdvf.M(p.sigma, p.tau, p.dim) ;
                    //                    hdvf.insert_matrices() ;
                    output_vtk(hdvf, complex) ;
                }
            }
        }
        else if (instr == std::string("W"))
        {
            int sigma, q ;
            std::cout << "Provide sigma / q (separate with space):" << std::endl ;
            std::cin >> sigma >> q ;
            bool found = false ;
            std::vector<Cell_pair> possW(hdvf.find_pairs_W(q, found, sigma)) ;
            std::cout << "-> possible pairings for W:" << std::endl ;
            for (int i = 0; i<possW.size(); ++i)
                std::cout << i << " : " << possW.at(i) << std::endl ;
            if (possW.size() > 0)
            {
                std::cout << "chose indice of pairing (out the range = quit): " ;
                std::cin >> ipair ;
                if ((ipair >= 0) && (ipair < possW.size()))
                {
                    Cell_pair p(possW.at(ipair)) ;
                    std::cout << "W(" << p << ")" << std::endl ;
                    hdvf.W(p.sigma, p.tau, p.dim) ;
                    //                    hdvf.insert_matrices() ;
                    output_vtk(hdvf, complex) ;
                }
            }
        }
        else if (instr == std::string("MW"))
        {
            int sigma, q ;
            std::cout << "Provide sigma / q (separate with space):" << std::endl ;
            std::cin >> sigma >> q ;
            bool found = false ;
            std::vector<Cell_pair> possMW(hdvf.find_pairs_MW(q, found, sigma)) ;
            std::cout << "-> possible pairings for MW:" << std::endl ;
            for (int i = 0; i<possMW.size(); ++i)
                std::cout << i << " : " << possMW.at(i) << std::endl ;
            if (possMW.size() > 0)
            {
                std::cout << "chose indice of pairing (out the range = quit): " ;
                std::cin >> ipair ;
                if ((ipair >= 0) && (ipair < possMW.size()))
                {
                    Cell_pair p(possMW.at(ipair)) ;
                    std::cout << "MW(" << p << ")" << std::endl ;
                    hdvf.MW(p.sigma, p.tau, p.dim) ;
                    //                    hdvf.insert_matrices() ;
                    output_vtk(hdvf, complex) ;
                }
            }
        }
        else if (instr == std::string("Q"))
            over = true ;
    }
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_HDVF_TOOLS_H
