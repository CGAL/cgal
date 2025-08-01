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

#include <iostream>

#include <CGAL/HDVF/tools_io.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

// IO tools

namespace CGAL {
namespace HDVF {

inline std::ostream & operator<< (std::ostream &out, const PairCell& p)
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
 
template <typename CoefType, typename ComplexType>
void interaction_loop(Hdvf<CoefType, ComplexType> &hdvf,
                      ComplexType &complex,
                      const std::function<void(Hdvf<CoefType, ComplexType> &hdvf, ComplexType &complex)> &output_vtk)
{
    bool over = false ;
    output_vtk(hdvf, complex) ;
    while (!over)
    {
        string instr ;
        int ipair ;
        cout << "Next instruction : M, W, MW or Q (to quit)" << endl ;
        cin >> instr ;
        if (instr == string("M"))
        {
            int sigma, q ;
            cout << "Provide sigma / q (separate with space):" << endl ;
            cin >> sigma >> q ;
            bool found = false ;
            vector<PairCell> possM(hdvf.find_pairs_M(q, found, sigma)) ;
            std::cout << "-> possible pairings for M:" << std::endl ;
            for (int i = 0; i<possM.size(); ++i)
                std::cout << i << " : " << possM.at(i) << std::endl ;
            if (possM.size() > 0)
            {
                cout << "chose indice of pairing (out the range = quit): " ;
                cin >> ipair ;
                if ((ipair >= 0) && (ipair < possM.size()))
                {
                    PairCell p(possM.at(ipair)) ;
                    std::cout << "M(" << p << ")" << std::endl ;
                    hdvf.M(p.sigma, p.tau, p.dim) ;
                    //                    hdvf.print_matrices() ;
                    output_vtk(hdvf, complex) ;
                }
            }
        }
        else if (instr == string("W"))
        {
            int sigma, q ;
            cout << "Provide sigma / q (separate with space):" << endl ;
            cin >> sigma >> q ;
            bool found = false ;
            vector<PairCell> possW(hdvf.find_pairs_W(q, found, sigma)) ;
            cout << "-> possible pairings for W:" << endl ;
            for (int i = 0; i<possW.size(); ++i)
                cout << i << " : " << possW.at(i) << endl ;
            if (possW.size() > 0)
            {
                cout << "chose indice of pairing (out the range = quit): " ;
                cin >> ipair ;
                if ((ipair >= 0) && (ipair < possW.size()))
                {
                    PairCell p(possW.at(ipair)) ;
                    cout << "W(" << p << ")" << endl ;
                    hdvf.W(p.sigma, p.tau, p.dim) ;
                    //                    hdvf.print_matrices() ;
                    output_vtk(hdvf, complex) ;
                }
            }
        }
        else if (instr == string("MW"))
        {
            int sigma, q ;
            cout << "Provide sigma / q (separate with space):" << endl ;
            cin >> sigma >> q ;
            bool found = false ;
            vector<PairCell> possMW(hdvf.find_pairs_MW(q, found, sigma)) ;
            cout << "-> possible pairings for MW:" << endl ;
            for (int i = 0; i<possMW.size(); ++i)
                cout << i << " : " << possMW.at(i) << endl ;
            if (possMW.size() > 0)
            {
                cout << "chose indice of pairing (out the range = quit): " ;
                cin >> ipair ;
                if ((ipair >= 0) && (ipair < possMW.size()))
                {
                    PairCell p(possMW.at(ipair)) ;
                    cout << "MW(" << p << ")" << endl ;
                    hdvf.MW(p.sigma, p.tau, p.dim) ;
                    //                    hdvf.print_matrices() ;
                    output_vtk(hdvf, complex) ;
                }
            }
        }
        else if (instr == string("Q"))
            over = true ;
    }
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_HDVF_TOOLS_H
