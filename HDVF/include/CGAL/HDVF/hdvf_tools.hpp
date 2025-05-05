#ifndef HDVF_TOOLS_HPP
#define HDVF_TOOLS_HPP

#include <iostream>

//#include "hdvf_common.hpp"
//#include "Simplex.hpp"
#include "tools_io.hpp"
//#include "AbstractSimpComplex.hpp"
//#include "CubComplex.hpp"
#include "hdvf.hpp"
#include "CGAL/OSM/OSM.hpp"

namespace CGAL {
namespace HDVF {

// IO tools

inline ostream & operator<< (ostream &out, const PairCell &p)
{
    out << "(" << p.sigma << "/" << p.tau << " - dim " << p.dim << ")" ;
    return out ;
}


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
            cout << "-> possible pairings for M:" << endl ;
            for (int i = 0; i<possM.size(); ++i)
                cout << i << " : " << possM.at(i) << endl ;
            if (possM.size() > 0)
            {
                cout << "chose indice of pairing (out the range = quit): " ;
                cin >> ipair ;
                if ((ipair >= 0) && (ipair < possM.size()))
                {
                    PairCell p(possM.at(ipair)) ;
                    cout << "M(" << p << ")" << endl ;
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

#endif // HDVF_TOOLS_HPP
