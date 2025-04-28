//
//  hdvf_duality.hpp
//  HDVF
//
//  Created by umenohana on 27/08/2024.
//

#ifndef HDVF_DUALITY_H
#define HDVF_DUALITY_H

#include <vector>
#include <cassert>
#include <iostream>
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/OSM/Bitboard.hpp"
#include "CGAL/HDVF/Abstract_simplicial_chain_complex.hpp"
#include "CGAL/HDVF/Simplex.hpp"
#include "CGAL/HDVF/Cubical_complex.hpp"
#include "CGAL/HDVF/hdvf.hpp"
#include "CGAL/HDVF/sub_chain_complex.hpp"
#include "CGAL/HDVF/SubSparseMatrix.hpp"


/**
 * \class Dual_HDVF
 * \brief Implementation of HDVF duality and associate operations.
 *
 * The Dual_HDVF class contains all functions to build HDVF for Alexander duality.
 *
 * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam _ComplexType The type of complex  (default is SimpComplex)
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 */

template<typename _CoefficientType, typename _ComplexType>
class Dual_HDVF : public HDVF<_CoefficientType, _ComplexType, OSM::Chain, OSM::SubSparseMatrix> {
private:
    // Matrices types
    typedef OSM::Chain<_CoefficientType, OSM::COLUMN> CChain;
    typedef OSM::Chain<_CoefficientType, OSM::ROW> RChain;
    
    // HDVF Type
    typedef HDVF<_CoefficientType, _ComplexType, OSM::Chain, OSM::SubSparseMatrix> HDVF_type ;
    
    // Complex L
    const _ComplexType& _L ;
    const int _hdvf_opt ;
    // Subcomplex K
    SubChainComplex<_CoefficientType,_ComplexType> _KCC, _subCC ;
public:
    // Critical cells of perfect HDVFs
    std::vector<std::vector<int> > critical_K, critical_L_K ;
    
    /**
     * \brief Build Dual_HDVF from a complex and a sub chain complex
     *
     * \warning The Dual_HDVF will check that the sub chain complex is compatible. TODO !!!
     *
     * \param[in] L The (full) complex.
     * \param[in] K The sub chain complex.
     *
     * \see \link OSM::Dual_HDVF \endlink
     *
     * \author Bac A.
     * \version 0.1.0
     * \date 28/08/2024
     */
    Dual_HDVF(const _ComplexType& L, SubChainComplex<_CoefficientType, _ComplexType>& K, int hdvf_opt = OPT_FULL) ;
    
    /** \brief find a valid PairCell for A in dimension q */
    virtual PairCell findPairA(int q, bool &found); // const;
    /** \brief find a valid PairCell containing tau for A in dimension q */
    virtual PairCell findPairA(int q, bool &found, int tau); // const;
    /** \brief find all the valid PairCell for A in dimension q */
    virtual std::vector<PairCell> findPairsA(int q, bool &found); // const;
    /** \brief find all the valid PairCell containing tau for A in dimension q */
    virtual std::vector<PairCell> findPairsA(int q, bool &found, int tau); // const;
    
    /** \brief Set _DD_col to _subCC */
    inline void screen_DD_col ()
    {
        for (int q=0; q<=_L.dim(); ++q)
        {
            this->_DD_col.at(q).set_sub(_subCC.get_BitBoard(q)) ;
        }
    }
    
    /** \brief Set _subCC to complex K and screen _DD_col accordingly */
    inline void set_mask_K ()
    {
        _subCC = _KCC ;
        screen_DD_col() ;
    }
    
    /** \brief Set _subCC to cocomplex L-K and screen _DD_col accordingly */
    inline void set_mask_L_K ()
    {
        _subCC = _KCC.complement() ;
        screen_DD_col() ;
    }
    
    /**
     * \brief Compute "perfect" HDVF (if they exist) over both sub complexes + pairing matrix
     *
     * \see \link OSM::Dual_HDVF \endlink
     *
     * \author Bac A.
     * \version 0.1.0
     * \date 28/08/2024
     */
    void computeDualPerfectHDVF() ;
    
    /**
     * \brief Compute "pairing" HDVF between K and L-K
     *
     * \see \link OSM::Dual_HDVF \endlink
     *
     * \author Bac A.
     * \version 0.1.0
     * \date 28/08/2024
     */
    std::vector<PairCell> computePairingHDVF() ;
    
    // Dual_HDVF getters
    // Method to get cells if with a given flag (P,S,C) for each dimension
    vector<vector<int> > get_flag (FlagType flag) const ;
    // Method to get cells with a given flag (P,S,C) for a given dimension
    vector<int> get_flag_dim (FlagType flag, int q) const ;
    // Method to get the flag (P,S,C) of a cell tau of dimension q
    
    // Dual_HDVF I/O
    
    virtual ostream& print_reduction(ostream& out = cout)
    {
        // Print K
        out << "----> K" << endl ;
        _subCC = _KCC ;
        screen_DD_col() ;
        print_reduction_sub(out) ;
        
        // Print L-K
        _subCC = _KCC.complement() ;
        screen_DD_col() ;
        print_reduction_sub(out) ;
        
        // Set back _subCC to K
        _subCC = _KCC ;
        
        return out ;
    }
    
    virtual ostream& print_reduction_sub(ostream& out = cout) // const;
    {
        // Print critical cells
        out << "----- critical cells:" << endl;
        for (int q = 0; q <= _L.dim(); ++q) {
            out << "--- dim " << q << endl;
            for (int i = 0; i < _L.nb_cells(q); ++i)
            {
                if ((this->_flag[q][i] == CRITICAL) && (_subCC.get_Bit(q, i))) {
                    out << i << " ";
                }
            }
            out << endl;
        }
        
        if (_hdvf_opt & (OPT_FULL | OPT_G))
        {
            // Print matrices g
            out << "----- g:" << endl;
            for (int q = 0; q <= _L.dim(); ++q) {
                out << "--- dim " << q << endl;
                for (int i = 0; i < _L.nb_cells(q); ++i) {
                    if ((this->_flag[q][i] == CRITICAL) && (_subCC.get_Bit(q, i))) {
                        out << "g(" << i << ") = (" << i << ")";
                        // Iterate over the ith column of _G_col
                        typename HDVF_type::CChain col(OSM::getColumn(this->_G_col.at(q), i)) ; // TODO cget
                        for (typename HDVF_type::CChain::const_iterator it_col = col.cbegin(); it_col != col.cend(); ++it_col) {
                            out << " + " << it_col->second << ".(" << it_col->first << ") + ";
                        }
                        out << endl;
                    }
                }
            }
        }
        
        if (_hdvf_opt & (OPT_FULL | OPT_F))
        {
            // Print matrices f*
            out << "----- f*:" << endl;
            for (int q = 0; q <= _L.dim(); ++q) {
                out << "--- dim " << q << endl;
                for (int i = 0; i < _L.nb_cells(q); ++i) {
                    if ((this->_flag[q][i] == CRITICAL) && (_subCC.get_Bit(q, i))) {
                        out << "f*(" << i << ") = (" << i << ")";
                        // Iterate over the ith row of _F_row
                        typename HDVF_type::RChain row(OSM::getRow(this->_F_row.at(q), i)) ; // TODO cget
                        for (typename HDVF_type::RChain::const_iterator it_row = row.cbegin(); it_row != row.cend(); ++it_row) {
                            out << " + " << it_row->second << ".(" << it_row->first << ") + ";
                        }
                        out << endl;
                    }
                }
            }
        }
        return out ;
    }
    
    // Method to print reduced boundary over the pairing
    
    ostream& print_bnd_pairing(ostream& out = cout)
    {
        SubChainComplex<_CoefficientType, _ComplexType> subPair(_L, false) ;
        for (int q=0; q<=_L.dim(); ++q)
        {
            for (int i=0; i<critical_K.at(q).size(); ++i)
                subPair.set_BitOn(q, critical_K.at(q).at(i)) ;
            for (int i=0; i<critical_L_K.at(q).size(); ++i)
                subPair.set_BitOn(q, critical_L_K.at(q).at(i)) ;
        }
        // Print corresponding submatrices _DD_col
        for (int q=1; q<=_L.dim(); ++q)
        {
            out << "--> dim " << q << " : q / q-1 cells" << endl ;
            out << "id " << q << " : " ;
            for (typename OSM::Bitboard::iterator it = subPair.get_BitBoard(q).begin(); it != subPair.get_BitBoard(q).end(); ++it)
                out << *it << " " ;
            out << endl ;
            out << "id " << q-1 << " : " ;
            for (typename OSM::Bitboard::iterator it = subPair.get_BitBoard(q-1).begin(); it != subPair.get_BitBoard(q-1).end(); ++it)
                out << *it << " " ;
            out << endl ;
            for (typename OSM::Bitboard::iterator it = subPair.get_BitBoard(q).begin(); it != subPair.get_BitBoard(q).end(); ++it)
            {
                for (typename OSM::Bitboard::iterator it2 = subPair.get_BitBoard(q-1).begin(); it2 != subPair.get_BitBoard(q-1).end(); ++it2)
                {
                    if (this->_DD_col.at(q).get_coef(*it2, *it) == 0)
                        out << ".\t" ;
                    else
                        out << this->_DD_col.at(q).get_coef(*it2, *it) << "\t" ;
                }
                out << endl ;
            }
        }
        return out ;
    }
    
    // Method to generate labels for visualisation
    // For PSC, no additional arguments, for FSTAR and G, specify the critical cell concerned
    vector<vector<int> > export_label (ExportType type, int cell=0, int dim=0) const
    {
        vector<vector<int> > labels(this->_K.dim()+1) ;
        if (type == PSC)
        {
            for (int q=0; q<=this->_K.dim(); ++q)
            {
                for (int i = 0; i<this->_K.nb_cells(q); ++i)
                {
                    if (_subCC.get_Bit(q, i)) // i belongs to _subCC
                    {
                        if (this->_flag.at(q).at(i) == PRIMARY)
                            labels.at(q).push_back(-1) ;
                        else if (this->_flag.at(q).at(i) == SECONDARY)
                            labels.at(q).push_back(1) ;
                        else
                            labels.at(q).push_back(0) ;
                    }
                    else // i does not belongs to _subCC
                        labels.at(q).push_back(2) ;
                }
            }
        }
        else
        {
            if ((type == FSTAR) && (_hdvf_opt & (OPT_FULL | OPT_F)))
            {
                for (int q=0; q<=this->_K.dim(); ++q)
                {
                    labels.at(q).resize(this->_K.nb_cells(q)) ;
                }
                
                if (dim < this->_K.dim())
                {
                    labels.at(dim).at(cell) = 1 ;
                    const typename HDVF_type::RChain& fstar_cell = OSM::cgetRow(this->_F_row.at(dim), cell) ;
                    for (typename HDVF_type::RChain::const_iterator it = fstar_cell.cbegin(); it != fstar_cell.cend(); ++it)
                    {
                        // Keep cells of the chain belonging to _subCC
                        if (_subCC.get_Bit(dim, it->first))
                        {
                            // Set the cofaces belonging to _subCC of it->first in dimension dim+1
                            typename HDVF_type::RChain cofaces(this->_K.cod(it->first,dim)) ;
                            for (typename HDVF_type::RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                            {
                                if (_subCC.get_Bit(dim+1, it2->first))
                                    labels.at(dim+1).at(it2->first) = 1 ; // set the flag to 1 only for FSTAR(cell) cells
                            }
                        }
                    }
                }
                else
                    throw "Error: cannot export f* for a cell of maximal dimension" ;
            }
            else if (_hdvf_opt & (OPT_FULL | OPT_G)) // G
            {
                //                cout << "Export : _subCC : " << _subCC << endl ;
                for (int q=0; q<=this->_K.dim(); ++q)
                {
                    labels.at(q).resize(this->_K.nb_cells(q)) ;
                }
                
                labels.at(dim).at(cell) = 1 ;
                const typename HDVF_type::CChain& g_cell = OSM::cgetColumn(this->_G_col.at(dim), cell) ;
                for (typename HDVF_type::CChain::const_iterator it = g_cell.cbegin(); it != g_cell.cend(); ++it)
                {
                    // Keep cells of the chain belonging to _subCC
                    if (_subCC.get_Bit(dim, it->first))
                    {
                        //                        cout << "g(" << cell << ", " << dim << ") : " << it->first << endl ;
                        labels.at(dim).at(it->first) = 1 ; // set the flag to 1 only for G(cell) cells
                    }
                    
                }
                
            }
        }
        return labels ;
    }
    
    
    // Method to generate PSC labels for visualisation
    virtual vector<vector<int> > export_labelsPSC () const
    {
        vector<vector<int> > labels(this->_K.dim()+1) ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            for (int i = 0; i<this->_K.nb_cells(q); ++i)
            {
                if (_subCC.get_Bit(q, i)) // i belongs to _subCC
                {
                    if (this->_flag.at(q).at(i) == PRIMARY)
                        labels.at(q).push_back(-1) ;
                    else if (this->_flag.at(q).at(i) == SECONDARY)
                        labels.at(q).push_back(1) ;
                    else
                        labels.at(q).push_back(0) ;
                }
                else // i does not belongs to _subCC
                    labels.at(q).push_back(2) ;
            }
        }
        return labels ;
    }
    
    // Method to export a G chain for visualisation
    // with cell indices in initial _K
    virtual CChain export_GChain (int cell, int dim) const
    {
        if ((dim<0) || (dim>this->_K.dim()))
            throw "Error : export_GChain with dim out of range" ;
        //        if (_K_to_per.at(dim).at(cell) > _t_dim.at(dim))
        //            throw "Error : export_GChain with 'future' cell wrt persistence" ;
        
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            // Get g(cell, dim) with per indices
            CChain g_cell(OSM::getColumn(this->_G_col.at(dim), cell)) ;
            // Add 1 to the cell
            g_cell[cell] = 1 ;
            // Keep cells of the chain belonging to _subCC
            CChain g_cell_sub(g_cell.dimension()) ;
            for (typename CChain::const_iterator it = g_cell.begin(); it != g_cell.end(); ++it)
            {
                
                if (_subCC.get_Bit(dim, it->first))
                {
                    g_cell_sub[it->first] = it->second ;
                }
            }
            return g_cell_sub ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
    
    // Method to export a FSTAR chain for visualisation
    // with cell indices in initial _K
    virtual CChain export_FSTARChain (int cell, int dim) const
    {
        if ((dim<0) || (dim>this->_K.dim()))
            throw "Error : export_FSTARChain with dim out of range" ;
        //        if (_K_to_per.at(dim).at(cell) > _t_dim.at(dim))
        //            throw "Error : export_GChain with 'future' cell wrt persistence" ;
        
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            RChain fstar_cell(OSM::getRow(this->_F_row.at(dim), cell)) ;
            // Add 1 to the cell
            fstar_cell[cell] = 1 ;
            // Compute the cofaces
            if (dim < this->_K.dim())
            {
                CChain fstar_cofaces(this->_K.nb_cells(dim+1)) ;
                for (typename RChain::const_iterator it = fstar_cell.cbegin(); it != fstar_cell.cend(); ++it)
                {
                    // Set the cofaces of it->first belonging to _subCC in dimension dim+1
                    RChain cofaces(this->_K.cod(it->first,dim)) ;
                    for (typename RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                    {
                        if (_subCC.get_Bit(dim+1, it2->first))
                            fstar_cofaces[it2->first] = 1 ;
                    }
                }
                return fstar_cofaces ;
            }
            else
                return CChain(0) ;
        }
    }
            
//            // Get fstar(cell, dim) with per indices
//            RChain fstar_cell(OSM::cgetRow(this->_F_row.at(dim), cell)) ;
//            // Add 1 to the cell
//            fstar_cell[cell] = 1 ;
//            // Keep cells of the chain belonging to _subCC
//            // Compute the chain (transpose of F_STAR column)
//            CChain fstar_cell_sub(fstar_cell.dimension()) ;
//            for (typename RChain::const_iterator it = fstar_cell.begin(); it != fstar_cell.end(); ++it)
//            {
//                // Keep cells of the chain belonging to _subCC
//                if (_subCC.get_Bit(dim, it->first))
//                {
//                    fstar_cell_sub[it->first] = it->second ;
//                }
//            }
//            return fstar_cell_sub ;
//        }
//        else
//            throw "Error : trying to export g_chain without proper HDVF option" ;
    
} ;

template<typename _CoefficientType, typename _ComplexType>
Dual_HDVF<_CoefficientType,_ComplexType>::Dual_HDVF(const _ComplexType& L, SubChainComplex<_CoefficientType, _ComplexType>& K, int hdvf_opt) :
HDVF<_CoefficientType, _ComplexType, OSM::Chain, OSM::SubSparseMatrix>(L,hdvf_opt), _L(L), _hdvf_opt(hdvf_opt), _KCC(K), _subCC(K) {}

/** \brief find a valid PairCell for A in dimension q */
template<typename _CoefficientType, typename _ComplexType>
PairCell Dual_HDVF<_CoefficientType,_ComplexType>::findPairA(int q, bool &found) // const
{
    found = false;
    PairCell p;
    
    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = this->_DD_col[q+1].begin(); (it_col != this->_DD_col[q+1].end() && !found); ++it_col)
    {
        const typename HDVF_type::CChain& col(OSM::cgetColumn(this->_DD_col[q+1], *it_col)) ;
        
        // Iterate through the entries of the column
        // Check that the row belongs to the subchaincomplex
        for (typename HDVF_type::CChain::const_iterator it = col.begin(); (it != col.end() && !found); ++it) {
            if (_subCC.get_Bit(q, it->first) && (abs(it->second) == 1)) {
                // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
                p.sigma = it->first;
                p.tau = *it_col;
                p.dim = q;
                found = true;
            }
        }
    }
    return p;
}

/** \brief find a valid PairCell containing tau for A in dimension q */
template<typename _CoefficientType, typename _ComplexType>
PairCell Dual_HDVF<_CoefficientType,_ComplexType>::findPairA(int q, bool &found, int tau) // const
{
    found = false;
    PairCell p ;
    // Check tau belongs to _subCC
    if (!_subCC.get_Bit(q, tau))
        throw("Dual_HDVF: searching for a cell tau outside _subCC") ;
    
    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    // and tau' belongs to _subCC
    const typename HDVF_type::CChain& tmp2(OSM::cgetColumn(this->_DD_col.at(q), tau)) ;
    for (typename HDVF_type::CChain::const_iterator it = tmp2.cbegin(); (it != tmp2.cend() && !found); ++it)
    {
        if (_subCC.get_Bit(q-1, it->first) && abs(it->second) == 1)
        {
            found = true ;
            p.sigma = it->first ;
            p.tau = tau ;
            p.dim = q-1 ;
        }
    }
    
    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    // and tau' belongs to _subCC
    typename HDVF_type::RChain tmp(OSM::getRow(this->_DD_col.at(q+1), tau)) ;
    for (typename HDVF_type::RChain::const_iterator it = tmp.cbegin(); (it != tmp.cend() && !found); ++it)
    {
        if (_subCC.get_Bit(q+1, it->first) && (abs(it->second) == 1))
        {
            found = true ;
            PairCell p ;
            p.sigma = tau ;
            p.tau = it->first ;
            p.dim = q ;
        }
    }
    return p;
}

/** \brief find all the valid PairCell for A in dimension q */
template<typename _CoefficientType, typename _ComplexType>
std::vector<PairCell> Dual_HDVF<_CoefficientType,_ComplexType>::findPairsA(int q, bool &found) // const
{
    vector<PairCell> pairs;
    found = false ;
    
    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = this->_DD_col[q+1].begin(); it_col != this->_DD_col[q+1].end(); ++it_col)
    {
        const typename HDVF_type::CChain& col(OSM::cgetColumn(this->_DD_col[q+1], *it_col)) ;
        
        // Iterate through the entries of the column
        for (typename HDVF_type::CChain::const_iterator it = col.begin(); it != col.end(); ++it) {
            if (_subCC.get_Bit(q, it->first) && ((it->second == 1) || (it->second == -1))) {
                // If an entry of _subCC with coefficient 1 or -1 is found, set the pair and mark as found
                PairCell p;
                p.sigma = it->first;
                p.tau = *it_col;
                p.dim = q;
                pairs.push_back(p) ;
                found = true;
            }
        }
    }
    return pairs;
}

/** \brief find all the valid PairCell containing tau for A in dimension q */
template<typename _CoefficientType, typename _ComplexType>
std::vector<PairCell> Dual_HDVF<_CoefficientType,_ComplexType>::findPairsA(int q, bool &found, int tau) // const
{
    found = false;
    std::vector<PairCell> pairs;
    // Check if tau belongs to _subCC
    if (!_subCC.get_Bit(q, tau))
        throw("Dual_HDVF: searching for a cell tau outside _subCC") ;
    
    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    // and tau' belongs to _subCC
    typename HDVF_type::RChain tmp(OSM::getRow(this->_DD_col.at(q+1), tau)) ;
    for (typename HDVF_type::RChain::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
    {
        if (_subCC.get_Bit(q+1, it->first) && (abs(it->second) == 1))
        {
            found = true ;
            PairCell p ;
            p.sigma = tau ;
            p.tau = it->first ;
            p.dim = q ;
            pairs.push_back(p) ;
        }
    }
    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    // and tau' belongs to _subCC
    const typename HDVF_type::CChain& tmp2(OSM::cgetColumn(this->_DD_col.at(q), tau)) ;
    for (typename HDVF_type::CChain::const_iterator it = tmp2.cbegin(); it != tmp2.cend(); ++it)
    {
        if (_subCC.get_Bit(q-1, it->first) && (abs(it->second) == 1))
        {
            found = true ;
            PairCell p ;
            p.sigma = it->first ;
            p.tau = tau ;
            p.dim = q-1 ;
            pairs.push_back(p) ;
        }
    }
    return pairs;
}

template<typename _CoefficientType, typename _ComplexType>
void Dual_HDVF<_CoefficientType,_ComplexType>::computeDualPerfectHDVF()
{
    std::cout << "==== Compute perfect HDVF over K" << endl ;
    // Set _subCC to K
    _subCC = _KCC ;
    // Restrict _DD_col accordingly
    screen_DD_col() ;
    // Compute perfect HDVF over K
    this->computePerfectHDVF() ;
    critical_K = this->get_flag(CRITICAL) ;
//    this->print_matrices() ;
//    this->print_reduction() ;

    std::cout << "==== Compute perfect HDVF over L-K" << endl ;
    // set _subCC to L-K
    _subCC = _KCC.complement() ;
    // Restrict _DD_col accordingly
    screen_DD_col() ;
    // Compute perfect HDVF over L-K
    this->computePerfectHDVF() ;
    critical_L_K = this->get_flag(CRITICAL) ;
//    this->print_matrices() ;
//    this->print_reduction() ;
    
//    std::cout << "==== Pairing matrix" << endl ;
//
//    print_bnd_pairing() ;
}

template<typename _CoefficientType, typename _ComplexType>
std::vector<PairCell> Dual_HDVF<_CoefficientType,_ComplexType>::computePairingHDVF()
{
    // TODO : check both HDVFs are perfect
    
    std::cout << "==== Compute pairing" << endl ;

    // Create a full SubChainComplex
    _subCC = SubChainComplex<_CoefficientType, _ComplexType>(_L) ;
    screen_DD_col() ;
    // Copy the HDVF before computing the pairing -> otherwise we loose it...
    std::vector<PairCell> pairing = this->computePerfectHDVF() ;
    return pairing ;
}


// Method to get cells if with a given flag (P,S,C) for each dimension
template<typename CoefficientType, typename ComplexType>
vector<vector<int> > Dual_HDVF<CoefficientType,ComplexType>::get_flag (FlagType flag) const
{
    vector<vector<int> > res(_L.dim()+1) ;
    for (int q=0; q<=_L.dim(); ++q)
    {
        for (int i=0; i<_L.nb_cells(q); ++i)
        {
            if (_subCC.get_Bit(q, i) && (this->_flag.at(q).at(i) == flag))
                res.at(q).push_back(i) ;
        }
    }
    return res ;
}

// Method to get cells with a given flag (P,S,C) for a given dimension
template<typename CoefficientType, typename ComplexType>
vector<int> Dual_HDVF<CoefficientType,ComplexType>::get_flag_dim (FlagType flag, int q) const
{
    vector<int> res ;
    for (int i=0; i<this->_K.nb_cells(q); ++i)
    {
        if (_subCC.get_Bit(q, i) && (this->_flag.at(q).at(i) == flag))
            res.push_back(i) ;
    }
    return res ;
}

#endif /* HDVF_DUALITY_H */
