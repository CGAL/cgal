/**
 * \file hdvf.hpp
 * \brief HDVF code library.
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 *
 * HDVF code library (constructors, operations, finders).
 */

#ifndef HDVF_H
#define HDVF_H

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/OSM/Bitboard.hpp"
#include "CGAL/HDVF/sub_chain_complex.hpp"

namespace CGAL {

/** \brief Enum for the HDVF label of cells. */
enum FlagType {
    PRIMARY,
    SECONDARY,
    CRITICAL,
    NONE // For duality and persistence
};

/** \brief HDVF options (partial computation of the reduction). */
// Compute only reduced boundary
const int OPT_BND = 0b0001;
// Compute reduced boundary and f
const int OPT_F = 0b0010;
// Compute reduced boundary and g
const int OPT_G = 0b0100;
// Compute full reduction
const int OPT_FULL = 0b1000;

/** \brief Enum for export: lets chose the type of export.
 *
 * \tparam PSC export a flag encoding the label of the cell: P (-1), S (1), C (0)
 * \tparam FSTAR exports a boolean flag for a given cohomology generator (true: the cell belongs to the cycle, false: the cell does not belong)
 * \tparam G exports a boolean flag for a given homology generator (true: the cell belongs to the cycle, false: the cell does not belong)
 */
enum ExportType {
    PSC,
    FSTAR,
    G
};


/** \brief PairCell: Structure to represent a pair of cells (for HDVF operations).
 * The PairCell follows the following conventions depending on the operation type:
 * - A: sigma is of dim q and tau is of dim q+1 ;
 * - R: sigma is of dim q and tau is of dim q+1 ;
 * - M: sigma is primary and tau is critical of dim q ;
 * - W: sigma is secondary and tau is critical of dim q ;
 * - MW: sigma is primary and tau is secondary of dim q.
 */
struct PairCell {
    int sigma;  // Index of the first cell
    int tau;    // Index of the second cell
    int dim;    // Dimension of cells: q/q+1 for A and R, q for other operations
};


/**
 * \class HDVF
 * \brief Implementation of HDVF and associate operations.
 *
 * The HDVF class contains all functions to build HDVF, create perfect HDVF (A and associate finders) and delineated generators through HDVF operations (R, M, W, MW and associate finders).
 *
 * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam _ComplexType The type of complex  (default is SimpComplex)
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 */


template<typename _CoefficientType, typename _ComplexType, template <typename, int> typename _ChainType = OSM::Chain, template <typename, int> typename _SparseMatrixType = OSM::SparseMatrix>
class HDVF {
public:
    /** \brief Typedefs for coefficient type and different chain and matrix types. */
    typedef _ChainType<_CoefficientType, OSM::COLUMN> CChain;
    typedef _ChainType<_CoefficientType, OSM::ROW> RChain;
    typedef _SparseMatrixType<_CoefficientType, OSM::COLUMN> CMatrix;
    typedef _SparseMatrixType<_CoefficientType, OSM::ROW> RMatrix;
    
protected:
    // Properties
    /** \brief Flags for the cells. */
    std::vector<std::vector<FlagType>> _flag;
    /** \brief Number of PRIMARY cells. */
    std::vector<int> _nb_P;
    /** \brief Number of SECONDARY cells. */
    std::vector<int> _nb_S;
    /** \brief Number of CRITICAL cells. */
    std::vector<int> _nb_C;
    /** \brief Row matrices for f. */
    std::vector<RMatrix> _F_row;
    /** \brief Column matrices for g. */
    std::vector<CMatrix> _G_col;
    /** \brief Column matrices for h. */
    std::vector<CMatrix> _H_col;
    /** \brief Column matrices for reduced boundary. */
    std::vector<CMatrix> _DD_col;
    
    /** \brief Reference to the underlying complex. */
    const _ComplexType& _K;
    
    /** \brief HDVF options for computation (computation of partial reduction). */
    int _hdvf_opt;
    
public:
    // Constructor - default : full reduction computed, no sub-complex
    HDVF(const _ComplexType& K, int hdvf_opt = OPT_FULL) ;
    // Copy constructor
    HDVF(const HDVF& hdvf) : _flag(hdvf._flag), _nb_P(hdvf._nb_P), _nb_S(hdvf._nb_S), _nb_C(hdvf._nb_C), _F_row(hdvf._F_row), _G_col(hdvf._G_col), _H_col(hdvf._H_col), _DD_col(hdvf._DD_col), _K(hdvf._K), _hdvf_opt(hdvf._hdvf_opt) { }
    // Destructor
    ~HDVF() { }
    
    // findPair functions
    
    /** \brief find a valid PairCell for A in dimension q */
    virtual PairCell findPairA(int q, bool &found); // const;
    /** \brief find a valid PairCell containing tau for A in dimension q */
    virtual PairCell findPairA(int q, bool &found, int tau); // const;
    /** \brief find all the valid PairCell for A in dimension q */
    virtual std::vector<PairCell> findPairsA(int q, bool &found); // const;
    /** \brief find all the valid PairCell containing tau for A in dimension q */
    virtual std::vector<PairCell> findPairsA(int q, bool &found, int tau); // const;
    
    /** \brief find a valid PairCell for M in dimension q */
    PairCell findPairM(int q, bool &found); // const;
    /** \brief find a valid PairCell containing tau for M in dimension q */
    PairCell findPairM(int q, bool &found, int tau); // const;
    /** \brief find all the valid PairCell for M in dimension q */
    std::vector<PairCell> findPairsM(int q, bool &found); // const;
    /** \brief find all the valid PairCell containing tau for M in dimension q */
    std::vector<PairCell> findPairsM(int q, bool &found, int tau); // const;
    
    /** \brief find a valid PairCell for W in dimension q */
    PairCell findPairW(int q, bool &found); // const;
    /** \brief find a valid PairCell containing tau for W in dimension q */
    PairCell findPairW(int q, bool &found, int tau); // const;
    /** \brief find all the valid PairCell for W in dimension q */
    std::vector<PairCell> findPairsW(int q, bool &found); // const;
    /** \brief find all the valid PairCell containing tau for W in dimension q */
    std::vector<PairCell> findPairsW(int q, bool &found, int tau); // const;
    
    /** \brief find a valid PairCell for MW in dimension q */
    PairCell findPairMW(int q, bool &found); // const;
    /** \brief find a valid PairCell containing tau for MW in dimension q */
    PairCell findPairMW(int q, bool &found, int tau); // const;
    /** \brief find all the valid PairCell for MW in dimension q */
    std::vector<PairCell> findPairsMW(int q, bool &found); // const;
    /** \brief find all the valid PairCell containing tau for MW in dimension q */
    std::vector<PairCell> findPairsMW(int q, bool &found, int tau); // const;
    
    // HDVF methods
    /** \brief HDVF operation A(gamma1, gamma2) */
    void A(int gamma1, int gamma2, int dim);
    void R(int pi, int sigma, int dim);
    void M(int pi, int gamma, int dim);
    void W(int sigma, int gamma, int dim);
    void MW(int pi, int sigma, int dim);
    
    /** \brief Compute a perfect HDVF
     From dimension q-1 to dimension 0
     As long as findPairA(q) returns a pair
     -> Perform operation A
     -> Add the pair to the output vector
     */
    std::vector<PairCell> computePerfectHDVF(bool verbose = false);
    
    /** \brief Compute a perfect HDVF
     From dimension q-1 to dimension 0
     As long as findPairsA(q) returns pairs
     -> Choose one randomly
     -> Perform operation A
     -> Add the pair to the output vector
     Slower than computePerfectHDVF
     */
    std::vector<PairCell> computeRandPerfectHDVF(bool verbose = false);
    
    // HDVF getters
    // Method to get cells if with a given flag (P,S,C) for each dimension
    virtual std::vector<std::vector<int> > get_flag (FlagType flag) const ;
    // Method to get cells with a given flag (P,S,C) for a given dimension
    virtual std::vector<int> get_flag_dim (FlagType flag, int q) const ;
    // Method to get the flag (P,S,C) of a cell tau of dimension q
    FlagType get_cell_flag (int q, int tau) const { return _flag.at(q).at(tau); }
    
    // Method to get options
    int get_hdvf_opts () const { return _hdvf_opt ; }
    
    // Method to print the matrices _F_row, _G_col, _H_col
    virtual std::ostream& print_matrices(std::ostream &out = std::cout) const;
    
    // Method to print the reduction
    virtual std::ostream& print_reduction(std::ostream &out = std::cout); // const;
    
    // Method to print A-pairs
    virtual std::ostream& print_pairs(const std::vector<PairCell>& pairs, std::ostream &out = std::cout);
    
    // Method to project a chain onto P, S, or C cells
    template<int ChainTypeFlag>
    _ChainType<_CoefficientType, ChainTypeFlag> projection(const _ChainType<_CoefficientType, ChainTypeFlag>& chain, FlagType flag, int dim) const {
        // Create a new chain to store the result
        // Better to initialize 'result' directly with the correct size and iterate over it
        _ChainType<_CoefficientType, ChainTypeFlag> result(chain);
        
        // Iterate over each element of the chain
        std::vector<int> tmp ;
        for (typename _ChainType<_CoefficientType, ChainTypeFlag>::const_iterator it = result.cbegin(); it != result.cend(); ++it)
        {
            int cell_index = it->first;
            _CoefficientType value = it->second;
            
            // Check the flag of the corresponding cell
            if (_flag[dim][cell_index] != flag) {
                // Set to 0
                tmp.push_back(cell_index) ;
                
            }
        }
        result /= tmp ;
        return result;
    }
    
    // [OLD] Method to generate labels for visualisation
    // Kept for compatibility
    // For PSC, no additional arguments, for FSTAR and G, specify the critical cell concerned
    virtual std::vector<std::vector<int> > export_label (ExportType type, int cell=0, int dim=0) const
    {
        std::vector<std::vector<int> > labels(_K.dim()+1) ;
        if (type == PSC)
        {
            for (int q=0; q<=_K.dim(); ++q)
            {
                for (int i = 0; i<_K.nb_cells(q); ++i)
                {
                    if (_flag.at(q).at(i) == PRIMARY)
                        labels.at(q).push_back(-1) ;
                    else if (_flag.at(q).at(i) == SECONDARY)
                        labels.at(q).push_back(1) ;
                    else
                        labels.at(q).push_back(0) ;
                }
            }
        }
        else
        {
            if ((type == FSTAR) && (_hdvf_opt & (OPT_FULL | OPT_F)))
            {
                for (int q=0; q<=_K.dim(); ++q)
                {
                    labels.at(q).resize(_K.nb_cells(q)) ;
                }
                
                if (dim < _K.dim())
                {
                    labels.at(dim).at(cell) = 1 ;
                    const RChain& fstar_cell = OSM::cgetRow(_F_row.at(dim), cell) ;
                    for (typename RChain::const_iterator it = fstar_cell.cbegin(); it != fstar_cell.cend(); ++it)
                    {
                        // Set the cofaces of it->first in dimension dim+1
                        RChain cofaces(_K.cod(it->first,dim)) ;
                        for (typename RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                            labels.at(dim+1).at(it2->first) = 1 ; // set the flag to 1 only for FSTAR(cell) cells
                    }
                }
                else
                    throw "Error: cannot export f* for a cell of maximal dimension" ;
            }
            else if (_hdvf_opt & (OPT_FULL | OPT_G)) // G
            {
                for (int q=0; q<=_K.dim(); ++q)
                {
                    labels.at(q).resize(_K.nb_cells(q)) ;
                }
                
                labels.at(dim).at(cell) = 1 ;
                const CChain& g_cell = OSM::cgetColumn(_G_col.at(dim), cell) ;
                for (typename CChain::const_iterator it = g_cell.cbegin(); it != g_cell.cend(); ++it)
                    labels.at(dim).at(it->first) = 1 ; // set the flag to 1 only for FSTAR(cell) cells
            }
        }
        return labels ;
    }
    
    // Method to generate PSC labels for visualisation
    virtual std::vector<std::vector<int> > export_labelsPSC () const
    {
        std::vector<std::vector<int> > labels(_K.dim()+1) ;
        for (int q=0; q<=_K.dim(); ++q)
        {
            for (int i = 0; i<_K.nb_cells(q); ++i)
            {
                if (_flag.at(q).at(i) == PRIMARY)
                    labels.at(q).push_back(-1) ;
                else if (_flag.at(q).at(i) == SECONDARY)
                    labels.at(q).push_back(1) ;
                else if (_flag.at(q).at(i) == CRITICAL)
                    labels.at(q).push_back(0) ;
                else // NONE
                    labels.at(q).push_back(2) ;
            }
        }
        return labels ;
    }
    
    // Method to export G chains for visualisation
    virtual CChain export_GChain (int cell, int dim) const
    {
        if ((dim<0) || (dim>_K.dim()))
            throw "Error : export_GChain with dim out of range" ;
        if (_hdvf_opt & (OPT_FULL | OPT_G))
        {
            CChain g_cell(OSM::getColumn(_G_col.at(dim), cell)) ;
            // Add 1 to the cell
            g_cell[cell] = 1 ;
            return g_cell ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
    
    // Method to export FSTAR chains for visualisation
    virtual CChain export_FSTARChain (int cell, int dim) const
    {
        if ((dim<0) || (dim>_K.dim()))
            throw "Error : export_GChain with dim out of range" ;
        if (_hdvf_opt & (OPT_FULL | OPT_F))
        {
            RChain fstar_cell(OSM::getRow(_F_row.at(dim), cell)) ;
            // Add 1 to the cell
            fstar_cell[cell] = 1 ;
            // Compute the cofaces
            if (dim < _K.dim())
            {
                CChain fstar_cofaces(_K.nb_cells(dim+1)) ;
                for (typename RChain::const_iterator it = fstar_cell.cbegin(); it != fstar_cell.cend(); ++it)
                {
                    // Set the cofaces of it->first in dimension dim+1
                    RChain cofaces(_K.cod(it->first,dim)) ;
                    for (typename RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                        fstar_cofaces[it2->first] = 1 ;
                }
                return fstar_cofaces ;
            }
            else
                return CChain(0) ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
protected:
    void progress_bar(int i, int n)
    {
        const int step(n/20) ;
        if ((i%step)==0)
        {
            const float percentage(float(i)/(n-1)) ;
            const char PBSTR[] = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" ;
            const int PBWIDTH(60) ;
            int val = (int) (percentage * 100);
            int lpad = (int) (percentage * PBWIDTH);
            int rpad = PBWIDTH - lpad;
            printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
            fflush(stdout);
        }
        if (i==(n-1))
            std::cout << std::endl ;
    }
    
};

// Constructor for the HDVF class
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::HDVF(const ComplexType& K, int hdvf_opt) : _K(K) {
    // Get the dimension of the simplicial complex
    int dim = _K.dim();
    std::cout << "----> Starting HDVF creation / dim " << dim << std::endl ;
    // HDVF options
    _hdvf_opt = hdvf_opt ;
    
    // Resize the _DD_col vector to hold dim+1 elements
    _DD_col.resize(dim + 1);
    
    // Resize the _F_row vector to hold dim+1 elements
    if (_hdvf_opt & (OPT_FULL | OPT_F))
        _F_row.resize(dim + 1);
    
    // Resize the _G_col vector to hold dim+1 elements
    if (_hdvf_opt & (OPT_FULL | OPT_G))
        _G_col.resize(dim + 1);
    
    // Resize the _H_col vector to hold dim+1 elements
    if (_hdvf_opt & OPT_FULL)
        _H_col.resize(dim + 1);
    
    // Resize flag and count vectors to hold dim+1 elements
    _flag.resize(dim + 1);
    _nb_P.resize(dim + 1);
    _nb_S.resize(dim + 1);
    _nb_C.resize(dim + 1);
    
    // Initialize matrices and counters
    
    for (int q = 0; q <= dim; q++) {
        // Initialize _F_row[q] as a row matrix with dimensions (dim(q) x dim(q))
        if (_hdvf_opt & (OPT_FULL | OPT_F))
            _F_row[q] = RMatrix(_K.nb_cells(q), _K.nb_cells(q));
        
        // Initialize _G_col[q] as a column matrix with dimensions (dim(q) x dim(q))
        if (_hdvf_opt & (OPT_FULL | OPT_G))
            _G_col[q] = CMatrix(_K.nb_cells(q), _K.nb_cells(q));
        
        // Initialize _H_col[q] as a column matrix with dimensions (dim(q+1) x dim(q))
        if (_hdvf_opt & OPT_FULL)
            _H_col[q] = CMatrix(_K.nb_cells(q + 1), _K.nb_cells(q));
        
        // Initialize the counters for PRIMARY, SECONDARY, and CRITICAL cells
        _nb_P[q] = 0;
        _nb_S[q] = 0;
        _nb_C[q] = _K.nb_cells(q);
    }
    
    // Initialize the flags for each dimension to CRITICAL
    for (int q = 0; q <= dim; q++) {
        _flag[q] = std::vector<FlagType>(_K.nb_cells(q), CRITICAL);
    }
    
    // Populate the DD matrices
    _DD_col.resize(_K.dim()+1) ;
    for (int q=0; q<=_K.dim(); ++q)
        _DD_col.at(q) = _K.get_bnd_matrix(q) ;
    std::cout << "------> End HDVF creation" << std::endl ;
}

// Method to print the matrices
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::ostream& HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::print_matrices(std::ostream& out) const {
    // Iterate through each dimension and print the corresponding matrices
    for (int q = 0; q <= _K.dim(); ++q) {
        out << "------- Dimension " << q << std::endl;
        
        out << "Matrices _DD_col:" << std::endl;
        out << _DD_col[q] << std::endl;
        
        if (_hdvf_opt & (OPT_FULL | OPT_F))
        {
            out << "Matrices _F_row:" << std::endl;
            out << _F_row[q] << std::endl;
        }
        
        if (_hdvf_opt & (OPT_FULL | OPT_G))
        {
            out << "Matrices _G_col:" << std::endl;
            out << _G_col[q] << std::endl;
        }
        
        if (_hdvf_opt & OPT_FULL)
        {
            out << "Matrices _H_col:" << std::endl;
            out << _H_col[q] << std::endl;
        }
    }
    return out ;
}


// Methods to find a pair of cells for A in dimension q
//  -> <d(cell2),cell1> = +- 1 with cell1, cell2 critical
// First version: returns a pair of dimensions q / q+1
// Second version: returns all the pairs containing sigma

/** \brief find a valid PairCell for A in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairA(int q, bool &found) // const
{
    found = false;
    PairCell p;
    
    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = _DD_col[q+1].begin(); (it_col != _DD_col[q+1].end() && !found); ++it_col)
    {
        const CChain& col(OSM::cgetColumn(_DD_col[q+1], *it_col)) ;
        
        // Iterate through the entries of the column
        for (typename CChain::const_iterator it = col.begin(); (it != col.end() && !found); ++it) {
            if ((it->second == 1) || (it->second == -1)) {
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
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairA(int q, bool &found, int tau) // const
{
    found = false;
    PairCell p ;
    
    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    const CChain& tmp2(OSM::cgetColumn(_DD_col.at(q), tau)) ;
    for (typename CChain::const_iterator it = tmp2.cbegin(); (it != tmp2.cend() && !found); ++it)
    {
        if (abs(it->second) == 1)
        {
            found = true ;
            p.sigma = it->first ;
            p.tau = tau ;
            p.dim = q-1 ;
        }
    }
    
    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    RChain tmp(OSM::getRow(_DD_col.at(q+1), tau)) ;
    for (typename RChain::const_iterator it = tmp.cbegin(); (it != tmp.cend() && !found); ++it)
    {
        if (abs(it->second) == 1)
        {
            found = true ;
            p.sigma = tau ;
            p.tau = it->first ;
            p.dim = q ;
        }
    }
    return p;
}

/** \brief find all the valid PairCell for A in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsA(int q, bool &found) // const
{
    std::vector<PairCell> pairs;
    found = false ;
    
    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = this->_DD_col[q+1].begin(); it_col != this->_DD_col[q+1].end(); ++it_col)
    {
        const CChain& col(OSM::cgetColumn(this->_DD_col[q+1], *it_col)) ;
        
        // Iterate through the entries of the column
        for (typename CChain::const_iterator it = col.begin(); it != col.end(); ++it) {
            if ((it->second == 1) || (it->second == -1)) {
                // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
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
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsA(int q, bool &found, int tau) // const
{
    found = false;
    std::vector<PairCell> pairs;
    
    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    RChain tmp(OSM::getRow(_DD_col.at(q+1), tau)) ;
    for (typename RChain::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
    {
        if (abs(it->second) == 1)
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
    const CChain& tmp2(OSM::cgetColumn(_DD_col.at(q), tau)) ;
    for (typename CChain::const_iterator it = tmp2.cbegin(); it != tmp2.cend(); ++it)
    {
        if (abs(it->second) == 1)
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

// Methods to find a pair of cells for M in dimension q
//  -> <f(cell1),cell2> = +- 1 with cell1 primary, cell2 critical
// First version: returns a pair of dimensions q
// Second version: returns all the pairs containing sigma

/** \brief find a valid PairCell for M in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairM(int q, bool &found) // const
{
    found = false;
    PairCell p;
    
    if (_hdvf_opt & OPT_F)
    {
        // Search for +-1 in _F - iterate over rows
        for (OSM::Bitboard::iterator it_row = _F_row[q].begin(); (it_row != _F_row[q].end() && !found); ++it_row)
        {
            const RChain &row(OSM::cgetRow(_F_row[q],*it_row));
            
            // Iterate through the entries of the row
            for (typename RChain::const_iterator it = row.begin(); (it != row.end() && !found); ++it) {
                if ((it->second == 1) || (it->second == -1)) {
                    // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
                    p.sigma = it->first; // primary cell
                    p.tau = *it_row; // critical cell
                    p.dim = q;
                    found = true;
                }
            }
        }
    }
    return p;
}

/** \brief find a valid PairCell containing tau for M in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairM(int q, bool &found, int tau) // const
{
    found = false;
    PairCell p;
    if (_flag.at(q).at(tau) == SECONDARY)
        return p ; // Empty / found false
    
    if (_hdvf_opt & OPT_F)
    {
        // If tau is primary, search for gamma such that <f(tau),gamma>=+-1
        if (_flag.at(q).at(tau) == PRIMARY)
        {
            for (OSM::Bitboard::iterator it_row = _F_row.at(q).begin(); (it_row != _F_row.at(q).end() && !found); ++it_row)
            {
                if (abs(_F_row.at(q).get_coef(*it_row,tau)) == 1)
                {
                    p.sigma = tau ; // primary cell
                    p.tau = *it_row ; // critical cell
                    p.dim = q ;
                    found = true;
                }
            }
        }
        // If tau is critical, search for pi such that <f(pi),tau>=+-1
        if (_flag.at(q).at(tau) == CRITICAL)
        {
            const RChain& row(OSM::getRow(_F_row.at(q), tau)) ;
            for (typename RChain::const_iterator it = row.cbegin(); (it != row.cend() && !found); ++it)
            {
                if (abs(it->second) == 1)
                {
                    p.sigma = it->first ; // primary cell
                    p.tau = tau ; // critical cell
                    p.dim = q ;
                    found = true;
                }
            }
        }
    }
    return p;
}

/** \brief find all the valid PairCell for M in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsM(int q, bool &found) // const
{
    found = false;
    std::vector<PairCell> pairs;
    if (_hdvf_opt & OPT_F)
    {
        // Search for +-1 in _F - iterate over rows
        for (OSM::Bitboard::iterator it_row = _F_row[q].begin(); it_row != _F_row[q].end() ; ++it_row)
        {
            const RChain &row(OSM::cgetRow(_F_row[q],*it_row));
            
            // Iterate through the entries of the row
            for (typename RChain::const_iterator it = row.begin(); it != row.end(); ++it) {
                if ((it->second == 1) || (it->second == -1)) {
                    // If an entry with coefficient 1 or -1 is found, set the pair and add it
                    PairCell p ;
                    p.sigma = it->first; // primary cell
                    p.tau = *it_row; // critical cell
                    p.dim = q;
                    pairs.push_back(p) ;
                }
            }
        }
    }
    found = !pairs.empty() ;
    return pairs;
}

/** \brief find all the valid PairCell containing tau for M in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsM(int q, bool &found, int tau) // const
{
    found = false;
    std::vector<PairCell> pairs;
    if (_flag.at(q).at(tau) == SECONDARY)
        return pairs ; // Empty / found false
    
    if (_hdvf_opt & OPT_F)
    {
        // If tau is primary, search for gamma such that <f(tau),gamma>=+-1
        if (_flag.at(q).at(tau) == PRIMARY)
        {
            for (OSM::Bitboard::iterator it_row = _F_row.at(q).begin(); it_row != _F_row.at(q).end(); ++it_row)
            {
                if (abs(_F_row.at(q).get_coef(*it_row,tau)) == 1)
                {
                    PairCell p ;
                    p.sigma = tau ; // primary cell
                    p.tau = *it_row ; // critical cell
                    p.dim = q ;
                    pairs.push_back(p) ;
                }
            }
        }
        // If tau is critical, search for pi such that <f(pi),tau>=+-1
        if (_flag.at(q).at(tau) == CRITICAL)
        {
            const RChain& row(OSM::getRow(_F_row.at(q), tau)) ;
            for (typename RChain::const_iterator it = row.cbegin(); it != row.cend(); ++it)
            {
                if (abs(it->second) == 1)
                {
                    PairCell p ;
                    p.sigma = it->first ; // primary cell
                    p.tau = tau ; // critical cell
                    p.dim = q ;
                    pairs.push_back(p) ;
                }
            }
        }
    }
    found = !pairs.empty() ;
    return pairs;
}

// Methods to find a pair of cells for W in dimension q
//  -> <g(cell2),cell1> = +- 1 with cell1 secondary, cell2 critical
// First version: returns a pair of dimensions q
// Second version: returns all the pairs containing sigma

/** \brief find a valid PairCell for W in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairW(int q, bool &found) // const
{
    found = false;
    PairCell p;
    
    if (_hdvf_opt & OPT_G)
    {
        // Search for +-1 in _G - iterate over cols
        for (OSM::Bitboard::iterator it_col = _G_col[q].begin(); (it_col != _F_row[q].end() && !found); ++it_col)
        {
            CChain &col = _G_col[q][*it_col];
            
            // Iterate through the entries of the col
            for (typename CChain::const_iterator it = col.begin(); (it != col.end() && !found); ++it) {
                if ((it->second == 1) || (it->second == -1)) {
                    // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
                    p.sigma = it->first; // secondary cell
                    p.tau = *it_col; // critical cell
                    p.dim = q;
                    found = true;
                }
            }
        }
    }
    return p;
}

/** \brief find a valid PairCell containing tau for W in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairW(int q, bool &found, int tau) // const
{
    found = false;
    PairCell p;
    if (_flag.at(q).at(tau) == PRIMARY)
        return p ; // Empty / found false
    
    if (_hdvf_opt & OPT_G)
    {
        // If tau is primary, search for gamma such that <g(gamma),tau>=+-1
        if (_flag.at(q).at(tau) == SECONDARY)
        {
            for (OSM::Bitboard::iterator it_col = _G_col.at(q).begin(); (it_col != _G_col.at(q).end() && !found); ++it_col)
            {
                if (abs(_G_col.at(q).get_coef(tau, *it_col)) == 1)
                {
                    p.sigma = tau ; // secondary cell
                    p.tau = *it_col ; // critical cell
                    p.dim = q ;
                    found = true;
                }
            }
        }
        // If tau is critical, search for sigma such that <g(tau),sigma>=+-1
        if (_flag.at(q).at(tau) == CRITICAL)
        {
            CChain col(OSM::getColumn(_G_col.at(q), tau)) ;
            for (typename CChain::const_iterator it = col.cbegin(); (it != col.cend() && !found); ++it)
            {
                if (abs(it->second) == 1)
                {
                    p.sigma = it->first ; // secondary cell
                    p.tau = tau ; // critical cell
                    p.dim = q ;
                    found = true;
                }
            }
        }
    }
    return p;
}

/** \brief find all the valid PairCell for W in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsW(int q, bool &found) // const
{
    found = false;
    std::vector<PairCell> pairs;
    
    if (_hdvf_opt & OPT_FULL)
    {
        // Search for +-1 in _G - iterate over cols
        for (OSM::Bitboard::iterator it_col = _G_col[q].begin(); it_col != _F_row[q].end(); ++it_col)
        {
            CChain &col = _G_col[q][*it_col];
            
            // Iterate through the entries of the col
            for (typename CChain::const_iterator it = col.begin(); it != col.end(); ++it) {
                if ((it->second == 1) || (it->second == -1)) {
                    // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
                    PairCell p ;
                    p.sigma = it->first; // secondary cell
                    p.tau = *it_col; // critical cell
                    p.dim = q;
                    pairs.push_back(p) ;
                }
            }
        }
    }
    found = !pairs.empty() ;
    return pairs;
}

/** \brief find all the valid PairCell containing tau for W in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsW(int q, bool &found, int tau) // const
{
    found = false;
    std::vector<PairCell> pairs;
    if (_flag.at(q).at(tau) == PRIMARY)
        return pairs ; // Empty / found false
    
    if (_hdvf_opt & OPT_FULL)
    {
        // If tau is primary, search for gamma such that <g(gamma),tau>=+-1
        if (_flag.at(q).at(tau) == SECONDARY)
        {
            for (OSM::Bitboard::iterator it_col = _G_col.at(q).begin(); it_col != _G_col.at(q).end(); ++it_col)
            {
                if (abs(_G_col.at(q).get_coef(tau, *it_col)) == 1)
                {
                    PairCell p ;
                    p.sigma = tau ; // secondary cell
                    p.tau = *it_col ; // critical cell
                    p.dim = q ;
                    pairs.push_back(p) ;
                }
            }
        }
        // If tau is critical, search for sigma such that <g(tau),sigma>=+-1
        if (_flag.at(q).at(tau) == CRITICAL)
        {
            CChain col(OSM::getColumn(_G_col.at(q), tau)) ;
            for (typename CChain::const_iterator it = col.cbegin(); it != col.cend(); ++it)
            {
                if (abs(it->second) == 1)
                {
                    PairCell p ;
                    p.sigma = it->first ; // secondary cell
                    p.tau = tau ; // critical cell
                    p.dim = q ;
                    pairs.push_back(p) ;
                }
            }
        }
    }
    found = !pairs.empty() ;
    return pairs;
}

// Method to find a pair of cells for MW in dimension q
//  -> <hd(cell1),cell2> = +- 1 and <dh(cell1),cell2> = +- 1 with cell1 primary, cell2 secondary
// First version: returns a pair of dimensions q
// Second version: returns all the pairs containing sigma

/** \brief find a valid PairCell for MW in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairMW(int q, bool &found) // const
{
    found = false;
    PairCell p;
    
    if (_hdvf_opt & OPT_FULL)
    {
        // pi and sigma must at least satisfy that col pi of H_q and row sigma of H_q-1 are non zero
        // iterate on H_q and H_q-1 accordingly
        for (OSM::Bitboard::iterator it_pi = _H_col[q].begin(); (it_pi != _H_col[q].end() && !found); ++it_pi)
        {
            const int pi = *it_pi ;
            CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
            CChain d_pi = _K.d(pi, q) ;
            CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
            for (int sigma = 0; (sigma < _H_col.at(q-1).dimensions().first && !found); ++sigma)
            {
                RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
                if (!H11q1.isNull())
                {
                    RChain cod_sigma = _K.cod(sigma, q) ;
                    RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
                    
                    // Compute xi and xi' to test the validity of MW
                    
                    const int xi = projS_cod_sigma * H11 ;
                    const int xip = H11q1 * projP_d_pi ;
                    found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                    if (found)
                    {
                        p.sigma = pi ; // primary cell
                        p.tau = sigma ; // secondary cell
                        p.dim = q ;
                    }
                }
            }
        }
    }
    return p;
}

/** \brief find a valid PairCell containing tau for MW in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
PairCell HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairMW(int q, bool &found, int tau) // const
{
    found = false;
    PairCell p;
    
    if (_flag.at(q).at(tau) == CRITICAL)
        return p ; // Empty / found false
    
    if (_hdvf_opt & OPT_FULL)
    {
        // If tau is primary (rename pi), search for a valid sigma
        if (_flag.at(q).at(tau) == PRIMARY)
        {
            const int pi(tau) ;
            // Col pi of H_q and proj_P(d(pi)) must at least be non empty
            CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
            CChain d_pi = _K.d(pi, q) ;
            CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
            
            if (H11.isNull() || projP_d_pi.isNull())
                return p ;
            
            // Search for sigma with col_sigma(H_q-1) non empty
            for (int sigma = 0; (sigma < _H_col.at(q-1).dimensions().first && !found); ++sigma)
            {
                RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
                if (!H11q1.isNull())
                {
                    // and proj_S(cod(sigma)) non empty
                    RChain cod_sigma = _K.cod(sigma, q) ;
                    RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
                    if (!projS_cod_sigma.isNull())
                    {
                        // test xi and xip
                        const int xi = projS_cod_sigma * H11 ;
                        const int xip = H11q1 * projP_d_pi ;
                        found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                        if (found)
                        {
                            p.sigma = pi ; // primary cell
                            p.tau = sigma ; // critical cell
                            p.dim = q ;
                        }
                    }
                }
            }
        }
        else // cell is secondary
        {
            // cout << "secondary" << endl ;
            const int sigma(tau) ;
            // Row sigma of H_q-1 and proj_S(cod(sigma)) must at least be non empty
            RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
            RChain cod_sigma = _K.cod(sigma, q) ;
            RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
            if (H11q1.isNull() || projS_cod_sigma.isNull())
                return p ;
            
            // Search for pi with col pi of H_q non empty
            for (OSM::Bitboard::iterator it_pi = _H_col[q].begin(); (it_pi != _H_col[q].end() && !found); ++it_pi)
            {
                const int pi(*it_pi) ;
                
                CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
                CChain d_pi = _K.d(pi, q) ;
                CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
                // proj_P of d(pi) must also be non empty
                if (!projP_d_pi.isNull())
                {
                    // test xi and xip
                    const int xi = projS_cod_sigma * H11 ;
                    const int xip = H11q1 * projP_d_pi ;
                    found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                    if (found)
                    {
                        p.sigma = pi ; // primary cell
                        p.tau = sigma ; // secondary cell
                        p.dim = q ;
                    }
                }
            }
        }
    }
    return p;
}

/** \brief find all the valid PairCell for MW in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsMW(int q, bool &found) // const
{
    found = false;
    std::vector<PairCell> pairs;
    
    if (_hdvf_opt & OPT_FULL)
    {
        // pi and sigma must at least satisfy that col pi of H_q and row sigma of H_q-1 are non zero
        // iterate on H_q and H_q-1 accordingly
        for (OSM::Bitboard::iterator it_pi = _H_col[q].begin(); it_pi != _H_col[q].end(); ++it_pi)
        {
            const int pi = *it_pi ;
            CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
            CChain d_pi = _K.d(pi, q) ;
            CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
            for (int sigma = 0; sigma < _H_col.at(q-1).dimensions().first; ++sigma)
            {
                RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
                if (!H11q1.isNull())
                {
                    RChain cod_sigma = _K.cod(sigma, q) ;
                    RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
                    
                    // Compute xi and xi' to test the validity of MW
                    
                    const int xi = projS_cod_sigma * H11 ;
                    const int xip = H11q1 * projP_d_pi ;
                    found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                    if (found)
                    {
                        PairCell p;
                        p.sigma = pi ; // primary cell
                        p.tau = sigma ; // secondary cell
                        p.dim = q ;
                        pairs.push_back(p) ;
                    }
                }
            }
        }
    }
    found = !pairs.empty() ;
    return pairs;
}

/** \brief find all the valid PairCell containing tau for MW in dimension q */
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::findPairsMW(int q, bool &found, int tau) // const
{
    found = false;
    std::vector<PairCell> pairs;
    if (_flag.at(q).at(tau) == CRITICAL)
        return pairs ; // Empty / found false
    
    if (_hdvf_opt & OPT_FULL)
    {
        // If tau is primary (rename pi), search for a valid sigma
        if (_flag.at(q).at(tau) == PRIMARY)
        {
            const int pi(tau) ;
            // Col pi of H_q and proj_P(d(pi)) must at least be non empty
            CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
            CChain d_pi = _K.d(pi, q) ;
            CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
            
            if (H11.isNull() || projP_d_pi.isNull())
                return pairs ;
            
            // Search for sigma with col_sigma(H_q-1) non empty
            for (int sigma = 0; sigma < _H_col.at(q-1).dimensions().first; ++sigma)
            {
                RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
                if (!H11q1.isNull())
                {
                    // and proj_S(cod(sigma)) non empty
                    RChain cod_sigma = _K.cod(sigma, q) ;
                    RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
                    if (!projS_cod_sigma.isNull())
                    {
                        // test xi and xip
                        const CoefficientType xi = projS_cod_sigma * H11 ;
                        const CoefficientType xip = H11q1 * projP_d_pi ;
                        found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                        if (found)
                        {
                            PairCell p ;
                            p.sigma = pi ; // primary cell
                            p.tau = sigma ; // critical cell
                            p.dim = q ;
                            pairs.push_back(p) ;
                        }
                    }
                }
            }
        }
        else // cell is secondary
        {
            // cout << "secondary" << endl ;
            const int sigma(tau) ;
            // Row sigma of H_q-1 and proj_S(cod(sigma)) must at least be non empty
            RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
            RChain cod_sigma = _K.cod(sigma, q) ;
            RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
            if (H11q1.isNull() || projS_cod_sigma.isNull())
                return pairs ;
            
            // Search for pi with col pi of H_q non empty
            for (OSM::Bitboard::iterator it_pi = _H_col[q].begin(); it_pi != _H_col[q].end(); ++it_pi)
            {
                const int pi(*it_pi) ;
                
                CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
                CChain d_pi = _K.d(pi, q) ;
                CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
                // proj_P of d(pi) must also be non empty
                if (!projP_d_pi.isNull())
                {
                    // test xi and xip
                    const CoefficientType xi = projS_cod_sigma * H11 ;
                    const CoefficientType xip = H11q1 * projP_d_pi ;
                    found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                    if (found)
                    {
                        PairCell p ;
                        p.sigma = pi ; // primary cell
                        p.tau = sigma ; // secondary cell
                        p.dim = q ;
                        pairs.push_back(p) ;
                    }
                }
            }
        }
    }
    found = !pairs.empty() ;
    return pairs;
}


// Method to perform operation A
// tau1 is in dimension q, tau2 is in dimension q+1
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
void HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::A(int tau1, int tau2, int q) {
    //----------------------------------------------- Submatrices of D ----------------------------------------------------
    
    // Output operation details to the console
    //    cout << "A of " << tau1 << "(dim " << q << ") / " << tau2 << "(dim " << q + 1 << ")" << endl;
    
    // Extract submatrices from _DD_col
    RChain D12(OSM::getRow(_DD_col.at(q+1),tau1)); // D12 is a row chain from _DD_col[q+1] at index tau1
    CChain D21(OSM::getColumn(_DD_col.at(q + 1),tau2)); // D21 is a column chain from _DD_col[q+1] at index tau2
    CoefficientType D11(D12[tau2]); // D11 is the coefficient at the intersection of tau2 in D12
    
    // Assert that D11 is either 1 or -1 (check invertibility)
    assert((D11 == 1) || (D11 == -1)); // !!!!! Test invertibility
    
    // Compute the inverse of D11 (which is itself, since D11 is 1 or -1)
    CoefficientType D11_inv = D11;
    
    // Perform operations to remove the row and column contributions
    D12 /= std::vector<int>({tau2}); // Remove tau2 column from D12
    D21 /= std::vector<int>({tau1}); // Remove tau1 row from D21
    
    // Delete rows and columns from _DD_col
    _DD_col[q + 1].delRow(tau1); // Remove row tau1 from _DD_col[q+1]
    _DD_col[q + 1].delColumn(tau2); // Remove column tau2 from _DD_col[q+1]
    
    //---------------------------------------------- Submatrices of F -----------------------------------------------------
    
    RChain F11 ;
    CChain G11 ;
    if (_hdvf_opt & (OPT_FULL | OPT_F))
    {
        // Extract the relevant submatrix from _F_row
        F11 = OSM::getRow(_F_row.at(q),tau1); // F11 is a row chain from _F_row[q] at index tau1
        
        // Delete the row tau1 from _F_row
        _F_row[q].delRow(tau1);
    }
    
    //--------------------------------------------- Submatrices of G ------------------------------------------------------
    
    if (_hdvf_opt & (OPT_FULL | OPT_G))
    {
        // Extract the relevant submatrix from _G_col
        G11 = OSM::getColumn(_G_col.at(q + 1),tau2); // G11 is a column chain from _G_col[q+1] at index tau2
        
        // Delete the column tau2 from _G_col
        _G_col[q + 1].delColumn(tau2);
    }
    
    //--------------------------------------------- Update matrices -------------------------------------------------------
    
    // ---- Update _F_row
    
    if (_hdvf_opt & (OPT_FULL | OPT_F))
    {
        // Update _F_row[q]
        // Subtract the product of (D21 * D11_inv) and F11 from _F_row[q]
        // Note: % operator returns a row matrix, so be careful with operations
        _F_row[q] -= (D21 * D11_inv) % F11;
        
        // Set the column tau1 of _F_row[q] to (D21 * (-D11_inv))
        OSM::setColumn(_F_row[q], tau1, D21 * (-D11_inv));
        
        // Remove the row tau2 from _F_row[q+1]
        _F_row[q + 1].delRow(tau2);
    }
    
    // ---- Update _G_col
    
    if (_hdvf_opt & (OPT_FULL | OPT_G))
    {
        // Update _G_col[q + 1]
        // Subtract the product of (G11 * D11_inv) and D12 from _G_col[q+1]
        _G_col[q + 1] -= (G11 * D11_inv) * D12;
        
        // Set the row tau2 of _G_col[q + 1] to (D12 * (-D11_inv))
        OSM::setRow(_G_col[q + 1], tau2, D12 * (-D11_inv));
        
        // Remove the column tau1 from _G_col[q]
        _G_col[q].delColumn(tau1);
    }
    
    // ---- Update _H_col
    
    if (_hdvf_opt & OPT_FULL)
    {
        // Update _H_col[q]
        // Compute the temporary matrix product G11 * F11
        CMatrix tmp = G11 * F11;
        
        // Add the product of tmp and D11_inv to _H_col[q]
        _H_col[q] += (tmp) * D11_inv;
        
        // Set the row tau2 of _H_col[q] to (F11 * D11_inv)
        OSM::setRow(_H_col[q], tau2, F11 * D11_inv);
        
        // Set the column tau1 of _H_col[q] to (G11 * D11_inv)
        OSM::setColumn(_H_col[q], tau1, G11 * D11_inv);
        
        // Set the coefficient at (tau2, tau1) in _H_col[q] to D11_inv
        _H_col[q].set_coef(tau2, tau1, D11_inv);
    }
    
    // ---- Update _DD_col
    
    // Update _DD_col
    _DD_col[q + 1] -= (D21 * D12) * D11_inv;
    
    // Remove columns and rows from _DD_col as necessary
    if (q > 0) {
        _DD_col[q].delColumn(tau1); // Remove column tau1 from _DD_col[q]
    }
    if (q + 2 <= _K.dim()) {
        _DD_col[q + 2].delRow(tau2); // Remove row tau2 from _DD_col[q+2]
    }
    
    // Update flags
    _flag[q][tau1] = PRIMARY; // Set the flag of tau1 in dimension q to PRIMARY
    --_nb_C.at(q) ;
    ++_nb_P.at(q) ;
    _flag[q + 1][tau2] = SECONDARY; // Set the flag of tau2 in dimension q+1 to SECONDARY
    --_nb_C.at(q+1) ;
    ++_nb_S.at(q+1) ;
    
    // -----------------------------------------------------------------------------------------------------------------------
}

// Method to perform operation R
// pi is in dimension q, sigma is in dimension q+1
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
void HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::R(int pi, int sigma, int q) {
    //----------------------------------------------- Submatrices of H ----------------------------------------------------
    
    if (_hdvf_opt & OPT_FULL)
    {
        // Output operation details to the console
        std::cout << "R of " << pi << "(dim " << q << ") / " << sigma << "(dim " << q + 1 << ")" << std::endl;
        
        // Extract the relevant row and column chains from _H_col
        RChain H12 = OSM::getRow(_H_col[q], sigma); // H12 is the row chain from _H_col[q] at index sigma
        CChain H21 = OSM::getColumn(_H_col[q], pi); // H21 is the column chain from _H_col[q] at index pi
        
        // Get the coefficient at the intersection of H12 and H21
        CoefficientType H11(H12[pi]); // H11 is the coefficient at row sigma and column pi
        
        // Assert that H11 is either 1 or -1 (check invertibility)
        assert((H11 == 1) || (H11 == -1)); // !!!!! Test invertibility
        CoefficientType H11_inv = H11; // Inverse of H11 (which is itself for 1 or -1)
        
        // Remove the contributions of pi from H12 and sigma from H21
        H12 /= std::vector<int>({pi}); // Remove column pi from H12
        H21 /= std::vector<int>({sigma}); // Remove row sigma from H21
        
        // Remove the corresponding row and column from _H_col
        _H_col[q].delRow(sigma); // Remove row sigma from _H_col[q]
        _H_col[q].delColumn(pi); // Remove column pi from _H_col[q]
        
        //---------------------------------------------- Submatrices of F -----------------------------------------------------
        
        // Extract the relevant column chain from _F_row
        CChain F11 = OSM::getColumn(_F_row[q], pi); // F11 is the column chain from _F_row[q] at index pi
        
        // Remove the column pi from _F_row
        _F_row[q].delColumn(pi);
        
        //--------------------------------------------- Submatrices of G ------------------------------------------------------
        
        // Extract the relevant row chain from _G_col
        RChain G11 = OSM::getRow(_G_col[q + 1], sigma); // G11 is the row chain from _G_col[q+1] at index sigma
        
        // Remove the row sigma from _G_col
        _G_col[q + 1].delRow(sigma);
        
        //--------------------------------------------- Update matrices -------------------------------------------------------
        
        // Update _F_row[q]
        // Subtract the product of (H12 * F11) and H11 from _F_row[q]
        // Note: % operator returns a row matrix, so be careful with operations
        _F_row[q] -= (F11 % H12) * H11_inv;
        
        // Set the row pi of _F_row[q] to (H12 * (-H11_inv))
        OSM::setRow(_F_row[q], pi, H12 * (-H11_inv));
        
        // Update _G_col[q + 1]
        // Subtract the product of (H21 * G11) and H11_inv from _G_col[q + 1]
        _G_col[q + 1] -= (H21 * G11) * H11_inv;
        
        // Set the column sigma of _G_col[q + 1] to (H21 * (-H11_inv))
        OSM::setColumn(_G_col[q + 1], sigma, H21 * (-H11_inv));
        
        // Update _DD_col
        // Compute the temporary matrix product F11 * G11
        _DD_col[q + 1] += (F11 * G11) * H11_inv;
        
        // Set the row pi in _DD_col[q + 1] to (G11 * H11_inv)
        OSM::setRow(_DD_col[q + 1], pi, G11 * H11_inv);
        
        // Set the column sigma in _DD_col[q + 1] to (F11 * H11_inv)
        OSM::setColumn(_DD_col[q + 1], sigma, F11 * H11_inv);
        
        // Set the coefficient at (pi, sigma) in _DD_col to H11_inv
        _DD_col[q + 1].set_coef(pi, sigma, H11_inv);
        
        // Update _H_col[q]
        // Subtract the product of (H21 * H12) and H11_inv from _H_col[q]
        _H_col[q] -= (H21 * H12) * H11_inv;
        
        // Perform additional updates
        
        // Extract boundary and coboundary chains
        CChain bnd_pi(_K.d(pi, q)); // Boundary of pi in dimension q
        RChain cobnd_sigma(_K.cod(sigma, q + 1)); // Coboundary of sigma in dimension q+1
        
        // Project the boundary and coboundary chains onto PRIMARY and SECONDARY flags
        CChain proj_P_pi(projection(bnd_pi, PRIMARY, q));
        RChain proj_S_sigma(projection(cobnd_sigma, SECONDARY, q + 1));
        
        if (q > 0) {
            // Update _DD_col[q] with projections and _F_row[q-1]
            CChain c1(_F_row[q - 1] * proj_P_pi + projection(bnd_pi, CRITICAL, q));
            OSM::setColumn(_DD_col[q], pi, c1);
            OSM::setColumn(_G_col[q], pi, _H_col[q - 1] * proj_P_pi);
        }
        
        // Update _F_row[q+1] with projection and _H_col[q+1]
        OSM::setRow(_F_row[q + 1], sigma, proj_S_sigma * _H_col[q + 1]);
        
        if (q + 2 <= _K.dim()) {
            // Update _DD_col[q+2] with projections
            RChain c4(projection(cobnd_sigma, CRITICAL, q + 1) + proj_S_sigma * _G_col[q + 2]);
            OSM::setRow(_DD_col[q + 2], sigma, c4);
        }
        
        // Update flags
        _flag[q][pi] = CRITICAL; // Set the flag of pi in dimension q to CRITICAL
        ++_nb_C.at(q) ;
        --_nb_P.at(q) ;
        _flag[q + 1][sigma] = CRITICAL; // Set the flag of sigma in dimension q+1 to CRITICAL
        ++_nb_C.at(q+1) ;
        --_nb_S.at(q+1) ;
    }
    else
        std::cout << "!!! R impossible with partial reduction options" << std::endl;
}

// Method to perform operation M
// pi is in dimension q, gamma is in dimension q
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
void HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::M(int pi, int gamma, int q) {
    //----------------------------------------------- Submatrices of F ----------------------------------------------------
    
    if (_hdvf_opt & OPT_FULL)
    {
        // Output the operation details to the console
        std::cout << "M_" << q << "(" << pi << "," << gamma << ")" << std::endl;
        
        if (q == _K.dim())
            throw("M operation in max dimension !!!") ;
        
        // Extract row and column chains from _F_row
        RChain F12(OSM::getRow(_F_row[q], gamma)); // F12 is the row chain from _F_row[q] at index gamma
        CChain F21(OSM::getColumn(_F_row[q], pi)); // F21 is the column chain from _F_row[q] at index pi
        
        // Get the coefficient at the intersection of F12 and F21
        CoefficientType F11(F12[pi]); // F11 is the coefficient at row gamma and column pi
        
        // Assert that F11 is either 1 or -1 (for invertibility)
        assert((F11 == 1) || (F11 == -1)); // !!!!! Test invertibility
        CoefficientType F11_inv = F11; // Inverse of F11 (which is itself for 1 or -1)
        
        // Remove the contributions of pi from F12 and gamma from F21
        F12 /= std::vector<int>({pi}); // Remove column pi from F12
        F21 /= std::vector<int>({gamma}); // Remove row gamma from F21
        
        // Remove the corresponding row and column from _F_row
        _F_row[q].delRow(gamma); // Remove row gamma from _F_row[q]
        _F_row[q].delColumn(pi); // Remove column pi from _F_row[q]
        
        //--------------------------------------------- Submatrices of G ------------------------------------------------------
        
        // Extract the relevant column chain from _G_col
        // G11_q is the column chain from _G_col[q] at index gamma
        //        CChain G11_q(OSM::getColumn(_G_col[q], gamma));
        _G_col[q].delColumn(gamma); // Remove column gamma from _G_col[q]
        
        //---------------------------------------------- Submatrices of H -----------------------------------------------------
        
        // Extract the relevant column chain from _H_col
        CChain H11(OSM::getColumn(_H_col[q], pi)); // H11 is the column chain from _H_col[q] at index pi
        
        // Remove the column pi from _H_col
        _H_col[q].delColumn(pi);
        
        //--------------------------------------------- Submatrices of DD_q+1 ------------------------------------------------------
        
        // For DD_q+1 and DD_q:
        // Extract the relevant row chains from _DD_col
        
        // DD_q+1 (corresponds to the row matrix of _DD_col)
        RChain D11(OSM::getRow(_DD_col[q+1], gamma)); // D11 is the row chain from _DD_col[q+1] at index gamma
        _DD_col[q + 1].delRow(gamma); // Remove row gamma from _DD_col[q + 1]
        
        //--------------------------------------------- Submatrices of DD ------------------------------------------------------
        
        // DD_q (corresponds to the column matrix of _DD_col)
        // CChain D11_q(OSM::getColumn(_DD_col[q], gamma));
        // D11_q is the column chain from _DD_col[q] at index gamma
        if (q > 0)
            _DD_col[q].delColumn(gamma); // Remove column gamma from _DD_col[q]
        
        //--------------------------------------------- Update matrices -------------------------------------------------------
        
        // Update _H_col
        // Subtract the product of (H11 * F11_inv) and F12 from _H_col[q]
        // Note: * operator is used for matrix column results
        _H_col[q] -= (H11 * F11_inv) * F12;
        
        // Set the column gamma of _H_col[q] to (-1 * H11) * F11_inv
        OSM::setColumn(_H_col[q], gamma, (-1 * H11) * (F11_inv));
        
        // Update _F_row
        // Add the product of (F21 * F11_inv) and F12 to _F_row[q]
        _F_row[q] -= (F21 * F11_inv) * F12;
        
        // Set the row pi of _F_row[q] to F11_inv * F12
        OSM::setRow(_F_row[q], pi, F11_inv * F12);
        
        // Set the column gamma of _F_row[q] to F21 * (-F11_inv)
        OSM::setColumn(_F_row[q], gamma, F21 * (-F11_inv));
        
        // Set the coefficient at (pi, gamma) in _F_row[q] to F11_inv
        _F_row[q].set_coef(pi, gamma, F11_inv);
        
        // Update _G_col
        // Add the product of (H11 * F11_inv) and D11 to _G_col[q+1]
        _G_col[q + 1] += (H11 * F11_inv) * D11;
        
        // Update _DD_col[q + 1]
        _DD_col[q + 1] -= (F21 * F11_inv) * D11; // Subtract the product from _DD_col[q + 1]
        OSM::setRow(_DD_col[q + 1], pi, D11 * F11_inv); // Set the row pi in _DD_col[q + 1]
        
        // Update _G_col (for dimension q) and _DD_col (for dimension q)
        if (q>0)
        {
            // Extract boundary chain and project it
            CChain c = _K.d(pi, q); // Boundary of pi in dimension q
            CChain projection_p(projection(c, PRIMARY, q-1)); // Project boundary chain to PRIMARY
            CChain projection_c = projection(c, CRITICAL, q-1); // Project boundary chain to CRITICAL
            
            // Set the column pi of _G_col[q] to (-1 * _H_col[q-1]) * projection_p
            OSM::setColumn(_G_col[q], pi, (CoefficientType(-1) * _H_col[q - 1]) * projection_p);
            
            // Update _DD_col
            // Extract projections and perform updates
            CChain tmp(_F_row[q - 1] * projection_p + projection_c); // Compute the product of _F_row[q - 1] and projection_p
            // Set the column pi of _DD_col[q] to cc
            OSM::setColumn(_DD_col[q], pi, tmp);
        }
        // else: if the dimension is 0, no change for _G_col[q] and _DD_col[q]
        
        // Update flags
        _flag[q][pi] = CRITICAL; // Set the flag of pi in dimension q to PRIMARY
        _flag[q][gamma] = PRIMARY; // Set the flag of gamma in dimension q to CRITICAL
    }
    else
        std::cout << "!!! M impossible with partial reduction options" << std::endl;
}

// Method to perform operation W
// gamma is in dimension q, sigma is in dimension q

template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
void HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::W(int sigma, int gamma, int q) {
    
    //----------------------------------------------- Submatrices of G ----------------------------------------------------
    
    if (_hdvf_opt & OPT_FULL)
    {
        // Output the operation details to the console
        std::cout << "W_" << q << "(" << sigma << "," << gamma << ")" << std::endl;
        
        if (q == 0)
            throw("W operation in dimension 0 !!!") ;
        
        // Extract row and column chains from _G_col
        RChain G12(OSM::getRow(_G_col[q], sigma)); // G12 is the row chain from _G_col[q] at index sigma
        CChain G21(OSM::getColumn(_G_col[q], gamma)); // G21 is the column chain from _G_col[q] at index gamma
        
        // Get the coefficient at the intersection of G12 and G21
        CoefficientType G11(G12[gamma]); // G11 is the coefficient at row sigma and column gamma
        
        // Assert that G11 is either 1 or -1 (for invertibility)
        assert((G11 == 1) || (G11 == -1)); // !!!!! Test invertibility
        CoefficientType G11_inv = G11; // Inverse of G11 (which is itself for 1 or -1)
        
        // Remove the contributions of gamma from G12 and sigma from G21
        G12 /= std::vector<int>({gamma}); // Remove column gamma from G12
        G21 /= std::vector<int>({sigma}); // Remove row sigma from G21
        
        // Remove the corresponding row and column from _G_col
        _G_col[q].delRow(sigma); // Remove row sigma from _G_col[q]
        _G_col[q].delColumn(gamma); // Remove column gamma from _G_col[q]
        
        //---------------------------------------------- Submatrices of F -----------------------------------------------------
        
        // Extract the row chain from _F_row
        //        RChain F11(OSM::getRow(_F_row[q], gamma)); // F11 is the row chain from _F_row[q] at index gamma
        
        // Remove the row gamma from _F_row
        _F_row[q].delRow(gamma);
        
        //--------------------------------------------- Submatrices of H ------------------------------------------------------
        
        // Extract the row chain from _H_col
        RChain H11(OSM::getRow(_H_col[q-1], sigma)); // H11 is the row chain from _H_col[q] at index sigma
        
        // Remove the row sigma from _H_col
        _H_col[q-1].delRow(sigma);
        
        //--------------------------------------------- Submatrices of DD_q+1 ------------------------------------------------------
        
        // Extract the row chain from _DD_col[q+1]
        // RChain D11_q_plus1(OSM::getRow(_DD_col[q+1], gamma)); // D11_q_plus1 is the row chain from _DD_col[q + 1] at index gamma
        
        // Remove the row gamma from _DD_col
        if (q < _K.dim())
            _DD_col[q + 1].delRow(gamma);
        
        //--------------------------------------------- Submatrices of DD_q ------------------------------------------------------
        
        // Extract the column chain from _DD_col[q]
        CChain D11_q(OSM::getColumn(_DD_col[q], gamma)); // D11_q is the column chain from _DD_col[q] at index gamma
        
        // Remove the column gamma from _DD_col
        _DD_col[q].delColumn(gamma);
        
        //--------------------------------------------- Update matrices -------------------------------------------------------
        
        // Update _H_col
        // Subtract the product of (G21 * G11_inv) and H11 from _H_col[q]
        // Note: % operator is used for matrix row results
        _H_col[q-1] -= (G21 * G11_inv) * H11;
        
        // Set the row gamma of _H_col[q] to (-1 * H11) * G11_inv
        OSM::setRow(_H_col[q-1], gamma, (-1 * H11) * (G11_inv));
        
        // Update _G_col
        // Subtract the product of G11_inv and (G21 * G12) from _G_col[q]
        _G_col[q] -= (G21 * G11_inv) * G12;
        
        // Set the row gamma of _G_col[q] to G11_inv * G12
        OSM::setRow(_G_col[q], gamma, -G11_inv * G12);
        
        // Set the column sigma of _G_col[q] to G11_inv * G21
        OSM::setColumn(_G_col[q], sigma, G21 * G11_inv);
        
        // Set the coefficient at (gamma, sigma) in _G_col[q] to G11_inv
        _G_col[q].set_coef(gamma, sigma, G11_inv);
        
        // Update _F_row (for dimension q-1)
        // Add the product of (D11_q * G11_inv) and H11 to _F_row[q - 1]
        _F_row[q - 1] += (D11_q * G11_inv) * H11;
        
        // Update _DD_col (for dimension q)
        _DD_col[q] -= (D11_q * G11_inv) * G12;// Subtract the product of (D11_q * G11_inv) and G12
        OSM::setColumn(_DD_col[q], sigma, D11_q * G11_inv) ;
        
        // Update _F_row (for dimension q) and _DD_col (for dimension q+1)
        if (q < _K.dim())
        {
            // Extract boundary chain and project it
            RChain c = _K.cod(sigma, q); // Boundary of sigma in dimension q
            RChain projection_s(projection(c, SECONDARY, q+1)); // Project boundary chain to SECONDARY
            RChain projection_c(projection(c, CRITICAL, q+1)); // Project boundary chain to SECONDARY
            
            // Set the row sigma of _F_row[q] to (-1 * projection_s * _H_col[q])
            OSM::setRow(_F_row[q], sigma, (-1 * projection_s) * _H_col[q]);
            
            // Set the row sigma of _DD_col[q + 1] to projection_s * _G_col[q + 1] + projection_c
            
            RChain tmp(projection_s * _G_col[q + 1] + projection_c) ;
            OSM::setRow(_DD_col[q + 1], sigma, tmp);
        }
        // else : if the dimension is maximal, no update of _DD_col[q+1] and _F_row[q]
        
        // Update flags
        _flag[q][gamma] = SECONDARY; // Set the flag of gamma in dimension q to SECONDARY
        _flag[q][sigma] = CRITICAL; // Set the flag of sigma in dimension q to CRITICAL
    }
    else
        std::cout << "!!! W impossible with partial reduction options" << std::endl;
}

// Method to perform operation MW
// gamma is in dimension q, sigma is in dimension q

template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
void HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::MW(int pi, int sigma, int q) {
    
    //----------------------------------------------- Submatrices of G ----------------------------------------------------
    
    if (_hdvf_opt & OPT_FULL)
    {
        // Output the operation details to the console
        std::cout << "MW_" << q << "(" << pi << "," << sigma << ")" << std::endl;
        
        if (q <= 0)
            throw("MW operation in dimension 0 !!!") ;
        if (q >= _K.dim())
            throw("MW operation in maximal dimension !!!") ;
        
        // In order to compute xi and xi', extract sub-matrices of H_q, H_q-1 and compute d(pi) and cod(sigma)
        
        // H_q extractions
        
        CChain H11(OSM::getColumn(_H_col.at(q), pi)) ;
        // H21 -> delete H11
        _H_col.at(q) /= std::vector<int>({pi}) ;
        
        // H_q-1 extractions
        
        RChain H11q1(OSM::getRow(_H_col.at(q-1), sigma)) ;
        // H21_q-1 -> delete H11q1
        _H_col.at(q-1).delRow(sigma) ;
        
        // d(pi)
        
        CChain d_pi = _K.d(pi, q) ;
        CChain projP_d_pi = projection(d_pi, PRIMARY, q-1) ;
        CChain projC_d_pi = projection(d_pi, CRITICAL, q-1) ;
        
        // cod(sigma)
        
        RChain cod_sigma = _K.cod(sigma, q) ;
        RChain projS_cod_sigma = projection(cod_sigma, SECONDARY, q+1) ;
        RChain projC_cod_sigma = projection(cod_sigma, CRITICAL, q+1) ;
        
        // Compute xi and xi' to test the validity of MW
        
        CoefficientType xi = projS_cod_sigma * H11 ;
        CoefficientType xip = H11q1 * projP_d_pi ;
        
        if (abs(xi) != 1)
            throw "MW impossible, xi non invertible" ;
        if (abs(xip) != 1)
            throw "MW impossible, xi' non invertible" ;
        
        // F_q extraction
        
        CChain F11(OSM::getColumn(_F_row.at(q), pi)) ;
        // F12 -> delete col F11
        _F_row.at(q).delColumn(pi) ;
        
        // G_q extractions
        
        RChain G11(OSM::getRow(_G_col.at(q), sigma)) ;
        // G21 -> dele row G11
        _G_col.at(q).delRow(sigma) ;
        
        // ----------- Update of the reduction
        
        // H_q
        
        RChain tmp1 = projS_cod_sigma * _H_col.at(q) ;
        
        _H_col.at(q) -= (H11 * xi) * tmp1 ;
        OSM::setColumn(_H_col.at(q), sigma, H11 * xi) ;
        
        // F_q
        
        _F_row.at(q) += (F11 * xi) * tmp1 ;
        OSM::setColumn(_F_row.at(q), sigma, F11 * (-xi)) ;
        
        // G_q+1 // note: G_q+1 is not be modified if the HDVF is perfect
        
        RChain tmp2(projS_cod_sigma * _G_col.at(q+1)) ;
        tmp2 += projC_cod_sigma ;
        _G_col.at(q+1) -= (H11 * xi) * tmp2 ;
        
        // DD_col_q+1 / DD_row_q
        
        _DD_col.at(q+1) -= (F11 * xi) * tmp2 ;
        
        // H_q-1
        
        CChain tmp3(_H_col.at(q-1) * projP_d_pi) ;
        
        _H_col.at(q-1) -= (tmp3 * xip) * H11q1 ;
        OSM::setRow(_H_col.at(q-1), pi, xip * H11q1) ;
        
        // G_q
        
        _G_col.at(q) -= (tmp3 * xip) * G11 ;
        OSM::setRow(_G_col.at(q), pi, xip * G11) ;
        
        // F_q-1 // note: F_q-1 is not be modified if the HDVF is perfect
        
        CChain tmp4(_F_row.at(q-1) * projP_d_pi) ;
        tmp4 += projC_d_pi ;
        _F_row.at(q-1) -= (tmp4 * xip) * H11q1 ;
        
        // DD_col_q
        
        _DD_col.at(q) += (tmp4 * xip) * G11 ;
        
        // Update flags
        _flag[q][pi] = SECONDARY; // Set the flag of gamma in dimension q to SECONDARY
        _flag[q][sigma] = PRIMARY; // Set the flag of sigma in dimension q to CRITICAL
    }
    else
        std::cout << "!!! MW impossible with partial reduction options" << std::endl;
}


// Method to compute a perfect HDVF (Higher Dimensional Vector Field)
// Returns a vector of PairCell objects representing the pairs found

template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::computePerfectHDVF(bool verbose) {
    std::vector<PairCell> pair_list; // Vector to store the list of pairs
    bool trouve = false; // Flag to indicate whether a pair was found
    int dim = _K.dim(); // Get the dimension of the complex K
    
    // Loop through dimensions from q-1 to 0
    for (int q = dim - 1; q >= 0; --q) {
        std::cout << "-> pairing cells of dimension " << q << " and " << q+1 << std::endl ;
        //        cout << _K.nb_cells(q) << " cells of dimension " << q << endl ;
        // Incorrect: the number of cells is the number of cols in _DD_col ... (duality)
        
        // Find a pair of cells in dimension q
        PairCell pair = findPairA(q, trouve);
        
        // While a pair is found
        while (trouve) {
            progress_bar(_K.nb_cells(q)-_nb_C.at(q), _K.nb_cells(q)) ;
            // Add the found pair to the list
            pair_list.push_back(pair);
            
            // Perform operation A with the found pair
            A(pair.sigma, pair.tau, q);
            if (verbose)
            {
                std::cout << "A : " << pair.sigma << " - " << pair.tau << " (dim " << pair.dim << ")" << std::endl ;
                print_matrices(std::cout) ;
            }
            
            // Find another pair of cells in dimension q
            pair = findPairA(q, trouve);
        }
    }
    
    // Return the list of pairs found
    return pair_list;
}

// Method to compute a random perfect HDVF (Higher Dimensional Vector Field)
// Returns a vector of PairCell objects representing the pairs found

template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<PairCell> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::computeRandPerfectHDVF(bool verbose) {
    std::vector<PairCell> pair_list; // Vector to store the list of pairs
    bool trouve = false; // Flag to indicate whether a pair was found
    int dim = _K.dim(); // Get the dimension of the complex K
    
    // Init random generator
    std::random_device dev;
    std::mt19937 rng(dev()); // Random
    
    // Loop through dimensions from q-1 to 0
    for (int q = dim - 1; q >= 0; --q) {
        std::cout << "-> pairing cells of dimension " << q << " and " << q+1 << std::endl ;
        // Incorrect: the number of cells is the number of cols in _DD_col ... (duality)
        
        std::vector<PairCell> pairs = findPairsA(q, trouve);
        
        PairCell pair ;
        
        // While a pair is found
        while (trouve)
        {
            // Add one of the pairs (randomly) to the list
            {
                // Pickup a random cell sigma
                std::uniform_int_distribution<std::mt19937::result_type> rand_dist(0,pairs.size()-1);
                int i(rand_dist(rng)) ;
                pair = pairs.at(i) ;
            }
            pair_list.push_back(pair);
            
            // Perform operation A with the chosen pair
            A(pair.sigma, pair.tau, pair.dim);
            if (verbose)
            {
                std::cout << "A : " << pair.sigma << " - " << pair.tau << " (dim " << pair.dim << ")" << std::endl ;
                print_matrices(std::cout) ;
            }
            
            // Compute possible pairings
            pairs = findPairsA(q, trouve);
        }
    }
    // Return the list of pairs found
    return pair_list;
}

// Method to get cells if with a given flag (P,S,C) for each dimension
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<std::vector<int> > HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::get_flag (FlagType flag) const
{
    std::vector<std::vector<int> > res(_K.dim()+1) ;
    for (int q=0; q<=_K.dim(); ++q)
    {
        for (int i=0; i<_K.nb_cells(q); ++i)
        {
            if (_flag.at(q).at(i) == flag)
                res.at(q).push_back(i) ;
        }
    }
    return res ;
}

// Method to get cells with a given flag (P,S,C) for a given dimension
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::vector<int> HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::get_flag_dim (FlagType flag, int q) const
{
    std::vector<int> res ;
    for (int i=0; i<_K.nb_cells(q); ++i)
    {
        if (_flag.at(q).at(i) == flag)
            res.push_back(i) ;
    }
    return res ;
}

// Method to print the current state of the reduction
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::ostream& HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::print_reduction(std::ostream& out) // const
{
    // Print PSC
    out << "----- flags of cells:" << std::endl;
    for (int q = 0; q <= _K.dim(); ++q) {
        out << "--- dim " << q << std::endl;
        for (int i = 0; i < _K.nb_cells(q); ++i)
        {
            const int flag(_flag.at(q).at(i)) ;
            if (flag == PRIMARY)
                out << i << " -> P" << std::endl ;
            else if (flag == SECONDARY)
                out << i << " -> S" << std::endl ;
            else
                out << i << " -> C" << std::endl ;
        }
        out << std::endl;
    }
    
    // Print critical cells
    out << "----- critical cells:" << std::endl;
    std::vector<std::vector<int> > critical(get_flag(CRITICAL)) ;
    for (int q = 0; q <= _K.dim(); ++q) {
        out << "--- dim " << q << std::endl;
        for (int i = 0; i < critical.at(q).size(); ++i)
        {
            out << critical.at(q).at(i) << " ";
        }
        out << std::endl;
    }
    
    if (_hdvf_opt & (OPT_FULL | OPT_G))
    {
        // Print matrices g
        out << "----- g:" << std::endl;
        for (int q = 0; q <= _K.dim(); ++q) {
            out << "--- dim " << q << std::endl;
            for (int i = 0; i < critical.at(q).size(); ++i)
            {
                const int id(critical.at(q).at(i)) ;
                out << "g(" << id << ") = (" << id << ")";
                // Iterate over the ith column of _G_col
                CChain col(OSM::getColumn(_G_col.at(q), id)) ; // TODO cget
                for (typename CChain::const_iterator it_col = col.cbegin(); it_col != col.cend(); ++it_col) {
                    out << " + " << it_col->second << ".(" << it_col->first << ") + ";
                }
                out << std::endl;
            }
        }
    }
    
    if (_hdvf_opt & (OPT_FULL | OPT_F))
    {
        // Print matrices f*
        out << "----- f*:" << std::endl;
        for (int q = 0; q <= _K.dim(); ++q) {
            out << "--- dim " << q << std::endl;
            for (int i = 0; i < critical.at(q).size(); ++i)
            {
                const int id(critical.at(q).at(i)) ;
                out << "f*(" << id << ") = (" << id << ")";
                // Iterate over the ith row of _F_row
                RChain row(OSM::getRow(_F_row.at(q), id)) ; // TODO cget
                for (typename RChain::const_iterator it_row = row.cbegin(); it_row != row.cend(); ++it_row) {
                    out << " + " << it_row->second << ".(" << it_row->first << ") + ";
                }
                out << std::endl;
            }
        }
    }
    return out ;
}

// Method to print A-pairs
//template<typename CoefficientType, typename ComplexType>
template<typename CoefficientType, typename ComplexType, template <typename, int> typename _ChainType, template <typename, int> typename _SparseMatrixType>
std::ostream& HDVF<CoefficientType, ComplexType, _ChainType, _SparseMatrixType>::print_pairs(const std::vector<PairCell>& pairs, std::ostream& out) // const
{
    for (const auto& pair : pairs) {
        out << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }
    return out ;
}

}

#endif // HDVF_H
