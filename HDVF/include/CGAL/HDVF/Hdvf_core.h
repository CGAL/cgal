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

#ifndef CGAL_HDVF_HDVF_CORE_H
#define CGAL_HDVF_HDVF_CORE_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <CGAL/OSM/OSM.h>
#include <CGAL/OSM/Bitboard.h>

namespace CGAL {
namespace HDVF {

/** \brief HDVF Enum for the label of cells. */
static enum FlagType {
    PRIMARY,
    SECONDARY,
    CRITICAL,
    NONE // For duality and persistence
};

/** \brief HDVF option (compute only reduced boundary). */
const int OPT_BND = 0b0001;
/** \brief HDVF option (compute only reduced boundary and f). */
const int OPT_F = 0b0010;
/** \brief HDVF option (compute only reduced boundary and g). */
const int OPT_G = 0b0100;
/** \brief HDVF option (compute full reduction). */
const int OPT_FULL = 0b1000;


/** \brief Structure to represent data for HDVF operations (pairs of cells).
 *
 * Cells are always sorted so that the dimension of `sigma` is lesser than the dimension of `tau`.
 */
struct Pair_cells {
    size_t sigma;  /// Index of the first cell
    size_t tau;    /// Index of the second cell
    int dim;    /// Dimension of cells: `dim`/`dim`+1 for A and R, `dim`/`dim` for other operations
};

/** \brief Overload of operator<< for Pair_cells type. */
inline std::ostream& operator<<(std::ostream &out, const std::vector<Pair_cells>& pairs)
{
    for (const auto& pair : pairs) {
        out << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }
    return out ;
}

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf_core` is the core implementation of homological discrete vector fields (HDVF for short).

 An enumeration `FlagType` is defined in the `HDVF` namespace and the `Hdvf_core` class maps each cell to one of the flags (namely `PRIMARY`, `SECONDARY`, `CRITICAL`). The NONE flag is used in child classes (such as `Hdvf_duality` or `Hdvf_reduced`) when computing relative homology on a sub-complex.
 The flag of each cell is stored in an appropriate structure and getters are provided to access to this information.

 The `Hdvf_core` class stores the associated reduction in sparse matrices: row-major for \f$f\f$, and column-major for \f$g\f$, \f$h\f$ and \f$\partial'\f$. Getters are provided to access this information. However, according to the chosen HDVF computation option (`OPT_BND`, `OPT_F`, `OPT_G`, `OPT_FULL`) the reduction can be computed only partially (and thus faster).

 The class provides perfect HDVF constuction operations: `compute_perfect_hdvf()` and `compute_rand_perfect_hdvf()`, which build perfect HDVFs by pairing iteratively critical cells through the `A()` operation.

 If the user wishes to build an HDVF using other criteria, several `find_pair_A()` functions are provided (searching for valid pairs of cells for `A`respecting various constraints). The `A` operation can be applied to any pair returned by these functions.

 Homology/cohomology generators are actually algebraic objects, namely chains. Methods `homology_chain()` and `cohomology_chain()` return the homology and cohomology generator chain associated to a given critical cell. VTK export functions output all the cells of such chains with non zero coefficients.


 \cgalModels{HDVF}

 \tparam CoefficientType a model of the `Ring` concept  providing the ring used to compute homology.
 \tparam ComplexType a model of the `AbstractChainComplex` concept, providing the type of abstract chain complex used.
 \tparam ChainType a model of the `SparseChain` concept (by default, `OSM::Sparse_chain`), providing the type of sparse chains used (should be coherent with `SparseMatrixType`).
 \tparam SparseMatrixType a model of the `SparseMatrix` concept (by default, `OSM::Sparse_matrix`), providing the type of sparse matrices used.
 */


template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType = OSM::Sparse_chain, template <typename, int> typename SparseMatrixType = OSM::Sparse_matrix>
class Hdvf_core {
public:
    /*!
     Type of column-major chains
     */
    typedef ChainType<CoefficientType, CGAL::OSM::COLUMN> Col_chain;

    /*!
     Type of row-major chains
     */
    typedef ChainType<CoefficientType, CGAL::OSM::ROW> Row_chain;

    /*!
     Type of column-major sparse matrices
     */
    typedef SparseMatrixType<CoefficientType, CGAL::OSM::COLUMN> Col_matrix;

    /*!
     Type of row-major sparse matrices
     */
    typedef SparseMatrixType<CoefficientType, CGAL::OSM::ROW> Row_matrix;

protected:
    /* \brief Flags of the cells.
     * _flag.at(q) contains the flags of cells of dimension q
     */
    std::vector<std::vector<FlagType>> _flag;
    /* \brief Number of `PRIMARY` cells. */
    std::vector<size_t> _nb_P;
    /* \brief Number of `SECONDARY` cells. */
    std::vector<size_t> _nb_S;
    /* \brief Number of `CRITICAL` cells. */
    std::vector<size_t> _nb_C;
    /* \brief Row matrices for f. */
    std::vector<Row_matrix> _F_row;
    /* \brief Column matrices for g. */
    std::vector<Col_matrix> _G_col;
    /* \brief Column matrices for h. */
    std::vector<Col_matrix> _H_col;
    /* \brief Column matrices for reduced boundary. */
    std::vector<Col_matrix> _DD_col;

    /* \brief Reference to the underlying complex. */
    const ComplexType& _K;

    /* \brief Hdvf_core options for computation (computation of partial reduction). */
    int _hdvf_opt;

public:
    /**
     * \brief Default constructor.
     *
     * Builds an "empty" HDVF_core associated to K (with all cells critical). By default, the HDVF option is set to OPT_FULL (full reduction computed).
     *
     * \param[in] K A chain complex (a model of `AbstractChainComplex`)
     * \param[in] hdvf_opt Option for HDVF computation (`OPT_BND`, `OPT_F`, `OPT_G` or `OPT_FULL`)
     */
    Hdvf_core(const ComplexType& K, int hdvf_opt = OPT_FULL) ;

    /*
     * \brief Constructor by copy.
     *
     * Builds a HDVF by copy from another, including options.
     *
     * \param[in] hdvf An initial HDVF.
     */
    Hdvf_core(const Hdvf_core& hdvf) : _flag(hdvf._flag), _nb_P(hdvf._nb_P), _nb_S(hdvf._nb_S), _nb_C(hdvf._nb_C), _F_row(hdvf._F_row), _G_col(hdvf._G_col), _H_col(hdvf._H_col), _DD_col(hdvf._DD_col), _K(hdvf._K), _hdvf_opt(hdvf._hdvf_opt) { }

    /*
     * \brief HDVF_core destructor. */
    ~Hdvf_core() { }

    /**
     * \brief Finds a valid Pair_cells of dimension q / q+1 for A.
     *
     * The function searches a pair of critical cells \f$(\gamma_1, \gamma2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Lower dimension of the pair.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual Pair_cells find_pair_A(int q, bool &found) const;

    /**
     * \brief Finds a valid Pair_cells for A containing `gamma` (a cell of dimension `q`)
     *
     * The function searches a cell \f$\gamma'\f$ such that one of the following conditions holds:
     * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
     * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
     *
     * \param[in] q Dimension of the cell `gamma`.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] gamma Index of a cell to pair.
     */
    virtual Pair_cells find_pair_A(int q, bool &found, size_t gamma) const;

    /**
     * \brief Finds *all* valid Pair_cells of dimension q / q+1 for A.
     *
     * The function searches all pairs of critical cells \f$(\gamma_1, \gamma2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Lower dimension of the pair.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual std::vector<Pair_cells> find_pairs_A(int q, bool &found) const;

    /**
     * \brief Finds *all* valid Pair_cells for A containing `gamma` (a cell of dimension `q`)
     *
     * The function searches all `CRITICAL` cells \f$\gamma'\f$ such that one of the following conditions holds:
     * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
     * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the cell `gamma`.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] gamma Index of a cell to pair.
     */
    virtual std::vector<Pair_cells> find_pairs_A(int q, bool &found, size_t gamma) const;

    /**
     * \brief A operation: pairs critical cells.
     *
     * A pair of critical cells \f$(\gamma_1, \gamma_2)\f$ of respective dimension q and q+1 is valid for A if \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ is invertible. After the `A()` operation, \f$\gamma_1\f$ becomes `PRIMARY`, \f$\gamma_2\f$ becomes `SECONDARY`. The A method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param[in] gamma1 First cell of the pair (dimension `q`)
     * \param[in] gamma2 Second cell of the pair (dimension `q+1`)
     * \param[in] q Dimension of the pair
     */
    void A(size_t gamma1, size_t gamma2, int q);

    /**
     * \brief Computes a perfect HDVF.
     *
     * As long as valid pairs for A exist, the function selects the first available pair (returned by `find_pair_A`()) and applies the corresponding `A()` operation.
     * If the `Ring` of coefficients is a field, this operation always produces a perfect HDVF (ie.\ the reduced boundary is null and the reduction provides homology and cohomology information).
     * Otherwise the operation produces a maximal HDVF with a residual boundary matrix over critical cells.
     *
     * If the HDVF is initially not trivial (some cells have already been paired), the function completes it into a perfect HDVF.
     *
     * \param[in] verbose If this parameter is `true`, all intermediate reductions are printed out.
     *
     * \return The vector of all `Pair_cells` paired with A.
     */
    std::vector<Pair_cells> compute_perfect_hdvf(bool verbose = false);

    /**
     * \brief Computes a random perfect HDVF.
     *
     * As long as valid pairs for A exist, the function selects a random pair (among pairs returned by `find_pairs_A()`) and applies the corresponding `A()` operation.
     * If the `Ring` of coefficients is a field, this operation always produces a perfect HDVF (that  is the reduced boundary is null and the reduction provides homology and cohomology information).
     *
     * If the HDVF is initially not trivial (some cells have already been paired), the function randomly completes it into a perfect HDVF.
     *
     * \warning This method is slower that `compute_perfect_hdvf()` (finding out all possible valid pairs requires additional time).
     *
     * \param[in] verbose If this  parameter is `true`, all intermediate reductions are printed out.
     *
     * \return The vector of all pairs of cells used for apply A.
     */
    std::vector<Pair_cells> compute_rand_perfect_hdvf(bool verbose = false);

    /**
     * \brief Tests if a HDVF is perfect.
     *
     * The function returns `true` is the reduced boundary matrix is null and `false` otherwise.
     */
    bool is_perfect_hdvf()
    {
        bool res = true ;
        int q = 0 ;
        while ((q<=_K.dim()) && res)
        {
            res = res && _DD_col.at(q).is_null() ;
            ++q;
        }
        return res ;
    }

    // Hdvf_core getters

    /**
     * \brief Gets cells with a given `flag` in any dimension.
     *
     * The function returns a vector containing, for each dimension, the vector of cells with a given `flag`.
     *
     * \param[in] flag Flag to select.
     */

    // !!! Why should it be virtual for duality?????

    virtual std::vector<std::vector<size_t> > flag (FlagType flag) const ;

    /**
     * \brief Gets cells with a given `flag` in dimension `q`.
     *
     * The function returns the vector of cells of dimension `q` with a given `flag`.
     *
     * \param[in] flag Flag to select.
     * \param[in] q Dimension visited.
     */
    virtual std::vector<size_t> flag_dim (FlagType flag, int q) const ;

    /*!
     * \brief Gets the flag of the cell `tau` in dimension `q`.
     *
     * \param[in] tau Index of the cell.
     * \param[in] q Dimension of the cell.
     */
    FlagType cell_flag (int q, size_t tau) const { return _flag.at(q).at(tau); }

    /**
     * \brief Gets HDVF computation option.
     */
    int hdvf_opts () const { return _hdvf_opt ; }

    /**
     * \brief Gets the row-major matrix of \f$f\f$ (from the reduction associated to the HDVF).
     */
    const Row_matrix& matrix_f (int q) const { return _F_row.at(q); }

    /**
     * \brief Gets the column-major matrix of \f$g\f$ (from the reduction associated to the HDVF).
     */
    const Col_matrix& matrix_g (int q) const { return _G_col.at(q); }

    /**
     * \brief Gets the column-major matrix of \f$h\f$ (from the reduction associated to the HDVF).
     */
    const Col_matrix& matrix_h (int q) const { return _H_col.at(q); }

    /**
     * \brief Gets the column-major matrix of \f$\partial'\f$, reduced boundary operator (from the reduction associated to the HDVF).
     */
    const Col_matrix& matrix_dd (int q) const { return _DD_col.at(q); }

    /**
     * \brief Prints the matrices of the reduction.
     *
     * Prints the matrices of the reduction (that is \f$f\f$, \f$g\f$, \f$h\f$, \f$\partial'\f$ the reduced boundary).
     * By default, outputs the complex to `std::cout`.
    */
    std::ostream& insert_matrices(std::ostream &out = std::cout) const;

    /**
     * \brief Prints the homology and cohomology reduction information.
     *
     * Prints \f$f^*\f$, \f$g\f$ \f$\partial'\f$ the reduced boundary over each critical cell.
     *
     * By default, outputs the complex to `std::cout`.
    */
    std::ostream& insert_reduction(std::ostream &out = std::cout) const;


    /**
     * \brief Exports primary/secondary/critical labels (in particular for vtk export)
     *
     * The method exports the labels of every cells in each dimension.
     *
     * \return A vector containing, for each dimension, the vector of labels by cell index.
     */
    virtual std::vector<std::vector<int> > psc_labels () const
    {
        std::vector<std::vector<int> > labels(_K.dim()+1) ;
        for (int q=0; q<=_K.dim(); ++q)
        {
            for (size_t i = 0; i<_K.nb_cells(q); ++i)
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

    /**
     * \brief Gets homology generators associated to `cell` (critical cell) of dimension  `q` (used by vtk export).
     *
     * The method exports the chain \f$g(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell` and dimension `q`.
     *
     * \return A column-major chain.
     */
    virtual Col_chain homology_chain (size_t cell, int q) const
    {
        if ((q<0) || (q>_K.dim()))
            throw "Error : homology_chain with dim out of range" ;
        if (_hdvf_opt & (OPT_FULL | OPT_G))
        {
            Col_chain g_cell(OSM::get_column(_G_col.at(q), cell)) ;
            // Add 1 to the cell
            g_cell.set_coef(cell, 1) ;
            return g_cell ;
        }
        else
            throw "Error : trying to export g_chain without proper Hdvf_core option" ;
    }

    /**
     * \brief Gets cohomology generators associated to `cell` (critical cell) of dimension  `q` (used by vtk export).
     *
     * The method exports the chain \f$f^\star(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell` and dimension `q`.
     *
     * \param[in] cell Index of the (critical) cell.
     * \param[in] dim Dimension of the (critical) cell.
     *
     * \return A column-major chain.
     */
    virtual Col_chain cohomology_chain (size_t cell, int dim) const
    {
        if ((dim<0) || (dim>_K.dim()))
            throw "Error : cohomology_chain with dim out of range" ;
        if (_hdvf_opt & (OPT_FULL | OPT_F))
        {
            Row_chain fstar_cell(OSM::get_row(_F_row.at(dim), cell)) ;
            // Add 1 to the cell
            fstar_cell.set_coef(cell, 1) ;

            return fstar_cell.transpose() ;

        }
        else
            throw "Error : trying to export fstar_chain without proper Hdvf_core option" ;
    }

    /**
     * \brief Saves a HDVF together with the associated reduction (f, g, h, d matrices)
     *
     * Save a HDVF to a `.hdvf` file, a simple text file format (see for a specification).
     */
    std::ostream& insert_hdvf_reduction(std::ostream& out) ;

    /**
     * \brief Loads a HDVF together with the associated reduction (f, g, h, d matrices)
     *
     * Load a HDVF and its reduction from a `.hdvf` file, a simple text file format (see for a specification).
     * \warning The underlying complex is not stored in the file!
     */
    std::istream& extract_hdvf_reduction(std::istream& out) ;

protected:
    /* \brief Project a chain onto a given flag
     * The methods cancels all the coefficients of the chain of dimension `q` that do not correspond to`flag`. This is actually an implementation of the projection operator onto the sub A-module generated by cells of flag `flag`.
     *
     * \param[in] chain The chain projected.
     * \param[in] flag The flag onto which the chain is projected.
     * \param[in] q Dimension of the chain
     *
     * \result Returns a copy of `chain` where only coefficients of cells of flag `flag` are kept (all other coefficients are cancelled).
     */
    template<int ChainTypeFlag>
    ChainType<CoefficientType, ChainTypeFlag> projection(const ChainType<CoefficientType, ChainTypeFlag>& chain, FlagType flag, int q) const {
        // Create a new chain to store the result
        // Better to initialize 'result' directly with the correct size and iterate over it
        ChainType<CoefficientType, ChainTypeFlag> result(chain);

        // Iterate over each element of the chain
        std::vector<size_t> tmp ;
        for (typename ChainType<CoefficientType, ChainTypeFlag>::const_iterator it = result.cbegin(); it != result.cend(); ++it)
        {
            size_t cell_index = it->first;
            CoefficientType value = it->second;

            // Check the flag of the corresponding cell
            if (_flag[q][cell_index] != flag) {
                // Mark for cancellation
                tmp.push_back(cell_index) ;

            }
        }
        // Cancel all coefficients of tmp
        result /= tmp ;
        return result;
    }

    /* \brief Display a text progress bar. */
    void progress_bar(size_t i, size_t n)
    {
        const size_t step(n/20) ;
        if ((i%step)==0)
        {
            const float percentage(float(i)/(n-1)) ;
            const char PBSTR[] = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" ;
            const size_t PBWIDTH(60) ;
            size_t val = (size_t) (percentage * 100);
            size_t lpad = (size_t) (percentage * PBWIDTH);
            size_t rpad = PBWIDTH - lpad;
            printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
            fflush(stdout);
        }
        if (i==(n-1))
            std::cout << std::endl ;
    }

};

// Constructor for the Hdvf_core class
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::Hdvf_core(const ComplexType& K, int hdvf_opt) : _K(K) {
    // Get the dimension of the simplicial complex
    int dim = _K.dim();
    std::cout << "----> Starting Hdvf_core creation / dim " << dim << std::endl ;
    // Hdvf_core options
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
            _F_row[q] = Row_matrix(_K.nb_cells(q), _K.nb_cells(q));

        // Initialize _G_col[q] as a column matrix with dimensions (dim(q) x dim(q))
        if (_hdvf_opt & (OPT_FULL | OPT_G))
            _G_col[q] = Col_matrix(_K.nb_cells(q), _K.nb_cells(q));

        // Initialize _H_col[q] as a column matrix with dimensions (dim(q+1) x dim(q))
        if (_hdvf_opt & OPT_FULL)
            _H_col[q] = Col_matrix(_K.nb_cells(q + 1), _K.nb_cells(q));

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
        _DD_col.at(q) = _K.boundary_matrix(q) ;
    std::cout << "------> End Hdvf_core creation" << std::endl ;
}

// Method to print the matrices
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::ostream& Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::insert_matrices(std::ostream& out) const {
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

// find a valid Pair_cells for A in dimension q
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
Pair_cells Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::find_pair_A(int q, bool &found) const
{
    found = false;
    Pair_cells p;

    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = _DD_col[q+1].begin(); (it_col != _DD_col[q+1].end() && !found); ++it_col)
    {
        const Col_chain& col(OSM::cget_column(_DD_col[q+1], *it_col)) ;

        // Iterate through the entries of the column
        for (typename Col_chain::const_iterator it = col.begin(); (it != col.end() && !found); ++it) {
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

// find a valid Pair_cells containing tau for A in dimension q
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
Pair_cells Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::find_pair_A(int q, bool &found, size_t gamma) const
{
    found = false;
    Pair_cells p ;

    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    const Col_chain& tmp2(OSM::cget_column(_DD_col.at(q), gamma)) ;
    for (typename Col_chain::const_iterator it = tmp2.cbegin(); (it != tmp2.cend() && !found); ++it)
    {
        if (abs(it->second) == 1)
        {
            found = true ;
            p.sigma = it->first ;
            p.tau = gamma ;
            p.dim = q-1 ;
        }
    }

    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    Row_chain tmp(OSM::get_row(_DD_col.at(q+1), gamma)) ;
    for (typename Row_chain::const_iterator it = tmp.cbegin(); (it != tmp.cend() && !found); ++it)
    {
        if (abs(it->second) == 1)
        {
            found = true ;
            p.sigma = gamma ;
            p.tau = it->first ;
            p.dim = q ;
        }
    }
    return p;
}

// find all the valid PairCells for A in dimension q
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::vector<Pair_cells> Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::find_pairs_A(int q, bool &found) const
{
    std::vector<Pair_cells> pairs;
    found = false ;

    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = this->_DD_col[q+1].begin(); it_col != this->_DD_col[q+1].end(); ++it_col)
    {
        const Col_chain& col(OSM::cget_column(this->_DD_col[q+1], *it_col)) ;

        // Iterate through the entries of the column
        for (typename Col_chain::const_iterator it = col.begin(); it != col.end(); ++it) {
            if ((it->second == 1) || (it->second == -1)) {
                // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
                Pair_cells p;
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

// find all the valid Pair_cells containing gamma for A in dimension q
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::vector<Pair_cells> Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::find_pairs_A(int q, bool &found, size_t gamma) const
{
    found = false;
    std::vector<Pair_cells> pairs;

    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    Row_chain tmp(OSM::get_row(_DD_col.at(q+1), gamma)) ;
    for (typename Row_chain::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
    {
        if (abs(it->second) == 1)
        {
            found = true ;
            Pair_cells p ;
            p.sigma = gamma ;
            p.tau = it->first ;
            p.dim = q ;
            pairs.push_back(p) ;
        }
    }
    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    const Col_chain& tmp2(OSM::cget_column(_DD_col.at(q), gamma)) ;
    for (typename Col_chain::const_iterator it = tmp2.cbegin(); it != tmp2.cend(); ++it)
    {
        if (abs(it->second) == 1)
        {
            found = true ;
            Pair_cells p ;
            p.sigma = it->first ;
            p.tau = gamma ;
            p.dim = q-1 ;
            pairs.push_back(p) ;
        }
    }
    return pairs;
}




// Method to perform operation A
// tau1 is in dimension q, tau2 is in dimension q+1
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
void Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::A(size_t tau1, size_t tau2, int q) {
    //----------------------------------------------- Submatrices of D ----------------------------------------------------

    // Output operation details to the console
    //    cout << "A of " << tau1 << "(dim " << q << ") / " << tau2 << "(dim " << q + 1 << ")" << endl;

    // Extract submatrices from _DD_col
    Row_chain D12(OSM::get_row(_DD_col.at(q+1),tau1)); // D12 is a row chain from _DD_col[q+1] at index tau1
    Col_chain D21(OSM::get_column(_DD_col.at(q + 1),tau2)); // D21 is a column chain from _DD_col[q+1] at index tau2
    const CoefficientType D11 = D12.get_coef(tau2); // D11 is the coefficient at the intersection of tau2 in D12

    // Assert that D11 is either 1 or -1 (check invertibility)
    assert((D11 == 1) || (D11 == -1)); // !!!!! Test invertibility

    // Compute the inverse of D11 (which is itself, since D11 is 1 or -1)
    CoefficientType D11_inv = D11;

    // Perform operations to remove the row and column contributions
    D12 /= std::vector<size_t>({tau2}); // Remove tau2 column from D12
    D21 /= std::vector<size_t>({tau1}); // Remove tau1 row from D21

    // Delete rows and columns from _DD_col
    del_row(_DD_col[q + 1], tau1); // Remove row tau1 from _DD_col[q+1]
    del_column(_DD_col[q + 1], tau2); // Remove column tau2 from _DD_col[q+1]

    //---------------------------------------------- Submatrices of F -----------------------------------------------------

    Row_chain F11 ;
    Col_chain G11 ;
    if (_hdvf_opt & (OPT_FULL | OPT_F))
    {
        // Extract the relevant submatrix from _F_row
        F11 = OSM::get_row(_F_row.at(q),tau1); // F11 is a row chain from _F_row[q] at index tau1

        // Delete the row tau1 from _F_row
        del_row(_F_row[q], tau1);
    }

    //--------------------------------------------- Submatrices of G ------------------------------------------------------

    if (_hdvf_opt & (OPT_FULL | OPT_G))
    {
        // Extract the relevant submatrix from _G_col
        G11 = OSM::get_column(_G_col.at(q + 1),tau2); // G11 is a column chain from _G_col[q+1] at index tau2

        // Delete the column tau2 from _G_col
        del_column(_G_col[q + 1], tau2);
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
        OSM::set_column(_F_row[q], tau1, D21 * (-D11_inv));

        // Remove the row tau2 from _F_row[q+1]
        del_row(_F_row[q + 1], tau2);
    }

    // ---- Update _G_col

    if (_hdvf_opt & (OPT_FULL | OPT_G))
    {
        // Update _G_col[q + 1]
        // Subtract the product of (G11 * D11_inv) and D12 from _G_col[q+1]
        _G_col[q + 1] -= (G11 * D11_inv) * D12;

        // Set the row tau2 of _G_col[q + 1] to (D12 * (-D11_inv))
        OSM::set_row(_G_col[q + 1], tau2, D12 * (-D11_inv));

        // Remove the column tau1 from _G_col[q]
        del_column(_G_col[q], tau1);
    }

    // ---- Update _H_col

    if (_hdvf_opt & OPT_FULL)
    {
        // Update _H_col[q]
        // Compute the temporary matrix product G11 * F11
        Col_matrix tmp = G11 * F11;

        // Add the product of tmp and D11_inv to _H_col[q]
        _H_col[q] += (tmp) * D11_inv;

        // Set the row tau2 of _H_col[q] to (F11 * D11_inv)
        OSM::set_row(_H_col[q], tau2, F11 * D11_inv);

        // Set the column tau1 of _H_col[q] to (G11 * D11_inv)
        OSM::set_column(_H_col[q], tau1, G11 * D11_inv);

        // Set the coefficient at (tau2, tau1) in _H_col[q] to D11_inv
        set_coef(_H_col[q], tau2, tau1, D11_inv);
    }

    // ---- Update _DD_col

    // Update _DD_col
    _DD_col[q + 1] -= (D21 * D12) * D11_inv;

    // Remove columns and rows from _DD_col as necessary
    if (q > 0) {
        del_column(_DD_col[q], tau1); // Remove column tau1 from _DD_col[q]
    }
    if (q + 2 <= _K.dim()) {
        del_row(_DD_col[q + 2], tau2); // Remove row tau2 from _DD_col[q+2]
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


// Method to compute a perfect Hdvf_core
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::vector<Pair_cells> Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::compute_perfect_hdvf(bool verbose) {
    std::vector<Pair_cells> pair_list; // Vector to store the list of pairs
    bool trouve = false; // Flag to indicate whether a pair was found
    int dim = _K.dim(); // Get the dimension of the complex K

    // Loop through dimensions from q-1 to 0
    for (int q = dim - 1; q >= 0; --q) {
        std::cout << std::endl << "-> pairing cells of dimension " << q << " and " << q+1 << std::endl ;

        // Find a pair of cells in dimension q
        Pair_cells pair = find_pair_A(q, trouve);

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
                insert_matrices(std::cout) ;
            }

            // Find another pair of cells in dimension q
            pair = find_pair_A(q, trouve);
        }
    }

    // Return the list of pairs found
    return pair_list;
}

// Method to compute a random perfect Hdvf_core
// Returns a vector of Pair_cells objects representing the pairs found
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::vector<Pair_cells> Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::compute_rand_perfect_hdvf(bool verbose) {
    std::vector<Pair_cells> pair_list; // Vector to store the list of pairs
    bool trouve = false; // Flag to indicate whether a pair was found
    int dim = _K.dim(); // Get the dimension of the complex K

    // Init random generator
    std::random_device dev;
    std::mt19937 rng(dev()); // Random

    // Loop through dimensions from q-1 to 0
    for (int q = dim - 1; q >= 0; --q) {
        std::cout << "-> pairing cells of dimension " << q << " and " << q+1 << std::endl ;
        // Incorrect: the number of cells is the number of cols in _DD_col ... (duality)

        std::vector<Pair_cells> pairs = find_pairs_A(q, trouve);

        Pair_cells pair ;

        // While a pair is found
        while (trouve)
        {
            // Add one of the pairs (randomly) to the list
            {
                // Pickup a random cell sigma
                std::uniform_int_distribution<std::mt19937::result_type> rand_dist(0,pairs.size()-1);
                size_t i(rand_dist(rng)) ;
                pair = pairs.at(i) ;
            }
            pair_list.push_back(pair);

            // Perform operation A with the chosen pair
            A(pair.sigma, pair.tau, pair.dim);
            if (verbose)
            {
                std::cout << "A : " << pair.sigma << " - " << pair.tau << " (dim " << pair.dim << ")" << std::endl ;
                insert_matrices(std::cout) ;
            }

            // Compute possible pairings
            pairs = find_pairs_A(q, trouve);
        }
    }
    // Return the list of pairs found
    return pair_list;
}

// Method to get cells if with a given flag (P,S,C) for each dimension
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::vector<std::vector<size_t> > Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::flag (FlagType flag) const
{
    std::vector<std::vector<size_t> > res(_K.dim()+1) ;
    for (int q=0; q<=_K.dim(); ++q)
    {
        for (size_t i=0; i<_K.nb_cells(q); ++i)
        {
            if (_flag.at(q).at(i) == flag)
                res.at(q).push_back(i) ;
        }
    }
    return res ;
}

// Method to get cells with a given flag (P,S,C) for a given dimension
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::vector<size_t> Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::flag_dim (FlagType flag, int q) const
{
    std::vector<size_t> res ;
    for (size_t i=0; i<_K.nb_cells(q); ++i)
    {
        if (_flag.at(q).at(i) == flag)
            res.push_back(i) ;
    }
    return res ;
}

// Method to print the current state of the reduction
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::ostream& Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::insert_reduction(std::ostream& out) const
{
    // Print PSC
    out << "----- flags of cells:" << std::endl;
    for (int q = 0; q <= _K.dim(); ++q) {
        out << "--- dim " << q << std::endl;
        for (size_t i = 0; i < _K.nb_cells(q); ++i)
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
    std::vector<std::vector<size_t> > critical(flag(CRITICAL)) ;
    for (int q = 0; q <= _K.dim(); ++q) {
        out << "--- dim " << q << std::endl;
        for (size_t i = 0; i < critical.at(q).size(); ++i)
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
            for (size_t i = 0; i < critical.at(q).size(); ++i)
            {
                const size_t id(critical.at(q).at(i)) ;
                out << "g(" << id << ") = (" << id << ")";
                // Iterate over the ith column of _G_col
                Col_chain col(OSM::get_column(_G_col.at(q), id)) ; // TODO cget
                for (typename Col_chain::const_iterator it_col = col.cbegin(); it_col != col.cend(); ++it_col) {
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
            for (size_t i = 0; i < critical.at(q).size(); ++i)
            {
                const size_t id(critical.at(q).at(i)) ;
                out << "f*(" << id << ") = (" << id << ")";
                // Iterate over the ith row of _F_row
                Row_chain row(OSM::get_row(_F_row.at(q), id)) ; // TODO cget
                for (typename Row_chain::const_iterator it_row = row.cbegin(); it_row != row.cend(); ++it_row) {
                    out << " + " << it_row->second << ".(" << it_row->first << ") + ";
                }
                out << std::endl;
            }
        }
    }
    return out ;
}

// Save HDVF and reduction
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::ostream& Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::insert_hdvf_reduction(std::ostream& out)
{
    // HDVF save type
    // 0: HDVF and reduction
    // 1: HDVF only
    out << 0 << std::endl ;
    // Dimension
    out << _K.dim() << std::endl ;
    // Number of cells in each dimension
    for (int q=0; q<=_K.dim(); ++q)
        out << _K.nb_cells(q) << " " ;
    out << std::endl ;
    // Flags
    // P : -1 / S : 1 / C : 0
    // Each dimension written on a row
    for (int q=0; q<=_K.dim(); ++q)
    {
        for (int i=0; i<_K.nb_cells(q); ++i)
        {
            if (_flag.at(q).at(i) == PRIMARY)
                out << -1 << " " ;
            else if (_flag.at(q).at(i) == SECONDARY)
                out << 1 << " " ;
            else // CRITICAL
                out << 0 << " " ;
        }
        out << std::endl ;
    }
    // F
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::write_matrix(_F_row.at(q), out) ;
    }
    // G
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::write_matrix(_G_col.at(q), out) ;
    }
    // H
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::write_matrix(_H_col.at(q), out) ;
    }
    // DD
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::write_matrix(_DD_col.at(q), out) ;
    }
    return out;
}

// Save HDVF and reduction
template<typename CoefficientType, typename ComplexType, template <typename, int> typename ChainType, template <typename, int> typename SparseMatrixType>
std::istream& Hdvf_core<CoefficientType, ComplexType, ChainType, SparseMatrixType>::extract_hdvf_reduction(std::istream& in)
{
    // Load and check HDVF save type
    int type ;
    in >> type ;
    if (type != 0)
    {
        std::cerr << "extract_hdvf_reduction error: trying to load a pure HDVF file..." << std::endl ;
        throw ("extract_hdvf_reduction error: trying to load a pure HDVF file...");
    }

    // Load and check dimension
    int d ;
    in >> d ;
    if (d != _K.dim())
    {
        std::cerr << "extract_hdvf_reduction error: dimension loaded incompatible with the dimension of the underlying complex" << std::endl ;
        throw ("extract_hdvf_reduction error: dimension loaded incompatible with the dimension of the underlying complex");
    }
    // Load and check number of cells
    int nb ;
    for (int q=0; q<=_K.dim(); ++q)
    {
        in >> nb ;
        if (nb != _K.nb_cells(q))
        {
            std::string mess("extract_hdvf_reduction error: incoherent number of cells in dimension ");
            mess += std::to_string(q);
            std::cerr << mess << std::endl ;
            throw (mess);
        }
    }
    // Load flags
    int flag ;
    for (int q=0; q<=_K.dim(); ++q)
    {
        for (int i=0; i<_K.nb_cells(q); ++i)
        {
            in >> flag ;
            if (flag == -1)
                _flag.at(q).at(i) = PRIMARY ;
            else if (flag == 1)
                _flag.at(q).at(i) = SECONDARY ;
            else
                _flag.at(q).at(i) = CRITICAL ;
        }
    }
    // Load reduction matrices
    // F
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::read_matrix(_F_row.at(q), in) ;
    }
    // G
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::read_matrix(_G_col.at(q), in) ;
    }
    // H
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::read_matrix(_H_col.at(q), in) ;
    }
    // DD
    for (int q=0; q<=_K.dim(); ++q)
    {
        OSM::read_matrix(_DD_col.at(q), in) ;
    }
    return in ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_HDVF_CORE_H
