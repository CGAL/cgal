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

#ifndef CGAL_HDVF_HDVF_DUALITY_H
#define CGAL_HDVF_HDVF_DUALITY_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <CGAL/OSM/OSM.h>
#include <CGAL/OSM/Bitboard.h>
#include <CGAL/HDVF/Hdvf_core.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/HDVF/Sub_sparse_matrix.h>

namespace CGAL {
namespace HDVF {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf_duality` is the implementation of homological discrete vector fields (HDVF for short) for Alexander duality computation.

 \warning The ring of coefficients provided should be a **field**.

 In dimension \f$n\f$, given a complex \f$L\f$ homeomorphic to \f$\mathcal S^n\f$ and a sub-complex \f$K\subseteq L\f$, Alexander duality states that for all \f$q\leqslant n\f$:
 \f\[\tilde H_q(K) \simeq \tilde H^{n-q-1}(L-K)\f\]
 where \f$\tilde H_q\f$ and \f$\tilde H^q\f$ denote reduced homology and cohomology groups.

 In [Gonzalez and al. 2025], the authors prove that, even if \f$L-K\f$ is not a sub-complex, it produces a valid chain complex (we call "co-complex" such a complementary of sub-complex). Hence, its homology/cohomology can be computed and for all \f$q\leqslant n\f$:
 \f\[\tilde H_q(K) \simeq \tilde H^{q+1}(L-K)\f\]

 HDVFs provide a fast and convenient mean to compute this isomorphism. In order to work with convenient finite complexes, the complex \f$L\f$ must be homeomorphic to a ball of dimension \f$n\f$ (thus \f$\mathcal S^n\f$ is actually homeomorphic to \f$L\f$ plus an infinite \f$n\f$-cell closing its boundary).

 Perfect HDVFs are first computed over \f$K\f$ and \f$L-K\f$ (providing corresponding relative homology) respectively and Alexander isomorphism gives rise to a pairing between critical cells in \f$K\f$ and \f$L-K\f$, that is a pairing between homology/cohomology generators in \f$K\f$ and \f$L-K\f$.

 The class provides HDVF constuction operations: `compute_perfect_hdvf()` and `compute_rand_perfect_hdvf()`, which build perfect HDVFs over \f$K\f$ and \f$L-K\f$ respectively.
 Then, `compute_alexander_pairing()` computes Alexander isomorphism (and provides a pairing between homology/cohomology generators in \f$K\f$ and \f$L-K\f$).


 <img src="HDVF_twirl_view1.png" align="center" width=35%/>
 <img src="HDVF_twirl_view2.png" align="center" width=30%/>

 Example of Alexander duality isomorphism. The twirl mesh is a subcomplex `K` of a larger complex `L` depicted in yellow, homeomorphic to the ball of dimension 3 (right - sectional view).

 \cgalFigureBegin{Duality_quartet,HDVF_twirl_quartet.png}
 Example of "homological quartet for the twirl model". <B>1:</B> Homology generators of the twirl \f$H_1(K)\f$, <B>2:</B> Cohomology generators of the twirl \f$H^1(K)\f$, <B>3:</B> Homology generators of the complementary of the twirl \f$H_1(L-K)\f$, <B>4:</B> Cohomology generators of the complementary of the twirl \f$H^1(L-K)\f$. Alexander isomorphism is represented through colours (paired generators have similar colours).
 \cgalFigureEnd

 Hence, each hole in \f$K\f$ gives rise to four generators (called its "homological quarted": its homology and cohomology generators in \f$K_q\f$ and the homology and cohomology generators paired with them in \f$(L-K)_{q+1}\f$).

 In order to compute relative homology, a sub chain complex mask is used to partially screen the complex `L` and thus restrict HDVF computation. This mask is called "current mask" (and can be set over `K` or `L-K`).

 \cgalModels{HDVF}

 \tparam CoefficientType a model of the `Ring` concept  providing the ring used to compute homology.
 \tparam ComplexType a model of the `AbstractChainComplex` concept, providing the type of abstract chain complex used.

 [Gonzalez and al. 2025] Gonzalez-Lorenzo, A., Bac, A. & Gazull, YS. A constructive approach of Alexander duality. J Appl. and Comput. Topology 9, 2 (2025).
 */

template<typename CoefficientType, typename ComplexType>
class Hdvf_duality : public Hdvf_core<CoefficientType, ComplexType, OSM::Sparse_chain, OSM::Sub_sparse_matrix> {
private:
     // Type of column-major chains
    typedef CGAL::OSM::Sparse_chain<CoefficientType, CGAL::OSM::COLUMN> Col_chain;
     // Type of row-major chains
    typedef CGAL::OSM::Sparse_chain<CoefficientType, CGAL::OSM::ROW> Row_chain;

    // Type of parent HDVF
    typedef Hdvf_core<CoefficientType, ComplexType, CGAL::OSM::Sparse_chain, CGAL::OSM::Sub_sparse_matrix> HDVF_type ;

    // Complex L
    const ComplexType& _L ;
    const int _hdvf_opt ;
    // Subcomplex K
    // _KCC is the Sub_chain_complex_mask describing the subcomplex K
    // _subCC is the Sub_chain_complex_mask describing the current subcomplex (K, L-K or remaining critical cells for pairing)
    Sub_chain_complex_mask<CoefficientType,ComplexType> _KCC, _subCC ;

    // Critical cells of perfect HDVFs (over K / L-K respectively)
    std::vector<std::vector<size_t> > _critical_K, _critical_L_K ;

public:
    /**
     * \brief Hdvf_duality constructor ( from a complex `L` and a sub-complex `K`)
     *
     * `L` is a complex of a given dimension \f$n\f$ homeomorphic to \f$\mathcal B^n\f$ and `K` is a sub-complex of `L` described by a bitboard (cells of `K` have a bit set to 1, cells of `K` have a bit set to 0).
     *
     * Initially, the sub chain complex mask is set to `K`.
     *
     * \param[in] L A complex of a given dimension \f$n\f$ homeomorphic to \f$\mathcal B^n\f$.
     * \param[in] K A sub complex of `L` encoded through a bitboard.
     * \param[in] hdvf_opt Option for HDVF computation (`OPT_BND`, `OPT_F`, `OPT_G` or `OPT_FULL`).
     */
    Hdvf_duality(const ComplexType& L, Sub_chain_complex_mask<CoefficientType, ComplexType>& K, int hdvf_opt = OPT_FULL) ;

    /**
     * \brief Finds a valid Pair_cells of dimension q / q+1 for A *in the current sub chain complex*.
     *
     * The function searches a pair of critical cells, *in the current sub chain complex*, \f$(\gamma_1, \gamma2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Lower dimension of the pair.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual Pair_cells find_pair_A(int q, bool &found) const;

    /**
     * \brief Finds a valid Pair_cells for A containing `gamma` *in the current sub chain complex* (a cell of dimension `q`)
     *
     * The function searches a cell \f$\gamma'\f$ *in the current sub chain complex* such that one of the following conditions holds:
     * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
     * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
     *
     * \param[in] q Dimension of the cell `gamma`.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] gamma Index of a cell to pair.
     */
    virtual Pair_cells find_pair_A(int q, bool &found, size_t gamma) const;

    /**
     * \brief Finds *all* valid Pair_cells of dimension q / q+1 *in the current sub chain complex* for A.
     *
     * The function searches all pairs of critical cells \f$(\gamma_1, \gamma2)\f$ *in the current sub chain complex* of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Lower dimension of the pair.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual std::vector<Pair_cells> find_pairs_A(int q, bool &found) const;

    /**
     * \brief Finds *all* valid Pair_cells for A containing `gamma` *in the current sub chain complex* (a cell of dimension `q`)
     *
     * The function searches all `CRITICAL` cells \f$\gamma'\f$ *in the current sub chain complex* such that one of the following conditions holds:
     * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
     * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the cell `gamma`.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] gamma Index of a cell to pair.
     */
    virtual std::vector<Pair_cells> find_pairs_A(int q, bool &found, size_t gamma) const;

    /** \brief Sets the current sub chain complex masks over `K`.
     *
     * Further HDVF computations will be restricted to `K` (ie.\ computation of reduced homology).
     */
    inline void set_mask_K ()
    {
        _subCC = _KCC ;
        _subCC.screen_matrices(this->_DD_col);
    }

    /** \brief Sets the current sub chain complex masks over `L-K`.
     *
     * Further HDVF computations will be restricted to `L-K` (ie.\ computation of reduced homology).
     */
    inline void set_mask_L_K ()
    {
        _subCC = _KCC.complement() ;
        _subCC.screen_matrices(this->_DD_col);
    }

    /** \brief Returns the value of the current sub chain complex mask.
     */
    Sub_chain_complex_mask<CoefficientType,ComplexType> get_current_mask()
    {
        return _subCC;
    }

    /**
     * \brief Computes a perfect HDVF over the current sub chain complex.
     *
     * As long as valid pairs for A exist in the current sub chain complex, the function selects the first available pair (returned by `find_pair_A()`) and applies the corresponding `A()` operation.
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
     * \brief Computes a random perfect HDVF over the current sub chain complex.
     *
     * As long as valid pairs for A exist in the current sub chain complex, the function selects a random pair (among pairs returned by `find_pairs_A()`) and applies the corresponding `A()` operation.
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
     * \brief Computes a "pairing" HDVF between K and L-K
     *
     * \warning Run `compute_perfect_hdvf()` first (to build perfect HDVFs over `K` and `L-K` respectively).
     *
     * The function computes a perfect HDVF over remaining critical cells. Each pair of cells inserted with the `A()` operation maps corresponding homology/cohomology generators in the Alexander isomorphism.
     *
     * \return The vector of paired critical cells (encoding Alexander isomorphism).
     */
    std::vector<Pair_cells> compute_pairing_hdvf() ;

    /**
     * \brief Computes a random "pairing" HDVF between K and L-K
     *
     * \warning Run `compute_perfect_hdvf()` first (to build perfect HDVFs over `K` and `L-K` respectively).
     *
     * The function computes a random perfect HDVF over remaining critical cells. Each pair of cells inserted with the `A() operation maps corresponding homology/cohomology generators in the Alexander isomorphism.
     *
     * \return The vector of paired critical cells (encoding Alexander isomorphism).
     */
    std::vector<Pair_cells> compute_rand_pairing_hdvf() ;

    // Hdvf_duality getters
    /**
     * \brief Gets cells with a given `flag` in any dimension *in the current sub chain complex*.
     *
     * The function returns a vector containing, for each dimension, the vector of cells with a given `flag`.
     *
     * \param[in] flag Flag to select.
     */
    vector<vector<size_t> > flag (FlagType flag) const ;

    /**
     * \brief Gets cells with a given `flag` in dimension `q` *in the current sub chain complex*.
     *
     * The function returns the vector of cells of dimension `q` with a given `flag`.
     *
     * \param[in] flag Flag to select.
     * \param[in] q Dimension visited.
     */
    vector<size_t> flag_dim (FlagType flag, int q) const ;

    // Hdvf_duality I/O

    /**
     * \brief Prints the homology and cohomology reduction information for `K` and `L-K`.
     *
     * Prints \f$f^*\f$, \f$g\f$ \f$\partial'\f$ the reduced boundary over each critical cell.
     *
     * By default, outputs the complex to `std::cout`.
    */
    virtual ostream& insert_reduction(ostream& out = cout)
    {
        // Print K
        out << "----> K" << endl ;
        _subCC = _KCC ;
        _subCC.screen_matrices(this->_DD_col);
        print_reduction_sub(out) ;

        // Print L-K
        _subCC = _KCC.complement() ;
        _subCC.screen_matrices(this->_DD_col);
        print_reduction_sub(out) ;

        // Set back _subCC to K
        _subCC = _KCC ;

        return out ;
    }

    /**
     * \brief Prints the homology and cohomology reduction information for the current such chain complex.
     *
     * Prints \f$f^*\f$, \f$g\f$ \f$\partial'\f$ the reduced boundary over each critical cell.
     *
     * By default, outputs the complex to `std::cout`.
    */
    virtual ostream& print_reduction_sub(ostream& out = cout) // const;
    {
        // Print critical cells
        out << "----- critical cells:" << endl;
        for (int q = 0; q <= _L.dim(); ++q) {
            out << "--- dim " << q << endl;
            for (size_t i = 0; i < _L.nb_cells(q); ++i)
            {
                if ((this->_flag[q][i] == CRITICAL) && (_subCC.get_bit(q, i))) {
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
                for (size_t i = 0; i < _L.nb_cells(q); ++i) {
                    if ((this->_flag[q][i] == CRITICAL) && (_subCC.get_bit(q, i))) {
                        out << "g(" << i << ") = (" << i << ")";
                        // Iterate over the ith column of _G_col
                        typename HDVF_type::Col_chain col(OSM::get_column(this->_G_col.at(q), i)) ; // TODO cget
                        for (typename HDVF_type::Col_chain::const_iterator it_col = col.cbegin(); it_col != col.cend(); ++it_col) {
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
                for (size_t i = 0; i < _L.nb_cells(q); ++i) {
                    if ((this->_flag[q][i] == CRITICAL) && (_subCC.get_bit(q, i))) {
                        out << "f*(" << i << ") = (" << i << ")";
                        // Iterate over the ith row of _F_row
                        typename HDVF_type::Row_chain row(OSM::get_row(this->_F_row.at(q), i)) ; // TODO cget
                        for (typename HDVF_type::Row_chain::const_iterator it_row = row.cbegin(); it_row != row.cend(); ++it_row) {
                            out << " + " << it_row->second << ".(" << it_row->first << ") + ";
                        }
                        out << endl;
                    }
                }
            }
        }
        return out ;
    }

    /**
     * \brief Prints the reduced boundary over critical cells of `K` and `L-K`.
     *
     * The method prints out the reduced boundary matrix in each dimension, restricted to critical cells of `K` and `L-K` (ie.\ the matrix used to compute Alexander pairing).
     *
     * \warning Call this method after `compute_perfect_hdvf()`.
     *
     * By default, outputs the complex to `std::cout`.
    */
    ostream& print_bnd_pairing(ostream& out = cout)
    {
        Sub_chain_complex_mask<CoefficientType, ComplexType> subPair(_L, false) ;
        for (int q=0; q<=_L.dim(); ++q)
        {
            for (size_t i=0; i<_critical_K.at(q).size(); ++i)
                subPair.set_bit_on(q, _critical_K.at(q).at(i)) ;
            for (size_t i=0; i<_critical_L_K.at(q).size(); ++i)
                subPair.set_bit_on(q, _critical_L_K.at(q).at(i)) ;
        }
        // Print corresponding submatrices _DD_col
        for (int q=1; q<=_L.dim(); ++q)
        {
            out << "--> dim " << q << " : q / q-1 cells" << endl ;
            out << "id " << q << " : " ;
            for (typename OSM::Bitboard::iterator it = subPair.get_bitboard(q).begin(); it != subPair.get_bitboard(q).end(); ++it)
                out << *it << " " ;
            out << endl ;
            out << "id " << q-1 << " : " ;
            for (typename OSM::Bitboard::iterator it = subPair.get_bitboard(q-1).begin(); it != subPair.get_bitboard(q-1).end(); ++it)
                out << *it << " " ;
            out << endl ;
            for (typename OSM::Bitboard::iterator it = subPair.get_bitboard(q).begin(); it != subPair.get_bitboard(q).end(); ++it)
            {
                for (typename OSM::Bitboard::iterator it2 = subPair.get_bitboard(q-1).begin(); it2 != subPair.get_bitboard(q-1).end(); ++it2)
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

    /**
     * \brief Exports primary/secondary/critical labels *of the current sub chain complex* for vtk export.
     *
     * The method exports the labels of every cells in each dimension.
     *
     * \return A vector containing, for each dimension, the vector of labels by cell index.
     */
    virtual vector<vector<int> > psc_labels () const
    {
        vector<vector<int> > labels(this->_K.dim()+1) ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            for (size_t i = 0; i<this->_K.nb_cells(q); ++i)
            {
                if (_subCC.get_bit(q, i)) // i belongs to _subCC
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

    /**
     * \brief Exports homology generators *of the current sub chain complex* associated to `cell` (critical cell) of dimension  `q` (used by vtk export).
     *
     * The method exports the chain \f$g(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell` and dimension `q`.
     *
     * \return A column-major chain.
     */
    virtual Col_chain homology_chain (size_t cell, int dim) const
    {
        if ((dim<0) || (dim>this->_K.dim()))
            throw "Error : homology_chain with dim out of range" ;
        if (!_subCC.get_bit(dim, cell))
            throw "Error : homology_chain for a cell out of current sub chain complex" ;

        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            // Get g(cell, dim) with per indices
            Col_chain g_cell(OSM::get_column(this->_G_col.at(dim), cell)) ;
            // Add 1 to the cell
            g_cell.set_coef(cell, 1) ;
            // Keep cells of the chain belonging to _subCC
            Col_chain g_cell_sub(g_cell.dimension()) ;
            for (typename Col_chain::const_iterator it = g_cell.begin(); it != g_cell.end(); ++it)
            {

                if (_subCC.get_bit(dim, it->first))
                {
                    g_cell_sub.set_coef(it->first, it->second) ;
                }
            }
            return g_cell_sub ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }

    /**
     * \brief Exports cohomology generators *of the current sub chain complex* associated to `cell` (critical cell) of dimension  `q` (used by vtk export).
     *
     * The method exports the chain \f$f^\star(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell` and dimension `q`.
     *
     * \return A column-major chain.
     */
    virtual Col_chain cohomology_chain (size_t cell, int dim) const
    {
        if ((dim<0) || (dim>this->_K.dim()))
            throw "Error : cohomology_chain with dim out of range" ;
        if (!_subCC.get_bit(dim, cell))
            throw "Error : cohomology_chain for a cell out of current sub chain complex" ;

        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            Row_chain fstar_cell(OSM::get_row(this->_F_row.at(dim), cell)) ;
            // Add 1 to the cell
            fstar_cell.set_coef(cell, 1) ;

            // Keep cells of the chain belonging to _subCC
            Col_chain fstar_cell_sub(fstar_cell.dimension()) ;
            for (typename Col_chain::const_iterator it = fstar_cell.begin(); it != fstar_cell.end(); ++it)
            {

                if (_subCC.get_bit(dim, it->first))
                {
                    fstar_cell_sub.set_coef(it->first, it->second) ;
                }
            }
            return fstar_cell_sub ;
        }
    }
} ;

// Constructor
template<typename CoefficientType, typename ComplexType>
Hdvf_duality<CoefficientType,ComplexType>::Hdvf_duality(const ComplexType& L, Sub_chain_complex_mask<CoefficientType, ComplexType>& K, int hdvf_opt) :
Hdvf_core<CoefficientType, ComplexType, OSM::Sparse_chain, OSM::Sub_sparse_matrix>(L,hdvf_opt), _L(L), _hdvf_opt(hdvf_opt), _KCC(K), _subCC(K) {}

// find a valid Pair_cells for A in dimension q
template<typename CoefficientType, typename ComplexType>
Pair_cells Hdvf_duality<CoefficientType,ComplexType>::find_pair_A(int q, bool &found) const
{
    found = false;
    Pair_cells p;

    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = this->_DD_col[q+1].begin(); (it_col != this->_DD_col[q+1].end() && !found); ++it_col)
    {
        const typename HDVF_type::Col_chain& col(OSM::cget_column(this->_DD_col[q+1], *it_col)) ;

        // Iterate through the entries of the column
        // Check that the row belongs to the subchaincomplex
        for (typename HDVF_type::Col_chain::const_iterator it = col.begin(); (it != col.end() && !found); ++it) {
            if (_subCC.get_bit(q, it->first) && (abs(it->second) == 1)) {
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
template<typename CoefficientType, typename ComplexType>
Pair_cells Hdvf_duality<CoefficientType,ComplexType>::find_pair_A(int q, bool &found, size_t tau) const
{
    found = false;
    Pair_cells p ;
    // Check tau belongs to _subCC
    if (!_subCC.get_bit(q, tau))
        throw("Hdvf_duality: searching for a cell tau outside _subCC") ;

    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    // and tau' belongs to _subCC
    const typename HDVF_type::Col_chain& tmp2(OSM::cget_column(this->_DD_col.at(q), tau)) ;
    for (typename HDVF_type::Col_chain::const_iterator it = tmp2.cbegin(); (it != tmp2.cend() && !found); ++it)
    {
        if (_subCC.get_bit(q-1, it->first) && abs(it->second) == 1)
        {
            found = true ;
            p.sigma = it->first ;
            p.tau = tau ;
            p.dim = q-1 ;
        }
    }

    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    // and tau' belongs to _subCC
    typename HDVF_type::Row_chain tmp(OSM::get_row(this->_DD_col.at(q+1), tau)) ;
    for (typename HDVF_type::Row_chain::const_iterator it = tmp.cbegin(); (it != tmp.cend() && !found); ++it)
    {
        if (_subCC.get_bit(q+1, it->first) && (abs(it->second) == 1))
        {
            found = true ;
            Pair_cells p ;
            p.sigma = tau ;
            p.tau = it->first ;
            p.dim = q ;
        }
    }
    return p;
}

// find all the valid Pair_cells for A in dimension q
template<typename CoefficientType, typename ComplexType>
std::vector<Pair_cells> Hdvf_duality<CoefficientType,ComplexType>::find_pairs_A(int q, bool &found) const
{
    vector<Pair_cells> pairs;
    found = false ;

    // Iterate through columns of _DD_col[q+1]
    for (OSM::Bitboard::iterator it_col = this->_DD_col[q+1].begin(); it_col != this->_DD_col[q+1].end(); ++it_col)
    {
        const typename HDVF_type::Col_chain& col(OSM::cget_column(this->_DD_col[q+1], *it_col)) ;

        // Iterate through the entries of the column
        for (typename HDVF_type::Col_chain::const_iterator it = col.begin(); it != col.end(); ++it) {
            if (_subCC.get_bit(q, it->first) && ((it->second == 1) || (it->second == -1))) {
                // If an entry of _subCC with coefficient 1 or -1 is found, set the pair and mark as found
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

// find all the valid Pair_cells containing tau for A in dimension q
template<typename CoefficientType, typename ComplexType>
std::vector<Pair_cells> Hdvf_duality<CoefficientType,ComplexType>::find_pairs_A(int q, bool &found, size_t tau) const
{
    found = false;
    std::vector<Pair_cells> pairs;
    // Check if tau belongs to _subCC
    if (!_subCC.get_bit(q, tau))
        throw("Hdvf_duality: searching for a cell tau outside _subCC") ;

    // Search for a q+1 cell tau' such that <_d(tau'),tau> invertible, ie <_cod(tau),tau'> invertible
    // and tau' belongs to _subCC
    typename HDVF_type::Row_chain tmp(OSM::get_row(this->_DD_col.at(q+1), tau)) ;
    for (typename HDVF_type::Row_chain::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
    {
        if (_subCC.get_bit(q+1, it->first) && (abs(it->second) == 1))
        {
            found = true ;
            Pair_cells p ;
            p.sigma = tau ;
            p.tau = it->first ;
            p.dim = q ;
            pairs.push_back(p) ;
        }
    }
    // Search for a q-1 cell tau' such that <_d(tau),tau'> invertible
    // and tau' belongs to _subCC
    const typename HDVF_type::Col_chain& tmp2(OSM::cget_column(this->_DD_col.at(q), tau)) ;
    for (typename HDVF_type::Col_chain::const_iterator it = tmp2.cbegin(); it != tmp2.cend(); ++it)
    {
        if (_subCC.get_bit(q-1, it->first) && (abs(it->second) == 1))
        {
            found = true ;
            Pair_cells p ;
            p.sigma = it->first ;
            p.tau = tau ;
            p.dim = q-1 ;
            pairs.push_back(p) ;
        }
    }
    return pairs;
}

// Compute dual perfect HDVFs (over K and L-K)
template<typename CoefficientType, typename ComplexType>
std::vector<Pair_cells> Hdvf_duality<CoefficientType,ComplexType>::compute_perfect_hdvf(bool verbose)
{
    std::cout << endl << "==== Compute perfect HDVF over K" << endl ;
    // Set _subCC to K
    _subCC = _KCC ;
    // Restrict _DD_col accordingly
    _subCC.screen_matrices(this->_DD_col);
    // Compute perfect HDVF over K
    std::vector<Pair_cells> tmp = HDVF_type::compute_perfect_hdvf(verbose) ;
    std::cout << tmp.size() << " cells paired" << endl ;
    _critical_K = flag(CRITICAL) ;

    std::cout << endl << "==== Compute perfect HDVF over L-K" << endl ;
    // set _subCC to L-K
    _subCC = _KCC.complement() ;
    // Restrict _DD_col accordingly
    _subCC.screen_matrices(this->_DD_col);
    // Compute perfect HDVF over L-K
    std::vector<Pair_cells> tmp2 = HDVF_type::compute_perfect_hdvf(verbose) ;
    std::cout << tmp2.size() << " cells paired" << endl ;
    _critical_L_K = flag(CRITICAL) ;

    // Return the vector of paired cells
    tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
    return tmp;
}

// Compute random dual perfect HDVFs (over K and L-K)
template<typename CoefficientType, typename ComplexType>
std::vector<Pair_cells> Hdvf_duality<CoefficientType,ComplexType>::compute_rand_perfect_hdvf(bool verbose)
{
    std::cout << endl << "==== Compute perfect HDVF over K" << endl ;
    // Set _subCC to K
    _subCC = _KCC ;
    // Restrict _DD_col accordingly
    _subCC.screen_matrices(this->_DD_col);
    // Compute perfect HDVF over K
    std::vector<Pair_cells> tmp = HDVF_type::compute_rand_perfect_hdvf(verbose) ;
    std::cout << tmp.size() << " cells paired" << endl ;
    _critical_K = flag(CRITICAL) ;

    std::cout << endl << "==== Compute perfect HDVF over L-K" << endl ;
    // set _subCC to L-K
    _subCC = _KCC.complement() ;
    // Restrict _DD_col accordingly
    _subCC.screen_matrices(this->_DD_col);
    // Compute perfect HDVF over L-K
    std::vector<Pair_cells> tmp2 = HDVF_type::compute_rand_perfect_hdvf(verbose) ;
    std::cout << tmp2.size() << " cells paired" << endl ;
    _critical_L_K = flag(CRITICAL) ;

    // Return the vector of paired cells
    tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
    return tmp;
}

// Compute Alexander isomorphism (A pairing between critical cells of K / L-K)
template<typename CoefficientType, typename ComplexType>
std::vector<Pair_cells> Hdvf_duality<CoefficientType,ComplexType>::compute_pairing_hdvf()
{
    // TODO : check both HDVFs are perfect

    std::cout << endl << "==== Compute pairing" << endl ;

    // Create a full Sub_chain_complex_mask
    _subCC = Sub_chain_complex_mask<CoefficientType, ComplexType>(_L) ;
    _subCC.screen_matrices(this->_DD_col);
    // If necessary, copy the HDVF before computing the pairing -> otherwise we loose it...
    std::vector<Pair_cells> pairing = HDVF_type::compute_perfect_hdvf() ;
    return pairing ;
}

// Compute random Alexander isomorphism (A pairing between critical cells of K / L-K)
template<typename CoefficientType, typename ComplexType>
std::vector<Pair_cells> Hdvf_duality<CoefficientType,ComplexType>::compute_rand_pairing_hdvf()
{
    // TODO : check both HDVFs are perfect

    std::cout << endl << "==== Compute pairing" << endl ;

    // Create a full Sub_chain_complex_mask
    _subCC = Sub_chain_complex_mask<CoefficientType, ComplexType>(_L) ;
    _subCC.screen_matrices(this->_DD_col);
    // If necessary, copy the HDVF before computing the pairing -> otherwise we loose it...
    std::vector<Pair_cells> pairing = HDVF_type::compute_rand_perfect_hdvf() ;
    return pairing ;
}

// Method to get cells of _subCC with a given flag for each dimension
template<typename CoefficientType, typename ComplexType>
vector<vector<size_t> > Hdvf_duality<CoefficientType,ComplexType>::flag (FlagType flag) const
{
    vector<vector<size_t> > res(_L.dim()+1) ;
    for (int q=0; q<=_L.dim(); ++q)
    {
        for (size_t i=0; i<_L.nb_cells(q); ++i)
        {
            if (_subCC.get_bit(q, i) && (this->_flag.at(q).at(i) == flag))
                res.at(q).push_back(i) ;
        }
    }
    return res ;
}

// Method to get cells of _subCC with a given flag for a given dimension
template<typename CoefficientType, typename ComplexType>
vector<size_t> Hdvf_duality<CoefficientType,ComplexType>::flag_dim (FlagType flag, int q) const
{
    vector<size_t> res ;
    for (size_t i=0; i<this->_K.nb_cells(q); ++i)
    {
        if (_subCC.get_bit(q, i) && (this->_flag.at(q).at(i) == flag))
            res.push_back(i) ;
    }
    return res ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_HDVF_DUALITY_H
