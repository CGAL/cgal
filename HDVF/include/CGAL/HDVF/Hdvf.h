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

#ifndef CGAL_HDVF_HDVF_H
#define CGAL_HDVF_HDVF_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <CGAL/HDVF/Hdvf_core.h>
#include <CGAL/OSM/OSM.h>
#include <CGAL/OSM/Bitboard.h>

namespace CGAL {
namespace Homological_discrete_vector_field {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf` implements homology and cohomology computation via homological discrete vector fields (%HDVF for short). It derives from `Hdvf_core` and shares all its data and methods.

 But besides construction operations and methods (using the `A()` operation), the `Hdvf` class implements four other HDVF operations: R, M, W and MW together with appropriate "find_pair()" functions. These operations change the %HDVF (that is change homology / cohomology generators) and thus provide a convenient tool to move inside the "space of homology/cohomology computations".

- `R()` operation is the "dual" of the A pairing operation (it cancels the pairing and turns back a `PRIMARY`/`SECONDARY` pair into a pair of `CRITICAL` cells)
- `M()` operation exchanges a `PRIMARY` \f$\pi\f$ and a `CRITICAL` cell \f$\gamma\f$ (under conditions) and modifies the homology generator associated to \f$\gamma\f$ (while preserving is associated cohomology generator)
- `W()` operation exchanges a `SECONDARY` \f$\sigma\f$ and a `CRITICAL` cell \f$\gamma\f$ (under conditions) and modifies the cohomology generator associated to \f$\gamma\f$ (while preserving is associated homology generator)
- `MW()` operation exchanges a `PRIMARY` \f$\pi\f$ and a `SECONDARY` cell \f$\sigma\f$ (under conditions). See the introduction to HDVF for more details on this operation.

Using appropriate combinations of such operations, one can change a HDVF until corresponding homology or cohomology generators meet a given basis or delineate a hole.

 Let us consider the following simple cubical complex and a perfect HDVF (top) together with the three corresponding homology generators (bottom, highlighted in pink):

 <img src="HDVF_op1.png" align="center" width=50%/><br>
 <img src="HDVF_op1_g1.png" align="center" width=30%/>
 <img src="HDVF_op1_g2.png" align="center" width=30%/>
 <img src="HDVF_op1_g3.png" align="center" width=30%/>

 A W operation between cells of Khalimsky coordinates \f$\sigma = (3,2)\f$ (`CRITICAL`) and \f$\tau=(5,2)\f$ (`SECONDARY`) produces the following HDVF (homology generators did not change; note that cohomology generators are modified):

 <img src="HDVF_op2_W7_11.png" align="center" width=50%/>

 <img src="HDVF_op2_g1.png" align="center" width=30%/>
 <img src="HDVF_op2_g2.png" align="center" width=30%/>
 <img src="HDVF_op2_g3.png" align="center" width=30%/>

 Then, a MW operation between cells of Khalimsky coordinates \f$\sigma' = (4,1)\f$ (`PRIMARY`) and \f$\tau=(3,2)\f$ (`SECONDARY`) produces the following HDVF where homology generators become "minimal":

 <img src="HDVF_op3_MW5_7.png" align="center" width=50%/>

 <img src="HDVF_op3_g1.png" align="center" width=30%/>
 <img src="HDVF_op3_g2.png" align="center" width=30%/>
 <img src="HDVF_op3_g3.png" align="center" width=30%/>

 Perfect HDVFs provide various topological results:
 - Betti numbers: are well known topological descriptors (providing the number of holes in each dimension). For a perfect HDVF, Betti numbers equal the number of critical cells (that can be obtained through `psc_flags()` with the `CRITICAL``PSC_flag` as argument).
 - Homology/cohomology generators are actually algebraic objects, namely chains. Methods `homology_chain()` and `cohomology_chain()` return the homology and cohomology generator chain associated to a given critical cell. VTK export functions output all the cells of such chains with non zero coefficients. Figures here above illustrate such homology and co-homology generators.
 - Homology/cohomology annotation: use `get_annotation()` or `get_coannotation()` to get the annotation/co-annotation of a cycle/co-cycle in the homology/cohomology basis (as a chain of critical cells).
 - Homology/cohomology comparison: use `are_same_cycles()` or `are_same_cocycles()` to check if two cycles (or co-cycles) belong to the same homology / cohomology class.

 \cgalFigureBegin{HDVFannotation,HDVF_annotations.png}
 Illustration of cycles annotation and comparison. <B>Left:</B> Generator \f$g(\sigma)\f$ associated to the critical cell \f$\sigma\f$; <B>Middle:</B> A cycle \f$\alpha\f$ ; <B>Right:</B> A second cycle \f$\beta\f$.

 For \f$\mathbb Z_2\f$ homology, the annotation of \f$\alpha\f$ is \f$\sigma\f$ (\f$\alpha\f$ equals \f$g(\sigma)\f$ up to a boundary), while the annotation of \f$\beta\f$ is \f$\sigma+\tau+\gamma\f$ (\f$\beta\f$ equals \f$g(\sigma)+g(\tau)+g(\gamma)\f$ up to a boundary).

 Therefore, `are_same_cycles` will return `true` for \f$\alpha\f$  and \f$g(\sigma)\f$,  but `false` for \f$\beta\f$  and \f$g(\sigma)\f$.
 \cgalFigureEnd

\cgalModels{HDVF}

\tparam ChainComplex a model of the `AbstractChainComplex` concept, providing the type of abstract chain complex used.
 */

template<typename ChainComplex>
class Hdvf : public Hdvf_core<ChainComplex, OSM::Sparse_chain, OSM::Sparse_matrix> {
public:
    /*! \brief Type of coefficients used to compute homology. */
    typedef ChainComplex::Coefficient_ring Coefficient_ring;
    
    /*!
     Type of parent Hdvf_core class.
     */
    typedef Hdvf_core<ChainComplex, OSM::Sparse_chain, OSM::Sparse_matrix> HDVF_parent ;

    /**
     * \brief Default constructor.
     *
     * Builds an" empty" HDVF associated to K (with all cells critical). By default, the HDVF option is set to OPT_FULL (full reduction computed).
     *
     * \param[in] K A chain complex (a model of `AbstractChainComplex`)
     * \param[in] hdvf_opt Option for HDVF computation (`OPT_BND`, `OPT_F`, `OPT_G` or `OPT_FULL`)
     */
    Hdvf(const ChainComplex& K, int hdvf_opt = OPT_FULL) ;

    /*
     * \brief Constructor by copy.
     *
     * Builds a HDVF by copy from another, including options.
     *
     * \param[in] hdvf An initial HDVF.
     */
    Hdvf(const Hdvf& hdvf) : Hdvf_core<ChainComplex, OSM::Sparse_chain, OSM::Sparse_matrix>(hdvf) { }

    /*
     * \brief HDVF destructor. */
    ~Hdvf() { }

    // findPair functions for M

    /**
     * \brief Finds a valid Cell_pair of dimension q for M.
     *
     * The function searches a pair of cells \f$(\pi, \gamma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\gamma\f$ `CRITICAL`, valid for M (ie.\ such that \f$\langle f(\pi), \gamma \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    Cell_pair find_pair_M(int q, bool &found) const;

    /**
     * \brief Finds a valid Cell_pair of dimension q for M cointaining `tau`.
     *
     * The function searches a pair of cells \f$(\pi, \gamma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\gamma\f$ `CRITICAL` (one of them is `tau`), valid for M (ie.\ such that \f$\langle f(\pi), \gamma \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] tau Cell of dimension `q` to pair.
     */
    Cell_pair find_pair_M(int q, bool &found, size_t tau) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q for M.
     *
     * The function searches all pairs of cells \f$(\pi, \gamma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\gamma\f$ `CRITICAL`, valid for M (ie.\ such that \f$\langle f(\pi), \gamma \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    std::vector<Cell_pair> find_pairs_M(int q, bool &found) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q for M cointaining `tau`.
     *
     * The function searches all pairs of cells \f$(\pi, \gamma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\gamma\f$ `CRITICAL` (one of them is `tau`), valid for M (ie.\ such that \f$\langle f(\pi), \gamma \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] tau Cell of dimension `q` to pair.
     */
    std::vector<Cell_pair> find_pairs_M(int q, bool &found, size_t tau) const;

    // findPair functions for W

    /**
     * \brief Finds a valid Cell_pair of dimension q for W.
     *
     * The function searches a pair of cells \f$(\sigma, \gamma)\f$ with \f$\sigma\f$ `SECONDARY` and \f$\gamma\f$ `CRITICAL`, valid for W (ie.\ such that \f$\langle g(\gamma), \sigma \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    Cell_pair find_pair_W(int q, bool &found) const;

    /**
     * \brief Finds a valid Cell_pair of dimension q for W cointaining `tau`.
     *
     * The function searches a pair of cells \f$(\sigma, \gamma)\f$ with \f$\sigma\f$ `SECONDARY` and \f$\gamma\f$ `CRITICAL` (one of them is `tau`), valid for W (ie.\ such that \f$\langle g(\gamma), \sigma \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] tau Cell of dimension `q` to pair.
     */
    Cell_pair find_pair_W(int q, bool &found, size_t tau) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q for W.
     *
     * The function searches all pairs of cells \f$(\sigma, \gamma)\f$ with \f$\sigma\f$ `SECONDARY` and \f$\gamma\f$ `CRITICAL`, valid for W (ie.\ such that \f$\langle g(\gamma), \sigma \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    std::vector<Cell_pair> find_pairs_W(int q, bool &found) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q for W cointaining `tau`.
     *
     * The function searches all pairs of cells \f$(\sigma, \gamma)\f$ with \f$\sigma\f$ `SECONDARY` and \f$\gamma\f$ `CRITICAL` (one of them is `tau`), valid for W (ie.\ such that \f$\langle g(\gamma), \sigma \rangle\f$ invertible). It returns the first valid pair found by iterators.
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] tau Cell of dimension `q` to pair.
     */
    std::vector<Cell_pair> find_pairs_W(int q, bool &found, size_t tau) const;

    // findPair functions for MW

    /**
     * \brief Finds a valid Cell_pair of dimension q for MW.
     *
     * The function searches a pair of cells \f$(\pi, \sigma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\sigma\f$ `SECONDARY`, valid for MW (ie.\ such that \f$\langle h_{q-1}\partial_q(\pi), \sigma \rangle\f$ invertible and \f$\langle \partial_{q+1} h_q(\sigma), \pi \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    Cell_pair find_pair_MW(int q, bool &found) const;

    /**
     * \brief Finds a valid Cell_pair of dimension q for MW cointaining `tau`.
     *
     * The function searches a pair of cells \f$(\pi, \sigma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\sigma\f$ `SECONDARY` (one of them is `tau`), valid for MW (ie.\ such that \f$\langle h_{q-1}\partial_q(\pi), \sigma \rangle\f$ invertible and \f$\langle \partial_{q+1} h_q(\sigma), \pi \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] tau Cell of dimension `q` to pair.
     */
    Cell_pair find_pair_MW(int q, bool &found, size_t tau) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q for MW.
     *
     * The function searches all pairs of cells \f$(\pi, \sigma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\sigma\f$ `SECONDARY`, valid for MW (ie.\ such that \f$\langle h_{q-1}\partial_q(\pi), \sigma \rangle\f$ invertible and \f$\langle \partial_{q+1} h_q(\sigma), \pi \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    std::vector<Cell_pair> find_pairs_MW(int q, bool &found) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q for W cointaining `tau`.
     *
     * The function searches all pairs of cells \f$(\pi, \sigma)\f$ with \f$\pi\f$ `PRIMARY` and \f$\sigma\f$ `SECONDARY` (one of them is `tau`), valid for MW (ie.\ such that \f$\langle h_{q-1}\partial_q(\pi), \sigma \rangle\f$ invertible and \f$\langle \partial_{q+1} h_q(\sigma), \pi \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * \param[in] q Dimension of the pair searched.
     * \param[in] found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     * \param[in] tau Cell of dimension `q` to pair.
     */
    std::vector<Cell_pair> find_pairs_MW(int q, bool &found, size_t tau) const;

    // Hdvf methods

    /**
     * \brief R operation (cancels a A operation).
     *
     * A pair of cells \f$(\pi, \sigma)\f$ of respective dimension q and q+1, with \f$\pi\f$ `PRIMARY` and \f$\sigma\f$ `SECONDARY`, is valid for R if \f$\langle h(\pi), \sigma \rangle\f$ is invertible. After the R operation, \f$\pi\f$ and \f$\sigma\f$ become `CRITICAL`. The R method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param[in] pi First cell of the pair (dimension `q`)
     * \param[in] sigma Second cell of the pair (dimension `q+1`)
     * \param[in] q Dimension of the pair
     */
    void R(size_t pi, size_t sigma, int q);

    /**
     * \brief M operation.
     *
     * A pair of cells \f$(\pi, \gamma)\f$ of dimension q, with \f$\pi\f$ `PRIMARY` and \f$\gamma\f$ `CRITICAL`, is valid for M if \f$\langle f(\pi), \gamma \rangle\f$ is invertible. After the M operation, \f$\pi\f$ becomes `CRITICAL` and \f$\gamma\f$ become `PRIMARY`. The M method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param[in] pi First cell of the pair (dimension `q`)
     * \param[in] gamma Second cell of the pair (dimension `q`)
     * \param[in] q Dimension of the pair
     */
    void M(size_t pi, size_t gamma, int q);

    /**
     * \brief W operation.
     *
     * A pair of cells \f$(\sigma, \gamma)\f$ of dimension q, with \f$\sigma\f$ `SECONDARY` and \f$\gamma\f$ `CRITICAL`, is valid for W if \f$\langle g(\gamma), \sigma \rangle\f$ is invertible. After the W operation, \f$\sigma\f$ becomes `CRITICAL` and \f$\gamma\f$ become `SECONDARY`. The W method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param[in] sigma First cell of the pair (dimension `q`)
     * \param[in] gamma Second cell of the pair (dimension `q`)
     * \param[in] q Dimension of the pair
     */
    void W(size_t sigma, size_t gamma, int q);

    /**
     * \brief MW operation.
     *
     * A pair of cells \f$(\pi, \sigma)\f$ of dimension q, with \f$\pi\f$ `PRIMARY` and \f$\sigma\f$ `SECONDARY`, is valid for MW if \f$\langle h_{q-1}\partial_q(\pi), \sigma \rangle\f$ is invertible and \f$\langle \partial_{q+1} h_q(\sigma), \pi \rangle\f$ is invertible. After the MW operation, \f$\pi\f$ becomes `SECONDARY` and \f$\sigma\f$ become `PRIMARY`. The MW method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param[in] pi First cell of the pair (dimension `q`)
     * \param[in] sigma Second cell of the pair (dimension `q`)
     * \param[in] q Dimension of the pair
     */
    void MW(size_t pi, size_t sigma, int q);

    /**
     * \brief Gets the annotation of a cycle in the homology basis.
     *
     * The method returns the image of a cycle by the morphism \f$f\f$, that is, a linear combination of critical cells (corresponding to the decomposition of the cycle in the homology basis).
     * If the annotation has a single non zero coefficient for a given critical cell \f$\sigma\f$, then the cycle belongs to the class of the homology generator \f$g(\sigma)\f$.
     *
     * \warning Will raise an error is the chain provided is not a cycle.
     * \warning The HDVF must be perfect.
     *
     * \param[in] chain The cycle to annotate in the homology basis.
     * \param[in] dim Dimension of the cycle.
     */
    typename HDVF_parent::Column_chain get_annotation(typename HDVF_parent::Column_chain chain, int dim) const
    {
        // Check that the chain is a cycle (must belong to the kernel of the boundary operator)
        typename HDVF_parent::Column_chain bnd(this->_DD_col.at(dim) * chain) ;
        if (!bnd.is_null())
            throw("get_annotation: the chain provided is not a cycle");

        // Compute the annotation
        return (this->_F_row.at(dim) * chain + this->projection(chain, CRITICAL, 1)) ;
    }

    /**
     * \brief Gets the co-annotation of a co-cycle in the cohomology basis.
     *
     * The method returns the image of a co-cycle by the morphism \f$g^*\f$, that is, a linear combination of critical cells (corresponding to the decomposition of the co-cycle in the cohomology basis).
     * If the annotation has a single non zero coefficient for a given critical cell \f$\sigma\f$, then the co-cycle belongs to the class of the cohomology generator \f$f^*(\sigma)\f$.
     *
     * \warning Will raise an error is the chain provided is not a co-cycle.
     * \warning The HDVF must be perfect.
     *
     * \param[in] chain The co-cycle to annotate in the homology basis.
     * \param[in] dim Dimension of the co-cycle.
     */
    typename HDVF_parent::Row_chain get_coannotation(typename HDVF_parent::Row_chain chain, int dim) const
    {
        // Check that the chain is a co-cycle (must belong to the kernel of the boundary operator)
        typename HDVF_parent::Row_chain bnd(chain * this->_DD_col.at(dim+1)) ;
        if (!bnd.is_null())
            throw("get_coannotation: the chain provided is not a co-cycle");

        // Compute the co-annotation
        return (chain * this->_G_col.at(dim) + this->projection(chain, CRITICAL, 1)) ;
    }

    /**
     * \brief Checks if two cycles belong to the same homology class.
     *
     * \warning Will raise an error is chains provided are not cycles.
     * \warning The HDVF must be perfect.
     *
     * \param[in] chain1 First cycle.
     * \param[in] chain2 Second cycle.
     * \param[in] dim Dimension of both cycles.
     */
    bool are_same_cycles (typename HDVF_parent::Column_chain chain1, typename HDVF_parent::Column_chain chain2, int dim)
    {
        // Check if both chains are cycles (must belong to the kernel of the boundary operator)
        typename HDVF_parent::Column_chain bnd1(this->_DD_col.at(dim) * chain1), bnd2(this->_DD_col.at(dim) * chain2) ;
        if (!bnd1.is_null())
            throw("get_annotation: chain1 is not a cycle");
        if (!bnd2.is_null())
            throw("get_annotation: chain2 is not a cycle");

        typename HDVF_parent::Column_chain annot1(get_annotation(chain1, dim)), annot2(get_annotation(chain2, dim)) ;
        return annot1 == annot2 ;
    }

    /**
     * \brief Checks if two co-cycles belong to the same cohomology class.
     *
     * \warning Will raise an error is chains provided are not co-cycles.
     * \warning The HDVF must be perfect.
     *
     * \param[in] chain1 First co-cycle.
     * \param[in] chain2 Second co-cycle.
     * \param[in] dim Dimension of both co-cycles.
     */
    bool are_same_cocycles (typename HDVF_parent::Row_chain chain1, typename HDVF_parent::Row_chain chain2, int dim)
    {
        // Check if both chains are cocycles (must belong to the kernel of the boundary^* operator)
        typename HDVF_parent::Row_chain cobnd1(chain1 * this->_DD_col.at(dim)), cobnd2(chain2 * this->_DD_col.at(dim)) ;
        if (!cobnd1.is_null())
            throw("get_coannotation: chain1 is not a co-cycle");
        if (!cobnd2.is_null())
            throw("get_coannotation: chain2 is not a co-cycle");

        typename HDVF_parent::Row_chain coannot1(get_coannotation(chain1, dim)), coannot2(get_coannotation(chain2, dim)) ;
        return coannot1 == coannot2 ;
    }
};

// Constructor for the Hdvf class
template<typename ChainComplex>
Hdvf<ChainComplex>::Hdvf(const ChainComplex& K, int hdvf_opt) : Hdvf_core<ChainComplex, OSM::Sparse_chain, OSM::Sparse_matrix>(K, hdvf_opt) { }



// Methods to find a pair of cells for M in dimension q
//  -> <f(cell1),cell2> = +- 1 with cell1 primary, cell2 critical
// First version: returns a pair of dimensions q
// Second version: returns all the pairs containing sigma

// \brief find a valid Cell_pair for M in dimension q

template<typename ChainComplex>
Cell_pair Hdvf<ChainComplex>::find_pair_M(int q, bool &found) const
{
    found = false;
    Cell_pair p;

    if (this->_hdvf_opt & (OPT_F | OPT_FULL))
    {
        // Search for +-1 in _F - iterate over rows
        for (OSM::Bitboard::iterator it_row = this->_F_row[q].begin(); (it_row != this->_F_row[q].end() && !found); ++it_row)
        {
            const typename HDVF_parent::Row_chain &row(OSM::cget_row(this->_F_row[q],*it_row));

            // Iterate through the entries of the row
            for (typename HDVF_parent::Row_chain::const_iterator it = row.begin(); (it != row.end() && !found); ++it) {
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

// find a valid Cell_pair containing tau for M in dimension q

template<typename ChainComplex>
Cell_pair Hdvf<ChainComplex>::find_pair_M(int q, bool &found, size_t tau) const
{
    found = false;
    Cell_pair p;
    if (this->_flag.at(q).at(tau) == SECONDARY)
        return p ; // Empty / found false

    if (this->_hdvf_opt & (OPT_F | OPT_FULL))
    {
        // If tau is primary, search for gamma such that <f(tau),gamma>=+-1
        if (this->_flag.at(q).at(tau) == PRIMARY)
        {
            // Get the column of tau
            const typename HDVF_parent::Column_chain& col(OSM::get_column(this->_F_row.at(q), tau)) ;
            for (typename HDVF_parent::Column_chain::const_iterator it = col.cbegin(); (it != col.cend() && !found); ++it)
            {
                if (abs(it->second) == 1)
                {
                    p.sigma = tau ;
                    p.tau = it->first ;
                    p.dim = q ;
                    found = true ;
                }
            }
        }
        // If tau is critical, search for pi such that <f(pi),tau>=+-1
        if (this->_flag.at(q).at(tau) == CRITICAL)
        {
            // Search along the row of tau
            const typename HDVF_parent::Row_chain& row(OSM::cget_row(this->_F_row.at(q), tau)) ;
            for (typename HDVF_parent::Row_chain::const_iterator it = row.cbegin(); (it != row.cend() && !found); ++it)
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

// find all the valid Cell_pair for M in dimension q
template<typename ChainComplex>
std::vector<Cell_pair> Hdvf<ChainComplex>::find_pairs_M(int q, bool &found) const
{
    found = false;
    std::vector<Cell_pair> pairs;
    if (this->_hdvf_opt & (OPT_F | OPT_FULL))
    {
        // Search for +-1 in _F - iterate over rows
        for (OSM::Bitboard::iterator it_row = this->_F_row[q].begin(); it_row != this->_F_row[q].end() ; ++it_row)
        {
            const typename HDVF_parent::Row_chain &row(OSM::cget_row(this->_F_row[q],*it_row));

            // Iterate through the entries of the row
            for (typename HDVF_parent::Row_chain::const_iterator it = row.begin(); it != row.end(); ++it) {
                if ((it->second == 1) || (it->second == -1)) {
                    // If an entry with coefficient 1 or -1 is found, set the pair and add it
                    Cell_pair p ;
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

// find all the valid Cell_pair containing tau for M in dimension q
template<typename ChainComplex>
std::vector<Cell_pair> Hdvf<ChainComplex>::find_pairs_M(int q, bool &found, size_t tau) const
{
    found = false;
    std::vector<Cell_pair> pairs;
    if (this->_flag.at(q).at(tau) == SECONDARY)
        return pairs ; // Empty / found false

    if (this->_hdvf_opt & (OPT_F | OPT_FULL))
    {
        // If tau is primary, search for gamma such that <f(tau),gamma>=+-1
        if (this->_flag.at(q).at(tau) == PRIMARY)
        {
            std::cout << "------PRIMARY" << std::endl ;
            // Get the column of tau
            const typename HDVF_parent::Column_chain& col(OSM::get_column(this->_F_row.at(q), tau)) ;
            for (typename HDVF_parent::Column_chain::const_iterator it = col.cbegin(); it != col.cend(); ++it)
            {
                if (abs(it->second) == 1)
                {
                    Cell_pair p ;
                    p.sigma = tau ;
                    p.tau = it->first ;
                    p.dim = q ;
                    found = true ;
                    pairs.push_back(p);
                }
            }
        }
        // If tau is critical, search for pi such that <f(pi),tau>=+-1
        if (this->_flag.at(q).at(tau) == CRITICAL)
        {
            std::cout << "------CRITICAL" << std::endl ;
            const typename HDVF_parent::Row_chain& row(OSM::get_row(this->_F_row.at(q), tau)) ;
            for (typename HDVF_parent::Row_chain::const_iterator it = row.cbegin(); it != row.cend(); ++it)
            {
                if (abs(it->second) == 1)
                {
                    Cell_pair p ;
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

// find a valid Cell_pair for W in dimension q
template<typename ChainComplex>
Cell_pair Hdvf<ChainComplex>::find_pair_W(int q, bool &found) const
{
    found = false;
    Cell_pair p;

    if (this->_hdvf_opt & (OPT_G | OPT_FULL))
    {
        // Search for +-1 in _G - iterate over cols
        for (OSM::Bitboard::iterator it_col = this->_G_col[q].begin(); (it_col != this->_F_row[q].end() && !found); ++it_col)
        {
            typename HDVF_parent::Column_chain &col(this->_G_col[q][*it_col]);

            // Iterate through the entries of the col
            for (typename HDVF_parent::Column_chain::const_iterator it = col.begin(); (it != col.end() && !found); ++it) {
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

// find a valid Cell_pair containing tau for W in dimension q
template<typename ChainComplex>
Cell_pair Hdvf<ChainComplex>::find_pair_W(int q, bool &found, size_t tau) const
{
    found = false;
    Cell_pair p;
    if (this->_flag.at(q).at(tau) == PRIMARY)
        return p ; // Empty / found false

    if (this->_hdvf_opt & (OPT_G | OPT_FULL))
    {
        // If tau is primary, search for gamma such that <g(gamma),tau>=+-1
        if (this->_flag.at(q).at(tau) == SECONDARY)
        {
            for (OSM::Bitboard::iterator it_col = this->_G_col.at(q).begin(); (it_col != this->_G_col.at(q).end() && !found); ++it_col)
            {
                if (abs(this->_G_col.at(q).get_coefficient(tau, *it_col)) == 1)
                {
                    p.sigma = tau ; // secondary cell
                    p.tau = *it_col ; // critical cell
                    p.dim = q ;
                    found = true;
                }
            }
        }
        // If tau is critical, search for sigma such that <g(tau),sigma>=+-1
        if (this->_flag.at(q).at(tau) == CRITICAL)
        {
            typename HDVF_parent::Column_chain col(OSM::get_column(this->_G_col.at(q), tau)) ;
            for (typename HDVF_parent::Column_chain::const_iterator it = col.cbegin(); (it != col.cend() && !found); ++it)
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

// find all the valid Cell_pair for W in dimension q
template<typename ChainComplex>
std::vector<Cell_pair> Hdvf<ChainComplex>::find_pairs_W(int q, bool &found) const
{
    found = false;
    std::vector<Cell_pair> pairs;

    if (this->_hdvf_opt & (OPT_G | OPT_FULL))
    {
        // Search for +-1 in _G - iterate over cols
        for (OSM::Bitboard::iterator it_col = this->_G_col[q].begin(); it_col != this->_F_row[q].end(); ++it_col)
        {
            typename HDVF_parent::Column_chain &col = this->_G_col[q][*it_col];

            // Iterate through the entries of the col
            for (typename HDVF_parent::Column_chain::const_iterator it = col.begin(); it != col.end(); ++it) {
                if ((it->second == 1) || (it->second == -1)) {
                    // If an entry with coefficient 1 or -1 is found, set the pair and mark as found
                    Cell_pair p ;
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

// find all the valid Cell_pair containing tau for W in dimension q
template<typename ChainComplex>
std::vector<Cell_pair> Hdvf<ChainComplex>::find_pairs_W(int q, bool &found, size_t tau) const
{
    found = false;
    std::vector<Cell_pair> pairs;
    if (this->_flag.at(q).at(tau) == PRIMARY)
        return pairs ; // Empty / found false

    if (this->_hdvf_opt & (OPT_G | OPT_FULL))
    {
        // If tau is primary, search for gamma such that <g(gamma),tau>=+-1
        if (this->_flag.at(q).at(tau) == SECONDARY)
        {
            for (OSM::Bitboard::iterator it_col = this->_G_col.at(q).begin(); it_col != this->_G_col.at(q).end(); ++it_col)
            {
                if (abs(get_coefficient(this->_G_col.at(q), tau, *it_col)) == 1)
                {
                    Cell_pair p ;
                    p.sigma = tau ; // secondary cell
                    p.tau = *it_col ; // critical cell
                    p.dim = q ;
                    pairs.push_back(p) ;
                }
            }
        }
        // If tau is critical, search for sigma such that <g(tau),sigma>=+-1
        if (this->_flag.at(q).at(tau) == CRITICAL)
        {
            typename HDVF_parent::Column_chain col(OSM::get_column(this->_G_col.at(q), tau)) ;
            for (typename HDVF_parent::Column_chain::const_iterator it = col.cbegin(); it != col.cend(); ++it)
            {
                if (abs(it->second) == 1)
                {
                    Cell_pair p ;
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

// find a valid Cell_pair for MW in dimension q
template<typename ChainComplex>
Cell_pair Hdvf<ChainComplex>::find_pair_MW(int q, bool &found) const
{
    found = false;
    Cell_pair p;

    if (this->_hdvf_opt & OPT_FULL)
    {
        // pi and sigma must at least satisfy that col pi of H_q and row sigma of H_q-1 are non zero
        // iterate on H_q and H_q-1 accordingly
        for (OSM::Bitboard::iterator it_pi = this->_H_col[q].begin(); (it_pi != this->_H_col[q].end() && !found); ++it_pi)
        {
            const size_t pi = *it_pi ;
            typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
            typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
            typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;
            for (size_t sigma = 0; (sigma < this->_H_col.at(q-1).dimensions().first && !found); ++sigma)
            {
                typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
                if (!H11q1.is_null())
                {
                    typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
                    typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;

                    // Compute xi and xi' to test the validity of MW

                    const Coefficient_ring xi = projS_cod_sigma * H11 ;
                    const Coefficient_ring xip = H11q1 * projP_d_pi ;
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

// find a valid Cell_pair containing tau for MW in dimension q
template<typename ChainComplex>
Cell_pair Hdvf<ChainComplex>::find_pair_MW(int q, bool &found, size_t tau) const
{
    found = false;
    Cell_pair p;

    if (this->_flag.at(q).at(tau) == CRITICAL)
        return p ; // Empty / found false

    if (this->_hdvf_opt & OPT_FULL)
    {
        // If tau is primary (rename pi), search for a valid sigma
        if (this->_flag.at(q).at(tau) == PRIMARY)
        {
            const size_t pi(tau) ;
            // Col pi of H_q and proj_P(d(pi)) must at least be non empty
            typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
            typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
            typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;

            if (H11.is_null() || projP_d_pi.is_null())
                return p ;

            // Search for sigma with col_sigma(H_q-1) non empty
            for (size_t sigma = 0; (sigma < this->_H_col.at(q-1).dimensions().first && !found); ++sigma)
            {
                typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
                if (!H11q1.is_null())
                {
                    // and proj_S(cod(sigma)) non empty
                    typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
                    typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;
                    if (!projS_cod_sigma.is_null())
                    {
                        // test xi and xip
                        const Coefficient_ring xi = projS_cod_sigma * H11 ;
                        const Coefficient_ring xip = H11q1 * projP_d_pi ;
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
            const size_t sigma(tau) ;
            // Row sigma of H_q-1 and proj_S(cod(sigma)) must at least be non empty
            typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
            typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
            typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;
            if (H11q1.is_null() || projS_cod_sigma.is_null())
                return p ;

            // Search for pi with col pi of H_q non empty
            for (OSM::Bitboard::iterator it_pi = this->_H_col[q].begin(); (it_pi != this->_H_col[q].end() && !found); ++it_pi)
            {
                const size_t pi(*it_pi) ;

                typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
                typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
                typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;
                // proj_P of d(pi) must also be non empty
                if (!projP_d_pi.is_null())
                {
                    // test xi and xip
                    const Coefficient_ring xi = projS_cod_sigma * H11 ;
                    const Coefficient_ring xip = H11q1 * projP_d_pi ;
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

// find all the valid Cell_pair for MW in dimension q
template<typename ChainComplex>
std::vector<Cell_pair> Hdvf<ChainComplex>::find_pairs_MW(int q, bool &found) const
{
    found = false;
    std::vector<Cell_pair> pairs;

    if (this->_hdvf_opt & OPT_FULL)
    {
        // pi and sigma must at least satisfy that col pi of H_q and row sigma of H_q-1 are non zero
        // iterate on H_q and H_q-1 accordingly
        for (OSM::Bitboard::iterator it_pi = this->_H_col[q].begin(); it_pi != this->_H_col[q].end(); ++it_pi)
        {
            const size_t pi = *it_pi ;
            typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
            typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
            typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;
            for (size_t sigma = 0; sigma < this->_H_col.at(q-1).dimensions().first; ++sigma)
            {
                typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
                if (!H11q1.is_null())
                {
                    typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
                    typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;

                    // Compute xi and xi' to test the validity of MW

                    const Coefficient_ring xi = projS_cod_sigma * H11 ;
                    const Coefficient_ring xip = H11q1 * projP_d_pi ;
                    found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                    if (found)
                    {
                        Cell_pair p;
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

// find all the valid Cell_pair containing tau for MW in dimension q
template<typename ChainComplex>
std::vector<Cell_pair> Hdvf<ChainComplex>::find_pairs_MW(int q, bool &found, size_t tau) const
{
    found = false;
    std::vector<Cell_pair> pairs;
    if (this->_flag.at(q).at(tau) == CRITICAL)
        return pairs ; // Empty / found false

    if (this->_hdvf_opt & OPT_FULL)
    {
        // If tau is primary (rename pi), search for a valid sigma
        if (this->_flag.at(q).at(tau) == PRIMARY)
        {
            const size_t pi(tau) ;
            // Col pi of H_q and proj_P(d(pi)) must at least be non empty
            typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
            typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
            typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;

            if (H11.is_null() || projP_d_pi.is_null())
                return pairs ;

            // Search for sigma with col_sigma(H_q-1) non empty
            for (size_t sigma = 0; sigma < this->_H_col.at(q-1).dimensions().first; ++sigma)
            {
                typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
                if (!H11q1.is_null())
                {
                    // and proj_S(cod(sigma)) non empty
                    typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
                    typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;
                    if (!projS_cod_sigma.is_null())
                    {
                        // test xi and xip
                        const Coefficient_ring xi = projS_cod_sigma * H11 ;
                        const Coefficient_ring xip = H11q1 * projP_d_pi ;
                        found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                        if (found)
                        {
                            Cell_pair p ;
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
            const size_t sigma(tau) ;
            // Row sigma of H_q-1 and proj_S(cod(sigma)) must at least be non empty
            typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
            typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
            typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;
            if (H11q1.is_null() || projS_cod_sigma.is_null())
                return pairs ;

            // Search for pi with col pi of H_q non empty
            for (OSM::Bitboard::iterator it_pi = this->_H_col[q].begin(); it_pi != this->_H_col[q].end(); ++it_pi)
            {
                const size_t pi(*it_pi) ;

                typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
                typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
                typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;
                // proj_P of d(pi) must also be non empty
                if (!projP_d_pi.is_null())
                {
                    // test xi and xip
                    const Coefficient_ring xi = projS_cod_sigma * H11 ;
                    const Coefficient_ring xip = H11q1 * projP_d_pi ;
                    found = ((abs(xi) == 1) && (abs(xip) == 1)) ;
                    if (found)
                    {
                        Cell_pair p ;
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


// Method to perform operation R
// pi is in dimension q, sigma is in dimension q+1
template<typename ChainComplex>
void Hdvf<ChainComplex>::R(size_t pi, size_t sigma, int q) {
    //----------------------------------------------- Submatrices of H ----------------------------------------------------

    if (this->_hdvf_opt & OPT_FULL)
    {
        // Output operation details to the console
        std::cout << "R of " << pi << "(dim " << q << ") / " << sigma << "(dim " << q + 1 << ")" << std::endl;

        // Extract the relevant row and column chains from this->_H_col
        typename HDVF_parent::Row_chain H12 = OSM::get_row(this->_H_col[q], sigma); // H12 is the row chain from this->_H_col[q] at index sigma
        typename HDVF_parent::Column_chain H21 = OSM::get_column(this->_H_col[q], pi); // H21 is the column chain from this->_H_col[q] at index pi

        // Get the coefficient at the intersection of H12 and H21
        Coefficient_ring H11(H12[pi]); // H11 is the coefficient at row sigma and column pi

        // Assert that H11 is either 1 or -1 (check invertibility)
        assert((H11 == 1) || (H11 == -1)); // !!!!! Test invertibility
        Coefficient_ring H11_inv = H11; // Inverse of H11 (which is itself for 1 or -1)

        // Remove the contributions of pi from H12 and sigma from H21
        H12 /= std::vector<size_t>({pi}); // Remove column pi from H12
        H21 /= std::vector<size_t>({sigma}); // Remove row sigma from H21

        // Remove the corresponding row and column from this->_H_col
        remove_row(this->_H_col[q], sigma); // Remove row sigma from this->_H_col[q]
        remove_column(this->_H_col[q], pi); // Remove column pi from this->_H_col[q]

        //---------------------------------------------- Submatrices of F -----------------------------------------------------

        // Extract the relevant column chain from this->_F_row
        typename HDVF_parent::Column_chain F11 = OSM::get_column(this->_F_row[q], pi); // F11 is the column chain from this->_F_row[q] at index pi

        // Remove the column pi from this->_F_row
        remove_column(this->_F_row[q], pi);

        //--------------------------------------------- Submatrices of G ------------------------------------------------------

        // Extract the relevant row chain from this->_G_col
        typename HDVF_parent::Row_chain G11 = OSM::get_row(this->_G_col[q + 1], sigma); // G11 is the row chain from this->_G_col[q+1] at index sigma

        // Remove the row sigma from this->_G_col
        remove_row(this->_G_col[q + 1], sigma);

        //--------------------------------------------- Update matrices -------------------------------------------------------

        // Update this->_F_row[q]
        // Subtract the product of (H12 * F11) and H11 from this->_F_row[q]
        // Note: % operator returns a row matrix, so be careful with operations
        this->_F_row[q] -= (F11 % H12) * H11_inv;

        // Set the row pi of this->_F_row[q] to (H12 * (-H11_inv))
        OSM::set_row(this->_F_row[q], pi, H12 * (-H11_inv));

        // Update this->_G_col[q + 1]
        // Subtract the product of (H21 * G11) and H11_inv from this->_G_col[q + 1]
        this->_G_col[q + 1] -= (H21 * G11) * H11_inv;

        // Set the column sigma of this->_G_col[q + 1] to (H21 * (-H11_inv))
        OSM::set_column(this->_G_col[q + 1], sigma, H21 * (-H11_inv));

        // Update this->_DD_col
        // Compute the temporary matrix product F11 * G11
        this->_DD_col[q + 1] += (F11 * G11) * H11_inv;

        // Set the row pi in this->_DD_col[q + 1] to (G11 * H11_inv)
        OSM::set_row(this->_DD_col[q + 1], pi, G11 * H11_inv);

        // Set the column sigma in this->_DD_col[q + 1] to (F11 * H11_inv)
        OSM::set_column(this->_DD_col[q + 1], sigma, F11 * H11_inv);

        // Set the coefficient at (pi, sigma) in this->_DD_col to H11_inv
        this->_DD_col[q + 1].set_coefficient(pi, sigma, H11_inv);

        // Update this->_H_col[q]
        // Subtract the product of (H21 * H12) and H11_inv from this->_H_col[q]
        this->_H_col[q] -= (H21 * H12) * H11_inv;

        // Perform additional updates

        // Extract boundary and coboundary chains
        typename HDVF_parent::Column_chain bnd_pi(this->_K.d(pi, q)); // Boundary of pi in dimension q
        typename HDVF_parent::Row_chain cobnd_sigma(this->_K.cod(sigma, q + 1)); // Coboundary of sigma in dimension q+1

        // Project the boundary and coboundary chains onto PRIMARY and SECONDARY flags
        typename HDVF_parent::Column_chain proj_P_pi(this->projection(bnd_pi, PRIMARY, q));
        typename HDVF_parent::Row_chain proj_S_sigma(this->projection(cobnd_sigma, SECONDARY, q + 1));

        if (q > 0) {
            // Update this->_DD_col[q] with projections and this->_F_row[q-1]
            typename HDVF_parent::Column_chain c1(this->_F_row[q - 1] * proj_P_pi + this->projection(bnd_pi, CRITICAL, q));
            OSM::set_column(this->_DD_col[q], pi, c1);
            OSM::set_column(this->_G_col[q], pi, this->_H_col[q - 1] * proj_P_pi);
        }

        // Update this->_F_row[q+1] with this->projection and this->_H_col[q+1]
        OSM::set_row(this->_F_row[q + 1], sigma, proj_S_sigma * this->_H_col[q + 1]);

        if (q + 2 <= this->_K.dimension()) {
            // Update this->_DD_col[q+2] with projections
            typename HDVF_parent::Row_chain c4(this->projection(cobnd_sigma, CRITICAL, q + 1) + proj_S_sigma * this->_G_col[q + 2]);
            OSM::set_row(this->_DD_col[q + 2], sigma, c4);
        }

        // Update flags
        this->_flag[q][pi] = CRITICAL; // Set the PSC_flag of pi in dimension q to CRITICAL
        ++this->_nb_C.at(q) ;
        --this->_nb_P.at(q) ;
        this->_flag[q + 1][sigma] = CRITICAL; // Set the PSC_flag of sigma in dimension q+1 to CRITICAL
        ++this->_nb_C.at(q+1) ;
        --this->_nb_S.at(q+1) ;
    }
    else
        std::cout << "!!! R impossible with partial reduction options" << std::endl;
}

// Method to perform operation M
// pi is in dimension q, gamma is in dimension q
template<typename ChainComplex>
void Hdvf<ChainComplex>::M(size_t pi, size_t gamma, int q) {
    //----------------------------------------------- Submatrices of F ----------------------------------------------------

    if (this->_hdvf_opt & OPT_FULL)
    {
        // Output the operation details to the console
        std::cout << "M_" << q << "(" << pi << "," << gamma << ")" << std::endl;

        if (q == this->_K.dimension())
            throw("M operation in max dimension !!!") ;

        // Extract row and column chains from this->_F_row
        typename HDVF_parent::Row_chain F12(OSM::get_row(this->_F_row[q], gamma)); // F12 is the row chain from this->_F_row[q] at index gamma
        typename HDVF_parent::Column_chain F21(OSM::get_column(this->_F_row[q], pi)); // F21 is the column chain from this->_F_row[q] at index pi

        // Get the coefficient at the intersection of F12 and F21
        const Coefficient_ring F11(F12.get_coefficient(pi)); // F11 is the coefficient at row gamma and column pi

        // Assert that F11 is either 1 or -1 (for invertibility)
        assert((F11 == 1) || (F11 == -1)); // !!!!! Test invertibility
        Coefficient_ring F11_inv = F11; // Inverse of F11 (which is itself for 1 or -1)

        // Remove the contributions of pi from F12 and gamma from F21
        F12 /= std::vector<size_t>({pi}); // Remove column pi from F12
        F21 /= std::vector<size_t>({gamma}); // Remove row gamma from F21

        // Remove the corresponding row and column from this->_F_row
        remove_row(this->_F_row[q], gamma); // Remove row gamma from this->_F_row[q]
        remove_column(this->_F_row[q], pi); // Remove column pi from this->_F_row[q]

        //--------------------------------------------- Submatrices of G ------------------------------------------------------

        // Extract the relevant column chain from this->_G_col
        remove_column(this->_G_col[q], gamma); // Remove column gamma from this->_G_col[q]

        //---------------------------------------------- Submatrices of H -----------------------------------------------------

        // Extract the relevant column chain from this->_H_col
        typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col[q], pi)); // H11 is the column chain from this->_H_col[q] at index pi

        // Remove the column pi from this->_H_col
        remove_column(this->_H_col[q], pi);

        //--------------------------------------------- Submatrices of DD_q+1 ------------------------------------------------------

        // For DD_q+1 and DD_q:
        // Extract the relevant row chains from this->_DD_col
        typename HDVF_parent::Row_chain D11(OSM::get_row(this->_DD_col[q+1], gamma)); // D11 is the row chain from this->_DD_col[q+1] at index gamma
        remove_row(this->_DD_col[q + 1], gamma); // Remove row gamma from this->_DD_col[q + 1]

        //--------------------------------------------- Submatrices of DD ------------------------------------------------------

        // DD_q (corresponds to the column matrix of this->_DD_col)
        if (q > 0)
            remove_column(this->_DD_col[q], gamma); // Remove column gamma from this->_DD_col[q]

        //--------------------------------------------- Update matrices -------------------------------------------------------

        // Update this->_H_col
        // Subtract the product of (H11 * F11_inv) and F12 from this->_H_col[q]
        // Note: * operator is used for matrix column results
        this->_H_col[q] -= (H11 * F11_inv) * F12;

        // Set the column gamma of this->_H_col[q] to (-1 * H11) * F11_inv
        OSM::set_column(this->_H_col[q], gamma, (-1 * H11) * (F11_inv));

        // Update this->_F_row
        // Add the product of (F21 * F11_inv) and F12 to this->_F_row[q]
        this->_F_row[q] -= (F21 * F11_inv) * F12;

        // Set the row pi of this->_F_row[q] to F11_inv * F12
        OSM::set_row(this->_F_row[q], pi, F11_inv * F12);

        // Set the column gamma of this->_F_row[q] to F21 * (-F11_inv)
        OSM::set_column(this->_F_row[q], gamma, F21 * (-F11_inv));

        // Set the coefficient at (pi, gamma) in this->_F_row[q] to F11_inv
        set_coefficient(this->_F_row[q], pi, gamma, F11_inv);

        // Update this->_G_col
        // Add the product of (H11 * F11_inv) and D11 to this->_G_col[q+1]
        this->_G_col[q + 1] += (H11 * F11_inv) * D11;

        // Update this->_DD_col[q + 1]
        this->_DD_col[q + 1] -= (F21 * F11_inv) * D11; // Subtract the product from this->_DD_col[q + 1]
        OSM::set_row(this->_DD_col[q + 1], pi, D11 * F11_inv); // Set the row pi in this->_DD_col[q + 1]

        // Update this->_G_col (for dimension q) and this->_DD_col (for dimension q)
        if (q>0)
        {
            // Extract boundary chain and project it
            typename HDVF_parent::Column_chain c = this->_K.d(pi, q); // Boundary of pi in dimension q
            typename HDVF_parent::Column_chain projection_p(this->projection(c, PRIMARY, q-1)); // Project boundary chain to PRIMARY
            typename HDVF_parent::Column_chain projection_c = this->projection(c, CRITICAL, q-1); // Project boundary chain to CRITICAL

            // Set the column pi of this->_G_col[q] to (-1 * this->_H_col[q-1]) * projection_p
            OSM::set_column(this->_G_col[q], pi, (Coefficient_ring(-1) * this->_H_col[q - 1]) * projection_p);

            // Update this->_DD_col
            // Extract projections and perform updates
            typename HDVF_parent::Column_chain tmp(this->_F_row[q - 1] * projection_p + projection_c); // Compute the product of this->_F_row[q - 1] and projection_p
            // Set the column pi of this->_DD_col[q] to cc
            OSM::set_column(this->_DD_col[q], pi, tmp);
        }
        // else: if the dimension is 0, no change for this->_G_col[q] and this->_DD_col[q]

        // Update flags
        this->_flag[q][pi] = CRITICAL; // Set the PSC_flag of pi in dimension q to PRIMARY
        this->_flag[q][gamma] = PRIMARY; // Set the PSC_flag of gamma in dimension q to CRITICAL
    }
    else
        std::cout << "!!! M impossible with partial reduction options" << std::endl;
}

// Method to perform operation W
// gamma is in dimension q, sigma is in dimension q
template<typename ChainComplex>
void Hdvf<ChainComplex>::W(size_t sigma, size_t gamma, int q) {
    //----------------------------------------------- Submatrices of G ----------------------------------------------------

    if (this->_hdvf_opt & OPT_FULL)
    {
        // Output the operation details to the console
        std::cout << "W_" << q << "(" << sigma << "," << gamma << ")" << std::endl;

        if (q == 0)
            throw("W operation in dimension 0 !!!") ;

        // Extract row and column chains from this->_G_col
        typename HDVF_parent::Row_chain G12(OSM::get_row(this->_G_col[q], sigma)); // G12 is the row chain from this->_G_col[q] at index sigma
        typename HDVF_parent::Column_chain G21(OSM::get_column(this->_G_col[q], gamma)); // G21 is the column chain from this->_G_col[q] at index gamma

        // Get the coefficient at the intersection of G12 and G21
        Coefficient_ring G11(G12.get_coefficient(gamma)); // G11 is the coefficient at row sigma and column gamma

        // Assert that G11 is either 1 or -1 (for invertibility)
        assert((G11 == 1) || (G11 == -1)); // !!!!! Test invertibility
        Coefficient_ring G11_inv = G11; // Inverse of G11 (which is itself for 1 or -1)

        // Remove the contributions of gamma from G12 and sigma from G21
        G12 /= std::vector<size_t>({gamma}); // Remove column gamma from G12
        G21 /= std::vector<size_t>({sigma}); // Remove row sigma from G21

        // Remove the corresponding row and column from this->_G_col
        remove_row(this->_G_col[q], sigma); // Remove row sigma from this->_G_col[q]
        remove_column(this->_G_col[q], gamma); // Remove column gamma from this->_G_col[q]

        //---------------------------------------------- Submatrices of F -----------------------------------------------------

        // Extract the row chain from this->_F_row

        // Remove the row gamma from this->_F_row
        remove_row(this->_F_row[q], gamma);

        //--------------------------------------------- Submatrices of H ------------------------------------------------------

        // Extract the row chain from this->_H_col
        typename HDVF_parent::Row_chain H11(OSM::get_row(this->_H_col[q-1], sigma)); // H11 is the row chain from this->_H_col[q] at index sigma

        // Remove the row sigma from this->_H_col
        remove_row(this->_H_col[q-1], sigma);

        //--------------------------------------------- Submatrices of DD_q+1 ------------------------------------------------------

        // Extract the row chain from this->_DD_col[q+1]

        // Remove the row gamma from this->_DD_col
        if (q < this->_K.dimension())
            remove_row(this->_DD_col[q + 1], gamma);

        //--------------------------------------------- Submatrices of DD_q ------------------------------------------------------

        // Extract the column chain from this->_DD_col[q]
        typename HDVF_parent::Column_chain D11_q(OSM::get_column(this->_DD_col[q], gamma)); // D11_q is the column chain from this->_DD_col[q] at index gamma

        // Remove the column gamma from this->_DD_col
        remove_column(this->_DD_col[q], gamma);

        //--------------------------------------------- Update matrices -------------------------------------------------------

        // Update this->_H_col
        // Subtract the product of (G21 * G11_inv) and H11 from this->_H_col[q]
        // Note: % operator is used for matrix row results
        this->_H_col[q-1] -= (G21 * G11_inv) * H11;

        // Set the row gamma of this->_H_col[q] to (-1 * H11) * G11_inv
        OSM::set_row(this->_H_col[q-1], gamma, (-1 * H11) * (G11_inv));

        // Update this->_G_col
        // Subtract the product of G11_inv and (G21 * G12) from this->_G_col[q]
        this->_G_col[q] -= (G21 * G11_inv) * G12;

        // Set the row gamma of this->_G_col[q] to G11_inv * G12
        OSM::set_row(this->_G_col[q], gamma, -G11_inv * G12);

        // Set the column sigma of this->_G_col[q] to G11_inv * G21
        OSM::set_column(this->_G_col[q], sigma, G21 * G11_inv);

        // Set the coefficient at (gamma, sigma) in this->_G_col[q] to G11_inv
        set_coefficient(this->_G_col[q], gamma, sigma, G11_inv);

        // Update this->_F_row (for dimension q-1)
        // Add the product of (D11_q * G11_inv) and H11 to this->_F_row[q - 1]
        this->_F_row[q - 1] += (D11_q * G11_inv) * H11;

        // Update this->_DD_col (for dimension q)
        this->_DD_col[q] -= (D11_q * G11_inv) * G12;// Subtract the product of (D11_q * G11_inv) and G12
        OSM::set_column(this->_DD_col[q], sigma, D11_q * G11_inv) ;

        // Update this->_F_row (for dimension q) and this->_DD_col (for dimension q+1)
        if (q < this->_K.dimension())
        {
            // Extract boundary chain and project it
            typename HDVF_parent::Row_chain c = this->_K.cod(sigma, q); // Boundary of sigma in dimension q
            typename HDVF_parent::Row_chain projection_s(this->projection(c, SECONDARY, q+1)); // Project boundary chain to SECONDARY
            typename HDVF_parent::Row_chain projection_c(this->projection(c, CRITICAL, q+1)); // Project boundary chain to SECONDARY

            // Set the row sigma of this->_F_row[q] to (-1 * projection_s * this->_H_col[q])
            OSM::set_row(this->_F_row[q], sigma, (-1 * projection_s) * this->_H_col[q]);

            // Set the row sigma of this->_DD_col[q + 1] to projection_s * this->_G_col[q + 1] + projection_c

            typename HDVF_parent::Row_chain tmp(projection_s * this->_G_col[q + 1] + projection_c) ;
            OSM::set_row(this->_DD_col[q + 1], sigma, tmp);
        }
        // else : if the dimension is maximal, no update of this->_DD_col[q+1] and this->_F_row[q]

        // Update flags
        this->_flag[q][gamma] = SECONDARY; // Set the PSC_flag of gamma in dimension q to SECONDARY
        this->_flag[q][sigma] = CRITICAL; // Set the PSC_flag of sigma in dimension q to CRITICAL
    }
    else
        std::cout << "!!! W impossible with partial reduction options" << std::endl;
}

// Method to perform operation MW
// gamma is in dimension q, sigma is in dimension q
template<typename ChainComplex>
void Hdvf<ChainComplex>::MW(size_t pi, size_t sigma, int q) {
    //----------------------------------------------- Submatrices of G ----------------------------------------------------

    if (this->_hdvf_opt & OPT_FULL)
    {
        // Output the operation details to the console
        std::cout << "MW_" << q << "(" << pi << "," << sigma << ")" << std::endl;

        if (q <= 0)
            throw("MW operation in dimension 0 !!!") ;
        if (q >= this->_K.dimension())
            throw("MW operation in maximal dimension !!!") ;

        // In order to compute xi and xi', extract sub-matrices of H_q, H_q-1 and compute d(pi) and cod(sigma)

        // H_q extractions

        typename HDVF_parent::Column_chain H11(OSM::get_column(this->_H_col.at(q), pi)) ;
        // H21 -> delete H11
        this->_H_col.at(q) /= std::vector<size_t>({pi}) ;

        // H_q-1 extractions

        typename HDVF_parent::Row_chain H11q1(OSM::get_row(this->_H_col.at(q-1), sigma)) ;
        // H21_q-1 -> delete H11q1
        remove_row(this->_H_col.at(q-1), sigma) ;

        // d(pi)

        typename HDVF_parent::Column_chain d_pi = this->_K.d(pi, q) ;
        typename HDVF_parent::Column_chain projP_d_pi = this->projection(d_pi, PRIMARY, q-1) ;
        typename HDVF_parent::Column_chain projC_d_pi = this->projection(d_pi, CRITICAL, q-1) ;

        // cod(sigma)

        typename HDVF_parent::Row_chain cod_sigma = this->_K.cod(sigma, q) ;
        typename HDVF_parent::Row_chain projS_cod_sigma = this->projection(cod_sigma, SECONDARY, q+1) ;
        typename HDVF_parent::Row_chain projC_cod_sigma = this->projection(cod_sigma, CRITICAL, q+1) ;

        // Compute xi and xi' to test the validity of MW

        Coefficient_ring xi = projS_cod_sigma * H11 ;
        Coefficient_ring xip = H11q1 * projP_d_pi ;

        if (abs(xi) != 1)
            throw "MW impossible, xi non invertible" ;
        if (abs(xip) != 1)
            throw "MW impossible, xi' non invertible" ;

        // F_q extraction

        typename HDVF_parent::Column_chain F11(OSM::get_column(this->_F_row.at(q), pi)) ;
        // F12 -> delete col F11
        remove_column(this->_F_row.at(q), pi) ;

        // G_q extractions

        typename HDVF_parent::Row_chain G11(OSM::get_row(this->_G_col.at(q), sigma)) ;
        // G21 -> dele row G11
        remove_row(this->_G_col.at(q), sigma) ;

        // ----------- Update of the reduction

        // H_q

        typename HDVF_parent::Row_chain tmp1 = projS_cod_sigma * this->_H_col.at(q) ;

        this->_H_col.at(q) -= (H11 * xi) * tmp1 ;
        OSM::set_column(this->_H_col.at(q), sigma, H11 * xi) ;

        // F_q

        this->_F_row.at(q) += (F11 * xi) * tmp1 ;
        OSM::set_column(this->_F_row.at(q), sigma, F11 * (-xi)) ;

        // G_q+1 // note: G_q+1 is not be modified if the Hdvf is perfect

        typename HDVF_parent::Row_chain tmp2(projS_cod_sigma * this->_G_col.at(q+1)) ;
        tmp2 += projC_cod_sigma ;
        this->_G_col.at(q+1) -= (H11 * xi) * tmp2 ;

        // DD_col_q+1 / DD_row_q

        this->_DD_col.at(q+1) -= (F11 * xi) * tmp2 ;

        // H_q-1

        typename HDVF_parent::Column_chain tmp3(this->_H_col.at(q-1) * projP_d_pi) ;

        this->_H_col.at(q-1) -= (tmp3 * xip) * H11q1 ;
        OSM::set_row(this->_H_col.at(q-1), pi, xip * H11q1) ;

        // G_q

        this->_G_col.at(q) -= (tmp3 * xip) * G11 ;
        OSM::set_row(this->_G_col.at(q), pi, xip * G11) ;

        // F_q-1 // note: F_q-1 is not be modified if the Hdvf is perfect

        typename HDVF_parent::Column_chain tmp4(this->_F_row.at(q-1) * projP_d_pi) ;
        tmp4 += projC_d_pi ;
        this->_F_row.at(q-1) -= (tmp4 * xip) * H11q1 ;

        // DD_col_q

        this->_DD_col.at(q) += (tmp4 * xip) * G11 ;

        // Update flags
        this->_flag[q][pi] = SECONDARY; // Set the PSC_flag of gamma in dimension q to SECONDARY
        this->_flag[q][sigma] = PRIMARY; // Set the PSC_flag of sigma in dimension q to PRIMARY
    }
    else
        std::cout << "!!! MW impossible with partial reduction options" << std::endl;
}

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_HDVF_H
