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

#ifndef CGAL_HDVF_PERSISTENCE_H
#define CGAL_HDVF_PERSISTENCE_H

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <functional>
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/HDVF/SubSparseMatrix.hpp"
#include "CGAL/HDVF/Hdvf_core.h"
#include "CGAL/HDVF/Filtration_lower_star.h"

namespace CGAL {
namespace HDVF {


/**
 * \class Persistent HDVF
 * \brief Implementation of persistent HDVF (inherits HDVF).
 *
 * The Hdvf_persistence class contains all functions to build persistent HDVF and create perfect persistent HDVF.
 *
 * \tparam CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam ComplexType The type of complex  (default is SimpComplex)
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 */

// Types for persistent homology / cohomology

/*! \brief Type to store persistent intervals filtration indices (birth/death indices)
 *
 * For "infinite" intervals borne at index i, the interval is set to (i,i-1)
 */
typedef std::pair<int, int> FiltrIndexPerInterval ;

/*! \brief Type for indexing uniquely a cell.
 * - First element of the pair: index of the cell.
 * - Second element of the pair: dimension of the cell.
 *
 * "Infinite cells" are defined as (-1,q+1) where q is the dimension of the chain complex.
 */
typedef std::pair<int, int> CellDim ;

/*! \brief Type for describing the pair of cells associated to a persistence interval:
 * - First element of the pair: cell entailing the birth of the hole.
 * - Second element of the pair: cell entailing the death of the hole.
 *
 * For infinite intervals, the "infinite cells" is defined as (-1,q+1) where q is the dimension of the chain complex.
 */
typedef std::pair<CellDim, CellDim> CellsPerInterval ;

/*! \brief Template for persistent intervals degrees (birth/death degrees)
 *
 * For "infinite" intervals borne at degree d, the interval is set to (d,d-1)
 */
template <typename DegType>
using DegreePerIntervalT =  std::pair<DegType,DegType> ;

/*! \brief Template for (full) persistent interval data:
 * - First element: persistent interval filtration indices
 * - Second element: persistent interval cells
 * - Third element: persistent interval degrees
 */
template <typename DegType>
using PerHoleT = std::tuple<FiltrIndexPerInterval, CellsPerInterval, DegreePerIntervalT<DegType> > ;

/*! \brief Overload of the `<<` operator to display persistent intervals (that is PerHoleT).
 *
 * Format:
 */
template <typename DegType>
ostream& operator<< (ostream& out, const PerHoleT<DegType>& hole)
{
    // time (cell, dim) -> time (cell, dim) / time to live
    const FiltrIndexPerInterval per_int(std::get<0>(hole)) ;
    const CellsPerInterval per_int_cells(std::get<1>(hole)) ;
    const DegreePerIntervalT<DegType> per_int_deg(std::get<2>(hole)) ;
    
    const DegType deg_duration(per_int_deg.second-per_int_deg.first) ;
    if (deg_duration >= 0) // finite interval
    {
        out << "[" << per_int.first << " (" << per_int_cells.first.first << ", " << per_int_cells.first.second << ") -> " ;
        out << per_int.second << " (" << per_int_cells.second.first << ", " << per_int_cells.second.second << ") / duration: " << deg_duration << "]" << endl ;
    }
    else
    {
        out << "[" << per_int.first << " (" << per_int_cells.first.first << ", " << per_int_cells.first.second << ") -> " ;
        out << "inf]" << endl ;
    }
    return out ;
}

/*!
 \ingroup PkgHDVFAlgorithmClasses
 
 The class `Hdvf_persistence` computes persistent homology using HDVFs. Hence, unlike other persistence algorithms, beside standard persistent intervals informations (birth/death indices, degrees, associated cells), `Hdvf_persistence` also provides homology and cohomology generators for persistent pairs. Intuitively, holes die when they are "filled" by a cell: associated homology and cohomology generators provide a representation of the hole and of the cells filling the hole.
 
 Given a `Filtration`, the `Hdvf_persistence`constructor basically builds a HDVF where indices of cells in the bases follow the filtration order (permutations between initial indices of cells in the chain complex and new indices given by the filtration are stored).
 Besides, `Hdvf_persistence` derives from `Hdvf_core` with the SparseMatrix parameter type set to `Sub_sparse_matrix`. Hence, `Hdvf_persistence` computes homology over a larger and larger sub-complex of `K` (encoded through a Bitboard mask of `Sub_sparse_matrix`) following the filtration `f`. At each step of the filtration, the new cell is paired (A operation) with the youngest cell valid for A.
 
 Persistent intervals can be exported to a file; moreover, persistent homology and cohomology generators can be exported to chains (and to vtk).
 
 \cgalModels{HDVF}
 
 \tparam CoefficientType a model of the `Ring` concept (by default, we use the `Z` model) providing the ring used to compute homology.
 \tparam ComplexType a model of the `AbstractChainComplex` concept, providing the type of abstract chain complex used.
 \tparam DegType a scalar data type used for the degrees of the filtration.
 \tparam FiltrationType a model of the `Filtration` concept, providing the filtration used to compute persistence.
 */

template<typename CoefficientType, typename ComplexType, typename DegType, typename FiltrationType >
class Hdvf_persistence : public Hdvf_core<CoefficientType,ComplexType, OSM::Chain, OSM::SubSparseMatrix>
{
public:
    // Matrices types
    /*!
     Type of column-major chains
     */
    typedef OSM::Chain<CoefficientType, OSM::COLUMN> CChain;
    
    /*!
     Type of row-major chains
     */
    typedef OSM::Chain<CoefficientType, OSM::ROW> RChain;
    
    /*!
     Type of column-major sparse matrices
     */
    typedef OSM::SubSparseMatrix<CoefficientType, OSM::COLUMN> CMatrix;
    
    /*!
     Type of row-major sparse matrices
     */
    typedef OSM::SubSparseMatrix<CoefficientType, OSM::ROW> RMatrix;
    
    /*! Type of parent HDVF class (Hdvf_core with appropriate template parameters)
     * The SparseMatrix model is set to Sub_sparse_matrix to activate (co)homology computation over a subcomplex.
     */
    typedef Hdvf_core<CoefficientType, ComplexType, OSM::Chain, OSM::SubSparseMatrix> HDVFParent ;
    
    /*! Type of filtrations used to compute persistence.
     */
    typedef FiltrationType Filtration;
    
    /*! Instanciation of `DegreePerIntervalT` (persistant intervals degrees) for `DegType`.
     */
    typedef DegreePerIntervalT<DegType> DegreePerInterval;
    
    /*! Instanciation of `PerHoleT` (full persistant intervals informations) for `DegType`.
     */
    typedef PerHoleT<DegType> PerHole;
    
    /*! Type for PRIMARY/SECONDARY/CRITICAL labels export
     *      Encoding used:
     *          - PRIMARY: -1
     *          - SECONDARY: +1
     *          - CRITICAL: 0
     */
    typedef std::vector<std::vector<int> > ExpLabels ;
    
    /*! Type  chains export (column-major chains)
     */
    typedef CChain ExpChain ;
    
    /*! For persistent diagram iterator (returned by *it)
     */
    typedef struct {
        PerHole hole ;
        ExpLabels labelsPSC ;
        ExpChain g_chain_sigma, g_chain_tau, fstar_chain_sigma, fstar_chain_tau ;
    } PerIntervalInformation ;
    
protected:
    /** \brief Reference to the filtration used for persistence */
    const Filtration &_f ;
    
    /** \brief Permutation between indices in the chain complex K and indices along the filtration
     * Indices along the filtration provide new indices for cells in each dimension.
     */
    std::vector<std::vector<int> > _K_to_per, _per_to_K ;
    
    /** \brief Vector of persistent pairs computed */
    std::vector<PerHole> _persist ;
    
    /** \brief Boolean determining weather or not export homology/cohomology generators associated to persistent pairs
     * - If _with_export is true, PSC labels and homology/cohomology generators are stored for each persistent pair of duration (that is, such as the difference between degrees of birth/death) strictly positive.
     * - If _with_export is false, only persistent intervals are stored.
     */
    bool _with_export ;
    
    /** \brief Vector of exported PSC labels */
    std::vector<ExpLabels> _export_labels ;
    
    /** \brief Vector of exported homology/cohomology generators */
    std::vector<std::pair<ExpChain, ExpChain> > _export_g, _export_fstar ;
    
    /** Current time (that is, current index) in the filtration */
    int _t ;
    /** Current time (that is, current index) along each dimension */
    std::vector<int> _t_dim ;
    
    /** Bitboard masks in each dimension
     * At a given filtration time, only cells already met have a bit set to 1
     */
    std::vector<OSM::Bitboard> _masks ;
    
    
private:
    // Hide find_pair_A methods of the Hdvf_core class
    // Hide A operation of the Hdvf_core class
    // This operation is redefined with a different prototype in Hdvf_persistence and set as private since persistence lets no choice for A pairing (rule of the "youngest")
    using HDVFParent::find_pair_A;
    using HDVFParent::find_pairs_A;
    using HDVFParent::A;
public:
    /**
     * \brief Hdvf_persistence default constructor
     *
     * Builds an "empty" HDVF_persistence (with all cells critical) associated to the chain complex `K` and the filtration `f`.
     * By default, the HDVF option is set to OPT_FULL (full reduction computed)
     *
     * \param[in] K A chain complex (a model of `AbstractChainComplex`).
     * \param[in] f A filtration (a model of `Filtration`).
     * \param[in] hdvf_opt Option for HDVF computation (`OPT_BND`, `OPT_F`, `OPT_G` or `OPT_FULL`)
     * \param[in] with_export Boolean option to activate or not the export of PSC labels and homology/cohomology generators for of persistent intervals of positive duration. This information is used by vtk exporters.
     */
    Hdvf_persistence(const ComplexType& K, const Filtration& f, int hdvf_opt = OPT_BND, bool with_export = false) ;
    
    /**
     * \brief Compute a perfect persistent HDVF.
     *
     * This method follows the filtration and considers cells one by one. For each of them, it searches the youngest possible cell valid for A (returned by `find_pair_A`), and applies the corresponding A operation.
     * By definition of persistent homology, the `Ring` of coefficients *must be* a field.
     *
     * \param[in] verbose If this parameter is `true`, all intermediate reductions are printed out.
     *
     * \returns The vector of all `PairCell` paired with A.
     */
    std::vector<PairCell> compute_perfect_hdvf(bool verbose = false)
    {
        bool found;
        PairCell pair;
        std::vector<PairCell> res ;
        for (int i=0; i < _f._filtration.size(); ++i)
        {
            this->progress_bar(i, _f._filtration.size()) ;
            found = false ;
            pair = step_persist(found, verbose) ;
            if (found)
                res.push_back(pair) ;
        }
        
        // Compute "infinite" holes
        vector<vector<int> > criticals(this->get_flag(CRITICAL)) ;
        for (int q=0; q < criticals.size(); ++q)
        {
            for (int i : criticals.at(q))
            {
                // i : persistence index
                const PairCell p = {i, -1, q} ;
                const int ki(_per_to_K.at(q).at(i)) ; // K index
                const CellDim c(ki,q) ;
                const int ti(_f._cell_to_t.at(c)) ;
                const DegType di(_f._deg.at(i)) ;
                FiltrIndexPerInterval per_int(ti,ti-1) ;
                CellDim inf(-1,q+1) ;
                CellsPerInterval per_int_cell(c,inf) ;
                DegreePerInterval per_deg_int(di,di-1) ;
                PerHole hole(per_int, per_int_cell, per_deg_int) ;
                _persist.push_back(hole) ;
                
                // If export is on, store export data
                if (_with_export)
                    export_hdvf_persistence_pair(p) ;
            }
        }
        return res;
    }
    
    /** \brief Get the "with_export" boolean flag.
     *  If the flag is `true`, homology/cohomology generators and corresponding PSC labels are exported for each persistent interval of positive duration.
     */
    bool with_export () { return _with_export ; }
    
    /** \brief Get a constant reference on the filtration
     */
    const Filtration& get_filtration() { return _f; }
    
    /** \brief Compute the (degree) duration of a persistent interval (ie. persistent hole)
     *
     * By definitions of "default" values, infinite intervals have a duration of -1.
     *
     * \param[in] hole Persistent interval considered for duration computation
     */
    DegType hole_duration (const PerHole hole) const
    {
        const DegreePerInterval per_int_deg(std::get<2>(hole)) ;
        return per_int_deg.second - per_int_deg.first ;
    }
    
    /** \brief Overload of operator<< for Hdvf_persistence.
     *
     * Prints out finite and infinite persistence intervals.
     *
     * \param[in] out Reference to an out stream.
     * \param[in] per_hdvf Constant reference on the Hdvf_persistence to print.
     */
    friend ostream& operator<< (ostream& out, const Hdvf_persistence& per_hdvf)
    {
        int i = 0 ;
        for (PerHole hole : per_hdvf._persist)
        {
            if (abs(per_hdvf.hole_duration(hole)) > 0)
                out << i << " --- duration : " << per_hdvf.hole_duration(hole) << " -- " << hole << endl ;
            ++i ;
        }
        return out ;
    }
    
    /** \brief Print informations related to the filtration.
     *
     * Prints out the filtration and associated permutations (_K_to_per and _per_to_K) between indices of cells in each dimension in the basis of K and indices along the filtration.
     *
     * \param[in] out Reference to an out stream.
     */
    ostream& print_hdvf_persistence_info (ostream& out)
    {
        out << "Filtration: " << _f << endl ;
        out << "_K_to_per and _per_to_K" << endl ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            out << "-> dim " << q << endl ;
            out << "index_per -(_per_to_K)-> index_K -(_K_to_per)-> index_per" << endl ;
            for (int i=0; i<this->_K.nb_cells(q); ++i)
            {
                const int id_K(_per_to_K.at(q).at(i)) ;
                out << i << " -> " << id_K << " -> " << _K_to_per.at(q).at(id_K) << endl ;
            }
        }
        return out ;
    }
    
    
    /**
     * \brief Export Primary/Secondary/Critical labels (e.g. for vtk export).
     *
     * The method exports the labels of every cells in each dimension.
     * Encoding used:
     * - PRIMARY: -1
     * - SECONDARY: +1
     * - CRITICAL: 0
     *
     * \returns A vector containing, for each dimension, the vector of labels by cell index.
     */
    virtual vector<vector<int> > export_psc_labels () const
    {
        vector<vector<int> > labels(this->_K.dim()+1) ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            for (int i = 0; i<this->_K.nb_cells(q); ++i)
            {
                const int id_per(_K_to_per.at(q).at(i)) ;
                if (id_per <=_t_dim.at(q))
                {
                    if (this->_flag.at(q).at(id_per) == PRIMARY)
                        labels.at(q).push_back(-1) ;
                    else if (this->_flag.at(q).at(id_per) == SECONDARY)
                        labels.at(q).push_back(1) ;
                    else if (this->_flag.at(q).at(id_per) == CRITICAL)
                        labels.at(q).push_back(0) ;
                    else // NONE
                        labels.at(q).push_back(2) ;
                }
                else
                    labels.at(q).push_back(2) ;
            }
        }
        return labels ;
    }
    
    /**
     * \brief Export homology generators associated to `cell` (critical cell) of dimension  `q` (e.g. forx@ vtk export).
     *
     * The method exports the chain \f$g(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell` and dimension `q`.
     *
     * \returns A column-major chain.
     */
    virtual CChain export_homology_chain (int cell, int q) const
    {
        if ((q<0) || (q>this->_K.dim()))
            throw "Error : export_homology_chain with dim out of range" ;
        //        if (_K_to_per.at(dim).at(cell) > _t_dim.at(dim))
        //            throw "Error : export_homology_chain with 'future' cell wrt persistence" ;
        
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            // Get g(cell, dim) with per indices
            CChain g_cell(OSM::get_column(this->_G_col.at(q), cell)) ;
            // Add 1 to the cell
            g_cell[cell] = 1 ;
            // Compute the chain with _K indices
            CChain g_cell_K(g_cell.dimension()) ;
            for (typename CChain::const_iterator it = g_cell.begin(); it != g_cell.end(); ++it)
            {
                const int i(_per_to_K.at(q).at(it->first)) ;
                g_cell_K[i] = it->second ;
            }
            return g_cell_K ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
    
    /**
     * \brief Export cohomology generators associated to `cell` (critical cell) of dimension  `q` (used by vtk export).
     *
     * The method exports the chain \f$f^\star(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell` and dimension `q`.
     *
     * \returns A column-major chain.
     */
    virtual CChain export_cohomology_chain (int cell, int q) const
    {
        if ((q<0) || (q>this->_K.dim()))
            throw "Error : export_homology_chain with dim out of range" ;
        //        if (_K_to_per.at(dim).at(cell) > _t_dim.at(dim))
        //            throw "Error : export_homology_chain with 'future' cell wrt persistence" ;
        
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            // Get fstar(cell, dim) with per indices
            RChain fstar_cell(OSM::get_row(this->_F_row.at(q), cell)) ;
            // Add 1 to the cell
            fstar_cell[cell] = 1 ;
            // Compute the cofaces of the chain with _K indices
            if (q < this->_K.dim())
            {
                CChain fstar_cofaces(this->_K.nb_cells(q+1)) ;
                for (typename RChain::const_iterator it = fstar_cell.begin(); it != fstar_cell.end(); ++it)
                {
                    // Set the cofaces of indices_K(it->first) in dimension dim+1
                    // belonging to _K(_t)
                    const int i(_per_to_K.at(q).at(it->first)) ;
                    RChain cofaces(this->_K.cod(i,q)) ;
                    for (typename RChain::const_iterator it2 =  cofaces.cbegin(); it2 != cofaces.cend(); ++it2)
                    {
                        const int id(it2->first) ;
                        if (_K_to_per.at(q+1).at(id) <=_t_dim.at(q+1))
                            fstar_cofaces[id] = 1 ;
                    }
                }
                return fstar_cofaces ;
            }
            else
                return CChain(0) ;
        }
        else
            throw "Error : trying to export g_chain without proper HDVF option" ;
    }
    
    /*! \brief Iterator over (finite) persistent intervals.
     *
     * Iterate over persistent intervals of finite degree duration.
     * If `discard_small` is true (which is the default), the iterator discards persistent intervals with a null degree duration (that is, small persistent holes).
     */
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = PerIntervalInformation;
        
        /*! \brief Iterator constructor
         *
         * \param[in] per_hdvf Constant reference over the Hdvf_persistence iterated.
         * \param[in] i The initial index.
         * \param[in] discard_small If `true` (default), only persistent intervals of (strictly) positive degree duration are iterated. Otherwise, all persistent intervals are iterated.
         */
        iterator(const Hdvf_persistence& per_hdvf, int i=0, bool discard_small = true) : _i(i), _per_hdvf(per_hdvf), _discard_small(discard_small)
        {
            if(_discard_small)
            {
                // Iterate only over holes of duration > 0
                while ((_i<_per_hdvf._persist.size()) && (_per_hdvf.hole_duration(_per_hdvf._persist.at(_i)) == 0))
                {
                    ++_i ;
                }
            }
        }
        
        /** \brief Returns the `discard_small` flag. */
        bool has_discard_small() { return _discard_small; }
        
        // Operators
        /*! \brief Iterator dereference
         *
         * \returns A `PerIntervalInformation` structure containing the information of the current persistence interval.
         */
        value_type operator*() const
        {
            PerIntervalInformation res ;
            res.hole = _per_hdvf._persist.at(_i) ;
            res.labelsPSC = _per_hdvf._export_labels.at(_i) ;
            if (_per_hdvf._hdvf_opt & (OPT_G | OPT_FULL))
            {
                res.g_chain_sigma = _per_hdvf._export_g.at(_i).first ;
                res.g_chain_tau = _per_hdvf._export_g.at(_i).second ;
            }
            if (_per_hdvf._hdvf_opt & (OPT_F | OPT_FULL))
            {
                res.fstar_chain_sigma = _per_hdvf._export_fstar.at(_i).first ;
                res.fstar_chain_tau = _per_hdvf._export_fstar.at(_i).second ;
            }
            return res ;
        }
        
        /**
         * \brief Prefix incrementation. Finds next persistent interval.
         *
         * If `discard_small` is true, the iterator searches next persistent interval with a (strictly) positive degree duration, otherwise, the iterator returns next persistent interval.
         *
         * \returns The reference to the current iterator.
         */
        iterator& operator++()
        {
            ++_i;
            if (_discard_small)
            {
                // Iterate only over holes of duration > 0
                while ((_i<_per_hdvf._persist.size()) && (_per_hdvf.hole_duration(_per_hdvf._persist.at(_i)) == 0))
                {
                    ++_i ;
                }
            }
            return *this;
        }
        
        /**
         * \brief Postfix incrementation. Finds the next not-null index.
         * \returns The pre-incremented iterator.
         */
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        
        /**
         * \brief Equality check.
         * \returns True if the indices are equal.
         */
        friend bool operator== (const iterator& a, const iterator& b) { return a._i == b._i; };
        
        /**
         * \brief Inequality check.
         * \returns True if the indices are different.
         */
        friend bool operator!= (const iterator& a, const iterator& b) { return a._i != b._i; };
        
    private:
        int _i ; // Index along _persist
        const Hdvf_persistence& _per_hdvf ; // per_hdvf iterated
        const bool _discard_small ;
    };
    
    /**
     * \brief Iterator to the beginning of persistent intervals.
     *
     * \param[in] discard_small If `true`, the iterator visits only persistent intervals of (strictly) positive degree duration. Otherwise, visit all persistent intervals.
     *
     * \returns The iterator to the beginning of the chains indices.
     */
    iterator begin(bool discard_small = true) { return iterator(*this, 0, discard_small) ; }
    
    /**
     * \brief Iterator to the ending of the chains indices.
     *
     * \returns The iterator to the ending of the chains indices.
     */
    iterator end() { return iterator(*this, _persist.size()) ; }

private:
    /** \brief Export current persistent pair informations.
     *
     * This method exports and saves PSC labels, homology and cohomology generators of the current persistent pair `p` to corresponding private vectors (_export_label, _export_g, _export_fstar).
     * The method is invoqued on the result of `find_pair_A` before pairing cells with the A operation.
     */
    void export_hdvf_persistence_pair(PairCell p)
    {
        // Export labels
        ExpLabels labels(this->export_psc_labels()) ;
        _export_labels.push_back(labels) ;
        // Export g (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            CChain chain_sigma(export_homology_chain(p.sigma, p.dim)) ;
            CChain chain_tau ;
            if (p.tau >= 0)
                chain_tau = export_homology_chain(p.tau, p.dim+1) ;
            _export_g.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
        // Export fstar (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            CChain chain_sigma(export_cohomology_chain(p.sigma, p.dim)) ;
            CChain chain_tau ;
            if (p.tau >= 0)
                chain_tau = export_cohomology_chain(p.tau, p.dim+1) ;
            _export_fstar.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
    }
    
    /** \brief Export empty persistent pair.
     * Export empty persistence information (for "small" discarded pairs).
     */
    void export_hdvf_persistence_pair()
    {
        // Export labels
        ExpLabels labels ;
        _export_labels.push_back(labels) ;
        // Export g (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_G))
        {
            CChain chain_sigma, chain_tau ;
            _export_g.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
        // Export fstar (according to options)
        if (this->_hdvf_opt & (OPT_FULL | OPT_F))
        {
            CChain chain_sigma, chain_tau ;
            _export_fstar.push_back(std::pair<CChain,CChain>(chain_sigma, chain_tau)) ;
        }
    }
    
    /**
     * \brief Find a valid PairCell for A for persistent homology.
     *
     * The function searches, at a given time \f$t\f$ in the filtration, the youngest critical cell \f$\gamma'\f$ forming a valid pair with the cell \f$\gamma\f$. Hence, \f$(\gamma', \gamma)\f$ valid pair is a valid pair
     * (ie. such that \f$\langle \partial(\gamma), \gamma' \rangle\f$ invertible).
     *
     * \param[in] found Reference to a boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    PairCell find_pair_A(bool &found) ;
    
    /**
     * \brief Step forward along the filtration.
     *
     * Searches a possible persistent pair for A with `find_pair_A`, apply A and update internal structures.
     */
    PairCell step_persist(bool& found, bool verbose = false) ;
} ;


template<typename CoefficientType, typename ComplexType, typename DegType, typename FiltrationType>
Hdvf_persistence<CoefficientType, ComplexType, DegType, FiltrationType>::Hdvf_persistence(const ComplexType& K, const FiltrationType& f, int hdvf_opt, bool with_export) : Hdvf_core<CoefficientType,ComplexType, OSM::Chain, OSM::SubSparseMatrix>(K,hdvf_opt), _f(f), _with_export(with_export), _t(-1)
{
    // Initialisation of _t_dim, _K_to_per and _per_to_K
    _t_dim.resize(this->_K.dim()+1, 0) ;
    _K_to_per.resize(this->_K.dim()+1) ;
    _per_to_K.resize(this->_K.dim()+1) ;
    for (int q=0; q<=this->_K.dim(); ++q)
    {
        _K_to_per.at(q).resize(this->_K.nb_cells(q)) ;
        _t_dim.at(q) = -1 ;
    }
    
    for (int i = 0; i<_f._filtration.size(); ++i)
    {
        const CellDim c(_f._filtration.at(i));
        const int q(c.second) ;
        const int ind_K_i(c.first) ;
        const int ind_per_i(_per_to_K.at(q).size()) ;
        _per_to_K.at(q).push_back(ind_K_i) ;
        _K_to_per.at(q).at(ind_K_i) = ind_per_i ;
    }
    
    // Init _masks
    _masks.resize(this->_K.dim()+1) ;
    for (int q=0; q<this->_K.dim()+1; ++q)
    {
        _masks.at(q) = OSM::Bitboard(this->_K.nb_cells(q)) ;
    }
    
    // Init boundary matrices
    vector<CMatrix> _DD_per(this->_K.dim()+1) ;
    
    // Copy _DD_col with filtration order (for dimensions q>0)
    for (int q = 0 ; q <= this->_K.dim(); ++q)
    {
        const std::pair<int, int> s(this->_DD_col.at(q).dimensions()) ;
        _DD_per.at(q) = CMatrix(s.first, s.second) ;
    }
    // Set empty mask for _DD_per[0]
    _DD_per.at(0).complement();
    
    for (int q = 1 ; q <= this->_K.dim(); ++q)
    {
        // Cross _DD_col.at(q) and set _DD_per.at(q) coefficients on the fly
        for (OSM::Bitboard::iterator it_col = this->_DD_col.at(q).begin(); it_col != this->_DD_col.at(q).end(); ++it_col)
        {
            const int j(*it_col) ;
            const CChain& col(OSM::cget_column(this->_DD_col.at(q), j)) ;
            for (typename CChain::const_iterator it = col.begin(); it != col.end(); ++it)
            {
                const int i(it->first) ;
                const CoefficientType v(it->second) ;
                // Cells in the _K basis : i(dim q-1) / j(dim q)
                // Convert to indices in the persistent order
                const int pi(_K_to_per.at(q-1).at(i)), pj(_K_to_per.at(q).at(j)) ;
                _DD_per.at(q).set_coef(pi, pj, v) ;
            }
        }
    }
    this->_DD_col = _DD_per ;
    
    // Init _DD_col mask (empty for all cells)
    for (int q = 1 ; q <= this->_K.dim(); ++q)
        this->_DD_col.at(q).set_sub(_masks.at(q)) ;
}

template<typename CoefficientType, typename ComplexType, typename DegType, typename FiltrationType>
PairCell Hdvf_persistence<CoefficientType, ComplexType, DegType, FiltrationType>::find_pair_A(bool &found)
{
    PairCell p ;
    // Get current cell (in the basis K)
    CellDim c(_f._filtration.at(_t)) ;
    const int q(c.second), sigma(_K_to_per.at(q).at(c.first))  ;
    // Search for pairing
    found = false;
    
    if (q >= 1)
    {
        // Compute bounded max while iterating over the Chain
        const CChain& tmp2(OSM::cget_column(this->_DD_col.at(q), sigma)) ;
        std::size_t tmax = _t_dim.at(q-1) ;
        std::size_t i ;
        for (typename CChain::const_iterator it = tmp2.cbegin(); it != tmp2.cend(); ++it)
        {
            if ((it->first <= tmax) && (abs(it->second) == 1)) // possible pairing
            {
                if (!found) // for first cell met
                {
                    found = true ;
                    i = it->first ;
                }
                else // update the i with the maximum cell index (persistance)
                {
                    if (it->first > i)
                        i = it->first ;
                }
            }
        }
        if (found)
        {
            p.sigma = i ;
            p.tau = sigma ;
            p.dim = q-1 ;
        }
    }
    return p;
}

template<typename CoefficientType, typename ComplexType, typename DegType, typename FiltrationType>
PairCell Hdvf_persistence<CoefficientType, ComplexType, DegType, FiltrationType>::step_persist(bool& found, bool verbose)
{
    // Compute next persistent pair
    
    ++_t ; // Step forward in the filtration
    const int q_current(_f._filtration.at(_t).second) ; // Get the dimension of the new current cell
    ++_t_dim.at(q_current) ; // Update time in the dimension of the current cell
    _masks.at(q_current).setOn(_t_dim.at(q_current)) ; // Update mask accordingly
    this->_DD_col.at(q_current).set_bitOn(_t_dim.at(q_current)) ; // Update _DD_col mask
    
    // Search for pairing
//    bool found ;
    PairCell p(find_pair_A(found)) ;
    if (found)
    {
        // Corresponding persistent interval
        const int q(p.dim) ;
        // indices of both cells in the _K basis
        const int ki(_per_to_K.at(q).at(p.sigma)), kj(_per_to_K.at(q+1).at(p.tau)) ;
        CellDim ci(ki, q), cj(kj, q+1) ; // cells of the interval - in the K basis
        int ti(_f._cell_to_t.at(ci)), tj(_f._cell_to_t.at(cj)) ; // times of the interval
        FiltrIndexPerInterval interval(ti, tj) ;
        CellsPerInterval interval_cells(ci, cj) ;
        DegreePerInterval interval_deg(_f._deg.at(ti), _f._deg.at(tj)) ;
        PerHole hole(interval, interval_cells, interval_deg) ;
        // Add this interval
        _persist.push_back(hole) ;
        
        // If export is on, store export data for significant persistant intervals
        if (_with_export)
        {
            if ((interval_deg.second-interval_deg.first)>0)
                export_hdvf_persistence_pair(p) ;
            else
                export_hdvf_persistence_pair() ;
        }

        
        // Prepare for next step
        this->A(p.sigma, p.tau, p.dim) ; // Update the reduction
        if (verbose)
        {
            std::cout << "A : " << p.sigma << " - " << p.tau << " (dim " << p.dim << ")" << std::endl ;
            this->print_matrices(std::cout) ;
        }
    }
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_PERSISTENCE_H
