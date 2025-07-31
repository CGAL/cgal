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

#ifndef CGAL_HDVF_FILTRATION_CORE_H
#define CGAL_HDVF_FILTRATION_CORE_H

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <functional>

namespace CGAL {
namespace HDVF {

/*!
 \ingroup PkgHDVFAlgorithmClasses
 
 The class `Filtration_core` implements data structures and methods required by the `CGAL::Filtration` concept.
 
 By definition, a filtration over a chain complex `K` is a sequence of sub complexes \f$\{K_t\,;\, t\in[d_1, \ldots, d_N]\}\f$, where \f$d_i\f$ are scalars, such that:
 
 - \f$\forall d\leqslant d'\f$, \f$K_d\subseteq K_{d'}\f$.
 - \f$|K_{d_{i+1}}| = |K_{d_{i}}|+1\f$.
    - As a consequence, any cell \f$\sigma\in K\f$, is "added" at a given *time* \f$i\f$ and the corresponding \f$d_i\f$ is called its *degree* (denoted by \f$\mathrm{deg}(\sigma)\f$).
 - \f$K_{d_1}\f$ is an empty complex.
 - \f$K_{d_N} = K\f$.
 - For any cell \f$\sigma\in K\f$, if \f$\tau\f$ is a proper face of \f$\sigma\f$, then \f$\mathrm{deg}(\tau)\leqslant \mathrm{deg}(\sigma)\f$.
 
 In consequence, a filtration associated to a complex `K`:
 - maps each cell to a scalar called its degree
 - orders the cells of `K` (each cell has a unique index along the filtration called *time* along the filtration) such that:
    - the index of a cell is larger that the indices of its faces,
    - the degree map is increasing along the filtration.
 
 The class `Filtration_core` provides elementary constructors and methods used in derived filtrations.

 \cgalModels{Filtration}
 
 \tparam CoefficientType a model of the `Ring` concept (ring used for homology computation).
 \tparam ComplexType a model of the `AbstractChainComplex` concept (type of the underlying chain complex).
 \tparam DegreeType the scalar type of degrees.
 */

template <typename CoefficientType, typename ComplexType, typename DegreeType>
class Filtration_core
{
public:
    /*! \brief Type for indexing uniquely a cell.
     * - First element of the pair: index of the cell.
     * - Second element of the pair: dimension of the cell.
     */
    typedef std::pair<size_t, int> CellDim ;
    
    /*! \brief Type of value returned by the iterator.
     */
    typedef struct {
        CellDim cell_dim ;
        DegreeType degree ;
    } FiltrationIterValue ;
    
protected:
    /** \brief Constant reference to the underlying chain complex. */
    const ComplexType& _K ;
    /** \brief Vector of cells of the filtration (full ordering of cells). */
    std::vector<CellDim> _filtration ;
    /** \brief Vector of degrees of cells along the filtration. */
    std::vector<DegreeType> _deg ;
    
    /** \brief Map from cells to their index in the filtration. */
    std::map<CellDim,size_t> _cell_to_t ;
    
    /*!
     Type of column-major sparse matrices
     */
    typedef OSM::Sparse_matrix<CoefficientType,OSM::COLUMN> CMatrix ;
    
    /*!
     Type of row-major sparse matrices
     */
    typedef OSM::Sparse_matrix<CoefficientType,OSM::ROW> RMatrix ;
    /*!
     Type of column-major chains
     */
    typedef OSM::Sparse_chain<CoefficientType,OSM::COLUMN> CChain ;
    
    /*!
     Type of row-major chains
     */
    typedef OSM::Sparse_chain<CoefficientType,OSM::ROW> RChain ;
public:
    /*! \brief Filtration_core default constructor
     *
     * Builds an "empty" filtration with `K` as underlying chain complex.
     *
     *\param[in] K A chain complex (a model of `AbstractChainComplex`), the underlying chain complex of the filtration.
     */
    Filtration_core(const ComplexType& K) : _K(K)
    {}
        
    /*! \brief Constructor from a vector of cells (ordering of cells) and an associated vector of degrees.
     *
     * The constructor check that the filtration is valid (a cell is introduced in the filtration after its faces and the degree vector is increasing) and throw an exception if not.
     * \param[in] K A chain complex (a model of `AbstractChainComplex`), the underlying chain complex of the filtration.
     * \param[in] filtration An ordering of the cells of `K` encoded as a vector of its cells.
     * \param[in] deg The (increasing) vector of cells degrees.
     */
    Filtration_core(const ComplexType& K, const std::vector<CellDim>& filtration, const std::vector<DegreeType>& deg) : _K(K), _filtration(filtration), _deg(deg)
    {
        if (!is_valid_filtration())
            throw ("Invalid filtration, Filtration_core constructor failed");
    }
    
    /**
     * \brief Constructor by copy.
     *
     * Builds a `Filtration_core` by copy from another.
     *
     * \param[in] f An initial `Filtration_core`.
     */
    Filtration_core(const Filtration_core& f) : _K(f._K), _filtration(f._filtration), _cell_to_t(f._cell_to_t) {}
    
    /** \brief Iterator over a filtration.
     *
     * Iterate the cells of the filtration by increasing index (and hence increasing degree).
     */
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = FiltrationIterValue;
        
        /** \brief Iterator constructor
         *
         * \param[in] f Constant reference on the `Filtration_core` iterated.
         * \param[in] i The initial index.
         */
        iterator(const Filtration_core& f, size_t i=0) : _i(i), _f(f) {}
        
        /*! \brief Iterator dereference
         *
         * \returns A `FiltrationIterValue` structure containing the information of the current cell and its degree.
         */
        value_type operator*() const
        {
            FiltrationIterValue res ;
            res.cell_dim = _f._filtration.at(_i) ;
            res.degree = _f._deg.at(_i) ;
            return res ;
        }
        
        /*! \brief Get the time (or index in the filtration) associated to current iterator.
         */
        size_t time () const { return _i ; }
        
        /*! \brief Get the cell (identified by its index and dimension) associated to current iterator.
         */
        CellDim cell_dim () const { return _f._filtration.at(_i); }
        
        /*! \brief Get degree associated to current iterator.
         */
        DegreeType degree () const { return _f._deg.at(_i); }
        
        /*!
         * \brief Prefix incrementation. Moves to next cell in the filtration.
         *
         * \returns The reference to the current iterator.
         */
        iterator& operator++()
        {
            ++_i;
            return *this;
        }
        
        /*!
         * \brief Postfix incrementation. Moves to next cell in the filtration.
         * \returns The pre-incremented iterator.
         */
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        
        /*!
         * \brief Equality check.
         * \returns True if the indices are equal.
         */
        friend bool operator== (const iterator& a, const iterator& b)  { return a._i == b._i; };
        
        /*!
         * \brief Inequality check.
         * \returns True if the indices are different.
         */
        friend bool operator!= (const iterator& a, const iterator& b)  { return a._i != b._i; };
        
    private:
        size_t _i ; // Index along _persist
        const Filtration_core& _f ; // Filtration_core iterated
    };
    
    /**
     * \brief Iterator to the beginning of the filtration.
     *
     * \return The iterator to the beginning of the filtration.
     */
    iterator begin() { return iterator(*this, 0) ; }
    
    /**
     * \brief Iterator to the ending of the filtration.
     *
     * \return The iterator to the ending of the filtration.
     */
    iterator end() { return iterator(*this, _filtration.size()) ; }
    
    // getters
    /*! \brief Get the filtration size.
     */
    size_t get_filtration_size () const { return _filtration.size();}
    
    /*! \brief Get the cell (that is cell index and dimension) at the index `i` of the filtration.
     */
    CellDim get_cell_dim (size_t i) const { return _filtration.at(i); }
    
    /*! \brief Get the degree of the `i`th element of the filtration.
     */
    DegreeType get_degree (size_t i) const { return _deg.at(i); }
    
    // Output filtration
    /*! \brief Overload of the `<<`operator for filtrations.
     */
    friend ostream & operator<<(ostream & out, const Filtration_core &f)
    {
        const size_t N(f._filtration.size()) ;
        for (size_t i=0; i<N; ++i)
            // Filtration_core
        {
            out << i << " -> (" << f._filtration.at(i).first << "- dim " << f._filtration.at(i).second << " , " << f._deg.at(i) << ") " << std::endl ;
        }
        return out ;
    }
    
    /**
     * \brief Export the filtration time indices.
     *
     * The method exports the time index of every cells in each dimension.
     *
     * \returns A vector containing, for each dimension, the vector of labels by cell index.
     */
    vector<vector<size_t> > export_filtration () const
    {
        vector<vector<size_t> > labels(_K.dim()+1) ;
        for (int q=0; q<=this->_K.dim(); ++q)
        {
            for (size_t i = 0; i<this->_K.nb_cells(q); ++i)
            {
                const CellDim cell(i,q) ;
                const size_t t(_cell_to_t.at(cell));
                labels.at(q).push_back(t) ;
            }
        }
        return labels ;
    }
    
    /*! \brief Check that a filtration is valid.
     * Checks that cells are ordered in increasing degrees and all cells have indices larger than their faces.
     * \returns `true` if the filtration is valid, `false` otherwise
     */
    bool is_valid_filtration() const ;
    
protected:
    /** \brief Sets _cell_to_t from the initial vector o cells. */
    void build_filtration_structure() ;
    
    /** Friend class: Hdvf_persistence. */
    template <typename CoefT, typename ComplexT, typename DegT, typename FiltrT>
    friend class Hdvf_persistence ;
};

template <typename CoefficientType, typename ComplexType, typename DegreeType>
void Filtration_core<CoefficientType, ComplexType, DegreeType>::build_filtration_structure()
{
    for (size_t i = 0; i<_filtration.size(); ++i)
    {
        const CellDim c(_filtration.at(i)) ;
        // c : filtration index i, index in the basis reordered by filtration : j
        _cell_to_t[c] = i ;
    }
}

template <typename CoefficientType, typename ComplexType, typename DegreeType>
bool Filtration_core<CoefficientType, ComplexType, DegreeType>::is_valid_filtration() const
{
    bool valid = true ;
    for (size_t i=0; i<_filtration.size() && valid; ++i)
    {
        if (i>0)
            valid = valid & (_deg.at(i) >= _deg.at(i-1)) ;
        CellDim c(_filtration.at(i)) ;
        cout << i << " -> " << c.first << " dim "<< c.second << endl ;
        const int q = c.second ;
        if (q>0)
        {
            CChain dc = _K.d(c.first, q) ;
            cout << "bnd : " << dc << endl ;
            for (typename CChain::iterator it = dc.begin(); it != dc.end() && valid; ++it)
            {
                // Faces of c
                const CellDim face(it->first,q-1) ;
                // Check if the face c belongs to the filtration
                auto it_face(_cell_to_t.find(face)) ;
                valid = valid & (it_face != _cell_to_t.end()) ;
                if (!valid)
                    cout << "face not found" << endl ;
                if (it_face != _cell_to_t.end())
                    valid = valid & (_cell_to_t[face] < i) ;
                if (!valid)
                    cout << "face " << it->first << " at time : " << _cell_to_t[face] << " with i : " << i << endl ;
            }
        }
        if (!valid)
            cout << "check failed : " << i << endl ;
    }
    return valid ;
}

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_FILTRATION_CORE_H
