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

#ifndef CGAL_HDVF_SUB_CHAIN_COMPLEX_MASK_H
#define CGAL_HDVF_SUB_CHAIN_COMPLEX_MASK_H

#include <CGAL/license/HDVF.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <memory>
#include <CGAL/OSM/OSM.h>
#include <CGAL/OSM/Bitboard.h>
#include <CGAL/HDVF/Sub_sparse_matrix.h>

namespace CGAL {
namespace HDVF {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Sub_chain_complex_mask` is a technical class implementing a sub chain complex. A sub chain complex \f$A\f$ of a chain complex \f$K\f$ is a subset \f$A\subseteq K\f$ such that the restricted boundary operator \f$\partial_A = \partial_K|_A\f$ still satisfies \f$\partial_A^2 = 0\f$.

 The `Sub_chain_complex_mask` class is used to compute  reduced homology. This class is based on a set of bitboard masks (one in each dimension) used to define sub chain complexes and their associated reduction encoded in sub-sparse matrices (`OSM::Sub_sparse_matrix` class). Technically, `Sub_chain_complex_mask` are used to partially screen chain complexes and chains in associated boundary matrices, and hence compute homology "locally".

 \warning For efficiency reasons, when a `Sub_chain_complex_mask` is used to screen sparse matrices  (with the `screen_matrices` method), screening is **only** performed on the major direction of matrices (thus column-major matrices are restricted over columns and row-major matrices are restricted over rows). Iterators are restricted accordingly. But chains themselves are not restricted.<br>
 However, if `A` is a proper sub-complex of `K` (that is, closed with respect to faces), chains automatically comply with the screening. Indeed, for any \f$q\f$-cell \f$\sigma\in A\f$ (thus the corresponding bit is on in the mask), all the faces of \f$\sigma\f$ also belong to \f$A\f$. Thus for any \f$q-1\f$-cell \f$\tau\f$ in the boundary of \f$\sigma\f$ (that is, such that \f$\langle\partial_k(\sigma),\tau\rangle\neq 0\f$), \f$\tau\f$ belongs to \f$A\f$ (and thus the corresponding bit in the mask is also on).

 \tparam CoefficientType a model of the `Ring` concept.
 \tparam ComplexType a model of the `AbstractChainComplex` concept (type of the chain complex screened by `Sub_chain_complex_mask`).
 */

template <typename CoefficientType, typename ComplexType>
class Sub_chain_complex_mask
{
protected:
    /** \brief Dimension of the underlying complex. */
    int _dim ;
    /** \brief Vector of bitboards coding masks in each dimension. */
    std::vector<OSM::Bitboard> _sub ;
    /** \brief Number of cells in the mask in each dimension. */
    std::vector<int> _nb_cells ;
    /** \brief Constant reference to the underlying complex. */
    const ComplexType& _K ;
    /** \brief Is full sub_complex.
     * This boolean flag is true is all bits in the bitboards are on. */
    bool _full ;
private:
    /** Load a set of cells and close the enumeration down.
     *
     * The method activates bits corresponding to the set of cells encoded by `cells` and recursively activates all bits corresponding to their faces.
     *
     * \param[in] cells A vector containing, in each dimension, a vector of cells indexes.
     */
    void down_closure (const std::vector<std::vector<int> > &cells)
    {
        std::vector<std::vector<int> > faces(cells.size()) ;
        bool rec_needed = false ;
        for (int q = 0; q < cells.size(); ++q)
        {
            for (int i : cells.at(q))
            {
                if (!_sub.at(q).isOn(i))
                {
                    _sub.at(q).setOn(i) ;
                    ++_nb_cells.at(q) ;
                    if (q>0) // Then faces must be considered
                    {
                        // Add all its faces to faces.at(q-1)
                        OSM::Sparse_chain<CoefficientType, OSM::COLUMN> bnd(_K.d(i,q)) ;
                        for (typename OSM::Sparse_chain<CoefficientType, OSM::COLUMN>::const_iterator it = bnd.cbegin(); it != bnd.cend(); ++it)
                        {
                            const int c(it->first) ;
                            if (!_sub.at(q-1).isOn(c))
                            {
                                faces.at(q-1).push_back(c) ;
                                rec_needed = true ;
                            }
                        }
                    }
                }
            }
        }
        // Call recursively
        if (rec_needed)
            down_closure(faces, _K);
    }

    /** Load a set of cells (without closure).
     *
     * The method activates bits corresponding to the set of cells encoded by `cells`.
     *
     * \param[in] cells A vector giving, for each dimension, a vector of cells indexes.
     */
    void set_cells (const std::vector<std::vector<int> > &cells)
    {
        for (int q = 0; q < cells.size(); ++q)
        {
            for (int i : cells.at(q))
            {
                _sub.at(q).setOn(i) ;
                ++_nb_cells.at(q) ;
            }
        }
    }

public:
    /** \brief Constructor from a complex.
     *
     * Build a `Sub_chain_complex_mask` associated to `K` with all bits set to 1 in the masks if `full` is true, and all bits set to 0 otherwise.
     *
     * \param[in] K A constant reference to the underlying complex.
     * \param[in] full Build full / empty masks (default: full).
     */
    Sub_chain_complex_mask(const ComplexType& K, bool full=true) : _K(K)
    {
        _dim = K.dimension() ;
        _sub.resize(_dim+1) ;
        _nb_cells.resize(_dim+1) ;
        // Create Bitboards
        for (int q=0; q<=K.dimension(); ++q)
        {
            if (full)
            {
                _sub.at(q) = OSM::Bitboard(K.number_of_cells(q),false) ; // all cells to 1
                _full = true ;
            }
            else
                _sub.at(q) = OSM::Bitboard(K.number_of_cells(q),true) ; // all cells to 0
            _nb_cells.at(q) = K.number_of_cells(q) ;
        }
    }


    /** \brief Constructor from an enumeration of cells.
     *
     * Build masks associated to the underlying complex `K` with all bits corresponding to `cells` (and their faces if `close` is true) set to 1.
     *
     * \param[in] K A constant reference to the underlying complex.
     * \param[in] cells A vector containing, in each dimension, a vector of cells indexes.
     * \param[in] close If this boolean is true, the faces of `cells` are also set to 1.
     */
    Sub_chain_complex_mask(const ComplexType& K, const std::vector<std::vector<int> > &cells, bool close = true) : _K(K)
    {
        _dim = K.dimension() ;
        _sub.resize(_dim+1) ;
        _nb_cells.resize(_dim+1) ;
        _full = true ;
        // Create Bitboards
        for (int q=0; q<=K.dimension(); ++q)
        {
            _sub.at(q) = OSM::Bitboard(K.number_of_cells(q)) ;
        }

        // Set bits
        if (close)
            down_closure(cells, K) ;
        else
            set_cells(cells, K) ;

        for (int q=0; q<=K.dimension(); ++q)
        {
            _full = (_full && (_nb_cells.at(q) == _sub.at(q).size())) ;
        }
    }

    /**
     * \brief Constructor by copy.
     *
     * Builds a `Sub_chain_complex_mask` by copy from another.
     *
     * \param[in] otherToCopy An initial `Sub_chain_complex_mask`.
     */
    Sub_chain_complex_mask(const Sub_chain_complex_mask& otherToCopy) : _K(otherToCopy._K)
    {
        _dim = otherToCopy._dim ;
        _nb_cells = otherToCopy._nb_cells ;
        _full = otherToCopy._full ;
        _sub = otherToCopy._sub ;
    }

    /** \brief Assignment operator
     *
     * \warning The operator argument must provide a sub chain complex mask over the *same underlying chain complex*. It not so, the assignment will throw an exception.
     *
     * \param[in] otherToCopy A `Sub_chain_complex_mask` copied into `this` (`otherToCopy` and `this` must have the same underlying chain complex).
     */
    Sub_chain_complex_mask& operator= (const Sub_chain_complex_mask& otherToCopy)
    {
        // Check that otherToCopy and current mask share the same chain complex
        if (_K.get_id() != otherToCopy._K.get_id() )
            throw("Error, operator= can only copy mask over the same chain complex.");

        // Perform assignment
        _dim = otherToCopy._dim ;
        _sub = otherToCopy._sub ;
        _nb_cells = otherToCopy._nb_cells ;
        _full = otherToCopy._full ;
        return *this;
    }

    /** \brief Returns the complement of the mask.
     *
     * The method return a new `Sub_chain_complex_mask` containing the complement of the current mask (0 and 1 bits in the mask are exchanged).
     */
    Sub_chain_complex_mask complement()
    {
        Sub_chain_complex_mask cSub(*this) ;
        // Complement bitboards
        for (int q=0; q<=_dim; ++q)
        {
            cSub._sub.at(q).bit_not() ;
            cSub._nb_cells.at(q) = _sub.at(q).size() - _nb_cells.at(q) ;
            cSub._full = ~_full ;
        }
        return cSub ;
    }

    /** \brief Gets a bit of the mask (bit i in dimension q). */
    inline bool get_bit(int q, int i) const
    {
        return _sub.at(q).isOn(i) ;
    }

    /** \brief Sets a bit to 1 (bit i in dimension q). */
    inline void set_bit_on(int q, int i)
    {
        if (!_sub.at(q).isOn(i))
        {
            _nb_cells[q]++ ;
            _sub.at(q).setOn(i) ;
        }
    }

    /** \brief Sets a bit to 0 (bit i in dimension q). */
    inline void set_bit_off(int q, int i)
    {
        if (_sub.at(q).isOn(i))
        {
            _nb_cells[q]--;
            _sub.at(q).setOff(i) ;
        }
    }

    /** \brief Gets the bitboards of the sub chain complex.
     *
     * Returns a constant reference to the vector of bitboards in each dimension.
     */
    inline const std::vector<OSM::Bitboard>& get_bitboard() const
    {
        return _sub ;
    }

    /** \brief Gets the bitboard of the sub chain complex in dimension q. */
    inline const OSM::Bitboard& get_bitboard(int q) const
    {
        return _sub.at(q) ;
    }

    /** \brief Screens a sequence of Sub_sparse_matrix (in each dimension).
     *
     * Given a sequence of matrices (vector of `Sub_sparse_matrices`) sets the masks of `Sub_sparse_matrices` in each dimension to the current `Sub_chain_complex_mask`.
     *
     * \warning For efficiency, screening is performed on the major direction of matrices (so along columns for column-major matrices and along row for row-major matrices).
     *
     * \param[in] matrices A vector of Sub_sparse_matrix (in each dimension).
     */
    template <typename CT, int CTF>
    void screen_matrices(std::vector<OSM::Sub_sparse_matrix<CT, CTF> >& matrices)
    {
        for (int q=0; q<=_K.dimension(); ++q)
        {
            matrices.at(q).set_sub(_sub.at(q)) ;
        }
    }

    /**
     * \brief Restricts a chain in a given dimension to the sub chain complex maks.
     *
     * Nullify all coefficients out of the mask.
     *
     * \param[in] chain The chain restricted.
     * \param[in] q Dimension of the chain.
     */
    template <typename CT, int CTF>
    inline void screen_chain(OSM::Sparse_chain<CT, CTF>& chain, int q)
    {
        vector<size_t> indices;
        for (typename OSM::Sparse_chain<CT, CTF>::iterator it = chain.begin(); it != chain.end(); ++it)
        {
            if (!get_bit(q, it->first))
                indices.push_back(it->first);
        }
        chain/=indices;
    }

    /*! \brief Overload of the `<<`operator for `Sub_chain_complex_mask`.
     */
    friend std::ostream & operator << (std::ostream & out, const Sub_chain_complex_mask & sub)
    {
        for (int q = 0; q <= sub._dim; ++q)
            out << "dim " << q << " : " << sub._sub.at(q) << std::endl ;
        return out ;
    }
};

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif // CGAL_HDVF_SUB_CHAIN_COMPLEX_MASK_H
