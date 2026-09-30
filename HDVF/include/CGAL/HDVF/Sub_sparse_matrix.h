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
//                  Kevin Fedyna <fedyna.kevin@gmail.com>

#ifndef CGAL_OSM_SUB_SPARSE_MATRIX_H
#define CGAL_OSM_SUB_SPARSE_MATRIX_H


#include <CGAL/license/HDVF.h>

#include <CGAL/OSM/Sparse_matrix.h>
#include <CGAL/OSM/Bitboard.h>

namespace CGAL {
namespace OSM {

// Forward declaration of Sub_sparse_matrix_core
#ifndef DOXYGEN_RUNNING
template <typename CoefficientRing, int StorageFormat, template <typename, int> typename SparseChainType>
class Sub_sparse_matrix_core;
#endif

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The structure `Sub_sparse_matrix` provides a user friendly interface for sub-sparse matrices built over a given `SparseChain`model (actually a Curryfication of the `Sub_parse_matrix_core` template).

 Given the template parameter `SparseChainType` (a model of the `SparseChain` concept), corresponding sparse chains and matrices templates are given by:
 \code
 Sub_sparse_matrix<SparseChainType>:: template Sparse_chain_type
 Sub_sparse_matrix<SparseChainType>:: template Sparse_matrix_type
 \endcode

 \cgalModels{SparseMatrix}

 \tparam SparseChainType a model of `SparseChain` used to store chains of the sparse matrix (default: `OSM::Sparse_chain`).

 \warning Check if the documentation of this structure is correct.
 */

template <template <typename, int> typename SparseChainType = OSM::Sparse_chain>
struct Sub_sparse_matrix {
    template <typename CT, int SF>
    using Sparse_chain_type = SparseChainType<CT,SF>;

    template <typename CT, int SF>
    using Sparse_matrix_type = Sub_sparse_matrix_core<CT,SF, SparseChainType>;
};

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Sub_sparse_matrix_core` is a technical class implementing the concept `SparseMatrix` together with a system of masks to partially screen matrices (and restrict computations to a subset of indices *along their major direction*). This class is used to compute reduced homology (and thus to compute persistent homology and Alexander duality).

 `Sub_sparse_matrix_core` inherits `Sparse_matrix_core` structure and basically adds two bitboards:
 - one describing indices of cells belonging to the subset of indices considered (let us denote it by \f$A\f$)
 - the second  providing indices of non-empty chains in \f$A\f$

 The class does not modify linear algebra operators, but focuses on adapted iterators, output operators and methods to adjust the "mask".

 \cgalModels{SparseMatrix}

 \tparam CoefficientRing a model of the `IntegralDomainWithoutDivision` concept, providing the ring used to compute homology.
 \tparam StorageFormat an integer constant encoding the storage format of matrices (`OSM::COLUMN` or `OSM::ROW`).
 \tparam SparseChainType a model of `SparseChain` used to store chains of the sparse matrix (default: `OSM::Sparse_chain`).
*/


template <typename CoefficientRing, int StorageFormat, template <typename, int> typename SparseChainType = OSM::Sparse_chain>
class Sub_sparse_matrix_core : public Sparse_matrix_core<CoefficientRing, StorageFormat, SparseChainType> {

    /** \brief Type of the chains associated to the matrix. */
    typedef Sparse_chain<CoefficientRing, StorageFormat> Matrix_chain;

    /**
     Type of the parent class.
     */
    typedef Sparse_matrix_core<CoefficientRing, StorageFormat, SparseChainType> Base;

protected:
    /* \brief A bitboard describing subchains restriction. */
    Bitboard _subChains;

    /* \brief A bitboard containing state of each chain (restricted to subchains). */
    Bitboard _subChainsStates;

public:
    /**
     * \brief Constructor with given rows/columns sizes and mask set to `full`.
     *
     * Constructor with sizes, initializes an empty `Sub_sparse_matrix_core` of a given size along rows/columns. The constructor sets the mask to `full`.
     *
     * \param rowCount The number of rows to preallocate (default 0).
     * \param columnCount The number of columns to preallocate (default 0).
     */
    Sub_sparse_matrix_core(size_t rowCount=0, size_t columnCount=0) : Sparse_matrix_core<CoefficientRing, StorageFormat, SparseChainType>(rowCount, columnCount)
    {
        if (StorageFormat == OSM::COLUMN)
            _subChains = OSM::Bitboard(columnCount,false) ;
        else
            _subChains = OSM::Bitboard(rowCount,false) ;
        _subChainsStates = this->_chainsStates & _subChains ;
    }

    /**
     * \brief Constructor with given rows/columns sizes and a mask.
     *
     * Create a new empty `Sub_sparse_matrix_core` of given size along rows/columns and a given mask.
     *
     * \param rowCount The number of rows to preallocate.
     * \param columnCount The number of columns to preallocate.
     * \param subChain Bitboard describing the subset of indices considered as a mask.
     */
    Sub_sparse_matrix_core(size_t rowCount, size_t columnCount, const Bitboard& subChain) : Sparse_matrix_core<CoefficientRing, StorageFormat, SparseChainType>(rowCount, columnCount), _subChains(subChain), _subChainsStates(this->_chainsStates & subChain)
    {
    }

    /** \brief Copy constructor.
     */
    Sub_sparse_matrix_core(const Sub_sparse_matrix_core& otherToCopy) : Sparse_matrix_core<CoefficientRing, StorageFormat, SparseChainType>(otherToCopy), _subChains(otherToCopy._subChains), _subChainsStates(otherToCopy._subChainsStates) {}

    /** \brief Copy constructor.
     */
    Sub_sparse_matrix_core(const Sparse_matrix_core<CoefficientRing,StorageFormat, SparseChainType>& otherToCopy) : Sparse_matrix_core<CoefficientRing, StorageFormat>(otherToCopy)
    {
        if (StorageFormat == OSM::COLUMN)
            _subChains = OSM::Bitboard(otherToCopy.dimensions().second,false) ;
        else
            _subChains = OSM::Bitboard(otherToCopy.dimensions().first,false) ;
        _subChainsStates = this->_chainsStates & _subChains ;
    }


    /** \brief Iterator to the beginning of the indices of non empty chains inside the mask. */
    inline Bitboard::iterator begin() const noexcept
    {
        return _subChainsStates.begin() ;
    }

    /** \brief Iterator past-the-end of the indices of non empty chains inside the mask. */
    inline Bitboard::iterator end() const noexcept
    {
        return _subChainsStates.end() ;
    }

    /** \brief Changes the indices subset mask.
     *
     * Set a new mask encoding a new subset of indices along the major dimension.
     */
    inline void set_sub (const Bitboard& new_subChains)
    {
        _subChains = new_subChains ;
        _subChainsStates = this->_chainsStates & _subChains ;
    }

    /** \brief Adds an index to the mask.
     *
     * Set the bit encoding a given index to 1 (ie.\ add the index in the mask).
     *
     * \param index Index to turn on in the mask.
     */
    inline void set_bit_on (size_t index)
    {
        _subChains.set_on(index) ;
        if (this->_chainsStates.is_on(index))
            _subChainsStates.set_on(index) ;
    }

    /** \brief Removes an index from the mask.
     *
     * Set the bit encoding a given index to 0 (ie.\ remove the index from the mask).
     *
     * \param index Index to turn off in the mask.
     */
    inline void set_bit_off (size_t index)
    {
        _subChains.set_off(index) ;
        _subChainsStates.set_off(index) ;
    }

    /** \brief Changes the mask to its complement. */
    inline void complement() { _subChains.bit_not() ; }

    /**
     * \brief Assignment.
     */
    inline Sub_sparse_matrix_core& operator=(const Sub_sparse_matrix_core &otherToCopy)
    {
        (dynamic_cast<Sparse_matrix_core<CoefficientRing,StorageFormat, SparseChainType>&>(*this)).operator=(otherToCopy) ;
        _subChains = otherToCopy._subChains ;
        _subChainsStates = otherToCopy._subChainsStates ;
        return *this ;
    }

    /**
     * \brief Displays a `Sub_sparse_matrix_core` in the output stream.
     *
     * Displays the sparse matrix as well as its mask.
     *
     * \param stream The output stream.
     * \param matrix The matrix to display.
     *
     * \return A reference to the modified stream.
     */
    template <typename _CR, int _SF, template <typename, int> typename _SCT>
    friend std::ostream& operator<<(std::ostream &stream, const Sub_sparse_matrix_core<_CR,_SF,_SCT> &matrix);

    /**
     * \brief Tests if a `Sub_sparse_matrix_core` is null.
     *
     * The function return `true` is the `Sub_sparse_matrix_core` is null (that is, all the chains in the mask are empty) and `false` otherwise.
     */
    bool is_null() const
    {
        return (_subChainsStates.begin() == _subChainsStates.end()) ;
    }

    /**
     * \brief Tests if a `Sub_sparse_matrix_core` with chains restricted to a `Bitboard` is null.
     *
     * The function return `true` is the restricted `Sub_sparse_matrix_core` is null (that is, all the chains in the mask, restricted to the `Bitboard` `b` are null) and `false` otherwise.
     */
    bool is_null(const Bitboard& b) const
    {
        bool res(true);
        // Test along chains
        for (Bitboard::iterator it = _subChainsStates.begin(); res && (it != _subChainsStates.end()); ++it) {
            // Checks that all non zero coefficients of the chain at index *it are out of b
            const Matrix_chain& chain(this->_chains.at(*it));
            for (typename Matrix_chain::const_iterator it2 = chain.begin(); it2 != chain.end(); ++it2) {
                res = !(b.is_on(it2->first));
                if (!res)
                    std::cout << "is_null: " << *it << " - " << it2->first << " b.is_on : " << b.is_on(it2->first) << std::endl;
            }
        }
        return res ;
    }
};

template <typename _CR, int _SF, template <typename, int> typename _SCT>
std::ostream& operator<<(std::ostream &stream, const Sub_sparse_matrix_core<_CR,_SF,_SCT> &matrix) {
    stream << static_cast<const Sparse_matrix_core<_CR,_SF,_SCT>&>(matrix) ;
    stream << matrix._subChains << std::endl ;
    return stream ;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SUB_SPARSE_MATRIX_H
