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

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Sub_sparse_matrix` is a technical class implementing the concept `SparseMatrix` together with a system of masks to partially screen matrices (and restrict computations to a subset of indices *along their major direction*). This class is used to compute reduced homology (and thus to compute persistent homology and Alexander duality).

 `Sub_sparse_matrix` inherits `Sparse_matrix` structure and basically adds two bitboards:
 - one describing indices of cells belonging to the subset of indices considered (let us denote it by \f$A\f$)
 - the second  providing indices of non-empty chains in \f$A\f$

 The class does not modify linear algebra operators, but focuses on adapted iterators, output operators and methods to adjust the "mask".

 \cgalModels{SparseMatrix}

 \tparam CoefficientRing a model of the `IntegralDomainWithoutDivision` concept, providing the ring used to compute homology.
 \tparam StorageFormat an integer constant encoding the storage format of matrices (`OSM::COLUMN` or `OSM::ROW`).
*/


template <typename CoefficientRing, int StorageFormat>
class Sub_sparse_matrix : public Sparse_matrix<CoefficientRing, StorageFormat> {

protected:
    /** \brief A bitboard describing subchains restriction. */
    Bitboard _subChains;

    /** \brief A bitboard containing state of each chain (restricted to subchains). */
    Bitboard _subChainsStates;

public:
    /**
     * \brief Default constructor of a new `Sub_sparse_matrix` (with given rows/columns sizes and mask set to `full`).
     *
     * Default constructor. Constructor with sizes, initialize an empty Sub_sparse_matrix of type `StorageFormat` with coefficients of type `CoefficientRing`, a given size along rows/columns. The constructor sets the mask to `full`.
     *
     * \param[in] rowCount The number of rows to preallocate (default 0).
     * \param[in] columnCount The number of columns to preallocate (default 0).
     */
    Sub_sparse_matrix(size_t rowCount=0, size_t columnCount=0) : Sparse_matrix<CoefficientRing, StorageFormat>(rowCount, columnCount)
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
     * Create a new empty `Sub_sparse_matrix` of type `StorageFormat` with coefficients of type `CoefficientRing`, a given size along rows/columns and a given mask.
     *
     * \param[in] rowCount The number of rows to preallocate.
     * \param[in] columnCount The number of columns to preallocate.
     * \param[in] subChain Bitboard describing the subset of indices considered as a mask.
     */
    Sub_sparse_matrix(size_t rowCount, size_t columnCount, const Bitboard& subChain) : Sparse_matrix<CoefficientRing, StorageFormat>(rowCount, columnCount), _subChains(subChain), _subChainsStates(this->_chainsStates & subChain)
    {
    }

    /** \brief Copy constructor from another `Sub_sparse_matrix`.
     *
     * Create a new empty `Sub_sparse_matrix` from another of type `StorageFormat` with coefficients of type `CoefficientRing`, a given size along rows/columns and a given mask.
     *
     * \param[in] otherToCopy `Sub_sparse_matrix` copied into `this`.
     */
    Sub_sparse_matrix(const Sub_sparse_matrix& otherToCopy) : Sparse_matrix<CoefficientRing, StorageFormat>(otherToCopy), _subChains(otherToCopy._subChains), _subChainsStates(otherToCopy._subChainsStates) {}

    /** \brief Copy constructor from `Sparse_matrix`.
     *
     * Create a new `Sub_sparse_matrix` from a `Sparse_matrix` object (with the same `StorageFormat`). Create a "full" mask.
     *
     * \param[in] otherToCopy The matrix copied.
     */
    Sub_sparse_matrix(const Sparse_matrix<CoefficientRing,StorageFormat>& otherToCopy) : Sparse_matrix<CoefficientRing, StorageFormat>(otherToCopy)
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

    /** \brief Iterator to the end of the of the indices of non empty chains inside the mask. */
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
     * \param[in] index Index to turn on in the mask.
     */
    inline void set_bit_on (size_t index)
    {
        _subChains.setOn(index) ;
        if (this->_chainsStates.isOn(index))
            _subChainsStates.setOn(index) ;
    }

    /** \brief Removes an index from the mask.
     *
     * Set the bit encoding a given index to 0 (ie.\ remove the index from the mask).
     *
     * \param[in] index Index to turn off in the mask.
     */
    inline void set_bit_off (size_t index)
    {
        _subChains.setOff(index) ;
        _subChainsStates.setOff(index) ;
    }

    /** \brief Changes the mask to its complement. */
    inline void complement() { _subChains.bit_not() ; }

    /**
     * \brief Assigns to other `Sub_sparse_matrix`.
     *
     * Assign to other matrix coefficient-wise, and copy the bitboard.
     *
     * \pre The matrices must have the same type.
     *
     * \param[in] otherToCopy The matrix we want to copy.
     *
     * \return The reference to the modified matrix.
     */
    inline Sub_sparse_matrix& operator=(const Sub_sparse_matrix &otherToCopy)
    {
        (dynamic_cast<Sparse_matrix<CoefficientRing,StorageFormat>&>(*this)).operator=(otherToCopy) ;
        _subChains = otherToCopy._subChains ;
        _subChainsStates = otherToCopy._subChainsStates ;
        return *this ;
    }

    /**
     * \brief Displays a `Sub_sparse_matrix` in the output stream.
     *
     * Displays the sparse matrix as well as its mask.
     *
     * \param[in] stream The output stream.
     * \param[in] matrix The matrix to display.
     *
     * \return A reference to the modified stream.
     */
    friend std::ostream& operator<<(std::ostream &stream, const Sub_sparse_matrix &matrix) {
        stream << static_cast<const Sparse_matrix<CoefficientRing, StorageFormat>&>(matrix) ;
        stream << matrix._subChains << std::endl ;
        return stream ;
    }
};

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SUB_SPARSE_MATRIX_H
