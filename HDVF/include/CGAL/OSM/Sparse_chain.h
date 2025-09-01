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

#ifndef CGAL_OSM_SPARSE_CHAIN_H
#define CGAL_OSM_SPARSE_CHAIN_H

#include <CGAL/license/HDVF.h>

#include <CGAL/OSM/__base.h>
#include <CGAL/OSM/Sparse_matrix.h>
#include <unordered_map>
#include <vector>
#include <iterator>
#include <iostream>

namespace CGAL {
namespace OSM {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Sparse_chain` implements the concept `SparseChain`, that is, sparse vectors (encoding homological chains) optimized for topological computations.
 Given a complex (for instance a simplicial, cubical or cellular complex) which cells in a given dimension \f$q\f$ are numbered
 \f$\{\sigma_i\,;\, i=1\ldots n_q\}\f$, homology with coefficients in a given ring \f$\mathbb K\f$ defines \f$q\f$-chains as formal linear combinations:
 \f\[\gamma = \sum_{i=1}^{n_q}\lambda_i\cdot \sigma_i\f\]
 with coefficients \f$\lambda_i\in\mathbb K\f$.
 In the basis \f$\{\sigma_i\,;\, i=1\ldots n_q\}\f$, the coordinates of \f$\gamma\f$ are given by the vector: \f$[\lambda_1,\ldots \lambda_{n_q}]\f$.

 Now, as chains considered in homology are boundaries of cells, most \f$\lambda_i\f$ coefficients are null.
 Hence `Sparse_chain` encodes such vectors in a "sparse" way, that is, storing only non zero coefficients through a map:
 \f\[i\mapsto \lambda_i\ \ \ \forall i=1\ldots n_q\text{ such that }\lambda_i\neq 0\f\]
 Moreover, as per any linear algebra vector, a `Sparse_chain` is either a column or row vector (the `ChainTypeFlag` parameter determines the type).

 The class `Sparse_chain` provides standard linear algebra operators and fast iterators and block operations (set, get and nullify) which are required to implement efficiently HDVFs.

 \cgalModels{SparseChain}

 \tparam CoefficientRing a model of the `Ring` concept, providing the ring used to compute homology.
 \tparam ChainTypeFlag an integer constant encoding the type of matrices (`OSM::COLUMN` or `OSM::ROW`).
*/

template <typename CoefficientRing, int ChainTypeFlag>
class Sparse_chain {

public:
    /*!
     Type of chains iterators.
     */
    typedef typename std::unordered_map<size_t, CoefficientRing>::iterator iterator;

    /*!
     Type of chains constant iterators.
     */
    typedef typename std::unordered_map<size_t, CoefficientRing>::const_iterator const_iterator;

    // Allow the Sparse_matrix class to access other templated Sparse_matrix and
    // Sparse_chain protected members.
    template <typename _CT, int _CTF>
    friend class Sparse_chain;

    template <typename _CT, int _CTF>
    friend class Sparse_matrix;

protected:
    /* \brief Type of data stored in the chain: map between indexes and coefficients. */
    typedef std::pair<size_t, CoefficientRing> pair;

    /* \brief The chain inner representation and storage of data. */
    std::unordered_map<size_t, CoefficientRing> _chainData;

    /* \brief The chain boundary. */
    size_t _upperBound;

public:
    /**
     * \brief Creates new empty Sparse_chain.
     *
     * Default constructor, initialize an empty Sparse_chain.
     * The default chain size is 0.
     */
    Sparse_chain() {
        _upperBound = 0;
        _chainData = std::unordered_map<size_t, CoefficientRing>();
    }

    /**
     * \brief Creates new empty SparseChain of given size.
     *
     * Constructor with size, initialize an empty Sparse_chain.
     *
     * \param[in] chainSize The size of the Sparse_chain.
     */
    Sparse_chain(const size_t chainSize) {
        _upperBound = chainSize;
        _chainData = std::unordered_map<size_t, CoefficientRing>();
    }



    /**
     * \brief Creates new SparseChain by copy.
     *
     * Copy constructor, initialize a SparseChain from an existing SparseChain.
     *
     * \pre The chains have the same `CoefficientRing` and `ChainTypeFlag`.

     * \param[in] otherToCopy The chain to copy.
     */
    Sparse_chain(const Sparse_chain &otherToCopy) {
        _upperBound = otherToCopy._upperBound;
        _chainData = otherToCopy._chainData;
    }

    /**
     * \brief Assigns to other chain.
     *
     * Assign to other chain coefficient-wise, equivalent to copying it.
     *
     * \pre The chains have the same coefficent type.
     *
     * \warning Chains must have the same `CoefficientRing`.
     *
     * \param[in] otherToCopy The chain we want to copy.
     */
    Sparse_chain& operator=(const Sparse_chain &otherToCopy) {
        _upperBound = otherToCopy._upperBound;
        _chainData = otherToCopy._chainData;

        return *this;
    }

    /**
     * \brief Size of the chain
     *
     * \return Size allocated for the chain.
     */
    size_t dimension() const { return _upperBound ; }

    /** \relates Sparse_chain
     *
     * \brief Displays a Sparse_chain in the output stream.
     *
     * \param[in] stream The output stream.
     * \param[in] chain The chain to display.
     *
     * \return A reference to the modified stream.
     */
    friend std::ostream& operator<<(std::ostream &stream, const Sparse_chain &chain) {
        stream << "[";
        for (const_iterator i = chain._chainData.begin() ; i != chain._chainData.end() ; ++i) {
            stream << i->first << ": " << i->second << ", ";
        }

        if (chain._chainData.size() > 0) {
            stream << "\b\b";
        }
        stream << "]";

        return stream;
    }

    /**
     * \brief Adds two chains together.
     *
     * Adds two chains together and return the result in a new matrix.
     *
     * \pre Chains must have the same `CoefficientRing` and the same `ChainTypeFlag`.
     *
     * \warning Will raise an error if the two chains are not the same `CoefficientRing`.
     * \warning Will raise a compilation error if the two chains don't have the same `ChainTypeFlag`.
     *
     * \param[in] other The other chain.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator+(const Sparse_chain &other) {
        Sparse_chain newChain = *this;
        newChain += other;

        return newChain;
    }

    /**
     * \brief Subtracts a chain from current chain.
     *
     * Subtract `other` chain from current chain and return the result in a new matrix.
     *
     * \pre Chains must have the same `CoefficientRing` and the same `ChainTypeFlag`.
     *
     * \warning Will raise an error if the two chains are not the same `CoefficientRing`.
     * \warning Will raise a compilation error if the two chains don't have the same `ChainTypeFlag`.
     *
     * \param[in] other The other chain.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator-(const Sparse_chain &other) {
        Sparse_chain newChain = *this;
        newChain -= other;

        return newChain;
    }

    /*! \relates Sparse_chain
     *
     * \brief Applies multiplication on each coefficient.
     *
     * \param[in] lambda The factor to apply.
     * \param[in] chain The  chain.
     *
     * \return A new chain representing the result.
     */
    template <int _CTF>
    friend Sparse_chain operator*(const CoefficientRing& lambda, const Sparse_chain<CoefficientRing, _CTF> &chain) {
        Sparse_chain newChain = chain;
        newChain *= lambda;

        return newChain;
    }

    /**
     * \brief Applies multiplication on each coefficient.
     *
     * \param[in] lambda The factor to apply.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator*(const CoefficientRing& lambda) {
        Sparse_chain newChain = *this;
        newChain *= lambda;

        return newChain;
    }

    /** \relates Sparse_chain
     *
     * \brief Performs matrix multiplication between two chains (COLUMN x ROW) and return a COLUMN matrix.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     *
     * \pre Chains must have the same `CoefficientRing`.
     *
     * \warning Will raise an error if chains do not have the same `CoefficientRing`.
     *
     * \param[in] column The column chain.
     * \param[in] row The row chain.
     *
     * \return The result of the matrix multiplication, column-based.
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_chain<_CT, COLUMN>& column, const Sparse_chain<_CT, ROW>& row);

    /** \relates Sparse_chain
     *
     * \brief Performs matrix multiplication between two chains (COLUMN x ROW) and return a ROW matrix.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre Chains must have the same `CoefficientRing`.
     *
     * \warning Will raise an error if chains do not have the same `CoefficientRing`.
     *
     * \param[in] column The column chain.
     * \param[in] row The row chain.
     *
     * \return The result of the matrix multiplication, row-based.
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_chain<_CT, COLUMN> &column, const Sparse_chain<_CT, ROW> &row);

    /** \relates Sparse_chain
     *
     * \brief Performs dot product between two chains (ROW x COLUMN).
     *
     * \pre Chains must have the same `CoefficientRing`.
     *
     * \warning Will raise an error if the chains do not have the same `CoefficientRing`.
     *
     * \param[in] row The row chain.
     * \param[in] column The column chain.
     *
     * \return The result of type CoefficientRing.
     */
    template <typename _CT>
    friend _CT operator*(const Sparse_chain<_CT, ROW> &row, const Sparse_chain<_CT, COLUMN> &column);

    /**
     * \brief Adds a chain to `this`.
     *
     * Add a chain to `this`.
     *
     * \pre Chains must have the same `CoefficientRing` and the same `ChainTypeFlag`.
     *
     * \warning Will raise an error if the two chains are not the same `CoefficientRing`.
     * \warning Will raise an error if the two chains don't have the same `ChainTypeFlag`.
     *
     * \param[in] other The other chain.
     *
     * \return The modified chain representing the result.
     */
    Sparse_chain& operator+=(const Sparse_chain &other) {
        if (this->_upperBound != other._upperBound) {
            throw std::runtime_error("Chains must be the same size.");
        }

        for (pair pair: other._chainData) {
            this->_chainData[pair.first] += pair.second;

            if (this->_chainData[pair.first] == 0) {
                this->_chainData.erase(pair.first);
            }
        }

        return *this;
    }

    /**
     * \brief Sustracts a chain to `this`.
     *
     * Subtract a chain to `this`.
     *
     * \pre Chains must have the same `CoefficientRing` and the same `ChainTypeFlag`.
     *
     * \warning Will raise an error if the two chains are not the same `CoefficientRing`.
     * \warning Will raise an error if the two chains don't have the same `ChainTypeFlag`.
     *
     * \param[in] other The other chain.
     *
     * \return The modified chain representing the result.
     */
    Sparse_chain& operator-=(const Sparse_chain &other) {
        if (this->_upperBound != other._upperBound) {
            throw std::runtime_error("Chains must be the same size.");
        }

        for (pair pair: other._chainData) {
            this->_chainData[pair.first] -= pair.second;

            if (this->_chainData[pair.first] == 0) {
                this->_chainData.erase(pair.first);
            }
        }

        return *this;
    }

    /**
     * \brief Applies multiplication on each coefficient of `this`.
     *
     * If `lambda` is null, this function comes to nullify the chain.
     *
     * \param[in] lambda The factor to apply.
     *
     * \return The modified chain representing the result.
     */
    Sparse_chain& operator*=(const CoefficientRing& lambda) {
        if (lambda == 0) {
            this->_chainData.clear();
            return *this;
        }

        for (pair pair: this->_chainData) {
            this->_chainData[pair.first] = pair.second * lambda;
        }

        return *this;
    }

    /**
     * \brief Gets the value of a coefficient of the chain.
     *
     * \warning The chain will perform boundary check.
     *
     * \param[in] index The coefficient index.
     *
     * \return The value of the coefficient.
     */
    CoefficientRing operator[](size_t index) const {
        if (index >= _upperBound) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_upperBound) + ".");
        }

        if (_chainData.find(index) == _chainData.end())
            return 0 ;
        else
            return _chainData.at(index);
    }

    /**
     * \brief Gets the value of a coefficient of the chain.
     *
     * \warning The chain will perform boundary check.
     *
     * \param[in] index The coefficient index.
     *
     * \return The value of the coefficient.
     */
    CoefficientRing get_coefficient(size_t index) const {
        if (index >= _upperBound) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_upperBound) + ".");
        }

        if (_chainData.find(index) == _chainData.end())
            return 0 ;
        else
            return _chainData.at(index);
    }

    /**
     * \brief Sets a given coefficient of the chain.
     *
     * Set the value of the coefficient in the chain at `index`.
     *
     * \warning The chain will perform boundary check.
     *
     * \param[in] index The coefficient index.
     * \param[in] d Value of the coefficient
     */
    void set_coefficient(size_t index, CoefficientRing d)
    {
        if (index >= _upperBound) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_upperBound) + ".");
        }
        if (d == 0)
            *this /= index ;
        else
            _chainData[index] = d ;
    }

    /**
     * \brief Checks if a coefficient is null.
     *
     * \param[in] index The index to check.
     * \return True if the data is null at given index.
     */
    const bool is_null(size_t index) const {
        return _chainData.find(index) == _chainData.end();
    }

    /**
     * \brief Checks if the chain is null.
     *
     * \return True if the chain is null.
     */
    const bool is_null() const {
        return _chainData.size() == 0;
    }

    /** \relates Sparse_chain
     *
     * \defgroup ChainChainComparison Compares two chains.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Compare two chains and return `true` if both chains equal (and `false` otherwise).
     *
     * @param chain The first chain.
     * @param other The second chain.
     * @return A boolean.
     * @{
     */

    /** \relates Sparse_chain
     *
     * \brief Comparison of two COLUMN chains.
     */
    template <typename _CT>
    friend bool operator==(const Sparse_chain<_CT, OSM::COLUMN>& chain, const Sparse_chain<_CT, OSM::COLUMN> &other);

    /** \relates Sparse_chain
     *
     * \brief Comparison of a COLUMN  and a ROW chain.
     */
    template <typename _CT>
    friend bool operator==(const Sparse_chain<_CT, OSM::COLUMN>& chain, const Sparse_chain<_CT, OSM::ROW> &other);

    /** \relates Sparse_chain
     *
     * \brief Comparison of a ROW and a COLUMN chain.
     */
    template <typename _CT>
    friend bool operator==(const Sparse_chain<_CT, OSM::ROW>& chain, const Sparse_chain<_CT, OSM::COLUMN> &other);

    /** \relates Sparse_chain
     *
     * \brief Comparison of two ROW chains.
     */
    template <typename _CT>
    friend bool operator==(const Sparse_chain<_CT, OSM::ROW>& chain, const Sparse_chain<_CT, OSM::ROW> &other);

    /** @} */

    /**
     * \brief Gets a subchain from the chain.
     *
     * Return a new chain where all coefficients of indices provided in the vector are removed.
     *
     * \note Will return a copy of the chain if `indices` is empty.
     *
     * \param[in] indices The indexes to remove.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator/(const std::vector<size_t> &indices) {
        Sparse_chain newChain = *this;
        newChain /= indices;
        return newChain;
    }

    /**
     * \brief Gets a subchain from the chain.
     *
     * Return a new chain where the coefficients at a given index is removed.
     *
     * \param[in] index The index to remove.
     */
    Sparse_chain operator/(size_t index) {
        Sparse_chain newChain = *this;
        newChain /= index;
        return newChain;
    }

    /**
     * \brief Restricts the chain to a sub-chain by removing indices.
     *
     * Removes all indices provided in the vector from the chain. Return a reference to the modified chain.
     *
     * \note Will not alter the chain if `indices` is empty.
     *
     * \param[in] indices The indices to remove.
     *
     * \return Return a reference to the modified chain.
     */
    Sparse_chain& operator/=(const std::vector<size_t> &indices) {
        for (size_t index : indices) {
            this->_chainData.erase(index);
        }

        return *this;
    }

    /**
     * \brief Restricts the chain to a sub-chain by removing a given index.
     *
     * Removes the index provided from the chain.
     *
     * \note Will not alter the chain if given vector is empty.
     *
     * \param[in] index The index to remove.
     *
     * \return Return a reference to the modified chain.
     */
    Sparse_chain& operator/=(const size_t index) {
        this->_chainData.erase(index);

        return *this;
    }

    /**
     * \brief Removes all coefficients from the chain.
     *
     * The function comes to set all coefficients to zero.
     */
    void nullify() {
        this->_chainData.clear();
    }

    /**
     * \brief Iterator to the beginning of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns an iterator to the first non zero index.
     */
    iterator begin() noexcept {
        return _chainData.begin();
    }

    /**
     * \brief Constant iterator to the beginning of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the first non zero index.
     */
    const_iterator begin() const noexcept {
        return _chainData.begin();
    }

    /**
     * \brief Constant iterator to the beginning of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the first non zero index.
     */
    const_iterator cbegin() const noexcept {
        return _chainData.cbegin();
    }

    /**
     * \brief Iterator to the end of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns an iterator to the ending of the chain.
     */
    iterator end() noexcept {
        return _chainData.end();
    }

    /**
     * \brief Constant iterator to the end of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the ending of the chain.
     */
    const_iterator end() const noexcept {
        return _chainData.end();
    }

    /**
     * \brief Constant iterator to the end of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the ending of the chain.
     */
    const_iterator cend() const noexcept {
        return _chainData.cend();
    }

    /**
     * \brief Transposes a Sparse_chain.
     *
     * The result is a chain with `ChainTypeFlag` switched between COLUMN and ROW.
     *
     * \return A new chain where the chain type flag is changed.
     */
    Sparse_chain<CoefficientRing, COLUMN + ROW - ChainTypeFlag> transpose() {
        Sparse_chain<CoefficientRing, COLUMN + ROW - ChainTypeFlag> chain;

        chain._upperBound = this->_upperBound;
        chain._chainData = this->_chainData;

        return chain;
    }

    /**
     * \brief Checks if chain is a column.
     *
     * \return true if chain is column-major, false otherwise.
     */
    bool is_column() const {
        return ChainTypeFlag == COLUMN;
    }

    /**
     * \brief Checks if chain is a row.
     *
     * \return true if chain is row-major, false otherwise.
     */
    bool is_row() const {
        return ChainTypeFlag == ROW;
    }

private:
    /*
     * \brief Get a reference on a coefficient of the chain.
     *
     * Used to set a coefficient in the chain.
     *
     * \warning The chain will perform boundary check.
     *
     * \param[in] index The coefficient index.
     *
     * \return The reference to the assigned coefficient.
     */
    CoefficientRing& operator[](const size_t index)  {
        if (index >= _upperBound) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_upperBound) + ".");
        }

        return _chainData[index];
    }
};


// COLUMN chain x ROW chain -> COLUMN matrix
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_chain<_CT, COLUMN> &column, const Sparse_chain<_CT, ROW> &row) {
    Sparse_matrix<_CT, COLUMN> matrix(column._upperBound, row._upperBound);

    for (std::pair<size_t, _CT> pair : row._chainData) {
        OSM::set_column(matrix,pair.first,column * pair.second) ;
    }

    return matrix;
}

// COLUMN chain x ROW chain -> ROW matrix
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_chain<_CT, COLUMN> &column, const Sparse_chain<_CT, ROW> &row) {
    Sparse_matrix<_CT, ROW> matrix(column._upperBound, row._upperBound);

    for (std::pair<size_t, _CT> pair : column._chainData) {
        OSM::set_row(matrix,pair.first,row * pair.second);
    }

    return matrix;
}

// Dot product (ROW chain x COLUMN chain)
template <typename CoefficientRing>
CoefficientRing operator*(const Sparse_chain<CoefficientRing, ROW> &row, const Sparse_chain<CoefficientRing, COLUMN> &column) {
    // Get indexes (avoid adding double indexes).
    std::unordered_map<size_t, int> indexes;
    for (std::pair<size_t, CoefficientRing> pair: row._chainData) {
        indexes[pair.first] = 1;
    }
    for (std::pair<size_t, CoefficientRing> pair: column._chainData) {
        indexes[pair.first] += 1;
    }

    // Perform dot product
    CoefficientRing result = CoefficientRing();
    for (std::pair<size_t, int> index: indexes) {
        if (index.second == 2) {
            result += row._chainData.at(index.first) * column._chainData.at(index.first);
        }
    }

    return result;
}

//// Get a subchain from the chain and assign.
//template <typename _CT, int _CTF>
//Sparse_chain<_CT, _CTF> operator/(const Sparse_chain<_CT, _CTF> &chain, const std::vector<size_t> &indices) {
//    Sparse_chain<_CT, _CTF> newChain = chain;
//    newChain /= indices;
//    return newChain;
//}
//
//// Get a subchain from the chain and assign.
//template <typename _CT, int _CTF>
//Sparse_chain<_CT, _CTF> operator/(const Sparse_chain<_CT, _CTF> &chain, size_t index) {
//    Sparse_chain<_CT, _CTF> newChain = chain;
//    newChain /= index;
//    return newChain;
//}

template <typename _CT>
bool operator==(const Sparse_chain<_CT, OSM::COLUMN>& chain, const Sparse_chain<_CT, OSM::COLUMN> &other)
{
    typedef Sparse_chain<_CT, OSM::COLUMN> ChainType;
    bool res = true ;
    // Check that chains have the same size
    res = res && (chain._upperBound == other._upperBound) ;
    // Check that each coefficient of chain also belongs to other
    for (typename ChainType::const_iterator it = chain.begin(); res && (it != chain.end()); ++it)
    {
        res = res && (it->second == other.get_coefficient(it->first)) ;
    }
    // Check that each coefficient of other also belongs to chain
    for (typename ChainType::const_iterator it = other.begin(); res && (it != other.end()); ++it)
    {
        res = res && (it->second == chain.get_coefficient(it->first)) ;
    }
    return res ;
}

template <typename _CT>
bool operator==(const Sparse_chain<_CT, OSM::ROW>& chain, const Sparse_chain<_CT, OSM::ROW> &other)
{
    typedef Sparse_chain<_CT, OSM::ROW> ChainType;
    bool res = true ;
    // Check that chains have the same size
    res = res && (chain._upperBound == other._upperBound) ;
    // Check that each coefficient of chain also belongs to other
    for (typename ChainType::const_iterator it = chain.begin(); (it != chain.end()) && res; ++it)
    {
        res = res && (it->second == other.get_coefficient(it->first)) ;
    }
    // Check that each coefficient of other also belongs to chain
    for (typename ChainType::const_iterator it = other.begin(); res && (it != other.end()); ++it)
    {
        res = res && (it->second == chain.get_coefficient(it->first)) ;
    }
    return res ;
}

template <typename _CT>
bool operator==(const Sparse_chain<_CT, OSM::COLUMN>& chain, const Sparse_chain<_CT, OSM::ROW> &other)
{
    return false;
}

template <typename _CT>
bool operator==(const Sparse_chain<_CT, OSM::ROW>& chain, const Sparse_chain<_CT, OSM::COLUMN> &other)
{
    return false;
}

template <typename _CT, int _CTF>
Sparse_chain<_CT, _CTF> operator*(const Sparse_chain<_CT, _CTF> &chain, const _CT& lambda) {
    Sparse_chain newChain = chain;
    newChain *= lambda;

    return newChain;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SPARSE_CHAIN_H
