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

#ifndef CGAL_OSM_SPARSE_MATRIX_H
#define CGAL_OSM_SPARSE_MATRIX_H

#include <CGAL/license/HDVF.h>

#include <stdint.h>
#include <cmath>
#include <unordered_set>
#include <iostream>
#include <fstream>

#include <CGAL/OSM/Sparse_chain.h>
#include <CGAL/OSM/Bitboard.h>

// DEBUG : matrix output for SparseMatrices / no DEBUG : chain output for SparseMatrices
//#define DEBUG

namespace CGAL {
namespace OSM {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Sparse_matrix` implements the concept `SparseMatrix`, that is, sparse matrices optimized for topological computations. It provides standard linear algebra operators and fast iterators and block operations (set, get and nullify) which are required to implement efficiently HDVFs.

 The implementation is based on mapped sparse matrices. Hence matrices of the `Sparse_matrix` class are either column of row major (the `ChainTypeFlag` parameter determines the type). A column-major (resp. row-major) `Sparse_matrix` is a vector of `Sparse_chains` which encode columns (res. rows). Moreover, in order to efficiently iterate over non empty columns (resp. rows) the `Bitboard` data structure implements the concept `SparseMatrix::NonZeroChainIndices`. A bitboard is basically a bucket of bits recording the indices of non empty chains. However, this data structure has been designed in order to efficiently remove or add indices, as well as  provide efficient iterators to visit non empty chains.

 For instance, let us consider  the \f$5\times 4\f$ matrix:
 \f[
 A = \left(\begin{array}{cccc}
 1 & \cdot & \cdot & \cdot \\
 -1 & \cdot & 2 & \cdot\\
 \cdot & \cdot & 1 & \cdot \\
 \cdot & \cdot & \cdot & \cdot \\
 \cdot & \cdot & \cdot & \cdot \\
 \end{array}\right)
 \f]
 where \f$\cdot\f$ means \f$0\f$.

 Figures below shows the data structures are created according to the chosen representation  (left: column-major, right: row-major):

 <img src="Sparse_matrix_example_col.pdf" align="center" width=20%/>
 <img src="Sparse_matrix_example_row.pdf" align="center" width=20%/>

 \cgalModels{SparseMatrix}

 \tparam CoefficientRing a model of the `Ring` concept, providing the ring used to compute homology.
 \tparam ChainTypeFlag an integer constant encoding the type of matrices (`OSM::COLUMN` or `OSM::ROW`).
*/

template <typename CoefficientRing, int ChainTypeFlag>
class Sparse_matrix {

public:

    /*!
     Type of chains associated to the matrix.
     */
    typedef Sparse_chain<CoefficientRing, ChainTypeFlag> Matrix_chain;

    // Allow the Sparse_matrix class to access other templated Sparse_matrix private members.
    template <typename _CT, int _CTF>
    friend class Sparse_matrix;

protected:
    /* \brief The inner chain storage. */
    std::vector<Sparse_chain<CoefficientRing, ChainTypeFlag>> _chains;

    /* \brief A bitboard containing state of each columns. */
    Bitboard _chainsStates;

    /* \brief The matrix size as a (row, column) pair. */
    std::pair<size_t, size_t> _size;

    /*
     * \brief Get a reference on a chain from a matrix and change its state (for assignment).
     *
     * Used only internally.
     *
     * \warning The operator changes the status of the chain and the matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The reference to the chain stored at given index.
     */
    Matrix_chain& operator[](const size_t _index) {
        if (ChainTypeFlag == COLUMN && _index >= _size.second) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.second) + ".");
        }
        if (ChainTypeFlag == ROW && _index >= _size.first) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.first) + ".");
        }

        _chainsStates |= _index;
        return _chains[_index];
    }

public:
    /**
     * \brief Default constructor (empty new `Sparse_matrix` object).
     *
     * Create an empty matrix of type `ChainTypeFlag` with coefficients of type `CoefficientRing`.
     * The default matrix size is 0x0.
     */
    Sparse_matrix() {
        _chains = std::vector<Matrix_chain>(0);
        _chainsStates = Bitboard(0);
        _size = {0, 0};
    }

    /**
     * \brief Constructor with given rows/columns sizes.
     *
     * Create a new empty Sparse_matrix object of type `ChainTypeFlag` with coefficients of type `CoefficientRing` and a given size along rows/columns.
     *
     * \param[in] rowCount The number of rows to preallocate.
     * \param[in] columnCount The number of columns to preallocate.
     */
    Sparse_matrix(const size_t rowCount, const size_t columnCount) {
        size_t mainSize = ChainTypeFlag == COLUMN ? columnCount : rowCount;
        size_t secondarySize = ChainTypeFlag == COLUMN ? rowCount : columnCount;

        _chains = std::vector<Sparse_chain<CoefficientRing, ChainTypeFlag>>(mainSize);
        for (size_t i = 0 ; i < mainSize ; i++) {
            _chains[i] = Sparse_chain<CoefficientRing, ChainTypeFlag>(secondarySize);
        }

        _chainsStates = Bitboard(mainSize);
        _size = {rowCount, columnCount};
    }

    /**
     * \brief Copy constructor.
     *
     * Create a new SparseMatrix from another SparseMatrix object (with possibly a different `ChainTypeFlag`). Initialize a SparseMatrix of same sizes, containing the same coefficients (but not necessarly of the same `ChainTypeFlag`).
     * If types are different, the constructor performs conversion.
     *
     * \param[in] otherToCopy The matrix copied.
     */
    template <int CTF>
    Sparse_matrix(const Sparse_matrix<CoefficientRing,CTF> &otherToCopy) {
        if (ChainTypeFlag == CTF)
        {
            _chainsStates = otherToCopy._chainsStates;
            _size = otherToCopy._size;
            // Copy of _chains as such
            _chains.resize(otherToCopy._chains.size()) ;
            for (size_t i = 0; i<otherToCopy._chains.size(); ++i)
            {
                const Sparse_chain<CoefficientRing, CTF>& tmp(otherToCopy._chains.at(i)) ;
                Sparse_chain<CoefficientRing,ChainTypeFlag> res(tmp.dimension()) ;
                for (typename Sparse_chain<CoefficientRing, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
                {
                    res[it->first] = it->second ;
                }
                _chains[i] = res ;
            }
        }
        else // Copy coefficient by coefficient
        {
            _size = otherToCopy._size;
            size_t vec_size = (ChainTypeFlag == OSM::COLUMN)?_size.second:_size.first ;
            size_t chain_size = (ChainTypeFlag == OSM::COLUMN)?_size.first:_size.second ;
            _chains.resize(vec_size) ;
            _chainsStates = OSM::Bitboard(vec_size) ;
            for (size_t i=0; i<vec_size; ++i)
            {
                _chains.at(i) = Sparse_chain<CoefficientRing,ChainTypeFlag>(chain_size) ;
            }

            for (size_t i = 0; i<otherToCopy._chains.size(); ++i)
            {
                const Sparse_chain<CoefficientRing, CTF>& tmp(otherToCopy._chains.at(i)) ;
                for (typename Sparse_chain<CoefficientRing, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
                {
                    _chains[it->first][i] = it->second ;
                    _chainsStates.setOn(it->first) ;
                }
            }
        }
    }

    /**
     * \brief Assigns to other matrix.
     *
     * Assign to other matrix coefficient-wise, equivalent to copying it.
     *
     * \pre The matrices must have the same type.
     *
     * \param[in] otherToCopy The matrix we want to copy.
     *
     * \return The reference to the modified matrix.
     */
    Sparse_matrix& operator=(const Sparse_matrix& otherToCopy)
    {
        _chainsStates = otherToCopy._chainsStates;
        _size = otherToCopy._size;
        _chains = otherToCopy._chains;

        return *this;
    }


    /**
     * \brief Cleans a SparseMatrix (set all coefficients to zero).
     *
     * Empty all structures of the sparse matrix.
     */
    void nullify() {
        std::vector<size_t> coefs ;
        for (OSM::Bitboard::iterator it = begin(); it != end(); ++it)
        {
            coefs.push_back(*it) ;
        }
        *this /= coefs ;
    }

    /**
     * \brief Tests if a SparseMatrix is null.
     *
     * The function return `true` is the SparseMatrix is null (that is, empty) and `false` otherwise.
    */
    bool is_null()
    {
        return (_chainsStates.begin() == _chainsStates.end()) ;
    }

    /**
     * \defgroup MatrixMatrixComparison Compares two matrices.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Compare two matrices and return `true` if both matrices equal (and `false` otherwise).
     *
     * @param matrix The first matrix.
     * @param other The second matrix.
     * @return A boolean.
     * @{
     */

    /** \brief Comparison of two COLUMN matrices. */
    template <typename _CT>
    friend bool operator==(const Sparse_matrix<_CT, OSM::COLUMN>& matrix, const Sparse_matrix<_CT, OSM::COLUMN> &other);

    /** \brief Comparison of a COLUMN  and a ROW matrix. */
    template <typename _CT>
    friend bool operator==(const Sparse_matrix<_CT, OSM::COLUMN>& matrix, const Sparse_matrix<_CT, OSM::ROW> &other);

    /** \brief Comparison of a ROW and a COLUMN matrix. */
    template <typename _CT>
    friend bool operator==(const Sparse_matrix<_CT, OSM::ROW>& matrix, const Sparse_matrix<_CT, OSM::COLUMN> &other);

    /** \brief Comparison of two ROW matrices. */
    template <typename _CT>
    friend bool operator==(const Sparse_matrix<_CT, OSM::ROW>& matrix, const Sparse_matrix<_CT, OSM::ROW> &other);

    /** @} */

    /**
     * \brief Displays a matrix in the output stream.
     *
     * \param[in] stream The output stream.
     * \param[in] matrix The matrix to display.
     *
     * \return A reference to the modified stream.
     */
    friend std::ostream& operator<<(std::ostream &stream, const Sparse_matrix &matrix) {
#ifndef DEBUG
        bool empty = true;

        stream << "[";
        for (size_t index : matrix._chainsStates) {
            empty = false;
            stream << index << ": " << matrix._chains[index] << ", ";
        }

        if (!empty) {
            stream << "\b\b";
        }
        stream << "]" << std::endl ;
#else
        for (size_t i=0; i<matrix._size.first; ++i)
        {
            for (size_t j=0; j<matrix._size.second; ++j)
            {
                CoefficientRing tmp = matrix.get_coefficient(i,j) ;

                if (tmp == 0)
                    stream << ".\t" ;
                else
                    stream << tmp << "\t" ;
            }
            stream << std::endl ;
        }
#endif

        return stream;
    }

    /** \relates Sparse_matrix
     *
     * \defgroup WriteMatrix Writes matrix to an output stream or a file.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Write a sparse matrix to an output stream or a file using the `.osm` file format.
     *
     * Output a sparse matrix  (the matrix can be reloaded using `read_matrix`).
     *
     * @{
     */

    /** \brief Writes a sparse COLUMN matrix to a stream. */
    template <typename _CT>
    friend std::ostream& write_matrix (const Sparse_matrix<_CT, OSM::COLUMN>& M, std::ostream& out);

    /** \brief Writes a sparse ROW matrix to a stream. */
    template <typename _CT>
    friend std::ostream& write_matrix (const Sparse_matrix<_CT, OSM::ROW>& M, std::ostream& out);

    /** \brief Writes a sparse COLUMN matrix to a file. */
    template <typename _CT>
    friend void write_matrix (const Sparse_matrix<_CT, OSM::COLUMN>& M, std::string filename);

    /** \brief Writes a sparse ROW matrix to a file. */
    template <typename _CT>
    friend void write_matrix (const Sparse_matrix<_CT, OSM::ROW>& M, std::string filename);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup ReadMatrix Reads matrix from an input stream.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Read a sparse matrix from an input stream or a file using the `.osm` file format.
     *
     * Read a sparse matrix (the input must respect the `.osm` file format).
     *
     * @{
     */

    /** \brief Reads a sparse COLUMN matrix from a stream. */
    template <typename _CT>
    friend std::istream& read_matrix (Sparse_matrix<_CT, OSM::COLUMN>& M, std::istream& in);

    /** \brief Reads a sparse ROW matrix from a stream. */
    template <typename _CT>
    friend std::istream& read_matrix (Sparse_matrix<_CT, OSM::ROW>& M, std::istream& in);

    /** \brief Reads a sparse COLUMN matrix from a file. */
    template <typename _CT>
    friend void read_matrix (Sparse_matrix<_CT, OSM::COLUMN>& M, std::string filename);

    /** \brief Reads a sparse ROW matrix from a file. */
    template <typename _CT>
    friend void read_matrix (Sparse_matrix<_CT, OSM::ROW>& M, std::string filename);

    /** @} */

    /**
     * \brief Adds two matrices together into a new matrix.
     *
     * Adds each coefficient of the matrices together and returns a new matrix (of the same type as `this`) representing the result (when possible, prefer `+=` for efficiency).
     *
     * \pre Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`.
     *
     * \warning Will raise an error if the other matrix is not the same `CoefficientRing`.
     *
     * \param[in] other The second matrix.
     *
     * \return A new matrix representing the result.
     */
    template <int _CTF>
    Sparse_matrix operator+(const Sparse_matrix<CoefficientRing, _CTF> &other) {
        Sparse_matrix newMatrix(*this);
        newMatrix += other;

        return newMatrix;
    }

    /**
     * \brief Substracts two matrices together into a new matrix.
     *
     * Substracts each coefficient of the matrix `other` and returns a new matrix (of the same type as `this`) representing the result (when possible, prefer `-=` for efficiency).
     *
     * \pre Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`.
     *
     * \warning Will raise an error if the other matrix is not the same `CoefficientRing`.
     *
     * \param[in] other The second matrix.
     *
     * \return A new matrix representing the result.
     */
    template <int _CTF>
    Sparse_matrix operator-(const Sparse_matrix<CoefficientRing, _CTF> &other) {
        Sparse_matrix newMatrix(*this);
        newMatrix -= other;

        return newMatrix;
    }

    /** \relates Sparse_matrix
     *
     * \brief Applies factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`.
     *
     * If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     *
     * \param[in] lambda The factor to apply.
     * \param[in] matrix The matrix.
     *
     * \return A new matrix representing the result.
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix operator*(const _CT& lambda, const Sparse_matrix<_CT, _CTF> &matrix) {
        Sparse_matrix newMatrix = matrix;
        newMatrix *= lambda;

        return newMatrix;
    }

    /**
     * \brief Applies factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`.
     *
     * If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     *
     * \param[in] lambda The factor to apply.
     *
     * \return A new matrix representing the result.
     */
    Sparse_matrix operator*(const CoefficientRing& lambda) {
        Sparse_matrix newMatrix = *this;
        newMatrix *= lambda;

        return newMatrix;
    }
    
    /** \relates Sparse_matrix
     *
     * \defgroup MatrixMatrixProdCol Matrices product (with column-based result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between matrices and returns a new column-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new column-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`. The product is optimized (using standard definition or block products) for each combination of `ChainTypeFlag`. However, efficiency depends on `ChainTypeFlag` (when possible, prefer row-major by column-major products).
     *
     * @param first The first matrix.
     * @param second The second matrix.
     * @return The result of the matrix multiplication, column-based.
     * @{
     */

    /** \brief Matrices product: COLUMN x COLUMN -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, COLUMN> &second);

    /** \brief Matrices product: ROW x COLUMN -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, COLUMN> &second);

    /** \brief Matrices product: COLUMN x ROW -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, ROW> &second);

    /** \brief Matrices product: ROW x ROW -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, ROW> &second);

    /** @} */

    /** \relates Sparse_matrix
     *
     *\defgroup MatrixChainProd Matrix / column-chain product (with column-chain result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between a sparse matrix (column or row major) and a column-chain. The function returns a new column-major chain.
     *
     * Perform standard linear algebra product between a matrix and a column-chain (ie.\ matrix / column vector product) and returns a new column-major chain. Both arguments must have the same `CoefficientRing` but the matrix can have any `ChainTypeFlag` (and the product is optimized for each of them).
     *
     * @param first The matrix.
     * @param second The column-major chain.
     * @return The result of the matrix multiplication, column-based.
     * @{
     */

    /** \brief Matrix/column chain product: COLUMN matrix x COLUMN chain -> COLUMN chain. */
    template <typename _CT>
    friend Sparse_chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_chain<_CT, COLUMN> &_second);

    /** \brief Matrix/column chain product: ROW matrix x COLUMN chain -> COLUMN chain. */
    template <typename _CT>
    friend Sparse_chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Sparse_chain<_CT, COLUMN> &_second);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup ChainMatrixProd Row-chain / matrix product (with row-chain result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between a row-chain and a sparse matrix (column or row major). The function returns a new row-major chain.
     *
     * Perform standard linear algebra product between a row-chain and a matrix (ie.\ row vector / matrix product) and returns a new row-major chain. Both arguments must have the same `CoefficientRing` but the matrix can have any `ChainTypeFlag` (and the product is optimized for each of them).
     *
     * @param first The row-major chain.
     * @param second The matrix.
     * @return The result of the matrix multiplication, row-based.
     * @{
     */

    /** \brief Row chain/matrix product: ROW chain x COLUMN matrix -> ROW chain. */
    template <typename _CT>
    friend Sparse_chain<_CT, ROW> operator*(const Sparse_chain<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second) ;

    /** \brief Row chain/matrix product: ROW chain x ROW matrix -> ROW chain. */
    template <typename _CT>
    friend Sparse_chain<_CT, ROW> operator*(const Sparse_chain<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second) ;

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup MatrixMatrixProdRow Matrices product (with row-based result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between matrices and returns a new row-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new row-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`. The product is optimized (using standard definition or block products) for each combination of `ChainTypeFlag`. However, efficiency depends on `ChainTypeFlag`.
     *
     * @param first The first matrix.
     * @param second The second matrix.
     * @return The result of the matrix multiplication, row-based.
     * @{
     */

    /** \brief Matrices product: COLUMN x COLUMN -> ROW */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, COLUMN> &_second);

    /** \brief Matrices product: ROW x COLUMN -> ROW */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second);

    /** \brief Matrices product: COLUMN x ROW -> ROW */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, ROW> &_second);

    /** \brief Matrices product: ROW x ROW -> ROW */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup MatrixMatrixAddAssign Sums matrices and assign.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform addition between matrices and assigns the result to `matrix`.
     *
     * Perform standard linear algebra addition. Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`. The addition is optimized for each combination of `ChainTypeFlag`.
     *
     * @param matrix The first matrix.
     * @param other The second matrix.
     * @return A reference to the modified matrix `matrix` containing the result.
     * @{
     */

    /** \brief Matrices sum and assign: COLUMN += COLUMN or ROW += ROW. */
    Sparse_matrix& operator+=(const Sparse_matrix &other) {
        if (this->_size != other._size) {
            throw std::runtime_error("Matrices must be the same _size.");
        }

        for (size_t index: other._chainsStates) {
            this->_chainsStates |= index;
            this->_chains[index] += other._chains[index];

            if (this->_chains[index].is_null()) {
                this->_chainsStates.setOff(index);
            }
        }

        return *this;
    }

    /** \brief Matrices sum and assign: COLUMN += ROW. */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator+=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, ROW> &other);

    /** \brief Matrices sum and assign: ROW += COLUMN. */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator+=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, COLUMN> &other);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup MatrixMatrixSubtractAssign Subtracts matrices and assign.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform subtraction between matrices and assigns the result to `matrix`.
     *
     * Perform standard linear algebra subtraction. Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`. The addition is optimized for each combination of `ChainTypeFlag`.
     *
     * @param matrix The first matrix.
     * @param other The second matrix.
     * @return A reference to the modified matrix `matrix` containing the result.
     * @{
     */

    /** \brief Matrices subtraction and assign: COLUMN -= COLUMN or ROW -= ROW. */
    Sparse_matrix& operator-=(const Sparse_matrix &other) {
        if (this->_size != other._size) {
            throw std::runtime_error("Matrices must be the same _size.");
        }

        for (size_t index: other._chainsStates) {
            this->_chainsStates |= index;
            this->_chains[index] -= other._chains[index];

            if (this->_chains[index].is_null()) {
                this->_chainsStates.setOff(index);
            }
        }

        return *this;
    }

    /** \brief Matrices subtraction and assign: COLUMN -= ROW. */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator-=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, ROW> &other);

    /** \brief Matrices subtraction and assign: ROW -= COLUMN. */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator-=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, COLUMN> &other);

    /** @} */

    /**
     * \brief Applies factor on each coefficients and assign.
     *
     * If `lambda` is 0, this comes to nullify the matrix.
     *
     * \param[in] lambda The factor to apply.
     *
     * \return The modified matrix representing the result.
     */
    Sparse_matrix& operator*=(const CoefficientRing& lambda) {
        if (lambda == 0) {
            this->nullify();
            return *this;
        }

        for (size_t index: this->_chainsStates) {
            this->_chains[index] *= lambda;
        }

        return *this;
    }

    /**
     * \brief Computes the negative of a matrix (unary operator).
     *
     * \return The resulting matrix.
     */
    Sparse_matrix operator-() {
        Sparse_matrix res(this->_size.first, this->_size.second) ;

        for (size_t index: this->_chainsStates)
        {
            const Matrix_chain& tmp_chain(this->_chains[index]) ;
            for (typename Matrix_chain::const_iterator it = tmp_chain.cbegin(); it != tmp_chain.cend(); ++it)
                res[index][it->first] = -it->second ;
        }
        return res ;
    }

    /** \relates Sparse_matrix
     *
     * \defgroup MatrixMatrixProdAssign Multiplies matrices and assign.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform product between matrices and assigns the result to `matrix`.
     *
     * Perform standard linear algebra product. Matrices must have the same `CoefficientRing` but can have different `ChainTypeFlag`. The product is optimized for each combination of `ChainTypeFlag`.
     *
     * @param matrix The first matrix.
     * @param other The second matrix.
     * @return A reference to the modified matrix `matrix` containing the result.
     * @{
     */

    /** \brief Matrices product and assign: COLUMN *= COLUMN. */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, COLUMN> &other);

    /** \brief Matrices product and assign: ROW *= ROW. */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, ROW> &other);

    /** \brief Matrices product and assign: COLUMN *= ROW. */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, ROW> &other);

    /** \brief Matrices product and assign: ROW *= COLUMN. */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, COLUMN> &other);

    /** @} */

    /**
     * \brief Gets the value of a chain from a const matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] index The chain index.
     *
     * \return The chain stored at given index.
     */
    Matrix_chain operator[](size_t index) const {
        if (ChainTypeFlag == COLUMN && index >= _size.second) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.second) + ".");
        }
        if (ChainTypeFlag == ROW && index >= _size.first) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.first) + ".");
        }

        return _chains[index];
    }

protected:
    // Protected method for set_coefficient
    void set_coefficient(const size_t i, const size_t j, const CoefficientRing d) {
        if (i >= _size.first) {
            throw std::runtime_error("Provided i index should be less than " + std::to_string(_size.first) + ".");
        }
        if (j >= _size.second) {
            throw std::runtime_error("Provided j index should be less than " + std::to_string(_size.second) + ".");
        }

        if (d != 0)
        {
            if (ChainTypeFlag == COLUMN)
            {
                (*this)[j][i] = d ;
            }
            else
                (*this)[i][j] = d ;
        }
        else
        {
            del_coefficient(i, j);
        }
    }
public:
    /** \relates Sparse_matrix
     *
     * \brief Sets a given coefficient in `matrix`.
     *
     * Assign the scalar `d` to the coefficient on row `i` and column `j`.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] matrix Reference on the matrix to modify.
     * \param[in] i The row index.
     * \param[in] j The column index.
     * \param[in] d The value.
     */
    template <typename _CT, int _CTF>
    friend void set_coefficient(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j, const _CT d);

protected:
    // Protected method for get_coefficient
    CoefficientRing get_coefficient(const size_t i, const size_t j) const {
        if (i >= _size.first) {
            throw std::runtime_error("Provided _i index should be less than " + std::to_string(_size.first) + ".");
        }
        if (j >= _size.second) {
            throw std::runtime_error("Provided _j index should be less than " + std::to_string(_size.second) + ".");
        }

        if (ChainTypeFlag == COLUMN)
            return (this->_chains)[j][i] ;
        else // ROW
            return (this->_chains)[i][j] ;
    }
public:
    /** \relates Sparse_matrix
     *
     * \brief Gets a given coefficient.
     *
     * Returns the coefficient on row `i` and column `j` of the matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] matrix Constant reference on the matrix.
     * \param[in] i The row index.
     * \param[in] j The column index.
     *
     * \return The value of the given coefficient.
     */
    template <typename _CT, int _CTF>
    friend _CT get_coefficient(const Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j);

    /** \relates Sparse_matrix
     *
     * \defgroup GetColumn Gets a column.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief Get the value of the column at a given `index` from the matrix (whatever the `ChainTypeFlag` of the matrix).
     *
     * \note For column-matrices, it is equivalent to `operator[]`, for row-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix considered.
     * @param index The coefficient index.
     * @return The column at given index.
     * @{
     */

    /** \brief Gets a column from a COLUMN matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index);

    /** \brief Gets a column from a ROW matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, ROW> &matrix,  size_t index);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup GetRow Gets a row.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief Get the value of the row at a given `index` from the matrix (whatever the `ChainTypeFlag` of the matrix).
     *
     * \note For row-matrices, it is equivalent to `operator[]`, for column-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix considered.
     * @param index The coefficient index.
     * @return The row at given index.
     * @{
     */

    /** \brief Gets a row from a COLUMN matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index);

    /** \brief Gets a row from a ROW matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, ROW> &matrix,  size_t index);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \brief Gets a const reference over a column from a column matrix.
     *
     * Constant time get.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] matrix The matrix considered.
     * \param[in] index The column index.
     *
     * \return A constant reference over the column stored at given index.
     */
    template <typename _CT>
    friend const Sparse_chain<_CT, COLUMN> & cget_column(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index);

    /** \relates Sparse_matrix
     *
     * \brief Gets a constant reference over a row from a row matrix
     *
     * Constant time get.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] matrix The matrix considered.
     * \param[in] index The row index.
     *
     * \return A const reference over the row stored at given index.
     */
    template <typename _CT>
    friend const Sparse_chain<_CT, ROW> & cget_row(const Sparse_matrix<_CT, ROW> &matrix, const size_t index);


    /** \relates Sparse_matrix
     *
     * \defgroup SetColumn Sets a column.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief Set the value of the column at a given `index` from the matrix to `chain` (whatever the `ChainTypeFlag` of the matrix).
     *
     * \note For column-matrices, it is equivalent to `operator[]` followed by an assignment, for row-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix.
     * @param index The column index.
     * @param chain The new column value.
     * @{
     */

    /** \brief Sets a column in a COLUMN matrix. */
    template <typename _CT>
    friend void set_column(Sparse_matrix<_CT, COLUMN> &matrix,  size_t index, const Sparse_chain<_CT, COLUMN> &chain);

    /** \brief Sets a column in a ROW matrix. */
    template <typename _CT>
    friend void set_column(Sparse_matrix<_CT, ROW> &matrix,  size_t index, const Sparse_chain<_CT, COLUMN> &chain);

    /** @} */

    /** \relates Sparse_matrix
     *
     * \defgroup SetRow Sets a row.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief Set the value of the row at a given `index` from the matrix to `chain` (whatever the `ChainTypeFlag` of the matrix).
     *
     * \note For row-matrices, it is equivalent to `operator[]` followed by an assignment, for column-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix.
     * @param index The row index.
     * @param chain The new row value.
     * @{
     */

    /** \brief Sets a row in a COLUMN matrix. */
    template <typename _CT>
    friend void set_row(Sparse_matrix<_CT, COLUMN> &matrix,  size_t index, const Sparse_chain<_CT, ROW> &chain);

    /** \brief Sets a row in a ROW matrix. */
    template <typename _CT>
    friend void set_row(Sparse_matrix<_CT, ROW> &matrix,  size_t index, const Sparse_chain<_CT, ROW> &chain);

    /** @} */

    /**
     * \defgroup GetBlockMatrix Gets a sub-matrix by removing chains in a matrix.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief Build and return a new matrix by copying a `matrix` and removing  all chains with indices in the `indices` vector (or at a given `index`). The result is thus a block of the initial matrix. When possible, for efficiency, prefer `\=`.
     *
     * \note Will return a copy of the matrix if given vector is empty.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The initial matrix.
     * @param indices/index The indice(s) of the chain(s) to remove.
     *
     * @return A new matrix containing the result.
     * @{
     */

    /** \brief Removes a set of chains from a copy of the matrix. */
    Sparse_matrix operator/(const std::vector<size_t> &_indices) {
        Sparse_matrix res(*this);
        res /= _indices;
        return res;
    }

    /** \brief Removes the chain at a given `index` from a copy of the matrix. */
    Sparse_matrix operator/(size_t index) {
        Sparse_matrix res(*this);
        res /= index;
        return res;
    }

    /** @} */

    /**
     * \defgroup GetBlockMatrixAssign Assign a sub-matrix by removing chains in a matrix.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief Remove  all chains with indices in the `indices` vector (or at a given `index`). The result is thus a block of the initial matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param indices/index The indice(s) of the chain(s) to remove.
     *
     * @return A reference on the modified matrix.
     * @{
     */

    /** \brief Removes a set of chains from a matrix. */
    Sparse_matrix& operator/=(const std::vector<size_t> &indexes) {
        for (size_t index : indexes) {
            *this /= index;
        }

        return *this;
    }

    /** \brief Removes the chain at a given `index` from a matrix. */
    Sparse_matrix& operator/=(const size_t index) {
        _chains[index].nullify();
        _chainsStates.setOff(index);

        return *this;
    }

    /** @} */

protected:
    // Protected version of del_column
    Sparse_matrix& del_column(size_t index) {
        std::vector<size_t> tmp_id{index} ;
        if (ChainTypeFlag == OSM::COLUMN)
        {
            (*this)/=tmp_id ;
        }
        else // OSM::ROW
        {
            for (size_t ind : _chainsStates)
            {
                Matrix_chain &tmp(_chains[ind]) ;
                tmp/=tmp_id ;
                // If the row index has become empty: update
                if (tmp.is_null())
                    *this /= std::vector<size_t>({ind}) ;
            }
        }

        return *this;
    }
public:
    /** \relates Sparse_matrix
     *
     * \brief Removes a column from the matrix.
     *
     * Removes column of index `index` whatever the `ChainTypeFlag` of the matrix. For column matrices, it just comes to the `\=` operator and for row matrices, it entails a traversal of the matrix.
     *
     * \param[in] matrix Reference on the matrix to modify.
     * \param[in] index The index to remove.
     *
     * \return The modified matrix representing the result.
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix<_CT, _CTF>& del_column(Sparse_matrix<_CT, _CTF>& matrix, size_t index);

protected:
    // Protected version of del_row
    Sparse_matrix& del_row(size_t index)
    {
        std::vector<size_t> tmp_id{index};
        if (ChainTypeFlag == OSM::ROW) {
            (*this) /= tmp_id;
        } else // OSM::COLUMN
        {
            for (size_t ind : _chainsStates) {
                Matrix_chain &tmp(_chains[ind]);
                tmp /= tmp_id;
                // If the column index has become empty: update
                if (tmp.is_null())
                    *this /= std::vector<size_t>({ind}) ;
            }
        }

        return *this;
    }

public:
    /** \relates Sparse_matrix
     *
     * \brief Removes a row from the matrix.
     *
     * Removes row of index `index` whatever the `ChainTypeFlag` of the matrix. For row matrices, it just comes to the `\=` operator and for column matrices, it entails a traversal of the matrix.
     *
     * \param[in] matrix Reference on the matrix to modify.
     * \param[in] index The index to remove.
     *
     * \return The modified matrix representing the result.
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix<_CT, _CTF>& del_row(Sparse_matrix<_CT, _CTF>& matrix, size_t index);

protected:
    // Protected version of del_coefficient
    Sparse_matrix& del_coefficient(size_t i, size_t j) {
        // OSM::COLUMN
        if (ChainTypeFlag == OSM::COLUMN) {
            std::vector<size_t> tmp_id({i}) ;
            Matrix_chain &tmp(_chains[j]);
            tmp /= tmp_id;
            if (tmp.is_null())
                _chainsStates.setOff(j) ;
        } else // OSM::ROW
        {
            std::vector<size_t> tmp_id({j}) ;
            Matrix_chain &tmp(_chains[i]);
            tmp /= tmp_id;
            if (tmp.is_null())
                _chainsStates.setOff(i) ;
        }
        return *this;
    }

public:
    /** \relates Sparse_matrix
     * 
     * \brief Removes a coefficient from the matrix.
     *
     * Removes coefficient at row `i` and column `j`.
     *
     * \param[in] matrix Reference on the matrix to modify.
     * \param[in] i Index of the row
     * \param[in] j Index of the column
     *
     * \return The modified matrix representing the result.
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix<_CT, _CTF>& del_coefficient(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j);

    /**
     * \brief Iterator to the index of the first non null chain.
     *
     * Return an iterator to the index of the first non null chain (the iterator visits indices of non null chains along the major dimension of the matrix).
     *
     * \return The iterator to the index of the first non null chain.
     */
    inline Bitboard::iterator begin() const noexcept { return _chainsStates.begin(); }

    /**
     * \brief Iterator to the ending of chains indices.
     *
     * \return The iterator to the ending of chains indices.
     */
    inline Bitboard::iterator end() const noexcept { return _chainsStates.end(); }

    /**
     * \brief Reverse iterator to the index of the last non null chain.
     *
     * Return a reverse iiterator to the index of the last non null chain (the iterator visits indices of non null chains, in decreading order, along the major dimension of the matrix).
     *
     * \return The reverse iterator to the index of the last non null chain.
     */
    inline Bitboard::reverse_iterator reverse_begin() noexcept { return _chainsStates.reverse_begin(); }
    inline Bitboard::reverse_iterator reverse_begin(size_t index) noexcept { return _chainsStates.reverse_begin(index); }

    /**
     * \brief Reverse iterator to the ending of chains indices.
     *
     * \return The reverse iterator to the ending of chains indices.
     */
    inline Bitboard::reverse_iterator reverse_end() noexcept { return _chainsStates.reverse_end(); }


    /**
     * \brief Transposes a matrix.
     *
     * \return A new matrix where the chain type flag has been swapped between COLUMN and ROW and data chains have been transposed.
     */
    Sparse_matrix<CoefficientRing, COLUMN + ROW - ChainTypeFlag> transpose() {
        Sparse_matrix<CoefficientRing, COLUMN + ROW - ChainTypeFlag> transposed(this->_size.second, this->_size.first);

        for (size_t index : this->_chainsStates) {
            transposed._chains[index] = this->_chains[index].transpose();
        }

        transposed._chainsStates = this->_chainsStates;

        return transposed;
    }

    /**
     * \brief Gets the matrix sizes.
     *
     * \return The matrix size as a row/column pair.
     */
    std::pair<size_t, size_t> dimensions() const {
        return this->_size;
    }
};

template <typename _CT, int _CTF>
Sparse_matrix<_CT, _CTF> operator*(const Sparse_matrix<_CT, _CTF> &matrix, const _CT& lambda){
    Sparse_matrix<_CT, _CTF> newMatrix = matrix;
     newMatrix *= lambda;

    return newMatrix;
}


// Matrix-matrix product
// COLUMN x COLUMN -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, COLUMN> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);

    // Perform col-col matrix multiplication with linear combination of columns.
    for (size_t index: second._chainsStates) {
        Sparse_chain<_CT, COLUMN> column(first._size.first);

        for (auto colRight: second._chains[index]) {
            if (first._chainsStates.isOn(colRight.first)) {
                column += colRight.second * first._chains[colRight.first];
            }
        }

        res[index] = column;
    }

    return res;
}

// Matrix-matrix product
// ROW x COLUMN -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, COLUMN> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);

    // Perform row-col matrix multiplication with dot products.
    for (size_t colRight: second._chainsStates) {
        Sparse_chain<_CT, COLUMN> column(first._size.first);

        for (size_t rowLeft: first._chainsStates) {
            _CT coef = first._chains[rowLeft] * second._chains[colRight];
            if (coef != 0) {
                column[rowLeft] = coef;
            }
        }

        res[colRight] = column;
    }

    return res;
}

// Matrix-matrix product
// COLUMN x ROW -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, ROW> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);

    // Perform row-col matrix multiplication with dot products.
    for (size_t colLeft: first._chainsStates) {
        res += first._chains[colLeft] * second._chains[colLeft];
    }

    return res;
}

// Matrix-matrix product
// ROW x ROW -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, ROW> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);

    // Perform row-col matrix multiplication with dot products.
    for (size_t i = 0 ; i < second._size.second ; i++) {
        Sparse_chain<_CT, COLUMN> column(first._size.first);

        for (size_t rowLeft: first._chainsStates) {
            _CT coef = first._chains[rowLeft] * get_column(second, i);
            if (coef != 0) {
                column[rowLeft] = coef;
            }
        }

        if (!column.is_null()) {
            res[i] = column;
        }
    }

    return res;
}

// Matrix - column chain product
// COLUMN matrix
template <typename _CT>
Sparse_chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_chain<_CT, COLUMN> &second)
{
    // Perform col-col matrix multiplication with linear combination of columns.
    Sparse_chain<_CT, COLUMN> column(first._size.first);

    for (typename Sparse_chain<_CT, COLUMN>::const_iterator it = second.begin(); it != second.end(); ++it)
    {
        column += it->second * first._chains[it->first];
    }

    return column;
}

// Matrix - column chain product
// ROW matrix
template <typename _CT>
Sparse_chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_chain<_CT, COLUMN> &second)
{
    // Perform row-col matrix multiplication with dots
    Sparse_chain<_CT, COLUMN> column(first._size.first);

    for (size_t index : first._chainsStates)
    {
        _CT tmp(first[index] * second) ;
        if (tmp != 0)
//            column[index] = tmp ;
            column.set_coefficient(index, tmp);
    }
    return column;
}

// Row chain - matrix product
// ROW matrix
template <typename _CT>
Sparse_chain<_CT, ROW> operator*(const Sparse_chain<_CT, ROW> &first, const Sparse_matrix<_CT, ROW> &second)
{
    // Perform row-row matrix multiplication with linear combination of rows.
    Sparse_chain<_CT, ROW> row(second._size.second);

    for (typename Sparse_chain<_CT, ROW>::const_iterator it = first.begin(); it != first.end(); ++it)
    {
        row += it->second * second._chains[it->first];
    }

    return row;
}

// Row chain - matrix product
// COLUMN matrix
template <typename _CT>
Sparse_chain<_CT, ROW> operator*(const Sparse_chain<_CT, ROW> &first, const Sparse_matrix<_CT, COLUMN> &second)
{
    // Perform row-col matrix multiplication with dots
    Sparse_chain<_CT, ROW> row(second._size.second);

    for (size_t index : second._chainsStates)
    {
        _CT tmp(first * second[index]) ;
        if (tmp != 0)
//            row[index] = tmp ;
            row.set_coefficient(index, tmp);
    }
    return row;
}

// Matrix-matrix product
// COLUMN x COLUMN -> ROW
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, COLUMN> &second) {
    Sparse_matrix<_CT, ROW> res(first._size.first, second._size.second);

    for (size_t i = 0 ; i < first._size.first ; i++) {
        Sparse_chain<_CT, ROW> row(second._size.second);

        for (size_t colRight: second._chainsStates) {
            _CT coef = get_row(first, i) * second._chains[colRight];
            if (coef != 0) {
                row[colRight] = coef;
            }
        }

        if (!row.is_null()) {
            res[i] = row;
        }
    }

    return res;
}

// Matrix-matrix product
// ROW x COLUMN -> ROW
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, COLUMN> &second) {
    Sparse_matrix<_CT, ROW> res(first._size.first, second._size.second);

    // Perform row-col matrix multiplication with dot products.
    for (size_t rowLeft: first._chainsStates) {
        Sparse_chain<_CT, ROW> row(second._size.second);

        for (size_t colRight: second._chainsStates) {
            _CT coef = first._chains[rowLeft] * second._chains[colRight];
            if (coef != 0) {
                row[colRight] = coef;
            }
        }

        res[rowLeft] = row;
    }

    return res;
}

// Matrix-matrix product
// COLUMN x ROW -> ROW
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, ROW> &second) {
    Sparse_matrix<_CT, ROW> res(first._size.first, second._size.second);

    // Perform row-col matrix multiplication with dot products.
    for (size_t colLeft: first._chainsStates) {
        res += first._chains[colLeft] % second._chains[colLeft];
    }

    return res;
}

// Matrix-matrix product
// ROW x ROW -> ROW
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, ROW> &second) {
    Sparse_matrix<_CT, ROW> res(first._size.first, second._size.second);

    // Perform row-row matrix multiplication with linear combination of rows.
    for (size_t index: first._chainsStates) {
        Sparse_chain<_CT, ROW> row(second._size.second);

        for (auto colRight: first._chains[index]) {
            if (first._chainsStates.isOn(colRight.first)) {
                row += colRight.second * second._chains[colRight.first];
            }
        }

        res[index] = row;
    }

    return res;
}

// Matrices sum and assign
// COLUMN += ROW
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator+=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, ROW> &other) {
    if (matrix._size != other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }

    for (size_t index = 0 ; index < other._size.second ; index++) {
        Sparse_chain<_CT, COLUMN> column = get_column(other, index);
        if (!column.is_null()) {
            matrix._chainsStates |= index;
            matrix._chains[index] += get_column(other, index);

            if (matrix._chains[index].is_null()) {
                matrix._chainsStates.setOff(index);
            }
        }
    }

    return matrix;
}

// Matrices sum and assign
// ROW += COLUMN
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator+=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, COLUMN> &other) {
    if (matrix._size != other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }

    for (size_t index = 0 ; index < other._size.first ; index++) {
        Sparse_chain<_CT, ROW> row = get_row(other, index);
        if (!row.is_null()) {
            matrix._chainsStates |= index;
            matrix._chains[index] += get_row(other, index);

            if (matrix._chains[index].is_null()) {
                matrix._chainsStates.setOff(index);
            }
        }
    }

    return matrix;
}

// Matrices subtraction and assign
// COLUMN -= ROW
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator-=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, ROW> &other) {
    if (matrix._size != other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }

    for (size_t index = 0 ; index < other._size.second ; index++) {
        Sparse_chain<_CT, COLUMN> column = get_column(other, index);
        if (!column.is_null()) {
            matrix._chainsStates |= index;
            matrix._chains[index] -= get_column(other, index);

            if (matrix._chains[index].is_null()) {
                matrix._chainsStates.setOff(index);
            }
        }
    }

    return matrix;
}

// Matrices subtraction and assign
// ROW -= COLUMN
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator-=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, COLUMN> &other) {
    if (matrix._size != other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }

    for (size_t index = 0 ; index < other._size.second ; index++) {
        Sparse_chain<_CT, ROW> row = get_row(other, index);
        if (!row.is_null()) {
            matrix._chainsStates |= index;
            matrix._chains[index] -= get_row(other, index);

            if (matrix._chains[index].is_null()) {
                matrix._chainsStates.setOff(index);
            }
        }
    }

    return matrix;
}

// Matrices product and assign
// COLUMN += COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, COLUMN> &other) {
    matrix = matrix * other;
    return matrix;
}

// Matrices product and assign
// ROW *= ROW
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, ROW> &other) {
    matrix = matrix % other;
    return matrix;
}

// Matrices product and assign
// COLUMN *= ROW
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &matrix, const Sparse_matrix<_CT, ROW> &other) {
    matrix = matrix * other;
    return matrix;
}

// Matrices product and assign
// ROW *= COLUMN
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &matrix, const Sparse_matrix<_CT, COLUMN> &other) {
    matrix = matrix % other;
    return matrix;
}

// Get column (in COLUMN matrix)
template <typename _CT>
Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index) {
    return matrix._chains[index];
}

// Get column (in ROW matrix)
template <typename _CT>
Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, ROW> &matrix,  size_t index) {
    Sparse_chain<_CT, COLUMN> column(matrix._size.first);
    if (matrix._size.first > 0)
    {
        for (size_t i : matrix._chainsStates) {
            if (!matrix._chains[i].is_null(index)) {
                column.set_coefficient(i, matrix._chains[i][index]);
            }
        }
    }

    return column;
}

// Get row (in COLUMN matrix)
template <typename _CT>
Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index) {
    Sparse_chain<_CT, ROW> row(matrix._size.second);
    if (matrix._size.second > 0)
    {
        for (size_t i : matrix._chainsStates) {
            if (!matrix._chains[i].is_null(index)) {
                row.set_coefficient(i, matrix._chains[i][index]);
            }
        }
    }

    return row;
}

// Get row (in COLUMN matrix)
template <typename _CT>
Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, ROW> &matrix,  size_t index) {
    return matrix._chains[index];
}

// Get constant reference over a column in a column-matrix
template <typename _CT>
const Sparse_chain<_CT, COLUMN> & cget_column(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index)
{
    if (index >= matrix._size.second) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(matrix._size.second) + ".");
    }
    return matrix._chains[index];
}

// Get constant reference over a row in a row-matrix
template <typename _CT>
const Sparse_chain<_CT, ROW> & cget_row(const Sparse_matrix<_CT, ROW> &matrix, const size_t index)
{
    if (index >= matrix._size.first) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(matrix._size.first) + ".");
    }
    return matrix._chains[index];
}

// Set column in a COLUMN matrix
template <typename _CT>
void set_column(Sparse_matrix<_CT, COLUMN> &matrix,  size_t index, const Sparse_chain<_CT, COLUMN> &chain) {
    if(matrix.dimensions().first != chain.dimension())
        throw std::runtime_error("set_column dimension error") ;
    matrix[index] = chain;
    if (chain.is_null())
        matrix._chainsStates.setOff(index) ;
}

// Set column in a ROW matrix
template <typename _CT>
void set_column(Sparse_matrix<_CT, ROW> &matrix,  size_t index, const Sparse_chain<_CT, COLUMN> &chain) {
    if(matrix.dimensions().first != chain.dimension())
        throw std::runtime_error("set_column dimension error") ;
    for (size_t i = 0 ; i < matrix._size.first ; i++) {
        if (!matrix._chains[i].is_null(index) && chain.is_null(i)) {
            matrix._chains[i] /= index;

            if (matrix._chains[i].is_null()) {
                matrix._chainsStates.setOff(i);
            }
        }

        if (!chain.is_null(i)) {
            matrix._chainsStates |= i;
            matrix._chains[i].set_coefficient(index, chain[i]);
        }
    }
}

// Set row in a COLUMN matrix
template <typename _CT>
void set_row(Sparse_matrix<_CT, COLUMN> &matrix,  size_t index, const Sparse_chain<_CT, ROW> &chain) {
    if(matrix.dimensions().second != chain.dimension())
        throw("set_column dimension error") ;
    for (size_t i = 0 ; i < matrix._size.second ; i++) {
        if (!matrix._chains[i].is_null(index) && chain.is_null(i)) {
            matrix._chains[i] /= index;

            if (matrix._chains[i].is_null()) {
                matrix._chainsStates.setOff(i);
            }
        }

        if (!chain.is_null(i)) {
            matrix._chainsStates |= i;
            matrix._chains[i].set_coefficient(index, chain[i]);
        }
    }
}

// Set row in a ROW matrix
template <typename _CT>
void set_row(Sparse_matrix<_CT, ROW> &matrix,  size_t index, const Sparse_chain<_CT, ROW> &chain) {
    if(matrix.dimensions().second != chain.dimension())
        throw("set_column dimension error") ;
    matrix[index] = chain;
    if (chain.is_null())
        matrix._chainsStates.setOff(index) ;
}

template <typename _CT, int _CTF>
inline void set_coefficient(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j, const _CT d)
{
    matrix.set_coefficient(i, j, d);
}

template <typename _CT, int _CTF>
inline _CT get_coefficient(const Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j)
{
    matrix.get_coefficient(i, j);
}

template <typename _CT, int _CTF>
inline Sparse_matrix<_CT, _CTF>& del_column(Sparse_matrix<_CT, _CTF>& matrix, size_t index)
{
    return matrix.del_column(index);
}

template <typename _CT, int _CTF>
inline Sparse_matrix<_CT, _CTF>& del_row(Sparse_matrix<_CT, _CTF>& matrix, size_t index)
{
    return matrix.del_row(index);
}

template <typename _CT, int _CTF>
inline Sparse_matrix<_CT, _CTF>& del_coefficient(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j)
{
    return matrix.del_coefficient(i, j);
}

template <typename _CT>
std::ostream& write_matrix (const Sparse_matrix<_CT, OSM::COLUMN>& M, std::ostream& out)
{
    typedef Sparse_chain<_CT, OSM::COLUMN> Column_chain;
    std::vector<size_t> vec_i, vec_j;
    std::vector<_CT> vec_val;
    // Matrix type : 0 for (COLUMN), 1 for (ROW)
    out << "0" << std::endl ;
    // Size : nb rows / nb cols
    out << M._size.first << " " << M._size.second << std::endl;
    // Get all coefficients
    for(OSM::Bitboard::iterator it = M.begin(); it != M.end(); ++it)
    {
        const Column_chain& col(OSM::cget_column(M, *it));
        // Iterate over the column
        for (typename Column_chain::const_iterator it_col = col.begin(); it_col != col.end(); ++it_col)
        {
            vec_j.push_back(*it) ;
            vec_i.push_back(it_col->first) ;
            vec_val.push_back(it_col->second) ;
        }
    }
    // Output the number of coefficients
    out << vec_i.size() << std::endl ;
    // Output all coefficients : i j val
    for (int n=0; n<vec_i.size(); ++n)
        out << vec_i.at(n) << " " << vec_j.at(n) << " " << vec_val.at(n) << std::endl ;
    return out ;
}

template <typename _CT>
std::ostream& write_matrix (const Sparse_matrix<_CT, OSM::ROW>& M, std::ostream& out)
{
    typedef Sparse_chain<_CT, OSM::ROW> Row_chain;
    std::vector<size_t> vec_i, vec_j;
    std::vector<_CT> vec_val;
    // Matrix type : 0 for (COLUMN), 1 for (ROW)
    out << "1" << std::endl ;
    // Size : nb rows / nb cols
    out << M._size.first << " " << M._size.second << std::endl;
    // Get all coefficients
    for(OSM::Bitboard::iterator it = M.begin(); it != M.end(); ++it)
    {
        const Row_chain& row(OSM::cget_row(M, *it));
        // Iterate over the column
        for (typename Row_chain::const_iterator it_row = row.begin(); it_row != row.end(); ++it_row)
        {
            vec_i.push_back(*it) ;
            vec_j.push_back(it_row->first) ;
            vec_val.push_back(it_row->second) ;
        }
    }
    // Output the number of coefficients
    out << vec_i.size() << std::endl ;
    // Output all coefficients : i j val
    for (int n=0; n<vec_i.size(); ++n)
        out << vec_i.at(n) << " " << vec_j.at(n) << " " << vec_val.at(n) << std::endl ;
    return out ;
}

template <typename _CT>
void write_matrix (const Sparse_matrix<_CT, OSM::COLUMN>& M, std::string filename)
{
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);
    if ( not out . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    CGAL::OSM::write_matrix(M, out) ;

    out.close();
}

template <typename _CT>
void write_matrix (const Sparse_matrix<_CT, OSM::ROW>& M, std::string filename)
{
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);
    if ( not out . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    CGAL::OSM::write_matrix(M, out) ;

    out.close();
}

template <typename _CT>
std::istream& read_matrix (Sparse_matrix<_CT, OSM::COLUMN>& M, std::istream& in)
{
    // Read and check type
    // Matrix type : 0 for (COLUMN), 1 for (ROW)
    int type ;
    in >> type ;
    if (type != 0)
        throw("read_matrix error: trying to load a ROW matrix representation into a COLUMN matrix");
    // Read and adjust size
    // Size : nb rows / nb cols
    size_t nrows, ncols ;
    in >> nrows >> ncols ;
    M = Sparse_matrix<_CT, OSM::COLUMN>(nrows, ncols) ;

    // Read number of coefficients
    size_t n ;
    in >> n ;
    // Read all coefficients and load them into the matrix
    size_t i, j ;
    _CT val ;
    for (size_t k=0; k<n; ++k)
    {
        in >> i >> j ;
        in >> val ;
        OSM::set_coefficient(M, i, j, val) ;
    }
    return in ;
}

template <typename _CT>
std::istream& read_matrix (Sparse_matrix<_CT, OSM::ROW>& M, std::istream& in)
{
    // Read and check type
    // Matrix type : 0 for (COLUMN), 1 for (ROW)
    int type ;
    in >> type ;
    if (type != 1)
        throw("read_matrix error: trying to load a COLUMN matrix representation into a ROW matrix");
    // Read and adjust size
    // Size : nb rows / nb cols
    size_t nrows, ncols ;
    in >> nrows >> ncols ;
    M = Sparse_matrix<_CT, OSM::ROW>(nrows, ncols) ;

    // Read number of coefficients
    size_t n ;
    in >> n ;
    // Read all coefficients and load them into the matrix
    size_t i, j ;
    _CT val ;
    for (size_t k=0; k<n; ++k)
    {
        in >> i >> j ;
        in >> val ;
        OSM::set_coefficient(M, i, j, val) ;
    }
    return in ;
}

template <typename _CT>
void read_matrix (Sparse_matrix<_CT, OSM::COLUMN>& M, std::string filename)
{
    std::ifstream in_file (filename);
    if ( not in_file . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    CGAL::OSM::read_matrix(M, in_file) ;

    in_file.close();
}

template <typename _CT>
void read_matrix (Sparse_matrix<_CT, OSM::ROW>& M, std::string filename)
{
    std::ifstream in_file (filename);
    if ( not in_file . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    CGAL::OSM::read_matrix(M, in_file) ;

    in_file.close();
}

template <typename _CT>
bool operator==(const Sparse_matrix<_CT, OSM::COLUMN>& matrix, const Sparse_matrix<_CT, OSM::COLUMN> &other)
{
    typedef Sparse_chain<_CT, OSM::COLUMN> SparseChainType;
    bool res = true ;
    // Checks that sizes are similar
    res = res && (matrix._size == other._size) ;
    // Checks that all the chains of matrix belong to other
    for (OSM::Bitboard::iterator it = matrix.begin() ; res && (it != matrix.end()); ++it)
    {
        const SparseChainType& chain1(OSM::cget_column(matrix, *it)) ;
        const SparseChainType& chain2(OSM::cget_column(other, *it)) ;
        res = res && (chain1 == chain2) ;
    }
    // Checks that all the chains of other belong to matrix
    for (OSM::Bitboard::iterator it = other.begin() ; res && (it != other.end()); ++it)
    {
        const SparseChainType& chain1(OSM::cget_column(matrix, *it)) ;
        const SparseChainType& chain2(OSM::cget_column(other, *it)) ;
        res = res && (chain1 == chain2) ;
    }
    return res ;
}

template <typename _CT>
bool operator==(const Sparse_matrix<_CT, OSM::ROW>& matrix, const Sparse_matrix<_CT, OSM::ROW> &other)
{
    typedef Sparse_chain<_CT, OSM::ROW> SparseChainType;
    bool res = true ;
    // Checks that sizes are similar
    res = res && (matrix._size == other._size) ;
    // Checks that all the chains of matrix belong to other
    for (OSM::Bitboard::iterator it = matrix.begin() ; it != matrix.end(); ++it)
    {
        const SparseChainType& chain1(OSM::cget_row(matrix, *it)) ;
        const SparseChainType& chain2(OSM::cget_row(other, *it)) ;
        res = res && (chain1 == chain2) ;
    }
    // Checks that all the chains of other belong to matrix
    for (OSM::Bitboard::iterator it = other.begin() ; it != other.end(); ++it)
    {
        const SparseChainType& chain1(OSM::cget_row(matrix, *it)) ;
        const SparseChainType& chain2(OSM::cget_row(other, *it)) ;
        res = res && (chain1 == chain2) ;
    }
    return res ;
}

template <typename _CT>
bool operator==(const Sparse_matrix<_CT, OSM::ROW>& matrix, const Sparse_matrix<_CT, OSM::COLUMN> &other)
{
    return false;
}

template <typename _CT>
bool operator==(const Sparse_matrix<_CT, OSM::COLUMN>& matrix, const Sparse_matrix<_CT, OSM::ROW> &other)
{
    return false;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SPARSE_MATRIX_H
