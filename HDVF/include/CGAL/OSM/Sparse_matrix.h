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


#include "Sparse_chain.h"
#include "Bitboard.h"
#include <stdint.h>
#include <cmath>
#include <unordered_set>

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
 
 \tparam CoefficientType a model of the `Ring` concept, providing the ring used to compute homology.
 \tparam ChainTypeFlag an integer constant encoding the type of matrices (`OSM::COLUMN` or `OSM::ROW`).
*/

template <typename CoefficientType, int ChainTypeFlag>
class Sparse_matrix {
    
public:
    
    /*!
     Type of chains associated to the matrix.
     */
    typedef Sparse_chain<CoefficientType, ChainTypeFlag> MatrixChain;
    
    // Allow the Sparse_matrix class to access other templated Sparse_matrix private members.
    template <typename _CT, int _CTF>
    friend class Sparse_matrix;
    
protected:
    /** \brief The inner chain storage. */
    std::vector<Sparse_chain<CoefficientType, ChainTypeFlag>> _chains;
    
    /** \brief A bitboard containing state of each columns. */
    Bitboard _chainsStates;
    
    /** \brief The matrix size as a (row, column) pair. */
    std::pair<size_t, size_t> _size;
    
    /**
     * \brief Get a reference on a chain from a matrix and change its state (for affectation).
     *
     * Used only internally.
     *
     * \warning The operator changes the status of the chain and the matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The reference to the chain stored at given index.
     */
    MatrixChain& operator[](const size_t _index) {
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
     * Create an empty Matrix of type `ChainTypeFlag` with coefficients of type `CoefficientType`.
     * The default matrix size is 0x0.
     */
    Sparse_matrix() {
        _chains = std::vector<MatrixChain>(0);
        _chainsStates = Bitboard(0);
        _size = {0, 0};
    }
    
    /**
     * \brief Constructor with given rows/columns sizes.
     *
     * Create a new empty Sparse_matrix object of type `ChainTypeFlag` with coefficients of type `CoefficientType` and a given size along rows/columns.
     *
     * \param[in] rowCount The number of rows to preallocate.
     * \param[in] columnCount The number of columns to preallocate.
     */
    Sparse_matrix(const size_t rowCount, const size_t columnCount) {
        size_t mainSize = ChainTypeFlag == COLUMN ? columnCount : rowCount;
        size_t secondarySize = ChainTypeFlag == COLUMN ? rowCount : columnCount;
        
        _chains = std::vector<Sparse_chain<CoefficientType, ChainTypeFlag>>(mainSize);
        for (size_t i = 0 ; i < mainSize ; i++) {
            _chains[i] = Sparse_chain<CoefficientType, ChainTypeFlag>(secondarySize);
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
    Sparse_matrix(const Sparse_matrix<CoefficientType,CTF> &otherToCopy) {
        if (ChainTypeFlag == CTF)
        {
            _chainsStates = otherToCopy._chainsStates;
            _size = otherToCopy._size;
            // Copy of _chains as such
            _chains.resize(otherToCopy._chains._size()) ;
            for (size_t i = 0; i<otherToCopy._chains._size(); ++i)
            {
                const Sparse_chain<CoefficientType, CTF>& tmp(otherToCopy._chains.at(i)) ;
                Sparse_chain<CoefficientType,ChainTypeFlag> res(tmp.dimension()) ;
                for (typename Sparse_chain<CoefficientType, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
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
                _chains.at(i) = Sparse_chain<CoefficientType,ChainTypeFlag>(chain_size) ;
            }
            
            for (size_t i = 0; i<otherToCopy._chains._size(); ++i)
            {
                const Sparse_chain<CoefficientType, CTF>& tmp(otherToCopy._chains.at(i)) ;
                for (typename Sparse_chain<CoefficientType, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
                {
                    _chains[it->first][i] = it->second ;
                    _chainsStates.setOn(it->first) ;
                }
            }
        }
    }
    
    /**
     * \brief Assign to other matrix.
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
     * \brief Clean a SparseMatrix (set all coefficients to zero).
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
     * \brief Test is a SparseMatrix is null.
     *
     * The function return `true` is the SparseMatrix is null (that is, empty) and `false` otherwise.
    */
    bool is_null()
    {
        return (_chainsStates.begin() == _chainsStates.end()) ;
    }
    
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
                CoefficientType tmp = matrix.get_coef(i,j) ;
                
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
    
    
    /**
     * \brief Adds two matrices together into a new matrix.
     *
     * Adds each coefficient of the matrices together and returns a new matrix (of the same type as `first`) representing the result (when possible, prefer `+=` for efficiency).
     *
     * \pre Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     *
     * \warning Will raise an error if the other matrix is not the same `CoefficientType`.
     *
     * \param[in] first The first matrix.
     * \param[in] second The second matrix.
     *
     * \return A new matrix representing the result.
     */
    template <int _CTF>
    friend Sparse_matrix operator+(const Sparse_matrix &first, const Sparse_matrix<CoefficientType, _CTF> &second) {
        Sparse_matrix newMatrix(first);
        newMatrix += second;
        
        return newMatrix;
    }
    
    /**
     * \brief Substracts two matrices together into a new matrix.
     *
     * Substracts each coefficient of the matrix together and returns a new matrix (of the same type as `first`) representing the result (when possible, prefer `-=` for efficiency).
     *
     * \pre Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     *
     * \warning Will raise an error if the other matrix is not the same `CoefficientType`.
     *
     * \param[in] first The first matrix.
     * \param[in] second The second matrix.
     *
     * \return A new matrix representing the result.
     */
    template <int _CTF>
    friend Sparse_matrix operator-(const Sparse_matrix &first, const Sparse_matrix<CoefficientType, _CTF> &second) {
        Sparse_matrix newMatrix = first;
        newMatrix -= second;
        
        return newMatrix;
    }
    
    /**
     * \brief Apply factor on each coefficients into a new matrix.
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
     * \brief Apply factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`.
     *
     * If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     *
     * \param[in] matrix The matrix.
     * \param[in] lambda The factor to apply.
     *
     * \return A new matrix representing the result.
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix operator*(const Sparse_matrix<_CT, _CTF> &matrix, const _CT& lambda) {
        Sparse_matrix newMatrix = matrix;
        newMatrix *= lambda;
        
        return newMatrix;
    }
    
    /**
     * \defgroup MatrixMatrixProdCol Matrices product (with column-based result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between matrices and returns a new column-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new column-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. The product is optimized (using standard definition or block products) for each combination of `ChainTypeFlag`. However, efficiency depends on `ChainTypeFlag` (when possible, prefer row-major by column-major products).
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
    
    /**
     * \defgroup MatrixChainProd Matrix / column-chain product (with column-chain result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between a sparse matrix (column or row major) and a column-chain. The function returns a new column-major chain.
     *
     * Perform standard linear algebra product between a matrix and a column-chain (ie. matrix / column vector product) and returns a new column-major chain. Both arguments must have the same `CoefficientType` but the matrix can have any `ChainTypeFlag` (and the product is optimized for each of them).
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
    
    /**
     * \defgroup ChainMatrixProd Row-chain / matrix product (with row-chain result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between a row-chain and a sparse matrix (column or row major). The function returns a new row-major chain.
     *
     * Perform standard linear algebra product between a row-chain and a matrix (ie. row vector / matrix product) and returns a new row-major chain. Both arguments must have the same `CoefficientType` but the matrix can have any `ChainTypeFlag` (and the product is optimized for each of them).
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
    
    /**
     * \defgroup MatrixMatrixProdRow Matrices product (with row-based result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between matrices and returns a new row-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new row-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. The product is optimized (using standard definition or block products) for each combination of `ChainTypeFlag`. However, efficiency depends on `ChainTypeFlag`.
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
    
    /**
     * \defgroup MatrixMatrixAddAssign Sum matrices and assign.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform addition between matrices and assigns the result to `matrix`.
     *
     * Perform standard linear algebra addition. Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. The addition is optimized for each combination of `ChainTypeFlag`.
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
    
    /**
     * \defgroup MatrixMatrixSubtractAssign Subtract matrices and assign.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform subtraction between matrices and assigns the result to `matrix`.
     *
     * Perform standard linear algebra subtraction. Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. The addition is optimized for each combination of `ChainTypeFlag`.
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
     * \brief Apply factor on each coefficients and assign.
     *
     * If `lambda` is 0, this comes to nullify the matrix.
     *
     * \param[in] lambda The factor to apply.
     *
     * \return The modified matrix representing the result.
     */
    Sparse_matrix& operator*=(const CoefficientType& lambda) {
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
     * \brief Compute the negative of a matrix (unary operator).
     *
     * \return The resulting matrix.
     */
    friend Sparse_matrix operator-(const Sparse_matrix& matrix) {
        Sparse_matrix res(matrix._size.first, matrix._size.second) ;
        
        for (size_t index: matrix._chainsStates)
        {
            const MatrixChain& tmp_chain(matrix._chains[index]) ;
            for (typename MatrixChain::const_iterator it = tmp_chain.cbegin(); it != tmp_chain.cend(); ++it)
                res[index][it->first] = -it->second ;
        }
        return res ;
    }
    
    /**
     * \defgroup MatrixMatrixProdAssign Multiplies matrices and assign.
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform product between matrices and assigns the result to `matrix`.
     *
     * Perform standard linear algebra product. Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. The product is optimized for each combination of `ChainTypeFlag`.
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
     * \brief Get the value of a chain from a const matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] index The chain index.
     *
     * \return The chain stored at given index.
     */
    MatrixChain operator[](size_t index) const {
        if (ChainTypeFlag == COLUMN && index >= _size.second) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.second) + ".");
        }
        if (ChainTypeFlag == ROW && index >= _size.first) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.first) + ".");
        }
        
        return _chains[index];
    }
    
protected:
    // Protected method for set_coef
    void set_coef(const size_t i, const size_t j, const CoefficientType d) {
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
            del_coef(i, j);
        }
    }
public:
    /**
     * \brief Set a given coefficient in `matrix`.
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
    friend void set_coef(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j, const _CT d);
    
protected:
    // Protected method for get_coef
    CoefficientType get_coef(const size_t i, const size_t j) const {
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
    /**
     * \brief Get a given coefficient.
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
    friend _CT get_coef(const Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j);
    
    /**
     * \defgroup GetColumn Get a column.
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
    
    /** \brief Get a column from a COLUMN matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index);
    
    /** \brief Get a column from a ROW matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, ROW> &matrix,  size_t index);
    
    /** @} */
    
    /**
     * \defgroup GetRow Get a row.
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
    
    /** \brief Get a row from a COLUMN matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, COLUMN> &matrix,  size_t index);
    
    /** \brief Get a row from a ROW matrix. */
    template <typename _CT>
    friend Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, ROW> &matrix,  size_t index);
    
    /** @} */
    
    /**
     * \brief Get a const reference over a column from a column matrix.
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
    
    /**
     * \brief Get a constant reference over a row from a row matrix
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
    
    
    /**
     * \defgroup SetColumn Set a column.
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
    
    /** \brief Set a column in a COLUMN matrix. */
    template <typename _CT>
    friend void set_column(Sparse_matrix<_CT, COLUMN> &matrix,  size_t index, const Sparse_chain<_CT, COLUMN> &chain);
    
    /** \brief Set a column in a ROW matrix. */
    template <typename _CT>
    friend void set_column(Sparse_matrix<_CT, ROW> &matrix,  size_t index, const Sparse_chain<_CT, COLUMN> &chain);
    
    /** @} */
    
    /**
     * \defgroup SetRow Set a row.
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
    
    /** \brief Set a row in a COLUMN matrix. */
    template <typename _CT>
    friend void set_row(Sparse_matrix<_CT, COLUMN> &matrix,  size_t index, const Sparse_chain<_CT, ROW> &chain);
    
    /** \brief Set a row in a ROW matrix. */
    template <typename _CT>
    friend void set_row(Sparse_matrix<_CT, ROW> &matrix,  size_t index, const Sparse_chain<_CT, ROW> &chain);
    
    /** @} */
    
    /**
     * \defgroup GetBlockMatrix Get a sub-matrix by removing chains in a matrix.
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

    /** \brief Remove a set of chains from a copy of a matrix. */
    friend Sparse_matrix operator/(const Sparse_matrix &matrix, const std::vector<size_t> &_indices) {
        Sparse_matrix res(matrix);
        res /= _indices;
        return res;
    }
    
    /** \brief Remove the chain at a given `index` from a copy of a matrix. */
    friend Sparse_matrix operator/(const Sparse_matrix &matrix, size_t index) {
        Sparse_matrix res(matrix);
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
    
    /** \brief Remove a set of chains from a matrix. */
    Sparse_matrix& operator/=(const std::vector<size_t> &indexes) {
        for (size_t index : indexes) {
            *this /= index;
        }
        
        return *this;
    }
    
    /** \brief Remove the chain at a given `index` from a matrix. */
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
                MatrixChain &tmp(_chains[ind]) ;
                tmp/=tmp_id ;
                // If the row index has become empty: update
                if (tmp.is_null())
                    *this /= std::vector<size_t>({ind}) ;
            }
        }
        
        return *this;
    }
public:
    /**
     * \brief Remove a column from the matrix.
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
                MatrixChain &tmp(_chains[ind]);
                tmp /= tmp_id;
                // If the column index has become empty: update
                if (tmp.is_null())
                    *this /= std::vector<size_t>({ind}) ;
            }
        }
        
        return *this;
    }
    
public:
    /**
     * \brief Remove a row from the matrix.
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
    // Protected version of del_coef
    Sparse_matrix& del_coef(size_t i, size_t j)
    {
        
        if (ChainTypeFlag == OSM::COLUMN) {
            std::vector<size_t> tmp_id({i}) ;
            MatrixChain &tmp(_chains[j]);
            tmp /= tmp_id;
            if (tmp.is_null())
                _chainsStates.setOff(j) ;
        } else // OSM::ROW
        {
            std::vector<size_t> tmp_id({j}) ;
            MatrixChain &tmp(_chains[i]);
            tmp /= tmp_id;
            if (tmp.is_null())
                _chainsStates.setOff(i) ;
        }
        
        return *this;
    }
    
public:
    /**
     * \brief Remove a coefficient from the matrix.
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
    friend Sparse_matrix<_CT, _CTF>& del_coef(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j);
    
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
     * \brief Transpose a matrix.
     *
     * \return A new matrix where the chain type flag has been swapped between COLUMN and ROW and data chains have been transposed.
     */
    Sparse_matrix<CoefficientType, COLUMN + ROW - ChainTypeFlag> transpose() {
        Sparse_matrix<CoefficientType, COLUMN + ROW - ChainTypeFlag> transposed(this->_size.second, this->_size.first);
        
        for (size_t index : this->_chainsStates) {
            transposed._chains[index] = this->_chains[index].transpose();
        }
        
        transposed._chainsStates = this->_chainsStates;
        
        return transposed;
    }
    
    /**
     * \brief Get the matrix sizes.
     *
     * \return The matrix size as a row/column pair.
     */
    std::pair<size_t, size_t> dimensions() const {
        return this->_size;
    }
};


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
Sparse_chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_chain<_CT, COLUMN> &_second)
{
    // Perform col-col matrix multiplication with linear combination of columns.
    Sparse_chain<_CT, COLUMN> column(_first._size.first);
    
    for (typename Sparse_chain<_CT, COLUMN>::const_iterator it = _second.begin(); it != _second.end(); ++it)
    {
        column += it->second * _first._chains[it->first];
    }
    
    return column;
}

// Matrix - column chain product
// ROW matrix
template <typename _CT>
Sparse_chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Sparse_chain<_CT, COLUMN> &_second)
{
    // Perform row-col matrix multiplication with dots
    Sparse_chain<_CT, COLUMN> column(_first._size.first);
    
    for (size_t index : _first._chainsStates)
    {
        _CT tmp(_first[index] * _second) ;
        if (tmp != 0)
            column[index] = tmp ;
    }
    return column;
}

// Row chain - matrix product
// ROW matrix
template <typename _CT>
Sparse_chain<_CT, ROW> operator*(const Sparse_chain<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second)
{
    // Perform row-row matrix multiplication with linear combination of rows.
    Sparse_chain<_CT, ROW> row(_second._size.second);
    
    for (typename Sparse_chain<_CT, ROW>::const_iterator it = _first.begin(); it != _first.end(); ++it)
    {
        row += it->second * _second._chains[it->first];
    }
    
    return row;
}

// Row chain - matrix product
// COLUMN matrix
template <typename _CT>
Sparse_chain<_CT, ROW> operator*(const Sparse_chain<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second)
{
    // Perform row-col matrix multiplication with dots
    Sparse_chain<_CT, ROW> row(_second._size.second);
    
    for (size_t index : _second._chainsStates)
    {
        _CT tmp(_first * _second[index]) ;
        if (tmp != 0)
            row[index] = tmp ;
    }
    return row;
}

// Matrix-matrix product
// COLUMN x COLUMN -> ROW
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, COLUMN> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    for (size_t i = 0 ; i < _first._size.first ; i++) {
        Sparse_chain<_CT, ROW> row(_second._size.second);
        
        for (size_t colRight: _second._chainsStates) {
            _CT coef = get_row(_first, i) * _second._chains[colRight];
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
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (size_t rowLeft: _first._chainsStates) {
        Sparse_chain<_CT, ROW> row(_second._size.second);
        
        for (size_t colRight: _second._chainsStates) {
            _CT coef = _first._chains[rowLeft] * _second._chains[colRight];
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
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, ROW> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (size_t colLeft: _first._chainsStates) {
        res += _first._chains[colLeft] % _second._chains[colLeft];
    }
    
    return res;
}

// Matrix-matrix product
// ROW x ROW -> ROW
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    // Perform row-row matrix multiplication with linear combination of rows.
    for (size_t index: _first._chainsStates) {
        Sparse_chain<_CT, ROW> row(_second._size.second);
        
        for (auto colRight: _first._chains[index]) {
            if (_first._chainsStates.isOn(colRight.first)) {
                row += colRight.second * _second._chains[colRight.first];
            }
        }
        
        res[index] = row;
    }
    
    return res;
}

// Matrices sum and assign
// COLUMN += ROW
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator+=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (size_t index = 0 ; index < _other._size.second ; index++) {
        Sparse_chain<_CT, COLUMN> column = get_column(_other, index);
        if (!column.is_null()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] += get_column(_other, index);
            
            if (_matrix._chains[index].is_null()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

// Matrices sum and assign
// ROW += COLUMN
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator+=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (size_t index = 0 ; index < _other._size.first ; index++) {
        Sparse_chain<_CT, ROW> row = get_row(_other, index);
        if (!row.is_null()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] += get_row(_other, index);
            
            if (_matrix._chains[index].is_null()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

// Matrices subtraction and assign
// COLUMN -= ROW
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator-=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (size_t index = 0 ; index < _other._size.second ; index++) {
        Sparse_chain<_CT, COLUMN> column = get_column(_other, index);
        if (!column.is_null()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] -= get_column(_other, index);
            
            if (_matrix._chains[index].is_null()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

// Matrices subtraction and assign
// ROW -= COLUMN
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator-=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (size_t index = 0 ; index < _other._size.second ; index++) {
        Sparse_chain<_CT, ROW> row = get_row(_other, index);
        if (!row.is_null()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] -= get_row(_other, index);
            
            if (_matrix._chains[index].is_null()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

// Matrices product and assign
// COLUMN += COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    _matrix = _matrix * _other;
    return _matrix;
}

// Matrices product and assign
// ROW += ROW
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    _matrix = _matrix % _other;
    return _matrix;
}

// Matrices product and assign
// COLUMN += ROW
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    _matrix = _matrix * _other;
    return _matrix;
}

// Matrices product and assign
// ROW += COLUMN
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    _matrix = _matrix % _other;
    return _matrix;
}

// Get column (in COLUMN matrix)
template <typename _CT>
Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, COLUMN> &_matrix,  size_t _index) {
    return _matrix._chains[_index];
}

// Get column (in ROW matrix)
template <typename _CT>
Sparse_chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, ROW> &_matrix,  size_t _index) {
    Sparse_chain<_CT, COLUMN> column(_matrix._size.first);
    if (_matrix._size.first > 0)
    {
        for (size_t i : _matrix._chainsStates) {
            if (!_matrix._chains[i].is_null(_index)) {
                column[i] = _matrix._chains[i][_index];
            }
        }
    }
    
    return column;
}

// Get row (in COLUMN matrix)
template <typename _CT>
Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, COLUMN> &_matrix,  size_t _index) {
    Sparse_chain<_CT, ROW> row(_matrix._size.second);
    if (_matrix._size.second > 0)
    {
        for (size_t i : _matrix._chainsStates) {
            if (!_matrix._chains[i].is_null(_index)) {
                row[i] = _matrix._chains[i][_index];
            }
        }
    }
    
    return row;
}

// Get row (in COLUMN matrix)
template <typename _CT>
Sparse_chain<_CT, ROW> get_row(const Sparse_matrix<_CT, ROW> &_matrix,  size_t _index) {
    return _matrix._chains[_index];
}

// Get constant reference over a column in a column-matrix
template <typename _CT>
const Sparse_chain<_CT, COLUMN> & cget_column(const Sparse_matrix<_CT, COLUMN> &_matrix,  size_t _index)
{
    if (_index >= _matrix._size.second) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(_matrix._size.second) + ".");
    }
    return _matrix._chains[_index];
}

// Get constant reference over a row in a row-matrix
template <typename _CT>
const Sparse_chain<_CT, ROW> & cget_row(const Sparse_matrix<_CT, ROW> &_matrix, const size_t _index)
{
    if (_index >= _matrix._size.first) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(_matrix._size.first) + ".");
    }
    return _matrix._chains[_index];
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
            matrix._chains[i][index] = chain[i];
        }
    }
}

// Set row in a COLUMN matrix
template <typename _CT>
void set_row(Sparse_matrix<_CT, COLUMN> &_matrix,  size_t _index, const Sparse_chain<_CT, ROW> &_chain) {
    if(_matrix.dimensions().second != _chain.dimension())
        throw("set_column dimension error") ;
    for (size_t i = 0 ; i < _matrix._size.second ; i++) {
        if (!_matrix._chains[i].is_null(_index) && _chain.is_null(i)) {
            _matrix._chains[i] /= _index;
            
            if (_matrix._chains[i].is_null()) {
                _matrix._chainsStates.setOff(i);
            }
        }
        
        if (!_chain.is_null(i)) {
            _matrix._chainsStates |= i;
            _matrix._chains[i][_index] = _chain[i];
        }
    }
}

// Set row in a ROW matrix
template <typename _CT>
void set_row(Sparse_matrix<_CT, ROW> &_matrix,  size_t _index, const Sparse_chain<_CT, ROW> &_chain) {
    if(_matrix.dimensions().second != _chain.dimension())
        throw("set_column dimension error") ;
    _matrix[_index] = _chain;
    if (_chain.is_null())
        _matrix._chainsStates.setOff(_index) ;
}

template <typename _CT, int _CTF>
void set_coef(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j, const _CT d)
{
    matrix.set_coef(i, j, d);
}

template <typename _CT, int _CTF>
_CT get_coef(const Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j)
{
    matrix.get_coef(i, j);
}

template <typename _CT, int _CTF>
Sparse_matrix<_CT, _CTF>& del_column(Sparse_matrix<_CT, _CTF>& matrix, size_t index)
{
    return matrix.del_column(index);
}

template <typename _CT, int _CTF>
Sparse_matrix<_CT, _CTF>& del_row(Sparse_matrix<_CT, _CTF>& matrix, size_t index)
{
    return matrix.del_row(index);
}

template <typename _CT, int _CTF>
Sparse_matrix<_CT, _CTF>& del_coef(Sparse_matrix<_CT, _CTF>& matrix, size_t i, size_t j)
{
    return matrix.del_coef(i, j);
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SPARSE_MATRIX_H
