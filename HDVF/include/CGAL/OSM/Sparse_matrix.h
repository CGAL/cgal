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


#include "Chain.h"
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
 
 The class `Sparse_matrix` implements the concept `SparseMatrix`, that is, sparse matrices optimized to topological computations. It provides standard linear algebra operators and fast iterators and block operations (set, get and nullify) which are required to implement efficiently HDVFs.
 
 The implementation is based on mapped sparse matrices. Hence matrices of the `Sparse_matrix` class are either column of row major (the `ChainTypeFlag` parameter determines the type). A column-major (resp. row-major) `Sparse_matrix` is a vector of `Sparse_chains` which encode columns (res. rows). Moreover, in order to efficiently iterate over non empty columns (resp. rows) the `Bitboard` data structure implements the concept `SparseMatrix::NonZeroChainIndices`. A bitboard is basically a bucket of bits recording the indices of non empty chains. However, this data structure has been designed in order to efficiently remove or add indices, as well as  provide efficient iterators to visit non empty chains.
 
 For instance, let us consider  the \f$3\times 4\f$ matrix:
 \f[
 A = \left(\begin{array}{cccc}
 1 & \cdot & \cdot & \cdot \\
 -1 & \cdot & 2 & \cdot\\
 \cdot & \cdot & 1 & \cdot \\
 \cdot & \cdot & \cdot & \cdot \\
 \end{array}\right)
 \f]
 where \f$\cdot\f$ means \f$0\f$.
 
 If the matrix is stored column-major, following data structures are created:
 
 <img src="Sparse_matrix_example_col.pdf" align="center" width=20%/>
 
 \cgalModels{SparseMatrix}
 
 \tparam CoefficientType a model of the `Ring` concept (by default, we use the `Z` model) providing the ring used to compute homology.
 \tparam ChainTypeFlag an integer constant encoding the type of matrices (`OSM::COLUMN` or `OSM::ROW`).
*/

template <typename CoefficientType, int ChainTypeFlag>
class Sparse_matrix {
    
public:
    
    /*!
     Type of chains associated to the matrix.
     */
    typedef Chain<CoefficientType, ChainTypeFlag> MatrixChain;
    
    // Allow the Sparse_matrix class to access other templated Sparse_matrix private members.
    template <typename _CT, int _CTF>
    friend class Sparse_matrix;
    
protected:
    /** \brief The inner chain storage. */
    std::vector<Chain<CoefficientType, ChainTypeFlag>> _chains;
    
    /** \brief A bitboard containing state of each columns. */
    Bitboard _chainsStates;
    
    /** \brief The matrix size as a (row, column) pair. */
    std::pair<int, int> _size;
    
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
    MatrixChain& operator[](const int _index) {
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
     * \brief Create an empty new SparseMatrix object.
     *
     * Default constructor, initialize an empty Matrix of type `ChainTypeFlag` with coefficients of type `CoefficientType`.
     * The default matrix size is 0x0.
     */
    Sparse_matrix() {
        _chains = std::vector<MatrixChain>(0);
        _chainsStates = Bitboard(0);
        _size = {0, 0};
    }
    
    /**
     * \brief Create a new SparseMatrix object with given rows/columns sizes.
     *
     * Constructor with sizes, initialize an empty Matrix of type `ChainTypeFlag` with coefficients of type `CoefficientType` and a given size along rows/columns.
     *
     * \param[in] rowCount The number of rows to preallocate.
     * \param[in] columnCount The number of columns to preallocate.
     */
    Sparse_matrix(const int rowCount, const int columnCount) {
        int mainSize = ChainTypeFlag == COLUMN ? columnCount : rowCount;
        int secondarySize = ChainTypeFlag == COLUMN ? rowCount : columnCount;
        
        _chains = std::vector<Chain<CoefficientType, ChainTypeFlag>>(mainSize);
        for (int i = 0 ; i < mainSize ; i++) {
            _chains[i] = Chain<CoefficientType, ChainTypeFlag>(secondarySize);
        }
        
        _chainsStates = Bitboard(mainSize);
        _size = {rowCount, columnCount};
    }
    
    /**
     * \brief Create a new SparseMatrix from another SparseMatrix object (with possibly a different `ChainTypeFlag`).
     *
     * Copy constructor, initialize a SparseMatrix of same sizes, containing the same coefficients (but not necessarly of the same `ChainTypeFlag`).
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
            for (int i = 0; i<otherToCopy._chains._size(); ++i)
            {
                const Chain<CoefficientType, CTF>& tmp(otherToCopy._chains.at(i)) ;
                Chain<CoefficientType,ChainTypeFlag> res(tmp.dimension()) ;
                for (typename Chain<CoefficientType, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
                {
                    res[it->first] = it->second ;
                }
                _chains[i] = res ;
            }
        }
        else // Copy coefficient by coefficient
        {
            _size = otherToCopy._size;
            int vec_size = (ChainTypeFlag == OSM::COLUMN)?_size.second:_size.first ;
            int chain_size = (ChainTypeFlag == OSM::COLUMN)?_size.first:_size.second ;
            _chains.resize(vec_size) ;
            _chainsStates = OSM::Bitboard(vec_size) ;
            for (int i=0; i<vec_size; ++i)
            {
                _chains.at(i) = Chain<CoefficientType,ChainTypeFlag>(chain_size) ;
            }
            
            for (int i = 0; i<otherToCopy._chains._size(); ++i)
            {
                const Chain<CoefficientType, CTF>& tmp(otherToCopy._chains.at(i)) ;
                for (typename Chain<CoefficientType, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
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
     * \param[in] _otherToCopy The matrix we want to copy.
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
        std::vector<int> coefs ;
        for (OSM::Bitboard::iterator it = begin(); it != end(); ++it)
        {
            coefs.push_back(*it) ;
        }
        *this /= coefs ;
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
        for (int index : matrix._chainsStates) {
            empty = false;
            stream << index << ": " << matrix._chains[index] << ", ";
        }
        
        if (!empty) {
            stream << "\b\b";
        }
        stream << "]" << std::endl ;
#else
        for (int i=0; i<matrix._size.first; ++i)
        {
            for (int j=0; j<matrix._size.second; ++j)
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
     * \defgroup MatrixProdCol Matrix product (with column-based result).
     * \ingroup PkgHDVFAlgorithmClasses
     * @brief  Perform multiplication between matrices and returns a new column-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new column-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. Efficiency of the product depends of `ChainTypeFlag` (when possible, prefer row-major by column-major products).
     *
     * @param first The first matrix.
     * @param second The second matrix.
     * @return The result of the matrix multiplication, column-based.
     * @{
     */
    
    //* COLUMN x COLUMN -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, COLUMN> &second);
    
    //* ROW x COLUMN -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, COLUMN> &second);
    
    //* COLUMN x ROW -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, ROW> &second);
    
    //* ROW x ROW -> COLUMN */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, ROW> &second);

    /** @} */
    
    
    /**
     * \brief Perform multiplication between a column-based matrix and a column-based chain.
     *
     * Generate a column-based chain from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The matrix.
     * \param[in] _second The chain.
     *
     * \return The result of the matrix multiplication, column-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Chain<_CT, COLUMN> &_second);
    
    /**
     * \brief Perform multiplication between a row-based matrix and a column-based chain.
     *
     * Generate a column-based chain from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The matrix.
     * \param[in] _second The chain.
     *
     * \return The result of the matrix multiplication, column-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Chain<_CT, COLUMN> &_second);
    
    /**
     * \brief Perform multiplication between a row-based chain and a row-based matrix.
     *
     * Generate a row-based chain from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The matrix.
     * \param[in] _second The chain.
     *
     * \return The result of the matrix multiplication, row-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, ROW> operator*(const Chain<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second) ;
    
    /**
     * \brief Perform multiplication between a row-based chain and a column-based matrix.
     *
     * Generate a row-based chain from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The matrix.
     * \param[in] _second The chain.
     *
     * \return The result of the matrix multiplication, row-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, ROW> operator*(const Chain<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second) ;
    
    /**
     * \brief Perform matrix multiplication between two matrices.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The first matrix.
     * \param[in] _second The second matrix.
     *
     * \return The result of the matrix multiplication, column-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, COLUMN> &_second);
    
    /**
     * \brief Perform matrix multiplication between two matrices.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The first matrix.
     * \param[in] _second The second matrix.
     *
     * \return The result of the matrix multiplication, column-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second);
    
    /**
     * \brief Perform matrix multiplication between two _chains.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The first matrix.
     * \param[in] _second The second matrix.
     *
     * \return The result of the matrix multiplication, column-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, ROW> &_second);
    
    /**
     * \brief Perform matrix multiplication between two _chains.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _first The first matrix.
     * \param[in] _second The second matrix.
     *
     * \return The result of the matrix multiplication, column-based.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second);
    
    /**
     * \brief Add a matrix and assign.
     *
     * Adds each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& operator+=(const Sparse_matrix &_other) {
        if (this->_size != _other._size) {
            throw std::runtime_error("Matrices must be the same _size.");
        }
        
        for (int index: _other._chainsStates) {
            this->_chainsStates |= index;
            this->_chains[index] += _other._chains[index];
            
            if (this->_chains[index].isNull()) {
                this->_chainsStates.setOff(index);
            }
        }
        
        return *this;
    }
    
    /**
     * \brief Add a matrix and assign.
     *
     * Adds each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator+=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other);
    
    /**
     * \brief Add a matrix and assign.
     *
     * Adds each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator+=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other);
    
    /**
     * \brief Substract a matrix and assign.
     *
     * Substracts each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& operator-=(const Sparse_matrix &_other) {
        if (this->_size != _other._size) {
            throw std::runtime_error("Matrices must be the same _size.");
        }
        
        for (int index: _other._chainsStates) {
            this->_chainsStates |= index;
            this->_chains[index] -= _other._chains[index];
            
            if (this->_chains[index].isNull()) {
                this->_chainsStates.setOff(index);
            }
        }
        
        return *this;
    }
    
    
    /**
     * \brief Substract a matrix and assign.
     *
     * Substracts each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator-=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other);
    
    /**
     * \brief Substract a matrix and assign.
     *
     * Substracts each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator-=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other);
    
    /**
     * \brief Apply factor on each coefficients and assign.
     *
     * \warning Will raise an error if _lambda is 0.
     *
     * \param[in] _lambda The factor to apply.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& operator*=(const CoefficientType& _lambda) {
        if (_lambda == 0) {
            std::pair<int, int> p(this->dimensions()) ;
            *this = Sparse_matrix(p.first, p.second);
            return *this;
        }
        
        for (int index: this->_chainsStates) {
            this->_chains[index] *= _lambda;
        }
        
        return *this;
    }
    
    /**
     * \brief Compute the negative of a matrix (unary operator).
     *
     * \return The resulting matrix.
     *
     * \see \link OSM::Sparse_matrix \endlink
     *
     * \author Bac A.
     * \version 0.1.0
     * \date 07/09/2024
     */
    friend Sparse_matrix operator-(const Sparse_matrix& _matrix) {
        Sparse_matrix res(_matrix._size.first, _matrix._size.second) ;
        
        for (int index: _matrix._chainsStates)
        {
            const MatrixChain& tmp_chain(_matrix._chains[index]) ;
            for (typename MatrixChain::const_iterator it = tmp_chain.cbegin(); it != tmp_chain.cend(); ++it)
                res[index][it->first] = -it->second ;
        }
        return res ;
    }
    
    /**
     * \brief Multiply a matrix and assign.
     *
     * Multiply each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other);
    
    /**
     * \brief Multiply a matrix and assign.
     *
     * Multiply each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, ROW> &_other);
    
    /**
     * \brief Multiply a matrix and assign.
     *
     * Multiply each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other);
    
    /**
     * \brief Multiply a matrix and assign.
     *
     * Multiply each coefficient of the matrix together.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two _chains are not the same coefficient type.
     *
     * \param[in] _matrix The matrix to reassign.
     * \param[in] _other The other matrix.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    template <typename _CT>
    friend Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other);
    
    /**
     * \brief Get a chain from a const matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The chain stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    MatrixChain operator[](const int _index) const {
        if (ChainTypeFlag == COLUMN && _index >= _size.second) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.second) + ".");
        }
        if (ChainTypeFlag == ROW && _index >= _size.first) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(_size.first) + ".");
        }
        
        return _chains[_index];
    }
    
    /**
     * \brief Set a given coefficient.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _i The row index.
     * \param[in] _j The column index.
     * \param[in] _d The value.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    void set_coef(const int _i, const int _j, const CoefficientType _d) {
        if (_i >= _size.first) {
            throw std::runtime_error("Provided _i index should be less than " + std::to_string(_size.first) + ".");
        }
        if (_j >= _size.second) {
            throw std::runtime_error("Provided _j index should be less than " + std::to_string(_size.second) + ".");
        }
        
        
        if (ChainTypeFlag == COLUMN)
        {
            (*this)[_j][_i] = _d ;
        }
        else
            (*this)[_i][_j] = _d ;
    }
    
    /**
     * \brief Get a given coefficient.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _i The row index.
     * \param[in] _j The column index.
     *
     * \return The value of the given coefficient.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    
    CoefficientType get_coef(const int _i, const int _j) const {
        if (_i >= _size.first) {
            //std::cout << "Provided _i index should be less than " << _size.first << "." << std::endl;
            throw std::runtime_error("Provided _i index should be less than " + std::to_string(_size.first) + ".");
        }
        if (_j >= _size.second) {
            //std::cout << "Provided _j index should be less than " << _size.second << "." << std::endl;
            throw std::runtime_error("Provided _j index should be less than " + std::to_string(_size.second) + ".");
        }
        
        
        if (ChainTypeFlag == COLUMN)
            return (*this)[_j][_i] ;
        else // ROW
            return (*this)[_i][_j] ;
    }
    
    /**
     * \brief Get a column from the matrix, even if matrix is a row-chain matrix.
     *
     * \note For column-chain matrixes, it is equivalent to operator[].
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The column stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, COLUMN> &_matrix, const int _index);
    
    /**
     * \brief Get a column from the matrix, even if matrix is a row-chain matrix.
     *
     * \note For column-chain matrixes, it is equivalent to operator[].
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The column stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, ROW> &_matrix, const int _index);
    
    /**
     * \brief Get a row from the matrix, even if matrix is a row-chain matrix.
     *
     * \note For row-chain matrixes, it is equivalent to operator[].
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The row stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, ROW> get_row(const Sparse_matrix<_CT, COLUMN> &_matrix, const int _index);
    
    /**
     * \brief Get a row from the matrix, even if matrix is a row-chain matrix.
     *
     * \note For row-chain matrixes, it is equivalent to operator[].
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The row stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend Chain<_CT, ROW> get_row(const Sparse_matrix<_CT, ROW> &_matrix, const int _index);
    
    /**
     * \brief Get a const reference over a column from a column matrix
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return A const reference over the column stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend const Chain<_CT, COLUMN> & cget_column(const Sparse_matrix<_CT, COLUMN> &_matrix, const int _index);
    
    /**
     * \brief Get a const reference over a row from a row matrix
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return A const reference over the row stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend const Chain<_CT, ROW> & cget_row(const Sparse_matrix<_CT, ROW> &_matrix, const int _index);
    
    
    /**
     * \brief Set a column from the matrix, even if matrix is a row-chain matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     * \param[in] _chain The new column.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend void set_column(Sparse_matrix<_CT, COLUMN> &_matrix, const int _index, const Chain<_CT, COLUMN> &_chain);
    
    /**
     * \brief Set a column from the matrix, even if matrix is a row-chain matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     * \param[in] _chain The new column.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend void set_column(Sparse_matrix<_CT, ROW> &_matrix, const int _index, const Chain<_CT, COLUMN> &_chain);
    
    /**
     * \brief Set a row from the matrix, even if matrix is a row-chain matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     * \param[in] _chain The new row.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend void set_row(Sparse_matrix<_CT, COLUMN> &_matrix, const int _index, const Chain<_CT, ROW> &_chain);
    
    /**
     * \brief Set a row from the matrix, even if matrix is a row-chain matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     * \param[in] _chain The new row.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    template <typename _CT>
    friend void set_row(Sparse_matrix<_CT, ROW> &_matrix, const int _index, const Chain<_CT, ROW> &_chain);
    
    /**
     * \brief Get a submatrix from the matrix.
     *
     * Removes all indexes provided in the vector from the matrix and returns it.
     *
     * \note Will return a copy of the matrix if given vector is empty.
     * \note The submatrix will be rectangular.
     *
     * \param[in] _matrix The matrix to process.
     * \param[in] _indexes The indexes to remove.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    friend Sparse_matrix operator/(const Sparse_matrix &_matrix, const std::vector<int> &_indexes) {
        Sparse_matrix res(_matrix);
        res /= _indexes;
        return res;
    }
    
    /**
     * \brief Get a submatrix from the matrix.
     *
     * Nullifies a chain along the major direction of a copy of the matrix returns it.
     *
     * \note Will return a copy of the matrix.
     * \note The submatrix will be rectangular.
     *
     * \param[in] _matrix The matrix to process.
     * \param[in] _index The index to remove.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    friend Sparse_matrix operator/(const Sparse_matrix &_matrix, const int _index) {
        Sparse_matrix res(_matrix);
        res /= _index;
        return res;
    }
    
    /**
     * \brief Get a submatrix from the matrix and assign.
     *
     * Removes all indexes provided in the vector from the matrix and returns it.
     *
     * \note Will not alter the matrix if given vector is empty.
     *
     * \param[in] _indexes The indexes to remove.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& operator/=(const std::vector<int> &_indexes) {
        for (int index : _indexes) {
            *this /= index;
        }
        
        return *this;
    }
    
    /**
     * \brief Get a submatrix from the matrix and assign.
     *
     * Removes index from the matrix and returns it.
     *
     * \note Will return a copy of the matrix if given vector is empty.
     * \note The submatrix will be rectangular.
     *
     * \param[in] _index The index to remove.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& operator/=(const int _index) {
        _chains[_index].nullify();
        _chainsStates.setOff(_index);
        
        return *this;
    }
    
    /**
     * \brief Nullifies a column from the matrix.
     *
     * Removes column of index `_index` whatever the `ChainTypeFlag` of the matrix. For column matrices, it just comes to the `\=` operator and for row matrices, it entails a traversal of the matrix.
     *
     * \param[in] _index The index to remove.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& del_column(int _index) {
        std::vector<int> tmp_id{_index} ;
        if (ChainTypeFlag == OSM::COLUMN)
        {
            (*this)/=tmp_id ;
        }
        else // OSM::ROW
        {
            for (int index : _chainsStates)
            {
                MatrixChain &tmp(_chains[index]) ;
                tmp/=tmp_id ;
                // If the row index has become empty: update
                if (tmp.isNull())
                    *this /= std::vector<int>({index}) ;
            }
        }
        
        return *this;
    }
    
    /**
     * \brief Nullifies a row from the matrix.
     *
     * Removes row of index `_index` whatever the `ChainTypeFlag` of the matrix. For row matrices, it just comes to the `\=` operator and for column matrices, it entails a traversal of the matrix.
     *
     * \param[in] _index The index to remove.
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& del_row(int _index)
    {
        std::vector<int> tmp_id{_index};
        if (ChainTypeFlag == OSM::ROW) {
            (*this) /= tmp_id;
        } else // OSM::COLUMN
        {
            for (int index : _chainsStates) {
                MatrixChain &tmp(_chains[index]);
                tmp /= tmp_id;
                // If the column index has become empty: update
                if (tmp.isNull())
                    *this /= std::vector<int>({index}) ;
            }
        }
        
        return *this;
    }
    
    /**
     * \brief Nullifies a coefficient from the matrix.
     *
     * Removes coefficient at row `i` and column `j` and update (if necessary) the status of corresponding _chains (null or not).
     *
     * \param[in] i Index of the row
     * \param[in] j Index of the column
     *
     * \return The modified matrix representing the result.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link std::vector \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& del_coef(int i, int j)
    {
        
        if (ChainTypeFlag == OSM::COLUMN) {
            std::vector<int> tmp_id({i}) ;
            MatrixChain &tmp(_chains[j]);
            tmp /= tmp_id;
            if (tmp.isNull())
                _chainsStates.setOff(j) ;
        } else // OSM::ROW
        {
            std::vector<int> tmp_id({j}) ;
            MatrixChain &tmp(_chains[i]);
            tmp /= tmp_id;
            if (tmp.isNull())
                _chainsStates.setOff(i) ;
        }
        
        return *this;
    }
    
    /**
     * \brief Iterator to the beginning of the _chains indices.
     *
     * \return The iterator to the beginning of the _chains indices.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * \see \link std::vector::iterator \endlink
     * \see \link std::vector::begin \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    inline Bitboard::iterator begin() const noexcept { return _chainsStates.begin(); }
    
    /**
     * \brief Iterator to the ending of the _chains indices.
     *
     * \return The iterator to the ending of the _chains indices.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * \see \link std::vector::iterator \endlink
     * \see \link std::vector::end \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    inline Bitboard::iterator end() const noexcept { return _chainsStates.end(); }
    
    /**
     * \brief Reverse iterator to the _chains indices.
     *
     * \return The reverse iterator to the _chains indices.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * \see \link std::vector::iterator \endlink
     * \see \link std::vector::begin \endlink
     *
     * \author Bac A.
     * \version 0.1.0
     * \date 08/04/2024
     */
    inline Bitboard::reverse_iterator reverse_begin() noexcept { return _chainsStates.reverse_begin(); }
    inline Bitboard::reverse_iterator reverse_begin(size_t index) noexcept { return _chainsStates.reverse_begin(index); }
    
    /**
     * \brief Reverse iterator to the ending of the _chains indices.
     *
     * \return The reverse iterator to the ending of the _chains indices.
     *
     * \see \link OSM::Sparse_matrix \endlink
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * \see \link std::vector::iterator \endlink
     * \see \link std::vector::end \endlink
     *
     * \author Bac A.
     * \version 0.1.0
     * \date 08/04/2024
     */
    inline Bitboard::reverse_iterator reverse_end() noexcept { return _chainsStates.reverse_end(); }
    
    
    /**
     * \brief Transpose a matrix.
     *
     * \return A new matrix where the chain type flag is changed.
     *
     * \see \link OSM::Sparse_matrix \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix<CoefficientType, COLUMN + ROW - ChainTypeFlag> transpose() {
        Sparse_matrix<CoefficientType, COLUMN + ROW - ChainTypeFlag> transposed(this->_size.second, this->_size.first);
        
        for (int index : this->_chainsStates) {
            transposed._chains[index] = this->_chains[index].transpose();
        }
        
        transposed._chainsStates = this->_chainsStates;
        
        return transposed;
    }
    
    /**
     * \brief Gets the matrix size.
     *
     * \return The matrix size as a Row/Column pair.
     *
     * \see \link OSM::Sparse_matrix \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 06/05/2024
     */
    std::pair<int, int> dimensions() const {
        return this->_size;
    }
};


/** @ingroup MatrixProdCol
 * @brief COLUMN x COLUMN -> COLUMN
 */
// COLUMN x COLUMN -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, COLUMN> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);
    
    // Perform col-col matrix multiplication with linear combination of columns.
    for (int index: second._chainsStates) {
        Chain<_CT, COLUMN> column(first._size.first);
        
        for (std::pair colRight: second._chains[index]) {
            if (first._chainsStates.isOn(colRight.first)) {
                column += colRight.second * first._chains[colRight.first];
            }
        }
        
        res[index] = column;
    }
    
    return res;
}

/** @ingroup MatrixProdCol
 * @brief ROW x COLUMN -> COLUMN
 */
// ROW x COLUMN -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, COLUMN> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int colRight: second._chainsStates) {
        Chain<_CT, COLUMN> column(first._size.first);
        
        for (int rowLeft: first._chainsStates) {
            _CT coef = first._chains[rowLeft] * second._chains[colRight];
            if (coef != 0) {
                column[rowLeft] = coef;
            }
        }
        
        res[colRight] = column;
    }
    
    return res;
}

/** @ingroup MatrixProdCol
 * @brief COLUMN x ROW -> COLUMN
 */
// COLUMN x ROW -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &first, const Sparse_matrix<_CT, ROW> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int colLeft: first._chainsStates) {
        res += first._chains[colLeft] * second._chains[colLeft];
    }
    
    return res;
}

/** @ingroup MatrixProdCol
 * @brief ROW x ROW -> COLUMN
 */
// ROW x ROW -> COLUMN
template <typename _CT>
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &first, const Sparse_matrix<_CT, ROW> &second) {
    Sparse_matrix<_CT, COLUMN> res(first._size.first, second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int i = 0 ; i < second._size.second ; i++) {
        Chain<_CT, COLUMN> column(first._size.first);
        
        for (int rowLeft: first._chainsStates) {
            _CT coef = first._chains[rowLeft] * get_column(second, i);
            if (coef != 0) {
                column[rowLeft] = coef;
            }
        }
        
        if (!column.isNull()) {
            res[i] = column;
        }
    }
    
    return res;
}

/**
 * \brief Perform multiplication between a column-based matrix and a column-based chain.
 *
 * Generate a column-based chain from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The matrix.
 * \param[in] _second The chain.
 *
 * \return The result of the matrix multiplication, column-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Chain<_CT, COLUMN> &_second)
{
    // Perform col-col matrix multiplication with linear combination of columns.
    Chain<_CT, COLUMN> column(_first._size.first);
    
    for (typename Chain<_CT, COLUMN>::const_iterator it = _second.begin(); it != _second.end(); ++it)
    {
        column += it->second * _first._chains[it->first];
    }
    
    return column;
}

/**
 * \brief Perform multiplication between a row-based matrix and a column-based chain.
 *
 * Generate a column-based chain from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The matrix.
 * \param[in] _second The chain.
 *
 * \return The result of the matrix multiplication, column-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Chain<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Chain<_CT, COLUMN> &_second)
{
    // Perform row-col matrix multiplication with dots
    Chain<_CT, COLUMN> column(_first._size.first);
    
    for (int index : _first._chainsStates)
    {
        _CT tmp(_first[index] * _second) ;
        if (tmp != 0)
            column[index] = tmp ;
    }
    return column;
}

/**
 * \brief Perform multiplication between a row-based chain and a row-based matrix.
 *
 * Generate a row-based chain from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The matrix.
 * \param[in] _second The chain.
 *
 * \return The result of the matrix multiplication, row-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Chain<_CT, ROW> operator*(const Chain<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second)
{
    // Perform row-row matrix multiplication with linear combination of rows.
    Chain<_CT, ROW> row(_second._size.second);
    
    for (typename Chain<_CT, ROW>::const_iterator it = _first.begin(); it != _first.end(); ++it)
    {
        row += it->second * _second._chains[it->first];
    }
    
    return row;
}

/**
 * \brief Perform multiplication between a row-based chain and a column-based matrix.
 *
 * Generate a row-based chain from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The matrix.
 * \param[in] _second The chain.
 *
 * \return The result of the matrix multiplication, row-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Chain<_CT, ROW> operator*(const Chain<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second)
{
    // Perform row-col matrix multiplication with dots
    Chain<_CT, ROW> row(_second._size.second);
    
    for (int index : _second._chainsStates)
    {
        _CT tmp(_first * _second[index]) ;
        if (tmp != 0)
            row[index] = tmp ;
    }
    return row;
}

/**
 * \brief Perform matrix multiplication between two _chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The first matrix.
 * \param[in] _second The second matrix.
 *
 * \return The result of the matrix multiplication, column-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, COLUMN> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    for (int i = 0 ; i < _first._size.first ; i++) {
        Chain<_CT, ROW> row(_second._size.second);
        
        for (int colRight: _second._chainsStates) {
            _CT coef = get_row(_first, i) * _second._chains[colRight];
            if (coef != 0) {
                row[colRight] = coef;
            }
        }
        
        if (!row.isNull()) {
            res[i] = row;
        }
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two _chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The first matrix.
 * \param[in] _second The second matrix.
 *
 * \return The result of the matrix multiplication, column-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int rowLeft: _first._chainsStates) {
        Chain<_CT, ROW> row(_second._size.second);
        
        for (int colRight: _second._chainsStates) {
            _CT coef = _first._chains[rowLeft] * _second._chains[colRight];
            if (coef != 0) {
                row[colRight] = coef;
            }
        }
        
        res[rowLeft] = row;
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two _chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The first matrix.
 * \param[in] _second The second matrix.
 *
 * \return The result of the matrix multiplication, column-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, ROW> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int colLeft: _first._chainsStates) {
        res += _first._chains[colLeft] % _second._chains[colLeft];
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two _chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _first The first matrix.
 * \param[in] _second The second matrix.
 *
 * \return The result of the matrix multiplication, column-based.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW> operator%(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second) {
    Sparse_matrix<_CT, ROW> res(_first._size.first, _second._size.second);
    
    // Perform row-row matrix multiplication with linear combination of rows.
    for (int index: _first._chainsStates) {
        Chain<_CT, ROW> row(_second._size.second);
        
        for (std::pair colRight: _first._chains[index]) {
            if (_first._chainsStates.isOn(colRight.first)) {
                row += colRight.second * _second._chains[colRight.first];
            }
        }
        
        res[index] = row;
    }
    
    return res;
}

/**
 * \brief Add a matrix and assign.
 *
 * Adds each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _matrix The matrix to reassign.
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator+=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other._size.second ; index++) {
        Chain column = get_column(_other, index);
        if (!column.isNull()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] += get_column(_other, index);
            
            if (_matrix._chains[index].isNull()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

/**
 * \brief Add a matrix and assign.
 *
 * Adds each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _matrix The matrix to reassign.
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator+=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other._size.first ; index++) {
        Chain row = get_row(_other, index);
        if (!row.isNull()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] += get_row(_other, index);
            
            if (_matrix._chains[index].isNull()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

/**
 * \brief Substract a matrix and assign.
 *
 * Substracts each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _matrix The matrix to reassign.
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator-=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other._size.second ; index++) {
        Chain column = get_column(_other, index);
        if (!column.isNull()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] -= get_column(_other, index);
            
            if (_matrix._chains[index].isNull()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

/**
 * \brief Substract a matrix and assign.
 *
 * Substracts each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _matrix The matrix to reassign.
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator-=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    if (_matrix._size != _other._size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other._size.second ; index++) {
        Chain row = get_row(_other, index);
        if (!row.isNull()) {
            _matrix._chainsStates |= index;
            _matrix._chains[index] -= get_row(_other, index);
            
            if (_matrix._chains[index].isNull()) {
                _matrix._chainsStates.setOff(index);
            }
        }
    }
    
    return _matrix;
}

/**
 * \brief Multiply a matrix and assign.
 *
 * Multiply each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Sparse_matrix \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    _matrix = _matrix * _other;
    return _matrix;
}

/**
 * \brief Multiply a matrix and assign.
 *
 * Multiply each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    _matrix = _matrix % _other;
    return _matrix;
}

/**
 * \brief Multiply a matrix and assign.
 *
 * Multiply each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, COLUMN>& operator*=(Sparse_matrix<_CT, COLUMN> &_matrix, const Sparse_matrix<_CT, ROW> &_other) {
    _matrix = _matrix * _other;
    return _matrix;
}

/**
 * \brief Multiply a matrix and assign.
 *
 * Multiply each coefficient of the matrix together.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two _chains are not the same coefficient type.
 *
 * \param[in] _other The other matrix.
 *
 * \return The modified matrix representing the result.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.3.0
 * \date 31/04/2024
 */
template <typename _CT>
Sparse_matrix<_CT, ROW>& operator*=(Sparse_matrix<_CT, ROW> &_matrix, const Sparse_matrix<_CT, COLUMN> &_other) {
    _matrix = _matrix % _other;
    return _matrix;
}

/**
 * \brief Get a column from the matrix, even if matrix is a row-chain matrix.
 *
 * \note For column-chain matrixes, it is equivalent to operator[].
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 *
 * \return The column stored at given index.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
Chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, COLUMN> &_matrix, const int _index) {
    return _matrix._chains[_index];
}

/**
 * \brief Get a column from the matrix, even if matrix is a row-chain matrix.
 *
 * \note For column-chain matrixes, it is equivalent to operator[].
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 *
 * \return The column stored at given index.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
Chain<_CT, COLUMN> get_column(const Sparse_matrix<_CT, ROW> &_matrix, const int _index) {
    Chain<_CT, COLUMN> column(_matrix._size.first);
    if (_matrix._size.first > 0)
    {
        for (int i : _matrix._chainsStates) {
            if (!_matrix._chains[i].isNull(_index)) {
                column[i] = _matrix._chains[i][_index];
            }
        }
    }
    
    return column;
}

/**
 * \brief Get a row from the matrix, even if matrix is a row-chain matrix.
 *
 * \note For row-chain matrixes, it is equivalent to operator[].
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 *
 * \return The row stored at given index.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
Chain<_CT, ROW> get_row(const Sparse_matrix<_CT, COLUMN> &_matrix, const int _index) {
    Chain<_CT, ROW> row(_matrix._size.second);
    if (_matrix._size.second > 0)
    {
        for (int i : _matrix._chainsStates) {
            if (!_matrix._chains[i].isNull(_index)) {
                row[i] = _matrix._chains[i][_index];
            }
        }
    }
    
    return row;
}

/**
 * \brief Get a row from the matrix, even if matrix is a row-chain matrix.
 *
 * \note For row-chain matrixes, it is equivalent to operator[].
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 *
 * \return The row stored at given index.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
Chain<_CT, ROW> get_row(const Sparse_matrix<_CT, ROW> &_matrix, const int _index) {
    return _matrix._chains[_index];
}

/**
 * \brief Get a const reference over a column from a column matrix
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 *
 * \return A const reference over the column stored at given index.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
const Chain<_CT, COLUMN> & cget_column(const Sparse_matrix<_CT, COLUMN> &_matrix, const int _index)
{
    if (_index >= _matrix._size.second) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(_matrix._size.second) + ".");
    }
    return _matrix._chains[_index];
}

/**
 * \brief Get a const reference over a row from a row matrix
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 *
 * \return A const reference over the row stored at given index.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
const Chain<_CT, ROW> & cget_row(const Sparse_matrix<_CT, ROW> &_matrix, const int _index)
{
    if (_index >= _matrix._size.first) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(_matrix._size.first) + ".");
    }
    return _matrix._chains[_index];
}

/**
 * \brief Set a column from the matrix, even if matrix is a row-chain matrix.
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 * \param[in] _chain The new column.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
void set_column(Sparse_matrix<_CT, COLUMN> &_matrix, const int _index, const Chain<_CT, COLUMN> &_chain) {
    if(_matrix.dimensions().first != _chain.dimension())
        throw std::runtime_error("set_column dimension error") ;
    _matrix[_index] = _chain;
    if (_chain.isNull())
        _matrix._chainsStates.setOff(_index) ;
}

/**
 * \brief Set a column from the matrix, even if matrix is a row-chain matrix.
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 * \param[in] _chain The new column.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
void set_column(Sparse_matrix<_CT, ROW> &_matrix, const int _index, const Chain<_CT, COLUMN> &_chain) {
    if(_matrix.dimensions().first != _chain.dimension())
        throw std::runtime_error("set_column dimension error") ;
    for (int i = 0 ; i < _matrix._size.first ; i++) {
        if (!_matrix._chains[i].isNull(_index) && _chain.isNull(i)) {
            _matrix._chains[i] /= _index;
            
            if (_matrix._chains[i].isNull()) {
                _matrix._chainsStates.setOff(i);
            }
        }
        
        if (!_chain.isNull(i)) {
            _matrix._chainsStates |= i;
            _matrix._chains[i][_index] = _chain[i];
        }
    }
}

/**
 * \brief Set a row from the matrix, even if matrix is a row-chain matrix.
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 * \param[in] _chain The new row.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
void set_row(Sparse_matrix<_CT, COLUMN> &_matrix, const int _index, const Chain<_CT, ROW> &_chain) {
    if(_matrix.dimensions().second != _chain.dimension())
        throw("set_column dimension error") ;
    for (int i = 0 ; i < _matrix._size.second ; i++) {
        if (!_matrix._chains[i].isNull(_index) && _chain.isNull(i)) {
            _matrix._chains[i] /= _index;
            
            if (_matrix._chains[i].isNull()) {
                _matrix._chainsStates.setOff(i);
            }
        }
        
        if (!_chain.isNull(i)) {
            _matrix._chainsStates |= i;
            _matrix._chains[i][_index] = _chain[i];
        }
    }
}

/**
 * \brief Set a row from the matrix, even if matrix is a row-chain matrix.
 *
 * \warning The matrix will perform boundary check.
 *
 * \param[in] _index The coefficient index.
 * \param[in] _chain The new row.
 *
 * \see \link OSM::Chain \endlink
 *
 * \author Fedyna K.
 * \version 0.1.0
 * \date 17/04/2024
 */
template <typename _CT>
void set_row(Sparse_matrix<_CT, ROW> &_matrix, const int _index, const Chain<_CT, ROW> &_chain) {
    if(_matrix.dimensions().second != _chain.dimension())
        throw("set_column dimension error") ;
    _matrix[_index] = _chain;
    if (_chain.isNull())
        _matrix._chainsStates.setOff(_index) ;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SPARSE_MATRIX_H
