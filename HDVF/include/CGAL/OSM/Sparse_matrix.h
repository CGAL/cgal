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
    std::vector<Chain<CoefficientType, ChainTypeFlag>> chains;
    
    /** \brief A bitboard containing state of each columns. */
    Bitboard chainsStates;
    
    /** \brief The matrix size as a (row, column) pair. */
    std::pair<int, int> size;
    
    /**
     * \brief Get a chain from a matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param[in] _index The coefficient index.
     *
     * \return The reference to the chain stored at given index.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    MatrixChain& operator[](const int _index) {
        if (ChainTypeFlag == COLUMN && _index >= size.second) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(size.second) + ".");
        }
        if (ChainTypeFlag == ROW && _index >= size.first) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(size.first) + ".");
        }
        
        chainsStates |= _index;
        return chains[_index];
    }
    
public:
    /**
     * \brief Create new Sparse_matrix object.
     *
     * Default constructor, initialize an empty Matrix as Z integers column chains.
     * The default matrix size is 128x128.
     *
     * \tparam CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix() {
        chains = std::vector<MatrixChain>(0);
        chainsStates = Bitboard(0);
        size = {0, 0};
    }
    
    /**
     * \brief Create new Sparse_matrix object.
     *
     * Constructor with size, initialize an empty Matrix as Z integers column chains.
     *
     * \tparam CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     * \param[in] _rowCount The number of rows to preallocate.
     * \param[in] _columnCount The number of columns to preallocate.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix(const int _rowCount, const int _columnCount) {
        int mainSize = ChainTypeFlag == COLUMN ? _columnCount : _rowCount;
        int secondarySize = ChainTypeFlag == COLUMN ? _rowCount : _columnCount;
        
        chains = std::vector<Chain<CoefficientType, ChainTypeFlag>>(mainSize);
        for (int i = 0 ; i < mainSize ; i++) {
            chains[i] = Chain<CoefficientType, ChainTypeFlag>(secondarySize);
        }
        
        chainsStates = Bitboard(mainSize);
        size = {_rowCount, _columnCount};
    }
    
    /**
     * \brief Create new Sparse_matrix from other Sparse_matrix object.
     *
     * Copy constructor, initialize a Sparse_matrix (not necessarly based on the same type of matrix).
     * If types are different, the constructor performs conversion.
     *
     * \tparam CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     * \param[in] _otherToCopy The matrix we want to copy.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int CTF>
    Sparse_matrix(const Sparse_matrix<CoefficientType,CTF> &_otherToCopy) {
        //        *this = _otherToCopy;
        if (ChainTypeFlag == CTF)
        {
            chainsStates = _otherToCopy.chainsStates;
            size = _otherToCopy.size;
            // Copy of chains as such
            chains.resize(_otherToCopy.chains.size()) ;
            for (int i = 0; i<_otherToCopy.chains.size(); ++i)
            {
                const Chain<CoefficientType, CTF>& tmp(_otherToCopy.chains.at(i)) ;
                Chain<CoefficientType,ChainTypeFlag> res(tmp.dimension()) ;
                for (typename Chain<CoefficientType, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
                {
                    res[it->first] = it->second ;
                }
                chains[i] = res ;
            }
        }
        else // Copy coefficient by coefficient
        {
            size = _otherToCopy.size;
            int vec_size = (ChainTypeFlag == OSM::COLUMN)?size.second:size.first ;
            int chain_size = (ChainTypeFlag == OSM::COLUMN)?size.first:size.second ;
            chains.resize(vec_size) ;
            chainsStates = OSM::Bitboard(vec_size) ;
            for (int i=0; i<vec_size; ++i)
            {
                chains.at(i) = Chain<CoefficientType,ChainTypeFlag>(chain_size) ;
            }
            
            for (int i = 0; i<_otherToCopy.chains.size(); ++i)
            {
                const Chain<CoefficientType, CTF>& tmp(_otherToCopy.chains.at(i)) ;
                for (typename Chain<CoefficientType, CTF>::const_iterator it = tmp.cbegin(); it != tmp.cend(); ++it)
                {
                    chains[it->first][i] = it->second ;
                    chainsStates.setOn(it->first) ;
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
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Sparse_matrix& operator=(const Sparse_matrix& _otherToCopy)
    {
        chainsStates = _otherToCopy.chainsStates;
        size = _otherToCopy.size;
        chains = _otherToCopy.chains;
        
        return *this;
    }
    
    
    /**
     * \brief Clean a Sparse_matrix (set everything to zero)
     *
     * Empty all structures of the sparse matrix.
     *
     * \warning Will raise an error if the other matrix is not the same chain type.
     *
     * \see \link OSM::Sparse_matrix \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
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
     * \brief Displays a Sparse_matrix in the output stream.
     *
     * \param[in] _stream The output stream.
     * \param[in] _matrix The matrix to display.
     *
     * \return A reference to the modified stream.
     *
     * \see \link OSM::Matrix \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    friend std::ostream& operator<<(std::ostream &_stream, const Sparse_matrix &_matrix) {
#ifndef DEBUG
        bool empty = true;
        
        _stream << "[";
        for (int index : _matrix.chainsStates) {
            empty = false;
            _stream << index << ": " << _matrix.chains[index] << ", ";
        }
        
        if (!empty) {
            _stream << "\b\b";
        }
        _stream << "]" << std::endl ;
#else
        for (int i=0; i<_matrix.size.first; ++i)
        {
            for (int j=0; j<_matrix.size.second; ++j)
            {
                
                //CoefficientType tmp = _matrix[j][i] ;  // get_coef(i,j) ;
                CoefficientType tmp = _matrix.get_coef(i,j) ;
                
                if (tmp == 0)
                    _stream << ".\t" ;
                else
                    _stream << tmp << "\t" ;
            }
            _stream << std::endl ;
        }
#endif
        
        return _stream;
    }
    
    
    /**
     * \brief Adds two matrices together.
     *
     * Adds each coefficient of the matrix together.
     *
     * \pre The matrixs have the same chain type.
     *
     * \warning Will raise an error if the other matrix is not the same chain type.
     *
     * \param[in] _first The first matrix.
     * \param[in] _second The second matrix.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int _CTF>
    friend Sparse_matrix operator+(const Sparse_matrix &_first, const Sparse_matrix<CoefficientType, _CTF> &_second) {
        Sparse_matrix newMatrix(_first);
        newMatrix += _second;
        
        return newMatrix;
    }
    
    /**
     * \brief Substracts two matrices together.
     *
     * Substracts each coefficient of the matrix together.
     *
     * \pre The matrixs have the same chain type.
     *
     * \warning Will raise an error if the other matrix is not the same chain type.
     *
     * \param[in] _first The first matrix.
     * \param[in] _second The second matrix.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int _CTF>
    friend Sparse_matrix operator-(const Sparse_matrix &_first, const Sparse_matrix<CoefficientType, _CTF> &_second) {
        Sparse_matrix newMatrix = _first;
        newMatrix -= _second;
        
        return newMatrix;
    }
    
    /**
     * \brief Apply factor on each coefficients.
     *
     * \warning Will raise an error if _lambda is 0.
     *
     * \param[in] _lambda The factor to apply.
     * \param[in] _matrix The matrix.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix operator*(const _CT& _lambda, const Sparse_matrix<_CT, _CTF> &_matrix) {
        Sparse_matrix newMatrix = _matrix;
        newMatrix *= _lambda;
        
        return newMatrix;
    }
    
    /**
     * \brief Apply factor on each coefficients.
     *
     * \warning Will raise an error if _lambda is 0.
     *
     * \param[in] _matrix The matrix.
     * \param[in] _lambda The factor to apply.
     *
     * \return A new matrix representing the result.
     *
     * \see \link OSM::Chain \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT, int _CTF>
    friend Sparse_matrix operator*(const Sparse_matrix<_CT, _CTF> &_matrix, const _CT& _lambda) {
        Sparse_matrix newMatrix = _matrix;
        newMatrix *= _lambda;
        
        return newMatrix;
    }
    
    /**
     * \brief Perform matrix multiplication between two matrices.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, COLUMN> &_second);
    
    /**
     * \brief Perform matrix multiplication between two matrices.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second);
    
    /**
     * \brief Perform matrix multiplication between two matrices.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, ROW> &_second);
    
    /**
     * \brief Perform matrix multiplication between two matrices.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
    friend Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second);
    
    /**
     * \brief Perform multiplication between a column-based matrix and a column-based chain.
     *
     * Generate a column-based chain from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \brief Perform matrix multiplication between two chains.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \brief Perform matrix multiplication between two chains.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre The matrix have the same coefficent type.
     *
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
        if (this->size != _other.size) {
            throw std::runtime_error("Matrices must be the same size.");
        }
        
        for (int index: _other.chainsStates) {
            this->chainsStates |= index;
            this->chains[index] += _other.chains[index];
            
            if (this->chains[index].isNull()) {
                this->chainsStates.setOff(index);
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
        if (this->size != _other.size) {
            throw std::runtime_error("Matrices must be the same size.");
        }
        
        for (int index: _other.chainsStates) {
            this->chainsStates |= index;
            this->chains[index] -= _other.chains[index];
            
            if (this->chains[index].isNull()) {
                this->chainsStates.setOff(index);
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
        
        for (int index: this->chainsStates) {
            this->chains[index] *= _lambda;
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
        Sparse_matrix res(_matrix.size.first, _matrix.size.second) ;
        
        for (int index: _matrix.chainsStates)
        {
            const MatrixChain& tmp_chain(_matrix.chains[index]) ;
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
     * \warning Will raise an error if the two chains are not the same coefficient type.
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
        if (ChainTypeFlag == COLUMN && _index >= size.second) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(size.second) + ".");
        }
        if (ChainTypeFlag == ROW && _index >= size.first) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(size.first) + ".");
        }
        
        return chains[_index];
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
        if (_i >= size.first) {
            throw std::runtime_error("Provided _i index should be less than " + std::to_string(size.first) + ".");
        }
        if (_j >= size.second) {
            throw std::runtime_error("Provided _j index should be less than " + std::to_string(size.second) + ".");
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
        if (_i >= size.first) {
            //std::cout << "Provided _i index should be less than " << size.first << "." << std::endl;
            throw std::runtime_error("Provided _i index should be less than " + std::to_string(size.first) + ".");
        }
        if (_j >= size.second) {
            //std::cout << "Provided _j index should be less than " << size.second << "." << std::endl;
            throw std::runtime_error("Provided _j index should be less than " + std::to_string(size.second) + ".");
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
        chains[_index].nullify();
        chainsStates.setOff(_index);
        
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
            for (int index : chainsStates)
            {
                MatrixChain &tmp(chains[index]) ;
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
            for (int index : chainsStates) {
                MatrixChain &tmp(chains[index]);
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
     * Removes coefficient at row `i` and column `j` and update (if necessary) the status of corresponding chains (null or not).
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
            MatrixChain &tmp(chains[j]);
            tmp /= tmp_id;
            if (tmp.isNull())
                chainsStates.setOff(j) ;
        } else // OSM::ROW
        {
            std::vector<int> tmp_id({j}) ;
            MatrixChain &tmp(chains[i]);
            tmp /= tmp_id;
            if (tmp.isNull())
                chainsStates.setOff(i) ;
        }
        
        return *this;
    }
    
    /**
     * \brief Iterator to the beginning of the chains indices.
     *
     * \return The iterator to the beginning of the chains indices.
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
    inline Bitboard::iterator begin() const noexcept { return chainsStates.begin(); }
    
    /**
     * \brief Iterator to the ending of the chains indices.
     *
     * \return The iterator to the ending of the chains indices.
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
    inline Bitboard::iterator end() const noexcept { return chainsStates.end(); }
    
    /**
     * \brief Reverse iterator to the chains indices.
     *
     * \return The reverse iterator to the chains indices.
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
    inline Bitboard::reverse_iterator reverse_begin() noexcept { return chainsStates.reverse_begin(); }
    inline Bitboard::reverse_iterator reverse_begin(size_t index) noexcept { return chainsStates.reverse_begin(index); }
    
    /**
     * \brief Reverse iterator to the ending of the chains indices.
     *
     * \return The reverse iterator to the ending of the chains indices.
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
    inline Bitboard::reverse_iterator reverse_end() noexcept { return chainsStates.reverse_end(); }
    
    
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
        Sparse_matrix<CoefficientType, COLUMN + ROW - ChainTypeFlag> transposed(this->size.second, this->size.first);
        
        for (int index : this->chainsStates) {
            transposed.chains[index] = this->chains[index].transpose();
        }
        
        transposed.chainsStates = this->chainsStates;
        
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
        return this->size;
    }
};



/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a column-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, COLUMN> &_second) {
    Sparse_matrix<_CT, COLUMN> res(_first.size.first, _second.size.second);
    
    // Perform col-col matrix multiplication with linear combination of columns.
    for (int index: _second.chainsStates) {
        Chain<_CT, COLUMN> column(_first.size.first);
        
        for (std::pair colRight: _second.chains[index]) {
            if (_first.chainsStates.isOn(colRight.first)) {
                column += colRight.second * _first.chains[colRight.first];
            }
        }
        
        res[index] = column;
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a column-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, COLUMN> &_second) {
    Sparse_matrix<_CT, COLUMN> res(_first.size.first, _second.size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int colRight: _second.chainsStates) {
        Chain<_CT, COLUMN> column(_first.size.first);
        
        for (int rowLeft: _first.chainsStates) {
            _CT coef = _first.chains[rowLeft] * _second.chains[colRight];
            if (coef != 0) {
                column[rowLeft] = coef;
            }
        }
        
        res[colRight] = column;
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a column-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, COLUMN> &_first, const Sparse_matrix<_CT, ROW> &_second) {
    Sparse_matrix<_CT, COLUMN> res(_first.size.first, _second.size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int colLeft: _first.chainsStates) {
        res += _first.chains[colLeft] * _second.chains[colLeft];
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a column-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
Sparse_matrix<_CT, COLUMN> operator*(const Sparse_matrix<_CT, ROW> &_first, const Sparse_matrix<_CT, ROW> &_second) {
    Sparse_matrix<_CT, COLUMN> res(_first.size.first, _second.size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int i = 0 ; i < _second.size.second ; i++) {
        Chain<_CT, COLUMN> column(_first.size.first);
        
        for (int rowLeft: _first.chainsStates) {
            _CT coef = _first.chains[rowLeft] * get_column(_second, i);
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Chain<_CT, COLUMN> column(_first.size.first);
    
    for (typename Chain<_CT, COLUMN>::const_iterator it = _second.begin(); it != _second.end(); ++it)
    {
        column += it->second * _first.chains[it->first];
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Chain<_CT, COLUMN> column(_first.size.first);
    
    for (int index : _first.chainsStates)
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Chain<_CT, ROW> row(_second.size.second);
    
    for (typename Chain<_CT, ROW>::const_iterator it = _first.begin(); it != _first.end(); ++it)
    {
        row += it->second * _second.chains[it->first];
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Chain<_CT, ROW> row(_second.size.second);
    
    for (int index : _second.chainsStates)
    {
        _CT tmp(_first * _second[index]) ;
        if (tmp != 0)
            row[index] = tmp ;
    }
    return row;
}

/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Sparse_matrix<_CT, ROW> res(_first.size.first, _second.size.second);
    
    for (int i = 0 ; i < _first.size.first ; i++) {
        Chain<_CT, ROW> row(_second.size.second);
        
        for (int colRight: _second.chainsStates) {
            _CT coef = get_row(_first, i) * _second.chains[colRight];
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
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Sparse_matrix<_CT, ROW> res(_first.size.first, _second.size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int rowLeft: _first.chainsStates) {
        Chain<_CT, ROW> row(_second.size.second);
        
        for (int colRight: _second.chainsStates) {
            _CT coef = _first.chains[rowLeft] * _second.chains[colRight];
            if (coef != 0) {
                row[colRight] = coef;
            }
        }
        
        res[rowLeft] = row;
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Sparse_matrix<_CT, ROW> res(_first.size.first, _second.size.second);
    
    // Perform row-col matrix multiplication with dot products.
    for (int colLeft: _first.chainsStates) {
        res += _first.chains[colLeft] % _second.chains[colLeft];
    }
    
    return res;
}

/**
 * \brief Perform matrix multiplication between two chains.
 *
 * Generate a row-based matrix from the matrix multiplication and return it.
 *
 * \pre The matrix have the same coefficent type.
 *
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    Sparse_matrix<_CT, ROW> res(_first.size.first, _second.size.second);
    
    // Perform row-row matrix multiplication with linear combination of rows.
    for (int index: _first.chainsStates) {
        Chain<_CT, ROW> row(_second.size.second);
        
        for (std::pair colRight: _first.chains[index]) {
            if (_first.chainsStates.isOn(colRight.first)) {
                row += colRight.second * _second.chains[colRight.first];
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    if (_matrix.size != _other.size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other.size.second ; index++) {
        Chain column = get_column(_other, index);
        if (!column.isNull()) {
            _matrix.chainsStates |= index;
            _matrix.chains[index] += get_column(_other, index);
            
            if (_matrix.chains[index].isNull()) {
                _matrix.chainsStates.setOff(index);
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    if (_matrix.size != _other.size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other.size.first ; index++) {
        Chain row = get_row(_other, index);
        if (!row.isNull()) {
            _matrix.chainsStates |= index;
            _matrix.chains[index] += get_row(_other, index);
            
            if (_matrix.chains[index].isNull()) {
                _matrix.chainsStates.setOff(index);
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    if (_matrix.size != _other.size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other.size.second ; index++) {
        Chain column = get_column(_other, index);
        if (!column.isNull()) {
            _matrix.chainsStates |= index;
            _matrix.chains[index] -= get_column(_other, index);
            
            if (_matrix.chains[index].isNull()) {
                _matrix.chainsStates.setOff(index);
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    if (_matrix.size != _other.size) {
        throw std::runtime_error("Matrices must be the same size.");
    }
    
    for (int index = 0 ; index < _other.size.second ; index++) {
        Chain row = get_row(_other, index);
        if (!row.isNull()) {
            _matrix.chainsStates |= index;
            _matrix.chains[index] -= get_row(_other, index);
            
            if (_matrix.chains[index].isNull()) {
                _matrix.chainsStates.setOff(index);
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
 * \warning Will raise an error if the two chains are not the same coefficient type.
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
    return _matrix.chains[_index];
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
    Chain<_CT, COLUMN> column(_matrix.size.first);
    if (_matrix.size.first > 0)
    {
        for (int i : _matrix.chainsStates) {
            if (!_matrix.chains[i].isNull(_index)) {
                column[i] = _matrix.chains[i][_index];
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
    Chain<_CT, ROW> row(_matrix.size.second);
    if (_matrix.size.second > 0)
    {
        for (int i : _matrix.chainsStates) {
            if (!_matrix.chains[i].isNull(_index)) {
                row[i] = _matrix.chains[i][_index];
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
    return _matrix.chains[_index];
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
    if (_index >= _matrix.size.second) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(_matrix.size.second) + ".");
    }
    return _matrix.chains[_index];
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
    if (_index >= _matrix.size.first) {
        throw std::runtime_error("Provided index should be less than " + std::to_string(_matrix.size.first) + ".");
    }
    return _matrix.chains[_index];
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
        _matrix.chainsStates.setOff(_index) ;
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
    for (int i = 0 ; i < _matrix.size.first ; i++) {
        if (!_matrix.chains[i].isNull(_index) && _chain.isNull(i)) {
            _matrix.chains[i] /= _index;
            
            if (_matrix.chains[i].isNull()) {
                _matrix.chainsStates.setOff(i);
            }
        }
        
        if (!_chain.isNull(i)) {
            _matrix.chainsStates |= i;
            _matrix.chains[i][_index] = _chain[i];
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
    for (int i = 0 ; i < _matrix.size.second ; i++) {
        if (!_matrix.chains[i].isNull(_index) && _chain.isNull(i)) {
            _matrix.chains[i] /= _index;
            
            if (_matrix.chains[i].isNull()) {
                _matrix.chainsStates.setOff(i);
            }
        }
        
        if (!_chain.isNull(i)) {
            _matrix.chainsStates |= i;
            _matrix.chains[i][_index] = _chain[i];
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
        _matrix.chainsStates.setOff(_index) ;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_SPARSE_MATRIX_H
