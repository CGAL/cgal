/**
 * \file Chain.hpp
 * \brief Namespace file for describing library.
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 * 
 * Define everything for the Chain class
 */

#ifndef __OSM_CHAIN__
#define __OSM_CHAIN__


#include "__base.hpp"
#include "SparseMatrix.hpp"
#include <unordered_map>
#include <vector>
#include <iterator>
#include <iostream>

namespace CGAL {
namespace OSM {

/**
 * \class Chain
 * \brief Allow to create row/column sparse vectors.
 * 
 * The Chain class contains all algebraic functions that are related to vectors.
 * 
 * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam _ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
 * 
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CoefficientType, int _ChainTypeFlag>
class Chain {
    
public:
    typedef std::pair<int, _CoefficientType> pair;
    typedef typename std::unordered_map<int, _CoefficientType>::iterator iterator;
    typedef typename std::unordered_map<int, _CoefficientType>::const_iterator const_iterator;
    
    // Add these lines to allow the SparseMatrix class to access other templated
    // SparsedMatrix private members.
    template <typename _CT, int _CTF>
    friend class Chain;
    template <typename _CT, int _CTF>
    friend class SparseMatrix;
protected:
    /** \brief The chain inner representation and storage of data. */
    std::unordered_map<int, _CoefficientType> chainData;
    
    /** \brief The chain coefficient type. */
    typedef _CoefficientType coefficientType;
    
    /** \brief The chain type flag. */
    int chainTypeFlag = _ChainTypeFlag;
    
    /** \brief The chain boundary. */
    int upperBound;
    
public:
    /**
     * \brief Create new Chain for SparseMatrix object.
     * 
     * Default constructor, initialize an empty Chain as a Z integers column chain.
     * The default chain size is 128.
     * 
     * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam _ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     * 
     * \see \link OSM::ZCoefficient \endlink
     * \see \link OSM::COLUMN \endlink
     * \see \link OSM::ROW \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain() {
        upperBound = 0;
        chainData = std::unordered_map<int, _CoefficientType>();
    }
    
    /**
     * \brief Create new Chain for SparseMatrix object.
     * 
     * Constructor with size, initialize an empty Chain as a Z integers column chain with boundary check.
     * 
     * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam _ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     * \param[in] _chainSize The upper bound of the Chain.
     * 
     * \see \link OSM::ZCoefficient \endlink
     * \see \link OSM::COLUMN \endlink
     * \see \link OSM::ROW \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain(const int _chainSize) {
        upperBound = _chainSize;
        chainData = std::unordered_map<int, _CoefficientType>();
    }
    
    /**
     * \brief Size of the chain
     *
     * \return Size of the chain (upperBound)
     *
     * \see \link OSM::ZCoefficient \endlink
     * \see \link OSM::COLUMN \endlink
     * \see \link OSM::ROW \endlink
     *
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    size_t dimension() const { return upperBound ; }
    
    /**
     * \brief Create new Chain for SparseMatrix object.
     * 
     * Copy constructor, initialize a Chain based on the same type of chain.
     * The resulting chain will be a copy of the passed chain.
     * 
     * \pre The chains have the same coefficent type.
     * 
     * \warning Will raise an error if the other chain is not the same coefficient type.
     * 
     * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam _ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     * \param[in] _otherToCopy The chain we want to copy.
     * 
     * \see \link OSM::ZCoefficient \endlink
     * \see \link OSM::COLUMN \endlink
     * \see \link OSM::ROW \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain(const Chain &_otherToCopy) {
        upperBound = _otherToCopy.upperBound;
        chainData = _otherToCopy.chainData;
        //        *this = _otherToCopy;
    }
    
    /**
     * \brief Assign to other chain.
     * 
     * Assign to other chain coefficient-wise, equivalent to copying it.
     * 
     * \pre The chains have the same coefficent type.
     * 
     * \warning Will raise an error if the other chain is not the same coefficient type.
     * 
     * \param[in] _otherToCopy The chain we want to copy.
     * 
     * \return The reference to the modified chain.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain& operator=(const Chain &_otherToCopy) {
        upperBound = _otherToCopy.upperBound;
        chainData = _otherToCopy.chainData;
        
        return *this;
    }
    
    /**
     * \brief Displays a Chain in the output stream.
     * 
     * \param[in] _stream The output stream.
     * \param[in] _chain The chain to display.
     * 
     * \return A reference to the modified stream.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    friend std::ostream& operator<<(std::ostream &_stream, const Chain &_chain) {
        _stream << "[";
        for (const_iterator i = _chain.chainData.begin() ; i != _chain.chainData.end() ; ++i) {
            _stream << i->first << ": " << i->second << ", ";
        }
        
        if (_chain.chainData.size() > 0) {
            _stream << "\b\b";
        }
        _stream << "]";
        
        return _stream;
    }
    
    /**
     * \brief Adds two chains together.
     * 
     * Adds each coefficient of the chain together.
     * 
     * \pre The chains have the same coefficent type and the same type flag.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * \warning Will raise an error if the two chains don't have the same type flag.
     * 
     * \param[in] _first The first chain.
     * \param[in] _second The second chain.
     * 
     * \return A new chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int _CTF>
    friend Chain operator+(const Chain &_first, const Chain<_CoefficientType, _CTF> &_second) {
        Chain newChain = _first;
        newChain += _second;
        
        return newChain;
    }
    
    /**
     * \brief Substract two chains together.
     * 
     * Substract each coefficient of the chain together.
     * 
     * \pre The chains have the same coefficent type and the same type flag.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * \warning Will raise an error if the two chains don't have the same type flag.
     * 
     * \param[in] _first The first chain.
     * \param[in] _second The second chain.
     * 
     * \return A new chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int _CTF>
    friend Chain operator-(const Chain &_first, const Chain<_CoefficientType, _CTF> &_second) {
        Chain newChain = _first;
        newChain -= _second;
        
        return newChain;
    }
    
    /**
     * \brief Apply factor on each coefficients.
     * 
     * \param[in] _lambda The factor to apply.
     * \param[in] _chain The second chain.
     * 
     * \return A new chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int _CTF>
    friend Chain operator*(const _CoefficientType& _lambda, const Chain<_CoefficientType, _CTF> &_chain) {
        Chain newChain = _chain;
        newChain *= _lambda;
        
        return newChain;
    }
    
    /**
     * \brief Apply factor on each coefficients.
     * 
     * \param[in] _chain The second chain.
     * \param[in] _lambda The factor to apply.
     * 
     * \return A new chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <int _CTF>
    friend Chain operator*(const Chain<_CoefficientType, _CTF> &_chain, const _CoefficientType& _lambda) {
        Chain newChain = _chain;
        newChain *= _lambda;
        
        return newChain;
    }
    
    /**
     * \brief Perform matrix multiplication between two chains.
     * 
     * Generate a column-based matrix from the matrix multiplication and return it.
     * 
     * \pre The chains have the same coefficent type.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * 
     * \param[in] _column The column chain.
     * \param[in] _row The row chain.
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
    friend SparseMatrix<_CT, COLUMN> operator*(const Chain<_CT, COLUMN> &_column, const Chain<_CT, ROW> &_row);
    
    /**
     * \brief Perform matrix multiplication between two chains.
     * 
     * Generate a row-based matrix from the matrix multiplication and return it.
     * 
     * \pre The chains have the same coefficent type.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * 
     * \param[in] _column The column chain.
     * \param[in] _row The row chain.
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
    friend SparseMatrix<_CT, ROW> operator%(const Chain<_CT, COLUMN> &_column, const Chain<_CT, ROW> &_row);
    
    /**
     * \brief Perform dot product between two chains.
     * 
     * \pre The chains have the same coefficent type.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * 
     * \param[in] _row The row chain.
     * \param[in] _column The column chain.
     * 
     * \return The result of type _CoefficientType.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT>
    friend _CT operator*(const Chain<_CT, ROW> &_row, const Chain<_CT, COLUMN> &_column);
    
    /**
     * \brief Add a chain and assign.
     * 
     * Adds each coefficient of the chain together.
     * 
     * \pre The chains have the same coefficent type and the same type flag.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * \warning Will raise an error if the two chains don't have the same type flag.
     * 
     * \param[in] _other The other chain.
     * 
     * \return The modified chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain& operator+=(const Chain &_other) {
        if (this->upperBound != _other.upperBound) {
            throw std::runtime_error("Chains must be the same size.");
        }
        
        for (pair pair: _other.chainData) {
            this->chainData[pair.first] += pair.second;
            
            if (this->chainData[pair.first] == 0) {
                this->chainData.erase(pair.first);
            }
        }
        
        return *this;
    }
    
    /**
     * \brief Substract a chain and assign.
     * 
     * Substract each coefficient of the chain together.
     * 
     * \pre The chains have the same coefficent type and the same type flag.
     * 
     * \warning Will raise an error if the two chains are not the same coefficient type.
     * \warning Will raise an error if the two chains don't have the same type flag.
     * 
     * \param[in] _other The other chain.
     * 
     * \return The modified chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain& operator-=(const Chain &_other) {
        if (this->upperBound != _other.upperBound) {
            throw std::runtime_error("Chains must be the same size.");
        }
        
        for (pair pair: _other.chainData) {
            this->chainData[pair.first] -= pair.second;
            
            if (this->chainData[pair.first] == 0) {
                this->chainData.erase(pair.first);
            }
        }
        
        return *this;
    }
    
    /**
     * \brief Apply factor on each coefficients and assign.
     * 
     * \param[in] _lambda The factor to apply.
     * 
     * \return The modified chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain& operator*=(const _CoefficientType& _lambda) {
        if (_lambda == 0) {
            this->chainData.clear();
            return *this;
        }
        
        for (pair pair: this->chainData) {
            this->chainData[pair.first] = pair.second * _lambda;
        }
        
        return *this;
    }
    
    /**
     * \brief Get a coefficient from the chain.
     * 
     * \warning The chain will perform boundary check.
     * 
     * \param[in] _index The coefficient index.
     * 
     * \return The coefficient stored in the chain.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    _CoefficientType operator[](const int _index) const {
        if (_index >= upperBound) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(upperBound) + ".");
        }
        
        return chainData.at(_index);
    }
    
    /**
     * \brief Set a coefficient in the chain.
     * 
     * \warning The chain will perform boundary check.
     * 
     * \param[in] _index The coefficient index.
     * 
     * \return The reference to the assigned coefficient.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    _CoefficientType& operator[](const int _index)  {
        if (_index >= upperBound) {
            throw std::runtime_error("Provided index should be less than " + std::to_string(upperBound) + ".");
        }
        
        return chainData[_index];
    }
    
    /**
     * \brief Checks if a coefficient is null. Should be used instead of operator[] for that case.
     * 
     * \param[in] _index The index to check.
     * \return True if the data is null at given index.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 29/05/2024
     */
    const bool isNull(const int _index) const {
        return chainData.find(_index) == chainData.end();
    }
    
    /**
     * \brief Checks if the chain is null.
     * 
     * \return True if the chain is null.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 29/05/2024
     */
    const bool isNull() const {
        return chainData.size() == 0;
    }
    
    /**
     * \brief Get a subchain from the chain.
     * 
     * Removes all indexes provided in the vector from the chain and returns it.
     * 
     * \note Will return a copy of the chain if given vector is empty.
     * 
     * \param[in] _chain The chain to process.
     * \param[in] _indexes The indexes to remove.
     * 
     * \return A new chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT, int _CTF>
    friend Chain<_CT, _CTF> operator/(const Chain<_CT, _CTF> &_chain, const std::vector<int> &_indexes);
    
    /**
     * \brief Get a subchain from the chain.
     * 
     * Removes the index provided from the chain and returns it.
     * 
     * \note Will return a copy of the chain if given vector is empty.
     * 
     * \param[in] _chain The chain to process.
     * \param[in] _index The index to remove.
     * 
     * \return A new chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CT, int _CTF>
    friend Chain<_CT, _CTF> operator/(const Chain<_CT, _CTF> &_chain, const int _indexes);
    
    /**
     * \brief Get a subchain from the chain and assign.
     * 
     * Removes all indexes provided in the vector from the chain and returns it.
     * 
     * \note Will not alter the chain if given vector is empty.
     * 
     * \param[in] _indexes The indexes to remove.
     * 
     * \return The modified chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain& operator/=(const std::vector<int> &_indexes) {
        for (int index : _indexes) {
            this->chainData.erase(index);
        }
        
        return *this;
    }
    
    /**
     * \brief Get a subchain from the chain and assign.
     * 
     * Removes the index provided from the chain and returns it.
     * 
     * \note Will not alter the chain if given vector is empty.
     * 
     * \param[in] _index The index to remove.
     * 
     * \return The modified chain representing the result.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::vector \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain& operator/=(const int _index) {
        this->chainData.erase(_index);
        
        return *this;
    }
    
    /**
     * \brief Set all coefficients to zero.
     * 
     * \author Fedyna K.
     * \version 0.3.0
     * \date 31/04/2024
     */
    void nullify() {
        this->chainData.clear();
    }
    
    /**
     * \brief Iterator to the beginning of the chain.
     * 
     * \warning The chain is stored unordered for speed reason.
     * 
     * \return The iterator to the beginning of the chain.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::unordered_map \endlink
     * \see \link std::unordered_map::iterator \endlink
     * \see \link std::unordered_map::begin \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    iterator begin() noexcept {
        return chainData.begin();
    }
    
    /**
     * \brief Constant iterator to the beginning of the chain.
     * 
     * \warning The chain is stored unordered for speed reason.
     * 
     * \return The constant iterator to the beginning of the chain.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::unordered_map \endlink
     * \see \link std::unordered_map::const_iterator \endlink
     * \see \link std::unordered_map::begin \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    const_iterator begin() const noexcept {
        return chainData.begin();
    }
    
    /**
     * \brief Constant iterator to the beginning of the chain.
     * 
     * \warning The chain is stored unordered for speed reason.
     * 
     * \return The constant iterator to the beginning of the chain.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::unordered_map \endlink
     * \see \link std::unordered_map::const_iterator \endlink
     * \see \link std::unordered_map::cbegin \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    const_iterator cbegin() const noexcept {
        return chainData.cbegin();
    }
    
    /**
     * \brief Iterator to the end of the chain.
     * 
     * \warning The chain is stored unordered for speed reason.
     * 
     * \return The iterator to the end of the chain.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::unordered_map \endlink
     * \see \link std::unordered_map::iterator \endlink
     * \see \link std::unordered_map::end \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    iterator end() noexcept {
        return chainData.end();
    }
    
    /**
     * \brief Constant iterator to the end of the chain.
     * 
     * \warning The chain is stored unordered for speed reason.
     * 
     * \return The constant iterator to the end of the chain.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::unordered_map \endlink
     * \see \link std::unordered_map::const_iterator \endlink
     * \see \link std::unordered_map::end \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    const_iterator end() const noexcept {
        return chainData.end();
    }
    
    /**
     * \brief Constant iterator to the end of the chain.
     * 
     * \warning The chain is stored unordered for speed reason.
     * 
     * \return The constant iterator to the end of the chain.
     * 
     * \see \link OSM::Chain \endlink
     * \see \link std::unordered_map \endlink
     * \see \link std::unordered_map::const_iterator \endlink
     * \see \link std::unordered_map::cend \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    const_iterator cend() const noexcept {
        return chainData.cend();
    }
    
    /**
     * \brief Transpose a Chain.
     * 
     * \return A new chain where the chain type flag is changed.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    Chain<_CoefficientType, COLUMN + ROW - _ChainTypeFlag> transpose() {
        Chain<_CoefficientType, COLUMN + ROW - _ChainTypeFlag> chain;
        
        chain.upperBound = this->upperBound;
        for (pair pair: this->chainData) {
            chain.chainData[pair.first] = pair.second;
        }
        
        return chain;
    }
    
    /**
     * \brief Checks if chain is a column.
     * 
     * \return true if chain is represented as a column, false otherwise.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    bool isColumn() const {
        return _ChainTypeFlag == COLUMN;
    }
    
    /**
     * \brief Checks if chain is a row.
     * 
     * \return true if chain is represented as a row, false otherwise.
     * 
     * \see \link OSM::Chain \endlink
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 17/04/2024
     */
    bool isRow() const {
        return _ChainTypeFlag == ROW;
    }
};


/**
 * \brief Perform matrix multiplication between two chains.
 * 
 * Generate a column-based matrix from the matrix multiplication and return it.
 * 
 * \pre The chains have the same coefficent type.
 * 
 * \warning Will raise an error if the two chains are not the same coefficient type.
 * 
 * \param[in] _column The column chain.
 * \param[in] _row The row chain.
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
SparseMatrix<_CT, COLUMN> operator*(const Chain<_CT, COLUMN> &_column, const Chain<_CT, ROW> &_row) {
    SparseMatrix<_CT, COLUMN> matrix(_column.upperBound, _row.upperBound);
    
    for (std::pair<int, _CT> pair : _row.chainData) {
        OSM::setColumn(matrix,pair.first,_column * pair.second) ;
    }
    
    return matrix;
}

/**
 * \brief Perform matrix multiplication between two chains.
 * 
 * Generate a row-based matrix from the matrix multiplication and return it.
 * 
 * \pre The chains have the same coefficent type.
 * 
 * \warning Will raise an error if the two chains are not the same coefficient type.
 * 
 * \param[in] _column The column chain.
 * \param[in] _row The row chain.
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
SparseMatrix<_CT, ROW> operator%(const Chain<_CT, COLUMN> &_column, const Chain<_CT, ROW> &_row) {
    SparseMatrix<_CT, ROW> matrix(_column.upperBound, _row.upperBound);
    
    for (std::pair<int, _CT> pair : _column.chainData) {
        OSM::setRow(matrix,pair.first,_row * pair.second);
    }
    
    return matrix;
}

/**
 * \brief Perform dot product between two chains.
 * 
 * \pre The chains have the same coefficent type.
 * 
 * \warning Will raise an error if the two chains are not the same coefficient type.
 * 
 * \param[in] _row The row chain.
 * \param[in] _column The column chain.
 * 
 * \return The result of type _CoefficientType.
 * 
 * \see \link OSM::Chain \endlink
 * 
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CoefficientType>
_CoefficientType operator*(const Chain<_CoefficientType, ROW> &_row, const Chain<_CoefficientType, COLUMN> &_column) {
    // Get indexes (avoid adding double indexes).
    std::unordered_map<int, int> indexes;
    for (std::pair<int, _CoefficientType> pair: _row.chainData) {
        indexes[pair.first] = 1;
    }
    for (std::pair<int, _CoefficientType> pair: _column.chainData) {
        indexes[pair.first] += 1;
    }
    
    // Perform dot product
    _CoefficientType result = _CoefficientType();
    for (std::pair<int, int> index: indexes) {
        if (index.second == 2) {
            result += _row.chainData.at(index.first) * _column.chainData.at(index.first);
        }
    }
    
    return result;
}

/**
 * \brief Get a subchain from the chain and assign.
 * 
 * Removes all indexes provided in the vector from the chain and returns it.
 * 
 * \note Will not alter the chain if given vector is empty.
 * 
 * \param[in] _indexes The indexes to remove.
 * 
 * \return The modified chain representing the result.
 * 
 * \see \link OSM::Chain \endlink
 * \see \link std::vector \endlink
 * 
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT, int _CTF>
Chain<_CT, _CTF> operator/(const Chain<_CT, _CTF> &_chain, const std::vector<int> &_indexes) {
    Chain newChain = _chain;
    newChain /= _indexes;
    return newChain;
}

/**
 * \brief Get a subchain from the chain and assign.
 * 
 * Removes the index provided from the chain and returns it.
 * 
 * \note Will not alter the chain if given vector is empty.
 * 
 * \param[in] _index The index to remove.
 * 
 * \return The modified chain representing the result.
 * 
 * \see \link OSM::Chain \endlink
 * \see \link std::vector \endlink
 * 
 * \author Fedyna K.
 * \version 0.1.0
 * \date 08/04/2024
 */
template <typename _CT, int _CTF>
Chain<_CT, _CTF> operator/(const Chain<_CT, _CTF> &_chain, const int _index) {
    Chain newChain = _chain;
    newChain /= _index;
    return newChain;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif
