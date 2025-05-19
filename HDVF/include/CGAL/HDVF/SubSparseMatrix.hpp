/**
 * \file SubSparseMatrix.hpp
 * \brief Namespace file for describing library.
 * \author Bac A.
 * \version 0.1.0
 * \date 02/10/2024
 *
 * Define everything for the SubSparseMatrix class
 */

#ifndef __SUB_SPARSE_MATRIX__
#define __SUB_SPARSE_MATRIX__


#include "CGAL/OSM/Sparse_matrix.h"
#include "CGAL/OSM/Bitboard.hpp"

namespace CGAL {
namespace OSM {

/**
 * \class SubSparseMatrix
 * \brief Vector<Map> implementation of sparse matrices with a bitboard submatrix mask.
 *
 * The SubSparseMatrix class contains functions related to OSM matrix restriction.
 *
 * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
 * \tparam _ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 02/10/2024
 */
template <typename _CoefficientType, int _ChainTypeFlag>
class SubSparseMatrix : public SparseMatrix<_CoefficientType, _ChainTypeFlag> {
    
protected:
    /** \brief A bitboard describing subchains restriction. */
    Bitboard _subChains;
    
    /** \brief A bitboard containing state of each chain (restricted to subchains). */
    Bitboard _subChainsStates;
    
public:
    /** \brief Constructor with sizes and Bitboard describing subchains. */
    SubSparseMatrix(const int _rowCount, const int _columnCount, const Bitboard& subChain) : SparseMatrix<_CoefficientType, _ChainTypeFlag>(_rowCount, _columnCount), _subChains(subChain), _subChainsStates(this->chainsStates & subChain)
    {
    }
    
    /** \brief Constructor with sizes (Bitboard describing subchains set to full Bitboard). */
    SubSparseMatrix(const int _rowCount=0, const int _columnCount=0) : SparseMatrix<_CoefficientType, _ChainTypeFlag>(_rowCount, _columnCount)
    {
        if (_ChainTypeFlag == OSM::COLUMN)
            _subChains = OSM::Bitboard(_columnCount,false) ;
        else
            _subChains = OSM::Bitboard(_rowCount,false) ;
        _subChainsStates = this->chainsStates & _subChains ;
    }
    
    /** \brief Copy constructor. */
    SubSparseMatrix(const SubSparseMatrix& otherToCopy) : SparseMatrix<_CoefficientType, _ChainTypeFlag>(otherToCopy), _subChains(otherToCopy._subChains), _subChainsStates(otherToCopy._subChainsStates) {}
    
    /** \brief Copy constructor from SparseMatrix. */
    SubSparseMatrix(const SparseMatrix<_CoefficientType,_ChainTypeFlag>& otherToCopy) : SparseMatrix<_CoefficientType, _ChainTypeFlag>(otherToCopy)
    {
        if (_ChainTypeFlag == OSM::COLUMN)
            _subChains = OSM::Bitboard(otherToCopy.dimensions().second,false) ;
        else
            _subChains = OSM::Bitboard(otherToCopy.dimensions().first,false) ;
        _subChainsStates = this->chainsStates & _subChains ;
    }
    
    
    /** \brief Iterator to the beginning of the subchains indices. */
    inline Bitboard::iterator begin() const noexcept
    {
        return _subChainsStates.begin() ;
    }
    
    /** \brief Iterator to the end of the subchains indices. */
    inline Bitboard::iterator end() const noexcept
    {
        return _subChainsStates.end() ; 
    }
    
    /** \brief Change subchains. */
    inline void set_sub (const Bitboard& new_subChains)
    {
        _subChains = new_subChains ;
        _subChainsStates = this->chainsStates & _subChains ;
    }
    
    /** \brief Set a bit on in subchains and update chainsStates accordingly. */
    inline void set_bitOn (int index)
    {
        _subChains.setOn(index) ;
        if (this->chainsStates.isOn(index))
            _subChainsStates.setOn(index) ;
    }
    
    /** \brief Set a bit off in subchains and update chainsStates accordingly. */
    inline void set_bitOff (int index)
    {
        _subChains.setOff(index) ;
        _subChainsStates.setOff(index) ;
    }
    
    /** \brief Change subchains to their complement. */
    inline void complement() { _subChains.bit_not() ; }
    
    /** \brief Operator=. */
    inline SubSparseMatrix& operator=(const SubSparseMatrix &other)
    {
        (dynamic_cast<SparseMatrix<_CoefficientType,_ChainTypeFlag>&>(*this)).operator=(other) ;
        _subChains = other._subChains ;
        _subChainsStates = other._subChainsStates ;
        return *this ;
    }
    
    /** \brief Operator<<. */
    friend std::ostream& operator<<(std::ostream &_stream, const SubSparseMatrix &_matrix) {
        _stream << static_cast<const SparseMatrix<_CoefficientType, _ChainTypeFlag>&>(_matrix) ;
        _stream << _matrix._subChains << std::endl ;
        return _stream ;
    }
};

} /* end namespace OSM */
} /* end namespace CGAL */

#endif
