/**
 * \file __base.hpp
 * \brief Defines all constants used by objects.
 * \author Fedyna K.
 * \version 0.1.0
 * \date 24/04/2024
 * 
 * Should not be included.
 */

#ifndef __OPTIMISED_SPARSED_MATRIX_BASE__
#define __OPTIMISED_SPARSED_MATRIX_BASE__

#include <ostream>


namespace OSM {
    /** \brief Chain type flag for column chain. */
    const int COLUMN = 0b01;
    /** \brief Chain type flag for row chain. */
    const int ROW    = 0b10;

    /** \brief The default type for signed integers. */
    typedef int ZCoefficient;

    /**
     * \class SparseMatrix
     * \brief Vector<Map> implementation of sparse matrices.
     * 
     * The SparseMatrix class contains all algebraic functions that are related to matrix.
     * 
     * \tparam _CoefficientType The chain's coefficient types (default is OSM::ZCoefficient)
     * \tparam _ChainTypeFlag The type of vector the chain is representing (default is OSM::COLUMN)
     * 
     * \author Fedyna K.
     * \version 0.1.0
     * \date 08/04/2024
     */
    template <typename _CoefficientType = OSM::ZCoefficient, int _ChainTypeFlag = OSM::COLUMN>
    class SparseMatrix;

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
    template <typename _CoefficientType = OSM::ZCoefficient, int _ChainTypeFlag = OSM::COLUMN>
    class Chain;

    /**
     * \class Bitboard
     * \brief Bitboards are long bitset with specific operations.
     * 
     * Bitboards are used with SparseMatrix to quickly find and set non empty columns.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    class Bitboard;
}

#endif