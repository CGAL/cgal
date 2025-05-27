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

namespace CGAL {
namespace OSM {

/** \brief Chain type flag for column chain. */
const int COLUMN = 0b01;
/** \brief Chain type flag for row chain. */
const int ROW    = 0b10;

/** \brief The default type for signed integers. */
typedef int ZCoefficient;

// Class Sparse_matrix
template <typename _CoefficientType = OSM::ZCoefficient, int _ChainTypeFlag = OSM::COLUMN>
class Sparse_matrix;

// Class Sparse_chain

template <typename _CoefficientType = OSM::ZCoefficient, int _ChainTypeFlag = OSM::COLUMN>
class Sparse_chain;

// Class Bitboard

class Bitboard;

} /* end namespace OSM */
} /* end namespace CGAL */

#endif
