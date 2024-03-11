/**
 * @file   algo/3d/LineInFacet.h
 * @author Gernot Walzl
 * @date   2012-04-26
 */

#ifndef ALGO_3D_LINEINFACET_H
#define ALGO_3D_LINEINFACET_H

#include "data/3d/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;

/**
 * Checks if the intersection point of the given line is inside the facet.
 * The algorithm is based on the ray casting algorithm
 * that solves the point-in-polygon problem.
 */
bool IsLineInFacet(FacetSPtr facet, Line3SPtr line);

} }

#endif /* ALGO_3D_LINEINFACET_H */
