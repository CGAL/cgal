// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//                 Andreas Fabri
// =============================================================================

#ifndef CGAL_FRECHET_DISTANCE_NEIGHBOR_SEARCH_H
#define CGAL_FRECHET_DISTANCE_NEIGHBOR_SEARCH_H

#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance.h>
#include <CGAL/basic.h>
#include <CGAL/Frechet_distance/internal/Neighbor_search.h>

#include <iterator>
#include <vector>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

namespace CGAL {
namespace Frechet_distance {

/*!
* \ingroup PkgFrechetDistanceRef
 * A data structure to store polylines with a function that enables to find those polylines which are closer than a distance bound to a query polyline.
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer` with `Traits::Point_d` as value type.
*/
template <class PointRange, class Traits>
class Neighbor_search
{
    using PT = Traits; //  Polyline_traits_2<Traits, double>;
    using FT = typename PT::FT;
    using Point = typename PT::Point_d;
    using Polyline = PointRange;
    using Polylines = std::vector<Polyline>;

public:

/*!
 * constructs a neighbor search data structure for `polylines`
 * \tparam ConcurrencyTag enables sequential versus parallel construction.
 * Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * \tparam PolylineRange must be a model of `RandomAccessRange` with value type `PointRange`
*/
    template <typename ConcurrencyTag = Sequential_tag, typename PolylineRange>
    Neighbor_search(const PolylineRange& polylines, ConcurrencyTag tag = ConcurrencyTag())
    : curves(polylines.begin(), polylines.end())
    {
      CGAL_USE(tag);
      kd_tree.insert(curves);
      kd_tree.template build<ConcurrencyTag>();
    }
#ifdef DOXYGEN_RUNNING
    using Point = Traits::Point_d;
#endif

     /*!  returns the indices of the inserted polylines that are closer than `distance` to the polyline `query`.
        * \param query the query polyline
        * \param distance the distance bound
        * \tparam ConcurrencyTag enables sequential versus parallel computation.
        * Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
        * \return a vector of indices of the inserted polylines that are closer than `distance` to the polyline `query`.
     */
    template <typename ConcurrencyTag = Sequential_tag>
    std::vector<std::size_t> get_close_curves(const PointRange& query, double distance)
    {
#ifndef CGAL_LINKED_WITH_TBB
        static_assert(!std::is_same_v<ConcurrencyTag, Parallel_tag>);
#endif
        return get_close_curves(query, distance, ConcurrencyTag());
    }

private:


    std::vector<std::size_t> get_close_curves(const PointRange& query, double distance, Sequential_tag);
#ifdef CGAL_LINKED_WITH_TBB
    std::vector<std::size_t> get_close_curves(const PointRange& query, double distance, Parallel_tag);
#endif

    Polylines curves;
    Frechet_distance::internal::FrechetKdTree<Traits> kd_tree;
};


template <class PointRange, class Traits>
#ifdef DOXYGEN_RUNNING
std::vector<std::size_t>
#else
auto
#endif
Neighbor_search<PointRange, Traits>::get_close_curves(
    const PointRange& curve, double distance, Sequential_tag) -> std::vector<std::size_t>
{
    auto result = kd_tree.search(curve, distance);

    auto predicate = [&](std::size_t id) {
      CGAL_assertion(id < curves.size());
        return  is_Frechet_distance_larger(
            curve, curves[id], distance, parameters::geom_traits(Traits()));
    };
 #ifdef CGAL_FRECHET_VERBOSE
    /std::cout << result.size() << " curves returned from tree" << std::endl;
#endif
    auto new_end = std::remove_if(result.begin(), result.end(), predicate);
    result.erase(new_end, result.end());

    return result;
}

#ifdef CGAL_LINKED_WITH_TBB
template <class PointRange, class Traits>
#ifdef DOXYGEN_RUNNING
std::vector<std::size_t>
#else
auto
#endif
Neighbor_search<PointRange, Traits>::get_close_curves(
    const PointRange& curve, double distance, Parallel_tag) -> std::vector<std::size_t>
{
    std::vector<std::size_t> result;
    result = kd_tree.search(curve, distance);
    std::vector<int> to_remove(result.size(), 0);

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, result.size()), [&](const tbb::blocked_range<std::size_t>& r) {
        for (std::size_t i = r.begin(); i != r.end(); ++i) {
            to_remove[i] = is_Frechet_distance_larger(curve, curves[result[i]], distance, parameters::geom_traits(Traits()));
        }
    });

    int compact = 0;
    for(std::size_t i = 0; i < to_remove.size(); ++i) {
        if(! to_remove[i]) {
            result[compact] = result[i];
            ++compact;
        }
    }
    result.resize(compact);
    return result;
}
#endif

}
}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_NEIGHBOR_SEARCH_H
