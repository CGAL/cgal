// Copyright (c) 2026 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Sebastien Loriot <sebastien.loriot@cgal.org>
//

#ifndef CGAL_TETRAHEDRON_SOUP_TO_LCC_H
#define CGAL_TETRAHEDRON_SOUP_TO_LCC_H

#include <boost/container_hash/hash.hpp>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>

namespace CGAL
{
/**
 * @ingroup PkgLinearCellComplexConstructions
 *
 * \brief imports a 3D tetrahedron soup into a linear cell complex.
 *
 * Creates a 3D linear cell complex in `lcc` from a set of points and a range
 * of tetrahedra defined by indices into `points`. The topological 3-combinatorial
 * map is constructed by sewing adjacent 3-volumes (tetrahedra) along matching 2-faces.
 *
 * @tparam LCC a model of the `LinearCellComplex` concept.
 * @tparam PointRange a model of `RandomAccessContainer` with `Point_3` as value type, being compatible with the point type of LCC.
 * @tparam TetraRange a model of `ConstRange` where each element is a container of 4 point indices, accessible using `operator[](int)`.
 *
 * @param points range of 3D points.
 * @param tetras range of tetrahedra, where each tetrahedron is represented by 4 indices corresponding to entries in `points`.
 * @param lcc the target linear cell complex.
 *
 * @return a descriptor to a dart of the created complex, or `LCC::null_descriptor` if the input is empty.
 *
 * @pre `LCC::dimension >= 3` and `LCC::ambient_dimension == 3`.
 * @pre Indices in `tetras` must be valid 0-based indices into `points`.
 *
 * @sa `CGAL::triangulation_3_to_lcc()`
 */
  template <class LCC, class PointRange, class TetraRange>
  typename LCC::Dart_descriptor
  tetrahedron_soup_to_lcc(const PointRange& points,
                          const TetraRange& tetras,
                          LCC& lcc)
  {
    static_assert( LCC::dimension>=3 && LCC::ambient_dimension==3 );

    if (points.empty() || tetras.empty())
      return LCC::null_descriptor;

    using Dart_descriptor = typename LCC::Dart_descriptor;

    std::vector<typename LCC::Vertex_attribute_descriptor> vertices;
    vertices.reserve(points.size());
    for(const auto& pt : points)
      vertices.push_back(lcc.create_vertex_attribute(pt));

    struct Array_hasher
    {
      std::size_t operator()(const std::array<std::size_t, 3>& a) const
      {
        std::size_t seed = 0;
        boost::hash_combine(seed, a[0]);
        boost::hash_combine(seed, a[1]);
        boost::hash_combine(seed, a[2]);

        return seed;
      }
    };

    std::unordered_map< std::array<std::size_t, 3>, Dart_descriptor, Array_hasher> face_map;

    Dart_descriptor dart = LCC::null_descriptor;
    for (const auto& t : tetras)
    {
      Dart_descriptor res = lcc.make_tetrahedron(vertices[t[0]],
                                                 vertices[t[1]],
                                                 vertices[t[2]],
                                                 vertices[t[3]]);

      if (dart==LCC::null_descriptor) dart = res;

      for (int i=0; i<4; ++i)
      {
        Dart_descriptor curr= LCC::null_descriptor;
        switch (i)
        {
          case 0: curr = lcc.template opposite<2>(lcc.next(res)); break;
          case 1: curr = lcc.template opposite<2>(lcc.previous(res)); break;
          case 2: curr = lcc.template opposite<2>(res); break;
          default: curr = res; break;
        }

        std::array<std::size_t, 3> f{t[(i+1)%4],t[(i+2)%4],t[(i+3)%4]};
        std::sort(f.begin(), f.end());
        auto insert_res = face_map.emplace(f, curr);
        if (!insert_res.second)
        {
          Dart_descriptor curr2=insert_res.first->second;
          while (lcc.vertex_attribute(curr2) !=
                 lcc.vertex_attribute(lcc.other_extremity(curr)) )
          {
            curr2 = lcc.next(curr2);
          }
          lcc.template topo_sew<3>(curr, lcc.other_orientation(curr2));
        }
      }
    }

    CGAL_assertion(dart!=LCC::null_descriptor);
    return dart;
  }
} // end of CGAL namespace


#endif // CGAL_TETRAHEDRON_SOUP_TO_LCC_H
