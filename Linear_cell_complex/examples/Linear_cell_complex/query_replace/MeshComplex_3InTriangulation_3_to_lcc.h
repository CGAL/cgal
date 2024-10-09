// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef MESH_COMPLEX_3_IN_TRIANGULATION_3_TO_LCC_H
#define MESH_COMPLEX_3_IN_TRIANGULATION_3_TO_LCC_H

#include <CGAL/assertions.h>
#include <unordered_map>
#include <CGAL/Weighted_point_3.h>

namespace internal 
{
template<typename Point>
struct Get_point
{
    static const Point& run(const Point& p)
    { return p; }
};

template<typename Kernel>
struct Get_point<CGAL::Weighted_point_3<Kernel> >
{
    static const typename Kernel::Point_3& run(const CGAL::Weighted_point_3<Kernel>& p)
    { return p.point(); }
};
}

/** Convert a given Triangulation_3 into a 3D linear cell complex.
 * @param alcc the used linear cell complex.
 * @param atr the Triangulation_3.
 * @param avol_to_dart a pointer to a std::map associating to each
 *        tetrahedron of atr a corresponding dart in alcc. Not used if nullptr.
 * @return A dart incident to the infinite vertex.
 */
template < class LCC, class Triangulation >
typename LCC::Dart_handle import_from_complex_3_in_triangulation_3
(LCC& alcc, const Triangulation &atr,
std::unordered_map<typename Triangulation::Cell_handle,
 typename LCC::Dart_handle >* avol_to_dart=nullptr)
{
  CGAL_static_assertion( LCC::dimension>=3 && LCC::ambient_dimension==3 );

  // Case of empty triangulations.
  if (atr.triangulation().number_of_vertices()==0) return LCC::null_handle;

  // Check the dimension.
  if (atr.triangulation().dimension()!=3) return LCC::null_handle;
  CGAL_assertion(atr.is_valid());

  typedef typename Triangulation::Vertex_handle TVertex_handle;
  typedef typename Triangulation::Cell_iterator TCell_iterator;

  // Create vertices in the map and associate in a map
  // TVertex_handle and vertices in the map.
  std::unordered_map< TVertex_handle, typename LCC::Vertex_attribute_handle > TV;
  for (auto itv=atr.triangulation().vertices_begin();
       itv!=atr.triangulation().vertices_end(); ++itv)
  {
    TV[itv]=alcc.create_vertex_attribute
            (internal::Get_point<typename Triangulation::Point>::run(itv->point()));
  }

  // Create the tetrahedron and create a map to link Cell_iterator
  // and tetrahedron.
  TCell_iterator it;
  std::unordered_map<typename Triangulation::Cell_handle, typename LCC::Dart_handle> TC;
  std::unordered_map<typename Triangulation::Cell_handle, typename LCC::Dart_handle>*
      mytc = (avol_to_dart==nullptr?&TC:avol_to_dart);

  typename LCC::Dart_handle res=LCC::null_handle, dart=LCC::null_handle;
  typename LCC::Dart_handle cur=LCC::null_handle, neighbor=LCC::null_handle;

  for (it=atr.cells_in_complex_begin(); it!=atr.cells_in_complex_end(); ++it)
  {
    if (it->vertex(0)!=atr.triangulation().infinite_vertex() &&
        it->vertex(1)!=atr.triangulation().infinite_vertex() &&
        it->vertex(2)!=atr.triangulation().infinite_vertex() &&
        it->vertex(3)!=atr.triangulation().infinite_vertex())
    {
      assert(TV.count(it->vertex(0))==1);
      assert(TV.count(it->vertex(1))==1);
      assert(TV.count(it->vertex(2))==1);
      assert(TV.count(it->vertex(3))==1);

      res=alcc.make_tetrahedron(TV[it->vertex(0)],
          TV[it->vertex(1)],
          TV[it->vertex(2)],
          TV[it->vertex(3)]);

      if ( dart==LCC::null_handle )
      { dart=res; }

      for (unsigned int i=0; i<4; ++i)
      {
        assert(&*(it->neighbor(i)->neighbor(atr.triangulation().mirror_index(it, i)))==&*it);
        assert(&*(it->neighbor(i))!=&*it);

        switch (i)
        {
          case 0: cur=alcc.template opposite<2>(alcc.next(res)); break;
          case 1: cur=alcc.template opposite<2>(alcc.previous(res)); break;
          case 2: cur=alcc.template opposite<2>(res); break;
          case 3: cur=res; break;
        }

        auto maptcell_it=mytc->find(it->neighbor(i));
        if (maptcell_it!=mytc->end())
        {
          switch (atr.triangulation().mirror_index(it,i) )
          {
            case 0: neighbor=alcc.template opposite<2>(alcc.next(maptcell_it->second));
              break;
            case 1: neighbor=alcc.template opposite<2>(alcc.previous(maptcell_it->second));
              break;
            case 2: neighbor=alcc.template opposite<2>(maptcell_it->second); break;
            case 3: neighbor=maptcell_it->second; break;
          }

          while (alcc.vertex_attribute(neighbor)!=
                 alcc.vertex_attribute(alcc.other_extremity(cur)))
          { neighbor = alcc.next(neighbor); }

          assert(alcc.vertex_attribute(alcc.other_extremity(cur))==
                 alcc.vertex_attribute(alcc.other_orientation(neighbor)));
          assert(alcc.vertex_attribute(alcc.other_extremity(alcc.next(cur)))==
                 alcc.vertex_attribute(alcc.other_orientation(alcc.previous(neighbor))));
          assert(alcc.vertex_attribute(alcc.other_extremity(alcc.previous(cur)))==
                 alcc.vertex_attribute(alcc.other_orientation(alcc.next(neighbor))));
          assert(alcc.template is_free<3>(cur));
          assert(alcc.template is_free<3>(alcc.other_orientation(neighbor)));

          alcc.template topo_sew<3>(cur, alcc.other_orientation(neighbor));
        }
      }
      (*mytc)[it]=res;
    }
  }
  CGAL_assertion(dart!=LCC::null_handle);
  alcc.correct_invalid_attributes(); // Necessary to correct problem of two
  // tetra adjacent by a vertex (non manifold case).
  return dart;
}

#endif // MESH_COMPLEX_3_IN_TRIANGULATION_3_TO_LCC_H
////////////////////////////////////////////////////////////////////////////////
