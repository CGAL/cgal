// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef CMAP_QUERY_REPLACE_GEOMETRY_H
#define CMAP_QUERY_REPLACE_GEOMETRY_H

#include <tuple>
#include <unordered_map>

#include <CGAL/Linear_cell_complex/lcc_barycentric_coord.h>
#include <CGAL/Linear_cell_complex/query_replace/lcc_pattern.h>

namespace CGAL::internal
{
template<class LCC>
using Dart_mapping=std::unordered_map<typename LCC::Dart_descriptor,
                                      typename LCC::Dart_descriptor>;

///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Point
compute_point_2D(LCC& lcc,
                 Dart_mapping<LCC>& links_from_pattern_to_face,
                 Dart_mapping<LCC>& pattern_to_global,
                 Barycentric_coord<LCC, 1>& m_barycentric_coords)
{
  typename LCC::Point p;
  typename LCC::Vector res=CGAL::NULL_VECTOR;
  typename LCC::Dart_descriptor cur, dh1, dh2;
  std::size_t nb=0;
  typename LCC::Dart_descriptor firstdh=
      links_from_pattern_to_face[pattern_to_global
                                     [std::get<0>(m_barycentric_coords.m_coords.front())]];
  typename LCC::Point bary2=lcc.template barycenter<2>(firstdh);
  for(std::tuple<typename LCC::Dart_descriptor, double, double, double>& e:
       m_barycentric_coords.m_coords)
  {
    CGAL_assertion(pattern_to_global.find(std::get<0>(e))!=pattern_to_global.end());
    CGAL_assertion(links_from_pattern_to_face.find(pattern_to_global[std::get<0>(e)])
           !=links_from_pattern_to_face.end());
    cur=links_from_pattern_to_face[pattern_to_global[std::get<0>(e)]];
    dh1=lcc.other_extremity(cur);
    if(compute_point_from_alpha_beta_gamma(lcc.point(cur), lcc.point(dh1),
                                            bary2, std::get<1>(e),
                                            std::get<2>(e), std::get<3>(e), p))
    {
      res+=typename LCC::Vector(p.x(), p.y(), p.z());
      ++nb;
    }
  }
  CGAL_assertion(nb>0);
  return typename LCC::Point(res.x()/nb, res.y()/nb, res.z()/nb);
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Point
compute_point_3D(LCC& lcc,
                 Dart_mapping<LCC>& links_from_pattern_to_volume,
                 Dart_mapping<LCC>& pattern_to_global,
                 Barycentric_coord<LCC, 3>& m_barycentric_coords)
{
  typename LCC::Point p;
  typename LCC::Vector res=CGAL::NULL_VECTOR;
  typename LCC::Dart_descriptor cur, dh1, dh2, dh3;
  std::size_t nb=0;
  for(std::tuple<typename LCC::Dart_descriptor, double, double, double, double>&
           e: m_barycentric_coords.m_coords)
  {
    CGAL_assertion(pattern_to_global.find(std::get<0>(e))!=pattern_to_global.end());
    CGAL_assertion(links_from_pattern_to_volume.find(pattern_to_global[std::get<0>(e)])
        !=links_from_pattern_to_volume.end());
    cur=links_from_pattern_to_volume[pattern_to_global[std::get<0>(e)]];
    dh1=lcc.template beta<0>(cur);
    dh2=lcc.other_extremity(cur);
    dh3=lcc.template beta<2,1,2>(cur);
    if(compute_point_from_alpha_beta_gamma_delta(lcc.point(cur), lcc.point(dh1),
                                                 lcc.point(dh2), lcc.point(dh3),
                                                 std::get<1>(e), std::get<2>(e),
                                                 std::get<3>(e), std::get<4>(e),
                                                 p))
    {
      res+=typename LCC::Vector(p.x(), p.y(), p.z());
      ++nb;
    }
  }
  CGAL_assertion(nb>0);
  return typename LCC::Point(res.x()/nb, res.y()/nb, res.z()/nb);
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Point
compute_point_3D_v2(LCC& lcc,
                    Dart_mapping<LCC>& links_from_pattern_to_volume,
                    Dart_mapping<LCC>& pattern_to_global,
                    Barycentric_coord<LCC, 3>& m_barycentric_coords)
{
  CGAL_assertion(!m_barycentric_coords.m_coords.empty());
  typename LCC::Point p;
  typename LCC::Vector res=CGAL::NULL_VECTOR;
  typename LCC::Dart_descriptor cur, dh1, dh2, dh3;
  std::size_t nb=0;
  typename LCC::Dart_descriptor firstdh=
      links_from_pattern_to_volume[pattern_to_global
      [std::get<0>(m_barycentric_coords.m_coords.front())]];
  typename LCC::Point bary3=lcc.template barycenter<3>(firstdh);
  for(std::tuple<typename LCC::Dart_descriptor, double, double, double, double>& e:
      m_barycentric_coords.m_coords)
  {
    CGAL_assertion(pattern_to_global.find(std::get<0>(e))!=pattern_to_global.end());
    CGAL_assertion(links_from_pattern_to_volume.find(pattern_to_global[std::get<0>(e)])
        !=links_from_pattern_to_volume.end());
    cur=links_from_pattern_to_volume[pattern_to_global[std::get<0>(e)]];
    dh1=lcc.other_extremity(cur);

    // TODO avoid to recompute barycenters several times (?)
    if(compute_point_from_alpha_beta_gamma_delta(lcc.point(cur), lcc.point(dh1),
                                                 lcc.template barycenter<2>(cur),
                                                 bary3,
                                                 std::get<1>(e), std::get<2>(e),
                                                 std::get<3>(e), std::get<4>(e),
                                                 p))
    {
      res+=typename LCC::Vector(p.x(), p.y(), p.z());
      ++nb;
    }
  }
  CGAL_assertion(nb>0);
  return typename LCC::Point(res.x()/nb, res.y()/nb, res.z()/nb);
}
///////////////////////////////////////////////////////////////////////////////
/// Transform the geometry of the fpattern according to the geometry of the
/// target.
template<typename LCC>
void transform_geometry_of_fpattern(LCC& lcc,
                                    Dart_mapping<LCC>& links_from_pattern_to_face,
                                    Dart_mapping<LCC>& pattern_to_global,
                                    Pattern<LCC, 1>& pattern)
{
  for(Barycentric_coord<LCC, 1>& inner: pattern.barycentric_coords())
  {
    CGAL_assertion(pattern_to_global.find(inner.m_dart)!=pattern_to_global.end());
    typename LCC::Dart_descriptor res=pattern_to_global[inner.m_dart];
    // TODO avoid to recompute barycenters several times (?)
    lcc.point(res)=compute_point_2D(lcc,
                                    links_from_pattern_to_face,
                                    pattern_to_global,
                                    inner);
  }
}
///////////////////////////////////////////////////////////////////////////////
/// Transform the geometry of the spattern according to the geometry of the
/// target. For now same method than transform_geometry_of_fpattern
template<typename LCC>
void transform_geometry_of_spattern(LCC& lcc,
                                    Dart_mapping<LCC>& links_from_pattern_to_face,
                                    Dart_mapping<LCC>& pattern_to_global,
                                    Pattern<LCC, 2>& pattern)
{
  for(Barycentric_coord<LCC, 1>& inner: pattern.barycentric_coords())
  {
    CGAL_assertion(pattern_to_global.find(inner.m_dart)!=pattern_to_global.end());
    typename LCC::Dart_descriptor res=pattern_to_global[inner.m_dart];
    // TODO avoid to recompute barycenters several times (?)
    lcc.point(res)=compute_point_2D(lcc,
                                    links_from_pattern_to_face,
                                    pattern_to_global,
                                    inner);
  }
}
////////////////////////////////////////////////////////////////////////////////
/// Transform the geometry of the vpattern according to the geometry of the
/// target. All darts to keep (between the query and the replace) are marked
/// with amark.
/// For now simple solution that does not work for any pattern. TODO better? TODO is this comment still true?
template<typename LCC>
void transform_geometry_of_vpattern(LCC& lcc, typename LCC::size_type /*amark*/,
                                    Dart_mapping<LCC>& links_from_pattern_to_volume,
                                    Dart_mapping<LCC>& pattern_to_global,
                                    Pattern<LCC, 3>& pattern)
{
/*  /// 1) Compute the oriented bounding box of surviving points
  typename LCC::size_type surviving_vertex=lcc.get_new_mark();
  std::unordered_set<typename LCC::Point> points;
  for(const auto& it: links_from_pattern_to_volume) // This is only "surviving" darts
  {
    CGAL_assertion(lcc.is_marked(it.first, amark));
    if(!lcc.is_marked(it.first, surviving_vertex))
    {
      points.insert(lcc.point(it.second)); // Survinv point in the match
      lcc.template mark_cell<0>(it.first, surviving_vertex);
    }
  }
  std::array<typename LCC::Point, 8> obb_points;
  CGAL::oriented_bounding_box(points, obb_points,
                              CGAL::parameters::use_convex_hull(true));

  / * std::cout<<"OBB: ";
  for(std::size_t i=0; i<8; ++i) {std::cout<<obb_points[i]<<"  ";}
  std::cout<<std::endl; * /

  /// 2) Compute new vertices coordinates
  std::vector<double> coords;
  for(const auto& it: pattern_to_global) // This is all darts of the copy of the pattern
  {
    if(!lcc.is_marked(it.second, surviving_vertex))
    {
      coords.clear();
      pattern.mv->operator()((const typename LCC::Point)(lcc.point(it.second)),
                             std::back_inserter(coords)); // ici on utilise les coordonnées barycentrique
      double x=0., y=0., z=0.;
      for(std::size_t i=0; i<8; ++i)
      {
        x+=obb_points[i].x()*coords[i];
        y+=obb_points[i].y()*coords[i];
        z+=obb_points[i].z()*coords[i];
      }
      // std::cout<<"Update: "<<lcc.point(it.second)<<"  ->   ";
      lcc.point(it.second)=typename LCC::Point(x, y, z);
      // std::cout<<lcc.point(it.second)<<std::endl;
      lcc.template mark_cell<0>(it.second, surviving_vertex);
    }
  }

  /// 3) Unmark vertices
  for(const auto& it: pattern_to_global)
  {
    if(lcc.is_marked(it.second, surviving_vertex))
    { lcc.template unmark_cell<0>(it.second, surviving_vertex); }
  }
  CGAL_assertion(lcc.is_whole_map_unmarked(surviving_vertex));
  lcc.free_mark(surviving_vertex);
*/
  // return;
  for(Barycentric_coord<LCC, 3>& inner: pattern.barycentric_coords())
  {
    CGAL_assertion(pattern_to_global.find(inner.m_dart)!=pattern_to_global.end());
    typename LCC::Dart_descriptor res=pattern_to_global[inner.m_dart];
    // std::cout<<"[transform_geometry_of_vpattern] "<<lcc.point()<<" -> before "
    //          <<lcc.point()<<" and  after ";
    /* lcc.point(res)=compute_point_3D(lcc,
                                    links_from_pattern_to_volume,
                                    pattern_to_global,
                                    inner); */
    lcc.point(res)=compute_point_3D_v2(lcc,
                                       links_from_pattern_to_volume,
                                       pattern_to_global,
                                       inner);
    // std::cout<<lcc.point()<<std::endl;
  }
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // CMAP_QUERY_REPLACE_GEOMETRY_H
