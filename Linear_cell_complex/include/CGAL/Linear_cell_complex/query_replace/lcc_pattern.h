// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_PATTERN_H
#define LCC_PATTERN_H

#include "lcc_barycentric_coord.h"
#include "cmap_signature.h"

#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates_3.h>
//#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/optimal_bounding_box.h>
//#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
//#include <CGAL/Surface_mesh.h>

#include <cassert>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace CGAL::internal
{
template<class LCC>
class Pattern_substituer;

template<class LCC>
using Dart_mapping=std::unordered_map<typename LCC::Dart_descriptor,
                                      typename LCC::Dart_descriptor>;

///////////////////////////////////////////////////////////////////////////////
/// type==1 for face, 2 for surface, 3 for volume
template<class LCC, unsigned int type>
class Pattern;
////////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern<LCC, 1> // Face pattern
{
  friend class Pattern_substituer<LCC>;

  using Dart_descriptor=typename LCC::Dart_descriptor;
  using size_type=typename LCC::size_type;
  using Point=typename LCC::Point;
  using Vector=typename LCC::Vector;
public:
  Pattern(): m_mark_to_preserve(LCC::INVALID_MARK)
  {}

  LCC& lcc()
  { return m_lcc; }

  size_type reserve_mark_to_preserve()
  {
    if(m_mark_to_preserve==LCC::INVALID_MARK)
    { m_mark_to_preserve=m_lcc.get_new_mark(); }
    return m_mark_to_preserve;
  }

  size_type mark_to_preserve() const
  { return m_mark_to_preserve; }

  std::vector<Barycentric_coord<LCC, 1>>& barycentric_coords()
  { return m_barycentric_coords; }

  /// Compute the barycentric coordinates of new vertices (vertices with all
  /// darts not marked surviving).
  /// These coordinates are stored in m_barycentric_coords vector.
  /// A dart is called surviving if it is shared before and after the
  /// "transformation".
  void compute_barycentric_coord(size_type surviving)
  {
    /// For replace patterns, a vertex having all its darts non surviving is
    /// a new vertex (it didn't exist in the query pattern). Thus, a vertex
    /// having at least one surviving dart is an old vertex (it existed in the
    /// query pattern).
    size_type old_vertex=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, surviving) && !m_lcc.is_marked(it, old_vertex))
      { m_lcc.template mark_cell<0>(it, old_vertex); }
    }

    typename LCC::Point bary2=CGAL::ORIGIN;
    std::vector<Dart_descriptor> surviving_darts;
    std::vector<Dart_descriptor> new_vertices;
    std::size_t nb1=0;
    auto vertex_treated=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, surviving))
      { surviving_darts.push_back(it); }
      if(!m_lcc.is_marked(it, vertex_treated))
      { // TODO we can improve the traversal of cells (regroup, use basic it...)
        m_lcc.template mark_cell<0>(it, vertex_treated);
        if(m_lcc.is_marked(it, old_vertex))
        { // Here is incident to a surviving vertex
          const Point& p=m_lcc.point(it);
          bary2=Point(bary2.x()+p.x(), bary2.y()+p.y(), bary2.z()+p.z());
          ++nb1;
          m_lcc.template unmark_cell<0>(it, old_vertex); // to speedup the freemark at end
        }
        else
        { // Here it is incident to a new vertex
          new_vertices.push_back(it);
        }
      }
    }
    m_lcc.free_mark(old_vertex);
    m_lcc.free_mark(vertex_treated);
    assert(nb1>0);
    bary2=Point(bary2.x()/nb1, bary2.y()/nb1, bary2.z()/nb1);

    nb1=0;
    m_barycentric_coords.resize(new_vertices.size());
    for(auto it: new_vertices)
    {
      m_barycentric_coords[nb1].m_coords.reserve(surviving_darts.size());
      m_barycentric_coords[nb1].m_dart=it;
      ++nb1;
    }

    double alpha, beta, gamma;
    for(auto& itd: surviving_darts)
    {
      const Point& p0=m_lcc.point(itd);
      const Point& p1=m_lcc.point(m_lcc.other_extremity(itd));
      nb1=0;
      for(auto& it: new_vertices)
      { // TODO we can test if tetra (p0, p1, p2, bary3) before to enter in the loop
        if(compute_alpha_beta_gamma_of_point(p0, p1, bary2,
                                              m_lcc.point(it),
                                              alpha, beta, gamma))
        {
          m_barycentric_coords[nb1].m_coords.push_back
              (std::make_tuple(itd, alpha, beta, gamma));
          ++nb1;
        }
      }
    }
  }

  /// Compute the barycentric coordinates of inner vertices (vertices having
  /// no dart 2-free). These coordinates are stored in m_barycentric_coords
  /// vector.
  void compute_barycentric_coord()
  {
    size_type surviving=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.template is_free<2>(it))
      { m_lcc.mark(it, surviving); }
    }

    compute_barycentric_coord(surviving);
    m_lcc.free_mark(surviving);
  }

  void display()
  {
    for(auto& it: m_barycentric_coords)
    { it.display(m_lcc); }
  }

protected:
  LCC m_lcc;
  typename LCC::size_type m_mark_to_preserve;
  /// For each inner point, its barycentric coordinates for each external point
  std::vector<Barycentric_coord<LCC, 1>> m_barycentric_coords;
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern<LCC, 2> // Surfacic pattern
{
  friend class Pattern_substituer<LCC>;

  using Dart_descriptor=typename LCC::Dart_descriptor;
  using size_type=typename LCC::size_type;
  using Point=typename LCC::Point;
  using Vector=typename LCC::Vector;

public:
  Pattern():
    m_mark_faceborder(m_lcc.get_new_mark()),
    m_mark_to_preserve(LCC::INVALID_MARK)
  {}

  LCC& lcc()
  { return m_lcc; }

  size_type reserve_mark_to_preserve()
  {
    if(m_mark_to_preserve==LCC::INVALID_MARK)
    { m_mark_to_preserve=m_lcc.get_new_mark(); }
    return m_mark_to_preserve;
  }

  size_type mark_to_preserve() const
  { return m_mark_to_preserve; }

  // same barycentric coords than for face => 1
  std::vector<Barycentric_coord<LCC, 1>>& barycentric_coords()
  { return m_barycentric_coords; }

  /// Compute the barycentric coordinates of new vertices (vertices with all
  /// darts not marked surviving).
  /// These coordinates are stored in m_barycentric_coords vector.
  /// A dart is called surviving if it is shared before and after the
  /// "transformation".
  void compute_barycentric_coord(size_type surviving)
  {
    /// For replace patterns, a vertex having all its darts non surviving is
    /// a new vertex (it didn't exist in the query pattern). Thus, a vertex
    /// having at least one surviving dart is an old vertex (it existed in the
    /// query pattern).
    size_type old_vertex=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, surviving) && !m_lcc.is_marked(it, old_vertex))
      { m_lcc.template mark_cell<0>(it, old_vertex); }
    }

    /*
    // Compute mean value coordinates for all interior points.
    std::cout << std::endl << "mean value coordinates (interior): " << std::endl << std::endl;

    for (const auto& query : interior_points) {
      coordinates.clear();
      CGAL::Barycentric_coordinates::mean_value_coordinates_2(
          star_shaped, query, std::back_inserter(coordinates), policy);

      // Output mean value coordinates.
      for (std::size_t i = 0; i < coordinates.size() - 1; ++i) {
        std::cout << coordinates[i] << ", ";
      }
      std::cout << coordinates[coordinates.size() - 1] << std::endl;
    }
    */

    typename LCC::Point bary2=CGAL::ORIGIN;
    std::vector<Dart_descriptor> surviving_darts;
    std::vector<Dart_descriptor> new_vertices;
    std::size_t nb1=0, new_index=0;
    Dart_descriptor cur, other;
    auto treated=m_lcc.get_new_mark();
    auto vertex_treated=m_lcc.get_new_mark();
    std::queue<Dart_descriptor> to_treat;
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(!m_lcc.is_marked(it, treated))
      {
        surviving_darts.clear();
        new_vertices.clear();
        nb1=0;
        bary2=CGAL::ORIGIN;

        // Here we iterate through the cc of faces inside a same cycle of edges
        // marked by m_mark_faceborder
        to_treat.push(it);
        m_lcc.mark(it, treated);
        while(!to_treat.empty())
        {
          cur=to_treat.front();
          to_treat.pop();

          if(m_lcc.is_marked(cur, m_mark_faceborder))
          {  surviving_darts.push_back(cur); }

          if(!m_lcc.is_marked(cur, vertex_treated))
          { // TODO we can improve the traversal of cells (regroup, use basic it...)
            m_lcc.template mark_cell<0>(cur, vertex_treated);
            if(m_lcc.is_marked(cur, old_vertex))
            { // Here is incident to a surviving vertex
              const Point& p=m_lcc.point(cur);
              bary2=Point(bary2.x()+p.x(), bary2.y()+p.y(), bary2.z()+p.z());
              ++nb1;
            }
            else
            { // Here it is incident to a new vertex
              new_vertices.push_back(cur);
            }
          }

          other=m_lcc.template beta<1>(cur);
          if(!m_lcc.is_marked(other, treated))
          {
            to_treat.push(other);
            m_lcc.mark(other, treated);
          }

          if(!m_lcc.is_marked(cur, m_mark_faceborder))
          {
            other=m_lcc.template beta<2>(cur);
            if(!m_lcc.is_marked(other, treated))
            {
              to_treat.push(other);
              m_lcc.mark(other, treated);
            }
          }
        }

        assert(nb1>0);
        // Now compute the barycentric coordinates of inner vertices
        if(new_vertices.size()>0)
        {
          bary2=Point(bary2.x()/nb1, bary2.y()/nb1, bary2.z()/nb1);

          new_index=m_barycentric_coords.size();
          nb1=new_index;
          m_barycentric_coords.resize(m_barycentric_coords.size()+
                                      new_vertices.size());
          for(auto iti: new_vertices)
          {
            m_barycentric_coords[nb1].m_coords.reserve(surviving_darts.size());
            m_barycentric_coords[nb1].m_dart=iti;
            ++nb1;
            if(m_lcc.is_marked(iti, vertex_treated))
            { m_lcc.template unmark_cell<0>(iti, vertex_treated);}
          }

          double alpha, beta, gamma;
          for(auto& itd: surviving_darts)
          {
            const Point& p0=m_lcc.point(itd);
            const Point& p1=m_lcc.point(m_lcc.other_extremity(itd));
            nb1=new_index;
            for(auto& iti: new_vertices)
            {
              if(compute_alpha_beta_gamma_of_point(p0, p1, bary2,
                                                    m_lcc.point(iti),
                                                   alpha, beta, gamma))
              {
                m_barycentric_coords[nb1].m_coords.push_back
                    (std::make_tuple(itd, alpha, beta, gamma));
                ++nb1;
              }
            }
            if(m_lcc.is_marked(itd, vertex_treated))
            { m_lcc.template unmark_cell<0>(itd, vertex_treated);}
          }
        }
        else
        {
          for(auto& itd: surviving_darts)
          {
            if(m_lcc.is_marked(itd, vertex_treated))
            { m_lcc.template unmark_cell<0>(itd, vertex_treated);}
          }
        }
        assert(m_lcc.is_whole_map_unmarked(vertex_treated));
      }
    }
    assert(m_lcc.is_whole_map_marked(treated));
    m_lcc.free_mark(treated);
    m_lcc.free_mark(vertex_treated);
 }

  /// Compute the barycentric coordinates of inner vertices (vertices having
  /// no dart marked with m_mark_faceborder). These coordinates are stored
  /// in m_barycentric_coords vector.
  void compute_barycentric_coord()
  {
    size_type surviving=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, m_mark_faceborder))
      { m_lcc.mark(it, surviving); }
    }

    compute_barycentric_coord(surviving);
    m_lcc.free_mark(surviving);
  }

  void display()
  {
    for(auto& it: m_barycentric_coords)
    { it.display(m_lcc); }
  }

  typename LCC::size_type mark_faceborder() const
  { return m_mark_faceborder; }

protected:
  LCC m_lcc;
  typename LCC::size_type m_mark_faceborder;
  typename LCC::size_type m_mark_to_preserve;
  /// For each inner point, its barycentric coordinates for each external point
  std::vector<Barycentric_coord<LCC, 1>> m_barycentric_coords;
};
////////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern<LCC, 3> // Volumic pattern
{
  friend class Pattern_substituer<LCC>;

  using Dart_descriptor=typename LCC::Dart_descriptor;
  using size_type=typename LCC::size_type;
  using Point=typename LCC::Point;
  using Vector=typename LCC::Vector;
public:
  Pattern(): m_mark_to_preserve(LCC::INVALID_MARK)
  {}

  LCC& lcc()
  { return m_lcc; }

  size_type reserve_mark_to_preserve()
  {
    if(m_mark_to_preserve==LCC::INVALID_MARK)
    { m_mark_to_preserve=m_lcc.get_new_mark(); }
    return m_mark_to_preserve;
  }

  size_type mark_to_preserve() const
  { return m_mark_to_preserve; }

  std::vector<Barycentric_coord<LCC, 3>>& barycentric_coords()
  { return m_barycentric_coords; }

  /// Compute the barycentric coordinates of new vertices (vertices with all
  /// darts not marked surviving).
  /// These coordinates are stored in m_barycentric_coords vector.
  /// A dart is called surviving if it is shared before and after the
  /// "transformation".
  void compute_barycentric_coord(size_type surviving)
  {
    /* namespace PMP = CGAL::Polygon_mesh_processing;
    bounding_mesh.clear();
    delete mv;

    /// 1) Compute the oriented bounding box of surviving points
    std::unordered_set<typename LCC::Point> points;
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
        it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, surviving))
      { points.insert(m_lcc.point(it)); }
    }
    std::array<Point, 8> obb_points;
    CGAL::oriented_bounding_box(points, obb_points,
                                CGAL::parameters::use_convex_hull(true));
    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                          obb_points[4], obb_points[5], obb_points[6], obb_points[7],
                          bounding_mesh);
    PMP::triangulate_faces(faces(bounding_mesh), bounding_mesh);

    / * std::cout<<"OBB: ";
    for(std::size_t i=0; i<8; ++i) {std::cout<<obb_points[i]<<"  ";}
    std::cout<<std::endl; * /

    /// 2) Compute the mean value coordinates_3 for this bounding box
    mv=new CGAL::Barycentric_coordinates::Mean_value_coordinates_3
      <CGAL::Surface_mesh<typename LCC::Point>,
        CGAL::Exact_predicates_inexact_constructions_kernel>
        (bounding_mesh); //,  CGAL::Barycentric_coordinates::Computation_policy_3::FAST);
*/

    size_type old_vertex=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, surviving) && !m_lcc.is_marked(it, old_vertex))
      { m_lcc.template mark_cell<0>(it, old_vertex); }
    }

    typename LCC::Point bary3=CGAL::ORIGIN;
    std::vector<Dart_descriptor> surviving_darts;
    std::vector<Dart_descriptor> new_vertices;
    std::size_t nb1=0;
    auto vertex_treated=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.is_marked(it, surviving))
      { surviving_darts.push_back(it); }
      if(!m_lcc.is_marked(it, vertex_treated))
      { // TODO we can improve the traversal of cells (regroup, use basic it...)
        m_lcc.template mark_cell<0>(it, vertex_treated);
        if(m_lcc.is_marked(it, old_vertex))
        { // Here is incident to a surviving vertex
          const Point& p=m_lcc.point(it);
          bary3=Point(bary3.x()+p.x(), bary3.y()+p.y(), bary3.z()+p.z());
          ++nb1;
          m_lcc.template unmark_cell<0>(it, old_vertex); // to speedup the freemark at end
        }
        else
        { // Here it is incident to a new vertex
          new_vertices.push_back(it);
        }
      }
    }
    assert(m_lcc.is_whole_map_unmarked(old_vertex));
    assert(m_lcc.is_whole_map_marked(vertex_treated));
    m_lcc.free_mark(old_vertex);
    m_lcc.free_mark(vertex_treated);
    assert(nb1>0);
    bary3=Point(bary3.x()/nb1, bary3.y()/nb1, bary3.z()/nb1);

    nb1=0;
    m_barycentric_coords.resize(new_vertices.size());
    for(auto it: new_vertices)
    {
      m_barycentric_coords[nb1].m_coords.reserve(surviving_darts.size());
      m_barycentric_coords[nb1].m_dart=it;
      ++nb1;
    }

    //double alpha, beta, gamma, delta;
    for(auto& itd: surviving_darts)
    {
      const Point& p0=m_lcc.point(itd);
      const Point& p1=m_lcc.point(m_lcc.other_extremity(itd));
      const Point p2=m_lcc.template barycenter<2>(itd);
      nb1=0;
      typename CGAL::Kernel_traits<Point>::Kernel::Tetrahedron_3
          t(p0, p1, p2, bary3);
      if(!t.is_degenerate())
      {
        std::vector<double> coords; coords.reserve(4);
        for(auto& it: new_vertices)
        {
          coords.clear();
          CGAL::Barycentric_coordinates::tetrahedron_coordinates_3
              (p0, p1, p2, bary3, m_lcc.point(it), std::back_inserter(coords));
          /*if(compute_alpha_beta_gamma_delta_of_point(p0, p1, p2, bary3,
                                                   m_lcc.point(it),
                                                   alpha, beta, gamma, delta))*/
          {
            m_barycentric_coords[nb1].m_coords.push_back
                //(std::make_tuple(itd, alpha, beta, gamma, delta));
                (std::make_tuple(itd, coords[0], coords[1], coords[2], coords[3]));
            ++nb1;
          }
        }
      }
    }
  }

  /// Compute the barycentric coordinates of inner vertices (vertices having
  /// no dart 3-free). These coordinates are stored in m_barycentric_coords
  /// vector. Valid only for query replace of one cell!!
  void compute_barycentric_coord()
  {
    size_type surviving=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end();
         it!=itend; ++it)
    {
      if(m_lcc.template is_free<3>(it))
      { m_lcc.mark(it, surviving); }
    }

    compute_barycentric_coord(surviving);
    m_lcc.free_mark(surviving);
  }

  void display()
  {
    for(auto& it: m_barycentric_coords)
    { it.display(m_lcc); }
  }

protected:
  LCC m_lcc;
  size_type m_mark_to_preserve;
  /// For each inner point, its barycentric coordinates for each external point
  std::vector<Barycentric_coord<LCC, 3>> m_barycentric_coords;

  /// Use CGAL (unofficial) package Barycentric_coordinates_3
  /* CGAL::Surface_mesh<typename LCC::Point> bounding_mesh;

public:
  CGAL::Barycentric_coordinates::Mean_value_coordinates_3
  <CGAL::Surface_mesh<typename LCC::Point>,
   CGAL::Exact_predicates_inexact_constructions_kernel>* mv=nullptr; */
};
///////////////////////////////////////////////////////////////////////////////
enum Pattern_type
{
  APPLY,
  DELETE,
  EXTEND,
  REPLACE
};
///////////////////////////////////////////////////////////////////////////////
// Class to store all info of a pattern for the query/replace of multiple cells.
// Type is 1 for face, 3 for volume
template<class LCC, unsigned int type>
class Pattern_for_multiple_cells
{
public:
  Pattern<LCC, type> replace_pattern;

  Signature query_word; // We only need the word for the query, no need to keep the whole lcc

  // keys are ids of darts to keep (id is the label in the signature),
  // values are their corresponding dart on replace pattern
  std::unordered_map<MyInt, typename LCC::Dart_descriptor> darts_to_keep;

  // for all automorphisms in the query pattern, store the id of the dart
  std::vector<MyInt> automorphisms;

  Pattern_type pattern_type; // REPLACE, ADD, REMOVE or APPLY

  std::string function_to_apply;
};
////////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_PATTERN_H
