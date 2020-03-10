// Copyright (c) 2019  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/Conforming_Delaunay_triangulation_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>

#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_data_structure_2.h>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/optional.hpp>

namespace CGAL {

template <typename Gt, typename Vb = Triangulation_vertex_base_3<Gt> >
class Constrained_Delaunay_triangulation_vertex_base_3
  : public Conforming_Delaunay_triangulation_vertex_base_3<Gt, Vb>
{
  using Base = Conforming_Delaunay_triangulation_vertex_base_3<Gt, Vb>;
public:
  bool original_point = false;

  // To get correct vertex type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
    typedef Constrained_Delaunay_triangulation_vertex_base_3 <Gt, Vb3> Other;
  };

  using Base::Base;

};

template <typename Gt, typename Cb = Triangulation_cell_base_3<Gt> >
class Constrained_Delaunay_triangulation_cell_base_3
  : public Cb
{
  using Base = Cb;
public:
  int face_id = -1;

  // To get correct cell type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Constrained_Delaunay_triangulation_cell_base_3 <Gt, Cb3> Other;
  };

  using Base::Base;
};

template <typename T_3>
class Constrained_Delaunay_triangulation_3 : public Conforming_Delaunay_triangulation_3<T_3> {
  using Conforming_Dt = Conforming_Delaunay_triangulation_3<T_3>;
  using Vertex_handle = typename T_3::Vertex_handle;
  using Cell_handle = typename T_3::Cell_handle;
  using Point_3 = typename T_3::Point;
  using Vector_3 = typename T_3::Geom_traits::Vector_3;
  using Locate_type = typename T_3::Locate_type;
  using Geom_traits = typename T_3::Geom_traits;

  using Face_index = int;

  struct CDT_2_types {
    using Projection_traits = Triangulation_2_projection_traits_3<Geom_traits>;
    static_assert(std::is_nothrow_move_constructible<Projection_traits>::value,
                  "move cstr is missing");

    struct Vertex_info {
      Vertex_handle vertex_handle_3d = {};
    };
    struct Face_info {
      int is_outside_the_face = 0;
    };
    using Vb = Triangulation_vertex_base_with_info_2<Vertex_info,
                                                     Projection_traits>;
    using Fb1 = Triangulation_face_base_with_info_2<Face_info,
                                                    Projection_traits>;
    using Fb = Constrained_triangulation_face_base_2<Projection_traits, Fb1>;
    using TDS = Triangulation_data_structure_2<Vb,Fb>;
    using Itag = Exact_predicates_tag;
    using CDT_base =
        Constrained_Delaunay_triangulation_2<Projection_traits, TDS, Itag>;
    using CDT = Constrained_triangulation_plus_2<CDT_base>;

    struct CDT_2_dual_color_map {
      using category = boost::read_write_property_map_tag;
      using reference = int&;
      using value_type = int;
      using key_type = typename CDT::Face_handle;

      friend reference get(CDT_2_dual_color_map, key_type fh) {
        return fh->info().is_outside_the_face;
      }
      friend void put(CDT_2_dual_color_map, key_type fh, int value) {
        fh->info().is_outside_the_face = value;
      }
    };
  }; // CDT_2_types
  using CDT_2 = typename CDT_2_types::CDT;
  using CDT_2_traits = typename CDT_2_types::Projection_traits;
  static_assert(std::is_nothrow_move_constructible<CDT_2>::value,
                "move cstr is missing");
  static_assert(std::is_nothrow_move_assignable<CDT_2>::value,
                "move assignment is missing");

  protected:
    using Constraint_hierarchy = typename Conforming_Dt::Constraint_hierarchy;
    using Constraint_id = typename Constraint_hierarchy::Constraint_id;
    using Subconstraint = typename Constraint_hierarchy::Subconstraint;

    class Insert_in_conflict_visitor {
      Constrained_Delaunay_triangulation_3<T_3> &self;
      typename Conforming_Dt::Insert_in_conflict_visitor conforming_dt_visitor;

    public:
      Insert_in_conflict_visitor(Constrained_Delaunay_triangulation_3 &self)
          : self(self), conforming_dt_visitor(self) {}

      template <class InputIterator>
      void process_cells_in_conflict(InputIterator cell_it, InputIterator end) {
        conforming_dt_visitor.process_cells_in_conflict(cell_it, end);
      }
      void reinsert_vertices(Vertex_handle) const {}
      Vertex_handle replace_vertex(Cell_handle c, int index,
                                   const Point_3 &) const {
        return c->vertex(index);
      }
      void hide_point(Cell_handle, const Point_3 &) const {}
  };

public:
  Vertex_handle insert(const Point_3 &p, Locate_type lt, Cell_handle c,
                               int li, int lj)
  {
    auto v = Conforming_Dt::private_insert(p, lt, c, li, lj, insert_in_conflict_visitor);
    Conforming_Dt::restore_Delaunay();
    return v;
  }

  Vertex_handle insert(const Point_3 &p, Cell_handle start = {}) {
    Locate_type lt;
    int li, lj;

    Cell_handle c = tr.locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
  }

  template <typename Polygon>
  CGAL_CPP20_REQUIRE_CLAUSE(
      std::ranges::common_range<Polygon>
      && (std::is_convertible_v<std::ranges::range_value_t<Polygon>, Point_3>))
  boost::optional<Face_index> insert_constrained_polygon(Polygon&& polygon) {
    std::vector<Vertex_handle> handles;
    handles.reserve(polygon.size());
    boost::optional<Cell_handle> hint;
    for(const auto& p : polygon) {
      handles.push_back(this->insert(p, hint.value_or(Cell_handle{})));
      hint = handles.back()->cell();
    }
    return insert_constrained_face(std::move(handles));
  }

  template <typename Vertex_handles>
  CGAL_CPP20_REQUIRE_CLAUSE(
      std::ranges::common_range<Vertex_handles>
      && (std::is_convertible_v<std::ranges::range_value_t<Vertex_handles>, Vertex_handle>))
  boost::optional<Face_index> insert_constrained_face(Vertex_handles&& vertex_handles) {
    using std::begin;
    using std::endl;
    const auto first_it = begin(vertex_handles);
    const auto vend =  end(vertex_handles);
    const auto size = std::distance(first_it, vend);
    if(size < 2) return {};
    if(size == 2) {
      this->insert_constrained_edge(*first_it, *std::next(first_it));
      return {};
    }
    const auto accumulated_normal = [&] {
      const auto last_it = std::next(first_it, size - 1);
      const auto &last_point = tr.point(*last_it);

      auto &&traits = tr.geom_traits();
      auto &&cross_product = traits.construct_cross_product_vector_3_object();
      auto &&vector = traits.construct_vector_3_object();
      auto &&sum_vector = traits.construct_sum_of_vectors_3_object();

      Vector_3 accumulated_normal = vector(CGAL::NULL_VECTOR);
      for (auto vit = first_it, next_it = std::next(first_it);
           next_it != last_it; ++vit, ++next_it) {
        accumulated_normal =
            sum_vector(accumulated_normal,
                       cross_product(vector(last_point, tr.point(*vit)),
                                     vector(last_point, tr.point(*next_it))));
      }
      return accumulated_normal;
    }();

    CGAL::Circulator_from_container<std::remove_reference_t<Vertex_handles>>
        circ{&vertex_handles}, next_circ{circ}, circ_end{circ};

    { // create and fill the 2D triangulation
      faces_to_preserve.emplace_back(CDT_2_traits{accumulated_normal});
      CDT_2& cdt_2 = faces_to_preserve.back();
      boost::optional<typename CDT_2::Vertex_handle> first;
      boost::optional<typename CDT_2::Vertex_handle> previous;
      do {
        auto vh_2 = cdt_2.insert(tr.point(*circ));
        vh_2->info().vertex_handle_3d = *circ;
        if(!first) first = vh_2;
        if(previous) {
          cdt_2.insert_constraint(*previous, vh_2);
          std::cerr << "insert constraint ("
                    << tr.point((*previous)->info().vertex_handle_3d)
                    << " , "
                    << tr.point((vh_2)->info().vertex_handle_3d)
                    << '\n';
        }
        previous = vh_2;
        ++circ;
        if(circ == circ_end) {
          cdt_2.insert_constraint(vh_2, *first);
          std::cerr << "insert constraint ("
                    << tr.point((vh_2)->info().vertex_handle_3d)
                    << " , "
                    << tr.point((*first)->info().vertex_handle_3d)
                    << '\n';
          break;
        }
      } while (circ != circ_end);
      { // Now, use BGL BFS algorithm to mark the faces reachable from
        // an infinite face as `is_outside_the_face`.
        // I use the fact that `white_color()` evaluates to 0, and
        // `black_color()` to 4 (and thus to a true Boolean).
        auto cdt_2_dual_graph = dual(cdt_2.tds());
        using pred_type = bool (*)(typename CDT_2::Edge);
        pred_type no_constrained_edge = [](typename CDT_2::Edge edge) mutable {
          return !edge.first->is_constrained(edge.second);
        };
        boost::filtered_graph<decltype(cdt_2_dual_graph),
                              pred_type>
            dual(cdt_2_dual_graph, no_constrained_edge);
        using Color_map = typename CDT_2_types::CDT_2_dual_color_map;
        boost::breadth_first_search(dual, cdt_2.infinite_vertex()->face(),
                                    boost::color_map(Color_map()));
      } // end of Boost BFS
#if CGAL_DEBUG_CDT_3
      int counter = 0;
      for(auto fh: cdt_2.finite_face_handles()) {
        if(!fh->info().is_outside_the_face) ++counter;
      }
      std::cerr << counter << " triangles(s) in the face\n";
#endif // CGAL_DEBUG_CDT_3
      for(auto fh: cdt_2.all_face_handles()) {
        if(fh->info().is_outside_the_face) continue;
        const auto v0 = fh->vertex(0)->info().vertex_handle_3d;
        const auto v1 = fh->vertex(1)->info().vertex_handle_3d;
        const auto v2 = fh->vertex(2)->info().vertex_handle_3d;
        Cell_handle c;
        int i, j, k;
        if(!tr.is_facet(v0, v1, v2, c, i, j, k)) {
          std::cerr << "Missing triangle: \n";
          std::cerr << "4"
                    << " " << tr.point(v0)
                    << " " << tr.point(v1)
                    << " " << tr.point(v2)
                    << " " << tr.point(v0) << '\n';
        }
      }
    } // end of the construction of the CDT_2
    CGAL_assertion( circ == circ_end && next_circ == circ);
    do {
      ++next_circ;
      auto constraint_id = this->insert_constrained_edge(*circ,*next_circ);
    } while(++circ != circ_end);

    return faces_to_preserve.size();
  }

  /// @{
  /// remove functions cannot be called
  void remove(Vertex_handle) = delete;
  void remove_cluster() = delete;
  /// @}

protected:
private:
  T_3 &tr = *this;
  Conforming_Dt &conforming_dt = *this;
  Insert_in_conflict_visitor insert_in_conflict_visitor = {*this};

  std::vector<CDT_2> faces_to_preserve;
};

} // end CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
