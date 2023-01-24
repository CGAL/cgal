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

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_3.h>

#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_data_structure_2.h>

#include <CGAL/Mesh_3/io_signature.h>

#include <CGAL/Conforming_Delaunay_triangulation_3.h>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/optional.hpp>
#include <boost/dynamic_bitset.hpp>

#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/small_vector.hpp>

#include <ranges>

namespace CGAL {

using CDT_3_face_index = int;

template <typename Gt, typename Vb = Triangulation_vertex_base_3<Gt> >
class Constrained_Delaunay_triangulation_vertex_base_3
  : public Conforming_Delaunay_triangulation_vertex_base_3<Gt, Vb>
{
  using Base = Conforming_Delaunay_triangulation_vertex_base_3<Gt, Vb>;
public:
  bool original_point = false;  // currently not used

  // To get correct vertex type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
    typedef Constrained_Delaunay_triangulation_vertex_base_3 <Gt, Vb3> Other;
  };

  using Base::Base;

  static std::string io_signature() {
    return Get_io_signature<Base>()();
  }
};

template <typename Gt, typename Cb = Triangulation_cell_base_3<Gt> >
class Constrained_Delaunay_triangulation_cell_base_3
  : public Cb
{
  using Base = Cb;
  std::array<CDT_3_face_index, 4> face_id = { -1, -1, -1, -1 };
  std::array<void*, 4> facet_2d = {nullptr, nullptr, nullptr, nullptr};

public:
  // To get correct cell type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Constrained_Delaunay_triangulation_cell_base_3 <Gt, Cb3> Other;
  };

  // Constructor
  using Base::Base;

  bool is_facet_constrained(int i) const { return face_id[unsigned(i)] >= 0; }

  template <typename Facet_handle>
  void set_facet_constraint(int i, CDT_3_face_index face_id,
                            Facet_handle facet_2d)
  {
    this->face_id[unsigned(i)] = face_id;
    this->facet_2d[unsigned(i)] = static_cast<void*>(std::addressof(*facet_2d));
  }

  CDT_3_face_index face_constraint_index(int i) const {
    return face_id[unsigned(i)];
  }

  template <typename CDT_2>
  auto face_2 (const CDT_2& cdt, int i) const {
    using Face = typename CDT_2::Face;
    auto ptr = static_cast<Face*>(facet_2d[unsigned(i)]);
    return cdt.tds().faces().iterator_to(*ptr);
  }

  static std::string io_signature() {
    return Get_io_signature<Base>()() + "+(" + Get_io_signature<int>()()
      + ")[4]";
  }

  friend std::ostream&
  operator<<(std::ostream& os,
             const Constrained_Delaunay_triangulation_cell_base_3& c)
  {
    os << static_cast<const Base&>(c);
    for( unsigned li = 0; li < 4; ++li ) {
      if(is_ascii(os)) {
        os << " " << c.face_id[li];
      } else {
        CGAL::write(os, c.face_id[li]);
      }
    }
    return os;
  }
  friend std::istream&
  operator>>(std::istream& is,
             Constrained_Delaunay_triangulation_cell_base_3& c)
  {
    is >> static_cast<Base&>(c);
    if(!is) return is;
    for( int li = 0; li < 4; ++li ) {
      int i;
      if(is_ascii(is)) {
        is >> i;
      } else {
        CGAL::read(is, i);
      }
      if(!is) return is;
      c.face_id[li] = i;
    }
    return is;
  }
};

template <typename T_3>
class Constrained_Delaunay_triangulation_3 : public Conforming_Delaunay_triangulation_3<T_3> {
public:
  using Conforming_Dt = Conforming_Delaunay_triangulation_3<T_3>;
  using Vertex_handle = typename T_3::Vertex_handle;
  using Cell_handle = typename T_3::Cell_handle;
  using Edge = typename T_3::Edge;
  using Point_3 = typename T_3::Point;
  using Vector_3 = typename T_3::Geom_traits::Vector_3;
  using Locate_type = typename T_3::Locate_type;
  using Geom_traits = typename T_3::Geom_traits;

  using Face_index = CDT_3_face_index;

  static std::string io_signature() {
  return Get_io_signature<Conforming_Dt>()();
}

private:
  struct CDT_2_types {
    using Projection_traits = Projection_traits_3<Geom_traits>;
    static_assert(std::is_nothrow_move_constructible<Projection_traits>::value,
                  "move cstr is missing");

    struct Vertex_info {
      Vertex_handle vertex_handle_3d = {};
    };

    using Color_value_type = std::int8_t;
    struct Face_info {
      Color_value_type is_outside_the_face = 0;
      Color_value_type is_in_region = 0;
      std::bitset<3> is_edge_also_in_3d_triangulation;
      bool missing_subface = true;
    };
    using Vb1 = Triangulation_vertex_base_with_info_2<Vertex_info,
                                                      Projection_traits>;
    using Vb = Base_with_time_stamp<Vb1>;
    using Fb1 = Triangulation_face_base_with_info_2<Face_info,
                                                    Projection_traits>;
    using Fb = Constrained_triangulation_face_base_2<Projection_traits, Fb1>;
    using TDS = Triangulation_data_structure_2<Vb,Fb>;
    using Itag = No_constraint_intersection_tag;
    using CDT_base =
        Constrained_Delaunay_triangulation_2<Projection_traits, TDS, Itag>;
    using CDT = CDT_base;

    template <Color_value_type Face_info::* member_ptr>
    struct CDT_2_dual_color_map {
      using category = boost::read_write_property_map_tag;
      using reference = Color_value_type&;
      using value_type = Color_value_type;
      using key_type = typename CDT::Face_handle;

      friend reference get(CDT_2_dual_color_map, key_type fh) {
        return fh->info().*member_ptr;
      }
      friend void put(CDT_2_dual_color_map, key_type fh, value_type value) {
        fh->info().*member_ptr = value;
      }
    };
    using Color_map_is_outside_the_face =
        CDT_2_dual_color_map<&Face_info::is_outside_the_face>;
    using Color_map_is_in_region =
        CDT_2_dual_color_map<&Face_info::is_in_region>;
  }; // CDT_2_types
  using CDT_2 = typename CDT_2_types::CDT;
  using CDT_2_traits = typename CDT_2_types::Projection_traits;
  using CDT_2_face_handle = typename CDT_2::Face_handle;
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
      void process_cells_in_conflict(const InputIterator cell_it_begin, const InputIterator end) {
        CGAL_assertion(self.dimension() >= 2);
        const int first_li = self.dimension() == 2 ? 3 : 0;
        for(auto cell_it = cell_it_begin; cell_it != end; ++cell_it) {
          auto c = *cell_it;
          for(int li = first_li; li < 4; ++li) {
            if(c->is_facet_constrained(li)) {
              const auto face_id = static_cast<std::size_t>(c->face_constraint_index(li));
              self.face_constraint_misses_subfaces.set(face_id);
              auto fh_2 = c->face_2(self.face_cdt_2[face_id], li);
#if CGAL_DEBUG_CDT_3
              std::cerr << "Add missing triangle (from visitor): \n";
              self.write_2d_triangle(std::cerr, fh_2);
#endif // CGAL_DEBUG_CDT_3

              fh_2->info().missing_subface = true;
            }
          }
        }
        conforming_dt_visitor.process_cells_in_conflict(cell_it_begin, end);
      }
      void after_insertion(Vertex_handle) {

      }

      void reinsert_vertices(Vertex_handle v) {
        after_insertion(v);
      }
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
    auto v = Conforming_Dt::insert_impl(p, lt, c, li, lj, insert_in_conflict_visitor);
    Conforming_Dt::restore_Delaunay(insert_in_conflict_visitor);
    return v;
  }

  Vertex_handle insert(const Point_3 &p, Cell_handle start = {}) {
    Locate_type lt;
    int li, lj;

    Cell_handle c = tr.locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
  }

  Constraint_id insert_constrained_edge(Vertex_handle va, Vertex_handle vb)
  {
    return this->insert_constrained_edge_impl(va, vb, insert_in_conflict_visitor);
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
    CGAL::Circulator_from_container<std::remove_reference_t<Vertex_handles>> circ{&vertex_handles};
    const auto circ_end{circ};
    auto& border = this->face_border.emplace_back();
    do {
      const auto va = *circ;
      ++circ;
      const auto vb = *circ;
      const auto c_id = this->constraint_from_extremities(va, vb);
      if(c_id != Constraint_id{}) {
        const bool constraint_c_id_is_reversed = true;
        border.push_back(Face_edge{c_id, constraint_c_id_is_reversed});
      } else {
        const auto c_id = this->insert_constrained_edge(va, vb);
        CGAL_assertion(c_id != Constraint_id{});
        border.push_back(Face_edge{c_id});
      }
    } while(circ != circ_end);

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
      if (accumulated_normal.z() < 0 ||
          (accumulated_normal.z() == 0 && accumulated_normal.y() > 0) ||
          (accumulated_normal.z() == 0 && accumulated_normal.y() == 0 &&
           accumulated_normal.x() > 0)
          )
      {
        accumulated_normal = - accumulated_normal;
      }
      return accumulated_normal;
    }();

    face_cdt_2.emplace_back(CDT_2_traits{accumulated_normal});
    face_constraint_misses_subfaces.resize(face_cdt_2.size());
    const auto polygon_contraint_id = static_cast<CDT_3_face_index>(face_cdt_2.size() - 1);

    // search_for_missing_subfaces(polygon_contraint_id);
    // restore_constrained_Delaunay();

    return polygon_contraint_id;
  }

private:
  void fill_cdt_2(CDT_3_face_index polygon_contraint_id)
  {
    CDT_2& cdt_2 = face_cdt_2[polygon_contraint_id];

    const auto handles = [this, polygon_contraint_id]() {
      std::vector<Vertex_handle> handles;
      for(const auto& face_edge : this->face_border[polygon_contraint_id]) {
        const auto c_id = face_edge.constraint_id;
        const bool reversed = face_edge.is_reverse;
        const auto v_begin = this->constraint_hierarchy.vertices_in_constraint_begin(c_id);
        const auto v_end = this->constraint_hierarchy.vertices_in_constraint_end(c_id);
        CGAL_assertion(std::distance(v_begin, v_end) >= 2);
        auto va = *v_begin;
        auto vb_it = v_end;
        --vb_it;
        auto vb = *vb_it;
        if(reversed) {
          using std::swap;
          swap(va, vb);
        }
        if(handles.empty()) {
          handles.push_back(va);
        } else {
          CGAL_assertion(handles.back() == va);
        }
        handles.push_back(vb);
      }
      CGAL_assertion(handles.front() == handles.back());
      handles.pop_back();
      return handles;
    }();

    CGAL::Circulator_from_container circ{&handles};
    const auto circ_end{circ};
    { // create and fill the 2D triangulation
      const auto first_2d  = cdt_2.insert(tr.point(*circ));
      first_2d->info().vertex_handle_3d = *circ;
      auto previous_2d = first_2d;
      do {
        const auto va = *circ;
        CGAL_assertion(previous_2d->info().vertex_handle_3d == va);
        ++circ;
        const auto vb = *circ;
        const auto c_id = this->constraint_from_extremities(va, vb);
        if(c_id != Constraint_id{}) {
          auto vit = this->constraint_hierarchy.vertices_in_constraint_begin(c_id);
          auto v_end = this->constraint_hierarchy.vertices_in_constraint_end(c_id);
          CGAL_assertion_code(const auto constraint_size = std::distance(vit, v_end);)
          if(vit != v_end) {
            const bool constraint_c_id_is_reversed = (*vit != va);
            CGAL_assertion(*vit == (constraint_c_id_is_reversed ? vb : va));
            if(++vit != v_end && vit != --v_end) {
              CGAL_assertion(constraint_size == std::distance(vit, v_end) + 2);
              CGAL_assertion(*v_end == (constraint_c_id_is_reversed ? va : vb));
              if(constraint_c_id_is_reversed) {
                using std::swap;
                swap(vit, v_end);
                --vit;
                --v_end;
              };
              while(vit != v_end) {
                auto vh_2d = cdt_2.insert(tr.point(*vit));
                vh_2d->info().vertex_handle_3d = *vit;
                std::cerr << "cdt_2.insert_constraint ("
                          << tr.point(previous_2d->info().vertex_handle_3d)
                          << " , "
                          << tr.point(vh_2d->info().vertex_handle_3d)
                          << ")\n";
                cdt_2.insert_constraint(previous_2d, vh_2d);
                previous_2d = vh_2d;
                if(constraint_c_id_is_reversed) {
                  --vit;
                } else {
                  ++vit;
                };
              }
            }
          }
        }

        auto vh_2d = circ == circ_end ? first_2d : cdt_2.insert(tr.point(vb));
        if(circ != circ_end) {
          vh_2d->info().vertex_handle_3d = vb;
        }
        std::cerr << "cdt_2.insert_constraint ("
                  << tr.point(previous_2d->info().vertex_handle_3d)
                  << " , "
                  << tr.point(vh_2d->info().vertex_handle_3d)
                  << ")\n";
        cdt_2.insert_constraint(previous_2d, vh_2d);
        previous_2d = vh_2d;
      } while (circ != circ_end);
      const auto cdt_2_dual_graph = dual(cdt_2.tds());
      { // Now, use BGL BFS algorithm to mark the faces reachable from
        // an infinite face as `is_outside_the_face`.
        // I use the fact that `white_color()` evaluates to 0, and
        // `black_color()` to 4 (and thus to a true Boolean).
        const boost::filtered_graph dual(cdt_2_dual_graph,
                                         [](typename CDT_2::Edge edge) {
                                           return false == edge.first->is_constrained(edge.second);
                                         });
        using Color_map = typename CDT_2_types::Color_map_is_outside_the_face;
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
    }  // end of the construction of the CDT_2
  }

  void search_for_missing_subfaces(CDT_3_face_index polygon_contraint_id)
  {
    const CDT_2& cdt_2 = face_cdt_2[polygon_contraint_id];

    for(const auto fh: cdt_2.all_face_handles())
    {
      if(fh->info().is_outside_the_face) continue;
      const auto v0 = fh->vertex(0)->info().vertex_handle_3d;
      const auto v1 = fh->vertex(1)->info().vertex_handle_3d;
      const auto v2 = fh->vertex(2)->info().vertex_handle_3d;
      Cell_handle c;
      int i, j, k;
      if(!tr.is_facet(v0, v1, v2, c, i, j, k)) {
        fh->info().missing_subface = true;
        face_constraint_misses_subfaces.set(static_cast<std::size_t>(polygon_contraint_id));
#if CGAL_DEBUG_CDT_3
        std::cerr << "Missing triangle: \n";
        write_triangle(std::cerr, v0, v1, v2);
#endif // CGAL_DEBUG_CDT_3
      } else {
        fh->info().missing_subface = false;
        const int facet_index = 6 - i - j - k;
        c->set_facet_constraint(facet_index, polygon_contraint_id, fh);
        if(tr.dimension() > 2) {
          const auto [n, n_index] = tr.mirror_facet({c, facet_index});
          n->set_facet_constraint(n_index, polygon_contraint_id, fh);
        }
      }
    }
  }

  static auto region(const CDT_2& cdt_2, CDT_2_face_handle fh)
  {
    std::vector<CDT_2_face_handle> fh_region;
    const auto cdt_2_dual_graph = dual(cdt_2.tds());
    const boost::filtered_graph dual(
        cdt_2_dual_graph,
        [](auto edge) {
          const auto face = edge.first;
          const auto i = unsigned(edge.second);
          return false == face->info().is_edge_also_in_3d_triangulation.test(i);
        },
        [](auto face_handle) { return false == face_handle->info().is_outside_the_face; });
    boost::breadth_first_search(dual, fh,
                                boost::color_map(typename CDT_2_types::Color_map_is_in_region())
                                    .visitor(boost::make_bfs_visitor(boost::write_property(
                                        boost::typed_identity_property_map<CDT_2_face_handle>(),
                                        std::back_inserter(fh_region), boost::on_finish_vertex()))));
    return fh_region;
  }

  auto brute_force_border_3_of_region(const std::vector<CDT_2_face_handle>& fh_region) {
    std::set<std::pair<Vertex_handle, Vertex_handle>> border_edges;
    for(auto fh: fh_region) {
      for(int i = 0; i < 3; ++i) {
        const auto va = fh->vertex(CDT_2::cw(i))->info().vertex_handle_3d;
        const auto vb = fh->vertex(CDT_2::ccw(i))->info().vertex_handle_3d;
        if(this->tds().is_edge(va, vb)) {
          border_edges.insert(CGAL::make_sorted_pair(va, vb));
        }
      }
    }
    return border_edges;
  }

  void restore_face(CDT_3_face_index i) {
    const CDT_2& cdt_2 = face_cdt_2[i];
#if CGAL_DEBUG_CDT_3
    std::cerr << "cdt_2 has " << cdt_2.number_of_vertices() << " vertices\n";
#endif // CGAL_DEBUG_CDT_3
    for(auto edge : cdt_2.finite_edges()) {
      const auto fh = edge.first;
      const auto i = edge.second;
      const auto va_3d = fh->vertex(cdt_2.cw(i))->info().vertex_handle_3d;
      const auto vb_3d = fh->vertex(cdt_2.ccw(i))->info().vertex_handle_3d;
      const bool is_3d = this->tds().is_edge(va_3d, vb_3d);
#if CGAL_DEBUG_CDT_3 && __has_include(<format>)
      std::cerr << std::format("Edge is 3D: {:6}  ({} , {})\n",
                                is_3d,
                                oformat(this->point(va_3d)),
                                oformat(this->point(vb_3d)));
#endif // CGAL_DEBUG_CDT_3
      CGAL_assertion(is_3d || !cdt_2.is_constrained(edge));
      fh->info().is_edge_also_in_3d_triangulation[unsigned(i)] = is_3d;
      const auto reverse_edge = cdt_2.mirror_edge(edge);
      reverse_edge.first->info().is_edge_also_in_3d_triangulation[unsigned(reverse_edge.second)] = is_3d;
    }
    for(CDT_2_face_handle fh : cdt_2.finite_face_handles()) {
      if(fh->info().is_outside_the_face) continue;
      if(false == fh->info().missing_subface) continue;
      auto fh_region = region(cdt_2, fh);
      assert(!fh_region.empty());
      assert(fh == fh_region[0]);
      auto border_edges = brute_force_border_3_of_region(fh_region);
#if CGAL_DEBUG_CDT_3
      std::cerr << "region size is: " << fh_region.size() << "\n";
      std::cerr << "region border size is: " << border_edges.size() << "\n";
      if(border_edges.size() < 3) {
        std::ofstream dump_region("dump_region_with_size_2.polylines.txt");
        dump_region.precision(17);
        write_region(dump_region, fh_region);
      }
#endif // CGAL_DEBUG_CDT_3
      const auto found_seg = [&]() -> std::optional<Edge> {
        for(auto fh_2d : fh_region) {
          CGAL_assertion(true == fh_2d->info().missing_subface);
          CGAL_assertion(false == fh_2d->info().is_outside_the_face);
          for(int index = 0; index < 3; ++index) {
            const auto va_3d = fh_2d->vertex(cdt_2.cw(index))->info().vertex_handle_3d;
            const auto vb_3d = fh_2d->vertex(cdt_2.ccw(index))->info().vertex_handle_3d;
            Cell_handle c;
            int i, j;
            const bool is_3d = this->tds().is_edge(va_3d, vb_3d, c, i, j);
            CGAL_assertion(fh_2d->info().is_edge_also_in_3d_triangulation[unsigned(index)] == is_3d);
            if(is_3d) {
              auto cell_circ = this->incident_cells(c, i, j), end = cell_circ;
              CGAL_assertion(cell_circ != nullptr);
              do {
                if(this->is_infinite(cell_circ)) {
                  continue;
                }
                const auto index_va = cell_circ->index(va_3d);
                const auto index_vb = cell_circ->index(vb_3d);
                const auto index_vc = this->next_around_edge(index_va, index_vb);
                const auto index_vd = this->next_around_edge(index_vb, index_va);
                const auto vc = cell_circ->vertex(index_vd);
                const auto vd = cell_circ->vertex(index_vc);
                CGAL_assertion(cell_circ->has_vertex(vc));
                CGAL_assertion(cell_circ->has_vertex(vd));
                const auto pc = this->point(vc);
                const auto pd = this->point(vd);
                const typename Geom_traits::Segment_3 seg{pc, pd};
                for(auto fh_2d : fh_region) {
                  const auto triangle = cdt_2.triangle(fh_2d);
                  if(do_intersect(seg, triangle)) {
                    std::cerr << "Segment " << seg << " intersects triangle " << triangle << "\n";
                    return { Edge{cell_circ, index_vc, index_vd} };
                  }
                }
              } while(++cell_circ != end);
            }
          }
        }
        return {};
      }();
      if(!found_seg) {
        std::cerr << "No segment found\n";
        {
          std::ofstream dump("dump.binary.cgal");
          CGAL::Mesh_3::save_binary_file(dump, *this);
          std::ofstream dump_region("dump_region.polylines.txt");
          dump_region.precision(17);
          write_region(dump_region, fh_region);
        }
      }
      CGAL_assertion(found_seg != std::nullopt);
    }
  }

public:
  void restore_constrained_Delaunay()
  {
    for(int i = 0, end = face_constraint_misses_subfaces.size(); i < end; ++i) {
      fill_cdt_2(i);
      search_for_missing_subfaces(i);
    }
    const auto npos = face_constraint_misses_subfaces.npos;
    auto i = face_constraint_misses_subfaces.find_first();
    while(i != npos) {
      restore_face(i);
      i = face_constraint_misses_subfaces.find_next(i);
    }
  }

  void write_region(std::ostream& out, auto region)
  {
    for(auto fh_2d : region) {
      write_2d_triangle(out, fh_2d);
    }
  }

  void write_triangle(std::ostream &out,
                      Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
  {
    out << "4"
        << " " << tr.point(v0) << " " << tr.point(v1) << " " << tr.point(v2)
        << " " << tr.point(v0) << '\n';
  }

  void write_2d_triangle(std::ostream &out, const CDT_2_face_handle fh)
  {
    const auto v0 = fh->vertex(0)->info().vertex_handle_3d;
    const auto v1 = fh->vertex(1)->info().vertex_handle_3d;
    const auto v2 = fh->vertex(2)->info().vertex_handle_3d;
    write_triangle(out, v0, v1, v2);
  }

  void write_missing_subfaces_file(std::ostream& out) {
    const auto npos = face_constraint_misses_subfaces.npos;
    auto i = face_constraint_misses_subfaces.find_first();
    while(i != npos) {
      const CDT_2& cdt = face_cdt_2[i];
      for(const auto fh: cdt.finite_face_handles()) {
        if (false == fh->info().is_outside_the_face &&
            true == fh->info().missing_subface)
        {
          write_2d_triangle(out, fh);
        }
      }
      i = face_constraint_misses_subfaces.find_next(i);
    }
  }

  /// @{
  /// remove functions cannot be called
  void remove(Vertex_handle) = delete;
  void remove_cluster() = delete;
  /// @}

protected:
  T_3 &tr = *this;
  Conforming_Dt &conforming_dt = *this;
  Insert_in_conflict_visitor insert_in_conflict_visitor = {*this};
  std::vector<CDT_2> face_cdt_2;
  struct Face_edge {
    Constraint_id constraint_id;
    bool is_reverse = false;
  };
  std::vector<std::vector<Face_edge>> face_border;
  boost::dynamic_bitset<> face_constraint_misses_subfaces;
};

} // end CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
