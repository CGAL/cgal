// Copyright (c) 2019-2023  GeometryFactory Sarl (France).
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

#include <CGAL/Triangulation_3/internal/CDT_3_config.h>

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_3.h>

#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_data_structure_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Compact_container.h>

#include <CGAL/Mesh_3/io_signature.h>

#include <CGAL/Conforming_Delaunay_triangulation_3.h>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/optional.hpp>
#include <boost/dynamic_bitset.hpp>

#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <unordered_map>
#include <ranges>
#if __has_include(<format>)
#  include <format>
#  include <concepts>
#elif CGAL_DEBUG_CDT_3
#  error "Compiler needs <format>"
#endif

namespace CGAL {

#if __cpp_lib_concepts >= 201806L
template <typename Polygon, typename Kernel>
concept Polygon_3 = std::ranges::common_range<Polygon>
      && (std::is_convertible_v<std::ranges::range_value_t<Polygon>,
                                typename Kernel::Point_3>);
#endif // concepts

using CDT_3_face_index = int;

template <typename Gt, typename Vb = Triangulation_vertex_base_3<Gt> >
class Constrained_Delaunay_triangulation_vertex_base_3
  : public Conforming_Delaunay_triangulation_vertex_base_3<Gt, Vb>
{
  using Base = Conforming_Delaunay_triangulation_vertex_base_3<Gt, Vb>;
public:
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
      if(IO::is_ascii(os)) {
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
      if(IO::is_ascii(is)) {
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

struct With_point_and_info_tag {};

template <class DSC, bool Const>
struct Output_rep<CGAL::internal::CC_iterator<DSC, Const>, With_point_and_info_tag>
  : public Output_rep<CGAL::internal::CC_iterator<DSC, Const>>
{
  using Base = Output_rep<CGAL::internal::CC_iterator<DSC, Const>>;
  using Base::Base;

  std::ostream& operator()(std::ostream& out) const {
    Base::operator()(out);
    if(this->it.operator->() != nullptr)
      return  out << (this->it->is_Steiner_vertex_on_edge() ? "(Steiner)" : "")
                  << (this->it->is_Steiner_vertex_in_face() ? "(Steiner in face)" : "")
                  << "= " << this->it->point();
    else
      return out;
  }
};


template <typename T_3>
class Constrained_Delaunay_triangulation_3 : public Conforming_Delaunay_triangulation_3<T_3> {
public:
  using Conforming_Dt = Conforming_Delaunay_triangulation_3<T_3>;
  using Vertex_handle = typename T_3::Vertex_handle;
  using Cell_handle = typename T_3::Cell_handle;
  using Edge = typename T_3::Edge;
  using Facet = typename T_3::Facet;
  using Point_3 = typename T_3::Point;
  using Segment_3 = typename T_3::Geom_traits::Segment_3;
  using Vector_3 = typename T_3::Geom_traits::Vector_3;
  using Locate_type = typename T_3::Locate_type;
  using Geom_traits = typename T_3::Geom_traits;

  using Face_index = CDT_3_face_index;

  static std::string io_signature() {
  return Get_io_signature<Conforming_Dt>()();
}

private:
  struct CDT_2_types {
    struct Projection_traits : public Projection_traits_3<Geom_traits> {
      // struct Side_of_oriented_circle_2 : public Geom_traits::Coplanar_side_of_bounded_circle_3 {
      //   using result_type = Oriented_side;
      //   template <typename... Arg> auto operator()(Arg&&... arg) const
      //   {
      //     return CGAL::enum_cast<Oriented_side>(
      //         Geom_traits::Coplanar_side_of_bounded_circle_3::operator()(std::forward<Arg>(arg)...));
      //   }
      // };
      // Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const { return {}; }

      // using Compare_xy_2 = typename Geom_traits::Compare_xyz_3;
      // Compare_xy_2 compare_xy_2_object() const { return {}; }

      using Projection_traits_3<Geom_traits>::Projection_traits_3; // inherit cstr
    };
    static_assert(std::is_nothrow_move_constructible<Projection_traits>::value,
                  "move cstr is missing");

    struct Vertex_info {
      Vertex_handle vertex_handle_3d = {};
    };

    using Color_value_type = std::int8_t;
    struct Face_info {
      Color_value_type is_outside_the_face = 0;
      Color_value_type is_in_region = 0;
      std::bitset<3> is_edge_also_in_3d_triangulation = 0;
      bool missing_subface = true;
    };
    using Vb1 = Triangulation_vertex_base_with_info_2<Vertex_info,
                                                      Projection_traits>;
    using Vb = Base_with_time_stamp<Vb1>;
    using Fb1 = Triangulation_face_base_with_info_2<Face_info,
                                                    Projection_traits>;
    using Fb = Constrained_triangulation_face_base_2<Projection_traits, Fb1>;
    using TDS = Triangulation_data_structure_2<Vb,Fb>;
    using Itag = No_constraint_intersection_requiring_constructions_tag;
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
  struct PLC_error { int face_index; int region_index; };
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
            std::cerr << "Add missing triangle (from visitor), face #F" << face_id << ": \n";
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

    void insert_Steiner_point_on_constraint(Constraint_id constraint,
                                            Vertex_handle va,
                                            Vertex_handle vb,
                                            Vertex_handle v_Steiner) const
    {
      const auto point = self.point(v_Steiner);
      if(!self.cdt_2_are_initialized) return;
      for(const auto [_, poly_id] : CGAL::make_range(self.constraint_to_faces.equal_range(constraint))) {
        auto& cdt_2 = self.face_cdt_2[poly_id];

        auto opt_edge = self.edge_of_cdt_2(cdt_2, va, vb);
        CGAL_assume(opt_edge != std::nullopt);
        CGAL_assertion(cdt_2.is_constrained(*opt_edge));
        auto [fh_2d, edge_index]= *opt_edge;
        const auto va_2d = fh_2d->vertex(cdt_2.cw(edge_index));
        const auto vb_2d = fh_2d->vertex(cdt_2.ccw(edge_index));

        const auto [mirror_fh_2d, mirror_edge_index] = cdt_2.mirror_edge(*opt_edge);

        const auto outside_on_right = fh_2d->info().is_outside_the_face;
        const auto outside_on_left = mirror_fh_2d->info().is_outside_the_face;
        const auto in_region_on_right = fh_2d->info().is_in_region;
        const auto in_region_on_left = mirror_fh_2d->info().is_in_region;

        fh_2d->set_constraint(edge_index, false);
        mirror_fh_2d->set_constraint(mirror_edge_index, false);
        const auto v_Steiner_2d = cdt_2.insert(point, fh_2d);
        v_Steiner_2d->info().vertex_handle_3d = v_Steiner;
        cdt_2.insert_constraint(va_2d, v_Steiner_2d);
        cdt_2.insert_constraint(vb_2d, v_Steiner_2d);

        // update the edge {fh_2d, edge_index}
        const bool is_edge = cdt_2.is_edge(va_2d, v_Steiner_2d, fh_2d, edge_index);
        CGAL_assume(is_edge);

        auto fc = cdt_2.incident_faces(v_Steiner_2d, fh_2d), fcbegin(fc);
        // circulators are counter-clockwise, so we start at the right of [va,v]
        do {
          fc->info().is_outside_the_face = outside_on_right;
          fc->info().is_in_region = in_region_on_right;
          ++fc;
        } while ( fc->vertex(cdt_2.ccw(fc->index(v_Steiner_2d))) != vb_2d );

        do {
          fc->info().is_outside_the_face = outside_on_left;
          fc->info().is_in_region = in_region_on_left;
        } while(++fc != fcbegin);
        fc = cdt_2.incident_faces(v_Steiner_2d, fh_2d);
        fcbegin  = fc;
        do {
          fc->info().missing_subface = true;
          const auto v_Steiner_index = fc->index(v_Steiner_2d);
          const auto other_edge = cdt_2.mirror_edge({fc, v_Steiner_index});
          fc->info().is_edge_also_in_3d_triangulation.set(other_edge.first->info().is_edge_also_in_3d_triangulation.test(other_edge.second));
        } while(++fc != fcbegin);

        self.face_constraint_misses_subfaces.set(poly_id);
      }
      conforming_dt_visitor.insert_Steiner_point_on_constraint(constraint, va, vb, v_Steiner);
    }
  };

public:
  Vertex_handle insert(const Point_3 &p, Locate_type lt, Cell_handle c,
                       int li, int lj, bool restore_Delaunay = true)
  {
    auto v = Conforming_Dt::insert_impl(p, lt, c, li, lj, insert_in_conflict_visitor);
    if(restore_Delaunay) {
      Conforming_Dt::restore_Delaunay(insert_in_conflict_visitor);
    }
    return v;
  }

  Vertex_handle insert(const Point_3 &p, Cell_handle start = {}, bool restore_Delaunay = true) {
    Locate_type lt;
    int li, lj;

    Cell_handle c = tr.locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj, restore_Delaunay);
  }

  Constraint_id insert_constrained_edge(Vertex_handle va, Vertex_handle vb, bool restore_Delaunay = true)
  {
    const auto id = this->insert_constrained_edge_impl(va, vb, insert_in_conflict_visitor);
    if(restore_Delaunay) {
      this->restore_Delaunay();
    }
    return id;
  }

  void restore_Delaunay() {
    Conforming_Dt::restore_Delaunay(insert_in_conflict_visitor);
  }

  bool is_constrained(Facet f) const {
    return f.first->is_facet_constrained(f.second);
  }

  void set_facet_constrained(Facet f, CDT_3_face_index polygon_contraint_id,
                             CDT_2_face_handle fh)
  {
    const auto [c, facet_index] = f;
    c->set_facet_constraint(facet_index, polygon_contraint_id, fh);
    if(tr.dimension() > 2) {
      const auto [n, n_index] = tr.mirror_facet({c, facet_index});
      n->set_facet_constraint(n_index, polygon_contraint_id, fh);
    }
    // for(int i = 0; i < 3; ++i) {
    //   const int index0 = tr.vertex_triple_index(facet_index, i);
    //   const int index1 = tr.vertex_triple_index(facet_index, tr.cw(i));
    //   const Edge edge{c, index0, index1};
    //   auto facet_circ = this->incident_facets(edge);
    //   const auto facet_circ_end = facet_circ;
    //   auto counter = 0u;
    //   do {
    //     if(is_constrained(*facet_circ)) {
    //       ++counter;
    //     }
    //   } while(++facet_circ != facet_circ_end);
    //   CGAL_assertion(counter <= 2);
    // }
  }

  template <CGAL_TYPE_CONSTRAINT(Polygon_3<Geom_traits>) Polygon>
  boost::optional<Face_index> register_new_constrained_polygon(Polygon&& polygon) {
    return insert_constrained_polygon(std::forward<Polygon>(polygon), false);
  }


template <CGAL_TYPE_CONSTRAINT(Polygon_3<Geom_traits>) Polygon>
  boost::optional<Face_index> insert_constrained_polygon(Polygon&& polygon, bool restore_Delaunay = true) {
    std::vector<Vertex_handle> handles;
    handles.reserve(polygon.size());
    boost::optional<Cell_handle> hint;
    for(const auto& p : polygon) {
      const auto v = this->insert(p, hint.value_or(Cell_handle{}), restore_Delaunay);
      handles.push_back(v);
      hint = v->cell();
    }
    return insert_constrained_face(std::move(handles), restore_Delaunay);
  }

  template <typename Vertex_handles>
  CGAL_CPP20_REQUIRE_CLAUSE(
      std::ranges::common_range<Vertex_handles>
      && (std::is_convertible_v<std::ranges::range_value_t<Vertex_handles>, Vertex_handle>))
  boost::optional<Face_index> insert_constrained_face(Vertex_handles&& vertex_handles, bool restore_Delaunay = true) {
#if CGAL_DEBUG_CDT_3 & 2
    std::cerr << "insert_constrained_face (" << std::size(vertex_handles) << " vertices):\n";
    for(auto v: vertex_handles) {
      std::cerr << "  - " << this->display_vert(v) << '\n';
    }
#endif // CGAL_DEBUG_CDT_3 & 2
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
    const auto polygon_contraint_id = static_cast<CDT_3_face_index>(this->face_border.size() - 1);
    do {
      const auto va = *circ;
      ++circ;
      const auto vb = *circ;
      const auto c_id = this->constraint_from_extremities(va, vb);
      if(c_id != Constraint_id{}) {
        const bool constraint_c_id_is_reversed = true;
        border.push_back(Face_edge{c_id, constraint_c_id_is_reversed});
        constraint_to_faces.emplace(c_id, polygon_contraint_id);
      } else {
        const auto c_id = this->insert_constrained_edge(va, vb, restore_Delaunay);
        CGAL_assertion(c_id != Constraint_id{});
        border.push_back(Face_edge{c_id});
        constraint_to_faces.emplace(c_id, polygon_contraint_id);
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
      if (accumulated_normal.x() < 0 ||
          (accumulated_normal.x() == 0 && accumulated_normal.y() < 0) ||
          (accumulated_normal.x() == 0 && accumulated_normal.y() == 0 &&
           accumulated_normal.z() < 0)
          )
      {
        accumulated_normal = - accumulated_normal;
      }
      return accumulated_normal;
    }();

    face_cdt_2.emplace_back(CDT_2_traits{accumulated_normal});
    face_constraint_misses_subfaces.resize(face_cdt_2.size());

    return polygon_contraint_id;
  }

  auto sequence_of_Steiner_vertices(Vertex_handle va, Vertex_handle vb) const
  {
    std::vector<Vertex_handle> steiner_vertices;
    const auto c_id = this->constraint_from_extremities(va, vb);
    if(c_id != Constraint_id{}) {
      auto vit = this->constraint_hierarchy.vertices_in_constraint_begin(c_id);
      auto v_end = this->constraint_hierarchy.vertices_in_constraint_end(c_id);
      CGAL_assertion_code(const auto constraint_size = std::distance(vit, v_end);) if(vit != v_end)
      {
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
            steiner_vertices.push_back(*vit);
            if(constraint_c_id_is_reversed) {
              --vit;
            } else {
              ++vit;
            };
          }
        }
      }
    } else {
      CGAL_error();
    }
    return steiner_vertices;
  }

private:
  void fill_cdt_2(CDT_2& cdt_2, CDT_3_face_index polygon_contraint_id)
  {
#if CGAL_DEBUG_CDT_3
    std::cerr << "Polygon #" << polygon_contraint_id << " normal is: " << cdt_2.geom_traits().normal() << '\n';
#endif // CGAL_DEBUG_CDT_3
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
#if CGAL_DEBUG_CDT_3
    std::cerr << "  points\n";
    do {
      std::cerr << "     " << tr.point(*circ) << '\n';
    } while(++circ != circ_end);
#endif

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
#if CGAL_DEBUG_CDT_3 & 4
                std::cerr << "cdt_2.insert_constraint ("
                          << tr.point(previous_2d->info().vertex_handle_3d)
                          << " , "
                          << tr.point(vh_2d->info().vertex_handle_3d)
                          << ")\n";
#endif // CGAL_DEBUG_CDT_3
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
#if CGAL_DEBUG_CDT_3 & 4
        std::cerr << "cdt_2.insert_constraint ("
                  << tr.point(previous_2d->info().vertex_handle_3d)
                  << " , "
                  << tr.point(vh_2d->info().vertex_handle_3d)
                  << ")\n";
#endif // CGAL_DEBUG_CDT_3
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
      for(const auto fh: cdt_2.finite_face_handles()) {
        if(!fh->info().is_outside_the_face) ++counter;
      }
      std::cerr << counter << " triangles(s) in the face\n";
#endif // CGAL_DEBUG_CDT_3
      if(Algebraic_structure_traits<typename Geom_traits::FT>::Is_exact::value &&
         !std::all_of(cdt_2.finite_face_handles().begin(), cdt_2.finite_face_handles().end(), [=](const auto fh) {
           const auto p0 = cdt_2.point(fh->vertex(0));
           const auto v1 = cdt_2.point(fh->vertex(1)) - p0;
           const auto v2 = cdt_2.point(fh->vertex(2)) - p0;
           return cross_product(cdt_2.geom_traits().normal(), cross_product(v1, v2)) == NULL_VECTOR;
         }))
      {
        std::cerr << std::string("Polygon #") + std::to_string(polygon_contraint_id) +
                         " is not coplanar.\n";
      }
    } // end of the construction of the CDT_2

    if(cdt_2.number_of_vertices() == 4) {
      // for polygons with 4 vertices, 2 triangles
      for(auto [ch, index]: cdt_2.finite_edges()) {
        if(!ch->is_constrained(index)) {
          // here the edge {ch, index} is the diagonal [ac] of the polygon [abcd]
          const auto vb = ch->vertex(index);
          const auto [ch2, index2] = cdt_2.mirror_edge({ch, index});
          const auto vd = ch2->vertex(index2);
          CGAL_assertion(!cdt_2.is_edge(vb, vd));
          const auto vb_3d = vb->info().vertex_handle_3d;
          const auto vd_3d = vd->info().vertex_handle_3d;
          if(tr.tds().is_edge(vb_3d, vd_3d)) {
            // let's insert the diagonal [bd] in the CDT_2
            cdt_2.insert_constraint(vb, vd);
#if CGAL_DEBUG_CDT_3 & 64
            std::cerr << "NOTE: CDT_2 has 4 vertices. Flip the diagonal\n";
#endif
          }
        }
        break;
      }
    }
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
        std::cerr << std::format("Missing triangle in polygon #{}:\n", polygon_contraint_id);
        write_triangle(std::cerr, v0, v1, v2);
#endif // CGAL_DEBUG_CDT_3
      } else {
        fh->info().missing_subface = false;
        const int facet_index = 6 - i - j - k;
        set_facet_constrained({c, facet_index}, polygon_contraint_id, fh);
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
    CGAL_assertion(!fh_region.empty());
    CGAL_assertion(fh == fh_region[0]);
    return fh_region;
  }

  auto brute_force_border_3_of_region(const std::vector<CDT_2_face_handle>& fh_region) {
    std::set<std::pair<Vertex_handle, Vertex_handle>> border_edges_set;
    for(const auto fh: fh_region) {
      for(int i = 0; i < 3; ++i) {
        const auto va = fh->vertex(CDT_2::cw(i))->info().vertex_handle_3d;
        const auto vb = fh->vertex(CDT_2::ccw(i))->info().vertex_handle_3d;
        if(this->tds().is_edge(va, vb)) {
          CGAL_assertion_code(const auto result =)
          border_edges_set.insert(CGAL::make_sorted_pair(va, vb));
          CGAL_assertion(result.second == true);
        }
      }
    }
    std::vector<Edge> border_edges;
    border_edges.reserve(border_edges_set.size());
    for(const auto& [va, vb]: border_edges_set) {
      Cell_handle c;
      int i, j;
      CGAL_assume(this->tds().is_edge(va, vb, c, i, j));
      border_edges.emplace_back(c, i, j);
    }
#if CGAL_DEBUG_CDT_3
    std::cerr << "region size is: " << fh_region.size() << "\n";
    std::cerr << "region border size is: " << border_edges.size() << "\n";
#endif // CGAL_DEBUG_CDT_3
    return border_edges;
  }

  struct Intersection_error : public std::runtime_error {
    using Seg = typename Geom_traits::Segment_3;
    using Tri = typename Geom_traits::Triangle_3;
    Intersection_error(Seg s, Tri t, std::string what) : std::runtime_error(what), segment(s), triangle(t) {}

    Seg segment;
    Tri triangle;
  };

  int does_edge_intersect_region(Cell_handle cell, int index_vc, int index_vd,
                                 const CDT_2& cdt_2, const auto& fh_region)
  {
    const auto vc = cell->vertex(index_vd);
    const auto vd = cell->vertex(index_vc);
    const auto pc = this->point(vc);
    const auto pd = this->point(vd);
    const typename Geom_traits::Segment_3 seg{pc, pd};
    for(const auto fh_2d : fh_region) {
      const auto t0 = cdt_2.point(fh_2d->vertex(0));
      const auto t1 = cdt_2.point(fh_2d->vertex(1));
      const auto t2 = cdt_2.point(fh_2d->vertex(2));

      const auto opc = CGAL::orientation(t0, t1, t2, pc);
      const auto opd = CGAL::orientation(t0, t1, t2, pd);
      if(opc == CGAL::COPLANAR || opd == CGAL::COPLANAR || opc == opd) {
        continue;
      } else {
        // otherwise the segment intersects the plane of the triangle
        if(CGAL::orientation(pc, pd, t0, t1) != opc &&
           CGAL::orientation(pc, pd, t1, t2) != opc &&
           CGAL::orientation(pc, pd, t2, t0) != opc)
        {
          return static_cast<int>(opc);
        }
      }
    }
    return 0;
  }

  // Given a region and a border edge of it, returns an edge in the link of the
  // border edge that intersects the region.
  // The returned edge has its first vertex above the region.
  std::optional<Edge> search_first_intersection(CDT_3_face_index /*face_index*/,
                                                const CDT_2& cdt_2,
                                                const auto& fh_region,
                                                const Edge border_edge)
  {
    const auto [c, i, j] = border_edge;
    const Vertex_handle va_3d = c->vertex(i);
    const Vertex_handle vb_3d = c->vertex(j);
    auto cell_circ = this->incident_cells(c, i, j), end = cell_circ;
#if CGAL_DEBUG_CDT_3 > 32
    std::ofstream dump_edges_around("dump_edges_around.polylines.txt");
    dump_edges_around.precision(17);
#endif // CGAL_DEBUG_CDT_3
    CGAL_assertion(cell_circ != nullptr);
    do {
      if(this->is_infinite(cell_circ)) {
        continue;
      }
      const auto index_va = cell_circ->index(va_3d);
      const auto index_vb = cell_circ->index(vb_3d);
      const auto index_vc = this->next_around_edge(index_va, index_vb);
      const auto index_vd = this->next_around_edge(index_vb, index_va);
#if CGAL_DEBUG_CDT_3 > 32
      write_segment(dump_edges_around, cell_circ->vertex(index_vc), cell_circ->vertex(index_vd));
#endif
      int cd_intersects_region = does_edge_intersect_region(cell_circ, index_vc, index_vd, cdt_2, fh_region);
      if(cd_intersects_region == 1) {
        return { Edge{cell_circ, index_vc, index_vd} };
      }
      if(cd_intersects_region == -1) {
        return { Edge{cell_circ, index_vd, index_vc} };
      }
    } while(++cell_circ != end);
    return {};
  }

  struct Next_region : std::logic_error {
    using std::logic_error::logic_error;
    CDT_2_face_handle fh_2d;
    // create a new region
    Next_region(const std::string& what, CDT_2_face_handle fh) : std::logic_error(what), fh_2d(fh) {}
  };

  static constexpr auto vertex_pair(Edge e) {
    const auto [c, i, j] = e;
    return std::pair<Vertex_handle, Vertex_handle>{c->vertex(i), c->vertex(j)};
  }

  auto construct_cavities(CDT_3_face_index face_index,
                          int region_count,
                          const CDT_2& cdt_2,
                          const auto& fh_region,
                          const auto& polygon_vertices,
                          Edge first_intersecting_edge)
  {
    // outputs
    struct Outputs
    {
      std::vector<Edge> intersecting_edges;
      std::set<Cell_handle> intersecting_cells;
      std::set<Vertex_handle> vertices_of_upper_cavity;
      std::set<Vertex_handle> vertices_of_lower_cavity;
      std::set<Facet> facets_of_upper_cavity;
      std::set<Facet> facets_of_lower_cavity;
    } outputs{
        {}, {}, {polygon_vertices.begin(), polygon_vertices.end()}, {polygon_vertices.begin(), polygon_vertices.end()},
        {}, {}};

    auto& [_, intersecting_cells, vertices_of_upper_cavity, vertices_of_lower_cavity,
           facets_of_upper_cavity, facets_of_lower_cavity] = outputs;

    auto& intersecting_edges = outputs.intersecting_edges;

    // marker for already visited elements
    std::set<Vertex_handle> visited_vertices;
    std::map<std::pair<Vertex_handle, Vertex_handle>, bool> visited_edges;
    std::set<Cell_handle> visited_cells;

    auto make__new_element_functor = [](auto& visited_set) {
      return [&visited_set](auto... e) {
        const auto [_, not_already_visited] = visited_set.emplace(e...);
        return not_already_visited;
      };
    };

    auto new_vertex = make__new_element_functor(visited_vertices);
    auto new_cell = make__new_element_functor(visited_cells);
    auto new_edge = [&](Vertex_handle v0, Vertex_handle v1, bool does_insersect) {
      return visited_edges.emplace(CGAL::make_sorted_pair(v0, v1), does_insersect);
    };

    intersecting_edges.push_back(first_intersecting_edge);
    const auto [v0, v1] = vertex_pair(first_intersecting_edge);
    (void)new_edge(v0, v1, true);
    for(std::size_t i = 0; i < intersecting_edges.size(); ++i) {
      const auto intersecting_edge = intersecting_edges[i];
      const auto [v_above, v_below] = vertex_pair(intersecting_edge);
#if CGAL_DEBUG_CDT_3
      std::cerr << std::format("restore_subface_region face index: {}, region #{}, Edge #{:6}: ({} , {})\n",
                                face_index, region_count, i,
                                IO::oformat(v_above, with_point),
                                IO::oformat(v_below, with_point));
#endif
      CGAL_assertion(false == polygon_vertices.contains(v_above));
      CGAL_assertion(false == polygon_vertices.contains(v_below));
      if(new_vertex(v_above)) {
        vertices_of_upper_cavity.insert(v_above);
      }
      if(new_vertex(v_below)) {
        vertices_of_lower_cavity.insert(v_below);
      }
      auto facet_circ = this->incident_facets(intersecting_edge);
      const auto facet_circ_end = facet_circ;
      do { // loop facets around [v_above, v_below]
        CGAL_assertion(false == this->is_infinite(*facet_circ));
        const auto cell = facet_circ->first;
        const auto facet_index = facet_circ->second;
        CGAL_assertion_msg(!cell->is_facet_constrained(facet_index), "intersecting polygons!");
        if(new_cell(cell)) {
          intersecting_cells.insert(cell);
        }
        const auto index_v_above = cell->index(v_above);
        const auto index_v_below = cell->index(v_below);
        const auto index_vc = 6 - index_v_above - index_v_below - facet_index;
        const auto vc = cell->vertex(index_vc);
        if(polygon_vertices.contains(vc)) continue; // intersecting edges cannot touch the border

        auto test_edge = [&](Vertex_handle v0, int index_v0, Vertex_handle v1, int index_v1, int expected) {
          auto [cached_value_it, not_visited] = new_edge(v0, v1, false);
          if(!not_visited) return cached_value_it->second;
          int v0v1_intersects_region = does_edge_intersect_region(cell, index_v0, index_v1, cdt_2, fh_region);
          if(v0v1_intersects_region != 0) {
            if(v0v1_intersects_region != expected) {
              throw PLC_error{face_index, region_count};
            }
            // report the edge with first vertex above the region
            if(v0v1_intersects_region < 0) {
              std::swap(index_v0, index_v1);
            }
            intersecting_edges.emplace_back(cell, index_v0, index_v1);
            cached_value_it->second = true;
            return true;
          }
          else {
            cached_value_it->second = false;
            return false;
          }
        };

        if(!test_edge(v_above, index_v_above, vc, index_vc, 1) &&
           !test_edge(v_below, index_v_below, vc, index_vc, -1))
        {
          dump_triangulation();
          dump_region(face_index, region_count, cdt_2);
          {
            std::ofstream out(std::string("dump_two_edges_") + std::to_string(face_index) + ".polylines.txt");
            write_segment(out, Edge{cell, index_v_above, index_vc});
            write_segment(out, Edge{cell, index_v_below, index_vc});
          }
          throw PLC_error{face_index, region_count};
        }
      } while(++facet_circ != facet_circ_end);
    }
    for(auto intersecting_edge: intersecting_edges) {
      const auto [v_above, v_below] = vertex_pair(intersecting_edge);

      auto cell_circ = this->incident_cells(intersecting_edge), end = cell_circ;
      CGAL_assume(cell_circ != nullptr);
      do {
        const Cell_handle cell = cell_circ;
        const auto index_v_above = cell->index(v_above);
        const auto index_v_below = cell->index(v_below);
        const auto cell_above = cell->neighbor(index_v_below);
        const auto cell_below = cell->neighbor(index_v_above);
        if(!intersecting_cells.contains(cell_above)) {
          facets_of_upper_cavity.emplace(cell_above, cell_above->index(cell));
        }
        if(!intersecting_cells.contains(cell_below)) {
          facets_of_lower_cavity.emplace(cell_below, cell_below->index(cell));
        }
      } while(++cell_circ != end);
    }
    return outputs;
  }

  template <typename Function>
  static void visit_convex_hull_of_triangulation(const Constrained_Delaunay_triangulation_3& tr, Function f)
  {
    const auto inf_vh = tr.infinite_vertex();
    tr.incident_cells(inf_vh, boost::make_function_output_iterator([&](Cell_handle c) {
                        const auto facet_index = c->index(inf_vh);
                        f(Facet{c, facet_index});
                        return true;
                      }));
  }

  static constexpr With_point_tag with_point{};
  static constexpr With_point_and_info_tag with_point_and_info{};

  void restore_subface_region(CDT_3_face_index face_index, int region_count,
                              const CDT_2& cdt_2, const auto& fh_region)
  {
    const auto border_edges = brute_force_border_3_of_region(fh_region);
    const auto polygon_vertices = [&]() {
      std::set<Vertex_handle> vertices;
      for(const auto& [c, i, j]: border_edges) {
        vertices.insert(c->vertex(i));
        vertices.insert(c->vertex(j));
      }
      return vertices;
    }();
#if CGAL_DEBUG_CDT_3
    std::cerr << "polygon_vertices.size() = " << polygon_vertices.size() << "\n";
    for(auto v : polygon_vertices) {
      std::cerr << std::format("  {}\n", IO::oformat(v, with_point));
    }
#endif
    const Edge first_border_edge{border_edges[0]};
    const auto found_edge_opt = search_first_intersection(face_index, cdt_2, fh_region, first_border_edge);

    [[maybe_unused]] auto debug_region_size_4 = [&] {
      if(polygon_vertices.size() == 4) {
        std::set<Vertex_handle> vertices;
        std::set<Vertex_handle> diagonal;
        for(auto fh : fh_region) {
          for(int i = 0; i < 3; ++i) {
            auto [it, new_vertex] = vertices.insert(fh->vertex(i)->info().vertex_handle_3d);
            if(!new_vertex) {
              diagonal.insert(*it);
            }
          }
        }
        std::set<Vertex_handle> other_diagonal;
        std::set_difference(polygon_vertices.begin(), polygon_vertices.end(),
                            diagonal.begin(), diagonal.end(),
                            std::inserter(other_diagonal, other_diagonal.begin()));
        CGAL_assertion(diagonal.size() == 2);
        CGAL_assertion(other_diagonal.size() == 2);
        std::cerr << std::format
            ("NOTE: diagonal: {:.6} {:.6}  {} in tr\n",
            IO::oformat(*diagonal.begin(), with_point),
            IO::oformat(*std::next(diagonal.begin()), with_point),
            tr.tds().is_edge(*diagonal.begin(), *std::next(diagonal.begin())) ? "IS" : "is NOT");
        std::cerr << std::format(
            "NOTE: the other diagonal: {:.6} {:.6}  {} in tr\n",
            IO::oformat(*other_diagonal.begin(), with_point),
            IO::oformat(*std::next(other_diagonal.begin()), with_point),
            tr.tds().is_edge(*other_diagonal.begin(), *std::next(other_diagonal.begin())) ? "IS" : "is NOT");
        if(cdt_2.geom_traits().side_of_oriented_circle_2_object()(
               (*polygon_vertices.begin())->point(), (*std::next(polygon_vertices.begin()))->point(),
               (*std::next(polygon_vertices.begin(), 2))->point(),
               (*std::next(polygon_vertices.begin(), 3))->point()) == CGAL::ZERO)
        {
          std::cerr << std::format(
              "NOTE: In polygon #{}, region {}, the 4 vertices are co-circular in the 2D triangulation\n",
              face_index, region_count);
        }
        if(CGAL::coplanar(
               (*polygon_vertices.begin())->point(),
               (*std::next(polygon_vertices.begin()))->point(),
               (*std::next(polygon_vertices.begin(), 2))->point(),
               (*std::next(polygon_vertices.begin(), 3))->point()))
        {
          std::cerr << std::format("NOTE: In polygon #{}, region {}, the 4 vertices are coplanar\n",
                                   face_index, region_count);
          if(CGAL::coplanar_side_of_bounded_circle(
               (*polygon_vertices.begin())->point(),
               (*std::next(polygon_vertices.begin()))->point(),
               (*std::next(polygon_vertices.begin(), 2))->point(),
               (*std::next(polygon_vertices.begin(), 3))->point()) == CGAL::ON_BOUNDARY)
          {
            std::cerr << std::format(
                "NOTE: In polygon #{}, region {}, the 4 vertices are co-circular in the 3D triangulation\n",
                face_index, region_count);
          }
        }
      }
    };
    if(!found_edge_opt) {
      // debug_region_size_4();
      // {
      //   Constrained_Delaunay_triangulation_3 new_tr;
      //   for(const auto v : polygon_vertices) {
      //     new_tr.insert(v->point());
      //   }
      //   std::cerr << "new_tr.dimension() = " << new_tr.dimension() << '\n';
      //   std::ofstream out(std::string("dump_polygon_") + std::to_string(face_index) + "_tr.off");
      //   out.precision(17);
      //   if(new_tr.dimension() == 2) {
      //     write_facets(out, new_tr, new_tr.finite_facets());
      //   }
      //   else {
      //     write_facets(out, new_tr, std::views::filter(new_tr.finite_facets(), [&](auto f) {
      //                    return new_tr.is_infinite(f.first) || new_tr.is_infinite(f.first->neighbor(f.second));
      //                  }));
      //   }
      // }
      // {
      //   dump_edge_link(std::string("dump_around_edge_") + std::to_string(face_index) + "_" +
      //                  std::to_string(region_count) + ".polylines.txt", first_border_edge);
      //   std::ofstream dump(std::string("dump_no_segment_found_") + std::to_string(face_index) + "_" +
      //                      std::to_string(region_count) + ".binary.cgal");
      //   CGAL::IO::save_binary_file(dump, *this);
      //   dump_region(face_index, region_count, cdt_2);
      // }
      throw Next_region{"No segment found", fh_region[0]};
    }
    CGAL_assertion(found_edge_opt != std::nullopt);

    const auto first_intersecting_edge = *found_edge_opt;
    auto cavities =
        construct_cavities(face_index, region_count, cdt_2, fh_region, polygon_vertices, first_intersecting_edge);
    auto& [intersecting_edges, intersecting_cells_vector, vertices_of_upper_cavity,
           vertices_of_lower_cavity, facets_of_upper_cavity, facets_of_lower_cavity] = cavities;

    const std::set<Cell_handle> intersecting_cells{intersecting_cells_vector.begin(), intersecting_cells_vector.end()};
    const std::set<Point_3> polygon_points = [&](){
      std::set<Point_3> polygon_points;
      for(auto vh : polygon_vertices) {
        polygon_points.insert(this->point(vh));
      }
      return polygon_points;
    }();

    auto is_facet_of_polygon = [&](const auto& tr, Facet f) {
      const auto [c, facet_index] = f;
      for(int i = 0; i < 3; ++i) {
        const auto vh = c->vertex(T_3::vertex_triple_index(facet_index, i));
        if(!polygon_points.contains(tr.point(vh))) {
          return false;
        }
      }
      return true;
    };

#if CGAL_DEBUG_CDT_3 & 64
    std::cerr << std::format("Cavity has {} cells and {} edges, "
                             "{} vertices in upper cavity and {} in lower, "
                             "{} facets in upper cavity and {} in lower\n",
                             intersecting_cells.size(),
                             intersecting_edges.size(),
                             vertices_of_upper_cavity.size(),
                             vertices_of_lower_cavity.size(),
                             facets_of_upper_cavity.size(),
                             facets_of_lower_cavity.size());
    if(intersecting_cells.size() > 3 || intersecting_edges.size() > 1) {
      std::cerr << "!! Interesting case !!\n";
      // dump_region(face_index, region_count, cdt_2);
      // {
      //   std::ofstream out(std::string("dump_intersecting_edges_") + std::to_string(face_index) + "_" +
      //                     std::to_string(region_count) + ".polylines.txt");
      //   out.precision(17);
      //   for(auto edge: intersecting_edges) {
      //     write_segment(out, edge);
      //   }
      // }
      // dump_facets_of_cavity(face_index, region_count, "lower", facets_of_lower_cavity);
      // dump_facets_of_cavity(face_index, region_count, "upper", facets_of_upper_cavity);
    }
#endif // CGAL_DEBUG_CDT_3
    typename T_3::Vertex_handle_unique_hash_map map_cavity_vertices_to_ambient_vertices;
    typename T_3::Vertex_handle_unique_hash_map map_lower_cavity_vertices_to_ambient_vertices;

#if CGAL_DEBUG_CDT_3 & 64
    std::cerr << "# upper cavity\n";
#endif // CGAL_DEBUG_CDT_3
    const auto [upper_cavity_triangulation, nb_of_add_vertices_upper] =
        triangulate_cavity(face_index, cdt_2, intersecting_cells, facets_of_upper_cavity,
                           map_cavity_vertices_to_ambient_vertices,
                           vertices_of_upper_cavity);
#if CGAL_DEBUG_CDT_3 & 64
    std::cerr << "# lower cavity\n";
#endif // CGAL_DEBUG_CDT_3
    const auto [lower_cavity_triangulation, nb_of_add_vertices_lower] =
        triangulate_cavity(face_index, cdt_2, intersecting_cells, facets_of_lower_cavity,
                           map_lower_cavity_vertices_to_ambient_vertices,
                           vertices_of_lower_cavity);

    CGAL_assertion(std::all_of(fh_region.begin(), fh_region.end(), [&](auto fh) {
      auto is_fh_facet_of = [this, fh](const auto& tr) -> std::optional<Facet> {
        const auto v0 = fh->vertex(0)->info().vertex_handle_3d;
        const auto v1 = fh->vertex(1)->info().vertex_handle_3d;
        const auto v2 = fh->vertex(2)->info().vertex_handle_3d;
        return this->vertex_triple_is_facet_of_other_triangulation(*this, v0, v1, v2, tr);
      };

      const bool fail_upper = !is_fh_facet_of(upper_cavity_triangulation);
      const bool fail_lower = !is_fh_facet_of(lower_cavity_triangulation);
      const bool test = !fail_upper && !fail_lower;
      if(fail_upper || fail_lower) {
        if(fail_upper) {
          std::cerr << "ERROR: Face is not a facet of the upper cavity\n";
        }
        if(fail_lower) {
          std::cerr << "ERROR: Face is not a facet of the lower cavity\n";
        }
        // debug_region_size_4();
        // dump_region(face_index, region_count, cdt_2);
        // dump_3d_triangulation(face_index, region_count, "lower", lower_cavity_triangulation);
        // dump_3d_triangulation(face_index, region_count, "upper", upper_cavity_triangulation);
        // auto dump_facets_of_cavity_border = [&](CDT_3_face_index face_index, int region_count, std::string type,
        //                                         const auto& cavity_triangulation) {
        //   std::ofstream out(std::string("dump_plane_facets_of_region_") + std::to_string(face_index) + "_" +
        //                     std::to_string(region_count) + "_" + type + ".off");
        //   std::ofstream other_out(std::string("dump_non_plane_facets_of_region_") + std::to_string(face_index) + "_" +
        //                     std::to_string(region_count) + "_" + type + ".off");
        //   out.precision(17);
        //   other_out.precision(17);

        //   std::vector<Facet> border_faces;
        //   std::vector<Facet> non_border_faces;
        //   visit_convex_hull_of_triangulation(cavity_triangulation,
        //       [&](Facet f) {
        //         if(is_facet_of_polygon(cavity_triangulation, f))
        //           border_faces.push_back(f);
        //         else
        //           non_border_faces.push_back(f);
        //       });
        //   CGAL_warning(!border_faces.empty());
        //   write_facets(out, cavity_triangulation, border_faces);
        //   write_facets(other_out, cavity_triangulation, non_border_faces);
        // };
        // dump_facets_of_cavity_border(face_index, region_count, "lower", lower_cavity_triangulation);
        // dump_facets_of_cavity_border(face_index, region_count, "upper", upper_cavity_triangulation);
        throw Next_region{"missing facet in polygon", fh_region[0]};
      }
      return test;
    }));

#if CGAL_DEBUG_CDT_3 & 64
    std::cerr << "# glu the upper triangulation of the cavity\n";

    if(nb_of_add_vertices_upper > 0 || nb_of_add_vertices_lower > 0)
    {
      std::cerr << std::format("!! Cavity has grown and not has {} cells and {} edges, "
                               "{} vertices in upper cavity and {} in lower, "
                               "{} facets in upper cavity and {} in lower\n",
                               intersecting_cells.size(),
                               intersecting_edges.size(),
                               vertices_of_upper_cavity.size(),
                               vertices_of_lower_cavity.size(),
                               facets_of_upper_cavity.size(),
                               facets_of_lower_cavity.size());
    }
#endif // CGAL_DEBUG_CDT_3

    typename T_3::Vertex_triple_Facet_map outer_map;
    auto add_to_outer_map = [&outer_map](typename T_3::Vertex_triple vt, Facet f,
                                         [[maybe_unused]] std::string_view extra = {}) {
      outer_map[vt] = f;
#if CGAL_DEBUG_CDT_3 & 64
      CGAL_assertion(vt.first != vt.second);
      CGAL_assertion(vt.first != vt.third);
      CGAL_assertion(vt.second != vt.third);
      std::cerr << std::format("outer map: Adding {}triple ({:.6}, {:.6}, {:.6})\n", extra,
                               IO::oformat(vt.first, with_point),
                               IO::oformat(vt.second, with_point),
                               IO::oformat(vt.third, with_point));
#endif // CGAL_DEBUG_CDT_3
    };
    auto fill_outer_map_of_cavity = [&](const auto&, const auto& facets) {
      for(auto f : facets) {
        typename T_3::Vertex_triple vt = this->make_vertex_triple(f);
        this->make_canonical_oriented_triple(vt);
        add_to_outer_map(vt, f);
      }
    };

    fill_outer_map_of_cavity(upper_cavity_triangulation, facets_of_upper_cavity);

    auto add_pseudo_cells_to_outer_map = [&](const auto& tr, bool is_upper_cavity) {
      std::vector<std::pair<Cell_handle, CDT_2_face_handle>> pseudo_cells;
      for(auto f : tr.finite_facets()) {
        if(!is_facet_of_polygon(tr, f))
          continue;
        const auto is_facet = facet_is_facet_of_cdt_2(tr, f, cdt_2);
        if(!is_facet) continue; // we might be in a sliver in the plane of the polygon
        const auto [fh_2d, reverse_orientation] = *is_facet;

        const auto vt_aux = this->make_vertex_triple(f);
        typename T_3::Vertex_triple vt(map_cavity_vertices_to_ambient_vertices[vt_aux.first],
                                       map_cavity_vertices_to_ambient_vertices[vt_aux.second],
                                       map_cavity_vertices_to_ambient_vertices[vt_aux.third]);
        this->make_canonical_oriented_triple(vt);
        if(reverse_orientation == is_upper_cavity) {
          std::swap(vt.second, vt.third);
        }
        auto new_cell = this->tds().create_cell();
        pseudo_cells.emplace_back(new_cell, fh_2d);
        new_cell->set_vertices(vt.first, vt.second, vt.third, this->infinite_vertex());
        CGAL_assertion(static_cast<bool>(facet_is_facet_of_cdt_2(*this, {new_cell, 3}, cdt_2)));
        add_to_outer_map(vt, {new_cell, 3}, "extra ");
      }
      return pseudo_cells;
    };
    const auto pseudo_cells = add_pseudo_cells_to_outer_map(upper_cavity_triangulation, true);

    auto inner_map_of_cavity = [&](const auto& tr) {
      typename T_3::Vertex_triple_Facet_map inner_map;
      auto add_facet_to_inner_map = [&](Facet f) {
        const auto vt_aux = this->make_vertex_triple(f);
        typename T_3::Vertex_triple vt(map_cavity_vertices_to_ambient_vertices[vt_aux.first],
                                       map_cavity_vertices_to_ambient_vertices[vt_aux.third],
                                       map_cavity_vertices_to_ambient_vertices[vt_aux.second]);
        this->make_canonical_oriented_triple(vt);
#if CGAL_DEBUG_CDT_3 & 64
        CGAL_assertion(vt.first != vt.second);
        CGAL_assertion(vt.first != vt.third);
        CGAL_assertion(vt.second != vt.third);
        std::cerr << std::format("inner map: Adding triple ({:.6}, {:.6}, {:.6})\n",
                                 IO::oformat(vt.first, with_point),
                                 IO::oformat(vt.second, with_point),
                                 IO::oformat(vt.third, with_point));
#endif // CGAL_DEBUG_CDT_3
        inner_map[vt] = f;
      };
      for(auto f : tr.finite_facets()) {
        add_facet_to_inner_map(f);
        add_facet_to_inner_map(this->mirror_facet(f));
      }
      return inner_map;
    };

    {
// #if CGAL_DEBUG_CDT_3 & 64
//       std::ofstream out("dump_upper_outer_map.off");
//       out.precision(17);
//       write_facets(out, *this, std::ranges::views::values(outer_map));
//       out.close();
// #endif // CGAL_DEBUG_CDT_3
      const auto upper_inner_map = inner_map_of_cavity(upper_cavity_triangulation);

      this->copy_triangulation_into_hole(map_cavity_vertices_to_ambient_vertices,
                                         std::move(outer_map),
                                         upper_inner_map,
                                         Emptyset_iterator{});
    }
#if CGAL_DEBUG_CDT_3 & 64
    std::cerr << "# glu the lower triangulation of the cavity\n";
#endif // CGAL_DEBUG_CDT_3

    map_cavity_vertices_to_ambient_vertices.clear();
    map_cavity_vertices_to_ambient_vertices = std::move(map_lower_cavity_vertices_to_ambient_vertices);

    outer_map.clear();
    std::vector<std::pair<Facet, CDT_2_face_handle>> new_constrained_facets;
    new_constrained_facets.reserve(pseudo_cells.size());
    for(const auto& [c, fh_2d] : pseudo_cells) {
      const Facet f = this->mirror_facet({c, 3});
      new_constrained_facets.emplace_back(f, fh_2d);
      CGAL_assertion(static_cast<bool>(facet_is_facet_of_cdt_2(*this, f, cdt_2)));
      auto vt = this->make_vertex_triple(f);
      this->make_canonical_oriented_triple(vt);
      add_to_outer_map(vt, f);
      this->tds().delete_cell(c);
    }
    fill_outer_map_of_cavity(lower_cavity_triangulation, facets_of_lower_cavity);
    {
      const auto lower_inner_map = inner_map_of_cavity(lower_cavity_triangulation);
#if CGAL_DEBUG_CDT_3 & 64
      std::cerr << "outer_map:\n";
      for(auto [vt, _] : outer_map) {
        std::cerr << std::format("  {:.6}, {:.6}, {:.6})\n",
                                 IO::oformat(vt.first,  with_point),
                                 IO::oformat(vt.second, with_point),
                                 IO::oformat(vt.third,  with_point));
      }
      std::ofstream out("dump_lower_outer_map.off");
      out.precision(17);
      write_facets(out, *this, std::ranges::views::values(outer_map));
      out.close();
#endif // CGAL_DEBUG_CDT_3
      this->copy_triangulation_into_hole(map_cavity_vertices_to_ambient_vertices, std::move(outer_map), lower_inner_map,
                                         Emptyset_iterator{});
    }
    for(auto c : intersecting_cells) {
      this->tds().delete_cell(c);
    }

    auto restore_markers = [&](Facet outside_facet) {
      const auto [outside_cell, outside_face_index] = outside_facet;
      const auto mirror_facet = this->mirror_facet(outside_facet);
      if(outside_cell->is_facet_constrained(outside_face_index)) {
        const auto poly_id = outside_cell->face_constraint_index(outside_face_index);
        const auto f2d = outside_cell->face_2(cdt_2, outside_face_index);
        set_facet_constrained(mirror_facet, poly_id, f2d);
      }
    };

    std::for_each(facets_of_lower_cavity.begin(), facets_of_lower_cavity.end(), restore_markers);
    std::for_each(facets_of_upper_cavity.begin(), facets_of_upper_cavity.end(), restore_markers);

    for(const auto& [f, f2d] : new_constrained_facets) {
      set_facet_constrained(f, face_index, f2d);
      f2d->info().missing_subface = false;
    }
    //CGAL_assertion(this->T_3::Tr_Base::is_valid(true));
  };

  struct Oriented_face_of_cdt_2 {
    CDT_2_face_handle fh;
    bool reversed_orientation = false;
  };

  template <typename Tr>
  static auto facet_is_facet_of_cdt_2(const Tr& tr, typename Tr::Facet f, const CDT_2& cdt_2)
      -> std::optional<Oriented_face_of_cdt_2>
  {
    const auto [c, facet_index] = f;
    const auto v0 = c->vertex(Tr::vertex_triple_index(facet_index, 0));
    const auto v1 = c->vertex(Tr::vertex_triple_index(facet_index, 1));
    const auto v2 = c->vertex(Tr::vertex_triple_index(facet_index, 2));

    auto v = [&, hint = CDT_2_face_handle{}](const auto& p) mutable {
      int i;
      typename CDT_2::Locate_type lt;
      const auto fh = cdt_2.locate(p, lt, i, hint);
      CGAL_assume(lt == CDT_2::VERTEX);
      hint = fh;
      return fh->vertex(i);
    };

    const auto cdt_2_v0 = v(tr.point(v0));
    const auto cdt_2_v1 = v(tr.point(v1));
    const auto cdt_2_v2 = v(tr.point(v2));

    CDT_2_face_handle fh;
    const bool is_face = cdt_2.is_face(cdt_2_v0, cdt_2_v1, cdt_2_v2, fh);
    if(is_face && fh->info().is_in_region != 0) {
      const int index_v0 = fh->index(cdt_2_v0);
      const bool reverse_orientation = (cdt_2_v2 == fh->vertex(T_3::ccw(index_v0)));
      return Oriented_face_of_cdt_2{fh, reverse_orientation};
    }
    else
      return std::nullopt;
  }

  auto edge_of_cdt_2(const CDT_2& cdt_2, const Vertex_handle va, const Vertex_handle vb) const
      -> std::optional<typename CDT_2::Edge>
  {
    auto v = [&, hint = CDT_2_face_handle{}](const auto& p) mutable {
      int i;
      typename CDT_2::Locate_type lt;
      const auto fh = cdt_2.locate(p, lt, i, hint);
      CGAL_assume(lt == CDT_2::VERTEX);
      hint = fh;
      return fh->vertex(i);
    };
    const auto cdt_2_v0 = v(this->point(va));
    const auto cdt_2_v1 = v(this->point(vb));
    CDT_2_face_handle fh;
    int edge_index;
    const bool is_edge = cdt_2.is_edge(cdt_2_v0, cdt_2_v1, fh, edge_index);
    if(is_edge) {
      typename CDT_2::Edge edge{fh, edge_index};
      // if(fh->vertex(cdt_2.cw(edge_index)) != cdt_2_v0) {
      //   edge = cdt_2.mirror_edge(edge);
      // }
      CGAL_assertion(edge.first->vertex(cdt_2.cw(edge.second)) == cdt_2_v0);
      CGAL_assertion(edge.first->vertex(cdt_2.ccw(edge.second)) == cdt_2_v1);
      return edge;
    }
    else
      return std::nullopt;
  }

  template <typename Tr1, typename Tr2, typename Vertex_handle1>
  static auto vertex_triple_is_facet_of_other_triangulation(
      const Tr1& tr, Vertex_handle1 v0, Vertex_handle1 v1, Vertex_handle1 v2, const Tr2& other_tr)
      -> std::optional<typename Tr2::Facet>
  {
    const auto p0 = tr.point(v0);
    const auto p1 = tr.point(v1);
    const auto p2 = tr.point(v2);
    auto v = [&, hint = typename Tr2::Cell_handle{}](const auto& p) mutable {
      int i, j;
      Locate_type lt;
      const auto c = other_tr.locate(p, lt, i, j, hint);
      CGAL_assume(lt == T_3::VERTEX);
      hint = c;
      return c->vertex(i);
    };
    typename Tr2::Cell_handle c;
    int i, j, k;
    const bool ok = other_tr.is_facet(v(p0), v(p1), v(p2), c, i, j, k);
    if(ok)
      return {typename Tr2::Facet(c, 6 - i - j - k)};
    else
      return {std::nullopt};
  };

  template <typename Vertex_map>
  auto triangulate_cavity([[maybe_unused]] CDT_3_face_index face_index,
                          const CDT_2& cdt_2,
                          std::set<Cell_handle> cells_of_cavity,
                          std::set<Facet>& facets_of_cavity_border,
                          Vertex_map& map_cavity_vertices_to_ambient_vertices,
                          std::set<Vertex_handle>& vertices_of_cavity) ///@TODO: not deterministic
  {
    struct {
      Constrained_Delaunay_triangulation_3 cavity_triangulation{};
      std::size_t number_of_added_vertices = 0;
    } result;
    auto& cavity_triangulation =  result.cavity_triangulation;
    CGAL::Unique_hash_map<Vertex_handle, Vertex_handle> vertex_map;

    auto insert_new_vertex = [&](Vertex_handle v, [[maybe_unused]] std::string_view extra = "") {
      const auto cavity_v = cavity_triangulation.insert(this->point(v));
      vertex_map[v] = cavity_v;
      map_cavity_vertices_to_ambient_vertices[cavity_v] = v;
#if CGAL_DEBUG_CDT_3 & 64
      std::cerr << std::format("inserted {}cavity vertex {:.6} -> {:.6}\n",
                               extra,
                               IO::oformat(cavity_v, with_point),
                               IO::oformat(v, with_point));
#endif
      ++result.number_of_added_vertices;
      return cavity_v;
    };

    for(const auto v : vertices_of_cavity) {
      insert_new_vertex(v);
    }

    std::vector<Facet> missing_faces;
    while(true) {
      missing_faces.clear();
      for(auto f : facets_of_cavity_border) {
        if(cells_of_cavity.contains(f.first)) {
          // internal facet, due to cavity growing
          continue;
        }
        const auto [v0, v1, v2] = this->make_vertex_triple(f);
        Cell_handle c;
        int i, j, k;
        if(!cavity_triangulation.is_facet(vertex_map[v0], vertex_map[v1], vertex_map[v2], c, i, j, k)) {
          missing_faces.push_back(f);
        }
      }
      if(missing_faces.empty()) {
        break;
      }
      for(auto [cell, facet_index] : missing_faces) {
        {
          if(cell->is_facet_constrained(facet_index)) {
            const auto polygon_contraint_id = cell->face_constraint_index(facet_index);
            auto f2d = cell->face_2(cdt_2, facet_index);
            f2d->info().missing_subface = true;
            CGAL_assertion(polygon_contraint_id != face_index);
            face_constraint_misses_subfaces.set(static_cast<std::size_t>(polygon_contraint_id));
          }
          auto [_, is_new_cell] = cells_of_cavity.insert(cell);
          if(!is_new_cell)
            continue;
        }
        const auto v3 = cell->vertex(facet_index);
        auto [_, v3_is_new_vertex] = vertices_of_cavity.insert(v3);
        if(v3_is_new_vertex) {
          insert_new_vertex(v3, "extra ");
        }
        for(int i = 0; i < 3; ++i) {
          Facet other_f{cell, this->vertex_triple_index(facet_index, i)};
          Facet mirror_f = this->mirror_facet(other_f);
          if(!cells_of_cavity.contains(mirror_f.first)) {
            facets_of_cavity_border.insert(mirror_f);
          }
        }
      }
    }
    CGAL_assertion(std::all_of(facets_of_cavity_border.begin(), facets_of_cavity_border.end(), [&](const auto& f) {
      const auto [v0, v1, v2] = this->make_vertex_triple(f);
      Cell_handle c;
      int i, j, k;
      return cavity_triangulation.is_facet(vertex_map[v0], vertex_map[v1], vertex_map[v2], c, i, j, k);
    }));
    return result;
  }

  bool restore_face(CDT_3_face_index face_index) {
    CDT_2& non_const_cdt_2 = face_cdt_2[face_index];
    const CDT_2& cdt_2 = non_const_cdt_2;
#if CGAL_DEBUG_CDT_3 && __has_include(<format>)
    std::cerr << std::format("restore_face({}): CDT_2 has {} vertices\n", face_index, cdt_2.number_of_vertices());
#endif // CGAL_DEBUG_CDT_3
    for(const auto& edge : cdt_2.finite_edges()) {
      const auto fh = edge.first;
      const auto i = edge.second;
      const auto va_3d = fh->vertex(cdt_2.cw(i))->info().vertex_handle_3d;
      const auto vb_3d = fh->vertex(cdt_2.ccw(i))->info().vertex_handle_3d;
      const bool is_3d = this->tds().is_edge(va_3d, vb_3d);
#if CGAL_DEBUG_CDT_3 & 128 && __has_include(<format>)
      std::cerr << std::format("Edge is 3D: {:6}  ({} , {})\n",
                                is_3d,
                                IO::oformat(this->point(va_3d)),
                                IO::oformat(this->point(vb_3d)));
#endif // CGAL_DEBUG_CDT_3
      CGAL_assertion(is_3d || !cdt_2.is_constrained(edge));
      fh->info().is_edge_also_in_3d_triangulation[unsigned(i)] = is_3d;
      const auto reverse_edge = cdt_2.mirror_edge(edge);
      reverse_edge.first->info().is_edge_also_in_3d_triangulation[unsigned(reverse_edge.second)] = is_3d;
    }
    std::set<CDT_2_face_handle> processed_faces;
    int region_count = 0;
    for(const CDT_2_face_handle fh : cdt_2.finite_face_handles()) {
      if(fh->info().is_outside_the_face) continue;
      CGAL_assertion(tr.is_facet(fh->vertex(0)->info().vertex_handle_3d,
                                 fh->vertex(1)->info().vertex_handle_3d,
                                 fh->vertex(2)->info().vertex_handle_3d) ||
                     (fh->info().missing_subface == true));
      if(false == fh->info().missing_subface) {
        continue;
      }
      Cell_handle c;
      int i, j, k;
      if(tr.is_facet(fh->vertex(0)->info().vertex_handle_3d,
                     fh->vertex(1)->info().vertex_handle_3d,
                     fh->vertex(2)->info().vertex_handle_3d, c, i, j, k))
      {
        const int facet_index = 6 - i - j - k;
        set_facet_constrained({c, facet_index}, face_index, fh);
        fh->info().missing_subface = false;
        continue;
      }
      if(processed_faces.contains(fh))
        continue;
      const auto fh_region = region(cdt_2, fh);
      processed_faces.insert(fh_region.begin(), fh_region.end());
      try {
        restore_subface_region(face_index, region_count++, cdt_2, fh_region);
      }
      catch(Next_region& e) {
        std::cerr << "ERROR: " << e.what() << " in sub-region " << (region_count - 1)
                  << " of face #F" << face_index << '\n';
#if CGAL_DEBUG_CDT_3 & 64 && __has_include(<format>)
        std::cerr << "  constrained edges are:\n";
        for(auto [c, index]: cdt_2.constrained_edges()) {
          const auto va = c->vertex(cdt_2.cw(index));
          const auto vb = c->vertex(cdt_2.ccw(index));
          const auto va_3d = va->info().vertex_handle_3d;
          const auto vb_3d = vb->info().vertex_handle_3d;
          std::cerr << std::format("    ({:.6} , {:.6})\n",
                                    IO::oformat(va_3d, with_point_and_info),
                                    IO::oformat(vb_3d, with_point_and_info));
        }
#endif // CGAL_DEBUG_CDT_3
        // return;
        const auto circ = CGAL::centroid(cdt_2.triangle(e.fh_2d));
        const auto other_fh = cdt_2.locate(circ, fh);
        if([&]() {
          for(int index = 0; index < 3; ++index) {
            if(!other_fh->is_constrained(index)) continue;
            const auto va = other_fh->vertex(cdt_2.cw(index));
            const auto vb = other_fh->vertex(cdt_2.ccw(index));
            const auto va_3d = va->info().vertex_handle_3d;
            const auto vb_3d = vb->info().vertex_handle_3d;
            const auto a = cdt_2.point(va);
            const auto b = cdt_2.point(vb);
            if(CGAL::angle(a, circ, b) == CGAL::OBTUSE) {
              const auto mid = CGAL::midpoint(a, b);
#if CGAL_DEBUG_CDT_3 & 64 && __has_include(<format>)
              std::cerr << std::format("Insert midpoint ({:.6}) of constrained edge ({:.6} , {:.6})\n",
                                        IO::oformat(mid),
                                        IO::oformat(va_3d, with_point_and_info),
                                        IO::oformat(vb_3d, with_point_and_info));
#endif // CGAL_DEBUG_CDT_3
              auto&& contexts = this->constraint_hierarchy.contexts_range(va_3d, vb_3d);
              CGAL_assertion(std::next(contexts.begin()) == contexts.end());
              const auto& context = *contexts.begin();
              const auto constraint_id = context.id();
              CGAL_assertion(constraint_id != Constraint_id{});
              this->insert_Steiner_point_on_subconstraint(mid, va_3d->cell(), {va_3d, vb_3d}, constraint_id,
                                                          insert_in_conflict_visitor);
              return false;
            }
          }
          return true;
        }()) {
#if CGAL_DEBUG_CDT_3 & 64 && __has_include(<format>)
          std::cerr << std::format("Inserting Steiner (circumcenter) point {} in non-coplanar face {}.\n",
                                    IO::oformat(circ),
                                    IO::oformat(cdt_2.triangle(e.fh_2d)));
#endif // CGAL_DEBUG_CDT_3
          const auto v = this->insert(circ);
          v->set_vertex_type(CDT_3_vertex_type::STEINER_IN_FACE);
          const auto v_2d = non_const_cdt_2.insert(circ, other_fh);
          v_2d->info().vertex_handle_3d = v;
        }
        return false;
      }
    }
    return true;
  }

public:
  void recheck_constrained_Delaunay() {
    for(int i = 0, end = face_constraint_misses_subfaces.size(); i < end; ++i) {
      search_for_missing_subfaces(i);
    }
  }

  void restore_constrained_Delaunay()
  {
    for(int i = 0, end = face_constraint_misses_subfaces.size(); i < end; ++i) {
      CDT_2& cdt_2 = face_cdt_2[i];
      fill_cdt_2(cdt_2, i);
      search_for_missing_subfaces(i);
    }
    cdt_2_are_initialized = true;
    const auto npos = face_constraint_misses_subfaces.npos;
    auto i = face_constraint_misses_subfaces.find_first();
    while(i != npos) {
      try {
        if(restore_face(i)) {
          face_constraint_misses_subfaces.reset(i);
        } else {
          std::cerr << "restore_face(" << i << ") incomplete, back to conforming...\n";
          Conforming_Dt::restore_Delaunay(insert_in_conflict_visitor);
          search_for_missing_subfaces(i);
        }
      }
      catch(PLC_error& e) {
        std::cerr << std::string("ERROR: PLC error with face #F") << std::to_string(e.face_index) + "\n";
      }
      i = face_constraint_misses_subfaces.find_first();
    }
  }

  static void write_region_to_OFF(std::ostream& out, const CDT_2& cdt_2) {
    out.precision(17);
    auto color_fn = [](CDT_2_face_handle fh_2d) -> CGAL::IO::Color {
      if(fh_2d->info().is_outside_the_face) return CGAL::IO::gray();
      if(fh_2d->info().is_in_region) return CGAL::IO::red();
      return CGAL::IO::green();
    };
    auto color_pmap = boost::make_function_property_map<CDT_2_face_handle>(color_fn);
    CGAL::IO::write_OFF(out, cdt_2, CGAL::parameters::face_color_map(color_pmap));
  }

  void write_region(std::ostream& out, const auto& region)
  {
    for(const auto fh_2d : region) {
      write_2d_triangle(out, fh_2d);
    }
  }

  void write_3d_triangulation_to_OFF(std::ostream& out, const Constrained_Delaunay_triangulation_3& tr) {
    write_facets(out, tr, tr.finite_facets());
  }

  void dump_3d_triangulation(CDT_3_face_index face_index,
                             int region_count,
                             std::string type,
                             const Constrained_Delaunay_triangulation_3& tr)
  {
    std::ofstream dump(std::string("dump_") + type + "_cavity_" + std::to_string(face_index) + "_" +
                       std::to_string(region_count) + ".off");
    dump.precision(17);
    write_3d_triangulation_to_OFF(dump, tr);
  }

  void dump_triangulation() const {
    std::ofstream dump("dump.binary.cgal");
    CGAL::IO::save_binary_file(dump, *this);
  }

  void dump_region(CDT_3_face_index face_index, int region_count, const CDT_2& cdt_2) {
    std::ofstream dump_region(std::string("dump_region_") + std::to_string(face_index) + "_" +
                              std::to_string(region_count) + ".off");
    write_region_to_OFF(dump_region, cdt_2);
  }

  void write_triangle(std::ostream &out,
                      Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
  {
    out.precision(17);
    out << "4"
        << " " << tr.point(v0) << " " << tr.point(v1) << " " << tr.point(v2)
        << " " << tr.point(v0) << '\n';
  }

  static void write_segment(std::ostream &out, Point_3 p0, Point_3 p1)
  {
    out.precision(17);
    out << "2" << " " << p0 << " " << p1 << '\n';
  }

  static void write_segment(std::ostream &out, Segment_3 seg) {
    write_segment(out, seg.source(), seg.target());
  }

  void write_segment(std::ostream& out, Vertex_handle v0, Vertex_handle v1)
  {
    write_segment(out, tr.point(v0), tr.point(v1));
  }

  void write_segment(std::ostream& out, Edge edge) {
    const auto [c, i, j] = edge;
    write_segment(out, c->vertex(i), c->vertex(j));
  }

  void dump_edge_link(std::string filename, Edge edge) {
    std::ofstream out(filename);
    out.precision(17);
    const auto [c, i, j] = edge;
    const Vertex_handle va = c->vertex(i);
    const Vertex_handle vb = c->vertex(j);
    auto cell_circ = this->incident_cells(edge), end = cell_circ;
    CGAL_assertion(cell_circ != nullptr);
    do {
      if(this->is_infinite(cell_circ)) {
        continue;
      }
      const auto index_va = cell_circ->index(va);
      const auto index_vb = cell_circ->index(vb);
      const auto index_vc = this->next_around_edge(index_va, index_vb);
      const auto index_vd = this->next_around_edge(index_vb, index_va);
      write_segment(out, cell_circ->vertex(index_vc), cell_circ->vertex(index_vd));
    } while(++cell_circ != end);
  }

  template <typename ...Args>
  void dump_segment(std::string filename, Args&& ...args)
  {
    std::ofstream out(filename);
    out.precision(17);
    write_segment(out, std::forward<Args>(args)...);
  }

  static auto export_facets_to_surface_mesh(const auto& tr, auto&& facets_range) {
    const auto size = std::distance(facets_range.begin(), facets_range.end());
    std::vector<typename Geom_traits::Point_3> points;
    points.reserve(size * 3);
    std::vector<std::array<std::size_t, 3>> facets;
    facets.reserve(size);

    for(std::size_t i = 0; const auto& [cell, facet_index] : facets_range) {
      const auto v0 = cell->vertex(T_3::vertex_triple_index(facet_index, 0));
      const auto v1 = cell->vertex(T_3::vertex_triple_index(facet_index, 1));
      const auto v2 = cell->vertex(T_3::vertex_triple_index(facet_index, 2));
      points.push_back(tr.point(v0));
      points.push_back(tr.point(v1));
      points.push_back(tr.point(v2));
      facets.push_back({i, i+1, i+2});
      i += 3;
    }
    CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, facets);
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, facets);
    Surface_mesh<Point_3> mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, facets, mesh);
    return mesh;
  }

  static void write_facets(std::ostream& out, const auto& tr, auto&& facets_range) {
    const auto mesh = export_facets_to_surface_mesh(tr, std::forward<decltype(facets_range)>(facets_range));
    CGAL::IO::write_OFF(out, mesh);
  }

  void dump_facets_of_cavity(CDT_3_face_index face_index, int region_count, std::string type, const auto& facets_range)
  {
    std::ofstream out(std::string("dump_facets_of_region_") + std::to_string(face_index) + "_" +
                      std::to_string(region_count) + "_" + type + ".off");
    out.precision(17);
    write_facets(out, *this, facets_range);
  }

  void write_2d_triangle(std::ostream &out, const CDT_2_face_handle fh)
  {
    const auto v0 = fh->vertex(0)->info().vertex_handle_3d;
    const auto v1 = fh->vertex(1)->info().vertex_handle_3d;
    const auto v2 = fh->vertex(2)->info().vertex_handle_3d;
    write_triangle(out, v0, v1, v2);
  }

  bool write_missing_subfaces_file(std::ostream& out) {
    const auto npos = face_constraint_misses_subfaces.npos;
    auto i = face_constraint_misses_subfaces.find_first();
    bool has_missing_subfaces = i != npos;
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
    return has_missing_subfaces;
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
  bool cdt_2_are_initialized = false;
  struct Face_edge {
    Constraint_id constraint_id;
    bool is_reverse = false;
  };
  std::vector<std::vector<Face_edge>> face_border;
  std::multimap<Constraint_id, CDT_3_face_index> constraint_to_faces;
  boost::dynamic_bitset<> face_constraint_misses_subfaces;
};

} // end CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
