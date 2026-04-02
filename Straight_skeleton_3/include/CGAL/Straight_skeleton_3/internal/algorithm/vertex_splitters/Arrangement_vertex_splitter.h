// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ARR_VERTEX_SPLITTER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ARR_VERTEX_SPLITTER_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/Straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <list>
#include <memory>
#include <vector>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

// These 2 visitors look generic, but they are not, it's really about a ValueType being a pointer
// and updating as to preserve the pointer information.
template <typename ValueType>
struct Range_updating_autoref_visitor
  : public CGAL::Polygon_mesh_processing::Autorefinement::Default_visitor,
    public CGAL::Polygon_mesh_processing::internal::Default_repair_PS_visitor
{
  using Autoref_visitor = CGAL::Polygon_mesh_processing::Autorefinement::Default_visitor;

  Range_updating_autoref_visitor(const std::vector<ValueType>& old_range,
                                 std::vector<ValueType>& new_range)
    : old_range_(old_range), new_range_(new_range) {
    new_range.reserve(old_range.size());
  }

  void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id) {
    // std::cout << "verbatim_triangle_copy " << tgt_id << " from " << src_id << std::endl;
    new_range_.resize(tgt_id + 1);
    new_range_[tgt_id] = old_range_[src_id];
  }

  void new_subtriangle (std::size_t tgt_id, std::size_t src_id) {
    // std::cout << "new_subtriangle " << tgt_id << " from " << src_id << std::endl;
    new_range_.resize(tgt_id + 1);
    new_range_[tgt_id] = old_range_[src_id];
  }

private:
  const std::vector<ValueType>& old_range_;
  std::vector<ValueType>& new_range_;
};

template <typename ValueType,
          typename BaseVisitor = CGAL::Polygon_mesh_processing::internal::Default_repair_PS_visitor>
struct Range_updating_repair_PS_visitor : public BaseVisitor
{
  Range_updating_repair_PS_visitor(std::vector<ValueType>& range,
                                    const BaseVisitor& base_visitor = BaseVisitor{})
    : BaseVisitor(base_visitor), range_(range) { }

  void swap(std::size_t pos_1, std::size_t pos_2) {
    std::cout << "swap(" << pos_1 << " " << pos_2 << ")" << std::endl;
    BaseVisitor::swap(pos_1, pos_2);
    CGAL_assertion(!(range_[pos_1]) || !(range_[pos_2]));
    std::swap(range_[pos_1], range_[pos_2]);
  }
  void duplicated_polygons(const std::vector<std::size_t>& duplicated_polygons) {
    ValueType ptr;
    // abusing ptr + knowing I'm keeping the first one
    for(std::size_t i=0; i<duplicated_polygons.size(); ++i) {
      if (range_[duplicated_polygons[i]]) {
        std::cout << "[] triangle " << duplicated_polygons[i] << " " << range_[duplicated_polygons[i]] << std::endl;
        if (ptr && range_[duplicated_polygons[i]]) {
          CGAL_assertion(ptr == range_[duplicated_polygons[i]]);
        }
        ptr = range_[duplicated_polygons[i]];
      }
    }
    range_[duplicated_polygons[0]] = ptr;
  }
  void resize(std::size_t new_size) {
    BaseVisitor::resize(new_size);
    range_.resize(new_size);
  }

private:
  std::vector<ValueType>& range_;
};

template <typename GeomTraits>
class Arr_vertex_splitter
  : public Abstract_vertex_splitter<GeomTraits>
{
  using Base = Abstract_vertex_splitter<GeomTraits>;
  using Arr_vertex_splitter_sptr = std::shared_ptr<Arr_vertex_splitter>;

private:
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  using Segment_3 = typename GeomTraits::Segment_3;
  using Ray_3 = typename GeomTraits::Ray_3;
  using Triangle_3 = typename GeomTraits::Triangle_3;
  using Plane_3 = typename GeomTraits::Plane_3;
  using Iso_cuboid_3 = typename GeomTraits::Iso_cuboid_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using Edge = typename Polyhedron::Edge;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using Skeleton_vertex_data = typename Polyhedron::Skeleton_vertex_data;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using Straight_skeleton_3 = CGAL::Straight_skeleton_3<GeomTraits>;

  using NodeWPtr = typename Straight_skeleton_3::NodeWPtr;
  using NodeSPtr = typename Straight_skeleton_3::NodeSPtr;

private:
  using Transformation = algorithm::Polyhedron_transformation<GeomTraits>;
  using Kernel_wrapper = kernel::Kernel_wrapper<GeomTraits>;

public:
  Arr_vertex_splitter()
  {
    this->type_ = Base::ARR_VERTEX_SPLITTER;
  }

  virtual ~Arr_vertex_splitter() { /*intentionally does nothing*/ }

  static Arr_vertex_splitter_sptr create()
  {
    return std::make_shared<Arr_vertex_splitter>();
  }

  static PolyhedronSPtr copy_vertex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    PolyhedronSPtr result = Polyhedron::create();
    VertexSPtr vertex_c = vertex->clone();
    result->add_vertex(vertex_c);
    FacetSPtr facet = vertex->first_facet();
    EdgeSPtr edge_first = vertex->find_edge(facet);
    EdgeSPtr edge;
    EdgeSPtr edge_prev;
    EdgeSPtr edge_prev_c;
    EdgeSPtr edge_first_c;
    FacetSPtr facet_next = facet;
    FacetSPtr facet_c;
    SkelFacetDataSPtr facet_c_data;
    while (edge != edge_first) {
      if (!edge) {
        edge = edge_first;
      }
      EdgeSPtr edge_c;
      VertexSPtr vertex_tgt_c;
      if (edge->source() == vertex) {
        vertex_tgt_c = edge->target()->clone();
      } else if (edge->target() == vertex) {
        vertex_tgt_c = edge->source()->clone();
      }
      result->add_vertex(vertex_tgt_c);
      edge_c = Edge::create(vertex_c, vertex_tgt_c);
      if (edge == edge_first) {
        edge_first_c = edge_c;
      }
      result->add_edge(edge_c);
      if (edge_prev_c) {
        facet_c = Facet::create();
        facet_c->set_plane(facet->get_plane());
        edge_prev_c->set_facet_L(facet_c);
        facet_c->add_edge(edge_prev_c);
        edge_c->set_facet_R(facet_c);
        facet_c->add_edge(edge_c);
        facet_c_data = Skeleton_facet_data::create(facet_c);
        facet_c_data->set_facet_origin(facet);
        if (facet->has_data()) {
          facet_c_data->set_speed(std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data())->get_speed());
        }
        result->add_facet(facet_c);
      }
      edge_prev_c = edge_c;
      edge_prev = edge;
      edge = edge->next(vertex);
      facet = facet_next;
      facet_next = edge->other(facet);
    }
    facet_c = Facet::create();
    facet_c->set_plane(facet->get_plane());
    edge_prev_c->set_facet_L(facet_c);
    facet_c->add_edge(edge_prev_c);
    edge_first_c->set_facet_R(facet_c);
    facet_c->add_edge(edge_first_c);
    facet_c_data = Skeleton_facet_data::create(facet_c);
    facet_c_data->set_facet_origin(facet);
    if (facet->has_data()) {
        facet_c_data->set_speed(std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data())->get_speed());
    }
    result->add_facet(facet_c);
    result->initialize_all_IDs();
    return result;
  }

  static PolyhedronSPtr apply(PolyhedronSPtr poly_split,
                              const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(poly_split);
    CGAL_SS3_DEBUG_SPTR(vertex);
    PolyhedronSPtr polyhedron = vertex->get_polyhedron();
    NodeWPtr node;
    if (vertex->has_data()) {
      SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
      node = data->get_wnode();
    }
    std::map<VertexSPtr, VertexSPtr> vertices;
    typename std::list<VertexSPtr>::iterator it_v = poly_split->vertices().begin();
    while (it_v != poly_split->vertices().end()) {
      VertexSPtr vertex_ps = *it_v++;
      if (vertex_ps->degree() > 1) {
        VertexSPtr vertex_vs = Vertex::create(vertex_ps->point());
        vertices[vertex_ps] = vertex_vs;
        polyhedron->add_vertex(vertex_vs);
        SkelVertexDataSPtr data = Skeleton_vertex_data::create(vertex_vs);
        data->set_wnode(node);
      }
    }

    typename std::list<EdgeSPtr>::iterator it_e = poly_split->edges().begin();
    while (it_e != poly_split->edges().end()) {
      EdgeSPtr edge_ps = *it_e++;
      VertexSPtr vertex_ps_src = edge_ps->source();
      VertexSPtr vertex_ps_tgt = edge_ps->target();
      FacetSPtr facet_ps_l = edge_ps->get_facet_L();
      FacetSPtr facet_ps_r = edge_ps->get_facet_R();
      if (vertex_ps_tgt->degree() == 1) {
        EdgeSPtr edge_vs;
        typename std::list<EdgeWPtr>::iterator it_ve = vertex->edges().begin();
        while (it_ve != vertex->edges().end()) {
          EdgeWPtr edge_wptr = *it_ve++;
          if (EdgeSPtr edge = edge_wptr.lock()) {
            if (edge->source()->point() == vertex_ps_tgt->point() ||
                edge->target()->point() == vertex_ps_tgt->point()) {
              edge_vs = edge;
              break;
            }
          }
        }
        VertexSPtr vertex_vs = vertices[vertex_ps_src];
        if (edge_vs->source() == vertex) {
          edge_vs->replace_vertex_src(vertex_vs);
        } else {
          edge_vs->replace_vertex_tgt(vertex_vs);
        }
        if (!edge_vs->get_facet_L()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_L()->add_vertex(vertex_vs);
        }
        if (!edge_vs->get_facet_R()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_R()->add_vertex(vertex_vs);
        }
      } else if (vertex_ps_src->degree() == 1) {
        EdgeSPtr edge_vs;
        typename std::list<EdgeWPtr>::iterator it_ve = vertex->edges().begin();
        while (it_ve != vertex->edges().end()) {
          EdgeWPtr edge_wptr = *it_ve++;
          if (EdgeSPtr edge = edge_wptr.lock()) {
            if (edge->source()->point() == vertex_ps_src->point() ||
                edge->target()->point() == vertex_ps_src->point()) {
              edge_vs = edge;
              break;
            }
          }
        }
        VertexSPtr vertex_vs = vertices[vertex_ps_tgt];
        if (edge_vs->source() == vertex) {
          edge_vs->replace_vertex_src(vertex_vs);
        } else {
          edge_vs->replace_vertex_tgt(vertex_vs);
        }
        if (!edge_vs->get_facet_L()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_L()->add_vertex(vertex_vs);
        }
        if (!edge_vs->get_facet_R()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_R()->add_vertex(vertex_vs);
        }
      } else {
        EdgeSPtr edge_vs = Edge::create(vertices[vertex_ps_src], vertices[vertex_ps_tgt]);
        SkelFacetDataSPtr data_l = std::dynamic_pointer_cast<Skeleton_facet_data>(facet_ps_l->get_data());
        SkelFacetDataSPtr data_r = std::dynamic_pointer_cast<Skeleton_facet_data>(facet_ps_r->get_data());
        FacetSPtr facet_vs_l = data_l->get_facet_origin();
        FacetSPtr facet_vs_r = data_r->get_facet_origin();
        edge_vs->set_facet_L(facet_vs_l);
        edge_vs->set_facet_R(facet_vs_r);
        facet_vs_l->add_edge(edge_vs);
        facet_vs_r->add_edge(edge_vs);
        polyhedron->add_edge(edge_vs);
      }
    }

    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facet->remove_vertex(vertex);
      }
    }
    polyhedron->remove_vertex(vertex);
    return polyhedron;
  }

public:
  // Create a sufficiently large bounding box containing all plane intersections
  static Iso_cuboid_3 compute_intersection_bbox(const VertexSPtr& vertex)
  {
    std::vector<Plane_3> all_planes; // @todo avoid copies
    for(const auto& facet_wptr : vertex->facets())
    {
      if(FacetSPtr facet = facet_wptr.lock())
      {
        all_planes.push_back(facet->get_plane());

        // Also add shifted planes
        if(facet->has_data())
        {
          auto skel_data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
          if(skel_data && skel_data->get_speed() > 0) {
            all_planes.push_back(Transformation::shift_plane(facet, -1));
          }
        }
      }
    }

    std::cout << "all planes: " << all_planes.size() << std::endl;

    std::ofstream inter_out("results/3-inter.xyz");
    inter_out.precision(17);

    // @todo obviously !slightly! suboptimal
    Bbox_3 bbox = vertex->point().bbox();
    for(const auto& pl_i : all_planes) {
      for(const auto& pl_j : all_planes) {
        for(const auto& pl_k : all_planes) {
          if(pl_i == pl_j || pl_i == pl_k || pl_j == pl_k) {
            continue;
          }

          std::optional<Point_3> opt_pt = Kernel_wrapper::intersection(pl_i, pl_j, pl_k);
          if (opt_pt.has_value()) { // might not exist because we can take a plane and a // plane
            bbox += opt_pt->bbox();
            inter_out << *(opt_pt) << "\n";
          }
        }
      }
    }

    bbox.scale(1.5);

    return Iso_cuboid_3(bbox);
  }

  static void get_clipped_plane_faces(const VertexSPtr vertex,
                                      const Iso_cuboid_3& bbox,
                                      std::vector<Point_3>& points,
                                      std::vector<std::vector<std::size_t> >& triangles,
                                      std::vector<FacetSPtr>& triangle_2_sptr)
  {
    using Vector_3 = typename GeomTraits::Vector_3;
    const Point_3& center = vertex->point();

    for(const auto& facet_wptr : vertex->facets())
    {
      if(FacetSPtr facet = facet_wptr.lock())
      {
        const Plane_3& plane = facet->get_plane();

        std::vector<Point_3> local_range;
        auto res = CGAL::intersection(bbox, plane);
        if (!res) {
          // Should not happen, as bbox is constructed to contain all intersections
          std::cerr << "no intersection between plane and bbox" << std::endl;
          CGAL_assertion(false);
          std::exit(1);
        } else if (const Triangle_3* itr = std::get_if<Triangle_3>(&*res)) {
          for (int i=0; i<3; ++i) {
            local_range.push_back((*itr)[i]);
          }
        } else if (const std::vector<Point_3>* ir = std::get_if<std::vector<Point_3> >(&*res)) {
          for (const Point_3& p : *ir) {
            local_range.push_back(p);
          }
        } else {
          std::cerr << "plane/bbox intersection is not a polygon" << std::endl;
          CGAL_assertion(false);
          std::exit(1);
        }

        // Ensure orientation: normal of local_range must match plane's orientation
        if(local_range.size() >= 3) {
          Vector_3 plane_normal = plane.orthogonal_vector();
          Vector_3 tri_normal = CGAL::cross_product(local_range[1] - local_range[0], local_range[2] - local_range[1]);
          if(tri_normal * plane_normal < 0) {
            std::reverse(local_range.begin(), local_range.end());
          }
        }

        // Insert extremities (projections) into local_range at the correct place
        const Point_3& prev_pt = vertex->prev(facet)->point();
        const Point_3& v_pt = vertex->point();
        const Point_3& next_pt = vertex->next(facet)->point();

        Vector_3 prev_dir = prev_pt - v_pt;
        Vector_3 next_dir = next_pt - v_pt;
        auto res_prev = CGAL::intersection(Ray_3{v_pt, prev_dir}, bbox);
        auto res_next = CGAL::intersection(Ray_3{v_pt, next_dir}, bbox);
        CGAL_assertion(res_prev && res_next);

        auto get_ray_bbox_extremity = [&](const auto& res, const Point_3& src) -> std::optional<Point_3> {
          if(const Segment_3* seg = std::get_if<Segment_3>(&*res)) {
            if(seg->source() == src)
              return seg->target();
            else if(seg->target() == src)
              return seg->source();
            else
              return std::nullopt;
          }
          return std::nullopt;
        };

        std::optional<Point_3> opt_prev = get_ray_bbox_extremity(res_prev, v_pt);
        std::optional<Point_3> opt_next = get_ray_bbox_extremity(res_next, v_pt);
        CGAL_assertion(opt_prev && opt_next);

        // check linearly to find in which segment the point belongs
        // (it belongs by construction)
        auto insert_in_order = [&local_range](const Point_3& new_p)
        {
          for(auto it = local_range.begin(); it != local_range.end(); ++it) {
            auto it_next = std::next(it);
            if(it_next == local_range.end()) {
              it_next = local_range.begin();
            }
            if(new_p != *it && new_p != *it_next &&
               CGAL::collinear(*it, new_p, *it_next) &&
               CGAL::collinear_are_strictly_ordered_along_line(*it, new_p, *it_next))
            {
              // std::cout << "insert " << new_p << " between " << *it << " and " << *it_next << std::endl;
              local_range.insert(it_next, new_p);
              break;
            }
          }
        };

        insert_in_order(*opt_prev);
        insert_in_order(*opt_next);

        // Triangulate by fanning from the center vertex (vertex->point())
        std::size_t center_idx = points.size();
        points.push_back(center);

        std::size_t base_idx = points.size();
        for(const Point_3& p : local_range) {
          points.push_back(p);
        }

        // --- Sector logic ---
        // Get prev/next points for this facet at the center vertex
        Vector_3 n = plane.orthogonal_vector();

        // Orientation planes: through center, normal is n x (prev-center) and n x (next-center)
        Vector_3 v_prev = prev_pt - center;
        Vector_3 v_next = next_pt - center;
        Vector_3 n_prev = CGAL::cross_product(n, v_prev);
        Vector_3 n_next = CGAL::cross_product(n, v_next);
        Plane_3 plane_prev(center, n_prev);
        Plane_3 plane_next(center, n_next);

        // Determine if angle at center is > 180°
        // If next is to the left of prev (in the facet's orientation), angle < 180°
        // If next is to the right of prev, angle > 180°
        // Use orientation of (center, prev, next) with normal n
        bool angle_gt_180 = (CGAL::orientation(center, next_pt, prev_pt, center + n) == CGAL::NEGATIVE);
        std::cout << "center = " << center << std::endl;
        std::cout << "prev_pt = " << prev_pt << std::endl;
        std::cout << "next_pt = " << next_pt << std::endl;
        std::cout << "center + n = " << center + n << std::endl;
        std::cout << "angle_gt_180 = " << angle_gt_180 << std::endl;

        for(std::size_t i=0; i<local_range.size(); ++i) {
          std::size_t i1 = base_idx + i;
          std::size_t i2 = base_idx + ((i+1)%local_range.size());

          // Compute midpoint of the two extremities
          const Point_3& p1 = points[i1];
          const Point_3& p2 = points[i2];
          Point_3 mid = CGAL::midpoint(p1, p2);
          std::cout << "test: " << mid << std::endl;

          // Test if mid is between the two orientation planes
          bool on_pos_side_prev = (plane_prev.oriented_side(mid) == CGAL::ON_NEGATIVE_SIDE);
          bool on_pos_side_next = (plane_next.oriented_side(mid) == CGAL::ON_POSITIVE_SIDE);
          std::cout << "on_pos_side_prev = " << on_pos_side_prev << std::endl;
          std::cout << "on_pos_side_next = " << on_pos_side_next << std::endl;

          bool in_sector = (on_pos_side_prev && on_pos_side_next);
          // If angle > 180°, invert logic: triangles inside are NOT part of the facet
          bool is_facet_triangle = angle_gt_180 ? !in_sector : in_sector;

          triangles.push_back({center_idx, i1, i2});
          triangle_2_sptr.push_back(is_facet_triangle ? facet : nullptr);
        }
      }
    }
  }

  static void get_clipped_shifted_plane_faces(const VertexSPtr vertex,
                                              const Iso_cuboid_3& bbox,
                                              std::vector<Point_3>& points,
                                              std::vector<std::vector<std::size_t> >& triangles,
                                              std::vector<FacetSPtr>& triangle_2_sptr)
  {
    for(const auto& facet_wptr : vertex->facets())
    {
      if(FacetSPtr facet = facet_wptr.lock())
      {
        const Plane_3 plane = Transformation::shift_plane(facet, -1);

        std::vector<Point_3> local_range;
        auto res = CGAL::intersection(bbox, plane);
        if (!res) {
          // that should not happen since we have every 3-plane intersection in the bbox
          std::cerr << "no inter plane/bbox" << std::endl;
          CGAL_assertion(false);
          std::exit(1);
        } else if (const Triangle_3* itr = std::get_if<Triangle_3>(&*res)) {
          std::cout << "triangle inter" << std::endl;
          for (int i=0; i<3; ++i) {
            local_range.push_back(itr->operator[](i));
          }
        } else if (const std::vector<Point_3>* ir = std::get_if<std::vector<Point_3> >(&*res)) {
          std::cout << "range inter (" << ir->size() << ")" << std::endl;
          for (const Point_3& p : *ir) {
            local_range.push_back(p);
          }
        } else {
          std::cerr << "plane/bbox inter shouldn't be a point or a segment" << std::endl;
          CGAL_assertion(false);
          std::exit(1);
        }

        // Check orientation: normal of local_range must match plane's orientation
        if(local_range.size() >= 3) {
          Vector_3 plane_normal = plane.orthogonal_vector();
          Vector_3 tri_normal = CGAL::cross_product(local_range[1] - local_range[0], local_range[2] - local_range[1]);
          if(tri_normal * plane_normal < 0) {
            std::reverse(local_range.begin(), local_range.end());
          }
        }

        // triangulate the range to get faces
        CGAL_assertion(local_range.size() >= 3);

        // Add all points to the global point list
        std::size_t base_idx = points.size();

        for(const Point_3& p : local_range) {
          points.push_back(p);
        }

        // Create triangles by fanning from the first vertex
        for(std::size_t i=1; i<local_range.size() - 1; ++i) {
          triangles.push_back({base_idx, base_idx + i, base_idx + i + 1});
          triangle_2_sptr.push_back(nullptr);
        }
      }
    }
  }

public:
  virtual bool split_vertex(const VertexSPtr& vertex)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    CGAL_SS3_DEBUG_SPTR(vertex);

    if(vertex->degree() <= 3)
      return true;

    vertex->sort();

    // soup to be used in the arrangement
    std::vector<Point_3> points;
    std::vector<std::vector<std::size_t> > triangles;
    std::vector<FacetSPtr> triangle_2_sptr;

    // Compute bounding box for intersections
    Iso_cuboid_3 bbox = compute_intersection_bbox(vertex);

    // Get the triangles from the base planes, with tagged faces
    get_clipped_plane_faces(vertex, bbox, points, triangles, triangle_2_sptr);

    std::cout << points.size() << " points, " << triangles.size() << " triangles [base]" << std::endl;

    // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
    {
      std::map<FacetSPtr, std::vector<std::vector<std::size_t> > > fsptr_to_faces;
      for (std::size_t i=0; i<triangles.size(); ++i) {
        FacetSPtr fsptr = triangle_2_sptr[i];
        if (fsptr) {
          fsptr_to_faces[fsptr].push_back(triangles[i]);
        }
      }
      std::size_t idx = 0;
      for (const auto& kv : fsptr_to_faces) {
        std::ostringstream oss;
        oss << "results/arr_base_fsptr_" << idx << ".off";
        CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
        ++idx;
      }
    }

    // Get the triangles for the shifted base planes
    get_clipped_shifted_plane_faces(vertex, bbox, points, triangles, triangle_2_sptr);

    std::cout << points.size() << " points, " << triangles.size() << " triangles [+shift]" << std::endl;

    CGAL::IO::write_OFF("results/arr_base_all_planes.off", points, triangles, CGAL::parameters::stream_precision(17));

    // Add the bounding box faces (@todo probably not necessary)
    std::size_t base_idx = points.size();
    points.push_back(Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin())); // v0
    points.push_back(Point_3(bbox.xmax(), bbox.ymin(), bbox.zmin())); // v1
    points.push_back(Point_3(bbox.xmax(), bbox.ymax(), bbox.zmin())); // v2
    points.push_back(Point_3(bbox.xmin(), bbox.ymax(), bbox.zmin())); // v3
    points.push_back(Point_3(bbox.xmin(), bbox.ymin(), bbox.zmax())); // v4
    points.push_back(Point_3(bbox.xmax(), bbox.ymin(), bbox.zmax())); // v5
    points.push_back(Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax())); // v6
    points.push_back(Point_3(bbox.xmin(), bbox.ymax(), bbox.zmax())); // v7

    // bottom face
    triangles.push_back({base_idx+0, base_idx+2, base_idx+1});
    triangles.push_back({base_idx+0, base_idx+3, base_idx+2});
    // top face
    triangles.push_back({base_idx+4, base_idx+5, base_idx+6});
    triangles.push_back({base_idx+4, base_idx+6, base_idx+7});
    // front face
    triangles.push_back({base_idx+0, base_idx+1, base_idx+5});
    triangles.push_back({base_idx+0, base_idx+5, base_idx+4});
    // back face
    triangles.push_back({base_idx+2, base_idx+3, base_idx+7});
    triangles.push_back({base_idx+2, base_idx+7, base_idx+6});
    // left face
    triangles.push_back({base_idx+0, base_idx+4, base_idx+7});
    triangles.push_back({base_idx+0, base_idx+7, base_idx+3});
    // right face
    triangles.push_back({base_idx+1, base_idx+2, base_idx+6});
    triangles.push_back({base_idx+1, base_idx+6, base_idx+5});

    triangle_2_sptr.resize(triangles.size(), nullptr);

    std::cout << points.size() << " points, " << triangles.size() << " triangles [+bbox]" << std::endl;

    for (std::size_t i=0; i<triangle_2_sptr.size(); ++i) {
      std::cout << "triangle " << i << " ==> ptr: " << triangle_2_sptr[i] << std::endl;
    }

    PMP::merge_duplicate_points_in_polygon_soup(points, triangles);
    CGAL::IO::write_OFF("results/arr_soup.off", points, triangles, CGAL::parameters::stream_precision(17));

    CGAL_assertion(triangles.size() == triangle_2_sptr.size());

    std::cout << "autorefining..." << std::endl;
    std::vector<FacetSPtr> updated_triangle_2_sptr;
    Range_updating_autoref_visitor<FacetSPtr> autoref_visitor(triangle_2_sptr, updated_triangle_2_sptr);

    PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::visitor(autoref_visitor));
    triangle_2_sptr = std::move(updated_triangle_2_sptr);
    CGAL_assertion(triangles.size() == triangle_2_sptr.size());

    CGAL::IO::write_OFF("results/arr_arr.off", points, triangles, CGAL::parameters::stream_precision(17));

    for (std::size_t i=0; i<triangles.size(); ++i) {
      std::cout << "[4] triangle " << i << " ptr: " << triangle_2_sptr[i] << std::endl;
    }

    // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
    {
      std::map<FacetSPtr, std::vector<std::vector<std::size_t> > > fsptr_to_faces;
      for (std::size_t i=0; i<triangles.size(); ++i) {
        FacetSPtr fsptr = triangle_2_sptr[i];
        if (fsptr) {
          fsptr_to_faces[fsptr].push_back(triangles[i]);
        }
      }
      std::size_t idx = 0;
      for (const auto& kv : fsptr_to_faces) {
        std::ostringstream oss;
        oss << "results/arr_arr_fsptr_" << idx << ".off";
        CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
        ++idx;
      }
    }

    std::cout << "repairing..." << std::endl;
    Range_updating_repair_PS_visitor<FacetSPtr> repair_ps_visitor(triangle_2_sptr);
    PMP::merge_duplicate_polygons_in_polygon_soup(points, triangles,
                                                  CGAL::parameters::visitor(repair_ps_visitor)
                                                                   .erase_all_duplicates(false) /*keep one*/
                                                                   .require_same_orientation(false)
                                                                   .verbose(true));
    CGAL::IO::write_OFF("results/arr_repaired.off", points, triangles, CGAL::parameters::stream_precision(17));

    // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
    {
      std::map<FacetSPtr, std::vector<std::vector<std::size_t> > > fsptr_to_faces;
      for (std::size_t i=0; i<triangles.size(); ++i) {
        FacetSPtr fsptr = triangle_2_sptr[i];
        if (fsptr) {
          fsptr_to_faces[fsptr].push_back(triangles[i]);
        }
      }
      std::size_t idx = 0;
      for (const auto& kv : fsptr_to_faces) {
        std::ostringstream oss;
        oss << "results/arr_repaired_fsptr_" << idx << ".off";
        CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
        ++idx;
      }
    }

    CGAL_assertion(triangles.size() == triangle_2_sptr.size());

    std::cout << "converting..." << std::endl;
    PMP::orient_polygon_soup(points, triangles);
    CGAL::IO::write_OFF("results/arr_oriented.off", points, triangles, CGAL::parameters::stream_precision(17));

    CGAL_assertion(triangles.size() == triangle_2_sptr.size());

    using Mesh = CGAL::Surface_mesh<Point_3>;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
    using face_iterator = typename boost::graph_traits<Mesh>::face_iterator;

    Mesh sm;
    auto fpm = sm.template add_property_map<face_descriptor, FacetSPtr>("f:fsptr").first;

    auto out = boost::make_function_output_iterator(
      [fpm, &triangle_2_sptr](const std::pair<std::size_t, face_descriptor>& pid_f) {
        std::cout << "add " << triangle_2_sptr[pid_f.first] << " to " << pid_f.second << std::endl;
        put(fpm, pid_f.second, triangle_2_sptr[pid_f.first]);
      });
    PMP::polygon_soup_to_polygon_mesh(points, triangles, sm,
                                      CGAL::parameters::polygon_to_face_output_iterator(out));
    CGAL::IO::write_OFF("results/arr_mesh.off", sm, CGAL::parameters::stream_precision(17));

    for (face_descriptor f : faces(sm)) {
      std::cout << "face " << f << " has shared ptr: " << get(fpm, f) << std::endl;
    }

    std::cout << "remesh planar stuff..." << std::endl;
    Mesh sm_remeshed;
    auto fpm_remeshed = sm.template add_property_map<face_descriptor, FacetSPtr>("f:fsptr").first;

    auto f_patch_pid_pm = sm.template add_property_map<face_descriptor, std::size_t>("f:patch_id").first;
    auto f_patch_pid_pm_remeshed = sm_remeshed.template add_property_map<face_descriptor, std::size_t>("f:patch_id").first;

    PMP::remesh_planar_patches(sm, sm_remeshed,
                               CGAL::parameters::face_patch_map(f_patch_pid_pm),
                               CGAL::parameters::face_patch_map(f_patch_pid_pm_remeshed)
                                                .do_not_triangulate_faces(true));
    CGAL::IO::write_OFF("results/arr_remeshed.off", sm_remeshed, CGAL::parameters::stream_precision(17));

    // now transport the info between the meshes
    std::vector<FacetSPtr> patch_id_2_sptr;
    for (face_descriptor f : faces(sm)) {
      const std::size_t pid = get(f_patch_pid_pm, f);
      const FacetSPtr& f_sptr = get(fpm, f);
      if (pid >= patch_id_2_sptr.size()) {
        patch_id_2_sptr.resize(pid + 1, nullptr);
      }
      std::cout << "f " << f << " on patch " << pid << " has sptr " << f_sptr << std::endl;
      if (patch_id_2_sptr[pid]) {
        std::cout << "  existing sptr " << patch_id_2_sptr[pid] << " already associated to patch" << std::endl;
        CGAL_assertion(patch_id_2_sptr[pid] == f_sptr);
      }
      patch_id_2_sptr[pid] = f_sptr;
    }

    for (face_descriptor f : faces(sm_remeshed)) {
      std::size_t pid = get(f_patch_pid_pm_remeshed, f);
      FacetSPtr f_sptr = patch_id_2_sptr[pid];
      put(fpm_remeshed, f, f_sptr);
    }

    for (face_descriptor f : faces(sm_remeshed)) {
      std::cout << "remeshed face " << f << " has shared ptr: " << get(fpm_remeshed, f) << std::endl;
    }

    // back to soup for clarity, we want to navigate into an arrangement of polygons
    points.clear();
    std::vector<std::vector<std::size_t> > polygons;
    PMP::polygon_mesh_to_polygon_soup(sm_remeshed, points, polygons);
    CGAL::IO::write_OFF("results/arr_polygons.off", points, triangles, CGAL::parameters::stream_precision(17));

    // @todo bit of an abuse here because we assume that pm_to_ps puts into the same order as faces(g)
    triangle_2_sptr.assign(faces(sm_remeshed).size(), nullptr);
    face_iterator fit = faces(sm_remeshed).begin();
    for (std::size_t i=0; i<faces(sm_remeshed).size(); ++i, ++fit) {
      triangle_2_sptr[i] = get(fpm_remeshed, *fit);
    }
    CGAL_assertion(polygons.size() == triangle_2_sptr.size());

    // dump only the polygons with a non nullptr
    {
      std::vector<std::vector<std::size_t> > input_polygons;
      for (std::size_t i=0; i<polygons.size(); ++i) {
        if (triangle_2_sptr[i]) {
          input_polygons.push_back(polygons[i]);
        }
      }

      std::cout << input_polygons.size() << " input polygons" << std::endl;
      CGAL::IO::write_OFF("results/arr_recovered_input_polygons.off", points, input_polygons, CGAL::parameters::stream_precision(17));
    }

    // Apply the split
    // apply(poly_split, vertex);

    std::exit(1);

    return true;
  }

  virtual std::string to_string() const
  {
    std::stringstream sstr;
    sstr << "Arr_vertex_splitter()";
    return sstr.str();
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ARR_VERTEX_SPLITTER_H */
