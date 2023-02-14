#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <CGAL/create_weighted_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/create_weighted_offset_polygons_2.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_weighted_straight_skeleton_2.h>

#include "print.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h> // could be avoided
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#include <boost/container/flat_map.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace SS = CGAL::CGAL_SS_i;
namespace PMP = CGAL::Polygon_mesh_processing;

// Kernel choice:
// EPICK: Robust and fast
// EPECK_with_sqrt: Exact and slow
// EPECK: More robust, and less slow than EPECK_with_sqrt

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
// using K = CGAL::Exact_predicates_exact_constructions_kernel; // @fixme need to implement ceil(Lazy_exact_NT first)
// using K = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;

using FT = K::FT;
using Point_2 = K::Point_2;
using Segment_2 = K::Segment_2;
using Line_2 = K::Line_2;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;

using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

using Offset_polygons = std::vector<boost::shared_ptr<Polygon_2> >;
using Offset_polygons_with_holes = std::vector<boost::shared_ptr<Polygon_with_holes_2> >;

using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
using Straight_skeleton_2_ptr = boost::shared_ptr<Straight_skeleton_2>;

using SS_Vertex_const_handle = Straight_skeleton_2::Vertex_const_handle;
using SS_Halfedge_const_handle = Straight_skeleton_2::Halfedge_const_handle;

using HDS = Straight_skeleton_2::Base;
using HDS_Vertex_const_handle = HDS::Vertex_const_handle;
using HDS_Halfedge_const_handle = HDS::Halfedge_const_handle;
using HDS_Face_handle = HDS::Face_handle;

// Standard CDT2 for the horizontal (z constant) faces
using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K>;
using Vbb = CGAL::Triangulation_vertex_base_2<K, Vb>;
using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
using TDS = CGAL::Triangulation_data_structure_2<Vb,Fb>;
using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
using CDT_Vertex_handle = CDT::Vertex_handle;
using CDT_Face_handle = CDT::Face_handle;

// Projection CDT2 for the lateral faces
using PK = CGAL::Projection_traits_3<K>;
using PVbb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, PK>;
using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
using PCDT_Vertex_handle = PCDT::Vertex_handle;
using PCDT_Face_handle = PCDT::Face_handle;

using Offset_builder_traits = CGAL::Polygon_offset_builder_traits_2<K>;

using Mesh = CGAL::Surface_mesh<Point_3>;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Below is to handle vertical slabs.

// @todo Maybe this postprocessing is not really necessary?...
#define CGAL_SLS_SNAP_TO_VERTICAL_SLABS
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS

// The purpose of this visitor is to snap back almost-vertical (see preprocessing_weights()) edges
// to actual vertical slabs.
Point_2 snap_point_to_contour_halfedge_plane(const Point_2& op,
                                             SS_Halfedge_const_handle ch)
{
  const SS_Vertex_const_handle sv = ch->opposite()->vertex();
  const SS_Vertex_const_handle tv = ch->vertex();

  if(sv->point().x() == tv->point().x())
  {
    // vertical edge
    // std::cout << "vertical edge, snapping " << op << " to " << sv->point().x() << " "  << op.y() << std::endl;
    return { sv->point().x(), op.y() };
  }
  else if(sv->point().y() == tv->point().y())
  {
    // horizontal edge
    // std::cout << "horizontal edge, snapping " << op << " to " << op.x() << " " << sv->point().y() << std::endl;
    return { op.x(), sv->point().y() };
  }
  else
  {
    // Project orthogonally onto the halfedge
    // @fixme, the correct projection should be along the direction of the other offset edge sharing this point
    Segment_2 s { sv->point(), tv->point() };
    boost::optional<Line_2> line = SS::compute_weighted_line_coeffC2(s, FT(1)); // the weight does not matter
    CGAL_assertion(bool(line)); // otherwise the skeleton would have failed already

    FT px, py;
    CGAL::line_project_pointC2(line->a(),line->b(),line->c(), op.x(),op.y(), px,py);
    // std::cout << "snapping " << op << " to " << px << " " << py << std::endl;
    return { px, py };
  }
};

void snap_skeleton_vertex(HDS_Halfedge_const_handle hds_h,
                          HDS_Halfedge_const_handle contour_h,
                          std::map<Point_2, Point_2>& snapped_positions)
{
  HDS_Vertex_const_handle hds_tv = hds_h->vertex();

  // this re-applies snapping towards contour_h even if the point was already snapped towards another contour
  auto insert_result = snapped_positions.emplace(hds_tv->point(), hds_tv->point());
  insert_result.first->second = snap_point_to_contour_halfedge_plane(insert_result.first->second, contour_h);

  // std::cout << "snap_skeleton_vertex(V" << hds_tv->id() << " pt: " << hds_h->vertex()->point() << ")"
  //           << " to " << insert_result.first->second << std::endl;
};

template <typename PointRange>
void apply_snapping(PointRange& points,
                    const std::map<Point_2, Point_2>& snapped_positions)
{
  for(Point_3& p3 : points)
  {
    auto it = snapped_positions.find({ p3.x(), p3.y() });
    if(it != snapped_positions.end())
      p3 = Point_3{it->second.x(), it->second.y(), p3.z()};
  }
}
#endif

class Skeleton_offset_correspondence_builder_visitor
  : public CGAL::Default_polygon_offset_builder_2_visitor<Offset_builder_traits, Straight_skeleton_2>
{
  using Base = CGAL::Default_polygon_offset_builder_2_visitor<Offset_builder_traits, Straight_skeleton_2>;

public:
  Skeleton_offset_correspondence_builder_visitor(const Straight_skeleton_2& ss,
                                                 std::unordered_map<HDS_Halfedge_const_handle, Point_2>& offset_points
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                                               , const FT vertical_weight
                                               , std::map<Point_2, Point_2>& snapped_positions
#endif
                                                 )
    : m_ss(ss)
    , m_offset_points(offset_points)
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    , m_vertical_weight(vertical_weight)
    , m_snapped_positions(snapped_positions)
#endif
  { }

public:
  void on_offset_contour_started() const
  {
    // std::cout << "~~ new contour ~~" << std::endl;
  }

  // can't modify the position yet because we need arrange_polygons() to still work properly
  void on_offset_point(const Point_2& op,
                       SS_Halfedge_const_handle hook) const
  {
    CGAL_assertion(hook->is_bisector());

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
   // @fixme technically, one could create a polygon thin-enough w.r.t. the max weight value such that
   // there is a skeleton vertex that wants to be snapped to two different sides...
    CGAL_assertion(m_snapped_positions.count(op) == 0);

    HDS_Halfedge_const_handle canonical_hook = (hook < hook->opposite()) ? hook : hook->opposite();

    SS_Halfedge_const_handle contour_h1 = hook->defining_contour_edge();
    CGAL_assertion(contour_h1->opposite()->is_border());
    SS_Halfedge_const_handle contour_h2 = hook->opposite()->defining_contour_edge();
    CGAL_assertion(contour_h2->opposite()->is_border());

    const bool is_h1_vertical = (contour_h1->weight() == m_vertical_weight);
    const bool is_h2_vertical = (contour_h2->weight() == m_vertical_weight);

    // this can happen when the offset is passing through vertices
    m_offset_points[canonical_hook] = op;

    // if both are vertical, it's the common vertex (which has to exist)
    if(is_h1_vertical && is_h2_vertical)
    {
      CGAL_assertion(contour_h1->vertex() == contour_h2->opposite()->vertex() ||
                     contour_h2->vertex() == contour_h1->opposite()->vertex());
      if(contour_h1->vertex() == contour_h2->opposite()->vertex())
        m_snapped_positions[op] = contour_h1->vertex()->point();
      else
        m_snapped_positions[op] = contour_h2->vertex()->point();
    }
    else if(is_h1_vertical)
    {
      m_snapped_positions[op] = snap_point_to_contour_halfedge_plane(op, contour_h1);
    }
    else if(is_h2_vertical)
    {
      m_snapped_positions[op] = snap_point_to_contour_halfedge_plane(op, contour_h2);
    }
#endif
  }

private:
  const Straight_skeleton_2& m_ss;
  std::unordered_map<HDS_Halfedge_const_handle, Point_2>& m_offset_points;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  const FT m_vertical_weight;
  std::map<Point_2, Point_2>& m_snapped_positions;
#endif
};

using Offset_builder = CGAL::Polygon_offset_builder_2<Straight_skeleton_2,
                                                      Offset_builder_traits,
                                                      Polygon_2,
                                                      Skeleton_offset_correspondence_builder_visitor>;

using HDS = Straight_skeleton_2::Base;

FT default_offset = std::numeric_limits<FT>::infinity();

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

bool read_dat_polygon(const char* filename,
                      Polygon_with_holes_2& p)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Could not read " << filename << std::endl;
    return false;
  }

  bool is_number_of_CC_in_input = false;
  if(CGAL::IO::internal::get_file_extension(filename) == "poly")
  {
    is_number_of_CC_in_input = true;
  }

  std::vector<Polygon_2> polys;

  auto read_polygon = [&in, &polys](int i) -> void
  {
    std::vector<Point_2> poly;

    int v_count = 0;
    in >> v_count;
    for(int j=0; j<v_count && in; ++j)
    {
      double x = 0., y = 0.;
      in >> x >> y;
      poly.push_back(Point_2(x, y));
    }

    if(poly.size() >= 3)
    {
      bool is_simple = CGAL::is_simple_2(poly.begin(), poly.end(), K());
      if(!is_simple)
        std::cerr << "Warning: input polygon not simple (hopefully it is strictly simple...)" << std::endl;

      CGAL::Orientation expected = (i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);

      const double area = CGAL::to_double(CGAL::polygon_area_2(poly.begin(), poly.end(), K()));
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

      if(orientation == expected)
        polys.emplace_back(poly.begin(), poly.end());
      else
        polys.emplace_back(poly.rbegin(), poly.rend());
    }
  };

  if(is_number_of_CC_in_input)
  {
    int ccb_count = 0;
    in >> ccb_count;
    for(int i=0; i<ccb_count && in; ++i)
      read_polygon(i);
  }
  else
  {
    int i = 0;
    while(in)
      read_polygon(i++);
  }

  if(polys.empty())
  {
    std::cerr << "Error: empty input?" << std::endl;
    return false;
  }

  std::cout <<"Polygon with border of size: " << polys[0].size() << std::endl;
  if(polys.size() > 1)
    std::cout << polys.size() - 1 << " hole(s)" << std::endl;

  p = Polygon_with_holes_2(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[i+1]);

  return true;
}

bool read_input_polygon(const char* filename,
                        Polygon_with_holes_2& p)
{
  std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "dat")
  {
    return read_dat_polygon(filename, p);
  }
  else
  {
    std::cerr << "Error: unknown file extension: " << ext << std::endl;
    return false;
  }
}

bool read_segment_speeds(const char* filename,
                        std::vector<std::vector<FT> >& weights)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Count not read " << filename << std::endl;
    return false;
  }

  std::vector<FT> border_weights;

  std::string line;
  while(getline(in, line))
  {
    if(line.empty())
    {
      weights.push_back(border_weights);
      border_weights.clear();
    }

    std::istringstream iss(line);
    double w;
    if(iss >> w)
      border_weights.push_back(w);
  }

  // in case the last line is not empty
  if(!border_weights.empty())
    weights.push_back(border_weights);

  return true;
}

Polygon_with_holes_2 generate_square_polygon()
{
  Polygon_2 poly;
  poly.push_back(Point_2(0,0));
  poly.push_back(Point_2(1,0));
  poly.push_back(Point_2(1,1));
  poly.push_back(Point_2(0,1));

  CGAL_assertion(poly.is_counterclockwise_oriented());
  return Polygon_with_holes_2{poly};
}

void generate_random_weights(const Polygon_with_holes_2& p,
                             const double min_weight,
                             const double max_weight,
                             const unsigned int seed,
                             std::vector<std::vector<FT> >& speeds)
{
  CGAL::Random rnd(seed);
  std::cout << "Seed is " << rnd.get_seed() << std::endl;

  CGAL_assertion(max_weight > 1);

  auto prev = [](const auto& it, const auto& container)
  {
    return it == container.begin() ? std::prev(container.end()) : std::prev(it);
  };

  auto next = [](const auto& it, const auto& container)
  {
    return (it == std::prev(container.end())) ? container.begin() : std::next(it);
  };

  auto generate_range_weights = [prev, next, min_weight, max_weight, &rnd](const auto& c)
  {
    using Container = typename std::remove_reference<decltype(c)>::type;
    using Iterator = typename Container::const_iterator;

    std::map<Iterator, std::size_t /*rnd weight*/> weight;

    // start somewhere not collinear
    Iterator start_it;
    for(Iterator it=c.begin(); it<c.end(); ++it)
    {
      // the edge is [prev_1 ; it], check for collinearity with the previous edge [prev_2; prev_1]
      auto prev_2 = prev(prev(it, c), c);
      auto prev_1 = prev(it, c);
      Segment_2 s0 {*prev_2, *prev_1}, s1 {*prev_1, *it};
      if(!SS::are_edges_orderly_collinear(s0, s1))
        start_it = it;
    }

    CGAL_assertion(start_it != Iterator()); // all collinear is impossible

    Iterator it=start_it, end=start_it;
    do
    {
      auto prev_2 = prev(prev(it, c), c);
      auto prev_1 = prev(it, c);
      Segment_2 s0 {*prev_2, *prev_1}, s1 {*prev_1, *it};
      if(SS::are_edges_orderly_collinear(s0, s1))
      {
        CGAL_assertion(weight.count(prev_1) != 0);
        weight[it] = weight[prev_1];
      }
      else
      {
        CGAL_assertion(weight.count(it) == 0);
        weight[it] = rnd.get_double(min_weight, max_weight);
      }

      it = next(it, c);
    }
    while(it != end);

    std::vector<FT> weights;
    for(auto it=c.begin(); it<c.end(); ++it)
      weights.push_back(weight[it]);

    return weights;
  };

  speeds.push_back(generate_range_weights(p.outer_boundary()));
  for(auto hit=p.holes_begin(); hit!=p.holes_end(); ++hit)
    speeds.push_back(generate_range_weights(*hit));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PointRange, typename FaceRange>
void construct_horizontal_faces(const Polygon_with_holes_2& p,
                                const FT altitude,
                                PointRange& points,
                                FaceRange& faces,
                                const bool invert_faces = false)
{
#ifdef CGAL_SLS_DEBUG_DRAW
  CGAL::draw(p);
#endif

  CDT cdt;
  cdt.insert_constraint(p.outer_boundary().begin(), p.outer_boundary().end(), true /*close*/);
  for(auto h_it=p.holes_begin(); h_it!=p.holes_end(); ++h_it)
    cdt.insert_constraint(h_it->begin(), h_it->end(), true /*close*/);

  std::size_t id = points.size(); // point ID offset (previous faces inserted their points)
  for(CDT_Vertex_handle vh : cdt.finite_vertex_handles())
  {
    points.emplace_back(cdt.point(vh).x(), cdt.point(vh).y(), altitude);
    vh->info() = id++;
  }

#ifdef CGAL_SLS_DEBUG_DRAW
  CGAL::draw(cdt);
#endif

  std::unordered_map<CDT_Face_handle, bool> in_domain_map;
  boost::associative_property_map< std::unordered_map<CDT_Face_handle, bool> > in_domain(in_domain_map);

  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  for(CDT_Face_handle f : cdt.finite_face_handles())
  {
    if(!get(in_domain, f))
      continue;

    // invert faces for the z=0 plane (bottom face)
    if(invert_faces)
      faces.push_back({f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info()});
    else
      faces.push_back({f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info()});
  }
}

template <typename PointRange, typename FaceRange>
void construct_horizontal_faces(const Offset_polygons_with_holes& p_ptrs,
                                const FT altitude,
                                PointRange& points,
                                FaceRange& faces)
{
  for(const auto& p_ptr : p_ptrs)
    construct_horizontal_faces(*p_ptr, altitude, points, faces);
}

template <typename SLSFacePoints, typename PointRange, typename FaceRange>
void triangulate_skeleton_face(SLSFacePoints& face_points,
                               const bool invert_faces,
                               PointRange& points,
                               FaceRange& faces)
{
  CGAL_precondition(face_points.size() >= 3);

  // shift once to ensure that face_points[0] and face_points[1] are at z=0 and thus the normal is correct
  std::rotate(face_points.rbegin(), face_points.rbegin() + 1, face_points.rend());
  CGAL_assertion(face_points[0][2] == 0 && face_points[1][2] == 0);

  const Vector_3 n = CGAL::cross_product(face_points[1] - face_points[0], face_points[2] - face_points[0]);
  PK traits(n);
  PCDT pcdt(traits);
  pcdt.insert_constraint(face_points.begin(), face_points.end(), true /*close*/);

  std::size_t id = points.size(); // point ID offset (previous faces inserted their points);
  for(PCDT_Vertex_handle vh : pcdt.finite_vertex_handles())
  {
    points.push_back(pcdt.point(vh));
    vh->info() = id++;
  }

#ifdef CGAL_SLS_DEBUG_DRAW
  CGAL::draw(pcdt);
#endif

  std::unordered_map<PCDT_Face_handle, bool> in_domain_map;
  boost::associative_property_map< std::unordered_map<PCDT_Face_handle, bool> > in_domain(in_domain_map);

  CGAL::mark_domain_in_triangulation(pcdt, in_domain);

  for(PCDT_Face_handle f : pcdt.finite_face_handles())
  {
    if(!get(in_domain, f))
      continue;

    // invert faces for exterior skeletons
    if(invert_faces)
      faces.push_back({f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info()});
    else
      faces.push_back({f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info()});
  }
}

// This version is for default offset, so just gather all the full faces
// @todo this doesn't not support holes in SLS faces
void construct_lateral_faces(const Straight_skeleton_2& ss,
                             std::vector<Point_3>& points,
                             std::vector<std::vector<std::size_t> >& faces,
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                             const FT& vertical_weight,
                             std::map<Point_2, Point_2>& snapped_positions,
#endif
                             const bool ignore_frame_faces = false,
                             const bool invert_faces = false)
{
  const HDS& hds = const_cast<const HDS&>(static_cast<const HDS&>(ss));

  std::size_t fc = 0;

  for(const HDS_Face_handle hds_f : CGAL::faces(hds))
  {
    std::vector<Point_3> face_points;

    // If they exist (exterior skeleton), the first four faces of the SLS correspond
    // to the outer frame, and should be ignored.
    if(ignore_frame_faces && fc++ < 4)
      continue;

    HDS_Halfedge_const_handle hds_h = hds_f->halfedge(), done = hds_h;
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    HDS_Halfedge_const_handle contour_h = hds_h->defining_contour_edge();
    CGAL_assertion(hds_h == contour_h);
    const bool is_vertical = (contour_h->weight() == vertical_weight);
#endif

    do
    {
      HDS_Vertex_const_handle hds_tv = hds_h->vertex();

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      // this computes the snapped position but does not change the geometry of the skeleton
      if(is_vertical && !hds_tv->is_contour())
        snap_skeleton_vertex(hds_h, contour_h, snapped_positions);
#endif

      face_points.emplace_back(hds_tv->point().x(), hds_tv->point().y(), hds_tv->time());

      hds_h = hds_h->next();
    }
    while(hds_h != done);

    if(face_points.size() < 3)
    {
      std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
      continue;
    }

    triangulate_skeleton_face(face_points, invert_faces, points, faces);
  }
}

// @todo this doesn't not support holes in SLS faces
void construct_lateral_faces(const Straight_skeleton_2& ss,
                             const Offset_builder& offset_builder,
                             const FT offset,
                             std::vector<Point_3>& points,
                             std::vector<std::vector<std::size_t> >& faces,
                             const std::unordered_map<HDS_Halfedge_const_handle, Point_2>& offset_points,
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                             const FT& vertical_weight,
                             std::map<Point_2, Point_2>& snapped_positions,
#endif
                             const bool ignore_frame_faces = false,
                             const bool invert_faces = false)
{
  CGAL_precondition(offset != default_offset);

  const HDS& hds = const_cast<const HDS&>(static_cast<const HDS&>(ss));

  std::size_t fc = 0;

  for(const HDS_Face_handle hds_f : CGAL::faces(hds))
  {
    // If they exist (exterior skeleton), the first four faces of the SLS correspond
    // to the outer frame, and should be ignored.
    if(ignore_frame_faces && fc++ < 4)
      continue;

    std::vector<Point_3> face_points;

    HDS_Halfedge_const_handle hds_h = hds_f->halfedge(), done = hds_h;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    HDS_Halfedge_const_handle contour_h = hds_h->defining_contour_edge();
    CGAL_assertion(hds_h == contour_h);
    const bool is_vertical = (contour_h->weight() == vertical_weight);
#endif

    do
    {
      HDS_Vertex_const_handle hds_sv = hds_h->opposite()->vertex();
      HDS_Vertex_const_handle hds_tv = hds_h->vertex();

      // Compare_offset_against_event_time compares offset to node->time(),
      // so when the offset is greater or equal than the node->time(), the node is the face
      auto compare_time_to_offset = [&](HDS_Vertex_const_handle node) -> CGAL::Comparison_result
      {
        if(node->is_contour())
          return CGAL::LARGER; // offset > 0 and contour nodes' time is 0
        else
          return offset_builder.Compare_offset_against_event_time(offset, node);
      };

      const CGAL::Comparison_result sc = compare_time_to_offset(hds_sv);
      const CGAL::Comparison_result tc = compare_time_to_offset(hds_tv);

      // if the offset is crossing at the source, it will be added when seen as a target
      // from the previous halfedge

      if(sc != tc && sc != CGAL::EQUAL && tc != CGAL::EQUAL)
      {
        // std::cout << "sc != tc" << std::endl;
        CGAL_assertion(sc != CGAL::EQUAL && tc != CGAL::EQUAL);

        HDS_Halfedge_const_handle hds_off_h = hds_h;
        if(hds_h->slope() == CGAL::NEGATIVE) // ensure same geometric point on both sides
          hds_off_h = hds_off_h->opposite();

        // The offset point must already been computed in the offset builder visitor
        auto off_p = offset_points.find(hds_off_h);
        CGAL_assertion(off_p != offset_points.end());

        face_points.emplace_back(off_p->second.x(), off_p->second.y(), offset);
      }

      if(tc != CGAL::SMALLER)
      {
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
        if(is_vertical && !hds_tv->is_contour())
          snap_skeleton_vertex(hds_h, contour_h, snapped_positions);
#endif

        const Point_2& off_p = hds_tv->point();
        face_points.emplace_back(off_p.x(), off_p.y(), hds_tv->time());
      }

      hds_h = hds_h->next();
    }
    while(hds_h != done);

    if(face_points.size() < 3)
    {
      std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
      continue;
    }

    triangulate_skeleton_face(face_points, invert_faces, points, faces);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

enum class Slope
{
  UNKNOWN = 0,
  INWARD,
  OUTWARD,
  VERTICAL
};

// handle vertical angles (inf speed)
// returns
// - whether the weights are positive or negative (inward / outward)
// - whether the input weights are valid
// - the weight of vertical slabs
template <typename WeightRange>
std::tuple<Slope, bool, FT> preprocess_weights(WeightRange& weights)
{
  CGAL_precondition(!weights.empty());

  Slope slope = Slope::UNKNOWN;

  FT max_value = 0; // non-inf, maximum absolute value
  for(auto& contour_weights : weights)
  {
    for(FT& w : contour_weights)
    {
      if(w == 0)
      {
        std::cerr << "Error: null weight (null angle) is not a valid input" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }

      // '0' means a vertical slab, aka 90° angle (see preprocess_angles())
      if(w == 0)
        continue;

      // determine whether weights indicate all inward or all outward
      if(slope == Slope::UNKNOWN)
      {
        // w is neither 0 nor inf here
        slope = (w > 0) ? Slope::INWARD : Slope::OUTWARD;
      }
      else if(slope == Slope::INWARD && w < 0)
      {
        std::cerr << "Error: mixing positive and negative weights is not yet supported" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }
      else if(slope == Slope::OUTWARD && w > 0)
      {
        std::cerr << "Error: mixing positive and negative weights is not yet supported" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }

      // if we are going outwards, it is just an interior skeleton with opposite weights
      w = CGAL::abs(w);
      if(w > max_value)
        max_value = w;
    }
  }

  if(slope == Slope::UNKNOWN)
  {
    std::cerr << "Warning: all edges vertical?" << std::endl;
    slope = Slope::VERTICAL;
  }

  // Take a weight which is a large % of the max value to ensure there's no ambiguity
  //
  // Since the max value might not be very close to 90°, take the max between of the large-% weight
  // and the weight corresponding to an angle of 89.9999999°
  const FT weight_of_89d9999999 = 572957787.3425436; // tan(89.9999999°)
  const FT scaled_max = (std::max)(weight_of_89d9999999, 1e3 * max_value);

  for(auto& contour_weights : weights)
  {
    for(FT& w : contour_weights)
    {
      if(w == FT(0))
        w = scaled_max;
    }
  }

  return {slope, true, scaled_max};
}

// convert angles (in degrees) to weights, and handle vertical angles
template <typename AngleRange>
std::tuple<Slope, bool, FT> preprocess_angles(AngleRange& angles)
{
  CGAL_precondition(!angles.empty());

  auto angle_to_weight = [](const FT angle) -> FT
  {
    CGAL_precondition(0 < angle && angle < 180);

    // @todo should this be an epsilon around 90°? As theta goes to 90°, tan(theta) goes to infinity
    // and thus we could get numerical issues (overlfows) if the kernel is not exact
    if(angle == 90)
      return 0;
    else
      return std::tan(CGAL::to_double(angle * CGAL_PI / 180));
  };

  for(auto& contour_angles : angles)
    for(FT& angle : contour_angles)
      angle = angle_to_weight(angle);

  return preprocess_weights(angles);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// this is roughly "CGAL::create_interior_weighted_skeleton_and_offset_polygons_with_holes_2()",
// but we want to know the intermediate straight skeleton to build the lateral faces of the 3D mesh
template <typename PointRange, typename FaceRange>
bool inward_construction(const Polygon_with_holes_2& pwh,
                         const std::vector<std::vector<FT> >& speeds,
                         const FT vertical_weight,
                         const FT offset,
                         PointRange& points,
                         FaceRange& faces)
{
  // Avoid recomputing offset points multiple times
  std::unordered_map<HDS_Halfedge_const_handle, Point_2> offset_points;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  // This is to deal with vertical slabs: we temporarily give a non-vertical slope to be able
  // to construct the SLS, and we then snap back the position to the vertical planes.
  // Note that points in non-vertical faces are also changed a bit
  std::map<Point_2, Point_2> snapped_positions;
#endif

  Straight_skeleton_2_ptr ss_ptr;

  if(offset == default_offset)
  {
    ss_ptr = CGAL::create_interior_weighted_straight_skeleton_2(
                      SS::vertices_begin(pwh.outer_boundary()),
                      SS::vertices_end(pwh.outer_boundary()),
                      pwh.holes_begin(), pwh.holes_end(),
                      speeds,
                      K());
  }
  else
  {
    ss_ptr = SS::create_partial_interior_weighted_straight_skeleton_2(
                    offset,
                    SS::vertices_begin(pwh.outer_boundary()),
                    SS::vertices_end(pwh.outer_boundary()),
                    pwh.holes_begin(), pwh.holes_end(),
                    speeds,
                    K());
  }

  if(!ss_ptr)
  {
    std::cerr << "Error: encountered an error during skeleton construction" << std::endl;
    return false;
  }

#ifdef CGAL_SLS_DEBUG_DRAW
  // print_straight_skeleton(*ss_ptr);
  CGAL::draw(*ss_ptr);
#endif

  if(offset == default_offset)
  {
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    construct_lateral_faces(*ss_ptr, points, faces, vertical_weight, snapped_positions);
#else
    construct_lateral_faces(*ss_ptr, points, faces);
#endif
  }
  else // offset != default_offset
  {
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    Skeleton_offset_correspondence_builder_visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
#else
    Skeleton_offset_correspondence_builder_visitor visitor(*ss_ptr, vertical_weight, offset_points);
#endif
    Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
    Offset_polygons raw_output;
    ob.construct_offset_contours(offset, std::back_inserter(raw_output));

    Offset_polygons_with_holes output = CGAL::arrange_offset_polygons_2<Polygon_with_holes_2>(raw_output);
    construct_horizontal_faces(output, offset, points, faces);

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    construct_lateral_faces(*ss_ptr, ob, offset, points, faces, offset_points, vertical_weight, snapped_positions);
#else
    construct_lateral_faces(*ss_ptr, ob, offset, points, faces, offset_points);
#endif
  }

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  apply_snapping(points, snapped_positions);
#endif

  return true;
}

template <typename PointRange, typename FaceRange>
bool outward_construction(const Polygon_with_holes_2& pwh,
                          const std::vector<std::vector<FT> >& speeds,
                          const FT vertical_weight,
                          const FT offset,
                          PointRange& points,
                          FaceRange& faces)
{
  CGAL_precondition(offset != default_offset); // was checked before, this is just a reminder

  // Avoid recomputing offset points multiple times
  std::unordered_map<HDS_Halfedge_const_handle, Point_2> offset_points;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  // This is to deal with vertical slabs: we temporarily give a non-vertical slope to be able
  // to construct the SLS, and we then snap back the position to the vertical planes.
  // Note that points in non-vertical faces are also changed a bit
  std::map<Point_2, Point_2> snapped_positions;
#endif

  Offset_polygons raw_output; // accumulates for both the outer boundary and the holes

  // the exterior of a polygon with holes is the exterior of its outer boundary,
  // and the interior of its inverted holes
  //
  // Start with the outer boundary
  {
    std::vector<std::vector<FT> > outer_speeds = { speeds[0] };
    Straight_skeleton_2_ptr ss_ptr = SS::create_partial_exterior_weighted_straight_skeleton_2(
                                        offset,
                                        SS::vertices_begin(pwh.outer_boundary()),
                                        SS::vertices_end(pwh.outer_boundary()),
                                        outer_speeds,
                                        K());

    if(!ss_ptr)
    {
      std::cerr << "Error: encountered an error during outer skeleton construction" << std::endl;
      return false;
    }

#ifdef CGAL_SLS_DEBUG_DRAW
    // print_straight_skeleton(*ss_ptr);
    CGAL::draw(*ss_ptr);
#endif

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    Skeleton_offset_correspondence_builder_visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
#else
    Skeleton_offset_correspondence_builder_visitor visitor(*ss_ptr, offset_points);
#endif
    Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
    ob.construct_offset_contours(offset, std::back_inserter(raw_output));

    // Manually filter the offset of the outer frame
    std::swap(raw_output[0], raw_output.back());
    raw_output.pop_back();

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    construct_lateral_faces(*ss_ptr, ob, offset, points, faces, offset_points,
                            vertical_weight, snapped_positions,
                            true /*ignore frame faces*/, true /*invert faces*/);
#else
    construct_lateral_faces(*ss_ptr, ob, offset, points, faces, offset_points,
                            true /*ignore frame faces*/, true /*invert faces*/);
#endif
  }

  // now, deal with the holes

  std::size_t hole_id = 1;
  for(auto hit=pwh.holes_begin(); hit!=pwh.holes_end(); ++hit, ++hole_id)
  {
    Polygon_2 hole = *hit; // intentional copy
    hole.reverse_orientation();

    // this is roughly "CGAL::create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2()",
    // but we want to know the intermediate straight skeleton to build the lateral faces of the 3D mesh

    std::vector<Polygon_2> no_holes;
    std::vector<std::vector<FT> > hole_speeds = { speeds[hole_id] };
    Straight_skeleton_2_ptr ss_ptr = SS::create_partial_interior_weighted_straight_skeleton_2(
                                        offset,
                                        SS::vertices_begin(hole), SS::vertices_end(hole),
                                        no_holes.begin(), no_holes.end(),
                                        hole_speeds,
                                        K());

    if(!ss_ptr)
    {
      std::cerr << "Error: encountered an error during skeleton construction" << std::endl;
      return EXIT_FAILURE;
    }

#ifdef CGAL_SLS_DEBUG_DRAW
    // print_straight_skeleton(*ss_ptr);
    CGAL::draw(*ss_ptr);
#endif

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    Skeleton_offset_correspondence_builder_visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
#else
    Skeleton_offset_correspondence_builder_visitor visitor(*ss_ptr, offset_points);
#endif
    Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
    ob.construct_offset_contours(offset, std::back_inserter(raw_output));

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    construct_lateral_faces(*ss_ptr, ob, offset, points, faces, offset_points,
                            vertical_weight, snapped_positions,
                            false /*no outer frame*/, true /*invert faces*/);
#else
    construct_lateral_faces(*ss_ptr, ob, offset, points, faces, offset_points,
                            false /*no outer frame*/, true /*invert faces*/);
#endif
  }

  // - the exterior offset of the outer boundary is built by creating an extra frame and turning
  // the outer boundary into a hole. Hence, it needs to be reversed back to proper orientation
  // - the exterior offset of the holes is built by reversing the holes and computing an internal
  // skeleton. Hence, the result also needs to be reversed.
  for(boost::shared_ptr<Polygon_2> ptr : raw_output)
    ptr->reverse_orientation();

  Offset_polygons_with_holes output = CGAL::arrange_offset_polygons_2<Polygon_with_holes_2>(raw_output);
  construct_horizontal_faces(output, offset, points, faces);

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  apply_snapping(points, snapped_positions);
#endif

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(!std::numeric_limits<FT>::has_infinity)
  {
    std::cerr << "Warning: number type has no infinity, using double's max()" << std::endl;
    default_offset = FT{std::numeric_limits<double>::max()};
  }

  const int argc_check = argc - 1;

  char* poly_filename = nullptr;
  char* speeds_filename = nullptr;

  FT offset = default_offset;
  bool use_angles = false; // whether the input is SLS edge weights, or taper angles
  bool flip_weights = false; // takes the opposite for weights, and the complement for angles

  // below is only used for random weight generation
  double min_weight = 1., max_weight = 10.;
  std::size_t seed = std::time(nullptr);

  for(int i = 1; i < argc; ++i)
  {
    if(!strcmp("-h", argv[i]) || !strcmp("--help", argv[i]) || !strcmp("-?", argv[i]))
    {
      std::cout << "Usage: " << argv[0] << "[options].\n"
        "Options:\n"
        "   -i <input_filename>: input polygon filename.\n"
        "   -t <value>: time (== height). Must be strictly positive.\n"
        "   -a <angles_filename>: angles. Format: one angle per line, a space to separate borders.\n"
        "   -w <weights_filename>: weights. Format: one weight per line, a space to separate borders.\n"
        " Note: -w and -a are exclusive.\n"
                << std::endl;

      return EXIT_FAILURE;
    } else if(!strcmp("-i", argv[i]) && i < argc_check) {
      poly_filename = argv[++i];
    } else if(!strcmp("-w", argv[i])) {
      if(speeds_filename != nullptr)
      {
        std::cerr << "Error: -w and -a are exclusive." << std::endl;
        return EXIT_FAILURE;
      }
      speeds_filename = argv[++i];
    } else if(!strcmp("-a", argv[i])) {
      if(speeds_filename != nullptr)
      {
        std::cerr << "Error: -w and -a are exclusive." << std::endl;
        return EXIT_FAILURE;
      }
      speeds_filename = argv[++i];
      use_angles = true;
    } else if(!strcmp("-t", argv[i]) && i < argc_check) {
      offset = std::stod(argv[++i]);
    } else if(!strcmp("-mw", argv[i]) && i < argc_check) {
      min_weight = std::stod(argv[++i]);
    } else if(!strcmp("-Mw", argv[i]) && i < argc_check) {
      max_weight = std::stod(argv[++i]);
    } else if(!strcmp("-f", argv[i]) && i < argc) {
      flip_weights = true;
    } else if(!strcmp("-s", argv[i]) && i < argc_check) {
      seed = std::stoi(argv[++i]);
    }
  }

  if(offset <= FT(0))
  {
    std::cerr << "Error: height/offset/time must be strictly positive; it is " << offset << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  Polygon_with_holes_2 pwh;
  if(poly_filename == nullptr)
  {
    pwh = generate_square_polygon();
  }
  else if(!read_input_polygon(poly_filename, pwh) || pwh.outer_boundary().is_empty())
  {
    std::cerr << "Error: failure during polygon read" << std::endl;
    return EXIT_FAILURE;
  }

#ifdef CGAL_SLS_OUTPUT_FILES
  std::ofstream out_poly("input.dat");
  out_poly.precision(17);
  out_poly << pwh;
  out_poly.close();
#endif

  // read segment speeds (angles or weights)
  std::vector<std::vector<FT> > speeds;
  if(speeds_filename == nullptr)
    generate_random_weights(pwh, min_weight, max_weight, seed, speeds);
  else
    read_segment_speeds(speeds_filename, speeds);

  if(flip_weights)
  {
    if(use_angles)
    {
      for(auto& contour_speeds : speeds)
        for(FT& a : contour_speeds)
          a = 180 - a;
    }
    else
    {
      for(auto& contour_speeds : speeds)
        for(FT& w : contour_speeds)
          w = -w;
    }
  }

  timer.stop();
  std::cout << "Reading input(s) took " << timer.time() << " s." << std::endl;

  // End of I/O, do some slope preprocessing and check the validity of the input(s)
  // -----------------------------------------------------------------------------------------------

  timer.reset();
  timer.start();

  Slope slope;
  bool valid_input;
  FT vertical_weight;

  if(use_angles)
    std::tie(slope, valid_input, vertical_weight) = preprocess_angles(speeds);
  else
    std::tie(slope, valid_input, vertical_weight) = preprocess_weights(speeds);

  if(!valid_input)
  {
    std::cerr << "Error: invalid input weights" << std::endl;
    return EXIT_FAILURE;
  }

  switch(slope)
  {
    case Slope::UNKNOWN: std::cout << "Slope is UNKNOWN??" << std::endl; break;
    case Slope::INWARD: std::cout << "Slope is INWARD" << std::endl; break;
    case Slope::OUTWARD: std::cout << "Slope is OUTWARD" << std::endl; break;
    case Slope::VERTICAL: std::cout << "Slope is VERTICAL" << std::endl; break;
  }

  if(slope != Slope::INWARD && offset == default_offset)
  {
    std::cerr << "Error: offset must be specified when using an outward (or vertical) slope" << std::endl;
    return EXIT_FAILURE;
  }

  // End of preprocessing, start the actual skeleton computation
  // -----------------------------------------------------------------------------------------------

  // build a soup, to be converted to a mesh afterwards
  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > faces;

  // just a reasonnable guess
  points.reserve(2 * pwh.outer_boundary().size());
  faces.reserve(2 * pwh.outer_boundary().size() + 2*pwh.number_of_holes());

  // bottom face (z=0)
  construct_horizontal_faces(pwh, 0 /*altitude*/, points, faces, true /*invert faces*/);

  bool res;
  if(slope != Slope::OUTWARD) // INWARD or VERTICAL
    res = inward_construction(pwh, speeds, vertical_weight, offset, points, faces);
  else
    res = outward_construction(pwh, speeds, vertical_weight, offset, points, faces);

  if(!res)
    return EXIT_FAILURE;

#ifdef CGAL_SLS_OUTPUT_FILES
  // This soup provides one connected component per edge of the input polygon
  CGAL::IO::write_polygon_soup("extruded_skeleton_soup.off", points, faces, CGAL::parameters::stream_precision(17));
#endif

  // Convert the triangle soup to a triangle mesh

  PMP::merge_duplicate_points_in_polygon_soup(points, faces);
  if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    PMP::orient_polygon_soup(points, faces);
  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));

  Mesh sm;
  PMP::polygon_soup_to_polygon_mesh(points, faces, sm);

  timer.stop();
  std::cout << "Offset computation took " << timer.time() << " s." << std::endl;

  CGAL::IO::write_polygon_mesh("extruded_skeleton.off", sm, CGAL::parameters::stream_precision(17));

  CGAL_assertion(is_valid_polygon_mesh(sm) && is_closed(sm));
  CGAL_warning(!PMP::does_self_intersect(sm)); // cannot be otherwise if there is non-manifoldness

  return EXIT_SUCCESS;
}