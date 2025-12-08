// Copyright (c) 2016-2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H

#include <algorithm>
#include <random>

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/clip_convex.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Exact_integer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace experimental {

template <class Kernel>
struct Three_point_cut_plane_traits
{
  using FT = typename Kernel::FT;
  // using Plane_3 = std::array<typename Kernel::Point_3, 3>;

  using Point_3 = typename Kernel::Point_3;
  struct Plane_3: public std::array<typename Kernel::Point_3, 3>{
    using Base = std::array<typename Kernel::Point_3, 3>;

    Plane_3(const Point_3 &a, const Point_3 &b, const Point_3 &c): Base({a, b, c}){}
    Plane_3(const std::array<Point_3, 3> &arr): Base(arr){}

    bool operator<(const Plane_3 &b) const{
      typename Kernel::Plane_3 pa((*this)[0], (*this)[1], (*this)[2]);
      typename Kernel::Plane_3 pb(b[0], b[1], b[2]);
      Comparison_result res = compare(pa.a(), pb.a());
      if(res == EQUAL)
        res = compare(pa.b(), pb.b());
      if(res == EQUAL)
        res = compare(pa.c(), pb.c());
      if(res == EQUAL)
        res = compare(pa.d(), pb.d());
      return res == SMALLER;
    };

    bool operator==(const Plane_3 &b) const{
      typename Kernel::Plane_3 pa((*this)[0], (*this)[1], (*this)[2]);
      typename Kernel::Plane_3 pb(b[0], b[1], b[2]);
      return pa==pb;
    }
  };
  using Vector_3 = typename Kernel::Vector_3;

  struct Does_not_support_CDT2{};

  struct Oriented_side_3
  {
    Oriented_side operator()(const Plane_3& plane, const Point_3& p)  const
    {
      return orientation(plane[0], plane[1], plane[2], p);
    }
  };

  struct Construct_plane_line_intersection_point_3
  {
    Point_3 operator()(const Plane_3& plane, const Point_3& p, const Point_3& q)
    {
      typename Kernel::Construct_plane_line_intersection_point_3 construction;
      return construction(plane[0], plane[1], plane[2], p, q);
    }
  };

  struct Construct_orthogonal_vector_3{
    Vector_3 operator()(const Plane_3& plane)
    {
      return typename Kernel::Plane_3(plane[0], plane[1], plane[2]).orthogonal_vector();
    }
  };

  struct Compute_squared_distance_3
  {
    using Compute_scalar_product_3 = typename Kernel::Compute_scalar_product_3;
    FT operator()(const Plane_3& plane, const Point_3& p)
    {
      typename Kernel::Plane_3 pl(plane[0], plane[1], plane[2]);
      return Compute_scalar_product_3()(Vector_3(ORIGIN, p), pl.orthogonal_vector())+pl.d();
    }
  };

  Oriented_side_3 oriented_side_3_object() const
  {
    return Oriented_side_3();
  }

  Construct_plane_line_intersection_point_3 construct_plane_line_intersection_point_3_object() const
  {
    return Construct_plane_line_intersection_point_3();
  }

  Construct_orthogonal_vector_3 construct_orthogonal_vector_3_object() const
  {
    return Construct_orthogonal_vector_3();
  }

  Compute_squared_distance_3 compute_squared_distance_3_object() const { return Compute_squared_distance_3(); }

#ifndef CGAL_PLANE_CLIP_DO_NOT_USE_BOX_INTERSECTION_D
// for does self-intersect
  using Segment_3 = typename Kernel::Segment_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Construct_segment_3 = typename Kernel::Construct_segment_3;
  using Construct_triangle_3 =typename  Kernel::Construct_triangle_3;
  using Do_intersect_3 = typename Kernel::Do_intersect_3;
  Construct_segment_3 construct_segment_3_object() const { return Construct_segment_3(); }
  Construct_triangle_3 construct_triangle_3_object() const { return Construct_triangle_3(); }
  Do_intersect_3 do_intersect_3_object() const { return Do_intersect_3(); }
#endif
};

template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel(const PolygonMesh& pm,
       const NamedParameters& np = parameters::default_values())
{
  // TODO: bench if with EPECK we can directly use Kernel::Plane_3
  // TODO: PolygonMesh as output is not correct since it is actually a PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  // using face_descriptor = typename BGT::face_descriptor;
  // using edge_descriptor = typename BGT::edge_descriptor;
  // using halfedge_descriptor = typename BGT::halfedge_descriptor;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using Point_3 = typename GT::Point_3;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Plane_3 = typename Three_point_cut_plane_traits<GT>::Plane_3;

  //TODO: what do we do with a mesh that is not closed?
  //TODO: what do we do if the input is not a triangle mesh?

  if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
    return PolygonMesh();



  CGAL::Bbox_3 bb3 = bbox(pm, np);
  PolygonMesh kernel;
  CGAL::make_hexahedron(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                        Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()),
                        kernel);
#ifdef CGAL_USE_OPTI_WITH_BBOX
  // We store the vertices that realized the bbox
  std::array<vertex_descriptor, 6> bbox_vertices;
  std::array<Point_3, 8> corners = CGAL::make_array(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()),
                                                    Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()),
                                                    Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()),
                                                    Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                                                    Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()),
                                                    Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()),
                                                    Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()),
                                                    Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()));
#endif

  Three_point_cut_plane_traits<GT> kgt;
#ifdef CGAL_USE_CONCAVE_FACE_FIRST
  auto compute_dihedral_angle=[](const PolygonMesh& mesh, const halfedge_descriptor h) {
    // Calcul de l'angle dihédral entre les deux faces à partir des quatre points
    return to_double(approximate_dihedral_angle(mesh.point(mesh.source(h)), mesh.point(mesh.target(h)),
                    mesh.point(mesh.target(mesh.next(h))), mesh.point(mesh.target(mesh.next(mesh.opposite(h))))));
  };

  std::vector<std::pair<face_descriptor, double>> faces_with_angles;
  for(face_descriptor f : faces(pm)){
    double max_dihedral_angle = 0.0;
    auto h = pm.halfedge(f);
    auto hf_circ = pm.halfedges_around_face(h);
    for (auto he : hf_circ) {
        if (!pm.is_border(he)) {
            max_dihedral_angle = (std::max)(max_dihedral_angle,compute_dihedral_angle(pm, he));
        }
        faces_with_angles.push_back({f, max_dihedral_angle});
    }
  }
  std::sort(faces_with_angles.begin(), faces_with_angles.end(),
        [](const std::pair<face_descriptor, double>& a, const std::pair<face_descriptor, double>& b) {
            return a.second > b.second;
        });
  for(auto pair : faces_with_angles)
  {
    auto h = halfedge(pair.first, pm);
    auto plane = make_array( get(vpm,source(h, pm)),
                             get(vpm,target(h, pm)),
                             get(vpm,target(next(h, pm), pm)) );
#else
  std::vector<Plane_3> planes;
  for (auto f : faces(pm)){
    auto h = halfedge(f, pm);
    planes.emplace_back(get(vpm,source(h, pm)),
                        get(vpm,target(h, pm)),
                        get(vpm,target(next(h, pm), pm)));
  }
  // std::sort(planes.begin(), planes.end());
  // std::unique(planes.begin(), planes.end());
  // std::random_shuffle(planes.begin(), planes.end());
  for(auto plane: planes)
  {
  // for (auto f : faces(pm))
  // {
  //   auto h = halfedge(f, pm);
  //   auto plane = make_array( get(vpm,source(h, pm)),
  //                            get(vpm,target(h, pm)),
  //                            get(vpm,target(next(h, pm), pm)) );
#endif

#ifdef CGAL_USE_OPTI_WITH_BBOX
    // TODO looking the sign of a, b, c of the plane, we can look only 2 orientations instead of possibly the eight one

    auto pred = kgt.oriented_side_3_object();
    int i=0;
    auto first_ori=pred(plane, corners[i]);
    while(++i!=8 && first_ori==ON_ORIENTED_BOUNDARY)
      first_ori=pred(plane, corners[i]);

    if (i==8) continue;
    bool all_the_same=true;
    for (;i<8;++i)
    {
      auto other_ori=pred(plane, corners[i]);
      if (other_ori!=ON_ORIENTED_BOUNDARY && other_ori!=first_ori)
      {
        all_the_same=false;
        break;
      }
    }

    if (all_the_same)
    {
      if (first_ori==ON_NEGATIVE_SIDE) continue;
      else
      {
        return PolygonMesh();
      }
    }
#endif
    clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(kgt).do_not_triangulate_faces(true).used_for_kernel(true));
#ifdef CGAL_USE_OPTI_WITH_BBOX
    // Update the bbox
// #define UPDATE_BBOX_AT(k, x, min, first, comp)
//     if(!kernel.is_valid(bbox_vertices[k])){
//       bb3.##x##min() = to_interval(get(vpm, boundaries[k]).##x()).##first;
//       auto it = boundaries.begin();
//       bbox_vertices[k] = *it;
//       for(++it; it!=boundaries.end(); ++it)
//         if(to_interval(get(vpm, *it).##x()).##first ##comp bb3.##x##min()){
//           bb3.##x##min() = to_interval(get(vpm, *it).x()).##first;
//           bbox_vertices[k] = *it;
//         }
//     }
//     UPDATE_BBOX(0,"x","min","first","<");
//     UPDATE_BBOX(2,"y","min","first","<");
    // UPDATE_BBOX(4,y,min,first,<);
    // UPDATE_BBOX(1,x,max,second,>);
    // UPDATE_BBOX(3,y,max,second,>);
    // UPDATE_BBOX(5,y,max,second,>);
#endif
    if (is_empty(kernel)) break;
  }

  return kernel;
};

template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel_using_chull(const PolygonMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<typename GT::Plane_3> planes;
  planes.reserve(faces(pm).size());


  for (auto f : faces(pm))
  {
    auto h = halfedge(f, pm);
    planes.emplace_back( get(vpm,source(h, pm)),
                         get(vpm,target(h, pm)),
                         get(vpm,target(next(h, pm), pm)) );
  }

  auto origin_opt = CGAL::halfspace_intersection_interior_point_3(planes.begin(), planes.end());

  PolygonMesh kernel;

  if (origin_opt.has_value())
    CGAL::halfspace_intersection_3(planes.begin(),
                                   planes.end(),
                                   kernel,
                                   origin_opt);

  return kernel;
}

template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel_using_chull_and_constructions(const PolygonMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<typename GT::Plane_3> planes;
  planes.reserve(faces(pm).size());


  for (auto f : faces(pm))
  {
    auto h = halfedge(f, pm);
    planes.emplace_back( get(vpm,source(h, pm)),
                         get(vpm,target(h, pm)),
                         get(vpm,target(next(h, pm), pm)) );
  }

  auto origin_opt = CGAL::halfspace_intersection_interior_point_3(planes.begin(), planes.end());

  PolygonMesh kernel;

  if (origin_opt.has_value())
    CGAL::halfspace_intersection_with_constructions_3(planes.begin(),
                                                      planes.end(),
                                                      kernel,
                                                      origin_opt);

  return kernel;
}



//__________________________________________________________________________________________________________________________________

template <class Kernel>
struct Plane_based_traits
{
  using FT = typename Kernel::FT;
  using RT = typename Kernel::RT;
  using Geometric_point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

  using plane_descriptor = std::size_t;
  // using Plane_3 = std::pair<typename Kernel::Plane_3, plane_descriptor>;
  struct Plane_3{
    Plane_3(typename Kernel::Plane_3 pl, plane_descriptor i):first(pl), second(i){}
    typename Kernel::Plane_3 first;
    plane_descriptor second;

    bool operator<(const Plane_3 &a) const{
      Comparison_result res = compare(first.a(), a.first.a());
      if(res == EQUAL)
        res = compare(first.b(), a.first.b());
      if(res == EQUAL)
        res = compare(first.c(), a.first.c());
      if(res == EQUAL)
        res = compare(first.d(), a.first.d());
      return res == SMALLER;
    }
    bool operator==(const Plane_3 &a) const{
      return first==a.first;
    }
  };

private:
  using plane_range_pointer = std::shared_ptr<std::vector<Plane_3>>;
  plane_range_pointer m_planes;
public:

  struct Construct_point_3;
  struct Point_3 : public Geometric_point_3{
    // std::array<plane_descriptor, 3> supports;
    // mutable std::set<plane_descriptor> other_coplanar;
    mutable std::vector<plane_descriptor> supports;
    using Base = Geometric_point_3;

    Point_3(){}
    Point_3(plane_descriptor a, plane_descriptor b, plane_descriptor c, Geometric_point_3 p): Base(p){
      supports.reserve(3);
      supports.push_back(a);
      supports.push_back(b);
      supports.push_back(c);
    }
    Point_3(plane_descriptor a, plane_descriptor b, plane_descriptor c, const std::vector<Plane_3> &planes):
          Point_3(a,b,c,Construct_point_3(&planes)(a,b,c)){}
    Point_3(const std::vector<plane_descriptor> &_supports, Geometric_point_3 p): Base(p), supports(_supports){}
    Point_3(const std::vector<plane_descriptor> &_supports, const std::vector<Plane_3> &planes):
          Point_3(_supports,Construct_point_3(&planes)(_supports)){}
  };

  struct Construct_point_3{
    plane_range_pointer m_planes;
    Construct_point_3(plane_range_pointer planes) : m_planes(planes){}
    Point_3 operator()(plane_descriptor a, plane_descriptor b, plane_descriptor c){
      const std::vector<Plane_3> &planes = *m_planes;
      auto res = intersection(planes[a].first, planes[b].first, planes[c].first);
      return Point_3(a, b, c, std::get<Geometric_point_3>(*res));
    }
    Point_3 operator()(const std::vector<plane_descriptor> &supports){
      const std::vector<Plane_3> &planes = *m_planes;
      auto res = intersection(planes[supports[0]].first, planes[supports[1]].first, planes[supports[2]].first);
      return Point_3(supports, std::get<Geometric_point_3>(*res));
    }

    Point_3 operator()(plane_descriptor a, plane_descriptor b, plane_descriptor c, Geometric_point_3 p){
      return Point_3(a, b, c, p);
    }
    Point_3 operator()(const std::vector<plane_descriptor> &supports, Geometric_point_3 p){
      return Point_3(supports, p);
    }
  };

  struct Does_not_support_CDT2{};

  struct Oriented_side_3
  {
    Oriented_side operator()(const Plane_3& plane, const Point_3& p)  const
    {
      for(plane_descriptor i: p.supports)
      if(i==plane.second)
        return COPLANAR;
      Oriented_side ori = plane.first.oriented_side(p);
      if(ori==COPLANAR)
        p.supports.push_back(plane.second);
      return ori;
    }
  };

  struct Compute_squared_distance_3
  {
    using Compute_scalar_product_3 = typename Kernel::Compute_scalar_product_3;
    FT operator()(const Plane_3& plane, const Point_3& p)
    {
      for(plane_descriptor i: p.supports)
      if(i==plane.second)
        return FT(0);
      FT sd = Compute_scalar_product_3()(Vector_3(ORIGIN, p), plane.first.orthogonal_vector())+plane.first.d();
      if(is_zero(sd))
        p.supports.push_back(plane.second);
      return sd;
    }
  };

  struct Construct_plane_line_intersection_point_3
  {
    plane_range_pointer m_planes;
    Construct_plane_line_intersection_point_3(plane_range_pointer planes) : m_planes(planes){}
    Point_3 operator()(const Plane_3& plane, const Point_3& p, const Point_3& q)
    {
      Construct_point_3 point_3(m_planes);
      std::vector<plane_descriptor> line_supports;
      line_supports.reserve(3);
      line_supports.push_back(plane.second);

      auto get_common_supports=[&](const Point_3& p, const Point_3 &q){
        for(plane_descriptor i: p.supports)
          for(plane_descriptor j: q.supports)
            if(i==j)
              line_supports.push_back(i);
        CGAL_assertion(line_supports.size()>=3);
      };
      get_common_supports(p, q);
      CGAL_assertion(((*m_planes)[line_supports[1]].first != (*m_planes)[line_supports[2]].first) &&
                     ((*m_planes)[line_supports[1]].first != (*m_planes)[line_supports[2]].first.opposite()));
      Point_3 res=point_3(line_supports);
      return res;
    }
  };

  template<class Mesh>
  Plane_based_traits(const Mesh &m){
    namespace mp = boost::multiprecision;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;

    auto to_int=[](const typename Mesh::Point &p){
      return Geometric_point_3(FT(p.x()),FT(p.y()),FT(p.z()));
    };
    auto to_int_plane=[&](face_descriptor f){
      auto pmap=boost::make_function_property_map<vertex_descriptor>([&](vertex_descriptor v){
        return to_int(m.point(v));
      });
      auto pl = compute_face_normal(f, m, parameters::vertex_point_map(pmap));
      return typename Kernel::Vector_3(pl.hx(),pl.hy(),pl.hz());
    };

    m_planes = std::make_shared<std::vector<Plane_3>>();
    m_planes->reserve(faces(m).size());

    // The ptr need to be initialized to get the functor
    Construct_plane_3 plane_3 = construct_plane_3_object();

    for(face_descriptor f : faces(m))
      plane_3(typename Kernel::Plane_3(to_int(m.point(m.target(m.halfedge(f)))),
                                       to_int_plane(f)));
  }

  struct Construct_plane_3{

    plane_range_pointer m_planes;
    Construct_plane_3(plane_range_pointer planes) : m_planes(planes){}
    Plane_3 operator()(const typename Kernel::Plane_3 &pl){
      namespace mp = boost::multiprecision;
      m_planes->emplace_back(pl, m_planes->size());
      return m_planes->back();
    }

    Plane_3 operator()(const RT &a, const RT &b, const RT &c, const RT &d){
      return (*this)(typename Kernel::Plane_3(a,b,c,d));
    }
  };

  Oriented_side_3 oriented_side_3_object() const { return Oriented_side_3(); }
  Construct_plane_3 construct_plane_3_object() { return Construct_plane_3(m_planes); }
  Construct_point_3 construct_point_3_object() { return Construct_point_3(m_planes); }
  Compute_squared_distance_3 compute_squared_distance_3_object() const { return Compute_squared_distance_3(); }
  Construct_plane_line_intersection_point_3 construct_plane_line_intersection_point_3_object() const {
    return Construct_plane_line_intersection_point_3(m_planes);
  }

  const std::vector<Plane_3>& planes(){
    return *m_planes;
  }

#ifndef CGAL_PLANE_CLIP_DO_NOT_USE_BOX_INTERSECTION_D
// for does self-intersect
  using Segment_3 = typename Kernel::Segment_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Construct_segment_3 = typename Kernel::Construct_segment_3;
  using Construct_triangle_3 =typename  Kernel::Construct_triangle_3;
  using Do_intersect_3 = typename Kernel::Do_intersect_3;
  Construct_segment_3 construct_segment_3_object() const { return Construct_segment_3(); }
  Construct_triangle_3 construct_triangle_3_object() const { return Construct_triangle_3(); }
  Do_intersect_3 do_intersect_3_object() const { return Do_intersect_3(); }
#endif
};

template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
Surface_mesh<typename Plane_based_traits<typename GetGeomTraits<PolygonMesh, NamedParameters>::type>::Point_3>
plane_based_kernel(const PolygonMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Kernel = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  using GT = Plane_based_traits<Kernel>;

  using Point_3 = typename GT::Point_3;
  using Plane_3 = typename GT::Plane_3;

  using Construct_plane_3 = typename GT::Construct_plane_3;
  using Construct_point_3 = typename GT::Construct_point_3;

  using InternMesh = Surface_mesh<Point_3>;

  GT gt(pm);
  std::vector<Plane_3> planes=gt.planes();
  std::sort(planes.begin(), planes.end());
  std::unique(planes.begin(), planes.end());
  std::shuffle(planes.begin(), planes.end(), std::default_random_engine());

  Construct_plane_3 plane_3 = gt.construct_plane_3_object();
  Construct_point_3 point_3 = gt.construct_point_3_object();

  if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
    return Surface_mesh<Point_3>();

  CGAL::Bbox_3 bb3 = bbox(pm);
  InternMesh kernel;
  Plane_3 xl=plane_3(1,0,0,int(-bb3.xmin()));
  Plane_3 yl=plane_3(0,1,0,int(-bb3.ymin()));
  Plane_3 zl=plane_3(0,0,1,int(-bb3.zmin()));
  Plane_3 xr=plane_3(1,0,0,int(-bb3.xmax()));
  Plane_3 yr=plane_3(0,1,0,int(-bb3.ymax()));
  Plane_3 zr=plane_3(0,0,1,int(-bb3.zmax()));
  CGAL::make_hexahedron(point_3(xl.second, yl.second, zl.second),
                        point_3(xl.second, yl.second, zr.second),
                        point_3(xl.second, yr.second, zr.second),
                        point_3(xl.second, yr.second, zl.second),
                        point_3(xr.second, yr.second, zl.second),
                        point_3(xr.second, yl.second, zl.second),
                        point_3(xr.second, yl.second, zr.second),
                        point_3(xr.second, yr.second, zr.second),
                        kernel);

                int k=0;
  for (auto plane: planes)
  {
    ++k;
    clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(gt).do_not_triangulate_faces(true).used_for_kernel(true));
    // std::ofstream("clip"+std::to_string(k%10)+".off") << kernel;
    if (is_empty(kernel)) break;
  }

  return kernel;
}


} } } // end of CGAL::Polygon_mesh_processing::experimental

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H