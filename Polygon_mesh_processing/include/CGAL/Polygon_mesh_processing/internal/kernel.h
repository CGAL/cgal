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

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/clip_convex.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Homogeneous.h>
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
  using Plane_3 = std::array<typename Kernel::Point_3, 3>;
  using Point_3 = typename Kernel::Point_3;
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

  struct Construct_point_on_plane_3{
    Point_3 operator()(const Plane_3& plane)
    {
      return plane[0];
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

  Construct_point_on_plane_3 construct_point_on_plane_3_object() const
  {
    return Construct_point_on_plane_3();
  }

  using Compute_scalar_product_3 = typename Kernel::Compute_scalar_product_3;
  Compute_scalar_product_3 compute_scalar_product_3_object() const { return Compute_scalar_product_3(); }

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


template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
TriangleMesh
kernel(const TriangleMesh& pm,
       const NamedParameters& np = parameters::default_values())
{
  // TODO: bench if with EPECK we can directly use Kernel::Plane_3
  // TODO: TriangleMesh as output is not correct since it is actually a PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using Point_3 = typename GT::Point_3;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  //TODO: what do we do with a mesh that is not closed?
  //TODO: what do we do if the input is not a triangle mesh?

  if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
    return TriangleMesh();



  CGAL::Bbox_3 bb3 = bbox(pm, np);
  TriangleMesh kernel;
  CGAL::make_hexahedron(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                        Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()),
                        kernel);

  Three_point_cut_plane_traits<GT> kgt;
#ifdef CGAL_USE_CONCAVE_FACE_FIRST
  auto compute_dihedral_angle=[](const TriangleMesh& mesh, const halfedge_descriptor h) {
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
  for (auto f : faces(pm))
  {
    auto h = halfedge(f, pm);
    auto plane = make_array( get(vpm,source(h, pm)),
                             get(vpm,target(h, pm)),
                             get(vpm,target(next(h, pm), pm)) );
#endif

#ifdef CGAL_USE_OPTI_WITH_BBOX
    auto pred = kgt.oriented_side_3_object();
    bb3 = bbox(kernel);
    std::array<Point_3, 8> corners = CGAL::make_array(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()),
                                                      Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()),
                                                      Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()),
                                                      Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                                                      Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()),
                                                      Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()),
                                                      Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()),
                                                      Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()));
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
        return TriangleMesh();
      }
    }
#endif
    clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(kgt).do_not_triangulate_faces(true).used_for_kernel(true));
    // clip(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(kgt).do_not_triangulate_faces(true).used_for_kernel(true));
    if (is_empty(kernel)) break;
  }

  return kernel;
};

template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
TriangleMesh
kernel_using_chull(const TriangleMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
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

  TriangleMesh kernel;

  if (origin_opt.has_value())
    CGAL::halfspace_intersection_3(planes.begin(),
                                   planes.end(),
                                   kernel,
                                   origin_opt);

  return kernel;
}

template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
TriangleMesh
kernel_using_chull_and_constructions(const TriangleMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
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

  TriangleMesh kernel;

  if (origin_opt.has_value())
    CGAL::halfspace_intersection_with_constructions_3(planes.begin(),
                                                      planes.end(),
                                                      kernel,
                                                      origin_opt);

  return kernel;
}



//__________________________________________________________________________________________________________________________________

using int256 = boost::multiprecision::int256_t;
using int512 = boost::multiprecision::int512_t;
using intexact = Exact_integer;

template <class Kernel>
struct Plane_based_traits
{
  using FT = typename Kernel::FT;
  using RT = typename Kernel::RT;
  using Geometric_point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

  using plane_descriptor = std::size_t;
  using Plane_3 = std::pair<typename Kernel::Plane_3, plane_descriptor>;

private:
  using plane_range_pointer = std::shared_ptr<std::vector<Plane_3>>;
  plane_range_pointer m_planes;
public:



  struct Point_3 : public Geometric_point_3{
    std::array<plane_descriptor, 3> supports;
    mutable std::set<plane_descriptor> other_coplanar;
    using Base = Geometric_point_3;

    Point_3(){}
    Point_3(plane_descriptor a, plane_descriptor b, plane_descriptor c, const std::vector<Plane_3> &planes):
          Base(*CGAL::Intersections::internal::intersection_point(planes[a].first, planes[b].first, planes[c].first, Kernel())),
          supports({a,b,c}){
    }

    Point_3(plane_descriptor a, plane_descriptor b, plane_descriptor c, Geometric_point_3 p): Base(p), supports({a,b,c}){}
  };

  struct Construct_point_3{
    plane_range_pointer m_planes;
    Construct_point_3(plane_range_pointer planes) : m_planes(planes){}
    Point_3 operator()(plane_descriptor a, plane_descriptor b, plane_descriptor c){
      const std::vector<Plane_3> &planes = *m_planes;
      auto res = CGAL::Intersections::internal::intersection_point(planes[a].first, planes[b].first, planes[c].first, Kernel());
      CGAL_assertion(res);
      return Point_3(a, b, c, *res);
    }

    Point_3 operator()(plane_descriptor a, plane_descriptor b, plane_descriptor c, Geometric_point_3 p){
      return Point_3(a, b, c, p);
    }
  };

  struct Does_not_support_CDT2{};

  struct Oriented_side_3
  {
    Oriented_side operator()(const Plane_3& plane, const Point_3& p)  const
    {
      if((p.supports[0]==plane.second) || (p.supports[1]==plane.second) || (p.supports[2]==plane.second))
        return COPLANAR;
      Oriented_side ori = plane.first.oriented_side(p);
      if(ori==COPLANAR)
        p.other_coplanar.emplace(plane.second);
      return ori;
    }
  };

  struct Construct_plane_line_intersection_point_3
  {
    plane_range_pointer m_planes;
    Construct_plane_line_intersection_point_3(plane_range_pointer planes) : m_planes(planes){}
    Point_3 operator()(const Plane_3& plane, const Point_3& p, const Point_3& q)
    {
      Construct_point_3 point_3(m_planes);
      const std::vector<Plane_3> &planes=*m_planes;

      auto get_common_supports=[&](const Point_3& p, const Point_3 &q){
        plane_descriptor first=-1;
        plane_descriptor second=-1;
        for(short int i=0; i!=3; ++i)
          for(short int j=0; j!=3; ++j)
            if(p.supports[i]==q.supports[j])
              if(first==-1){
                first=p.supports[i];
                break;
              } else {
                second=p.supports[i];
                return std::make_pair(first, second);
              }

        for(plane_descriptor pd: p.other_coplanar)
          for(short int j=0; j!=3; ++j)
            if(pd==q.supports[j])
              if(first==-1){
                first=pd;
                break;
              } else if(planes[first].first!=planes[pd].first){
                second=pd;
                return std::make_pair(first, second);
              }

        for(plane_descriptor pd: q.other_coplanar)
          for(short int j=0; j!=3; ++j)
            if(pd==p.supports[j])
              if(first==-1){
                first=pd;
                break;
              } else if(planes[first].first!=planes[pd].first){
                second=pd;
                return std::make_pair(first, second);
              }

        for(plane_descriptor pd: p.other_coplanar)
          for(plane_descriptor qd: q.other_coplanar)
            if(pd==qd)
              if(first==-1){
                first=pd;
                break;
              } else if(planes[first].first!=planes[pd].first){
                second=pd;
                return std::make_pair(first, second);
              }
        // The two points do not shair a common support
        std::cout << p.supports[0] << " " << p.supports[1] << " " << p.supports[2] << std::endl;
        std::cout << q.supports[0] << " " << q.supports[1] << " " << q.supports[2] << std::endl;
        CGAL_assertion_code(std::cout << "The two points do not share two common supporting planes" << std::endl;)
        CGAL_assertion(0);
        return std::make_pair(first, second);
      };

      std::pair<plane_descriptor, plane_descriptor> line_supports=get_common_supports(p, q);
      assert((planes[line_supports.first].first != planes[line_supports.second].first) &&
             (planes[line_supports.first].first != planes[line_supports.second].first.opposite()));
      Point_3 res=point_3(plane.second, line_supports.first, line_supports.second);

      namespace mp = boost::multiprecision;
      CGAL_assertion_code(int256 max2E195=mp::pow(int256(2),195);)
      CGAL_assertion_code(int256 max2E169=mp::pow(int256(2),169);)
      // CGAL_assertion((mp::abs(p.hx())<=max2E195) && (mp::abs(p.hy())<=max2E195) && (mp::abs(p.hz())<=max2E195) && (mp::abs(p.hw())<=max2E169));
      return res;
    }
  };

  template<class Mesh>
  Plane_based_traits(const Mesh &m){
    namespace mp = boost::multiprecision;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;

    auto to_int=[](const typename Mesh::Point &p){
      return Geometric_point_3(int(p.x()),int(p.y()),int(p.z()));
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

  Oriented_side_3 oriented_side_3_object() const
  {
    return Oriented_side_3();
  }

  Construct_plane_line_intersection_point_3 construct_plane_line_intersection_point_3_object()
  {
    return Construct_plane_line_intersection_point_3(m_planes);
  }

  struct Construct_plane_3{

    plane_range_pointer m_planes;
    Construct_plane_3(plane_range_pointer planes) : m_planes(planes){}
    Plane_3 operator()(const typename Kernel::Plane_3 &pl){
      namespace mp = boost::multiprecision;
      m_planes->emplace_back(pl, m_planes->size());
      CGAL_assertion_code(int256 max2E55=mp::pow(int256(2),55);)
      CGAL_assertion_code(int256 max2E82=mp::pow(int256(2),82);)
      // CGAL_assertion((mp::abs(pl.a())<=max2E55) && (mp::abs(pl.b())<=max2E55) && (mp::abs(pl.c())<=max2E82) && (mp::abs(pl.d())<=max2E82));
      return m_planes->back();
    }

    Plane_3 operator()(const RT &a, const RT &b, const RT &c, const RT &d){
      return (*this)(typename Kernel::Plane_3(a,b,c,d));
    }
  };

  Construct_plane_3 construct_plane_3_object()
  {
    return Construct_plane_3(m_planes);
  }

  Construct_point_3 construct_point_3_object()
  {
    return Construct_point_3(m_planes);
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

template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
Surface_mesh<Plane_based_traits<Homogeneous<int256>>::Point_3>
trettner_kernel(const TriangleMesh& pm,
                const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = Plane_based_traits<Homogeneous<int256>>;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));

  using Point_3 = typename GT::Point_3;
  using Plane_3 = typename GT::Plane_3;

  using Construct_plane_3 = GT::Construct_plane_3;
  using Construct_point_3 = GT::Construct_point_3;

  using InternMesh = Surface_mesh<Point_3>;

  GT gt(pm);
  const std::vector<Plane_3> &planes=gt.planes();

  Construct_plane_3 plane_3 = gt.construct_plane_3_object();
  Construct_point_3 point_3 = gt.construct_point_3_object();

  if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
    return Surface_mesh<Point_3>();

  CGAL::Bbox_3 bb3 = bbox(pm, np);
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

  for (auto plane: planes)
  {
    // clip(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(gt).do_not_triangulate_faces(true).used_for_kernel(true));
    if (is_empty(kernel)) break;
  }

  return kernel;
}

template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
Surface_mesh<Plane_based_traits<Exact_predicates_exact_constructions_kernel>::Point_3>
trettner_epeck_kernel(const TriangleMesh& pm,
                const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // using GT = Plane_based_traits<Homogeneous<int256>>;
  using GT = Plane_based_traits<Exact_predicates_exact_constructions_kernel>;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));

  using Point_3 = typename GT::Point_3;
  using Plane_3 = typename GT::Plane_3;

  using Construct_plane_3 = GT::Construct_plane_3;
  using Construct_point_3 = GT::Construct_point_3;

  using InternMesh = Surface_mesh<Point_3>;

  GT gt(pm);
  const std::vector<Plane_3> &planes=gt.planes();

  Construct_plane_3 plane_3 = gt.construct_plane_3_object();
  Construct_point_3 point_3 = gt.construct_point_3_object();

  if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
    return Surface_mesh<Point_3>();

  CGAL::Bbox_3 bb3 = bbox(pm, np);
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

  for (auto plane: planes)
  {
    // clip(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(gt).do_not_triangulate_faces(true).used_for_kernel(true));
    if (is_empty(kernel)) break;
  }

  return kernel;
}


} } } // end of CGAL::Polygon_mesh_processing::experimental

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H