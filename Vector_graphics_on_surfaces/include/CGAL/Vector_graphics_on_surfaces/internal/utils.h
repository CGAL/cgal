// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_INTERNAL_UTILS_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_INTERNAL_UTILS_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

namespace CGAL {
namespace Vector_graphics_on_surfaces {
namespace internal {

template <class GT, class TriangleMesh, class VertexPointMap>
std::tuple<bool, std::array<typename GT::FT,3>>
point_in_triangle(const VertexPointMap &vpm,
                  const TriangleMesh &mesh,
                  const typename boost::graph_traits<TriangleMesh>::face_descriptor& face,
                  const typename GT::Point_3 point,
                  typename GT::FT tol=1e-5)
{
  // http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  // pag.48
  std::array<typename GT::FT,3> b = make_array(0.,0.,0.);
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h=halfedge(face,mesh);
  typename GT::Point_3 v0 = get(vpm,source(h,mesh));
  typename GT::Point_3 v1 = get(vpm,target(h,mesh));
  typename GT::Point_3 v2 = get(vpm,target(next(h,mesh),mesh));

  typename GT::Vector_3 u = v1 - v0, v = v2 - v0, w = point - v0;
  typename GT::FT d00 = u.squared_length(), d01 = u*v, d11 = v.squared_length(), d20 = w*u,
                 d21 = w*v, d = d00 * d11 - d01 * d01;

  if (d == 0)
    return {false, make_array(0.,0.,0.)};

  b[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(b[2]));
  b[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(b[1]));
  b[0] = 1 - b[1] - b[2];
  assert(!isnan(b[0]));

  for (auto i = 0; i < 3; ++i) {
    if (b[i] < -tol || b[i] > 1.0 + tol)
      return {false, make_array(0.,0.,0.)};
  }

  return {true, b};
}

// TODO: recode using CGAL code?
template <class GT>
typename GT::Vector_2
intersect_circles(const typename GT::Vector_2 &c2, typename GT::FT R2,
                  const typename GT::Vector_2 &c1, typename GT::FT R1)
{
  auto R = (c2 - c1).squared_length();
  assert(R > 0);
  auto invR = typename GT::FT(1) / R;
  typename GT::Vector_2 result = c1+c2;

  result = result + (c2 - c1) * ((R1 - R2) * invR);
  auto A = 2 * (R1 + R2) * invR;
  auto B = (R1 - R2) * invR;
  auto s = A - B * B - 1;
  assert(s >= 0);
  result = result + typename GT::Vector_2(c2.y() - c1.y(), c1.x() - c2.x()) * sqrt(s);
  return result / 2.;
}

template <class GT, class VertexPointMap, class TriangleMesh>
std::array<typename GT::Vector_2, 3>
init_flat_triangle( const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor& h,
                    const VertexPointMap &vpm, const TriangleMesh &mesh)
{
  using Vector_2 = typename GT::Vector_2;
  auto triangle_vertices = make_array(source(h, mesh), target(h, mesh), target(next(h, mesh), mesh));

  std::array<typename GT::Vector_2, 3> tr2d;
  tr2d[0] = Vector_2(0, 0);
  tr2d[1] =
      Vector_2(0, sqrt(squared_distance(get(vpm, triangle_vertices[0]),
                                        get(vpm, triangle_vertices[1]))));
  auto rx = squared_distance(get(vpm, triangle_vertices[0]),
                             get(vpm, triangle_vertices[2]));
  auto ry = squared_distance(get(vpm, triangle_vertices[1]),
                             get(vpm, triangle_vertices[2]));
  tr2d[2] = intersect_circles<GT>(tr2d[0], rx, tr2d[1], ry);

  return tr2d;
}

//case k=0 assume that flat_tid has been flattened putting x-axis aligned with h_edge
template <class GT, class VertexPointMap, class TriangleMesh>
std::array<typename GT::Vector_2, 3>
unfold_face(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor& h_edge,
            const VertexPointMap &vpm, const TriangleMesh &mesh,
            const std::array<typename GT::Vector_2, 3>& flat_tid,const int k=0)
{
  const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_opp = opposite(h_edge, mesh);

  const typename boost::graph_traits<TriangleMesh>::vertex_descriptor v = target(next(h_opp,mesh),mesh);
  const typename boost::graph_traits<TriangleMesh>::vertex_descriptor a = source(h_edge, mesh);
  const typename boost::graph_traits<TriangleMesh>::vertex_descriptor b = target(h_edge,mesh);
  typename GT::FT r0 = squared_distance(get(vpm,v), get(vpm,a));
  typename GT::FT r1 = squared_distance(get(vpm,v), get(vpm,b));

  typename GT::Vector_2 v2 = intersect_circles<GT>(flat_tid[(k+1)%3], r1, flat_tid[k], r0);

  const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_ref_opp = halfedge(face(h_opp,mesh),mesh);

  std::array<typename GT::Vector_2, 3> res;

  if(h_ref_opp==h_opp)
  {
    res[0]=flat_tid[(k+1)%3];
    res[1]=flat_tid[k];
    res[2]=v2;
  } else if(next(h_ref_opp,mesh)==h_opp)
  {
     res[0]=v2;
     res[1]=flat_tid[(k+1)%3];
     res[2]=flat_tid[k];
  }else
  {
     assert(prev(h_ref_opp,mesh)==h_opp);
     res[0]=flat_tid[k];
     res[1]=v2;
     res[2]=flat_tid[(k+1)%3];
  }

  return res;
}

template <class GT, class VertexPointMap, class TriangleMesh>
std::array<typename GT::Vector_2, 2>
unfold_face(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_curr,
            typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_next,
            const VertexPointMap &vpm, const TriangleMesh &mesh,
            const std::array<typename GT::Vector_2, 2>& flat_tid)
{
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_next_opp = opposite(h_next, mesh);
  CGAL_assertion(face( h_curr, mesh) == face(h_next_opp, mesh));


  typename boost::graph_traits<TriangleMesh>::vertex_descriptor v = target(next(h_curr,mesh),mesh);
  typename boost::graph_traits<TriangleMesh>::vertex_descriptor a = target(h_curr,mesh);
  typename boost::graph_traits<TriangleMesh>::vertex_descriptor b = source(h_curr, mesh);
  typename GT::FT r0 = squared_distance(get(vpm,v), get(vpm,a));
  typename GT::FT r1 = squared_distance(get(vpm,v), get(vpm,b));

  typename GT::Vector_2 v2 = intersect_circles<GT>(flat_tid[1], r1, flat_tid[0], r0);


  std::array<typename GT::Vector_2, 2> res;
  if(next(h_curr, mesh) == h_next_opp)
  {
     res[0]=flat_tid[0];
     //res[2]=flat_tid[1];
     res[1]=v2;
#ifdef CGAL_DEBUG_BSURF
     std::cout << "4 " << res[0]  << " 0 "
               << res[1]  << " 0 "
               << flat_tid[1] << " 0 "
               << res[0]  << " 0\n";
#endif
  }
  else
  {
    CGAL_assertion(prev(h_curr, mesh) == h_next_opp);
    res[0]=v2;
    res[1]=flat_tid[1];
    //res[2]=flat_tid[0];
#ifdef CGAL_DEBUG_BSURF
    std::cout << "4 " << res[0]  << " 0 "
              << res[1]  << " 0 "
              << flat_tid[0] << " 0 "
              << res[0]  << " 0\n";
#endif
  }

  return res;
}

//TODO get rid of tol or at least get it dependent on the input?
template<class GT, class FaceLocation>
std::tuple<bool,int>
point_is_vert(const FaceLocation& p,const typename GT::FT& tol=1e-5)
{
   auto bary=p.second;
   if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol)
       return {true, 0};
   if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol)
      return {true, 1};
   if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol)
      return {true, 2};

 return {false, -1};
}

//TODO get rid of tol or at least get it dependent on the input?
template<class GT, class FaceLocation>
std::tuple<bool, int>
point_is_edge(const FaceLocation& p,
              const typename GT::FT& tol=1e-5)
{
  auto bary=p.second;
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol)
    return {true, 0};
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol)
    return {true, 1};
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol)
    return {true, 2};

  return {false, -1};
}


//https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
template <class GT>
std::pair<typename GT::FT, typename GT::FT>
intersect(const typename GT::Vector_2 &direction, const typename GT::Vector_2 &left,
          const typename GT::Vector_2 &right,const typename GT::Vector_2& origin)
{
  typename GT::Vector_2 v1 = origin-left;
  typename GT::Vector_2 v2 = right - left;
  typename GT::Vector_2 v3(-direction.y(), direction.x());
  typename GT::FT t0 = (v2.x()*v1.y()-v2.y()*v1.x()) / (v2*v3);
  typename GT::FT t1 = v1*v3/ ( v2*v3 );
  return std::make_pair(t0, t1);
};

template <class GT>
std::tuple<int,typename GT::FT>
segment_in_tri(const typename GT::Vector_2& p,
               const std::array<typename GT::Vector_2, 3>& tri,
               const typename GT::Vector_2& dir,const int offset)
{
  //rotated the triangle in order to test intersection at meaningful edges before
  std::array<typename GT::Vector_2, 3> rotated_tri=tri;
  // TODO rotated_tri.begin() + offset ?
  for(int k=0;k<offset;++k)
    std::rotate(rotated_tri.begin(), rotated_tri.begin() + 1, rotated_tri.end());
#ifdef CGAL_DEBUG_BSURF
  std::cout << "  offset " << offset << "\n";
  std::cout << "  tri " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
  std::cout << "  p " << p << "\n";
  std::cout << "  p+dir " << p+dir << "\n";
  std::cout << "  " << rotated_tri[0] << " " << rotated_tri[1] << " " << rotated_tri[2] << "\n";
#endif
  for (auto i = 0; i < 3; ++i)
  {
    auto [t0, t1] = intersect<GT>(dir, rotated_tri[(i+1)%3], rotated_tri[(i+2) % 3],p);
#ifdef CGAL_DEBUG_BSURF
    std::cout  << "  t0/t1 " << t0 << " / " << t1 << "\n";
#endif
    //TODO: replace intersection with CGAL code
    if (t0 > 0 && t1 >= -1e-4 && t1 <= 1 + 1e-4)
    {
#ifdef CGAL_DEBUG_BSURF
        int h=((i+1)%3 + offset)%3;
        Vector_2 intersection_point=(1-t1)*tri[h]+t1*tri[(h+1)%3];
        std::cout<<"  2D intersection point: "<< intersection_point << std::endl;
#endif
        return {((i+1)%3 + offset)%3, std::clamp(t1, 0., 1.)}; //return the offset w.r.t h_ref
    }
  }
  CGAL_assertion(false);
  return {-1,-1};
}

// return an angle in degrees
template <class GT, class TriangleMesh, class VertexPointMap>
typename GT::FT
get_total_angle(const typename boost::graph_traits<TriangleMesh>::vertex_descriptor& vid,
                const TriangleMesh& mesh,
                const VertexPointMap &vpm)
{
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_ref= halfedge(vid,mesh);
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_start=h_ref;
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_next = next(h_start,mesh);

  typename GT::FT theta=approximate_angle(get(vpm,source(h_start,mesh))-get(vpm,target(h_start,mesh)),
                                          get(vpm,target(h_next,mesh))-get(vpm,target(h_start,mesh)));
  h_start=opposite(h_next,mesh);
  h_next=next(h_start,mesh);
  while(h_start!=h_ref)
  {
    theta+=approximate_angle(get(vpm,source(h_start,mesh))-get(vpm,target(h_start,mesh)),
                 get(vpm,target(h_next,mesh))-get(vpm,target(h_start,mesh)));
    h_start=opposite(h_next,mesh);
    h_next=next(h_start,mesh);
  }

  return theta;
}

template <class GT, class TriangleMesh, class VertexPointMap>
std::tuple<typename GT::Vector_2,typename boost::graph_traits<TriangleMesh>::face_descriptor, typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>
polthier_condition_at_vert(const TriangleMesh& mesh,
                           const VertexPointMap &vpm,
                           const typename boost::graph_traits<TriangleMesh>::vertex_descriptor& vid,
                           const typename boost::graph_traits<TriangleMesh>::face_descriptor& tid,
                           const typename GT::FT& init_angle)
{
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_2 = typename GT::Vector_2;

  auto get_vid_offset=[&mesh](const halfedge_descriptor& h_ref,const vertex_descriptor& vid)
  {
    if(source(h_ref,mesh)==vid)
      return 0;

    if(target(h_ref,mesh)==vid)
      return 1;

    if(target(next(h_ref,mesh),mesh)==vid)
      return 2;

    std::cerr<<"Error! Halfedges are in different faces"<<std::endl;

    CGAL_assertion(false);

    return -1;
  };
  //TODO use interval for robustness + snap if ambiguous
  FT total_angle=get_total_angle<GT>(vid,mesh,vpm);
  FT theta = 0.5 * total_angle; //in degrees
  halfedge_descriptor h=halfedge(tid,mesh);

  while(target(h,mesh)!=vid)
    h=next(h,mesh);


  Point_3 vert = get(vpm,vid);
  Point_3 vert_adj=get(vpm,source(h,mesh));
  FT acc = init_angle;

  FT prev_angle = acc;
#ifdef CGAL_DEBUG_BSURF
  std::cout<<"initial h "<< edge(h,mesh)<<std::endl;
  std::cout<<"acc "<< acc<<std::endl;
  std::cout<<"theta "<< theta<<std::endl;
#endif
  while (acc < theta) {
    h=prev(opposite(h,mesh),mesh);
    prev_angle = acc;
    Point_3 next_vert_adj=get(vpm,source(h,mesh));
    acc += approximate_angle(vert_adj - vert,next_vert_adj - vert);
    vert_adj=next_vert_adj;
  }
  auto offset = theta - prev_angle;
  //moving to radians
  offset*=CGAL_PI/180.;
  Point_3 prev_vert_adj=get(vpm,target(next(h,mesh),mesh));

  FT l = sqrt(squared_distance(prev_vert_adj,vert));
  FT phi = approximate_angle(vert - prev_vert_adj, vert_adj - prev_vert_adj);
  phi*=CGAL_PI/180.;
  FT x = l * std::sin(offset) / std::sin(CGAL_PI - phi - offset);
  FT alpha = x / sqrt(squared_distance(vert_adj, prev_vert_adj));
  halfedge_descriptor prev_h=prev(h,mesh);

  std::array<Vector_2,3> flat_tid = init_flat_triangle<GT>(halfedge(face(h,mesh),mesh),vpm,mesh);
  int kv=get_vid_offset(halfedge(face(h,mesh),mesh),vid);


  Vector_2 q = (1 - alpha) * flat_tid[(kv+1)%3] + alpha * flat_tid[(kv+2)%3];
  Vector_2 new_dir = q - flat_tid[kv];
#ifdef CGAL_DEBUG_BSURF
  Vector_2 flat_p=flat_tid[kv];
  std::cout<<"final h "<< edge(h,mesh)<<std::endl;
  std::cout<<" prev_vert_adj "<< prev_vert_adj<<std::endl;
  std::cout<<" vert_adj "<< vert_adj<<std::endl;
  std::cout<< "------ outgoing tid"<<std::endl;
  std::cout << "  " << flat_tid[0] << " " << flat_tid[1] << " " << flat_tid[2] << "\n";

  std::cout<< "------ flat_p"<<std::endl;
  std::cout << "  " << flat_p<< "\n";
  std::cout<< "------ outgoing dir"<<std::endl;
  std::cout << "  " << new_dir[0]+flat_p[0] << " " << new_dir[1]+flat_p[1] << "\n";
  std::cout<<"outgoing face "<<face(prev_h,mesh)<<std::endl;
  std::cout<< "h_curr "<< edge(h,mesh)<<std::endl;
#endif
  return {new_dir,face(prev_h,mesh), h};
}

//h_ref is the reference halfedge of the face we are in, h_edge is the halfedge along which we want to unfold
template <class GT, class TriangleMesh, class VertexPointMap>
typename GT::Vector_2
compute_new_dir(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_ref,
                const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_edge,
                const std::array<typename GT::FT,3>& prev_coords,
                const std::array<typename GT::FT,3>& curr_coords,
                const VertexPointMap& vpm,
                const TriangleMesh& mesh)
{
  using Vector_2 = typename GT::Vector_2;
  int k=0;

  if(h_edge==prev(h_ref,mesh))
    k=2;
  else if(h_edge==next(h_ref,mesh))
    k=1;
  else
    assert(h_edge==h_ref);

  std::array<Vector_2,3> flat_curr = init_flat_triangle<GT>(h_ref,vpm,mesh);
  std::array<Vector_2,3> flat_prev = unfold_face<GT>(h_edge,vpm,mesh,flat_curr,k);

  Vector_2 prev_flat_p=prev_coords[0]*flat_prev[0]+prev_coords[1]*flat_prev[1]+prev_coords[2]*flat_prev[2];
  Vector_2 curr_flat_p=curr_coords[0]*flat_curr[0]+curr_coords[1]*flat_curr[1]+curr_coords[2]*flat_curr[2];

#ifdef CGAL_DEBUG_BSURF
  std::cout<<"   k is "<<k<<"\n";
  std::cout<<"   Flat_prev: ";
  std::cout << "  " << flat_prev[0] << " " << flat_prev[1] << " " << flat_prev[2] << "\n";
  std::cout<<"   Flat_curr: ";
  std::cout << "  " << flat_curr[0] << " " << flat_curr[1] << " " << flat_curr[2] << "\n";
  //Unfold face take care of the orientation so flat_t1[1] - flat_t1[0] is our old
  //reference frame
  std::cout<<"   h_ref_curr"<<std::endl;
  std::cout << "  " << source(h_ref,mesh) << " " << target(h_ref,mesh) << " " << target(next(h_ref,mesh),mesh) << "\n";
  halfedge_descriptor h_ref_prev=halfedge(face(opposite(h_edge,mesh),mesh),mesh);
  std::cout<<"   h_ref_prev"<<std::endl;
  std::cout << "  " << source(h_ref_prev,mesh) << " " << target(h_ref_prev,mesh) << " " << target(next(h_ref_prev,mesh),mesh) << "\n";
  std::cout<<"     prev_Flat_p "<<prev_flat_p.x()<<" "<< prev_flat_p.y()<<" )\n";
  std::cout<<"     curr_Flat_p "<<curr_flat_p.x()<<" "<< curr_flat_p.y()<<" )\n";
#endif

  return curr_flat_p-prev_flat_p;

}

} } } // end of CGAL::Vector_graphics_on_surfaces::internal namespace

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_INTERNAL_UTILS_H
