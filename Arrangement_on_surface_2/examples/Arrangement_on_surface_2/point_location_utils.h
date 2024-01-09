#include <cassert>
#include <CGAL/Arr_point_location_result.h>

//-----------------------------------------------------------------------------
// Print the result of a point-location query.
//
template <typename Arrangement_>
void print_point_location
(const typename Arrangement_::Point_2& q,
 typename CGAL::Arr_point_location_result<Arrangement_>::Type obj)
{
  using Arrangement_2 = Arrangement_;
  using Vertex_const_handle = typename Arrangement_2::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement_2::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement_2::Face_const_handle;

  const Vertex_const_handle* v;
  const Halfedge_const_handle* e;
  const Face_const_handle* f;

  std::cout << "The point (" << q << ") is located ";
  if ((f = std::get_if<Face_const_handle>(&obj)))          // inside a face
    std::cout << "inside "
              << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
              << " face." << std::endl;
  else if ((e = std::get_if<Halfedge_const_handle>(&obj))) // on an edge
    std::cout << "on an edge: " << (*e)->curve() << std::endl;
  else if ((v = std::get_if<Vertex_const_handle>(&obj)))   // on a vertex
    std::cout << "on " << (((*v)->is_isolated()) ? "an isolated" : "a")
              << " vertex: " << (*v)->point() << std::endl;
  else CGAL_error_msg("Invalid object.");
}

//-----------------------------------------------------------------------------
// Perform a point-location query and print the result.
//
template <typename PointLocation>
void locate_point(const PointLocation& pl,
                  const typename PointLocation::Arrangement_2::Point_2& q)
{
  // Perform the point-location query.
  using Point_location = PointLocation;
  using Arrangement_2 = typename Point_location::Arrangement_2;
  typename CGAL::Arr_point_location_result<Arrangement_2>::Type obj =
    pl.locate(q);

  // Print the result.
  print_point_location<Arrangement_2>(q, obj);
}

//-----------------------------------------------------------------------------
// Perform a vertical ray-shooting query and print the result.
//
template <typename VerticalRayShooting>
void shoot_vertical_ray(const VerticalRayShooting& vrs,
                        const typename
                        VerticalRayShooting::Arrangement_2::Point_2& q)
{
  using Vertical_ray_shooting = VerticalRayShooting;

  // Perform the point-location query.
  typename Vertical_ray_shooting::result_type obj = vrs.ray_shoot_up(q);

  // Print the result.
  using Arrangement_2 = typename Vertical_ray_shooting::Arrangement_2;
  using Vertex_const_handle = typename Arrangement_2::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement_2::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement_2::Face_const_handle;

  const Vertex_const_handle* v;
  const Halfedge_const_handle* e;
  const Face_const_handle* f;

  std::cout << "Shooting up from (" << q << ") : hit ";

  if ((v = std::get_if<Vertex_const_handle>(&obj)))         // hit a vertex
    std::cout << (((*v)->is_isolated()) ? "an isolated" : "a")
              << " vertex: " << (*v)->point() << std::endl;
  else if ((e = std::get_if<Halfedge_const_handle>(&obj)) ) // hit an edge
    std::cout << "an edge: " << (*e)->curve() << std::endl;
  else if ((f = std::get_if<Face_const_handle>(&obj))) {    // hit nothing
    assert((*f)->is_unbounded());
    std::cout << "nothing." << std::endl;
  }
  else CGAL_error_msg("Invalid object.");
}

//-----------------------------------------------------------------------------
// Construct the arrangement of segments needed for the point-location and
// the vertical ray-shooting examples.
// The function assumes that the arrangement is of line segments with integer
// coordinates.
//
template <typename Arrangement_>
void construct_segments_arr(Arrangement_& arr)
{
  using Arrangement_2 = Arrangement_;
  using Point_2 = typename Arrangement_2::Point_2;
  using Segment_2 = typename Arrangement_2::X_monotone_curve_2;
  using Halfedge_handle = typename Arrangement_2::Halfedge_handle;

  Point_2 p0(3,2), p1(0,3), p2(2,5), p3(4,5), p4(6,3), p5(3,0);
  Segment_2 s1(p1, p2), s2(p2, p3), s3(p3, p4), s4(p4, p5), s5(p5, p1);

  arr.insert_in_face_interior(p0, arr.unbounded_face());

  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Halfedge_handle e2 = arr.insert_from_left_vertex(s2, e1->target());
  Halfedge_handle e3 = arr.insert_from_left_vertex(s3, e2->target());
  Halfedge_handle e4 = arr.insert_from_right_vertex(s4, e3->target());
  arr.insert_at_vertices(s5, e4->target(), e1->source());
}
