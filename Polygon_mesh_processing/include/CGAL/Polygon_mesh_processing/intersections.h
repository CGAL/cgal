#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTIONS_H
#include <CGAL/license/Polygon_mesh_processing/predicate.h>


#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Kernel/global_functions_3.h>

#include <vector>
#include <exception>
#include <boost/foreach.hpp>
#include <boost/range.hpp>

#include <boost/function_output_iterator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL{
namespace internal {

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class VertexPointMap>
struct Intersect_facets
{
  // wrapper to check whether anything is inserted to output iterator
  struct Output_iterator_with_bool
  {
    Output_iterator_with_bool(OutputIterator* out, bool* intersected)
      : m_iterator(out), m_intersected(intersected) { }

    template<class T>
    void operator()(const T& t) {
      *m_intersected = true;
      *(*m_iterator)++ = t;
    }

    OutputIterator* m_iterator;
    bool* m_intersected;
  };


  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;

  // members
  const TM& m_tmesh1;
  const VertexPointMap m_vpmap1;
  const TM& m_tmesh2;
  const VertexPointMap m_vpmap2;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_facets(const TM& mesh1, const TM& mesh2,
                   OutputIterator it,
                   VertexPointMap vpmap1, VertexPointMap vpmap2,
                   const Kernel& kernel)
    :
      m_tmesh1(mesh1),
      m_vpmap1(vpmap1),
      m_tmesh2(mesh2),
      m_vpmap2(vpmap2),
      m_iterator(it),
      m_intersected(false),
      m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tmesh1);
    halfedge_descriptor g = halfedge(c->info(), m_tmesh2);


    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap1, target(h,m_tmesh1)),
                                    get(m_vpmap1, target(next(h,m_tmesh1),m_tmesh1)),
                                    get(m_vpmap1, target(next(next(h,m_tmesh1),m_tmesh1),m_tmesh1)));

    Triangle t2 = triangle_functor( get(m_vpmap2, target(g,m_tmesh2)),
                                    get(m_vpmap2, target(next(g,m_tmesh2),m_tmesh2)),
                                    get(m_vpmap2, target(next(next(g,m_tmesh2),m_tmesh2),m_tmesh2)));
    if(do_intersect_3_functor(t1, t2)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_facets

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class VertexPointMap>
struct Intersect_facet_polyline
{
  // wrapper to check whether anything is inserted to output iterator
  struct Output_iterator_with_bool
  {
    Output_iterator_with_bool(OutputIterator* out, bool* intersected)
      : m_iterator(out), m_intersected(intersected) { }

    template<class T>
    void operator()(const T& t) {
      *m_intersected = true;
      *(*m_iterator)++ = t;
    }

    OutputIterator* m_iterator;
    bool* m_intersected;
  };


  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point;


  // members
  const TM& m_tmesh;
  const std::vector<face_descriptor>& faces;
  const VertexPointMap m_vpmap;
  const std::vector<Point>& polyline;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_facet_polyline(const TM& mesh,
                           const std::vector<face_descriptor>& faces,
                           const std::vector<Point>& polyline,
                           OutputIterator it,
                           VertexPointMap vpmap,
                           const Kernel& kernel)
    :
      m_tmesh(mesh),
      faces(faces),
      polyline(polyline),
      m_vpmap(vpmap),
      m_iterator(it),
      m_intersected(false),
      m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
      segment_functor(kernel.construct_segment_3_object()),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(faces[b->info()], m_tmesh);


    // check for geometric intersection
    Triangle t = triangle_functor( get(m_vpmap, target(h,m_tmesh)),
                                   get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                   get(m_vpmap, target(next(next(h,m_tmesh),m_tmesh),m_tmesh)));

    Segment s = segment_functor(polyline[c->info()], polyline[c->info() + 1]);
    if(do_intersect_3_functor(t, s)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_facet_polyline

template<class Polyline,
         class Kernel,
         class Box,
         class OutputIterator>
struct Intersect_polylines
{
  // wrapper to check whether anything is inserted to output iterator
  struct Output_iterator_with_bool
  {
    Output_iterator_with_bool(OutputIterator* out, bool* intersected)
      : m_iterator(out), m_intersected(intersected) { }

    template<class T>
    void operator()(const T& t) {
      *m_intersected = true;
      *(*m_iterator)++ = t;
    }

    OutputIterator* m_iterator;
    bool* m_intersected;
  };


  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Point_3 Point;


  // members
  const std::vector<Point>& polyline1;
  const std::vector<Point>& polyline2;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_polylines(const Polyline& polyline1,
                      const Polyline& polyline2,
                      OutputIterator it,
                      const Kernel& kernel)
    :
      polyline1(polyline1),
      polyline2(polyline2),
      m_iterator(it),
      m_intersected(false),
      m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
      segment_functor(kernel.construct_segment_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {


    // check for geometric intersection

    Segment s1 = segment_functor(polyline1[b->info()], polyline1[b->info() + 1]);
    Segment s2 = segment_functor(polyline2[c->info()], polyline2[c->info() + 1]);
    if(do_intersect_3_functor(s1, s2)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_polylines

struct Throw_at_output {
  class Throw_at_output_exception: public std::exception
  { };

  template<class T>
  void operator()(const T& /* t */) const {
    throw Throw_at_output_exception();
  }
};
}// namespace internal

namespace Polygon_mesh_processing {

template <class TriangleMesh
          , class FaceRange
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections( const FaceRange& face_range1,
               const FaceRange& face_range2,
               const TriangleMesh& mesh1,
               const TriangleMesh& mesh2,
               OutputIterator out,
               const NamedParameters& np1,
               const NamedParameters& np2)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh1));
  CGAL_precondition(CGAL::is_triangle_mesh(mesh2));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor> Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range1), boost::end(face_range1) )
        );
  boxes2.reserve(
        std::distance( boost::begin(face_range2), boost::end(face_range2) )
        );

  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;

  VertexPointMap vpmap1 = boost::choose_param(get_param(np1, internal_np::vertex_point),
                                              get_const_property_map(boost::vertex_point, mesh1));
  VertexPointMap vpmap2 = boost::choose_param(get_param(np2, internal_np::vertex_point),
                                              get_const_property_map(boost::vertex_point, mesh2));

  BOOST_FOREACH(face_descriptor f, face_range1)
  {
    boxes1.push_back(Box( get(vpmap1, target(halfedge(f,mesh1),mesh1)).bbox()
                          + get(vpmap1, target(next(halfedge(f, mesh1), mesh1), mesh1)).bbox()
                          + get(vpmap1, target(next(next(halfedge(f, mesh1), mesh1), mesh1), mesh1)).bbox(),
                          f));
  }

  BOOST_FOREACH(face_descriptor f, face_range2)
  {
    boxes2.push_back(Box( get(vpmap2, target(halfedge(f,mesh2),mesh2)).bbox()
                          + get(vpmap2, target(next(halfedge(f, mesh2), mesh2), mesh2)).bbox()
                          + get(vpmap2, target(next(next(halfedge(f, mesh2), mesh2), mesh2), mesh2)).bbox(),
                          f));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr;
  std::vector<const Box*> box2_ptr;
  box1_ptr.reserve(num_faces(mesh1));
  box2_ptr.reserve(num_faces(mesh2));

  BOOST_FOREACH(Box& b, boxes1)
      box1_ptr.push_back(&b);
  BOOST_FOREACH(Box& b, boxes2)
      box2_ptr.push_back(&b);

  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;

  CGAL::internal::Intersect_facets<TM,
      GeomTraits,
      Box,
      OutputIterator,
      VertexPointMap>
      intersect_facets(mesh1, mesh2,
                       out,
                       vpmap1, vpmap2,
                       GeomTraits());

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           intersect_facets,cutoff);
  return intersect_facets.m_iterator;
}


template <class TriangleMesh
          , class FaceRange
          , class Polyline
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections( const FaceRange& face_range,
               const Polyline& polyline,
               const TriangleMesh& mesh,
               OutputIterator out,
               const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;

  VertexPointMap vpmap = boost::choose_param(get_param(np, internal_np::vertex_point),
                                             get_const_property_map(boost::vertex_point, mesh));
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  std::vector<face_descriptor> faces;
  faces.reserve(std::distance( boost::begin(face_range), boost::end(face_range) ));

  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t> Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range), boost::end(face_range) )
        );

  boxes2.reserve(
        std::distance( boost::begin(polyline), boost::end(polyline) ) - 1
        );


  BOOST_FOREACH(face_descriptor f, face_range)
  {
    faces.push_back(f);
    boxes1.push_back(Box( get(vpmap, target(halfedge(f,mesh),mesh)).bbox()
                          + get(vpmap, target(next(halfedge(f, mesh), mesh), mesh)).bbox()
                          + get(vpmap, target(next(next(halfedge(f, mesh), mesh), mesh), mesh)).bbox(),
                          faces.size()-1));
  }

  for(std::size_t i =0; i< polyline.size()-1; ++i)
  {
    Point p1 = polyline[i];
    Point p2 = polyline[i+1];
    boxes2.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr;
  std::vector<const Box*> box2_ptr;
  box1_ptr.reserve(boxes1.size());
  box2_ptr.reserve(boxes2.size());

  BOOST_FOREACH(Box& b, boxes1)
      box1_ptr.push_back(&b);
  BOOST_FOREACH(Box& b, boxes2)
      box2_ptr.push_back(&b);

  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;

  CGAL::internal::Intersect_facet_polyline<TM,
      GeomTraits,
      Box,
      OutputIterator,
      VertexPointMap>
      intersect_facet_polyline(mesh,
                               faces,
                               polyline,
                               out,
                               vpmap,
                               GeomTraits());

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           intersect_facet_polyline, cutoff);
  return intersect_facet_polyline.m_iterator;
}


template < class Polyline
           , class OutputIterator
           , class Kernel
           >
OutputIterator
intersections( const Polyline& polyline1,
               const Polyline& polyline2,
               OutputIterator out,
               const Kernel& )
{
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t> Box;
  typedef typename Kernel::Point_3 Point;
  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(polyline1), boost::end(polyline1) ) - 1
        );

  boxes2.reserve(
        std::distance( boost::begin(polyline2), boost::end(polyline2) ) - 1
        );

  for(std::size_t i =0; i< polyline1.size()-1; ++i)
  {
    Point p1 = polyline1[i];
    Point p2 = polyline1[i+1];
    boxes1.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  for(std::size_t i =0; i< polyline2.size()-1; ++i)
  {
    Point p1 = polyline2[i];
    Point p2 = polyline2[i+1];
    boxes2.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr;
  std::vector<const Box*> box2_ptr;
  box1_ptr.reserve(boxes1.size());
  box2_ptr.reserve(boxes2.size());

  BOOST_FOREACH(Box& b, boxes1)
      box1_ptr.push_back(&b);
  BOOST_FOREACH(Box& b, boxes2)
      box2_ptr.push_back(&b);

  // compute intersections filtered out by boxes

  CGAL::internal::Intersect_polylines<Polyline,
      Kernel,
      Box,
      OutputIterator>
      intersect_polylines(polyline1,
                          polyline2,
                          out,
                          Kernel());

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           intersect_polylines, cutoff);
  return intersect_polylines.m_iterator;
}


template <class TriangleMesh
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections(const TriangleMesh& mesh1,
              const TriangleMesh& mesh2,
              OutputIterator out,
              const NamedParameters& np1,
              const NamedParameters& np2)
{
  return intersections(faces(mesh1), faces(mesh2),
                       mesh1, mesh2, out, np1, np2);
}

template <class TriangleMesh
          , class Polyline
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections(const TriangleMesh& mesh,
              const Polyline& polyline,
              OutputIterator out,
              const NamedParameters& np)
{
  return intersections(faces(mesh), polyline, mesh, out, np);
}


template <class TriangleMesh
          , class CGAL_PMP_NP_TEMPLATE_PARAMETERS
          >
bool do_intersect(const TriangleMesh& mesh1
                  , const TriangleMesh& mesh2
                  , const CGAL_PMP_NP_CLASS& np1
                  , const CGAL_PMP_NP_CLASS& np2)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh1));
  CGAL_precondition(CGAL::is_triangle_mesh(mesh2));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    intersections(mesh1,mesh2, OutputIterator(), np1, np2);
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& )
  { return true; }

  return false;
}


template <class TriangleMesh
          , class Polyline
          , class CGAL_PMP_NP_TEMPLATE_PARAMETERS
          >
bool do_intersect(const TriangleMesh& mesh
                  , const Polyline& polyline
                  , const CGAL_PMP_NP_CLASS& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    intersections(mesh,polyline, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& )
  { return true; }

  return false;
}

template < class Polyline
           , class Kernel
           >
bool do_intersect( const Polyline& polyline1
                  , const Polyline& polyline2,
                   const Kernel&)
{
  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    intersections(polyline1,polyline2, OutputIterator(), Kernel());
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& )
  { return true; }

  return false;
}

}//end PMP
} //end CGAL
#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTIONS_H
