#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/boost/graph/graph_traits_Arrangement_2.h>

#include <boost/foreach.hpp>

#include <list>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

typedef boost::graph_traits<Arrangement_2>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Arrangement_2>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Arrangement_2>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Arrangement_2>::face_descriptor face_descriptor;
typedef boost::graph_traits<Arrangement_2>::vertex_iterator vertex_iterator;

int main ()
{
  // Construct the arrangement of five intersecting segments.
  Arrangement_2           arr;
  std::list<Segment_2>    segments;

  segments.push_back (Segment_2 (Point_2(0, 0), Point_2(1, 0)));
  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(3, 0)));
  segments.push_back (Segment_2 (Point_2(0, 0), Point_2(1, 2)));  
  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(1, 2)));
  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(2, 1)));
  segments.push_back (Segment_2 (Point_2(2, 1), Point_2(1, 2)));
  segments.push_back (Segment_2 (Point_2(2, 1), Point_2(3, 2)));
  segments.push_back (Segment_2 (Point_2(3, 0), Point_2(3, 2)));
  segments.push_back (Segment_2 (Point_2(1, 2), Point_2(3, 2)));

  insert (arr, segments.begin(), segments.end());

  BOOST_FOREACH(vertex_descriptor vd, vertices(arr)){
    std::cerr << vd->point() << std::endl;
  }

  vertex_iterator vi = *(vertices(arr).first);
  vertex_descriptor vd = *vi;

  halfedge_descriptor hd = halfedge(vd, arr);
  halfedge_descriptor hdo = opposite(hd, arr);

  face_descriptor fd = face(hd,arr);

  hd = halfedge(fd,arr);

  hd = next(hd,arr);
  hd = prev(hd,arr);

  BOOST_FOREACH(vertex_descriptor vd, vertices(arr)){
    std::cerr << vd->point() << std::endl;
    std::size_t d = degree(vd,arr), i=0;
    std::cerr << "degree = " << d << std::endl;
    BOOST_FOREACH(edge_descriptor ed, out_edges(vd,arr)){
      ++i;
      assert(source(ed,arr) == vd);
    }
    assert(i == d);
  }

  int i = 0;
  BOOST_FOREACH(face_descriptor fd, faces(arr)){
    std::cerr << "face degree = "<< degree(fd,arr) << std::endl;
    ++i;
  }
  std::cerr << i << " faces" << std::endl;

  i = 0;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(arr)){
    ++i;
  }
  std::cerr << i << " halfedges" << std::endl;


  BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(*(halfedges(arr).first), arr)){
    std::cerr << target(hd,arr)->point() <<  std::endl;
  }

  i = 0;
  BOOST_FOREACH(edge_descriptor ed, edges(arr)){
    halfedge_descriptor hd = halfedge(ed,arr);
    halfedge_descriptor ohd = opposite(hd,arr);
    assert(ed == edge(ohd,arr));
    ++i;
  }

  std::cerr << i << " edges" << std::endl;


  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << num_vertices(arr)
            << ",  E = " << num_edges(arr) 
            << ",  F = " << num_faces(arr) << std::endl;

  return 0;
}
