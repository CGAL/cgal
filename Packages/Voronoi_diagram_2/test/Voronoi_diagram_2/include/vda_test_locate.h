#ifndef VDA_TEST_LOCATE_H
#define VDA_TEST_LOCATE_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <vector>


template<class VDA, class Projector, class QStream, class OStream>
void test_locate(const VDA& vda, const Projector& project,
		 QStream& isq, OStream& os = std::cout)
{
  std::cout << std::endl;
  std::cout << "is Delaunay graph valid? "
	    << (vda.dual().is_valid() ? "yes" : "no") << std::endl;
  std::cout << std::endl;

  os << "Dimension of Delaunay graph: " << vda.dual().dimension()
     << std::endl << std::endl;

  os << "Vertices of the Delaunay graph:" << std::endl;
  typename VDA::Dual_graph::Finite_vertices_iterator vit;
  for (vit = vda.dual().finite_vertices_begin();
       vit != vda.dual().finite_vertices_end(); ++vit) {
    os << project(vit) << std::endl;
  }
  os << std::endl;

  typedef typename VDA::Voronoi_traits::Point_2  Point_2;
  Point_2 p;
  std::vector<Point_2>  vecp;

  while ( isq >> p ) {
    vecp.push_back(p);
  }

  test_locate_dg(vda, project, vecp, os);
  test_locate_vd(vda, vecp, os);
}

template<class VDA, class Projector, class Point_vector, class OStream>
void test_locate_dg(const VDA& vda, const Projector& project,
		    const Point_vector& vecp, OStream& os)
{
  typedef typename VDA::Voronoi_traits                Voronoi_traits;
  typedef typename Voronoi_traits::Point_locator      Point_locator;
  
  typedef typename Point_locator::Object              Object;
  typedef typename Point_locator::Assign              Assign;
  typedef typename Point_locator::Vertex_handle       Vertex_handle;
  typedef typename Point_locator::Face_handle         Face_handle;
  typedef typename Point_locator::Edge                Edge;

  Vertex_handle v;
  Face_handle   f;
  Edge          e;

  Point_locator locate = vda.voronoi_traits().point_locator_object();
  Assign assign = locate.assign_object();
  Object o;

  os << "Query sites and location feature in dual graph:" << std::endl;
  for (unsigned int i = 0; i < vecp.size(); ++i) {
    os << vecp[i] << "\t --> \t" << std::flush;
    o = locate(vda.dual(), vecp[i]);
    if ( assign(v, o) ) {
      os << "FACE";
    } else if ( assign(e, o) ) {
      os << "EDGE";
    } else if ( assign(f, o) ) {
      os << "VERTEX";
    } else {
      os << " *** NOT READY YET *** ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  if ( vda.dual().dimension() == 1 ) {
    os << "# of vertices: " << vda.dual().number_of_vertices() << std::endl;
    os << "# of edges: " << vda.dual().tds().number_of_edges() << std::endl;
    os << "# of faces: " << vda.dual().number_of_faces() << std::endl;
    os << std::endl << std::endl;

    typename VDA::Dual_graph::All_edges_iterator eit;
    for (eit = vda.dual().all_edges_begin();
	 eit != vda.dual().all_edges_end(); ++eit) {
      typename VDA::Dual_graph::Edge e = *eit;
      typename VDA::Dual_graph::Vertex_handle v1, v2;
      v1 = e.first->vertex((e.second+1)%3);
      v2 = e.first->vertex((e.second+2)%3);
      if ( vda.dual().is_infinite(v1) ) {
	os << "inf";
      } else {
	os << project( v1 );
      }
      os << " - ";
      if ( vda.dual().is_infinite(v2) ) {
	os << "inf";
      } else {
	os << project( v2 );
      }
      os << std::endl;
    }
    os << std::endl;
  }
}

template<class VDA, class Point_vector, class OStream>
void test_locate_vd(const VDA& vda, const Point_vector& vecp, OStream& os)
{
  typedef typename VDA::Voronoi_traits                Voronoi_traits;
  typedef typename Voronoi_traits::Point_locator      Point_locator;
  
  typedef typename Voronoi_traits::Object             Object;
  typedef typename Voronoi_traits::Assign             Assign;
  typedef typename VDA::Vertex_handle                 Vertex_handle;
  typedef typename VDA::Face_handle                   Face_handle;
  typedef typename VDA::Halfedge_handle               Halfedge_handle;

  Vertex_handle      v;
  Face_handle        f;
  Halfedge_handle    e;

  Assign assign = vda.voronoi_traits().assign_object();
  Object o;

  os << "Query sites and location feature in dual graph:" << std::endl;
  for (unsigned int i = 0; i < vecp.size(); ++i) {
    os << vecp[i] << "\t --> \t" << std::flush;
    o = vda.locate(vecp[i]);
    if ( assign(e, o) ) {
      os << "VORONOI EDGE"; 
    } else if ( assign(v, o) ) {
      os << "VORONOI VERTEX";
    } else if ( assign(f, o) ) {
      os << "VORONOI FACE";
    } else {
      os << " *** NOT READY YET *** ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


#endif // VDA_TEST_LOCATE_H
