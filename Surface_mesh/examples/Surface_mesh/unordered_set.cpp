#include <iostream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

template <typename T>
struct Unordered
{
  Unordered()
  {
    std::hash<T> fct;
    T t;
    fct(t);
  }
};


typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;

typedef CGAL::Surface_mesh<Point>       Surface_mesh;
typedef CGAL::Triangulation_2<Kernel>   Triangulation_2;
typedef CGAL::Polyhedron_3<Kernel>      Polyhedron_3;

int main()
{
  {
    std::cout << "Surface_mesh::.._index"<< std::endl;
    Unordered<Surface_mesh::Vertex_index> U1;
    Unordered<Surface_mesh::Halfedge_index> U2;
    Unordered<Surface_mesh::Edge_index> U3;
    Unordered<Surface_mesh::Face_index> U4;

    std::cout << "\ngraph_traits<Surface_mesh>::_descriptor"<< std::endl;
    Unordered<boost::graph_traits<Surface_mesh>::vertex_descriptor> U5;
    Unordered<boost::graph_traits<Surface_mesh>::halfedge_descriptor> U6;
    Unordered<boost::graph_traits<Surface_mesh>::edge_descriptor> U7;
    Unordered<boost::graph_traits<Surface_mesh>::face_descriptor> U8;
  }
  {
    std::cout << "\nTriangulation_2::.._handle"<< std::endl;
    Unordered<Triangulation_2::Vertex_handle> U1;
    Unordered<Triangulation_2::Face_handle> U2;

    std::cout << "\ngraph_traits<Triangulation_2>::_descriptor"<< std::endl;
    Unordered<boost::graph_traits<Triangulation_2>::vertex_descriptor> U3;
    Unordered<boost::graph_traits<Triangulation_2>::halfedge_descriptor> U4;
    Unordered<boost::graph_traits<Triangulation_2>::edge_descriptor> U5;
    Unordered<boost::graph_traits<Triangulation_2>::face_descriptor> U6;
  }
  {
    std::cout << "\nPolyhedron_3::.._handle"<< std::endl;
    Unordered<Polyhedron_3::Vertex_handle> U1;
    Unordered<Polyhedron_3::Halfedge_handle> U2;
    Unordered<Polyhedron_3::Face_handle> U3;

    std::cout << "\ngraph_traits<Polyhedron_3>::_descriptor"<< std::endl;
    Unordered<boost::graph_traits<Polyhedron_3>::vertex_descriptor> U5;
    Unordered<boost::graph_traits<Polyhedron_3>::halfedge_descriptor> U6;
    Unordered<boost::graph_traits<Polyhedron_3>::edge_descriptor> U7;
    Unordered<boost::graph_traits<Polyhedron_3>::face_descriptor> U8;
  }
  return 0;
}
