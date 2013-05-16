#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
typedef K::Point_3 Point_3;


int main()
{
  std::cout.precision(20);
  Polyhedron P;
  std::cin >> P;

  for(Halfedge_iterator it = P.halfedges_begin(); it != P.halfedges_end(); ++it){
    if(it->is_border()){
      std::cerr << "a border edge" << std::endl;
      Halfedge_around_facet_circulator circ(it), done(circ);
      do{
        std::cout << circ->vertex()->point() << std::endl;
      }while (++circ != done);
        std::cout << circ->vertex()->point() << std::endl;
      return 0;
    }
  } 
  return 0;
}
