/*
  agc course project,
  spherical arrangements of none intersecting arcs of great circles on a sphere

  an example for traversing incoming halfedges of vertices

  author: Kapelushnik Lior
*/

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Spherical_map.h>
#include <CGAL/Sphere_traits.h>

typedef CGAL::Gmpq                                    	Number_type;
typedef CGAL::Cartesian<Number_type>                  	Kernel;
typedef CGAL::Sphere_traits<Kernel>			Traits_2;
typedef Kernel::Direction_3			        Direction_3;
typedef CGAL::SphereTopologicalMap<Kernel>		SphereType;
typedef CGAL::Spherical_map<SphereType, Traits_2>	SphereMap;
typedef Traits_2::Curve_2				Curve_2;

int main() {
  SphereMap sMap; // the spherical arrangement
  std::cout << "inserting curve (1,0,0) - (1,1,1)" << std::endl;
  sMap.insert(Curve_2(Direction_3(1,0,0), Direction_3(1,1,1)));
  std::cout << "inserting curve (1,1,1) - (0,1,1)" << std::endl;
  sMap.insert(Curve_2(Direction_3(1,1,1), Direction_3(0,1,1)));
  std::cout << "inserting curve (1,0,0) - (0,-2,-1)" << std::endl;
  sMap.insert(Curve_2(Direction_3(1,0,0), Direction_3(0,-2,-1)));
  std::cout << "inserting curve (1,0,0) - (1,0,1)" << std::endl;
  sMap.insert(Curve_2(Direction_3(1,0,0), Direction_3(1,0,1)));
  std::cout << "inserting curve (1,0,1) - (1,1,1)" << std::endl;
  sMap.insert(Curve_2(Direction_3(1,0,1), Direction_3(1,1,1)));

  SphereMap::Vertex_iterator vit;
  // loop over vertices and print incident halfedges
  for (vit = sMap.vertices_begin(); vit != sMap.vertices_end(); ++vit) {
    std::cout << "now at vertex with direction: " << vit->direction() <<
      std::endl;

    SphereMap::Halfedge_around_vertex_circulator havc;
    havc = vit->incident_halfedges();
    SphereMap::Halfedge_around_vertex_circulator circEn=havc;
    do {
      std::cout << "incident halfedge: " << havc->curve() <<
	std::endl;
      ++havc;
    } while (havc != circEn);

  }

  std::cout <<sMap;

  return 0;
}
