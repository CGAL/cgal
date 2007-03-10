/*
  agc course project,
  spherical arrangements of none intersecting arcs of great circles on a sphere

  an example for traversing spherical faces boundaries

  author: Kapelushnik Lior
*/

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

#include "CGAL/Spherical_map.h"
#include <CGAL/Sphere_traits.h>

typedef CGAL::Gmpq                                    	Number_type;
typedef CGAL::Cartesian<Number_type>                  	Kernel;
typedef CGAL::Sphere_traits<Kernel>			Traits_2;
typedef Kernel::Direction_3			        Direction_3;
typedef CGAL::SphereTopologicalMap<Kernel>		SphereType;
typedef CGAL::Spherical_map<SphereType, Traits_2>	SphereMap;

/*
 display one spherical face
 */
void displayFace(SphereMap::Face_handle &face) {
  std::cout << "displaying a face curves: " << std::endl;
  std::cout << "outer ccb:" << std::endl;
  SphereMap::Ccb_halfedge_circulator hit;
  if (face->does_outer_ccb_exist()) {
    // show outer ccb if exists
    hit = face->outer_ccb();
    SphereMap::Ccb_halfedge_circulator circEn=hit;
    do {
      std::cout << "edge: (" << hit->source()->direction() << ") - (" <<
	hit->target()->direction() << ")" << std::endl;
      ++hit;
    } while (hit !=circEn);
  } else {
    std::cout << "no outer ccb" << std::endl;
  }
  SphereMap::Holes_iterator holesIt;
  for (holesIt = face->holes_begin(); holesIt != face->holes_end(); ++holesIt) {
    std::cout << "a hole: " << std::endl;
    SphereMap::Ccb_halfedge_circulator hit;
    hit = *holesIt;
    SphereMap::Ccb_halfedge_circulator circEn=hit;
    do {
      std::cout << "edge: (" << hit->source()->direction() << ") - (" <<
	hit->target()->direction() << ")" << std::endl;
      ++hit;
    } while (hit !=circEn);
  }

}

int main() {
  SphereMap sMap;

  std::cout << "displaying the empty cube" << std::endl;
  SphereMap::Face_iterator fit2 = sMap.faces_begin();
  displayFace(fit2);
  std::cout << "finished displaying empty cube" << std::endl;

  SphereMap::X_monotone_curve_2 curv1(Direction_3(1,0,0), Direction_3(0,1,0));
  SphereMap::X_monotone_curve_2 curv2(Direction_3(0,1,0), Direction_3(-1,0,0));
  SphereMap::X_monotone_curve_2 curv3(Direction_3(-1,0,0), Direction_3(0,-1,0));
  SphereMap::X_monotone_curve_2 curv4(Direction_3(0,-1,0), Direction_3(1,0,0));

  std::cout << "inserting curve (1,0,0) - (0,1,0)" << std::endl;
  sMap.insert(curv1);
  std::cout << "inserting curve (0,1,0) - (-1,0,0)" << std::endl;
  sMap.insert(curv2);
  std::cout << "inserting curve (-1,0,0) - (0,-1,0)" << std::endl;
  sMap.insert(curv3);
  std::cout << "inserting curve (0,-1,0) - (1,0,0)" << std::endl;
  sMap.insert(curv4);

  std::cout << "displaying current map faces" << std::endl;
  // display faces
  SphereMap::Face_iterator fit;
  for (fit = sMap.faces_begin(); fit != sMap.faces_end(); ++fit) {
	displayFace(fit);
  }

  SphereMap::X_monotone_curve_2 curv5(Direction_3(1,-1,1), Direction_3(1,1,1));
  SphereMap::X_monotone_curve_2 curv6(Direction_3(-1,0,1), Direction_3(1,1,1));
  SphereMap::X_monotone_curve_2 curv7(Direction_3(-1,0,1), Direction_3(1,-1,1));
  std::cout << "inserting curve (1,-1,1) - (1,1,1)" << std::endl;
  sMap.insert(curv5);
  std::cout << "inserting curve (-1,0,1) - (1,1,1)" << std::endl;
  sMap.insert(curv6);
  std::cout << "inserting curve (-1,0,1) - (1,-1,1)" << std::endl;
  sMap.insert(curv7);

  std::cout << "displaying current map faces" << std::endl;
  for (fit = sMap.faces_begin(); fit != sMap.faces_end(); ++fit) {
    displayFace(fit);
  }

  return 0;
}
