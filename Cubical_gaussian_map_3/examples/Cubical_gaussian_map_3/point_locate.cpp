/*
  agc course project,
  spherical arrangements of none intersecting arcs of great circles on a sphere

  an example for point location of spherical objects

  author: Kapelushnik Lior
*/

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Cubical_gaussian_map_3.h>
#include <CGAL/Sphere_traits.h>
#include "CGAL/Spherical_map.h"
#include <CGAL/Sphere_point_location.h>

typedef CGAL::Gmpq                                    	Number_type;
typedef CGAL::Cartesian<Number_type>                  	Kernel;
typedef CGAL::Sphere_traits<Kernel>                     Traits_2;
typedef Kernel::Direction_3			        Direction_3;
typedef CGAL::SphereTopologicalMap<Kernel>		SphereType;
typedef CGAL::Cubical_gaussian_map_3<Kernel,CGAL::Spherical_cgm_arr_dcel> CGM;
typedef CGAL::Spherical_map<SphereType, Traits_2>	SphereMap;


//typedef Sphere_naive_point_location<SphereMap>		PointLocation;
typedef CGAL::Sphere_walk_along_line_point_location<SphereMap>	PointLocation;

int main() {

  SphereMap sMap;

  std::cout << "inserting curve (1,1,0) - (1,1,1)" << std::endl;
  SphereMap::X_monotone_curve_2 curvy(Direction_3(1,1,0), Direction_3(1,1,1));
  sMap.insert(curvy);
  std::cout << "inserting curve (1,1,1) - (0,0,1)" << std::endl;
  SphereMap::X_monotone_curve_2 curvy2(Direction_3(1,1,1), Direction_3(0,0,1));
  sMap.insert(curvy2);
  std::cout << "inserting curve (0,2,0) - (1,0,0)" << std::endl;
  SphereMap::X_monotone_curve_2 curvy3(Direction_3(0,2,0), Direction_3(1,0,0));
  sMap.insert(curvy3);
  std::cout << "inserting curve (0,0,1) - (1,0,0)" << std::endl;
  SphereMap::X_monotone_curve_2 curvy4(Direction_3(0,0,1), Direction_3(1,0,0));
  sMap.insert(curvy4);
  std::cout << "inserting curve (0,1,0) - (0,0,1)" << std::endl;
  SphereMap::X_monotone_curve_2 curvy5(Direction_3(0,1,0), Direction_3(0,0,1));
  sMap.insert(curvy5);

  std::cout << sMap; // print map

  std::cout << std::endl;

  PointLocation pl(&sMap);
  std::cout << "point locating direction (0,1,0)" << std::endl;
  Direction_3 dirr(0,1,0);
  CGAL::Object obj = pl.locate(dirr);
  SphereMap::Vertex_handle foundVer;
  if (CGAL::assign(foundVer, obj)) {
    std::cout << "found a vertex with direction: " << foundVer->direction() <<
      std::endl;
  }
  std::cout << "end locating (0,1,0)" << std::endl;

  std::cout << "point locating direction (-1,11,-1)" << std::endl;
  Direction_3 dirr2k(-1,11,-1);
  obj = pl.locate(dirr2k);
  SphereMap::Face_handle foundFace;
  if (CGAL::assign(foundFace, obj)) {
    std::cout << "found a face" << std::endl;
    std::cout << "the face outer ccb" << std::endl;
    SphereMap::Ccb_halfedge_circulator heif;
    heif = foundFace->outer_ccb();
    SphereMap::Ccb_halfedge_circulator circEn = heif;
    do {
      std::cout << heif->curve() << std::endl;
      ++heif;
    } while (heif != circEn);
  }
  std::cout << "end locating (-1,11,-1)" << std::endl;

  std::cout << "point locating direction (0,1,1)" << std::endl;
  Direction_3 dirr4k(0,1,1);
  obj = pl.locate(dirr4k);
  SphereMap::Halfedge_handle foundHalf;
  if (CGAL::assign(foundHalf, obj)) {
    std::cout << "found a halfedge with curve: " << foundHalf->curve() <<
      std::endl;
  }
  std::cout << "end locating (0,1,1)" << std::endl;

  std::cout << "point locating direction (1/2,1,1)" << std::endl;
  Direction_3 dirr7k(0.5,1,1);
  obj = pl.locate(dirr7k);

  if (CGAL::assign(foundFace, obj)) {
    std::cout << "found a face" << std::endl;
    std::cout << "the face outer ccb" << std::endl;
    SphereMap::Ccb_halfedge_circulator heif;
    heif = foundFace->outer_ccb();
    SphereMap::Ccb_halfedge_circulator circEn = heif;
    do {
      std::cout << heif->curve() << std::endl;
      ++heif;
    } while (heif != circEn);
  }
  std::cout << "end locating (1/2,1,1)" << std::endl;
  return 0;
}
