#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Quotient.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/cartesian_homogeneous_conversion.h>
#include "visual_hull_creator.h"
#include <list>

typedef leda_integer NT;
typedef CGAL::Quotient<NT> CNT;
typedef CGAL::Cartesian<CNT> CKernel;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef visual_hull_creator<Nef_polyhedron> VHC;
typedef CKernel::FT FT;
typedef CKernel::Point_3 CPoint;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;

Point_3 read_point(std::ifstream& in) {
  double x,y,z;
  in >> x;
  in >> y;
  in >> z;
  CPoint p(x,y,z);
  return quotient_cartesian_to_homogeneous(p);
}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);

  std::ifstream in(argv[1]);

  VHC vhc(read_point(in), read_point(in));
  Nef_polyhedron N(Nef_polyhedron::COMPLETE);

  CGAL::Timer t;

  int ncameras;
  in >> ncameras;
  for(int cam=0; cam<ncameras; ++cam) {

    Point_3 camera(read_point(in));
    vhc.add_camera(camera);
    
    int npolygons;
    in >> npolygons;
    
    std::list<std::list<Point_3> > polygon_list;
    for(int poly=0; poly<npolygons; ++poly) {

      int npoints;
      in >> npoints;
      
      std::list<Point_3> input_points;
      for(int pnt=0; pnt<npoints; ++pnt)
	input_points.push_back(read_point(in));
      polygon_list.push_back(input_points);
    }

    std::list<std::list<Point_3> >::iterator li;
    for(li=polygon_list.begin(); li!=polygon_list.end(); ++li) {
      std::list<Point_3>::iterator pi(li->begin()), pimin(pi), pi_next,pi_prev;
      for(; pi!=li->end(); ++pi) {
	if(CGAL::lexicographically_xyz_smaller(*pi,*pimin))
	  pimin=pi;
      }
      pi_next=pi_prev=pimin;
      ++pi_next;
      if(pi_next==li->end()) pi_next=li->begin();
      if(pi_prev==li->begin()) pi_prev=li->end();
      --pi_prev;
      if(CGAL::orientation(*pi_prev,*pimin,*pi_next,camera) 
	 == CGAL::POSITIVE)
	li->reverse();
    }

    t.start();

    for(li=polygon_list.begin(); li!=polygon_list.end(); ++li)
      if(li==polygon_list.begin()) {
	vhc.add_outer_cycle_to_camera(li->begin(), li->end());
      } else {
	vhc.add_inner_cycle_to_camera(li->begin(), li->end());
      }

    for(li=polygon_list.begin(); li!=polygon_list.end(); ++li)
      if(li==polygon_list.begin()) {
	vhc.create_outer_cycles_opposites(li->begin(), li->end());
      } else {
	vhc.create_inner_cycles_opposites(li->begin(), li->end());
      }   
	
    vhc.recompute_scene();
    
    t.stop();
  }


  std::cerr << "Runtime Visual Hull :" << t.time() << std::endl;

  vhc.print_off();
}
