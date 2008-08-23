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
typedef Nef_polyhedron::SNC_structure SNC_structure;
typedef CGAL::visual_hull_creator<SNC_structure> VHC;
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
  std::list<Nef_polyhedron> N_list;
  
  CGAL::Timer t;

  Point_3 room_min = read_point(in);
  Point_3 room_max = read_point(in);

  int ncameras;
  in >> ncameras;
  for(int cam=0; cam<ncameras; ++cam) {

    Point_3 camera(read_point(in));
    
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

    Nef_polyhedron N;
    VHC vhc(room_min, room_max, camera, polygon_list);
    N.delegate(vhc);
    N_list.push_back(N);
    t.stop();

    std::cerr << "intermediate time " << t.time() << std::endl;
  }

  t.start();
  Nef_polyhedron N1,N2;
  while(N_list.size() > 1) {
    N1 = N_list.front();
    N_list.pop_front();
    N2 = N_list.front();
    N_list.pop_front();
    N1.is_valid();
    N2.is_valid();
    N_list.push_back(N1*N2);
  }
  t.stop();

  std::cerr << "Runtime Visual Hull :" << t.time() << std::endl;
  
  Nef_polyhedron result(N_list.front());
  //  std::cerr << result;

  /*
  QApplication a(argc,argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
  */
}
