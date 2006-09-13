#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <CGAL/Random.h>
#include <iomanip>


struct Do_work {
  void operator()(/*Qt_examiner_viewer_2 *q*/){
    CGAL::Bbox_3 box(std::numeric_limits<double>::max(),
		     std::numeric_limits<double>::max(),
		     std::numeric_limits<double>::max(),
		     -std::numeric_limits<double>::max(),
		     -std::numeric_limits<double>::max(),
		     -std::numeric_limits<double>::max());

    typedef Arrangement_of_spheres_traits_3::Geometric_traits K;
    typedef CGAL::Simple_cartesian<double> DK;
    std::vector<K::Sphere_3> spheres;
    std::ifstream in("data/random_point_locate.spheres");
    while (true){
      char buf[1000];
      in.getline(buf, 1000);
      if (!in) break;
      std::istringstream iss(buf);
      DK::Sphere_3 s;
      iss >> s;
      if (!iss) {
	std::cerr << "Can't parse line " << buf << std::endl;
      } else {
	spheres.push_back(K::Sphere_3(K::Point_3(s.center().x(), 
						 s.center().y(),
						 s.center().z()), 
				      s.squared_radius()));
	box= box+ spheres.back().bbox();
      }
      //std::cout << spheres.back() << std::endl;
    }


    std::cout << "Read " << spheres.size() << " spheres." << std::endl;
    std::cout << "Bounding box is from " << box.zmin() << " to " << box.zmax()
	      << std::endl;
    //

    CGAL::Random rand;

    double z= rand.get_double(box.zmin(), box.zmax());
    std::cout << std::setprecision(15);
    std::cout << "z is " << z << std::endl;

    Slice::T tr(spheres.begin(), spheres.end());
    Slice slice(tr);
    slice.initialize_at(z);

    for (unsigned int nit=0; nit != 100; ++nit){
      double x,y, radius;
      x= rand.get_double(box.xmin(), box.xmax());
      y= rand.get_double(box.ymin(), box.ymax());
      y= rand.get_double(0,4);
      K::Sphere_3 s(K::Point_3(x, y, z+radius), radius);
      K::FT lp[3];
      lp[plane_coordinate(0).index()]=x;
      lp[plane_coordinate(1).index()]=y;
      lp[sweep_coordinate().index()]=z;
      Arrangement_of_spheres_traits_3::Sphere_point_3 sp(s, 
							 K::Line_3(K::Point_3(lp[0], lp[1], lp[2]),
								   sweep_vector<K::Vector_3>()));
      std::cout << "Inserting " << s << std::endl;
      Slice::T::Key k=slice.debug_new_sphere(s);
      slice.insert_sphere(sp, k);
    }
  }
};


int main(int argc, char *argv[]){
  
  Do_work dw;
  dw();

  //Qt_debug_viewer_2<Do_work> qtd(dw, argc, argv);
  
  return EXIT_SUCCESS;
}
