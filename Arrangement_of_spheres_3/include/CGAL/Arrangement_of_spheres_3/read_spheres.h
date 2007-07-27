#ifndef CGAL_AOS3_READ_CIRCLES_H
#define CGAL_AOS3_READ_CIRCLES_H

#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Simple_cartesian.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

template <class K, bool VERBOSE>
void read_spheres(std::istream &in, std::vector<typename K::Sphere_3> &out) {
  typedef CGAL::Simple_cartesian<double> DK;
  CGAL::Bbox_3 bbox(std::numeric_limits<double>::max(),
		    std::numeric_limits<double>::max(),
		    std::numeric_limits<double>::max(),
		    -std::numeric_limits<double>::max(),
		    -std::numeric_limits<double>::max(),
		    -std::numeric_limits<double>::max());
		    
  while (true){
    char buf[1000];
    in.getline(buf, 1000);
    if (!in) break;
    {
      std::istringstream iss(buf);
      if (buf[0]=='#') continue;
      typename K::Sphere_3 ns;
      iss >> ns;
      if (iss) {
	bbox= bbox+ns.bbox();	
	out.push_back(ns);
      } else {
	std::istringstream iss2(buf);
	typename DK::Sphere_3 ns;
	iss2 >> ns;
	if (iss2) {
	  bbox= bbox+ns.bbox();
	  out.push_back(typename K::Sphere_3(typename K::Point_3(ns.center().x(), ns.center(). y(),
								 ns.center().z()), ns.squared_radius()));
	} else {
	  std::cerr << "Can't parse line " << buf << std::endl;
	  break;
	}
      }
    }
  }
  if (VERBOSE) {
    std::cout << "Bounding box is " << bbox << std::endl;
  }
}

CGAL_AOS3_END_INTERNAL_NAMESPACE

#endif
