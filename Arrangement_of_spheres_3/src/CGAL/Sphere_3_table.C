#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table.h>
#include <CGAL/Arrangement_of_spheres_3/coordinates.h>

void Sphere_3_table::initialize_1() {
  di_ = geometric_traits_object().intersect_3_object();
  spheres_.push_back(Sphere_3());
  spheres_.push_back(Sphere_3());
  spheres_.push_back(Sphere_3());
}


void Sphere_3_table::initialize_2() {
  bbox_= CGAL::Bbox_3(std::numeric_limits<double>::max(),
		      std::numeric_limits<double>::max(),
		      std::numeric_limits<double>::max(),
		      -std::numeric_limits<double>::max(),
		      -std::numeric_limits<double>::max(),
		      -std::numeric_limits<double>::max());
  for (unsigned int i=3; i< spheres_.size(); ++i){
    bbox_= bbox_+ spheres_[i].bbox();
  }
  NT inf=2*std::max(bbox_.xmax(),
		    std::max(std::abs(bbox_.xmin()),
			     std::max(bbox_.ymax(),
				      std::max(std::abs(bbox_.ymin()),
					       std::max(bbox_.zmax(),
							std::abs(bbox_.zmin()))))));
  spheres_[1]=Sphere_3(Point_3(-inf, -inf, -inf), 0);
  spheres_[2]=Sphere_3(Point_3(inf, inf, inf), 0);
}


/* 
   Constructions--------------------------------------------------------
*/





Sphere_3_table::Point_3 
Sphere_3_table::center(Key a) const {
  return sphere(a).center();
}


Sphere_3_table::Plane_3 
Sphere_3_table::separating_plane(Key a, Key b) const {
  Vector_3 d(center(b)-center(a));
  NT nv[3];
  nv[plane_coordinate(0).index()]=-d[plane_coordinate(1).index()];
  nv[plane_coordinate(1).index()]= d[plane_coordinate(0).index()];
  nv[sweep_coordinate().index()]=0;
  Vector_3 n(nv[0], nv[1], nv[2]);
  Plane_3 plane(center(b), n);
  std::cout << "The plane for " << a << " and " << b << " is " << plane << std::endl;
  return plane;
}



Sphere_3_table::Plane_3 
Sphere_3_table::equipower_plane(Key a, Key b) const {
  CGAL_precondition(a != b);
  Vector_3 n=2*(center(a)-center(b));
  NT d= discriminant(b) - discriminant(a);
  return Plane_3(n[0], n[1], n[2], d);
}

Sphere_3_table::Point_3 
Sphere_3_table::equipower_point(Key a, Key b) const {
  Plane_3 eqp= equipower_plane(a,b);
  Line_3 l(center(a), (center(a)-center(b)));
  CGAL::Object o= di_(eqp, l);
  Point_3 pt;
  if (CGAL::assign(pt, o)) {
    return pt;
  } else {
    CGAL_assertion(center(a) == center(b));
    //Line_3 l(center(a), Vector_3(0,0,1));
    // d^2= r^2
    throw Equal_centers_exception();
  }
}

Sphere_3_table::NT 
Sphere_3_table::discriminant(Key i) const {
  Vector_3 v= center(i)-CGAL::ORIGIN;
  return v*v - sphere(i).squared_radius();
}



Sphere_3_table::Sphere_3 Sphere_3_table::sphere(Key ind) const {
  return spheres_[ind.internal_index()];
}


