#ifndef CGAL_AOS3_SPHERE_3_TABLE_IMPL_H
#define CGAL_AOS3_SPHERE_3_TABLE_IMPL_H


#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table.h>




CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

template <class Traits_t>
void Sphere_3_table<Traits_t>::initialize_1() {
  di_ = geom_traits_object().intersect_3_object();
  spheres_.push_back(Sphere_3());
  spheres_.push_back(Sphere_3());
  spheres_.push_back(Sphere_3());
}

template <class Traits_t>
void Sphere_3_table<Traits_t>::initialize_2() {
  bbox_= CGAL::Bbox_3(std::numeric_limits<double>::max(),
		      std::numeric_limits<double>::max(),
		      std::numeric_limits<double>::max(),
		      -std::numeric_limits<double>::max(),
		      -std::numeric_limits<double>::max(),
		      -std::numeric_limits<double>::max());
  for (unsigned int i=3; i< spheres_.size(); ++i){
    bbox_= bbox_+ spheres_[i].bbox();
  }
  if (spheres_.size() <=3 ) {
    bbox_= CGAL::Bbox_3(0,0,0,0,0,0);
  }
  FT inf=16+2*std::max(bbox_.xmax(),
		    std::max(std::abs(bbox_.xmin()),
			     std::max(bbox_.ymax(),
				      std::max(std::abs(bbox_.ymin()),
					       std::max(bbox_.zmax(),
							std::abs(bbox_.zmin()))))));
  spheres_[1]=Sphere_3(Point_3(-inf, -inf, -inf), 0);
  spheres_[2]=Sphere_3(Point_3(inf, inf, inf), 0);
}


template <class Traits_t>
std::ostream &Sphere_3_table<Traits_t>::write(std::ostream &out) const {
  out << "Spheres:\n";
  out << Key(Key::BL) << ": " << operator[](Key(Key::BL)) << "\n";
  out << Key(Key::TR) << ": " << operator[](Key(Key::TR)) << "\n";
  for (Sphere_key_const_iterator it= sphere_keys_begin(); it != sphere_keys_end(); ++it) {
    out << *it << ": " << operator[](*it) << "\n";
  }
  out << "Bbox: " << bbox_ << std::endl;
  return out;
}

/* 
   Constructions--------------------------------------------------------
*/




template <class Traits_t>
typename Sphere_3_table<Traits_t>::Point_3 
Sphere_3_table<Traits_t>::center(Key a) const {
  return sphere(a).center();
}

template <class Traits_t>
typename Sphere_3_table<Traits_t>::FT 
Sphere_3_table<Traits_t>::squared_radius(Key a) const {
  return sphere(a).squared_radius();
}

template <class Traits_t>
typename Sphere_3_table<Traits_t>::Plane_3 
Sphere_3_table<Traits_t>::separating_plane(Key a, Key b) const {
  //Vector_3 d(center(b)-center(a));
  //FT nv[3];
  //nv[plane_coordinate(0).index()]=d[plane_coordinate(1).index()];
  //nv[plane_coordinate(1).index()]=-d[plane_coordinate(0).index()];
  //nv[sweep_coordinate().index()]=0;
  //Vector_3 n(-nv[0], -nv[1], -nv[2]);
  Plane_3 plane(center(a), center(b), center(a)-sweep_vector<Vector_3>());
  std::cout << "The plane for " << a << " and " << b << " is " << plane 
	    <<  "( " << plane.orthogonal_vector() << ")" << std::endl;
  return plane;
}


template <class Traits_t>
typename Sphere_3_table<Traits_t>::Plane_3 
Sphere_3_table<Traits_t>::equipower_plane(Key a, Key b) const {
  CGAL_precondition(a != b);
  Vector_3 n=2*(center(a)-center(b));
  FT d= discriminant(b) - discriminant(a);
  return Plane_3(n[0], n[1], n[2], d);
}

template <class Traits_t>
typename Sphere_3_table<Traits_t>::Point_3 
Sphere_3_table<Traits_t>::equipower_point(Key a, Key b) const {
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

template <class Traits_t>
typename Sphere_3_table<Traits_t>::Line_3 
Sphere_3_table<Traits_t>::equipower_line(Key a, Key b) const {
  Plane_3 eqp= equipower_plane(a,b);
  Plane_3 op(center(a), center(b), center(a)+sweep_vector<Vector_3>());
  Object o= di_(eqp, op);
  Line_3 l;
  if (assign(l,o)) {
    return l;
  } else {
    CGAL_assertion(0);
    return l;
  }
}

template <class Traits_t>
typename Sphere_3_table<Traits_t>::FT 
Sphere_3_table<Traits_t>::discriminant(Key i) const {
  Vector_3 v= center(i)-CGAL::ORIGIN;
  return v*v - sphere(i).squared_radius();
}


template <class Traits_t>
typename Sphere_3_table<Traits_t>::Sphere_3 Sphere_3_table<Traits_t>::sphere(Key ind) const {
  return spheres_[ind.internal_index()];
}

CGAL_AOS3_END_INTERNAL_NAMESPACE

#endif
