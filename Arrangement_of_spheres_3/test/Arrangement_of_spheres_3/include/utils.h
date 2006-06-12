NT squared_depth(Point_3 p, const NT &z){
  return CGAL::square(p.z()-z);
}

Circle intersect(Sphere s, const NT &z){
  NT r2= s.squared_radius()- squared_depth(s.center(), z);
  CGAL_assertion(r2>=0);
  Circle c(Point_2(s.center().x(), s.center().y()), r2);
  return c;
}


Point_2 project(Point_3 s) {
  return Point_2(s.x(), s.y());
}



Vector_2 project(Vector_3 s) {
  return Vector_2(s.x(), s.y());
}

Point_3 unproject(Point_2 s, NT z) {
  return Point_3(s.x(), s.y(), z);
}

Line_2 project(Line_3 s) {
  return Line_2(project(s.point()), project(s.to_vector()));
} 

bool has_overlap(const Sphere &s, const NT &z) {
  NT dz2= squared_depth(s.center(), z);
  if (dz2 < s.squared_radius()) return true;
  else return false;
}



struct Compare_spheres {
  bool operator()(const Sphere &a, const Sphere &b) const {
    double ra= std::sqrt(CGAL::to_double(a.squared_radius()));
    double rb= std::sqrt(CGAL::to_double(b.squared_radius()));
    return a.center().z() - NT(ra)  < b.center().z() - NT(rb); 
  }
};

template <int D, int C>
void clip(const Sphere &s, const Line_3 &l,
	  std::pair<Sphere_line_intersection<C>, Sphere_line_intersection<C> > &seg) {
  //std::cout << "Segment is (" << seg.first << ", " << seg.second << ")" << std::endl;
 
  Sphere_line_intersection<C> i0(s, l, true);
  Sphere_line_intersection<C> i1(s, l, false);
  CGAL_assertion(i0 <= i1);
  if (i0.is_finite()) {
    CGAL::Root_of_traits<NT>::RootOf_2 ri0= i0.exact_coordinate(C);
    CGAL::Root_of_traits<NT>::RootOf_2 ri1= i1.exact_coordinate(C);
    
    CGAL_assertion(ri0 <= ri1);
  }
  //std::cout << "Intersections are " << i0 << " and " << i1 << std::endl;
  if (i0 > i1) std::swap(i0, i1);
  bool b0f= i0 < seg.first;
  bool b0s= i0 < seg.second;
  bool b1f= i1 < seg.first;
  bool b1s= i1 < seg.second;
  if (D==0) {
    CGAL_assertion(seg.first < seg.second);
    if (b0f && !b1f && b1s){
      seg.second= i1;
    } else if (!b0f && b0s) {
      seg.second=i0;
    }
  } else {
    CGAL_assertion(seg.second < seg.first);
    if (b1f && !b1s){
      seg.second= i1;
    } else if (!b1f && !b0s && b0f) {
      seg.second=i0;
    }
  }
}



template <int CA, int CB>
void clip_seg(NT z,
	      const Line_3 &la,
	      const std::pair<Sphere_line_intersection<CB>, Sphere_line_intersection<CB> > &b,
	      const Line_3 &lb,
	      std::pair<Sphere_line_intersection<CA>, Sphere_line_intersection<CA> > &a){
  Line_2 la2(project(la.point()), project(la.to_vector()));
  Line_2 lb2(project(lb.point()), project(lb.to_vector()));
  CGAL::Object oi= CGAL::intersection(la2, lb2);
  Point_2 pi;
  if (CGAL::assign(pi, oi)){
    Sphere_line_intersection<CA> ii(Point_3(pi.x(), pi.y(), z));
    if (ii > a.first && ii < a.second || ii > a.second && ii < a.first){
      Sphere_line_intersection<CB> iib(Point_3(pi.x(), pi.y(), z));
      if (iib > b.first && iib < b.second || iib > b.second && iib < b.first){
	a.second=ii;
      } else {
	/*std::cout << CGAL::to_double(iib) << " is out of range b (" << CGAL::to_double(b.first) 
	  << ", " << CGAL::to_double(b.second) << ")" << std::endl;*/
      }
      // check b range too
    } else {
      /*std::cout << CGAL::to_double(ii) << " is out of range a (" << CGAL::to_double(a.first) 
	<< ", " << CGAL::to_double(a.second) << ")" << std::endl;*/
    }
  } else {
    /*std::cout << "No intersection point for lines " << la << " and " << lb << std::endl;
      std::cout << "As " << la2 << " and " << lb2 << std::endl;*/
  }
}
