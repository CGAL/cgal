#ifndef CGAL_ARR_SPHERE_UTILITIES_H
#define CGAL_ARR_SPHERE_UTILITIES_H
#include <CGAL/Root_of_traits.h>

template <class K>
inline double approximate_intersection(const typename K::Sphere_3 &s, const typename K::Line_3 &l,  int C) {
  bool first= CGAL::sign(l.to_vector()[C]) == CGAL::POSITIVE;
  typename K::Vector_3 lp= l.point()-CGAL::ORIGIN;
  typename K::Vector_3 vc=s.center()-CGAL::ORIGIN; 
  typename K::Vector_3 lv=l.to_vector(); 
  double a=CGAL::to_double(lv*lv);
  double b=CGAL::to_double(2*lv*(lp-vc)); //-2*lv* vc + 2*lv*lp;
  double c=CGAL::to_double(lp*lp + vc*vc-s.squared_radius()-2*lp*vc);

  double disc= b*b-4*a*c;
  CGAL_assertion(disc >= -.1*b*b);
  if (disc <0) disc=0;

  double boa= b/(2.0*a);
  double discoa= std::sqrt(disc)/(2.0*a);
  double t0=-boa+discoa;
  double t1=-boa-discoa;

  double p0= CGAL::to_double(lp[C])+CGAL::to_double(lv[C])*t0;
  double p1= CGAL::to_double(lp[C])+t1*CGAL::to_double(lv[C]);
  //double v;
  if (first) return std::min(p0,p1);
  else return std::max(p0,p1);
}

template <class K, int C>
inline double approximate_intersection(const typename K::Sphere_3 &s, const typename K::Line_3 &l) {
  return approximate_intersection<K>(s,l, C);
}

template <class K>
inline typename CGAL::Root_of_traits<typename K::FT>::RootOf_2 exact_intersection(const typename K::Sphere_3 s, const typename K::Line_3 l, int C) {

  //sort on correct coordinate
  bool first= CGAL::sign(l.to_vector()[C]) == CGAL::POSITIVE;

  typename K::Vector_3 lp= l.point()-CGAL::ORIGIN;
  typename K::Vector_3 vc=s.center()-CGAL::ORIGIN; 
  typename K::Vector_3 lv=l.to_vector(); 
  typename K::FT a=lv*lv;
  typename K::FT b=2*lv*(lp-vc); //-2*lv* vc + 2*lv*lp;
  typename K::FT c=lp*lp + vc*vc-s.squared_radius()-2*lp*vc;

  typename K::FT disc= b*b-4*a*c;
  CGAL_assertion(disc >= 0);
  CGAL_precondition(a!=0);

  bool rt= first && CGAL::sign(lv[C]) == CGAL::POSITIVE 
    || !first && CGAL::sign(lv[C]) == CGAL::NEGATIVE;
  typename CGAL::Root_of_traits<typename K::FT>::RootOf_2 t=CGAL::make_root_of_2(a,b,c,rt);

   
  if (first) {
    CGAL_assertion(t*lv[C] <= lv[C] *CGAL::make_root_of_2(a,b,c, !rt));
  } else {
    CGAL_assertion(t*lv[C] >= lv[C] *CGAL::make_root_of_2(a,b,c, !rt)); 
  }
  typename CGAL::Root_of_traits<typename K::FT>::RootOf_2 p0= lp[C] + lv[C]*t;
  //std::cout << p0 << " is " << CGAL::to_double(p0) << std::endl;
  return p0;
}


template <class K, int C>
inline typename CGAL::Root_of_traits<typename K::FT>::RootOf_2 exact_intersection(const typename K::Sphere_3 &s, const typename K::Line_3 &l) {
  return exact_intersection(s, l, C);
}
  
  
// might specialize for x,y since they are always axis aligned lines I think
template <class K>
inline typename K::FT coord_on_line(const typename K::Line_3 &l, typename K::FT v, int C, int CC){
  typename K::FT t=(v-l.point()[C])/l.to_vector()[C];
  return  l.point()[CC]+t*l.to_vector()[CC];
}

template <class K, int C, int CC>
inline typename K::FT coord_on_line(const typename K::Line_3 &l, typename K::FT v){
  return coord_on_line<K>(l, v, C, CC);
}

// might specialize for x,y since they are always axis aligned lines I think
template <class K>
inline typename K::Point_3 point_on_line(const typename K::Line_3 &l, typename K::FT v, int C){
  CGAL_precondition(l.to_vector()[C] != 0);
  typename K::FT t=(v-l.point()[C])/l.to_vector()[C];
  return  l.point()+t*l.to_vector();
}

template <class K, int C>
inline typename K::Point_3 point_on_line(const typename K::Line_3 &l, typename K::FT v){
  return point_on_line<K>(l, v, C);
}

/*static inline CGAL::Bounded_side side_of_bounded_sphere(const Sphere &s, const Point_3 &p){
  Vector_3 d= (s.center()-p);
  return CGAL::Bounded_side(CGAL::sign(s.squared_radius()- d*d));
  }*/

template <class K>
inline typename K::FT closest_coord(const typename K::Point_3 &p, const typename K::Line_3 &l, int C) {
  typename K::FT t= (l.to_vector()*(p-CGAL::ORIGIN) - l.to_vector()*(l.point()-CGAL::ORIGIN))/(l.to_vector()*l.to_vector());
  return l.point()[C]+t*l.to_vector()[C];
}

template <class K, int C>
inline typename K::FT closest_coord(const typename K::Point_3 &p, const typename K::Line_3 &l) {
  return closest_coord<K>(p,l, C);
}


/*template <class C>
  static CGAL::Bounded_side bounded_side(Sphere_3 s,
  const Sphere_point_3 &p) {
  Sphere_point_3 ls(s, p.line(), true);
  if (!ls.is_finite()) 
  return CGAL::ON_UNBOUNDED_SIDE;
  Sphere_point_3 rs(s, p.line(), false);
  if (ls > rs) std::swap(ls,rs);
  if (ls < p && rs > p) return CGAL::ON_BOUNDED_SIDE;
  else if (ls > p || rs < p ) return CGAL::ON_UNBOUNDED_SIDE;
  else return CGAL::ON_BOUNDARY;
  }*/

#endif
