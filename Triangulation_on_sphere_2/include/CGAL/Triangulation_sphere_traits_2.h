
#ifndef CGAL_TRIANGULATION_SPHERE_TRAITS_2_H
#define CGAL_TRIANGULATION_SPHERE_TRAITS_2_H


namespace CGAL { 
// REQUIRED PREDICATES
// colinear between 3D
//deux points p,q  dans un plan avec o, on veut savoir si q et a gauche ou a droite de op


template < typename K >
class Orientation_sphere_1
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Comparison_result        Comparison_result;

  Orientation_sphere_1(const Point_2& sphere);

  Comparison_result operator()(const Point_2& p, const Point_2& q) const
  {
    return coplanar_orientation(_sphere,p,q);
  }
	
 Comparison_result operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
 {
	return coplanar_orientation(_sphere,p,q,r);
	 //return coplanar_orientation(p,q,r,_sphere);
 }

	Comparison_result operator()(const Point_2& p, const Point_2& q, const Point_2& r,const Point_2& s) const
	{
		return coplanar_orientation(p,q,r,s);
		//return coplanar_orientation(p,q,r,_sphere);
	}

protected :
 Point_2  _sphere;
};
template < typename K >
Orientation_sphere_1<K>::
Orientation_sphere_1(const Point_2& sphere)
: _sphere(sphere)
{}

template < typename K >
class Orientation_sphere_2
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Comparison_result        Comparison_result;

  typedef Comparison_result   result_type;
  
  Orientation_sphere_2(const Point_2& sphere);

  Comparison_result operator()(const Point_2& p,
			       const Point_2& q,
			       const Point_2& test) const
  {
    return orientation(_sphere,p,q,test);
  }
	
	Comparison_result operator()(const Point_2& p, const Point_2& q,
								 const Point_2& r, const Point_2 & s) const
	{
		return orientation(p,q,r,s);
	}
	
	

protected :
 Point_2  _sphere;
};
template < typename K >
Orientation_sphere_2<K>::
Orientation_sphere_2(const Point_2& sphere)
: _sphere(sphere)
{}




template < typename K >
class Coradial_sphere_2
{
public:
  typedef typename K::Point_2                  Point_2;

  Coradial_sphere_2(const Point_2& sphere);

  bool operator()(const Point_2& p, const Point_2 q) const
  {
    return collinear(_sphere,p,q) &&
      ( are_ordered_along_line(_sphere,p,q) || are_ordered_along_line(_sphere,q,p) );
  }

protected :
 Point_2  _sphere;
};
template < typename K >
Coradial_sphere_2<K>::
Coradial_sphere_2(const Point_2& sphere)
: _sphere(sphere)
{}


template < typename K >
class Inside_cone_2
{
public:
  typedef typename K::Point_2                  Point_2;

  Inside_cone_2(const Point_2& sphere);

  bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
    if( collinear(_sphere,p,r)||collinear(_sphere,q,r)||orientation(_sphere,p,q,r)!=COLLINEAR)
      return false;
    if( collinear(_sphere,p,q) )
      return true;
   return coplanar_orientation(_sphere,p,q,r) == ( POSITIVE==coplanar_orientation(_sphere,q,p,r) );
  }

protected :
 Point_2  _sphere;
};
template < typename K >
Inside_cone_2<K>::
Inside_cone_2(const Point_2& sphere)
: _sphere(sphere)
{}


template < class R >
class Triangulation_sphere_traits_2 : public R 
{
 public :
  typedef Triangulation_sphere_traits_2<R>    Traits;   
  typedef R                                   K;
  typedef typename K::Point_3            Point_2;
  typedef typename K::Segment_3          Segment_2;
  typedef typename K::Triangle_3         Triangle_2;

  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;

  typedef typename K::Compare_x_3   Compare_x_3;
  typedef typename K::Compare_y_3   Compare_y_3;
	typedef typename K::Compare_z_3 Compare_z_3;

  typedef CGAL::Orientation_sphere_2<Traits> Orientation_2;
  typedef CGAL::Orientation_sphere_1<Traits> Orientation_1;
  typedef CGAL::Coradial_sphere_2<Traits>    Coradial_sphere_2;
  typedef CGAL::Inside_cone_2<Traits>        Inside_cone_2;
	
  
  Orientation_1
  orientation_1_object() const {
    return Orientation_1(_sphere);
  }

  Orientation_2
  orientation_2_object() const {
    return Orientation_2(_sphere);
  }

  Coradial_sphere_2
  coradial_sphere_2_object() const {
    return Coradial_sphere_2(_sphere);
  }

  Inside_cone_2
  inside_cone_2_object() const {
    return Inside_cone_2(_sphere);
  }
	
 
 

  
 public:
  Triangulation_sphere_traits_2(const Point_2& sphere=Point_2(0,0,0));  

 protected :
  Point_2 _sphere;
};

template < class R >
Triangulation_sphere_traits_2<R> ::
Triangulation_sphere_traits_2(const Point_2& sphere)
: _sphere(sphere)
{}

 


} //namespace CGAL

#endif // CGAL_TRIANGULATION_SPHERE_TRAITS_2_H

