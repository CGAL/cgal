#ifndef CGAL_MESH_AREA_TRAITS_2_H
#define CGAL_MESH_AREA_TRAITS_2_H

#include <CGAL/Mesh_default_traits_2.h>

namespace CGAL {

template <class K>
class Mesh_area_traits_2 : public Mesh_default_traits_2<K>
{
  double areabound;

public:
  typedef Mesh_default_traits_2<K> Base;

  Mesh_area_traits_2(const double aspect_bound = 0.125, 
		     const double area_bound = 0)
    : Base(aspect_bound), areabound(area_bound) {};

  inline
  double area_bound() const { return areabound; };

  inline
  void set_area_bound(const double ab) { areabound = ab; };

  class Is_bad: public Base::Is_bad
  {
    const double AB;
  public:
    typedef typename K::Point_2 Point_2;
    typedef typename K::Triangle_2 Triangle_2;

    Is_bad(const double aspect_bound,
	   const double area_bound)
      : Base::Is_bad(aspect_bound), AB(area_bound) {};

    bool operator()(const typename Base::Point_2& a,
		    const typename Base::Point_2& b,
		    const typename Base::Point_2& c) const
    {
      if(Base::Is_bad::operator()(a,b,c))
	return true;
      if(AB==0)
	return false;
      else
	{
	  typename K::Compute_area_2 area =
	    K().compute_area_2_object();
	  return area(Triangle_2(a, b, c)) > AB;
	}
    };
  };

  Is_bad is_bad_object() const
  { return Is_bad(bound(), area_bound()); }
};

} //end namespace

#endif
