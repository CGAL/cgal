#ifndef CGAL_MESH_SIZE_TRAITS_2_H
#define CGAL_MESH_SIZE_TRAITS_2_H

#include <CGAL/Mesh_default_traits_2.h>

namespace CGAL {

template <class K>
class Mesh_size_traits_2 : public Mesh_default_traits_2<K>
{
  double sizebound;

public:
  typedef Mesh_default_traits_2<K> Base;

  Mesh_size_traits_2(const double aspect_bound = 0.125, 
		     const double size_bound = 0)
    : Base(aspect_bound), sizebound(size_bound) {};

  inline
  double size_bound() const { return sizebound; };

  inline
  void set_size_bound(const double sb) { sizebound = sb; };

  class Is_bad: public Base::Is_bad
  {
    const double SB;
  public:
    typedef typename Base::Is_bad::Point_2 Point_2;
    typedef typename Base::Is_bad::Traits Traits;

    Is_bad(const double aspect_bound,
	   const double size_bound)
      : Base::Is_bad(aspect_bound), SB(size_bound) {};

    bool operator()(const typename Base::Point_2& a,
		    const typename Base::Point_2& b,
		    const typename Base::Point_2& c) const
    {
      if(Base::Is_bad::operator()(a,b,c))
	return true;
      if(SB==0)
	return false;
      else
	{
	  typename Traits::Compute_squared_distance_2 distance = 
	    Traits().compute_squared_distance_2_object();
	  const double 
	    dab = to_double(distance(a,b)), 
	    dbc = to_double(distance(b,c)),
	    dca = to_double(distance(c,a));

	  return ((dab > SB*SB) || (dbc > SB*SB) || (dca > SB*SB));
	}
    };
  };

  Is_bad is_bad_object() const
  { return Is_bad(bound(), size_bound()); }
};

} //end namespace

#endif
