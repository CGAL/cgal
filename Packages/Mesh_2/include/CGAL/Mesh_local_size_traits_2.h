#include <CGAL/Mesh_size_traits_2.h>

namespace CGAL {

template <class K>
class Mesh_local_size_traits_2 : public Mesh_size_traits_2<K>
{
public:
  typedef Mesh_size_traits_2<K> Base;
  typedef typename Base::Base PreviousBase;
  typedef typename K::FT FT;
  typedef typename K::Point_2 Point;

private:
  bool local;
  Point _p;

public:
  Mesh_local_size_traits_2(const double aspect_bound = 0.125, 
			   const double size_bound = 0,
			   const bool is_local_size = false,
			   const Point p = Point())
    : Base(aspect_bound, size_bound), local(is_local_size), _p(p) {};

  inline
  Point point() const { return _p; };

  inline 
  void set_point(const Point p) { _p = p; };

  inline
  bool is_local_size() const { return local; };

  inline 
  void set_local_size(bool local_size) { local = local_size; };

  class Is_bad: public Base::Is_bad
  {
  public:
    typedef typename Base::Is_bad Baseclass;
    typedef typename Baseclass::Point_2 Point_2;
    typedef typename Baseclass::Traits Traits;
    typedef typename Base::Base PreviousBase;
    typedef typename PreviousBase::Is_bad PreviousBaseclass;

  private:
    const bool local;
    const Point_2 p;

  public:
    Is_bad(const double aspect_bound,
	   const double size_bound,
	   const bool l,
	   const Point_2 _p)
      : Base::Is_bad(aspect_bound, size_bound), local(l), p(_p) {};

    bool operator()(const Point_2& a,
		    const Point_2& b,
		    const Point_2& c) const
    {
      if(!local)
	return Base::Is_bad::operator()(a,b,c);
      else
	{
	  typename Traits::Orientation_2 orient = 
	    Traits().orientation_2_object();
	  
	  if(PreviousBaseclass::operator()(a,b,c))
	    return true;
	  Orientation 
	    o1 = orient(a,b,p),
	    o2 = orient(b,c,p),
	    o3 = orient(c,a,p);
	  
	  if((o1==o2) && (o2==o3))
	    return Base::Is_bad::operator()(a,b,c);
	  else
	    return false;
	};
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(bound(), size_bound(), local, point()); };
};

}; //end namespace
