namespace CGAL {

template <class K>
class Mesh_default_traits_2 : public K
{
  double B;
public:

  Mesh_default_traits_2(const double bound = 0.125) : B(bound) {};

  inline
  double bound() const { return B; };

  inline 
  void set_bound(const double bound) { B = bound; };

  class Is_bad;
  Is_bad is_bad_object() const
    { return Is_bad(B); }

  class Compute_squared_minimum_sine_2;
  Compute_squared_minimum_sine_2 
  compute_squared_minimum_sine_2_object() const
    { return Compute_squared_minimum_sine_2(); }
};

template <class K>
class Mesh_default_traits_2<K>::Is_bad 
{
  const double B;
public:
  typedef typename K::Point_2 Point_2;
  typedef Mesh_default_traits_2<K> Traits;
  
  Is_bad(const double bound) : B(bound) {};
  
  bool operator()(const Point_2& a,
		  const Point_2& b,
		  const Point_2& c) const
    {
      Compute_squared_minimum_sine_2 squared_sine = 
	Traits().compute_squared_minimum_sine_2_object();
      
      if( squared_sine(a,b,c) >= B )
	return false;
      else
	return true;
    };
};

template <class K>
class Mesh_default_traits_2<K>::Compute_squared_minimum_sine_2 {
public:
  typedef typename K::Point_2 Point_2;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Compute_area_2 Compute_area_2;
  typedef typename K::Compute_squared_distance_2
  Compute_squared_distance_2;
  typedef typename K::FT FT;

  double operator()(const Point_2& pa,
		    const Point_2& pb,
		    const Point_2& pc) const
    {
      K k;
      Compute_area_2 area_2 = k.compute_area_2_object();
      Compute_squared_distance_2 squared_distance = 
	k.compute_squared_distance_2_object();

      Triangle_2 t = k.construct_triangle_2_object()(pa,pb,pc);
      double area = 2*to_double(area_2(t));
      area=area*area;

      double
	a = squared_distance(pb, pc),
	b = squared_distance(pc, pa),
	c = squared_distance(pa, pb);

      if(a<b)
	if(a<c)
	  return area/(b*c);
	else
	  return area/(a*b);
      else
	if(b<c)
	  return area/(a*c);
	else
	  return area/(a*b);
    }
};

}; // end namespace CGAL
