#ifndef DELAUNAY_MESH_SIZING_FIELD_CRITERIA_2_H
#define DELAUNAY_MESH_SIZING_FIELD_CRITERIA_2_H

#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>

#include <CGAL/Lipschitz_sizing_field_2.h>


namespace CGAL
{

template <class CDT, class SF>
class Lipschitz_sizing_field_criteria_2
  : public virtual Delaunay_mesh_criteria_2<CDT>
{

public:
    
  typedef SF Sizing_field;
    
protected:
    
  Sizing_field* sizing_field;

  typedef typename CDT::Geom_traits Geom_traits;
  Geom_traits traits;
    
public:
    
  typedef Delaunay_mesh_criteria_2<CDT> Base;
    
  Lipschitz_sizing_field_criteria_2(const double aspect_bound = 0.125,
				    Sizing_field* sf = 0,
                                    const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound), sizing_field(sf), traits(traits)
  {}

  Lipschitz_sizing_field_criteria_2& operator =(const Lipschitz_sizing_field_criteria_2<CDT,SF>& c)
  {
    if(&c == this) return *this;
    this->sizing_field = c.sizing_field;
    this->traits = c.traits;
    return *this;
  }

  inline const Sizing_field* sizing_field_object()
  {
    return sizing_field;
  }

  struct Quality : public std::pair<double, double>
  {
    typedef std::pair<double, double> Base;
      
    Quality() : Base() {};
    Quality(double _sine, double _size) : Base(_sine, _size) {};
      
    const double& size() const { return second; }
    const double& sine() const { return first; }
      
    // q1<q2 means q1 is prioritised over q2
    // ( q1 == *this, q2 == q )
    bool operator<(const Quality& q) const
    {
      if( size() > 1 )
	if( q.size() > 1 )
	  return ( size() > q.size() );
	else
	  return true; // *this is big but not q
      else
	if( q.size() >  1 )
	  return false; // q is big but not *this
      return( sine() < q.sine() );
    }
  };
    
  class Is_bad : public Base::Is_bad
  {
  protected:

    Sizing_field* sizing_field;

  public:

    typedef typename Base::Is_bad::Point_2 Point_2;

    Is_bad(const double aspect_bound, Sizing_field* sf, 
           const Geom_traits& traits)
      : Base::Is_bad(aspect_bound, traits), sizing_field(sf)
    {
    }

    Mesh_2::Face_badness operator ()(const typename CDT::Face_handle& fh, Quality& q)
    {
      typedef typename CDT::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2 Compute_squared_distance_2;
	    
      Geom_traits traits; /** @warning traits with data!! */
	    
      Compute_squared_distance_2 squared_distance = traits.compute_squared_distance_2_object();

      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();
	    
      double
	a = CGAL::to_double(squared_distance(pb, pc)),
	b = CGAL::to_double(squared_distance(pc, pa)),
	c = CGAL::to_double(squared_distance(pa, pb));
	    
      double max_length; // squared max edge length
      double second_max_length;
	    
      if(a < b)
	{
	  if(b < c) 
	    {
	      max_length = c;
	      second_max_length = b;
	    }
	  else 
	    { // c<=b
	      max_length = b;
	      second_max_length = ( a < c ? c : a );
	    }
	}
      else // b<=a
	{
	  if(a < c) 
	    {
	      max_length = c;
	      second_max_length = a;
	    }
	  else 
	    { // c<=a
	      max_length = a;
	      second_max_length = ( b < c ? c : b );
	    }
	}

      q.second = 0;

      // compute the sizing field at the centroid of the vertices
      const Point_2& pg = centroid(pa, pb, pc);
      double squared_size_bound = sizing_field->query(pg);
      squared_size_bound *= squared_size_bound;

      if( squared_size_bound != 0 )
	{
	  q.second = max_length / squared_size_bound;
	  // normalized by size bound to deal
	  // with size field
	  if( q.size() > 1 )
	    {
	      q.first = 1; // (do not compute sine)
	      //std::cout << "imperatively bad!!\n";
	      return Mesh_2::IMPERATIVELY_BAD;
	    }
	}

      Compute_area_2 area_2 = traits.compute_area_2_object();
	    
      double area = 2*CGAL::to_double(area_2(pa, pb, pc));
	    
      q.first = (area * area) / (max_length * second_max_length); // (sine)
	    
      if( q.sine() < this->B )
	{
	  //std::cout << "bad!\n";
	  return Mesh_2::BAD;
	}
      else
	{
	  //std::cout << "not bad.\n";
	  return Mesh_2::NOT_BAD;
	}
    }
	
    Mesh_2::Face_badness operator ()(const Quality q) const
    {
      if(q.size() > 1)
	return Mesh_2::IMPERATIVELY_BAD;
      if(q.sine() < this->B)
	return Mesh_2::BAD;
      else
	return Mesh_2::NOT_BAD;
    }

  };

  Is_bad is_bad_object() const
  { 
    return Is_bad(this->bound(), sizing_field, traits);
  }

};

}

#endif
