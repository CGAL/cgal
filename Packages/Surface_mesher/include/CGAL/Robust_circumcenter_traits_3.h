// Author(s)     :  Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
// This traits override Construct_circumcenter_3 
// to get a robust ones when called on slivers

#ifndef CGAL_ROBUST_CIRCUMCENTER_TRAITS_3_H
#define CGAL_ROBUST_CIRCUMCENTER_TRAITS_3_H

#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

CGAL_BEGIN_NAMESPACE 


template < class K >
class Robust_construct_circumcenter_3
{
  typedef Exact_predicates_exact_constructions_kernel EK;
  typedef typename K::Point_3         Point_3;

  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 4 >   Arity;

  Point_3 operator() (const Point_3 & p,
		      const Point_3 & q,
		      const Point_3 & r,
		      const Point_3 & s) const
    {
   
      Cartesian_converter<K,EK> to_exact;
      Cartesian_converter<EK,K, To_double<EK::FT> > back_from_exact;
      //Cartesian_converter<EK,K > back_from_exact;
      EK::Construct_circumcenter_3 
	exact_circumcenter = EK().construct_circumcenter_3_object();

      return back_from_exact(exact_circumcenter( to_exact(p),
						 to_exact(q),
						 to_exact(r),
						 to_exact(s)));

    }

  Point_3 operator() (const Point_3 & p,
		      const Point_3 & q,
		      const Point_3 & r) const
    {
      Cartesian_converter<K,EK> to_exact;
      Cartesian_converter<EK,K, To_double<EK::FT> > back_from_exact; 
      EK::Construct_circumcenter_3 
	exact_circumcenter = EK().construct_circumcenter_3_object();

      return back_from_exact(exact_circumcenter( to_exact(p),
						 to_exact(q),
						 to_exact(r)));
    }
};




template < class K>
class Robust_circumcenter_traits_3
  : public K
{
 public:
  typedef CGAL::Robust_construct_circumcenter_3<K> Construct_circumcenter_3;

  Construct_circumcenter_3 
  construct_circumcenter_3_object() const
  { return Construct_circumcenter_3(); }
  
};



CGAL_END_NAMESPACE

#endif //CGAL_ROBUST_CIRCUMCENTER_TRAITS_3_H
