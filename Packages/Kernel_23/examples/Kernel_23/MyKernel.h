
#ifndef MYKERNEL_H
#define MYKERNEL_H

#include <CGAL/Cartesian/Cartesian_base.h>
#include <CGAL/Simple_Handle_for.h>


#include "MyPointC2.h"
#include "MysegmentC2.h"



  template <typename K>
  class MyConstruct_point_2
  {
    typedef typename K::RT         RT;
    typedef typename K::Point_2    Point_2;
  public:
    typedef Point_2          result_type;
    typedef CGAL::Arity_tag< 1 >   Arity;

    Point_2
    operator()() const
    { return Point_2(); }

    Point_2
    operator()(CGAL::Origin o) const
    { return Point_2(0,0); }

    Point_2
    operator()(const RT& x, const RT& y) const
    { return Point_2(x, y); }

    
    // We need this one, as such a functor is in the Filtered_kernel
    Point_2
    operator()(const RT& x, const RT& y, const RT& w) const
    { 
      if(w != 1){
	return Point_2(x/w, y/w); 
      } else {
	return Point_2(x,y);
      }
    }
  };


// Taken from include/CGAL/Cartesian/Cartesian_base.h

template < typename K_ >
struct MyCartesian_base : public CGAL::Cartesian_base< K_ >
{
    typedef K_                           Kernel;
    typedef MyPointC2<Kernel>            Point_2;
    typedef MySegmentC2<Kernel>          Segment_2;
    typedef MyConstruct_point_2<Kernel>  Construct_point_2;
};






// Taken from include/CGAL/Cartesian.h

template < typename FT_, typename Kernel >
struct MyCartesian_base_no_ref_count
  : public MyCartesian_base< Kernel >
{
    typedef FT_                                           RT;
    typedef FT_                                           FT;

    // The mecanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef CGAL::Simple_Handle_for<T>   type; };

    template < typename Kernel2 >
    struct Base { typedef MyCartesian_base_no_ref_count<FT_, Kernel2>  Type; };

    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}
};



template < typename FT_ >
struct MyKernel
  : public MyCartesian_base_no_ref_count<FT_, MyKernel<FT_> >
{
public:
  typedef MyCartesian_base_no_ref_count<FT_, MyKernel<FT_> > Kernel_base;
};



CGAL_ITERATOR_TRAITS_POINTER_SPEC_TEMPLATE(GF::MyKernel)



#endif // MYKERNEL_H
