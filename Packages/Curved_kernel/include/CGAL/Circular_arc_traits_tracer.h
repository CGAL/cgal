// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Circular_arc_traits_tracer.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_TRACER_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_TRACER_H

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Profile_counter.h>

namespace CGAL {

/// Wrapper around a traits class for Arrangement which just forwards the
/// calls to its nested traits, and prints debug info for each argument and
/// return value.  Compiling with -DCGAL_PROFILE dumps the number of calls
/// to each predicate at the end of the program.

// to get better profiles for the traits interface :
#ifndef CGAL_NO_INLINE
#  define CGAL_NO_INLINE __attribute__((__noinline__))
//#define CGAL_NO_INLINE
#endif

// FIXME : should it be the responsibility of the tracer
//         to add the Debug_id stuff ?
template < typename Traits, int debug_level = 2 >
class Circular_arc_traits_tracer
 : public Traits
{

  typedef Circular_arc_traits_tracer<Traits, debug_level> Self;

// debug_level :
// 0 = no debug messages
// 1 = important messages
// 2 = maybe important messages
// 3 = all calls

  //  Traits traits;

  mutable int nested_level; /// for detection of nested calls
  mutable int call_number;  /// identifies the call number

// FIXME : transform these macros into member functions ?
#define CGAL_DEBUG(level, args) if (level <= debug_level) \
   { std::cout << "Call #" << t->call_number++ << " " << __PRETTY_FUNCTION__ \
               << "  ( " << args << " )"; \
     if (t->nested_level != 0) \
       std::cout << "   [ nested level " << t->nested_level << " ]"; \
     ++t->nested_level; \
     std::cout << std::endl; }

#define CGAL_DEBUG_RET(level, args) if (level <= debug_level) \
   { std::cout << "         returns : " << args; \
     --t->nested_level; \
     if (t->nested_level != 0) \
       std::cout << "   [ nested level " << t->nested_level << " ]"; \
     std::cout << std::endl; }

public:

  typedef typename Traits::Kernel                Kernel;
  typedef typename Traits::Curve_2               Curve_2;
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;

  typedef typename Traits::Point                 Point;
  typedef typename Traits::Point_2               Point_2;

  typedef typename Kernel::Circle_2                 Circle;
  typedef typename Kernel::Circular_arc_2           Circular_arc;

  typedef typename Traits::Has_left_category     Has_left_category;
  typedef typename Traits::Has_merge_category    Has_merge_category;
 
  const Traits & traits () const
    { return static_cast<const Traits &>(*this); }

  Circular_arc_traits_tracer(const Traits &t = Traits())
    : Traits(t), nested_level(0), call_number(0) {}

  struct Compare_x_2 {

    Compare_x_2(const Self *tt) :t(tt) {}

    const Self *t;

    CGAL::Comparison_result
    CGAL_NO_INLINE
    operator()(const Point_2 &p, const Point_2 &q) const
    {
      CGAL_PROFILER(__FUNCTION__);
      CGAL_DEBUG(2, p.id() << ", " << q.id());
      CGAL::Comparison_result ret = t->traits().compare_x_2_object()(p, q);
      CGAL_DEBUG_RET(2, ret);
      return ret;
    }
  };

  Compare_x_2
  compare_x_2_object() const
  { return Compare_x_2(this); }

  struct Compare_y_2 {

    Compare_y_2(const Self *tt) :t(tt) {}

    const Self *t;

    CGAL::Comparison_result
    CGAL_NO_INLINE
    operator()(const Point_2 &p, const Point_2 &q) const
    {
      CGAL_PROFILER(__FUNCTION__);
      CGAL_DEBUG(2, p.id() << ", " << q.id());
      CGAL::Comparison_result ret = t->traits().compare_y_2_object()(p, q);
      CGAL_DEBUG_RET(2, ret);
      return ret;
    }
  };

  Compare_y_2
  compare_y_2_object() const
  { return Compare_y_2(this); }

  struct Equal_2 {

    Equal_2(const Self *tt) :t(tt) {}

    const Self *t;

    bool
    CGAL_NO_INLINE
    operator()(const Point_2 &p, const Point_2 &q) const
    {
      CGAL_PROFILER(__FUNCTION__);
      CGAL_DEBUG(2, p.id() << ", " << q.id());
      bool ret = t->traits().equal_2_object()(p, q);
      CGAL_DEBUG_RET(2, ret);
      return ret;
    }

    bool
    CGAL_NO_INLINE
    operator()(const Circular_arc &a1, const Circular_arc &a2) const
    {
      CGAL_PROFILER(__FUNCTION__);
      CGAL_DEBUG(2, a1.id() << ", " << a2.id());
      bool ret = t->traits().equal_2_object()(a1, a2);
      CGAL_DEBUG_RET(2, ret);
      return ret;
    }
  };

  Equal_2
  equal_2_object() const
  { return Equal_2(this); }

  struct Intersect_2 {

    Intersect_2(const Self *tt) :t(tt) {}

    const Self *t;

    template < class OutputIterator >
    OutputIterator
    CGAL_NO_INLINE
    operator()(const Circle & c1, const Circle & c2, 
	       OutputIterator res) const
    {
      CGAL_PROFILER(__FUNCTION__);
      CGAL_DEBUG(2, c1.id() << ", " << c2.id());
      std::vector<CGAL::Object> ret;
      t->traits().intersect_2_object()(c1, c2,
                                                   std::back_inserter(ret));
      
      if (2 <= debug_level)
	{
	  std::cout << "         returns : ";
          for(int i=0; i<ret.size(); ++i) {
	    //Point_2 p;
	    //Circular_arc c;
	    //std::pair<Point_2, int> pp;
	    if (const Point_2 *p = CGAL::object_cast<Point_2>(&ret[i])) {
	      std::cout << *p << "  ";
	    }
	    else if (const Circular_arc *c = CGAL::object_cast<Circular_arc>(&ret[i])) {
	      std::cout << *c << "  ";
	    }
	    else if (const std::pair<Point_2, int> *pp = CGAL::object_cast<std::pair<Point_2, int> >(&ret[i])) {
	      std::cout << " pair ( " << pp->first << " ; " <<
		pp->second << " )  ";
	    }
	    else {
	      std::cout << " [ something else ? ] ";
	    }
	  }
	  --nested_level;
	  if (nested_level != 0)
	    std::cout << "   [ nested level " << nested_level << " ]";
	  std::cout << std::endl;
	}
      
      return std::copy(ret.begin(), ret.end(), res);
    }

    template < class OutputIterator >
    OutputIterator
    CGAL_NO_INLINE
    operator()(const Circular_arc & c1, const Circular_arc & c2, 
	       OutputIterator res) const
    {
      CGAL_PROFILER(__FUNCTION__);
      CGAL_DEBUG(2, c1.id() << ", " << c2.id());
      std::vector<CGAL::Object> ret;
      t->traits().intersect_2_object()(c1, c2, std::back_inserter(ret));
      
      if (2 <= debug_level)
	{
	  std::cout << "         returns : ";
          for(unsigned i=0; i<ret.size(); ++i) {
	   // Point_2 p;
	   // Circular_arc c;
	   // std::pair<Point_2, unsigned> pp;
	    if (const Point_2 *p = CGAL::object_cast<Point_2>(&ret[i])) {
	      std::cout << *p << "  ";
	    }
	    else if (const Circular_arc *c = CGAL::object_cast<Circular_arc>(&ret[i])) {
	      std::cout << *c << "  ";
	    }
	    else if (const std::pair<Point_2, unsigned> *pp = CGAL::object_cast<std::pair<Point_2, unsigned> > (&ret[i])) {
	      std::cout << " pair ( " << pp->first << " ; " <<
		pp->second << " )  ";
	    }
	    else {
	      std::cout << " [ something else ? ] ";
	    }
	  }
	  --t->nested_level;
	  if (t->nested_level != 0)
	    std::cout << "   [ nested level " << t->nested_level << " ]";
	  std::cout << std::endl;
	}
      
      return std::copy(ret.begin(), ret.end(), res);
    }
  };

  Intersect_2
  intersect_2_object() const
  { return Intersect_2(this); }

};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_TRACER_H
