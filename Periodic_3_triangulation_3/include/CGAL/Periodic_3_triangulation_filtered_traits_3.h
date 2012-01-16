// Copyright (c) 2004,2006-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>


#ifndef CGAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H

#include <string>
#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

namespace CGAL {

// This template class is a wrapper that implements the filtering for any
// predicate (dynamic filters with IA).

// TODO :
// - each predicate in the default kernel should define a tag that says if it
//   wants to be filtered or not (=> all homogeneous predicate define this
//   tag).  We could even test-suite that automatically.  It makes a strong
//   new requirement on the kernel though...
//   Could be done with a traits mechanism ?
//   A default template could use the current IA, but other tags or whatever
//   could specify no filtering at all, or static filtering...
// - same thing for constructions => virtual operator() ?
// - similarly, constructions should have a tag saying if they can throw or
//   not, or we let all this up to the compiler optimizer to figure out ?
// - Some caching could be done at the Point_2 level.


template <class EP, class AP, class C2E, class C2A, bool Protection = true>
class Filtered_periodic_predicate
{
  EP  ep;
  AP  ap;
  C2E c2e;
  C2A c2a;

  typedef typename AP::result_type  Ares;

protected:
  typename AP::Iso_cuboid_3 * _domain;

public:

  typedef AP    Approximate_predicate;
  typedef EP    Exact_predicate;
  typedef C2E   To_exact_converter;
  typedef C2A   To_approximate_converter;

  typedef typename EP::result_type  result_type;
  // AP::result_type must be convertible to EP::result_type.

  Filtered_periodic_predicate()
  {}

  // These constructors are used for constructive predicates.
  // You should try to avoid constructive predicates, as they will construct
  // the exact values systematically (in the ctor), rather than lazily.
  template <class OE, class OA> 
  Filtered_periodic_predicate(const OE * oe, const OA * oa) 
    : ep(oe), ap(oa) 
  {} 

  //  template <class O>
  //Filtered_periodic_predicate(const O * o)
  //  : ep(&c2e(*o)), ap(&c2a(*o)) {}

  //  template <class O1, class O2>
  //Filtered_periodic_predicate(const O1 &o1, const O2 &o2)
  //  : ep(c2e(o1), c2e(o2)), ap(c2a(o1), c2a(o2))
  //{}

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <typename... Args>
  result_type
  operator()(const Args&... args) const;
#else

  template <class A1>
  result_type
  operator()(const A1 &a1) const;

  template <class A1, class A2>
  result_type
  operator()(const A1 &a1, const A2 &a2) const;

  template <class A1, class A2, class A3>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const;

  template <class A1, class A2, class A3, class A4>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const;

  template <class A1, class A2, class A3, class A4, class A5>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9, class A10>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9, const A10 &a10) const;

  // Idem for more than 10 arguments.  Do it on demand.

#endif
};

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <typename... Args>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const Args&... args) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(args)...);
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(args)...);
}

#else

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  
	  Ares res = ap(c2a(a1));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7), c2a(a8));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7),
              c2e(a8));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7), c2a(a8), c2a(a9));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7),
              c2e(a8), c2e(a9));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9, class A10>
typename Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_periodic_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9, const A10 &a10) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7), c2a(a8), c2a(a9), c2a(a10));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7),
              c2e(a8), c2e(a9), c2e(a10));
}

#endif

// The Offset_converter is parametrized by a usual kernel converter,
// and adds the conversions for Offsets.
template < typename Converter >
struct Offset_converter_3
  : public Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Periodic_3_triangulation_traits_base_3<Source_kernel>
                   ::Offset  Source_off;
  typedef typename Periodic_3_triangulation_traits_base_3<Source_kernel>
                   ::Point_3  Source_pt;

  typedef typename Periodic_3_triangulation_traits_base_3<Target_kernel>
                   ::Offset  Target_off;
  typedef typename Periodic_3_triangulation_traits_base_3<Target_kernel>
                   ::Point_3  Target_pt;


  using Converter::operator();

  Target_off
  operator()(const Source_off &off) const
  {
    return off;
  }
};

// The argument is supposed to be a Filtered_kernel like kernel.
template < typename K, typename Off >
class Periodic_3_triangulation_filtered_traits_base_3
  : public Periodic_3_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_triangulation_traits_base_3<K, Off> Base;

  // Exact traits is based on the exact kernel.
  typedef Periodic_3_triangulation_traits_3<typename K::Exact_kernel,
                                            Off>
                                                   Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_3_triangulation_traits_3<typename K::Approximate_kernel,
                                            Off>
                                                   Filtering_traits;
private:
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef typename C2E::Target_kernel::Iso_cuboid_3 Exact_iso_cuboid_3;
  typedef typename C2F::Target_kernel::Iso_cuboid_3 Approximate_iso_cuboid_3;
 
public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  void set_domain(const Iso_cuboid_3& domain) {
    C2E c2e;
    C2F c2f;
    this->_domain = domain;
    this->_domain_e = c2e(this->_domain);
    this->_domain_f = c2f(this->_domain);
  }

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Compare_xyz_3,
            typename Filtering_traits::Compare_xyz_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Compare_xyz_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Coplanar_orientation_3,
            typename Filtering_traits::Coplanar_orientation_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Coplanar_orientation_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Orientation_3,
            typename Filtering_traits::Orientation_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Orientation_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Coplanar_side_of_bounded_circle_3,
            typename Filtering_traits::Coplanar_side_of_bounded_circle_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Coplanar_side_of_bounded_circle_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Side_of_oriented_sphere_3,
            typename Filtering_traits::Side_of_oriented_sphere_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Side_of_oriented_sphere_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Compare_distance_3,
            typename Filtering_traits::Compare_distance_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Compare_distance_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Side_of_bounded_sphere_3,
            typename Filtering_traits::Side_of_bounded_sphere_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Side_of_bounded_sphere_3;


  Compare_xyz_3 compare_xyz_3_object() const
  { return Compare_xyz_3(&_domain_e,&_domain_f);}

  Coplanar_orientation_3 coplanar_orientation_3_object() const
  { return Coplanar_orientation_3(&_domain_e,&_domain_f); }

  Orientation_3 orientation_3_object() const
  { return Orientation_3(&_domain_e,&_domain_f);}

  Coplanar_side_of_bounded_circle_3
  coplanar_side_of_bounded_circle_3_object() const 
  { return Coplanar_side_of_bounded_circle_3(&_domain_e,&_domain_f); }

  Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const
  { return Side_of_oriented_sphere_3(&_domain_e,&_domain_f);}

  Compare_distance_3 compare_distance_3_object() const
  { return Compare_distance_3(&_domain_e,&_domain_f);}

  Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object() const
  { return Side_of_bounded_sphere_3(&_domain_e,&_domain_f);}

  // The following are inherited since they are constructions :
  // Construct_segment_3
  // Construct_triangle_3
  // Construct_tetrahedron_3
  // Construct_circumcenter_3

 protected:
  Exact_iso_cuboid_3 _domain_e;
  Approximate_iso_cuboid_3 _domain_f;
};

} //namespace CGAL

#include <CGAL/Periodic_3_triangulation_statically_filtered_traits_3.h>

namespace CGAL {

template < typename K, typename Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_triangulation_filtered_traits_3
  : public Periodic_3_triangulation_statically_filtered_traits_3<
  Periodic_3_triangulation_filtered_traits_base_3<K, Off> > {
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H
