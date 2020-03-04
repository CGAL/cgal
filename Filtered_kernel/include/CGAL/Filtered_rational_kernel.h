// Copyright (c) 2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Andreas Fabri

#ifndef CGAL_FILTERED_RATIONAL_KERNEL_H
#define CGAL_FILTERED_RATIONAL_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <CGAL/Filtered_rational.h>
#include <CGAL/Static_filtered_predicate.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_2.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_3.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Cartesian_converter.h>

#include <type_traits>

namespace CGAL {

template <typename T1, typename T2>
struct Approximate_exact_pair
  : public std::pair<T1,T2>
{
  typedef std::pair<T1,T2> Base;

  typedef T1 Approximate_type;
  
  Approximate_exact_pair()
  {}

  Approximate_exact_pair(const T1& t1, const T2& t2)
    : Base(t1,t2)
  {}

  Approximate_exact_pair(const std::pair<T1,T2>& p)
    : Base(p)
  {}
  
  const T1& approx() const
  {
    return this->first;
  }

  const T2& exact() const
  {
    return this->second;
  }
};

namespace mpl {

BOOST_MPL_HAS_XXX_TRAIT_DEF(first_type)
BOOST_MPL_HAS_XXX_TRAIT_DEF(second_type)

template <class T>
struct is_pair_
{
  enum { value = has_first_type<T>::value && has_second_type<T>::value };
};

} // mpl namespace

namespace frk {
// Small utility to manipulate pairs for kernel objects, and
// simple things for bool, Sign...  Object is yet another case...
template < typename T1, typename T2, typename EK, typename FRK >
struct Pairify
{
  typedef typename Type_mapper<T2,EK,FRK>::type result_type;
  result_type operator()(const T1 &t1, const T2 &t2) const
  {
    typedef typename Type_mapper<T2,EK,FRK>::type T;
    return T(std::make_pair(t1, t2)); }
};

template <typename EK, typename FRK>
struct Pairify <bool, bool, EK, FRK>
{
  typedef bool   result_type;
  result_type operator()(const bool &t1, const bool &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <Sign, Sign, EK, FRK>
{
  typedef Sign   result_type;
  result_type operator()(const Sign &t1, const Sign &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <Bounded_side, Bounded_side, EK, FRK>
{
  typedef Bounded_side   result_type;
  result_type operator()(const Bounded_side &t1, const Bounded_side &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <Angle, Angle, EK, FRK>
{
  typedef Angle   result_type;
  result_type operator()(const Angle &t1, const Angle &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <const Interval_nt<false>&,const Gmpq&, EK, FRK>
{
  typedef Filtered_rational<Interval_nt<false>,Gmpq> result_type;

  result_type operator()(const Interval_nt<false>& a, const Gmpq& e) const
  {
    return Filtered_rational<Interval_nt<false>, Gmpq>(a,e);
  }
};

} // namespace frk
  
#define SPEC_NT_GETTER(X) \
template <> \
struct Getter<X>\
{ \
  typedef X first_type; \
  typedef X second_type; \
}; \
inline const X& approx(const X& x) { return x; }\
inline const X& exact(const X& x) { return x; }

template <class A>
struct Getter
{
  typedef typename Getter<typename A::Rep>::first_type first_type;
  typedef typename Getter<typename A::Rep>::second_type second_type;
};

SPEC_NT_GETTER(int)
SPEC_NT_GETTER(unsigned int)
SPEC_NT_GETTER(double)
SPEC_NT_GETTER(Origin)
SPEC_NT_GETTER(Null_vector)
SPEC_NT_GETTER(Sign)
SPEC_NT_GETTER(Object)
SPEC_NT_GETTER(Return_base_tag)

#undef SPEC_NT_GETTER

template<>
struct Getter<Bbox_2>
{
    typedef Bbox_2 first_type;
    typedef Bbox_2 second_type;
};

template<>
struct Getter<Bbox_3>
{
    typedef Bbox_3 first_type;
    typedef Bbox_3 second_type;
};


template <class A1, class A2>
struct Getter<std::pair<A1, A2> >
{
  typedef A1 first_type;
  typedef A2 second_type;
};

template <class A1, class A2>
struct Getter<Approximate_exact_pair<A1, A2> >
{
  typedef A1 first_type;
  typedef A2 second_type;
};

template <class A_FT, class E_FT>
struct Getter<Filtered_rational<A_FT,E_FT>>
{
  typedef A_FT first_type;
  typedef E_FT second_type;
};

template <class A1, class A2>
const A1&
approx(const Approximate_exact_pair<A1,A2>& p)
{
  return p.first;
}

template <class A1, class A2>
const A1&
approx(const Filtered_rational<A1,A2>& p)
{
  return p.n1();
}

template <class A>
const typename A::Rep::first_type&
approx(const A& a, typename boost::disable_if< mpl::is_pair_<A> >::type* = nullptr)
{
  return approx(a.rep());
}
  
template <class A1, class A2>
const A2&
exact(const std::pair<A1,A2>& p)
{
  return p.second;
}
  
template <class A1, class A2>
const A2&
exact(const Filtered_rational<A1,A2>& p)
{
  return p.n2();
}

template <class A>
const typename A::Rep::second_type&
exact(const A& a, typename boost::disable_if< mpl::is_pair_<A> >::type* = nullptr)
{
  return exact(a.rep());
}

template <class AP, class EP>
class Filtered_rational_predicate
{
  AP ap;
  EP ep;

public:
  // CGAL_static_assertion((std::is_same<typename AP::result_type, typename EP::result_type>::value));
  typedef typename EP::result_type result_type;
  
public:
  Filtered_rational_predicate(const AP &pap = AP(), const EP &pep = EP()) : ap(pap), ep(pep) { }

  template <class ... A>
  typename CGAL::cpp11::result_of<EP(typename Getter<A>::second_type...)>::type
  operator()(const A&... a) const
  {
    typedef typename CGAL::cpp11::result_of<AP(typename Getter<A>::first_type...)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<EP(typename Getter<A>::second_type...)>::type result_type_2;

    CGAL::Interval_nt<false>::Protector p;

    try
    {
      result_type_1 res1 = ap(approx(a)...);
      if (is_certain(res1))
        return make_certain(res1);
    }
    catch(Uncertain_conversion_exception&) { }

    result_type_2 res2 = ep(exact(a)...);
    return res2;
  }
};

template <class AP, class EP, class AK, class EK, class FRK>
class Filtered_rational_construction
{
  AP ap;
  EP ep;

  CGAL::Cartesian_converter<EK, AK> e2a;

public:
  Filtered_rational_construction(const AP &pap = AP(), const EP &pep = EP()) : ap(pap), ep(pep) { }

  template <class T>
  struct result;

  template<typename F, typename ... A>
  struct result<F(A...)>
  {
    typedef typename cpp11::result_of<AP(typename Getter<A>::first_type...)>::type R1;
    typedef typename cpp11::result_of<EP(typename Getter<A>::second_type...)>::type R2;
    typedef typename frk::Pairify<R1,R2,EK,FRK>::result_type type;
  };

  // In case the exact result type is a reference
  template <typename AT, bool b = std::is_lvalue_reference<AT>::value>
  struct Approx
  {
    CGAL::Cartesian_converter<EK, AK> e2a;
    AP ap;

    Approx(const CGAL::Cartesian_converter<EK, AK>& e2a, const AP& ap) : e2a(e2a), ap(ap) { }

    template <typename ERT, typename ... A>
    const AT& operator()(const ERT&, const A&... a) const { return ap(a ...); }
  };

  // In case we have to generate the approximation from the result of the exact construction
  template <typename AT>
  struct Approx<AT, false>
  {
    CGAL::Cartesian_converter<EK, AK> e2a;
    AP ap;

    Approx(const CGAL::Cartesian_converter<EK, AK>& e2a, const AP& ap) : e2a(e2a), ap(ap) { }

    template <typename ERT, typename ... A>
    AT operator()(const ERT& ert, const A&... ) const { return e2a(ert); }
  };

  // TODO: I think the result_of is simply using AP::result_type because arguments are not valid (pairs...)
  template <class ... A>
  typename frk::Pairify<typename CGAL::cpp11::result_of<AP(typename Getter<A>::first_type...)>::type,
  typename CGAL::cpp11::result_of<EP(typename Getter<A>::second_type...)>::type, EK, FRK>::result_type
  operator()(const A&... a) const
  {
    typedef typename CGAL::cpp11::result_of<AP(typename Getter<A>::first_type...)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<EP(typename Getter<A>::second_type...)>::type result_type_2;
    result_type_2 res2 = ep(exact(a)...);

    return frk::Pairify<result_type_1, result_type_2, EK, FRK>()(Approx<result_type_1>(e2a, ap)
                                                                 (res2, approx(a)...), res2);
  }
};

namespace FRK {
template < typename K1, typename K2 >
struct Approx_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  const typename Getter<T>::first_type&
  operator()(const T&t) const
  { return approx(t); }

  
  const Null_vector&
  operator()(const Null_vector& n) const
  { return n; }
  
  const Bbox_2&
  operator()(const Bbox_2& b) const
  { return b; }
  
  const Bbox_3&
  operator()(const Bbox_3& b) const
  { return b; }
};

template < typename K1, typename K2 >
struct Exact_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  const typename Getter<T>::second_type&
  operator()(const T&t) const
  { return exact(t); }

  const Null_vector&
  operator()(const Null_vector& n) const
  { return n; }
  
  const Bbox_2&
  operator()(const Bbox_2& b) const
  { return b; }
  
  const Bbox_3&
  operator()(const Bbox_3& b) const
  { return b; }
};

}// namespace FRK

  
template < class AK, class EK, class Kernel_ >
class Filtered_rational_kernel_generic_base
{
protected:
  AK ak;
  EK ek;

public:

  typedef bool                      Boolean;
  typedef CGAL::Sign                Sign;
  typedef CGAL::Comparison_result   Comparison_result;
  typedef CGAL::Orientation         Orientation;
  typedef CGAL::Oriented_side       Oriented_side;
  typedef CGAL::Bounded_side        Bounded_side;
  typedef CGAL::Angle               Angle;

  typedef CGAL::Object Object_2;
  typedef CGAL::Object Object_3;

  typedef Kernel_ Kernel;
  typedef AK     Approximate_kernel;
  typedef EK     Exact_kernel;

  typedef FRK::Approx_converter<Kernel, Approximate_kernel>   C2F;
  typedef FRK::Exact_converter<Kernel, Exact_kernel>    C2E;

  template < typename T >
  struct Ambient_dimension {
    typedef typename T::Ambient_dimension type;
  };

  typedef typename EK::Kernel_tag                       Kernel_tag;
  typedef typename EK::Rep_tag                          Rep_tag;

  enum { Has_filtered_predicates = true };
  enum { Has_static_filters = false };
  typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

  typedef Filtered_rational<typename AK::FT, typename EK::FT> FT;
  typedef FT RT;

  typedef Cartesian_coordinate_iterator_2<Kernel> Cartesian_const_iterator_2;
  typedef Cartesian_coordinate_iterator_3<Kernel> Cartesian_const_iterator_3;

  typedef CGAL::Aff_transformationC2<Kernel_> Aff_transformation_2;
  typedef CGAL::Aff_transformationC3<Kernel_> Aff_transformation_3;

  // Kernel objects are defined as pairs, with primitives run in parallel.
#define CGAL_frk_pair(X) typedef Approximate_exact_pair<typename AK::X, typename EK::X> X;


  // TODO : Object_[23] are subtil : should probably be Object(pair<...>).
  // Or should Assign_[23] be used, and that's it ?
  // In any case, Assign will have to be treated separately because it
  // takes its first argument by non-const reference.
  // Maybe Primitive_checker should provide a variant with non-const ref...

  CGAL_frk_pair(Point_2)
  CGAL_frk_pair(Weighted_point_2)
  CGAL_frk_pair(Vector_2)
  CGAL_frk_pair(Direction_2)
  CGAL_frk_pair(Line_2)
  CGAL_frk_pair(Ray_2)
  CGAL_frk_pair(Segment_2)
  CGAL_frk_pair(Triangle_2)
  CGAL_frk_pair(Iso_rectangle_2)
  CGAL_frk_pair(Circle_2)
  CGAL_frk_pair(Conic_2)


  CGAL_frk_pair(Point_3)
  CGAL_frk_pair(Weighted_point_3)
  CGAL_frk_pair(Plane_3)
  CGAL_frk_pair(Vector_3)
  CGAL_frk_pair(Direction_3)
  CGAL_frk_pair(Line_3)
  CGAL_frk_pair(Ray_3)
  CGAL_frk_pair(Segment_3)
  CGAL_frk_pair(Triangle_3)
  CGAL_frk_pair(Tetrahedron_3)
  CGAL_frk_pair(Iso_cuboid_3)
  CGAL_frk_pair(Circle_3)
  CGAL_frk_pair(Sphere_3)

#undef CGAL_frk_pair

#define CGAL_Kernel_pred(P, Pf) \
  typedef Filtered_rational_predicate<typename AK::P, typename EK::P> P; \
  P Pf() const { return P(ak.Pf(), ek.Pf()); }

#define CGAL_Kernel_cons(C, Cf) \
  typedef Filtered_rational_construction<typename AK::C, typename EK::C, AK, EK, Kernel_> C; \
  C Cf() const { return C(ak.Cf(), ek.Cf()); }


  public:

  #include <CGAL/Kernel/interface_macros.h>
};

template < class AK, class EK, class Kernel_ >
class Filtered_rational_kernel_base
    : public Filtered_rational_kernel_generic_base<AK,EK,Kernel_>
{
public:
  typedef Filtered_rational_kernel_base<AK,EK,Kernel_> Self;
  using typename Filtered_rational_kernel_generic_base<AK,EK,Kernel_>::Cartesian_const_iterator_2;
  using typename Filtered_rational_kernel_generic_base<AK,EK,Kernel_>::Cartesian_const_iterator_3;

  class Construct_object_2
  {
    typedef typename Kernel_::Object_2   Object_2;
  public:
    typedef Object_2         result_type;

    template <class Cls>
    Object_2
    operator()( const Cls& c) const
    { return make_object(c); }
  };

  Construct_object_2 construct_object_2_object() const
  {
    return Construct_object_2();
  }

  class Construct_object_3
  {
    typedef typename Kernel_::Object_3   Object_3;
  public:
    typedef Object_3         result_type;

    template <class Cls>
    Object_3
    operator()( const Cls& c) const
    { return make_object(c); }
  };


  Construct_object_3 construct_object_3_object() const
  {
    return Construct_object_3();
  }


  class Assign_2
  {
    typedef typename Kernel_::Object_2  Object_2;
  public:
    typedef bool                  result_type;

    template <class T>
    result_type
    operator()(T& t, const Object_2& o) const
    { return assign(t, o); }
  };

  Assign_2 assign_2_object() const
  {
    return Assign_2();
  }

  class Assign_3
  {
    typedef typename Kernel_::Object_3        Object_3;
  public:
    typedef bool                        result_type;

    template <class T>
    result_type
    operator()(T& t, const Object_3& o) const
    { return assign(t, o); }
  };

  Assign_3 assign_3_object() const
  {
    return Assign_3();
  }

  class Construct_bbox_2 {
  public:
    typedef Bbox_2 result_type;

    template <typename T1, typename T2>
    Bbox_2 operator()(const std::pair<T1,T2>& p) const
    {
      // AF: Or do we want to construct the BBox from Gmpq ?
      return typename AK::Construct_bbox_2()(p.first);
    }

  };

  Construct_bbox_2 construct_bbox_2_object() const
  {
    return Construct_bbox_2();
  }

  class Construct_bbox_3 {
  public:
    typedef Bbox_2 result_type;

    template <typename T1, typename T2>
    Bbox_3 operator()(const std::pair<T1,T2>& p) const
    {
      // AF: Or do we want to construct the BBox from Gmpq ?
      return typename AK::Construct_bbox_3()(p.first);
    }
  };

  Construct_bbox_3 construct_bbox_3_object() const
  {
    return Construct_bbox_3();
  }

  template <typename OptionalVariant>
  class Make_optional_variant
    : public boost::static_visitor<>
  {
    OptionalVariant& ov;
  public:

    Make_optional_variant(OptionalVariant& ov)
      : ov(ov)
    {}

    template <typename ET>
    void operator()( const ET & et ) const
    {
      typedef typename Type_mapper<ET, EK, Kernel_>::type T;
      Cartesian_converter<EK,AK> e2a;
      e2a(et);
      ov = OptionalVariant(T(std::make_pair(e2a(et),et)));
    }

    template <typename ET>
    void operator()( const std::vector<ET>& vec) const
    {
      typedef typename Type_mapper<ET, EK, Kernel_>::type T;
      std::vector<T> resvec;
      Cartesian_converter<EK,AK> e2a;
      for(const ET& et : vec){
        resvec.push_back(T(std::make_pair(e2a(et),et)));
      }
      ov = OptionalVariant(resvec);
    }
  };


  class Intersect_2 {
  public:

    template<typename>
    struct result;

    template<typename F, typename A, typename B>
    struct result<F(A,B)> {
      typedef typename Intersection_traits<Kernel_, A, B>::result_type type;
    };

    template <typename T1, typename T2>
    typename Intersection_traits<Kernel_,T1,T2>::result_type
    operator()(const T1& s1, const T2& s2) const
    {
      typedef typename Type_mapper<T1,Kernel_,EK>::type EKT1;
      typedef typename Type_mapper<T2,Kernel_,EK>::type EKT2;

      typedef  typename Intersection_traits<EK, EKT1, EKT2>::result_type Exact_optional_variant;
      typedef typename Exact_optional_variant::value_type Exact_variant;

      Exact_optional_variant  eres = typename EK::Intersect_2()(s1.second,s2.second);

      if(! eres){
        return boost::none;
      }
      Exact_variant ev = *eres;

      typedef typename Intersection_traits<Kernel_,T1,T2>::result_type result_type;
      result_type res;
      boost::apply_visitor( Make_optional_variant<result_type>(res), ev );

      return res;
    }

  };

  Intersect_2 intersect_2_object() const
  {
    return Intersect_2();
  }

  class Intersect_3 {
  public:

    template<typename>
    struct result;

    template<typename F, typename A, typename B>
    struct result<F(A,B)> {
      typedef typename Intersection_traits<Kernel_, A, B>::result_type type;
    };

    template <typename T1, typename T2>
    typename Intersection_traits<Kernel_,T1,T2>::result_type
    operator()(const T1& s1, const T2& s2) const
    {
      typedef typename Type_mapper<T1,Kernel_,EK>::type EKT1;
      typedef typename Type_mapper<T2,Kernel_,EK>::type EKT2;

      typedef  typename Intersection_traits<EK, EKT1, EKT2>::result_type Exact_optional_variant;
      typedef typename Exact_optional_variant::value_type Exact_variant;

      Exact_optional_variant  eres = typename EK::Intersect_3()(s1.second,s2.second);

      if(! eres){
        return boost::none;
      }
      Exact_variant ev = *eres;

      typedef typename Intersection_traits<Kernel_,T1,T2>::result_type result_type;
      result_type res;
      boost::apply_visitor( Make_optional_variant<result_type>(res), ev );

      return res;
    }

  };

  Intersect_3 intersect_3_object() const
  {
    return Intersect_3();
  }

  class Construct_cartesian_const_iterator_2 {
  public:

    typedef Cartesian_const_iterator_2 result_type;
    
    template <typename PV>
    Cartesian_const_iterator_2 operator()(const PV& pv) const
    {
      return Cartesian_const_iterator_2(&pv);
    }

    template <typename PV>
    Cartesian_const_iterator_2 operator()(const PV& pv, int i) const
    {
      return Cartesian_const_iterator_2(&pv, i);
    }
  };

  Construct_cartesian_const_iterator_2 construct_cartesian_const_iterator_2_object() const
  {
    return Construct_cartesian_const_iterator_2();
  }

  class Construct_cartesian_const_iterator_3 {
  public:
    
    typedef Cartesian_const_iterator_3 result_type;
    
    template <typename PV>
    Cartesian_const_iterator_3 operator()(const PV& pv) const
    {
      return Cartesian_const_iterator_3(&pv);
    }

    template <typename PV>
    Cartesian_const_iterator_3 operator()(const PV& pv, int i) const
    {
      return Cartesian_const_iterator_3(&pv, i);
    }
  };

  Construct_cartesian_const_iterator_3 construct_cartesian_const_iterator_3_object() const
  {
    return Construct_cartesian_const_iterator_3();
  }
};

template < class AK, class EK >
class Filtered_rational_kernel_without_type_equality
  : public Filtered_rational_kernel_base<AK,EK,Filtered_rational_kernel_without_type_equality<AK,EK>>
{};
  
template < class AK, class EK >
class Filtered_rational_kernel
  : public Static_filters_base_without_filtered_wrapping<
             Type_equality_wrapper<Filtered_rational_kernel_base<AK,EK, Filtered_rational_kernel<AK,EK>>,
                                   Filtered_rational_kernel<AK,EK> > >
{
};
  
} //namespace CGAL

#endif // CGAL_FILTERED_RATIONAL_KERNEL_H
