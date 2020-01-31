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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_FILTERED_RATIONAL_KERNEL_H
#define CGAL_FILTERED_RATIONAL_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Filtered_rational.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_2.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_3.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Cartesian_converter.h>

namespace CGAL {

  namespace mpl
{

BOOST_MPL_HAS_XXX_TRAIT_DEF(first_type)
BOOST_MPL_HAS_XXX_TRAIT_DEF(second_type)

template <class T> struct is_pair_ {
  enum { value = has_first_type<T>::value && has_second_type<T>::value }; 
};

} // mpl namespace
  
// Small utility to manipulate pairs for kernel objects, and
// simple things for bool, Sign...  Object is yet another case...
  template < typename T1, typename T2, typename EK, typename FRK >
struct Pairify {
  typedef typename Type_mapper<T2,EK,FRK>::type result_type;
  result_type operator()(const T1 &t1, const T2 &t2) const
  {
    typedef typename Type_mapper<T2,EK,FRK>::type T;
    return T(std::make_pair(t1, t2)); }
};

template <typename EK, typename FRK>
struct Pairify <bool, bool, EK, FRK> {
  typedef bool   result_type;
  result_type operator()(const bool &t1, const bool &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <Sign, Sign, EK, FRK> {
  typedef Sign   result_type;
  result_type operator()(const Sign &t1, const Sign &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <Bounded_side, Bounded_side, EK, FRK> {
  typedef Bounded_side   result_type;
  result_type operator()(const Bounded_side &t1, const Bounded_side &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <Angle, Angle, EK, FRK> {
  typedef Angle   result_type;
  result_type operator()(const Angle &t1, const Angle &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <typename EK, typename FRK>
struct Pairify <const Interval_nt<false>&,const Gmpq&, EK, FRK> {
  
  typedef Filtered_rational<Interval_nt<false>,Gmpq> result_type;
  
  result_type operator()(const Interval_nt<false>& a, const Gmpq& e) const
  {
    return Filtered_rational<Interval_nt<false>, Gmpq>(a,e);
  }
};

#define SPEC_NT_GETTER(X) \
template <> \
struct Getter<X>\
{ \
  typedef X first_type; \
  typedef X second_type; \
}; \
const X& get_first(const X& x) { return x; }\
const X& get_second(const X& x) { return x; }


template <class A>
struct Getter;

SPEC_NT_GETTER(int)
SPEC_NT_GETTER(unsigned int)
SPEC_NT_GETTER(double)
SPEC_NT_GETTER(Origin)
SPEC_NT_GETTER(Null_vector)
SPEC_NT_GETTER(Sign)
SPEC_NT_GETTER(Object)
SPEC_NT_GETTER(Return_base_tag)

#undef SPEC_NT_GETTER

template <class A1, class A2>
struct Getter<std::pair<A1,A2>>
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

template <class A>
struct Getter
  : Getter<typename A::Rep>
{};

template <class A1, class A2>
const A1&
get_first(const std::pair<A1,A2>& p)
{
  return p.first;
}

template <class A1, class A2>
const A1&
get_first(const Filtered_rational<A1,A2>& p)
{
  return p.n1();
}

template <class A>
const typename A::Rep::first_type&
get_first(const A& a, typename boost::disable_if< mpl::is_pair_<A> >::type* = nullptr)
{
  return get_first(a.rep());
}

template <class A1, class A2>
const A2&
get_second(const std::pair<A1,A2>& p)
{
  return p.second;
}

template <class A1, class A2>
const A2&
get_second(const Filtered_rational<A1,A2>& p)
{
  return p.n2();
}
  
template <class A>
const typename A::Rep::second_type&
get_second(const A& a, typename boost::disable_if< mpl::is_pair_<A> >::type* = nullptr)
{
  return get_second(a.rep());
}

template <class P1, class P2>
class Predicate_wrapper
{
  P1  p1;
  P2  p2;

public:
  Predicate_wrapper(const P1 &pp1 = P1(), const P2 &pp2 = P2())
    : p1(pp1), p2(pp2)
  { }

  template <class ... A>
  typename CGAL::cpp11::result_of<P2(const typename Getter<A>::second_type&...)>::type
  operator()(const A&... a) const
  {
    typedef typename CGAL::cpp11::result_of<P1(const typename Getter<A>::first_type&...)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(const typename Getter<A>::second_type&...)>::type result_type_2;
    CGAL::Interval_nt<false>::Protector p;
    try{
      result_type_1 res1 = p1(get_first(a)...);
      if (is_certain(res1))
        return make_certain(res1);
    }
    catch(Uncertain_conversion_exception&)
    {}
    result_type_2 res2 = p2(get_second(a)...);
    return res2;
  }
};

template <class P1, class P2, class AK, class EK, class FRK>
class Construction_wrapper
{
  P1  p1;
  P2  p2;

  CGAL::Cartesian_converter<EK, AK> to_k1;

public:
  Construction_wrapper(const P1 &pp1 = P1(), const P2 &pp2 = P2())
    : p1(pp1), p2(pp2)
  { }

  template <class T>
  struct result;

  template<typename F, typename ... A>
  struct result<F(A...)> {
    typedef typename cpp11::result_of<P1(const typename Getter<A>::first_type&...)>::type R1;
    typedef typename cpp11::result_of<P2(const typename Getter<A>::second_type&...)>::type R2;
    typedef typename Pairify<R1,R2,EK,FRK>::result_type type;
  };
  
  // TODO: I think the result_of is simply using P1::result_type because arguments are not valid (pairs...)
  template <class ... A>
  typename Pairify<typename CGAL::cpp11::result_of<P1(const typename Getter<A>::first_type&...)>::type,
                   typename CGAL::cpp11::result_of<P2(const typename Getter<A>::second_type&...)>::type,
                   EK,FRK>::result_type
  operator()(const A&... a) const
  {
    typedef typename CGAL::cpp11::result_of<P1(const typename Getter<A>::first_type&...)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(const typename Getter<A>::second_type&...)>::type result_type_2;
    result_type_2 res2 = p2(get_second(a)...);
    return Pairify<result_type_1, result_type_2,EK,FRK>()(to_k1(res2), res2);
  }


  // this is the overload for functors such as Construct_vertex_2
  template <typename AT, typename ET>
  typename Pairify<typename CGAL::cpp11::result_of<P1(const AT&,int)>::type,
                   typename CGAL::cpp11::result_of<P2(const ET&,int)>::type,
                   EK,FRK>::result_type
  operator()(const std::pair<AT,ET>& p, int i) const
  {
    typedef typename CGAL::cpp11::result_of<P1(const AT&,int)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(const ET&,int)>::type result_type_2;
    result_type_2 res2 = p2(get_second(p),i);
    typedef typename Type_mapper<ET,EK,FRK>::type T;
    return Pairify<result_type_1,result_type_2,EK,FRK>()(to_k1(res2), res2);
  }
};


template < class AK, class EK, class Kernel_ >
class Filtered_rational_kernel_generic_base
{
protected:
  AK k1;
  EK k2;

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
  typedef AK     Kernel1;
  typedef EK     Kernel2;

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
#define CGAL_kc_pair(X) typedef std::pair<typename AK::X, typename EK::X> X;


  // TODO : Object_[23] are subtil : should probably be Object(pair<...>).
  // Or should Assign_[23] be used, and that's it ?
  // In any case, Assign will have to be treated separately because it
  // takes its first argument by non-const reference.
  // Maybe Primitive_checker should provide a variant with non-const ref...

  //CGAL_kc_pair(Object_2)
  //CGAL_kc_pair(Object_3)

  CGAL_kc_pair(Point_2)
  CGAL_kc_pair(Weighted_point_2)
  CGAL_kc_pair(Vector_2)
  CGAL_kc_pair(Direction_2)
  CGAL_kc_pair(Line_2)
  CGAL_kc_pair(Ray_2)
  CGAL_kc_pair(Segment_2)
  CGAL_kc_pair(Triangle_2)
  CGAL_kc_pair(Iso_rectangle_2)
  CGAL_kc_pair(Circle_2)
  CGAL_kc_pair(Conic_2)


  CGAL_kc_pair(Point_3)
  CGAL_kc_pair(Weighted_point_3)
  CGAL_kc_pair(Plane_3)
  CGAL_kc_pair(Vector_3)
  CGAL_kc_pair(Direction_3)
  CGAL_kc_pair(Line_3)
  CGAL_kc_pair(Ray_3)
  CGAL_kc_pair(Segment_3)
  CGAL_kc_pair(Triangle_3)
  CGAL_kc_pair(Tetrahedron_3)
  CGAL_kc_pair(Iso_cuboid_3)
  CGAL_kc_pair(Circle_3)
  CGAL_kc_pair(Sphere_3)

#undef CGAL_kc_pair

#define CGAL_Kernel_pred(X, Y) \
  typedef Predicate_wrapper<typename AK::X, typename EK::X> X; \
  X Y() const { return X(k1.Y(), k2.Y()); }

#define CGAL_Kernel_cons(X, Y) \
  typedef Construction_wrapper<typename AK::X, typename EK::X, AK, EK, Kernel_> X; \
  X Y() const { return X(k1.Y(), k2.Y()); }


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
      Cartesian_converter<EK,AK> to_k1;
      to_k1(et);
      ov = OptionalVariant(T(std::make_pair(to_k1(et),et)));
    }

    template <typename ET>
    void operator()( const std::vector<ET>& vec) const
    {
      typedef typename Type_mapper<ET, EK, Kernel_>::type T;
      std::vector<T> resvec;
      Cartesian_converter<EK,AK> to_k1;
      for(const ET& et : vec){
        resvec.push_back(T(std::make_pair(to_k1(et),et)));
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
  : public Type_equality_wrapper<Filtered_rational_kernel_base<AK,EK, Filtered_rational_kernel<AK,EK>>,
                                 Filtered_rational_kernel<AK,EK> >
{
};
  
} //namespace CGAL

#endif // CGAL_FILTERED_RATIONAL_KERNEL_H
