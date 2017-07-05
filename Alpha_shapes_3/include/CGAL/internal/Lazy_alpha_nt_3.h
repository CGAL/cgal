// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sébastien Loriot <sebastien.loriot@geometryfactory.com>
//                 Mael Rouxel-Labbé

#ifndef CGAL_INTERNAL_LAZY_ALPHA_NT_3_H
#define CGAL_INTERNAL_LAZY_ALPHA_NT_3_H

#include <CGAL/license/Alpha_shapes_3.h>

#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/Has_conversion.h>

#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/optional.hpp>

#include <iostream>

namespace CGAL {

namespace internal{

// check whether Cartesian_converter can do the following conversions
//  -- Input_traits::(Weighted_)Point_3 to K2::(Weighted_)Point_3
//  -- Input_traits::(Weighted_)Point_3 to K3::(Weighted_)Point_3
//
template < class Input_traits, class Kernel_approx, class Kernel_exact,
           class Weighted_tag >
class Is_traits_point_convertible
{
  typedef typename Kernel_traits<typename Input_traits::Point_3>::Kernel   Kernel_input;

  typedef typename Input_traits::Point_3                                   K1P;
  typedef typename Kernel_approx::Point_3                                  K2P;
  typedef typename Kernel_exact::Point_3                                   K3P;

public:
  static const bool value
    = (Has_conversion<Kernel_input, Kernel_approx, K1P, K2P>::value &&
       Has_conversion<Kernel_input, Kernel_exact, K1P, K3P>::value);
};

template < class Input_traits, class Kernel_approx, class Kernel_exact >
class Is_traits_point_convertible<Input_traits, Kernel_approx, Kernel_exact,
                                  ::CGAL::Tag_true /* Weighted_tag */>
{
  typedef typename Kernel_traits<typename Input_traits::Point_3>::Kernel   Kernel_input;

  typedef typename Input_traits::Weighted_point_3                          K1WP;
  typedef typename Kernel_approx::Weighted_point_3                         K2WP;
  typedef typename Kernel_exact::Weighted_point_3                          K3WP;

public:
  static const bool value
    = (Has_conversion<Kernel_input, Kernel_approx, K1WP, K2WP>::value &&
       Has_conversion<Kernel_input, Kernel_exact, K1WP, K3WP>::value);
};

template <class T>
struct Input_points_for_lazy_alpha_nt_3
{
  int nbpts;
  const T* p0;
  const T* p1;
  const T* p2;
  const T* p3;
};

//non-weighted case  
template <class Weighted_tag,class Input_traits,class Kernel_input,class Kernel_approx,class Kernel_exact>
struct Types_for_alpha_nt_3
{
//Converter types
  typedef CGAL::Cartesian_converter<Kernel_input,Kernel_approx>    To_approx;
  typedef CGAL::Cartesian_converter<Kernel_input,Kernel_exact>     To_exact;
//Point types
  typedef typename Kernel_approx::Point_3                          Approx_point;
  typedef typename Kernel_exact::Point_3                           Exact_point;
  typedef typename Input_traits::Point_3                           Input_point;
//Constructions 
  typedef typename Kernel_approx::Compute_squared_radius_3         Approx_squared_radius;
  typedef typename Kernel_exact::Compute_squared_radius_3          Exact_squared_radius;
};
  
  
//weighted case
template <class Input_traits,class Kernel_input,class Kernel_approx,class Kernel_exact>
struct Types_for_alpha_nt_3< ::CGAL::Tag_true,Input_traits,Kernel_input,Kernel_approx,Kernel_exact>
{
//Converter types
  typedef CGAL::Cartesian_converter<Kernel_input,Kernel_approx>   To_approx;
  typedef CGAL::Cartesian_converter<Kernel_input,Kernel_exact>    To_exact;
//Point types
  typedef typename Kernel_approx::Weighted_point_3 Approx_point;
  typedef typename Kernel_exact::Weighted_point_3  Exact_point;
  typedef typename Input_traits::Weighted_point_3  Input_point;
//Constructions 
  typedef typename Kernel_approx::Compute_squared_radius_smallest_orthogonal_sphere_3           Approx_squared_radius;
  typedef typename Kernel_exact::Compute_squared_radius_smallest_orthogonal_sphere_3            Exact_squared_radius;
};


template<class Input_traits, bool mode, class Weighted_tag>
class Lazy_alpha_nt_3{
//NT & kernels
  typedef CGAL::Interval_nt<mode>                                                               NT_approx;
  //Gmpq or Quotient<MP_float>
  typedef Exact_field_selector<double>::Type                                                    NT_exact;
  typedef CGAL::Simple_cartesian<NT_approx>                                                     Kernel_approx;
  typedef CGAL::Simple_cartesian<NT_exact>                                                      Kernel_exact;
  typedef typename Kernel_traits<typename Input_traits::Point_3>::Kernel   Kernel_input;

//Helper class for weighted and non-weighted case  
  typedef Types_for_alpha_nt_3<Weighted_tag,Input_traits,Kernel_input,Kernel_approx,Kernel_exact> Types;  
  
//Converters
  typedef typename Types::To_approx                                                             To_approx;
  typedef typename Types::To_exact                                                              To_exact;
 
//Constructions class
  typedef typename Types::Approx_squared_radius                                                 Approx_squared_radius;
  typedef typename Types::Exact_squared_radius                                                  Exact_squared_radius;
  
//Point
  typedef typename Types::Approx_point                                                          Approx_point;
  typedef typename Types::Exact_point                                                           Exact_point;
  typedef typename Types::Input_point                                                           Input_point;
//Convertion functions
  Approx_point to_approx(const Input_point& wp) const
  {
    // The traits class' Point_3 must be convertible using the Cartesian converter
    CGAL_static_assertion((Is_traits_point_convertible<
                            Input_traits, Kernel_approx, Kernel_exact, Weighted_tag>::value));

    To_approx converter;
    return converter(wp);
  }
  
  Exact_point to_exact(const Input_point& wp) const
  {
    // The traits class' Point_3 must be convertible using the Cartesian converter
    CGAL_static_assertion((Is_traits_point_convertible<
                            Input_traits, Kernel_approx, Kernel_exact, Weighted_tag>::value));

    To_exact converter;
    return converter(wp);
  }

//members  
  //the members can be updated when calling method exact()
  mutable boost::optional<NT_exact> exact_;
  mutable NT_approx approx_;

//private functions
  typedef Input_points_for_lazy_alpha_nt_3<Input_point> Data_vector;
  Data_vector input_points;

  const Data_vector& data() const{ return input_points;}
  Data_vector& data(){ return input_points;}

public:

  typedef NT_exact               Exact_nt;
  typedef NT_approx              Approximate_nt;

  void update_exact() const{
    switch (data().nbpts){
      case 1:
        exact_ = Exact_squared_radius()( to_exact(*data().p0) );
      break;
      case 2:
        exact_ = Exact_squared_radius()( to_exact(*data().p0),to_exact(*data().p1) );
      break;
      case 3:
        exact_ = Exact_squared_radius()( to_exact(*data().p0),to_exact(*data().p1),to_exact(*data().p2) );
      break;
      case 4:
        exact_ = Exact_squared_radius()( to_exact(*data().p0),to_exact(*data().p1),to_exact(*data().p2),to_exact(*data().p3) );
      break;
      default:
        CGAL_assertion(false);
    }
  }
  
  void set_approx(){
    switch (data().nbpts){
      case 1:
        approx_ = Approx_squared_radius()( to_approx(*data().p0) );
      break;
      case 2:
        approx_ = Approx_squared_radius()( to_approx(*data().p0),to_approx(*data().p1) );
      break;
      case 3:
        approx_ = Approx_squared_radius()( to_approx(*data().p0),to_approx(*data().p1),to_approx(*data().p2) );
      break;
      case 4:
        approx_ = Approx_squared_radius()( to_approx(*data().p0),to_approx(*data().p1),to_approx(*data().p2),to_approx(*data().p3) );
      break;
      default:
        CGAL_assertion(false);
    }    
  }

  const NT_exact& exact() const {
    if (exact_ == boost::none){
      update_exact();
      approx_=to_interval(*exact_);
    }
    return *exact_;
  }

  const NT_approx& approx() const{
    return approx_;
  }
//Constructors  
  Lazy_alpha_nt_3()
   : exact_(Exact_nt(0)),approx_(0)
  {
    data().nbpts=0;
    data().p0=NULL;
    data().p1=NULL;
    data().p2=NULL;
    data().p3=NULL;
  }
  
  Lazy_alpha_nt_3(double d)
   : exact_(Exact_nt(d)),approx_(d)
  {
    data().nbpts=0;
    data().p0=NULL;
    data().p1=NULL;
    data().p2=NULL;
    data().p3=NULL;
  }
  
  Lazy_alpha_nt_3(const Input_point& wp1)
  {
    data().nbpts=1;
    data().p0=&wp1;
    data().p1=NULL;
    data().p2=NULL;
    data().p3=NULL;
    set_approx();
  }

  Lazy_alpha_nt_3(const Input_point& wp1,
                  const Input_point& wp2)
  {
    data().nbpts=2;
    data().p0=&wp1;
    data().p1=&wp2;
    data().p2=NULL;
    data().p3=NULL;
    set_approx();
  }

  Lazy_alpha_nt_3(const Input_point& wp1,
                  const Input_point& wp2,
                  const Input_point& wp3)
  {
    data().nbpts=3;
    data().p0=&wp1;
    data().p1=&wp2;
    data().p2=&wp3;
    data().p3=NULL;
    set_approx();
  }

  Lazy_alpha_nt_3(const Input_point& wp1,
                  const Input_point& wp2,
                  const Input_point& wp3,
                  const Input_point& wp4)
  {
    data().nbpts=4;
    data().p0=&wp1;
    data().p1=&wp2;
    data().p2=&wp3;
    data().p3=&wp4;
    set_approx();
  }
    
  #define CGAL_LANT_COMPARE_FUNCTIONS(CMP) \
  bool \
  operator CMP (const Lazy_alpha_nt_3<Input_traits,mode,Weighted_tag> &other) const \
  { \
    Uncertain<bool> res = this->approx() CMP other.approx(); \
    if (res.is_certain()) \
      return res; \
    else \
      return this->exact() CMP other.exact(); \
  } \

  CGAL_LANT_COMPARE_FUNCTIONS(<)
  CGAL_LANT_COMPARE_FUNCTIONS(>)
  CGAL_LANT_COMPARE_FUNCTIONS(>=)
  CGAL_LANT_COMPARE_FUNCTIONS(<=)
  CGAL_LANT_COMPARE_FUNCTIONS(==)
  CGAL_LANT_COMPARE_FUNCTIONS(!=)

  #undef CGAL_LANT_COMPARE_FUNCTIONS  
};

template<class Input_traits, bool mode, class Weighted_tag>
std::ostream&
operator<< (std::ostream& os,const Lazy_alpha_nt_3<Input_traits,mode,Weighted_tag>& a){
  return os << ::CGAL::to_double(a.approx());
}
  
//small class to select predicate in weighted and unweighted case
template <class GeomTraits,class Weighted_tag>
struct iCompute_squared_radius_3;

template <class GeomTraits>
struct iCompute_squared_radius_3<GeomTraits,Tag_false>
{
  template <class As>
  typename GeomTraits::Compute_squared_radius_3
  operator()(const As& as) const{
    return static_cast<const typename As::Triangulation&>(as).geom_traits().compute_squared_radius_3_object();
  }
};

template <class GeomTraits>
struct iCompute_squared_radius_3<GeomTraits,Tag_true>
{
  template <class As>
  typename GeomTraits::Compute_squared_radius_smallest_orthogonal_sphere_3
  operator()(const As& as) const{
    return static_cast<const typename As::Triangulation&>(as).geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object();
  }
};

template <class Type_of_alpha,class Point>
struct Lazy_compute_squared_radius_3 {
  Type_of_alpha operator() (const Point& p, 
                 const Point& q , 
                 const Point& r, 
                 const Point& s)
  {return Type_of_alpha(p,q,r,s);}

  Type_of_alpha operator() ( const Point& p, 
                  const Point& q , 
                  const Point& r)
  {return Type_of_alpha(p,q,r); }

  Type_of_alpha operator() (const Point& p, 
                 const Point& q )
  {return Type_of_alpha(p,q); }

  Type_of_alpha operator() (const Point& p) 
  {return Type_of_alpha(p);}
};


template <class GeomTraits,class ExactAlphaComparisonTag,class Weighted_tag>
struct Alpha_nt_selector_impl_3;

template <class GeomTraits,class Weighted_tag>
struct Alpha_nt_selector_impl_3<GeomTraits,Tag_false,Weighted_tag>
{
  typedef typename GeomTraits::FT Type_of_alpha;
  typedef iCompute_squared_radius_3<GeomTraits,Weighted_tag> Compute_squared_radius_3;
};

template <class GeomTraits,class Weighted_tag>
struct Alpha_nt_selector_impl_3<GeomTraits,Tag_true,Weighted_tag>
{
  typedef Lazy_alpha_nt_3<GeomTraits,true,Tag_false> Type_of_alpha;
  typedef Lazy_compute_squared_radius_3<Type_of_alpha,typename GeomTraits::Point_3> Functor;
  struct Compute_squared_radius_3{
    template<class As>
    Functor operator()(const As&){return Functor();}    
  };
};

template <class GeomTraits>
struct Alpha_nt_selector_impl_3<GeomTraits,Tag_true,Tag_true>
{
  typedef Lazy_alpha_nt_3<GeomTraits,true,Tag_true> Type_of_alpha;
  typedef Lazy_compute_squared_radius_3<Type_of_alpha,typename GeomTraits::Weighted_point_3> Functor;
  struct Compute_squared_radius_3{
    template<class As>
    Functor operator()(const As&){return Functor();}    
  };
};

template <class GeomTraits,class ExactAlphaComparisonTag,class Weighted_tag>
struct Alpha_nt_selector_3
  : public Alpha_nt_selector_impl_3<
             GeomTraits,
             // If the base traits is already exact then we don't need to do anything,
             // and we can simply directly use the traits class
             Boolean_tag<boost::is_floating_point<typename GeomTraits::FT>::value &&
                         ExactAlphaComparisonTag::value >,
             Weighted_tag>
{ };

} //namespace internal

template<class Input_traits, bool mode, class Weighted_tag>
double to_double(const internal::Lazy_alpha_nt_3<Input_traits, mode, Weighted_tag>& a)
{
  return to_double(a.approx());
}

} //namespace CGAL

#endif //CGAL_INTERNAL_LAZY_ALPHA_NT_3_H

