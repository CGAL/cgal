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
// $URL$
// $Id$
// 
//
// Author(s)     : Sébastien Loriot <sebastien.loriot@geometryfactory.com>

#ifndef CGAL_INTERNAL_LAZY_ALPHA_NT_3_H
#define CGAL_INTERNAL_LAZY_ALPHA_NT_3_H

#include <CGAL/assertions.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <iostream>

namespace CGAL {

namespace internal{

//non-weighted case  
template <class Weighted_tag,class Input_traits,class Kernel_input,class Kernel_approx,class Kernel_exact>
struct Types_for_alpha_nt_3
{
//Converter types
  typedef CGAL::Cartesian_converter<Kernel_input,Kernel_approx>    To_approx;
  typedef CGAL::Cartesian_converter<Kernel_input,Kernel_exact>     To_exact;
//Traits types
  
//Point types
  typedef typename Kernel_approx::Point_3                          Approx_point;
  typedef typename Kernel_exact::Point_3                           Exact_point;
  typedef typename Kernel_input::Point_3                           Input_point;
//Constructions 
  typedef typename Kernel_approx::Compute_squared_radius_3         Approx_squared_radius;
  typedef typename Kernel_exact::Compute_squared_radius_3          Exact_squared_radius;
};
  
  
//weighted case
template <class Input_traits,class Kernel_input,class Kernel_approx,class Kernel_exact>
struct Types_for_alpha_nt_3< ::CGAL::Tag_true,Input_traits,Kernel_input,Kernel_approx,Kernel_exact>
{
//Converter types
  typedef CGAL::Weighted_converter_3< CGAL::Cartesian_converter<Kernel_input,Kernel_approx> >   To_approx;
  typedef CGAL::Weighted_converter_3< CGAL::Cartesian_converter<Kernel_input,Kernel_exact> >    To_exact;
//Traits types
  typedef ::CGAL::Regular_triangulation_euclidean_traits_3<Kernel_approx>                       Approx_traits;
  typedef ::CGAL::Regular_triangulation_euclidean_traits_3<Kernel_exact>                        Exact_traits;
//Point types
  typedef typename Approx_traits::Weighted_point Approx_point;
  typedef typename Exact_traits::Weighted_point  Exact_point;
  typedef typename Input_traits::Weighted_point  Input_point;
//Constructions 
  typedef typename Approx_traits::Compute_squared_radius_smallest_orthogonal_sphere_3           Approx_squared_radius;
  typedef typename Exact_traits::Compute_squared_radius_smallest_orthogonal_sphere_3            Exact_squared_radius; 
};


template<class Input_traits, class Kernel_input, bool mode, class Weighted_tag>
class Lazy_alpha_nt_3{
//NT & kernels
  typedef CGAL::Interval_nt<mode>                                                               NT_approx;
  //Gmpq or Quotient<MP_float>
  typedef Exact_field_selector<double>::Type                                                    NT_exact;
  typedef CGAL::Simple_cartesian<NT_approx>                                                     Kernel_approx;
  typedef CGAL::Simple_cartesian<NT_exact>                                                      Kernel_exact;
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
  Approx_point to_approx(const Input_point& wp) const {
    To_approx converter;
    return converter(wp);
  }
  
  Exact_point to_exact(const Input_point& wp) const {
    To_exact converter;
    return converter(wp);
  }


//members  
  unsigned nb_pt;
  //the members can be updated when calling method exact()
  mutable bool updated;
  mutable NT_exact exact_;
  mutable NT_approx approx_;
  typedef std::vector<const Input_point*> Data_vector;
  boost::shared_ptr<Data_vector> inputs_ptr;

//private functions  
  const Data_vector& data() const{ return *inputs_ptr;}

  Data_vector& 
  data(){ return *inputs_ptr;}  
  
  
public:

  typedef NT_exact               Exact_nt;
  typedef NT_approx              Approximate_nt;

  void update_exact() const{
    switch (nb_pt){
      case 1:
        exact_ = Exact_squared_radius()( to_exact(*data()[0]) );
      break;
      case 2:
        exact_ = Exact_squared_radius()( to_exact(*data()[0]),to_exact(*data()[1]) );
      break;
      case 3:
        exact_ = Exact_squared_radius()( to_exact(*data()[0]),to_exact(*data()[1]),to_exact(*data()[2]) );
      break;
      case 4:
        exact_ = Exact_squared_radius()( to_exact(*data()[0]),to_exact(*data()[1]),to_exact(*data()[2]),to_exact(*data()[3]) );
      break;
      default:
        CGAL_assertion(false);
    }
    updated=true;
  }
  
  void set_approx(){
    switch (nb_pt){
      case 1:
        approx_ = Approx_squared_radius()( to_approx(*data()[0]) );
      break;
      case 2:
        approx_ = Approx_squared_radius()( to_approx(*data()[0]),to_approx(*data()[1]) );
      break;
      case 3:
        approx_ = Approx_squared_radius()( to_approx(*data()[0]),to_approx(*data()[1]),to_approx(*data()[2]) );
      break;
      case 4:
        approx_ = Approx_squared_radius()( to_approx(*data()[0]),to_approx(*data()[1]),to_approx(*data()[2]),to_approx(*data()[3]) );
      break;
      default:
        CGAL_assertion(false);
    }    
  }

  const NT_exact& exact() const {
    if (!updated){
      update_exact();
      approx_=to_interval(exact_);
    }
    return exact_;
  }

  const NT_approx& approx() const{
    return approx_;
  }
//Constructors  
  Lazy_alpha_nt_3():nb_pt(0),updated(true),exact_(0),approx_(0){}
  
  Lazy_alpha_nt_3(double d):nb_pt(0),updated(true),exact_(d),approx_(d){}
  
  Lazy_alpha_nt_3(const Input_point& wp1):nb_pt(1),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    set_approx();
  }

  Lazy_alpha_nt_3(const Input_point& wp1,
           const Input_point& wp2):nb_pt(2),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    data().push_back(&wp2);
    set_approx();
  }

  Lazy_alpha_nt_3(const Input_point& wp1,
           const Input_point& wp2,
           const Input_point& wp3):nb_pt(3),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    data().push_back(&wp2);
    data().push_back(&wp3);
    set_approx();
  }

  Lazy_alpha_nt_3(const Input_point& wp1,
           const Input_point& wp2,
           const Input_point& wp3,
           const Input_point& wp4):nb_pt(4),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    data().push_back(&wp2);
    data().push_back(&wp3);
    data().push_back(&wp4);
    set_approx();
  }
    
  #define CGAL_LANT_COMPARE_FUNCTIONS(CMP) \
  bool \
  operator CMP (const Lazy_alpha_nt_3<Input_traits,Kernel_input,mode,Weighted_tag> &other) const \
  { \
    try{ \
      return this->approx() CMP other.approx(); \
    } \
    catch(CGAL::Uncertain_conversion_exception&){ \
      return this->exact() CMP other.exact(); \
    } \
  } \

  CGAL_LANT_COMPARE_FUNCTIONS(<)
  CGAL_LANT_COMPARE_FUNCTIONS(>)
  CGAL_LANT_COMPARE_FUNCTIONS(>=)
  CGAL_LANT_COMPARE_FUNCTIONS(<=)
  CGAL_LANT_COMPARE_FUNCTIONS(==)
  CGAL_LANT_COMPARE_FUNCTIONS(!=)

  #undef CGAL_LANT_COMPARE_FUNCTIONS  
};

template<class Input_traits, class Kernel_input, bool mode, class Weighted_tag>
std::ostream&
operator<< (std::ostream& os,const Lazy_alpha_nt_3<Input_traits,Kernel_input,mode,Weighted_tag>& a){
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
  typedef Lazy_alpha_nt_3<GeomTraits,GeomTraits,true,Tag_false> Type_of_alpha;
  typedef Lazy_compute_squared_radius_3<Type_of_alpha,typename GeomTraits::Point_3> Functor;
  struct Compute_squared_radius_3{
    template<class As>
    Functor operator()(const As&){return Functor();}    
  };
};

template <class GeomTraits>
struct Alpha_nt_selector_impl_3<GeomTraits,Tag_true,Tag_true>
{
  typedef Lazy_alpha_nt_3<GeomTraits,typename GeomTraits::Kernel,true,Tag_true> Type_of_alpha;
  typedef Lazy_compute_squared_radius_3<Type_of_alpha,typename GeomTraits::Weighted_point> Functor;
  struct Compute_squared_radius_3{
    template<class As>
    Functor operator()(const As&){return Functor();}    
  };
};

template <class GeomTraits,class ExactAlphaComparisonTag,class Weighted_tag>
struct Alpha_nt_selector_3:
  public Alpha_nt_selector_impl_3<GeomTraits,
              Boolean_tag< boost::is_floating_point<typename GeomTraits::FT>::value && ExactAlphaComparisonTag::value >,
              Weighted_tag>
{};


} //namespace internal

} //namespace CGAL

#endif //CGAL_INTERNAL_LAZY_ALPHA_NT_3_H

