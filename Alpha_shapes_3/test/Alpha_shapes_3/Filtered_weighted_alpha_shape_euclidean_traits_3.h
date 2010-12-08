// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh:// $
// $Id: $
// 
//
// Author(s)     : Sébastien Loriot <Sebastien.Loriot@sophia.inria.fr>

#ifndef CGAL_FILTERED_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_FILTERED_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <boost/shared_ptr.hpp>

namespace CGAL {

namespace internal{

template<template <class T,class W,bool UseFilteredPredicates> class Traits, class Kernel_input, bool mode>
class Alpha_nt{
//NT & kernels
  typedef CGAL::Interval_nt<mode> NT_approx;
  typedef CGAL::Gmpq NT_exact;
  typedef CGAL::Simple_cartesian<NT_approx> Kernel_approx;
  typedef CGAL::Simple_cartesian<NT_exact>  Kernel_exact;
  
//Converters
  typedef CGAL::Weighted_converter_3< CGAL::Cartesian_converter<Kernel_input,Kernel_approx> >   To_approx;
  typedef CGAL::Weighted_converter_3< CGAL::Cartesian_converter<Kernel_input,Kernel_exact> >    To_exact;
  
//Traits
  typedef Traits<Kernel_input,typename Kernel_input::FT,Kernel_input::Has_filtered_predicates>  Input_traits;
  typedef Traits<Kernel_approx,NT_approx,Kernel_approx::Has_filtered_predicates>                Approx_traits;
  typedef Traits<Kernel_exact,NT_exact,Kernel_exact::Has_filtered_predicates>                   Exact_traits;
  
//Constructions class
  typedef typename Approx_traits::Compute_squared_radius_smallest_orthogonal_sphere_3 Approx_squared_radius;
  typedef typename Exact_traits::Compute_squared_radius_smallest_orthogonal_sphere_3  Exact_squared_radius;
  
//Convertion functions
  static 
  typename Approx_traits::Weighted_point
  to_approx(const typename Input_traits::Weighted_point& wp) {
    static To_approx converter;
    return converter(wp);
  }
  
  static  
  typename Exact_traits::Weighted_point
  to_exact(const typename Input_traits::Weighted_point& wp) {
    static To_exact converter;
    return converter(wp);
  }  


//members  
  unsigned nb_pt;
  //the members can be updated when calling method exact()
  mutable bool updated;
  mutable NT_exact exact_;
  mutable NT_approx approx_;
  typedef std::vector<const typename Input_traits::Weighted_point*> Data_vector;
  boost::shared_ptr<Data_vector> inputs_ptr;

//private functions  
  const Data_vector& data() const{ return *inputs_ptr;}

  Data_vector& 
  data(){ return *inputs_ptr;}  
  
  
public:

  void update_exact() const{
    switch (nb_pt){
      case 1:
        exact_ = -NT_exact( data()[0]->weight() );
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
        assert(false);
    }
    updated=true;
  }
  
  void set_approx(){
    switch (nb_pt){
      case 1:
        approx_ = - NT_approx( data()[0]->weight() );
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
        assert(false);
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
  Alpha_nt():nb_pt(0),updated(true),exact_(0),approx_(0){}
  
  Alpha_nt(const double& d):nb_pt(0),updated(true),exact_(d),approx_(d){}
  
  Alpha_nt(const typename Input_traits::Weighted_point_3& wp1):nb_pt(1),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    set_approx();
  }

  Alpha_nt(const typename Input_traits::Weighted_point_3& wp1,
           const typename Input_traits::Weighted_point_3& wp2):nb_pt(2),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    data().push_back(&wp2);
    set_approx();
  }

  Alpha_nt(const typename Input_traits::Weighted_point_3& wp1,
           const typename Input_traits::Weighted_point_3& wp2,
           const typename Input_traits::Weighted_point_3& wp3):nb_pt(3),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    data().push_back(&wp2);
    data().push_back(&wp3);
    set_approx();
  }

  Alpha_nt(const typename Input_traits::Weighted_point_3& wp1,
           const typename Input_traits::Weighted_point_3& wp2,
           const typename Input_traits::Weighted_point_3& wp3,
           const typename Input_traits::Weighted_point_3& wp4):nb_pt(4),updated(false),inputs_ptr(new Data_vector())
  {
    data().reserve(nb_pt);
    data().push_back(&wp1);
    data().push_back(&wp2);
    data().push_back(&wp3);
    data().push_back(&wp4);
    set_approx();
  }
  
};
  

unsigned& nb_call_total(){
  static unsigned n=0;
  return n;
}

unsigned& nb_failure(){
  static unsigned n=0;
  return n;
}

#define COMPARE_FUNCTIONS(CMP) \
template<template <class T,class W,bool b> class Traits, class Kernel_input,bool mode> \
inline \
bool \
operator CMP (const Alpha_nt<Traits,Kernel_input,mode> &a, const Alpha_nt<Traits,Kernel_input,mode> &b) \
{ \
 ++nb_call_total();\
  try{ \
    return a.approx() CMP b.approx(); \
  } \
  catch(CGAL::Uncertain_conversion_exception& e){ \
    ++nb_failure();\
    return a.exact() CMP b.exact(); \
  } \
} \
\
template<template <class T,class W,bool b> class Traits, class Kernel_input,bool mode> \
inline \
bool \
operator CMP (const Alpha_nt<Traits,Kernel_input,mode> &a, const double &b) \
{ \
  try{ \
    return a.approx() CMP b; \
  } \
  catch(CGAL::Uncertain_conversion_exception& e){ \
    return a.exact() CMP b; \
  } \
} \
\
template<template <class T,class W,bool b> class Traits, class Kernel_input,bool mode> \
inline \
bool \
operator CMP (const double& a, const Alpha_nt<Traits,Kernel_input,mode> &b) \
{ \
  try{ \
    return a CMP b.approx(); \
  } \
  catch(CGAL::Uncertain_conversion_exception& e){ \
    return a CMP b.exact(); \
  } \
} \

COMPARE_FUNCTIONS(<)
COMPARE_FUNCTIONS(>)
COMPARE_FUNCTIONS(>=)
COMPARE_FUNCTIONS(<=)
COMPARE_FUNCTIONS(==)
COMPARE_FUNCTIONS(!=)

} //namespace internal
   
//------------------ Traits class -------------------------------------

template <class K,bool mode>
class Filtered_weighted_alpha_shape_euclidean_traits_3: public 
Regular_triangulation_euclidean_traits_3<K>
{
  typedef internal::Alpha_nt<CGAL::Regular_triangulation_euclidean_traits_3,K,mode> Alpha_nt;
public:
  typedef Regular_triangulation_euclidean_traits_3<K> Base;
  typedef typename Base::Side_of_bounded_orthogonal_sphere_3 
                                       Side_of_bounded_sphere_3;
  
  typedef Alpha_nt                          FT;


 class  Compute_squared_radius_3 {
    typedef typename Base::Weighted_point   Weighted_point_3;
  public:
    FT operator() (const Weighted_point_3& p, 
                   const Weighted_point_3& q , 
                   const Weighted_point_3& r, 
                   const Weighted_point_3& s)
    {return FT(p,q,r,s);}

    FT operator() ( const Weighted_point_3& p, 
                    const Weighted_point_3& q , 
                    const Weighted_point_3& r)
    {return FT(p,q,r); }

    FT operator() (const Weighted_point_3& p, 
                   const Weighted_point_3& q )
    {return FT(p,q); }

    FT operator() (const Weighted_point_3& p) 
    {return FT(p);}
  };
 
  

  //---------------------------------------------------------------------

  Compute_squared_radius_3 
  compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_3();
    }
  //---------------------------------------------------------------------

  Side_of_bounded_sphere_3 
  side_of_bounded_sphere_3_object() const
    {
      return Side_of_bounded_sphere_3();
    }
};

} //namespace CGAL

#endif //CGAL_FILTERED_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

