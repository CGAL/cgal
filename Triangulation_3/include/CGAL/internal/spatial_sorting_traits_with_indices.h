// Copyright (c) 2010   GeometryFactory (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_INTERNAL_SPATIAL_SORTING_TRAITS_WITH_INDICES
#define CGAL_INTERNAL_SPATIAL_SORTING_TRAITS_WITH_INDICES

#include <boost/mpl/has_xxx.hpp>
#include <CGAL/Triangulation_vertex_base_with_info_3.h> 

namespace CGAL {

namespace internal{

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_typedef_Info,Info,false)  
 
template <class T,bool has_info=Has_typedef_Info<T>::value>
 struct Info_check{
 struct type{};
};

template <class T>
struct Info_check<T,true>{
 typedef typename T::Info type;
};

template <class T,bool use_reference = (sizeof(T) > sizeof(int*))>
struct Arg_type_selection;

template <class T>
struct Arg_type_selection<T,false>{
  typedef T type;
};

template <class T>
struct Arg_type_selection<T,true>{
  typedef const T& type;
};

template <class T>
class Vector_property_map{
  const std::vector<T>& data;
public:
  typedef std::size_t key_type;
  Vector_property_map(const std::vector<T>& input):data(input){}
  
  const T& operator[](key_type i) const{
    return data[i];
  }
};

template<class Base_traits,class PointPropertyMap>
class Spatial_sort_traits_with_property_map_3:public Base_traits{
  PointPropertyMap accessor_;
public:
  Spatial_sort_traits_with_property_map_3(const PointPropertyMap& accessor,Base_traits base=Base_traits())
  :Base_traits(base),accessor_(accessor){}

  typedef Base_traits Gt;
  typedef typename PointPropertyMap::key_type Point_3;
  typedef typename Arg_type_selection<Point_3>::type Arg_type;

  struct Less_x_3 : public Base_traits::Less_x_3{
    Less_x_3(const PointPropertyMap& accessor,const typename Base_traits::Less_x_3& base):
      Base_traits::Less_x_3(base),accessor_(accessor){}
    const PointPropertyMap& accessor_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_x_3*>(this)->operator()(accessor_[p],accessor_[q]);
    }
  };

  struct Less_y_3 : public Base_traits::Less_y_3{
    Less_y_3(const PointPropertyMap& accessor,const typename Base_traits::Less_y_3& base):
      Base_traits::Less_y_3(base),accessor_(accessor){}
    const PointPropertyMap& accessor_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_y_3*>(this)->operator()(accessor_[p],accessor_[q]);
    }
  };

  struct Less_z_3 : public Base_traits::Less_z_3{
    Less_z_3(const PointPropertyMap& accessor,const typename Base_traits::Less_z_3& base):
      Base_traits::Less_z_3(base),accessor_(accessor){}
    const PointPropertyMap& accessor_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_z_3*>(this)->operator()(accessor_[p],accessor_[q]);
    }
  };

  Less_x_3 less_x_3_object () const {return Less_x_3(accessor_,static_cast<const Gt*>(this)->less_x_3_object() );}
  Less_y_3 less_y_3_object () const {return Less_y_3(accessor_,static_cast<const Gt*>(this)->less_y_3_object() );}
  Less_z_3 less_z_3_object () const {return Less_z_3(accessor_,static_cast<const Gt*>(this)->less_z_3_object() );}
};

template<class Base_traits,class PointPropertyMap>
class Spatial_sort_traits_with_property_map_2:public Base_traits{
  PointPropertyMap accessor_;
public:
  Spatial_sort_traits_with_property_map_2(const PointPropertyMap& accessor,Base_traits base=Base_traits())
  :Base_traits(base),accessor_(accessor){}

  typedef Base_traits Gt;
  typedef typename PointPropertyMap::key_type Point_2;
  typedef typename Arg_type_selection<Point_2>::type Arg_type;

  struct Less_x_2 : public Base_traits::Less_x_2{
    Less_x_2(const PointPropertyMap& accessor,const typename Base_traits::Less_x_2& base):
      Base_traits::Less_x_2(base),accessor_(accessor){}
    const PointPropertyMap& accessor_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_x_2*>(this)->operator()(accessor_[p],accessor_[q]);
    }
  };

  struct Less_y_2 : public Base_traits::Less_y_2{
    Less_y_2(const PointPropertyMap& accessor,const typename Base_traits::Less_y_2& base):
      Base_traits::Less_y_2(base),accessor_(accessor){}
    const PointPropertyMap& accessor_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_y_2*>(this)->operator()(accessor_[p],accessor_[q]);
    }
  };

  Less_x_2 less_x_2_object () const {return Less_x_2(accessor_,static_cast<const Gt*>(this)->less_x_2_object() );}
  Less_y_2 less_y_2_object () const {return Less_y_2(accessor_,static_cast<const Gt*>(this)->less_y_2_object() );}

};
} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_SPATIAL_SORTING_TRAITS_WITH_INDICES
