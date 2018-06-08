// Copyright (c) 2018  GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_VERTEX_TO_POINT_TRAITS_ADAPTER_3_H
#define CGAL_VERTEX_TO_POINT_TRAITS_ADAPTER_3_H
#include <boost/call_traits.hpp>
#include <CGAL/property_map.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/convex_hull_3.h>


namespace CGAL {

template<class Base_traits,class VertexPointMap>
class Vertex_to_point_traits_adapter
    :public Base_traits
{
  VertexPointMap vpm_;
  
public:
  Vertex_to_point_traits_adapter(const VertexPointMap& vpmap, Base_traits base=Base_traits())
    :Base_traits(base), vpm_(vpmap)
  {}
  //typedef typename boost::property_traits<VertexPointMap>::key_type Vertex;
  //typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename boost::property_traits<VertexPointMap>::key_type Vertex;
  typedef Vertex Point_3;
  class Compute_x_3:public Base_traits::Compute_x_3
  {
    VertexPointMap vpm_;
    const typename Base_traits::Compute_x_3& base;
    typedef typename Base_traits::FT             FT;
  public:
    Compute_x_3(const VertexPointMap& map,const typename Base_traits::Compute_x_3& base):
      Base_traits::Compute_x_3(base),vpm_(map),base(base){}
    typedef const FT&                  result_type;
    
    result_type
    operator()(const Point_3& v) const
    {
      return base(get(vpm_, v));
    }
    
  };
  Compute_x_3 compute_x_3_object () const {return Compute_x_3(vpm_,static_cast<const Base_traits*>(this)->compute_x_3_object() );}
  
  class Compute_y_3:public Base_traits::Compute_y_3
  {
    VertexPointMap vpm_;
    const typename Base_traits::Compute_y_3& base;
    typedef typename Base_traits::FT             FT;
  public:
    Compute_y_3(const VertexPointMap& map,const typename Base_traits::Compute_y_3& base):
      Base_traits::Compute_y_3(base),vpm_(map), base(base){}
    typedef const FT&                  result_type;
    
    result_type
    operator()(const Point_3& v) const
    {
      return base(get(vpm_, v));
    }
    
  };
  Compute_y_3 compute_y_3_object () const {return Compute_y_3(vpm_,static_cast<const Base_traits*>(this)->compute_y_3_object() );}
  
  class Compute_z_3:public Base_traits::Compute_z_3
  {
    VertexPointMap vpm_;
    const typename Base_traits::Compute_z_3& base;
    typedef typename Base_traits::FT             FT;
  public:
    Compute_z_3(const VertexPointMap& map,const typename Base_traits::Compute_z_3& base):
      Base_traits::Compute_z_3(base),vpm_(map), base(base){}
    typedef const FT&                  result_type;
    
    result_type
    operator()(const Point_3& v) const
    {
      return base(get(vpm_, v));
    }
  };
  Compute_z_3 compute_z_3_object () const {return Compute_z_3(vpm_,static_cast<const Base_traits*>(this)->compute_z_3_object() );}
  
  class Equal_3:public Base_traits::Equal_3
  {
    VertexPointMap vpm_;
    const typename Base_traits::Equal_3& base;
  public:
    Equal_3(const VertexPointMap& map,const typename Base_traits::Equal_3& base):
      Base_traits::Equal_3(base),vpm_(map), base(base){}
    typedef bool       result_type;
    
    result_type
    operator()(const Point_3 &p, const Point_3 &q) const
    {
      return base(get(vpm_, p), get(vpm_, q));
    }
  };
  Equal_3 equal_3_object () const {return Equal_3(vpm_,static_cast<const Base_traits*>(this)->equal_3_object() );}
  
  class Collinear_3:public Base_traits::Collinear_3
  {
    VertexPointMap vpm_;
    const typename Base_traits::Collinear_3& base;
  public:
    Collinear_3(const VertexPointMap& map,const typename Base_traits::Collinear_3& base):
      Base_traits::Collinear_3(base),vpm_(map), base(base){}
    typedef bool    result_type;
    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      
      return base(get(vpm_,p), get(vpm_,q), get(vpm_,r));
    }
  };
  Collinear_3 collinear_3_object () const {return Collinear_3(vpm_,static_cast<const Base_traits*>(this)->collinear_3_object() );}
  
  class Coplanar_3:public Base_traits::Coplanar_3
  {
    VertexPointMap vpm_;
    typename Base_traits::Orientation_3 o;
    const typename Base_traits::Coplanar_3& base;
  public:
    Coplanar_3(const VertexPointMap& map,const typename Base_traits::Coplanar_3& base):
      Base_traits::Coplanar_3(base),vpm_(map),base(base){}
    typedef bool    result_type;
    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& s) const
    {
      return base(get(vpm_,p), get(vpm_,q), get(vpm_,r), get(vpm_,s));
    }
  };
  Coplanar_3 coplanar_3_object () const {return Coplanar_3(vpm_,static_cast<const Base_traits*>(this)->coplanar_3_object() );}
  
  class Less_distance_to_point_3:public Base_traits::Less_distance_to_point_3
  {
    VertexPointMap vpm_;
    const typename Base_traits::Less_distance_to_point_3& base;
  public:
    Less_distance_to_point_3(const VertexPointMap& map,const typename Base_traits::Less_distance_to_point_3& base):
      Base_traits::Less_distance_to_point_3(base),vpm_(map),base(base){}
    typedef bool       result_type;
    
    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return base(get(vpm_,p), get(vpm_,q), get(vpm_,r));
    }
  };
  Less_distance_to_point_3 less_distance_to_point_3_object() const 
  {return Less_distance_to_point_3(vpm_,static_cast<const Base_traits*>(this)->less_distance_to_point_3_object() );}
  
  class Less_signed_distance_to_plane_3
      :public Base_traits::Less_signed_distance_to_plane_3
  {  
    VertexPointMap vpm_;
    const typename Base_traits::Less_signed_distance_to_plane_3& base;
  public:
    typedef typename Base_traits::Plane_3 Plane_3;
    
    typedef bool             result_type;
    
    Less_signed_distance_to_plane_3(
        const VertexPointMap& map,
        const typename Base_traits::Less_signed_distance_to_plane_3& base):
      Base_traits::Less_signed_distance_to_plane_3(base),vpm_(map), base(base){}
    
    bool
    operator()( const Plane_3& h, const Vertex& p, const Vertex& q) const
    {
      return base(h, get(vpm_,p), get(vpm_,q));
    }
  };
  Less_signed_distance_to_plane_3 less_signed_distance_to_plane_3_object() const 
  {return Less_signed_distance_to_plane_3(
          vpm_,static_cast<const Base_traits*>(this)->less_signed_distance_to_plane_3_object() );}
  
  class Construct_plane_3:public Base_traits::Construct_plane_3
  {  
    VertexPointMap vpm_;
    const typename Base_traits::Construct_plane_3& base;
  public:
    Construct_plane_3(const VertexPointMap& map, const typename Base_traits::Construct_plane_3& base):
      Base_traits::Construct_plane_3(base),vpm_(map), base(base){}
    typename Base_traits::Plane_3 operator()(const Point_3& p, const Point_3& q, const Point_3& r)const 
    {
      return base(get(vpm_,p),get(vpm_,q),get(vpm_,r));
    }
  };
  Construct_plane_3 construct_plane_3_object() const
  {return Construct_plane_3(vpm_,static_cast<const Base_traits*>(this)->construct_plane_3_object());}
  
  class Has_on_positive_side_3:public Base_traits::Has_on_positive_side_3
  {  
    VertexPointMap vpm_;
    const typename Base_traits::Has_on_positive_side_3& base;
  public:
    Has_on_positive_side_3(const VertexPointMap& map,const typename Base_traits::Has_on_positive_side_3& base):
      Base_traits::Has_on_positive_side_3(base),vpm_(map), base(base){}
    
    typedef typename Base_traits::Plane_3          Plane_3;
  public:
    typedef bool          result_type;
    
    result_type
    operator()( const Plane_3& pl, const Point_3& p) const
    { 
      return base(pl, get(vpm_, p));
    }
  };
  Has_on_positive_side_3 has_on_positive_side_3_object() const {return Has_on_positive_side_3(
          vpm_,static_cast<const Base_traits*>(this)->has_on_positive_side_3_object() );}
  
  template<class Base_proj_traits>
  class Proj_traits_3:public Base_proj_traits
  {
    VertexPointMap vpm_;
    typedef Base_proj_traits Btt;
  public:
    Proj_traits_3(const VertexPointMap& map,const Btt& base):
      Base_proj_traits(base),vpm_(map){}
    typedef Point_3 Point_2;
    
    class Equal_2:public Btt::Equal_2
    {
      VertexPointMap vpm_;
      const typename Btt::Equal_2& base;
    public:
      Equal_2(const VertexPointMap& map,const typename Btt::Equal_2& base):
        Btt::Equal_2(base),vpm_(map), base(base){}
    public:      
      bool operator()(Point_2 p, Point_2 q) const
      { 
        return base(get(vpm_, p), get(vpm_, q));
      }
    };
    Equal_2 equal_2_object () const {return Equal_2(vpm_,static_cast<const Btt*>(this)->equal_2_object() );}
    
    class Less_xy_2:public Btt::Less_xy_2
    {
      VertexPointMap vpm_;
      const typename Btt::Less_xy_2& base;
    public:
      Less_xy_2(const VertexPointMap& map,const typename Btt::Less_xy_2& base):
        Btt::Less_xy_2(base),vpm_(map), base(base){}
    public:      
      bool operator()(Point_2 p, Point_2 q) const
      { 
        return base(get(vpm_, p), get(vpm_, q));
      }
    };
    Less_xy_2 less_xy_2_object ()const{return Less_xy_2(vpm_,static_cast<const Btt*>(this)->less_xy_2_object() );}
    
    class Less_yx_2:public Btt::Less_yx_2
    {
      VertexPointMap vpm_;
      const typename Btt::Less_yx_2& base;
    public:
      Less_yx_2(const VertexPointMap& map,const typename Btt::Less_yx_2& base):
        Btt::Less_yx_2(base),vpm_(map), base(base){}
    public:      
      bool operator()(Point_2 p, Point_2 q) const
      { 
        return base(get(vpm_, p), get(vpm_, q));
      }
    };
    Less_yx_2 less_yx_2_object ()const{return Less_yx_2(vpm_,static_cast<const Btt*>(this)->less_yx_2_object() );}
    
    class Less_signed_distance_to_line_2:public Btt::Less_signed_distance_to_line_2
    {
      VertexPointMap vpm_;
      const typename Btt::Less_signed_distance_to_line_2& base;
    public:
      Less_signed_distance_to_line_2(const VertexPointMap& map,const typename Btt::Less_signed_distance_to_line_2& base):
        Btt::Less_signed_distance_to_line_2(base),vpm_(map), base(base){}
    public:      
      bool operator()(Point_2 p, Point_2 q, Point_2 r,Point_2 s) const
      { 
        return base(get(vpm_, p), get(vpm_, q), get(vpm_, r), get(vpm_, s));
      }
    };
    Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object ()const
    {return Less_signed_distance_to_line_2(vpm_,static_cast<const Btt*>(this)->Less_signed_distance_to_line_2() );}
    
    class Less_rotate_ccw_2:public Btt::Less_rotate_ccw_2
    {
      VertexPointMap vpm_;
      const typename Btt::Less_rotate_ccw_2& base;
    public:
      Less_rotate_ccw_2(const VertexPointMap& map,const typename Btt::Less_rotate_ccw_2& base):
        Btt::Less_rotate_ccw_2(base),vpm_(map), base(base){}
    public:      
      bool operator()(Point_2 e, Point_2 p,Point_2 q) const
      { 
        return base(get(vpm_, e), get(vpm_, p), get(vpm_, q));
      }
    };
    Less_rotate_ccw_2 less_rotate_ccw_2_object ()const
    {return Less_rotate_ccw_2(vpm_,static_cast<const Btt*>(this)->less_rotate_ccw_2_object() );}
    
    class Left_turn_2:public Btt::Left_turn_2
    {
      VertexPointMap vpm_;                                                           
      const typename Btt::Left_turn_2& base;                                   
    public:                                                                          
      Left_turn_2(const VertexPointMap& map,const typename Btt::Left_turn_2& base):
        Btt::Left_turn_2(base),vpm_(map), base(base){}                         
    public:                                                                          
      bool operator()(Point_2 p, Point_2 q, Point_2 r) const                                    
      {                                                                              
        return base(get(vpm_, p), get(vpm_, q), get(vpm_, r));
      }
    };
    Left_turn_2 left_turn_2_object ()const{return Left_turn_2(vpm_,static_cast<const Btt*>(this)->left_turn_2_object() );}
    
    class Orientation_2:public Btt::Orientation_2
    {
      VertexPointMap vpm_;                                                           
      const typename Btt::Orientation_2& base;                                   
    public:                                                                          
      Orientation_2(const VertexPointMap& map,const typename Btt::Orientation_2& base):
        Btt::Orientation_2(base),vpm_(map), base(base){}                         
      
      typename CGAL::Orientation operator()(Point_2 e,Point_2 p, Point_2 q) const                                    
      {                                                                              
        return base(get(vpm_, e), get(vpm_, p), get(vpm_, q));                                     
      }
    };
    Orientation_2 orientation_2_object ()const{return Orientation_2(vpm_,static_cast<const Btt*>(this)->orientation_2_object() );}
  };
  
  typedef Proj_traits_3<typename Base_traits::Traits_xy_3> Traits_xy_3;
  typedef Proj_traits_3<typename Base_traits::Traits_yz_3> Traits_yz_3;
  typedef Proj_traits_3<typename Base_traits::Traits_xz_3> Traits_xz_3;
  
  Traits_xy_3 construct_traits_xy_3_object()const
  {return Traits_xy_3(vpm_, static_cast<const Base_traits*>(this)->construct_traits_xy_3_object());}
  Traits_yz_3 construct_traits_yz_3_object()const
  {return Traits_yz_3(vpm_, static_cast<const Base_traits*>(this)->construct_traits_yz_3_object());}
  Traits_xz_3 construct_traits_xz_3_object()const
  {return Traits_xz_3(vpm_, static_cast<const Base_traits*>(this)->construct_traits_xz_3_object());}
  
};

}//end CGAL

#endif // CGAL_VERTEX_TO_POINT_TRAITS_ADAPTER_3_H
