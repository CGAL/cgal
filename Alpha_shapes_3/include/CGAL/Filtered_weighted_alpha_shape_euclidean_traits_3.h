// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sébastien Loriot <sebastien.loriot@geometryfactory.com>

#ifndef CGAL_FILTERED_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_FILTERED_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

#include <CGAL/Filtered_alpha_shape_euclidean_traits_3.h>

namespace CGAL{

template <class K,bool mode=true>
class Filtered_weighted_alpha_shape_euclidean_traits_3: public 
Regular_triangulation_euclidean_traits_3<K>
{
  typedef internal::Alpha_nt<Regular_triangulation_euclidean_traits_3<K>,K,mode,Tag_true> Alpha_nt;
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

}//namespace CGAL

#endif //CGAL_FILTERED_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 
