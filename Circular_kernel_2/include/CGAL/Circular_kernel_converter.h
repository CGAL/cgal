// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CCIRCULAR_KERNEL_CONVERTER_H
#define CGAL_CCIRCULAR_KERNEL_CONVERTER_H 

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Algebraic_kernel_converter.h>
#include <CGAL/Object.h>

// TODO :
// - we should have a better default than Cartesian_converter.

namespace CGAL {

template < class C1, class C2,
  class LK_converter = Cartesian_converter<C1, C2>,
  //typename C1::Linear_kernel, typename C2::Linear_kernel>,
  class AK_converter = Algebraic_kernel_converter<typename C1::Algebraic_kernel,
                                                  typename C2::Algebraic_kernel > >
class Circular_kernel_converter
  : public LK_converter
{
public:

	typedef C1                       		      Source_kernel;
	typedef C2                          		      Target_kernel;
//	typedef typename C1::Linear_kernel  		      L1;
//	typedef typename C2::Linear_kernel  		      L2;
	typedef LK_converter                		      Linear_kernel_converter;
	typedef AK_converter                		      Algebraic_kernel_converter;
	typedef typename Linear_kernel_converter::Number_type_converter       RT_type_converter;
	typedef typename Algebraic_kernel_converter::Root_of_type_converter   Root_of_type_converter;




#ifdef _MSC_VER
   bool operator()(bool b) const { return Linear_kernel_converter::operator()(b); }

    Bounded_side operator()(Bounded_side bs) const {
      return Linear_kernel_converter::operator()(bs);
    }

    Comparison_result operator()(Comparison_result cr) const {
      return Linear_kernel_converter::operator()(cr);
    }

    Angle operator()(Angle a) const { return Base::operator()(a); }

    Origin
    operator()(const Origin& o) const
    {
        return o;
    }

    Null_vector
    operator()(const Null_vector& n) const
    {
        return n;
    }

    Bbox_2
    operator()(const Bbox_2& b) const
    {
        return b;
    }

    Bbox_3
    operator()(const Bbox_3& b) const
    {
        return b;
    }

    typename C2::FT
    operator()(const typename C1::FT &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Point_2
    operator()(const typename C1::Point_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Vector_2
    operator()(const typename C1::Vector_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Direction_2
    operator()(const typename C1::Direction_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Segment_2
    operator()(const typename C1::Segment_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Line_2
    operator()(const typename C1::Line_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Ray_2
    operator()(const typename C1::Ray_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }

    typename C2::Circle_2
    operator()(const typename C1::Circle_2 &a) const
    {
      return Linear_kernel_converter::operator()(a);
    }
#else
        using LK_converter::operator();
#endif


	typename C2::Circular_arc_point_2
	operator()(const typename C1::Circular_arc_point_2 &a) const
	{
		return typename C2::Circular_arc_point_2( typename C2::Circular_arc_point_2::Root_for_circles_2_2(
							     Root_of_type_converter()( a.x() ),
							     Root_of_type_converter()( a.y() )
	   	  )
	   );
	}


	typename C2::Circular_arc_2
	operator()(const typename C1::Circular_arc_2 &a) const
	{
		return typename C2::Circular_arc_2(operator()(a.supporting_circle()),
		                     	           operator()(a.source()),
		                     	           operator()(a.target()));	     
	}


	typename C2::Line_arc_2
	operator()(const typename C1::Line_arc_2 &a) const
	{
		return typename C2::Line_arc_2 (     operator()( a.supporting_line() ),
		                     		     operator()( a.source() ),
		                     		     operator()( a.target() ) );
	}
	

	  typename C2::Object_2
    operator()(const typename C1::Object_2 &obj) const
    {

      if (const typename C1::Circular_arc_2 * ptr = object_cast<typename C1::Circular_arc_2>(&obj)) {
        return make_object(operator()(*ptr));
      } else if (const typename C1::Circular_arc_point_2 * ptr = 
                 object_cast<typename C1::Circular_arc_point_2>(&obj)) {
        return make_object(operator()(*ptr));
      } else if (const std::pair<typename C1::Circular_arc_point_2,unsigned int> * ptr = 
                 object_cast<std::pair<typename C1::Circular_arc_point_2,unsigned int> >(&obj)) {
        return make_object(std::make_pair(operator()(ptr->first),ptr->second));
      } else if (const typename C1::Line_arc_2 * ptr = 
                 object_cast<typename C1::Line_arc_2>(&obj)) {
        return make_object(operator()(*ptr));
      }
      CGAL_error_msg("CircularK_converter is unable to determine what is wrapped in the Object");
      return Object();
	
    }
	

	
  std::vector<typename C2::Object_2>
  operator()(const std::vector<typename C1::Object_2>& v) const
  {
    std::vector<typename C2::Object_2> res;
    res.reserve(v.size());
    for(unsigned int i = 0; i < v.size(); i++){
      res.push_back(operator()(v[i]));
    }
    return res;
  }
	

        
  std::pair<typename C2::Line_arc_2,typename C2::Line_arc_2>
  operator()(const std::pair<typename C1::Line_arc_2,typename C1::Line_arc_2> &a) const
  {
    return std::make_pair (operator()( a.first ),
		           operator()( a.second ));
  }

  std::pair<typename C2::Circular_arc_2,typename C2::Circular_arc_2>
  operator()(const std::pair<typename C1::Circular_arc_2,typename C1::Circular_arc_2> &a) const
  {
    return std::make_pair (operator()( a.first ),
		           operator()( a.second ));
  }


}; 

} //namespace CGAL

#endif // CGAL_CCIRCULAR_KERNEL_CONVERTER_H
