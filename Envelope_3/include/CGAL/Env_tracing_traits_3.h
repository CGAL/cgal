// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
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
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*!
  \file   Env_tracing_traits_3.h
  \brief  The file is used to trace envelope_3 trace classes
  \todo   Revise this file.
  
*/


#ifndef CGAL_ENV_TRACE_TRAITS_H
#define CGAL_ENV_TRACE_TRAITS_H

#include <CGAL/license/Envelope_3.h>


/*! \file
 * This is the file containing a traits class used to trace other voronoi 
 * diagram traits.
 */

#include <list>
#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/functions_on_enums.h>
#include <CGAL/number_utils.h> 
#include <CGAL/Envelope_3/Envelope_base.h>

namespace CGAL {

/*!
 * \class The traits class
 */
template <class Traits_>
class Env_tracing_traits_3 : public Traits_
{
public:
  typedef Traits_                                         Base;
  
  typedef typename Base::Multiplicity                     Multiplicity;
  typedef typename Base::Point_2                          Point_2;
  typedef typename Base::Curve_2                          Curve_2;
  typedef typename Base::X_monotone_curve_2               X_monotone_curve_2;

  typedef typename Base::Boundary_category                Boundary_category;
    
  typedef typename Base::Xy_monotone_surface_3            Xy_monotone_surface_3;
  typedef typename Base::Surface_3                              Surface_3;
  
  typedef std::pair<X_monotone_curve_2, Multiplicity>     Intersection_curve;
    
  class Make_xy_monotone_3
  {
  public:
        
    template <class OutputIterator>
    OutputIterator operator()(const Surface_3& s,
                              bool is_lower,
                              OutputIterator o) const
    {
      Base base;
      return base.make_xy_monotone_3_object() (s, is_lower, o);
    }
  };
    
  Make_xy_monotone_3 make_xy_monotone_3_object() const
  {
    return Make_xy_monotone_3();
  }
    
  class Compare_z_at_xy_3
  {
  public:
        
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Base base;
      std::cerr << "Compare_z_at_xy_3:" << std::endl;
      std::cerr << "Point: " << p << std::endl;
      std::cerr << "Surface1: " << h1 << std::endl;
      std::cerr << "Surface2: " << h2 << std::endl;
      Comparison_result res = base.compare_z_at_xy_3_object() (p, h1, h2);
      std::cerr << "Result: " << res << std::endl;
      return res;
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Base base;
      std::cerr << "Compare_z_at_xy_3:" << std::endl;
      std::cerr << "Curve: " << cv << std::endl;
      std::cerr << "Surface1: " << h1 << std::endl;
      std::cerr << "Surface2: " << h2 << std::endl;
      Comparison_result res = base.compare_z_at_xy_3_object() (cv, h1, h2);
      std::cerr << "Result: " << res << std::endl;
      return res;
    }
        
        
    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const

    {
      Base base;
      std::cerr << "Compare_z_at_xy_3:" << std::endl;
      std::cerr << "Surface1: " << h1 << std::endl;
      std::cerr << "Surface2: " << h2 << std::endl;
      Comparison_result res = base.compare_z_at_xy_3_object() (h1, h2);
      std::cerr << "Result: " << res << std::endl;
      return res;
    }
    
  };
  
  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    return Compare_z_at_xy_3();
  }

  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Base base;
      std::cerr << "Compare_z_at_xy_above_3:" << std::endl;
      std::cerr << "Curve: " << cv << std::endl;
      std::cerr << "Surface1: " << h1 << std::endl;
      std::cerr << "Surface2: " << h2 << std::endl;
      Comparison_result res =
        base.compare_z_at_xy_above_3_object() (cv, h1, h2);
      std::cerr << "Result: " << res << std::endl;
      return res;
    }
  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }

  class Compare_z_at_xy_below_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Base base;
      std::cerr << "Compare_z_at_xy_below_3:" << std::endl;
      std::cerr << "Curve: " << cv << std::endl;
      std::cerr << "Surface1: " << h1 << std::endl;
      std::cerr << "Surface2: " << h2 << std::endl;
      Comparison_result res =
        base.compare_z_at_xy_below_3_object() (cv, h1, h2);
      std::cerr << "Result: " << res << std::endl;
      return res;
    }

  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const
  {
    return Compare_z_at_xy_below_3();
  }


  class Construct_projected_boundary_2
  {
  public:

    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const
    {
      Base base;
      std::cerr << "Construct_projected_boundary_2: JUST FIRST" << std::endl;
      std::cerr << "Surface: " << s << std::endl;
      std::list<CGAL::Object> l;
      base.construct_projected_boundary_2_object() (s, std::back_inserter(l));

      if (l.size() > 0)
      {
        std::pair<X_monotone_curve_2, CGAL::Oriented_side> i;
        if (CGAL::assign(i, l.front()))
          std::cerr << "First: " << i.first << std::endl;
        else
          std::cerr << "First intersection is a point" << std::endl;
      }
          
      std::copy(l.begin(), l.end(), o);
      return o;
    }
  };

  Construct_projected_boundary_2 
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2();
  }


  class Construct_projected_intersections_2
  {
  public:

    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const
    {
      Base base;
      std::cerr << "Construct_projected_intersections_2: JUST FIRST"
                << std::endl;
      std::cerr << "Surface1: " << s1 << std::endl;
      std::cerr << "Surface2: " << s2 << std::endl;
      std::list<CGAL::Object> l;
      base.construct_projected_intersections_2_object() (s1, s2,
                                                         std::back_inserter(l));
          
      if (l.size() > 0)
      {
        Intersection_curve i;
        if (CGAL::assign(i, l.front()))
          std::cerr << "First: " << i.first << std::endl;
        else
          std::cerr << "First intersection is not a point" << std::endl;
      }
          
      std::copy(l.begin(), l.end(), o);
      return o;
    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  {
    return Construct_projected_intersections_2();
  }

};

} //namespace CGAL

#endif
