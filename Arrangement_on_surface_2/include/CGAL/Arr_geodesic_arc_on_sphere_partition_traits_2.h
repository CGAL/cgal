// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ophir Setter         <ophir.setter@cs.au.ac.il>


#ifndef CGAL_ARR_GEODESIC_ARC_ON_SPHERE_PARTITION_TRAITS_2_H
#define CGAL_ARR_GEODESIC_ARC_ON_SPHERE_PARTITION_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * The partition traits class for geodesic arcs on the sphere enables 
 * the partition of geodesic polygons to convex polygons. It models the 
 * concept YMonotonePartitionTraits_2.
 * This partition of geodesic polygons is garenteed to work only for polygons
 * that are contained in a hemisphere and that do not intersect one of the 
 * boundaries. 
 * For larger polygons there is a chance that at least one steiner point may 
 * have to be added; see manuscript by Prof. Dan Halperin from 2008.
 * PAY ATTENTION TO THE FACT THAT WE REVERSE THE ROLES OF X AND Y SO IT WILL
 * BE EASIER TO IMPLEMENT A MODEL FOR YMonotonePartitionTraits_2 (actually
 * implementing XMonotonePartitionTraits_2)
 */


#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#include <vector>

namespace CGAL {


template <class T_Kernel, class Container_P =
          std::vector<typename Arr_geodesic_arc_on_sphere_traits_2< T_Kernel >::Point_2> >
class Arr_geodesic_arc_on_sphere_partition_traits_2 
  : public Arr_geodesic_arc_on_sphere_traits_2< T_Kernel >
{
private:
  typedef Arr_geodesic_arc_on_sphere_partition_traits_2< T_Kernel > Self;
  typedef Arr_geodesic_arc_on_sphere_traits_2< T_Kernel >           Base;
 
public:
 
  /// \name PartitionTraits_2 concept
  //@{
 
  typedef typename Base::Point_2                                    Point_2;
 
  /*! Class to represent polygon. Contain Vertex_const_iterator and access 
   *  functions.
   */
  class Polygon_2 : public Container_P
  {
  public:
    typedef typename Container_P::const_iterator         Vertex_const_iterator;
   
    Vertex_const_iterator vertices_begin() const
    { return this->begin();} 
   
    Vertex_const_iterator vertices_end() const
    { return this->end();}
  };

  /*! A functor that compares two points lexigoraphically: by x, then by y.
   */
  class Less_xy_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
   
  public:
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Less_xy_2(const Self * traits) : m_traits(traits) {}
   
    /*! Compare two points lexigoraphically: by x, then by y.
     *  We actually reversing the order, so x <--> y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return true - y(p1) < y(p2);
     *         true - y(p1) = y(p2) and x(p1) < x(p2);
     *         false - otherwise.
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
     
      // compare y and then x (reverse the order of x and y).
      Comparison_result res = m_traits->compare_y(p1, p2);
      if (res == EQUAL) 
        return m_traits->compare_x(p1, p2) == SMALLER;
      return res == SMALLER;
      //     return m_traits->compare_xy(p1, p2) == SMALLER;
    }
  };

  Less_xy_2 less_xy_2_object() const { return Less_xy_2(this); }

  /*! A functor that compares two points lexigoraphically: by y, then by x.
   */
  class Less_yx_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
      
  public:

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Less_yx_2(const Self * traits) : m_traits(traits) {}


    /*! Compare two points lexigoraphically: by y, then by x.
     *  We actually reversing the order, so x <--> y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return true - x(p1) < x(p2);
     *         true - x(p1) = x(p2) and y(p1) < y(p2);
     *         false - otherwise.
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
     
      // compare x and then y (reverse the order of x and y).
      return m_traits->compare_xy(p1, p2) == SMALLER;
    }
  };

  Less_yx_2 less_yx_2_object() const { return Less_yx_2(this); }
 

  /*! A functor that checks orientation of three points.
   */
  class Orientation_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
   
  public:

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Orientation_2 (const Self * traits) : m_traits(traits) {}

    /*! Checks the orientation between three points.
     *  We actually reversing the order, so x <--> y.
     * \param p 
     * \param q
     * \param r 
     * \return CGAL::LEFT_TURN - r lies to the "left" of the oriented line l 
     *         defined by p and q (do not forget that x <--> y.
     *         CGAL::RIGHT_TURN - r lies to the "right" of l
     *         CGAL::COLLINEAR - r lies on l. 
     * \pre p does not lie on the boundary.
     * \pre q does not lie on the boundary.
     * \pre r does not lie on the boundary.
     * \pre p and q are not antipodal (from the req. that the whole polygon
     *      is contained in a hemisphere.
     */
    CGAL::Orientation operator()(const Point_2 &p,
                                 const Point_2 &q,
                                 const Point_2 &r) const
    {
      CGAL_precondition(p.is_no_boundary());
      CGAL_precondition(q.is_no_boundary());
      CGAL_precondition(r.is_no_boundary());
     
      // the orientation is determined by the relative position of r with 
      // respect to the plane that contains p and q.
      typename Base::Vector_3 normal = 
        m_traits->construct_cross_product_vector_3_object() (p.vector(),
                                                             q.vector());

      Oriented_side res = CGAL::sign(normal * r.vector());

      return (res == ON_NEGATIVE_SIDE) ? RIGHT_TURN : 
        ((res == ON_POSITIVE_SIDE) ? LEFT_TURN : COLLINEAR);
    }
  };

  Orientation_2 orientation_2_object() const { return Orientation_2(this); }

  /*! A functor that checks if three points create a left turn.
   *  See Orientation_2 above.
   */
  class Left_turn_2
  {
  protected:
    const Self * m_traits;
   
  public:
    Left_turn_2 (const Self * traits) : m_traits(traits) {}

    bool operator()(const Point_2 &p,
                    const Point_2 &q,
                    const Point_2 &r) const
    {
      return m_traits->orientation_2_object()(p, q, r) == LEFT_TURN;
    }
  };

  Left_turn_2 left_turn_2_object() const { return Left_turn_2(this); }

  /*! As x switches parts with y (x <--> y) the compare_x_2 from base can
   *  be used here as compare_y_2.
   * */
  typedef typename Base::Compare_x_2              Compare_y_2;
 
  Compare_y_2 compare_y_2_object() const {return Base::compare_x_2_object(); }


  /*! A functor that compares two points by x coordinate.
   */
  class Compare_x_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
   
  public:
 
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_2(const Self * traits) : m_traits(traits) {}
   

    /*! Compare two points by y coordinate.
     *  We actually reversing the order, so x <--> y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2);
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
     
      // compare y and then x (reverse the order of x and y).
      return m_traits->compare_y(p1, p2);
    }
  };

  Compare_x_2 compare_x_2_object() const { return Compare_x_2(this); }
 
  //@}
  
  /// \name ConvexPartitionIsValidTraits_2 concept
  /// For now, we have a stub implementation. I am not sure how easy it is
  /// to create the real implementation.
  //@{

  struct Is_convex_2 
  {
    template <typename T> Is_convex_2(T) {}
    
    template<class InputIterator>
    bool operator ()(InputIterator first, InputIterator beyond) const
    { return true; }
  };

  Is_convex_2 is_convex_2_object() const { return Is_convex_2(); }

  struct Is_valid 
  {
    template <typename T> Is_valid(T) {}

    template<class InputIterator>
    bool operator ()(InputIterator first, InputIterator beyond) const
    { return true; }
  };

  // Is_valid is_valid_object() const { return Is_valid(); }
  
  //@}
 
};


} //namespace CGAL

#endif  // CGAL_ARR_GEODESIC_ARC_ON_SPHERE_PARTITION_TRAITS_2_H


// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ophir Setter         <ophir.setter@cs.au.ac.il>


#ifndef CGAL_ARR_GEODESIC_ARC_ON_SPHERE_PARTITION_TRAITS_2_H
#define CGAL_ARR_GEODESIC_ARC_ON_SPHERE_PARTITION_TRAITS_2_H

/*! \file
 * The partition traits class for geodesic arcs on the sphere enables 
 * the partition of geodesic polygons to convex polygons. It models the 
 * concept YMonotonePartitionTraits_2.
 * This partition of geodesic polygons is garenteed to work only for polygons
 * that are contained in a hemisphere and that do not intersect one of the 
 * boundaries. 
 * For larger polygons there is a chance that at least one steiner point may 
 * have to be added; see manuscript by Prof. Dan Halperin from 2008.
 * PAY ATTENTION TO THE FACT THAT WE REVERSE THE ROLES OF X AND Y SO IT WILL
 * BE EASIER TO IMPLEMENT A MODEL FOR YMonotonePartitionTraits_2 (actually
 * implementing XMonotonePartitionTraits_2)
 */


#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#include <vector>

namespace CGAL {


template <class T_Kernel, class Container_P = std::vector<
                            typename Arr_geodesic_arc_on_sphere_traits_2< T_Kernel >::Point_2> >
class Arr_geodesic_arc_on_sphere_partition_traits_2 
  : public Arr_geodesic_arc_on_sphere_traits_2< T_Kernel >
{
private:
  typedef Arr_geodesic_arc_on_sphere_partition_traits_2< T_Kernel > Self;
  typedef Arr_geodesic_arc_on_sphere_traits_2< T_Kernel >           Base;
 
public:
 
  /// \name PartitionTraits_2 concept
  //@{
 
  typedef typename Base::Point_2                                    Point_2;
 
  /*! Class to represent polygon. Contain Vertex_const_iterator and access 
   *  functions.
   */
  class Polygon_2 : public Container_P
  {
  public:
    typedef typename Container_P::const_iterator                     Vertex_const_iterator;
   
    Vertex_const_iterator vertices_begin() const
    { return this->begin();} 
   
    Vertex_const_iterator vertices_end() const
    { return this->end();}
  };
 
  /*! A functor that compares two points lexigoraphically: by x, then by y.
   */
  class Less_xy_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
   
  public:

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Less_xy_2(const Self * traits) : m_traits(traits) {}
   
    /*! Compare two points lexigoraphically: by x, then by y.
     *  We actually reversing the order, so x <--> y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return true - y(p1) < y(p2);
     *         true - y(p1) = y(p2) and x(p1) < x(p2);
     *         false - otherwise.
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
     
      // compare y and then x (reverse the order of x and y).
      Comparison_result res = m_traits->compare_y(p1, p2);
      if (res == EQUAL) 
        return compare_x(p1, p2) == SMALLER;
      return res == SMALLER;
      //     return m_traits->compare_xy(p1, p2) == SMALLER;
    }
  };

  Less_xy_2 less_xy_2_object() const { return Less_xy_2(this); }

  /*! A functor that compares two points lexigoraphically: by y, then by x.
   */
  class Less_yx_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
     
  public:

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Less_yx_2(const Self * traits) : m_traits(traits) {}

    /*! Compare two points lexigoraphically: by y, then by x.
     *  We actually reversing the order, so x <--> y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return true - x(p1) < x(p2);
     *         true - x(p1) = x(p2) and y(p1) < y(p2);
     *         false - otherwise.
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
     
      // compare x and then y (reverse the order of x and y).
      return m_traits->compare_xy(p1, p2) == SMALLER;
    }
  };

  Less_yx_2 less_yx_2_object() const { return Less_yx_2(this); }
 

  /*! A functor that checks orientation of three points.
   */
  class Orientation_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
   
  public:
   
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Orientation_2 (const Self * traits) : m_traits(traits) {}

    /*! Checks the orientation between three points.
     *  We actually reversing the order, so x <--> y.
     * \param p 
     * \param q
     * \param r 
     * \return CGAL::LEFT_TURN - r lies to the "left" of the oriented line l 
     *         defined by p and q (do not forget that x <--> y.
     *         CGAL::RIGHT_TURN - r lies to the "right" of l
     *         CGAL::COLLINEAR - r lies on l. 
     * \pre p does not lie on the boundary.
     * \pre q does not lie on the boundary.
     * \pre r does not lie on the boundary.
     * \pre p and q are not antipodal (from the req. that the whole polygon
     *      is contained in a hemisphere.
     */
    CGAL::Orientation operator()(const Point_2 &p,
                                 const Point_2 &q,
                                 const Point_2 &r) const
    {
      CGAL_precondition(p.is_no_boundary());
      CGAL_precondition(q.is_no_boundary());
      CGAL_precondition(r.is_no_boundary());
     
      // the orientation is determined by the relative position of r with 
      // respect to the plane that contains p and q.
      typename Base::Vector_3 normal = 
        m_traits->construct_cross_product_vector_3_object() (p.vector(),
                                                             q.vector());

      Oriented_side res = m_traits->oriented_side(normal, r);

      return (res == ON_NEGATIVE_SIDE) ? RIGHT_TURN : 
        ((res == ON_POSITIVE_SIDE) ? LEFT_TURN : COLLINEAR);
    }
  };

  Orientation_2 orientation_2_object() const { return Orientation_2(this); }

  /*! A functor that checks if three points create a left turn.
   *  See Orientation_2 above.
   */
  class Left_turn_2
  {
  protected:
    const Self * m_traits;
   
  public:

    Left_turn_2 (const Self * traits) : m_traits(traits) {}

    bool operator()(const Point_2 &p,
                    const Point_2 &q,
                    const Point_2 &r) const
    {
      return m_traits->orientation_2_object()(p, q, r) == LEFT_TURN;
    }
  };

  Left_turn_2 left_turn_2_object() const { return Left_turn_2(this); }

  /*! As x switches parts with y (x <--> y) the compare_x_2 from base can
   *  be used here as compare_y_2.
   */
  typedef typename Base::Compare_x_2              Compare_y_2;
 
  Compare_y_2 compare_y_2_object() const {return Base::compare_x_2_object(); }

  /*! A functor that compares two points by x coordinate.
   */
  class Compare_x_2
  {
  protected:
    /*! The traits (in case it has state) */
    const Self * m_traits;
      
  public:

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_2(const Self * traits) : m_traits(traits) {}

    /*! Compare two points by y coordinate.
     *  We actually reversing the order, so x <--> y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2);
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
     
      // compare y and then x (reverse the order of x and y).
      return m_traits->compare_y(p1, p2);
    }
  };

  Compare_x_2 compare_x_2_object() const { return Compare_x_2(this); }
 
  //@}
  
  /// \name ConvexPartitionIsValidTraits_2 concept
  /// For now, we have a stub implementation. I am not sure how easy it is
  /// to create the real implementation.
  //@{

  class Is_convex_2 
  {
    template<class InputIterator>
    bool operator ()(InputIterator first, InputIterator beyond) const
    { return true; }
  };

  Is_convex_2 is_convex_2_object() const { return Is_convex_2(); }
 
  //@}
 
};


} //namespace CGAL

#endif  // CGAL_ARR_GEODESIC_ARC_ON_SPHERE_PARTITION_TRAITS_2_H
