// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Planar_map_2/Pm_segment_slim_traits.h
// package       : Planar_map (5.80)
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// author(s)     : Efi Fogel <efif@post.tau.ac.il>
//                 Shai Hirsch <shaihi@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin halperin<@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_PM_SEGMENT_UTILITIES_2_H
#define CGAL_PM_SEGMENT_UTILITIES_2_H

CGAL_BEGIN_NAMESPACE

/*! Construct a direction in the same direction as a given segment which source
 * is the given point
 * \pre the given point must be one of the segment end-points
 */
template <class Kernel_, class Segment_2_>
struct Construct_direction_at_endpoint_2 : 
public Kernel_::Construct_direction_2
{
  typedef Kernel_                                Kernel;
  typedef Segment_2_                             Segment_2;
    
  typedef typename Kernel::Construct_direction_2 Base;
  typedef typename Kernel::Direction_2           Direction_2;
  typedef typename Kernel::Point_2               Point_2;
  typedef typename Kernel::Equal_2               Equal_2;

  Construct_direction_at_endpoint_2() {}
  Construct_direction_at_endpoint_2(const Kernel & in_kernel) 
    : m_kernel(in_kernel) {}
    
  Direction_2 operator()(const Segment_2 & cv, const Point_2 & point) const
  { 
    Equal_2 is_equal = m_kernel.equal_2_object();
    if (is_equal(m_kernel.construct_vertex_2_object()(cv, 0), point)) {
      const Direction_2 & d = Base::operator()(cv);
      return m_kernel.construct_opposite_direction_2_object()(d);
    }
    return Base::operator()(cv);
  }

private:
  Kernel m_kernel;
};

/*!
 */
template<class Kernel_>
struct Is_vertex_for_segments_2 
{
public:
  typedef Kernel_                     Kernel;
  typedef typename Kernel::Point_2    Point_2;
  typedef typename Kernel::Segment_2  Segment_2;

  Is_vertex_for_segments_2() {}
  Is_vertex_for_segments_2(const Kernel & in_kernel) : m_kernel(in_kernel) {}

  bool
  operator()(const Point_2 & point, Segment_2 segment)
  {
    typedef typename Kernel::Has_on_2           Has_on_2;
    typedef typename Kernel::Equal_2            Equal_2;
    typedef typename Kernel::Construct_vertex_2 Construct_vertex_2;

    Has_on_2 has_on = m_kernel.has_on_2_object();

    // if point is on segment
    if ( has_on(segment, point) )
    {
      Equal_2            equal = m_kernel.equal_2_object();
      Construct_vertex_2 construct_vertex = 
                                 m_kernel.construct_vertex_2_object();

      const Point_2 & source = construct_vertex(segment, 0);
      // is point is segment's source
      if ( equal(source, point) ) 
	return true;
      else
      {
	const Point_2 & target = construct_vertex(segment, 1);
	// is point is segment's target
	if ( equal(target, point) ) 
	  return true; 
	else 
	  return false;
      }
    }

    return false;
  }

private:
  Kernel m_kernel;
};

/*! 
 */
template <class Kernel_, class Segment_2_>
struct Counterclockwise_in_between_for_segments_2 :
public Kernel_::Counterclockwise_in_between_2
{
  typedef Kernel_                                        Kernel;
  typedef Segment_2_                                     Segment_2;

  typedef typename Kernel::Counterclockwise_in_between_2 Base;
  typedef typename Kernel::Point_2                       Point_2;
  typedef typename Kernel::Direction_2                   Direction_2;
  typedef CGAL::Construct_direction_at_endpoint_2<Kernel, Segment_2>
                                                         Construct_direction_2;
  typedef struct Is_vertex_for_segments_2<Kernel>        Is_vertex_2;

  Counterclockwise_in_between_for_segments_2() {}
  Counterclockwise_in_between_for_segments_2(const Kernel & in_kernel) 
    : m_kernel(in_kernel) {}

  bool
  operator()(const Point_2 & point,
             const Segment_2 & cv, 
	     const Segment_2 & first, 
	     const Segment_2 & second) const
  {
    // Preconditions:
    CGAL_precondition_code(Is_vertex_2 is_vertex;);
    CGAL_precondition_msg( is_vertex(point, cv),
			   "point should be an endpoint of cv.");
    CGAL_precondition_msg( is_vertex(point, first),
			   "point should be an endpoint of first.");
    CGAL_precondition_msg( is_vertex(point, second),
			   "point should be an endpoint of second.");

    Construct_direction_2 construct_direction_2 = 
      construct_direction_2_object();

    const Direction_2 & d  = construct_direction_2(cv,     point);
    const Direction_2 & d1 = construct_direction_2(first,  point);
    const Direction_2 & d2 = construct_direction_2(second, point);

    return m_kernel.counterclockwise_in_between_2_object()(d, d1, d2);
  }

private:
  inline Construct_direction_2 construct_direction_2_object() const
  { return Construct_direction_2(); }
private:
  Kernel m_kernel;
};

// SHAI: Turn this into a functor
template <class Segment_2>
Segment_2 curve_flip(const Segment_2 &cv)
{
  typedef typename Segment_2::R Kernel;
  return Kernel().construct_opposite_segment_2_object()(cv);
}

CGAL_END_NAMESPACE

#endif  // CGAL_PM_SEGMENT_UTILITIES_2_H
