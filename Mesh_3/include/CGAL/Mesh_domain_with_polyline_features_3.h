// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H
#define CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H

#include <CGAL/iterator.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <vector>
#include <set>
#include <map>
#include <boost/next_prior.hpp> // for boost::prior and boost::next
#include <boost/variant.hpp>

namespace CGAL {


namespace internal {
namespace Mesh_3 {

template <typename Kernel>
class Polyline
{
  typedef typename Kernel::Point_3  Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::FT       FT;

  typedef std::vector<Point_3>      Data;
  
public:
  typedef typename Data::const_iterator const_iterator;
  
  Polyline() {}
  ~Polyline() {}
  
  /// Add a point at the end of the polyline
  void add_point(const Point_3& p)
  {
    CGAL_assertion( points_.empty() || p != end_point() );
    points_.push_back(p);
  }
  
  /// Returns the starting point of the polyline
  const Point_3& start_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.front();
  }
  
  /// Returns the ending point of the polyline
  const Point_3& end_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.back();          
  }
  
  /// Returns true if the polyline is not degenerated
  bool is_valid() const
  {
    return points_.size() > 1;
  }
  
  /// Returns true if polyline is a cycle
  bool is_cycle() const
  {
    return start_point() == end_point();
  }
  
  /// Returns the length of the polyline
  FT length() const
  {
    //TODO: cache result
    FT result (0);
    const_iterator it = points_.begin();
    const_iterator previous = it++;
    
    for ( const_iterator end = points_.end() ; it != end ; ++it, ++previous )
    {
      result += distance(*previous,*it);
    }
    
    return result;
  }
  
  /// Returns signed geodesic distance between \c p and \c q
  FT geodesic_distance(const Point_3& p, const Point_3& q,
                       bool /*treat_cycle*/=true) const
  {
    CGAL_precondition(is_valid());
    
    // Locate p & q on polyline
    const_iterator pit = locate(p);
    const_iterator qit = locate(q,false);
    
    // Compute geodesic distance
    FT result = (pit <= qit) ? geodesic_distance(p,q,pit,qit)
                             : -geodesic_distance(q,p,qit,pit);
 
    // Treat cycles: return a positive value
    if ( is_cycle() && (p==q || result < FT(0)) )
    {
      result = length() + result;
    }

    return result;
  }
  
  
  /// Returns a point at geodesic distance \c distance from p along the
  /// polyline. The polyline is oriented from starting point to end point.
  /// The distance could be negative.
  Point_3 point_at(const Point_3& p, FT distance) const
  {
    // use first point of the polyline instead of p
    distance += geodesic_distance(start_point(),p,false);
    
    // If polyline is a cycle, ensure that distance is given from start_point()
    if ( is_cycle() )
    {
      if ( distance < FT(0) ) { distance += length(); }
      else if ( distance > length() ) { distance -= length(); }
    }
    
    CGAL_assertion( distance > FT(0) );
    CGAL_assertion( distance < length() );
    
    // Initialize iterators
    const_iterator pit = points_.begin();
    const_iterator previous = pit++;
    
    // Iterate to find which segment contains the point we want to construct
    FT segment_length = this->distance(*previous,*pit);
    while ( distance > segment_length )
    {
      distance -= segment_length;
      
      // Increment iterators and update length
      ++previous;
      ++pit;
      CGAL_assertion(pit != points_.end());
      
      segment_length = this->distance(*previous,*pit);
    }
    
    // return point at distance from current segment source
    typedef typename Kernel::Vector_3 Vector_3;
    Vector_3 v (*previous, *pit);
    
    return (*previous) + (distance / CGAL::sqrt(v.squared_length())) * v;
  }
  
  bool are_ordered_along(const Point_3& p, const Point_3& q) const
  {
    CGAL_precondition(!is_cycle());
    
    // Locate p & q on polyline
    const_iterator pit = locate(p);
    const_iterator qit = locate(q,true);
    
    // Points are not located on the same segment
    if ( pit != qit ) { return (pit <= qit); }
    
    // pit == qit, then we have to sort p&q along (pit,pit+1)
    return ( compare_distance(*pit,p,q) != CGAL::LARGER );
  }
  
private:
  const_iterator first_segment_source() const
  {
    CGAL_precondition(is_valid());
    return points_.begin();
  }
  
  const_iterator last_segment_source() const
  {
    CGAL_precondition(is_valid());
    return (points_.end() - 2);
  }
  
  FT geodesic_distance(const Point_3& p, const Point_3& q,
                       const_iterator pit, const_iterator qit) const
  {
    CGAL_precondition(std::distance(pit,qit) >= 0);
    
    // If p and q are in the same segment of the polyline
    if ( pit == qit )
    {
      FT result = distance(p,q);
      
      // Find the closest point to *pit
      if ( compare_distance(*pit,p,q) != CGAL::LARGER )
      { return result; }
      else
      { return -result; }
    }
    
    // p is inside [pit,pit+1], pit+1 != qit, q is inside [qit,qit+1]
    FT result = distance(p,*(pit+1));
    result += distance(*qit,q);
    
    // Add segments between pit+1 and qit to result
    for ( const_iterator it = (pit+1) ; it != qit ; ++it )
    {
      result += distance(*it,*(it+1));
    }
    
    return result;
  }

  /// Returns an iterator on the starting point of the segment of the 
  /// polyline which contains p
  /// if end_point_first is true, then --end is returned instead of begin
  /// if p is the starting point of a cycle.
  const_iterator locate(const Point_3& p, bool end_point_first=false) const
  {
    CGAL_precondition(is_valid());
    
    // First look if p is one of the points of the polyline
    const_iterator result = std::find(points_.begin(), points_.end(), p);
    if ( result != points_.end() )
    {
      if ( result != points_.begin() )
      { return --result; }
      else
      { 
        // Treat cycles
        if ( end_point_first && p == end_point() )
        { return last_segment_source(); }
        else
        { return result; }
      }
    }

    CGAL_assertion(result == points_.end());

    // Get result by projecting p on the polyline
    const_iterator it = points_.begin();
    const_iterator previous = it;
    Segment_3 nearest_segment;
    const_iterator nearest_vertex = it;
    result = nearest_vertex;
    bool nearest_is_a_segment = false;
    
    while ( ++it != points_.end() )
    {
      Segment_3 seg (*previous, *it);

      if(nearest_is_a_segment)
      {
        if(compare_distance(p, seg, nearest_segment) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          result = previous;
        }
        if(compare_distance(p, *it, nearest_segment) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          nearest_is_a_segment = false;
          result = it;
        }
      }
      else {
        if(compare_distance(p, *it, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          result = it;
        }
        if(compare_distance(p, seg, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          nearest_is_a_segment = true;
          result = previous;
        }
      }
      previous = it;
    } // end the while loop on the vertices of the polyline


    if(result == points_.begin()) {
      return (end_point_first && !nearest_is_a_segment) ? last_segment_source() : points_.begin();
    } else {
      return result;
    }
  }
  
  // FT squared_distance(const Point_3& p, const Point_3& q) const
  // {
  //   typename Kernel::Compute_squared_distance_3 sq_distance =
  //     Kernel().compute_squared_distance_3_object();
  //   return sq_distance(p,q);
  // }
  
  FT distance(const Point_3& p, const Point_3& q) const
  {
    return CGAL::sqrt(squared_distance(p, q));
  }

  Angle angle(const Point_3& p, 
              const Point_3& angle_vertex_point,
              const Point_3& q) const 
  {
    typename Kernel::Angle_3 compute_angle =  Kernel().angle_3_object();
    return compute_angle(p,angle_vertex_point,q);
  }
  
  template <typename T1, typename T2>
  CGAL::Sign compare_distance(const Point_3& p,
                              const T1& obj1,
                              const T2& obj2) const
  {
    typename Kernel::Compare_distance_3 compare_distance =
      Kernel().compare_distance_3_object();
    return compare_distance(p,obj1,obj2);
  }

private:
  Data points_;
}; // end class Polyline
  

}
}





/**
 * @class Mesh_domain_with_polyline_features_3
 *
 *
 */
template < typename MeshDomain >
class Mesh_domain_with_polyline_features_3
  : public MeshDomain
{
  typedef MeshDomain Base;

public:
  // Index types
  typedef typename Base::Index    Index;
  typedef int                     Curve_segment_index;
  typedef int                     Corner_index;

  typedef typename Base::R         Gt;
  typedef Gt                       R;
  typedef typename Base::Point_3   Point_3;
  typedef typename Gt::FT          FT;
  
  typedef CGAL::Tag_true           Has_features;

  /// Constructors
  /// Call the base class constructor
  Mesh_domain_with_polyline_features_3()
    : Base()
    , current_curve_index_(1) {}
  
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <typename ... T>
  Mesh_domain_with_polyline_features_3(const T& ...t)
    : Base(t...)
    , current_curve_index_(1) {}
#else
  template <typename T1>
  Mesh_domain_with_polyline_features_3(const T1& o1)
    : Base(o1),
      current_curve_index_(1) {}

  template <typename T1, typename T2>
  Mesh_domain_with_polyline_features_3(const T1& o1, const T2& o2)
    : Base(o1, o2),
      current_curve_index_(1) {}

  template <typename T1, typename T2, typename T3>
  Mesh_domain_with_polyline_features_3(const T1& o1, const T2& o2, 
                                       const T3& o3)
    : Base(o1, o2, o3),
      current_curve_index_(1) {}
#endif

  /// Destructor
  ~Mesh_domain_with_polyline_features_3() {}

  /// OutputIterator value type is std::pair<Corner_index, Point_3>
  template <typename OutputIterator>
  OutputIterator get_corners(OutputIterator out) const;
  
  /// OutputIterator value type is CGAL::cpp11::tuple<Curve_segment_index,
  /// std::pair<Point_3,Index>, std::pair<Point_3,Index> >
  template <typename OutputIterator>
  OutputIterator get_curve_segments(OutputIterator out) const;

  /// Returns the geodesic distance between points p and q of curve 
  /// \c curve_index
  FT geodesic_distance(const Point_3& p, const Point_3& q,
                       const Curve_segment_index& curve_index) const;

  /// Construct a point on curve \c curve_index at geodesic distance \c distance
  /// of \c starting_point
  Point_3
  construct_point_on_curve_segment(const Point_3& starting_point,
                                   const Curve_segment_index& curve_index,
                                   FT distance) const;
  
  /// Returns the sign of the orientation of p,q,r along curve segment
  /// of index \c index
  CGAL::Sign distance_sign_along_cycle(const Point_3& p,
                                       const Point_3& q,
                                       const Point_3& r,
                                       const Curve_segment_index& index) const;
  
  /// Returns true if curve \c curve_index is a cycle
  bool is_cycle(const Point_3&, const Curve_segment_index& index) const;
  
  /// Returns an Index from a Curve_segment_index
  Index index_from_curve_segment_index(const Curve_segment_index& index) const
  { return Index(index); }

  /// Returns an Curve_segment_index from an Index
  Curve_segment_index curve_segment_index(const Index& index) const
  { return boost::get<Curve_segment_index>(index); }
  
  /// Returns an Index from a Corner_index
  Index index_from_corner_index(const Corner_index& index) const
  { return Index(index); }
  
  /// Returns an Corner_index from an Index
  Corner_index corner_index(const Index& index) const
  { return boost::get<Corner_index>(index); }

  /// Insert a bunch of edges into domain
  ///   + InputIterator type should have begin() and end() function
  ///   + InputIterator::iterator value type must be Point_3
  //    + IndicesOutputIterator is an output iterator of value_type equal
  ///   to Curve_segment_index
  template <typename InputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  add_features(InputIterator first, InputIterator last,
               IndicesOutputIterator out /*= CGAL::Emptyset_iterator()*/);

  template <typename InputIterator>
  void 
  add_features(InputIterator first, InputIterator last)
  { add_features(first, last, CGAL::Emptyset_iterator()); }

  /// Insert one edge into domain
  /// InputIterator value type is Point_3
  template <typename InputIterator>
  Curve_segment_index insert_edge(InputIterator first, InputIterator last);

private:
  /// Returns the sign of the geodesic distance between \c p and \c q
  /// Precondition: index is not a cycle
  CGAL::Sign distance_sign(const Point_3& p, const Point_3& q,
                           const Curve_segment_index& index) const;  
  
  /// Returns Index associated to p (p must be the coordinates of a corner
  /// point)
  Index point_corner_index(const Point_3& p) const;
  
private:
  typedef std::map<Point_3,Corner_index> Corners;

  typedef internal::Mesh_3::Polyline<Gt> Polyline;
  typedef std::map<Curve_segment_index, Polyline> Edges;

  Corners corners_;
  Edges edges_;
  Curve_segment_index current_curve_index_;

private:
  // Disabled copy constructor & assignment operator
  typedef Mesh_domain_with_polyline_features_3 Self;
  Mesh_domain_with_polyline_features_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Mesh_domain_with_polyline_features_3



template <class MD_>
template <typename OutputIterator>
OutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_corners(OutputIterator out) const
{
  for ( typename Corners::const_iterator
       cit = corners_.begin(), end = corners_.end() ; cit != end ; ++cit )
  {
    *out++ = std::make_pair(cit->second,cit->first);
  }
  
  return out;
}

template <class MD_>
template <typename OutputIterator>
OutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_curve_segments(OutputIterator out) const
{
  for ( typename Edges::const_iterator
       eit = edges_.begin(), end = edges_.end() ; eit != end ; ++eit )
  {
    CGAL_assertion( eit->second.is_valid() );
    
    const Point_3& p = eit->second.start_point();
    const Point_3& q = eit->second.end_point();
    
    Index p_index, q_index;
    if ( ! eit->second.is_cycle() )
    {
      p_index = point_corner_index(p);
      q_index = point_corner_index(q);
    }
    else
    {
      p_index = index_from_curve_segment_index(eit->first);
      q_index = p_index;
    }
    
    *out++ = CGAL::cpp11::make_tuple(eit->first,
                                     std::make_pair(p,p_index),
                                     std::make_pair(q,q_index));
  }
  
  return out;
}
  

template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Index
Mesh_domain_with_polyline_features_3<MD_>::
point_corner_index(const Point_3& p) const
{
  typename Corners::const_iterator p_index_it = corners_.find(p);
  if ( p_index_it == corners_.end() )
  {
    CGAL_assertion(false);
    return Index();
  }
  
  return p_index_it->second;
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::FT
Mesh_domain_with_polyline_features_3<MD_>::
geodesic_distance(const Point_3& p, const Point_3& q,
                  const Curve_segment_index& curve_index) const
{
  // Get corresponding polyline
  typename Edges::const_iterator eit = edges_.find(curve_index);
  CGAL_assertion(eit != edges_.end());
  
  // Compute geodesic_distance
  return eit->second.geodesic_distance(p,q);
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Point_3
Mesh_domain_with_polyline_features_3<MD_>::
construct_point_on_curve_segment(const Point_3& starting_point,
                                 const Curve_segment_index& curve_index,
                                 FT distance) const
{
  // Get corresponding polyline
  typename Edges::const_iterator eit = edges_.find(curve_index);
  CGAL_assertion(eit != edges_.end());
  
  // Return point at geodesic_distance distance from starting_point
  return eit->second.point_at(starting_point,distance);
}



template <class MD_>
template <typename InputIterator, typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
add_features(InputIterator first, InputIterator last,
             IndicesOutputIterator indices_out)
{
  // Insert one edge for each element
  while ( first != last )
  {
    *indices_out++ = insert_edge(first->begin(), first->end());
    ++first;
  }
  return indices_out;
}



template <class MD_>
template <typename InputIterator>
typename Mesh_domain_with_polyline_features_3<MD_>::Curve_segment_index
Mesh_domain_with_polyline_features_3<MD_>::
insert_edge(InputIterator first, InputIterator last)
{
  CGAL_assertion(std::distance(first,last) > 1);
  
  const Curve_segment_index curve_index = current_curve_index_++;
  const Corner_index corner_index = 0;
  
  // Fill corners
  if ( *first != *boost::prior(last) )
  {
    corners_.insert(std::make_pair(*first,corner_index));
    corners_.insert(std::make_pair(*boost::prior(last),corner_index));
  }
  
  // Create a new polyline
  std::pair<typename Edges::iterator,bool> insertion =
    edges_.insert(std::make_pair(curve_index,Polyline()));
  
  // Fill polyline with data
  while ( first != last )
  {
    insertion.first->second.add_point(*first++);
  }
  return curve_index;
}


template <class MD_>
CGAL::Sign
Mesh_domain_with_polyline_features_3<MD_>::
distance_sign(const Point_3& p, const Point_3& q,
              const Curve_segment_index& index) const
{
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());
  CGAL_precondition( ! eit->second.is_cycle() );
  
  if ( p == q )
    return CGAL::ZERO;
  else if ( eit->second.are_ordered_along(p,q) )
    return CGAL::POSITIVE;
  else
    return CGAL::NEGATIVE;
}

template <class MD_>
CGAL::Sign
Mesh_domain_with_polyline_features_3<MD_>::
distance_sign_along_cycle(const Point_3& p,
                          const Point_3& q,
                          const Point_3& r,
                          const Curve_segment_index& index) const
{
  // Find edge
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());
  
  // If eit is not a cycle, then the orientation corresponds to the sign
  // of the distance
  if ( ! eit->second.is_cycle() )
  {
    return distance_sign(p,r,index);
  }
  
  // If p and r are the same point, it correspond to a complete loop on a cycle
  if ( p == r ) { return CGAL::POSITIVE; }
  
  // We are on a cycle without any clue (p==q). Return the shortest path as
  // orientation.
  if ( p == q )
  {
    FT pr = eit->second.geodesic_distance(p,r);
    FT rp = eit->second.geodesic_distance(r,p);
    if ( pr < rp ) { return CGAL::POSITIVE; }
    else { return CGAL::NEGATIVE; }
  }
  
  // If pq or pr is negative, edge is not a cycle, thus geodesic_distance
  // gives the answer.
  FT pq = eit->second.geodesic_distance(p,q);
  FT pr = eit->second.geodesic_distance(p,r);
  CGAL_assertion(pq > FT(0));
  CGAL_assertion(pr > FT(0));
  
  // Compare pq and pr 
  if ( pq <= pr ) { return CGAL::POSITIVE; }
  else { return CGAL::NEGATIVE; }
}

template <class MD_>
bool
Mesh_domain_with_polyline_features_3<MD_>::
is_cycle(const Point_3&, const Curve_segment_index& index) const
{
  // Find edge
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());
  
  return eit->second.is_cycle();
}


} //namespace CGAL


#endif // CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H
