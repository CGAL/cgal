// ======================================================================
//
// Copyright (c) 2003 GeometryFactory
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Level_interval.h
// package       : Interval_skip_list
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : GeometryFactory (<Andreas.Fabri@geometryfactory.com>)
//
// ======================================================================

#ifndef CGAL_LEVEL_INTERVAL_H
#define CGAL_LEVEL_INTERVAL_H

#include <iostream>


namespace CGAL {

  template <class FaceHandle>
  class Level_interval
  {
  public:
    typedef typename FaceHandle::value_type Face;
    typedef typename Face::Vertex Vertex;
    typedef typename Vertex::Point Point;
    typedef Kernel_traits<Point>::Kernel K;
    typedef typename K::FT Value;


  private:
    FaceHandle fh_; 
    Value inf_;
    Value sup_;  // left and right boundary values
  public:

    Level_interval(){}
    Level_interval(FaceHandle fh);
    const Value& inf() const {return inf_;}
    const Value& sup() const {return sup_;}
    FaceHandle face_handle() const { return fh_;}
    bool contains(const Value& V) const;

    // true iff this contains (l,r)
    bool contains_interval(const Value& l, const Value& r) const;  

    bool operator==(const Level_interval& I) const 
    {
      // there is no need to compare inf and sup, as these are derived from the face
      return face_handle() == I.face_handle();
    }

    bool operator!=(const Level_interval& I) const 
    {
      return face_handle() != I.face_handle();
    }
  };



  template <class V>
  std::ostream& operator<<(std::ostream& os, 
			   const Level_interval<V>& i)
  {
    os << i.face_handle()->vertex(0)->point() << ", " << 
      i.face_handle()->vertex(1)->point() << ", " << 
      i.face_handle()->vertex(2)->point() << std::endl; 
    return os;
  }


  template <class FaceHandle>
  Level_interval<FaceHandle>::Level_interval(FaceHandle fh)
    : fh_(fh), inf_(fh->vertex(0)->point().z()), sup_(inf_)
  {
    double z = fh->vertex(1)->point().z();
    sup_= (z>sup_)? z : sup_;
    inf_ = (z<inf_) ? z : inf_;
    z = fh->vertex(2)->point().z();
    sup_ = (z>sup_)? z : sup_;
    inf_ = (z<inf_) ? z : inf_;
  }


  template <class FaceHandle>
  bool
  Level_interval<FaceHandle>::contains_interval(const Value& i, 
					       const Value& s) const
    // true iff this contains (l,r)
  {
    return( (inf() <= i) && (sup() >= s) );
  }


  template <class FaceHandle>
  bool
  Level_interval<FaceHandle>::contains(const Value& v) const
  {
    // return true if this contains V, false otherwise
    if((v >= inf()) && (v <= sup()))
      return true;
    else
      return false;
  }




} // namespace CGAL

#endif // CGAL_LEVEL_INTERVAL_H







