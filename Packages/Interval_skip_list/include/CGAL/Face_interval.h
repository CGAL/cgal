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
// file          : include/CGAL/Face_interval.h
// package       : Interval_skip_list
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : GeometryFactory (<Andreas.Fabri@geometryfactory.com>)
//
// ======================================================================

#ifndef CGAL_FACE_INTERVAL_H
#define CGAL_FACE_INTERVAL_H

#include <iostream>


namespace CGAL {

  template <class FaceHandle>
  class Face_interval
  {
  public:
    // TODO: derive type from Ft of coordinates of points in vertex
    typedef double Value;

  private:
    FaceHandle fh_; 
    Value inf_;
    Value sup_;  // left and right boundary values
  public:

    Face_interval(){}
    Face_interval(FaceHandle fh);
    const Value& inf() const {return inf_;}
    const Value& sup() const {return sup_;}
    FaceHandle face_handle() const { return fh_;}
    bool contains(const Value& V) const;

    // true iff this contains (l,r)
    bool contains_interval(const Value& l, const Value& r) const;  

    bool operator==(const Face_interval& I) const 
    {
      return ( (inf() == I.inf()) && (sup() == I.sup()) && 
	       (face_handle() == I.face_handle()) );
    }

    bool operator!=(const Face_interval& I) const 
    {
      return ! (*this == I);
    }
  };



  template <class V>
  std::ostream& operator<<(std::ostream& os, 
			   const Face_interval<V>& i)
  {
    os << i.face_handle()->vertex(0)->point() << ", " << 
      i.face_handle()->vertex(1)->point() << ", " << 
      i.face_handle()->vertex(2)->point() << std::endl; 
    return os;
  }


  template <class FaceHandle>
  Face_interval<FaceHandle>::Face_interval(FaceHandle fh)
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
  Face_interval<FaceHandle>::contains_interval(const Value& i, 
					       const Value& s) const
    // true iff this contains (l,r)
  {
    return( (inf() <= i) && (sup() >= s) );
  }


  template <class FaceHandle>
  bool
  Face_interval<FaceHandle>::contains(const Value& v) const
  {
    // return true if this contains V, false otherwise
    if((v >= inf()) && (v <= sup()))
      return true;
    else
      return false;
  }




} // namespace CGAL

#endif // CGAL_FACE_INTERVAL_H







