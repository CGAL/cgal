// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Segment_Voronoi_diagram_kernel_wrapper_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H


#include <CGAL/Segment_Voronoi_diagram_site_2.h>



CGAL_BEGIN_NAMESPACE

template<class Kernel_base_2>
class Segment_Voronoi_diagram_kernel_wrapper_2
  : public Kernel_base_2
{
private:
  typedef Segment_Voronoi_diagram_kernel_wrapper_2<Kernel_base_2> Self;

public:
  //  typedef CGAL::Segment_Voronoi_diagram_site_2<Self>  Site_2;
  typedef CGAL::Segment_Voronoi_diagram_site_2<Kernel_base_2>  Site_2;
  //  typedef typename Site_2::Point_2    Point_2;
  //  typedef typename Site_2::Segment_2  Segment_2;

  //  typedef typename Kernel_base_2::Point_2    Point_2;
  //  typedef typename Kernel_base_2::Segment_2  Segment_2;
};




template<class K1, class K2, class Converter >
class Svd_cartesian_converter
  : public Converter
{
private:
  //  typedef Segment_Voronoi_diagram_kernel_wrapper_2<K1> K1W;
  //  typedef Segment_Voronoi_diagram_kernel_wrapper_2<K2> K2W;

  typedef typename K2::Site_2     K2_Site_2;
  typedef typename K2::Point_2    K2_Point_2;
  typedef typename K2::Segment_2  K2_Segment_2;

  typedef Converter               Base;

public:
  K2_Site_2
  operator()(const typename K1::Site_2& t) const
  {
    if ( t.is_point() ) {
      return K2_Site_2( Base::operator()(t.point()) );
    }

    return K2_Site_2( Base::operator()(t.segment()) );
  }


  K2_Point_2
  operator()(const typename K1::Point_2& p) const
  {
    return Base::operator()(p);
  }


  K2_Segment_2
  operator()(const typename K1::Segment_2& s) const
  {
    return  Base::operator()(s);
  }


  Sign
  operator()(const Sign& s) const
  {
    return s;
  }

};


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H
