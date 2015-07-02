// Copyright (c) 2015  Università della Svizzera italiana.
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
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_CONSTRUCTIONS_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_CONSTRUCTIONS_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>
#include <CGAL/enum.h>

#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>

#include <CGAL/Polychain_2.h>
#include <CGAL/intersections.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {


//***********************************************************************
//***********************************************************************
//                            CONSTRUCTIONS
//***********************************************************************
//***********************************************************************


//-----------------------------------------------------------------------
//                  Segment Delaunay graph site
//-----------------------------------------------------------------------
template<class Site,class ITag> class Construct_sdg_site_2;

template<class Site>
class Construct_sdg_site_2<Site,Tag_true>
{
public:
  typedef Site                             Site_2;
  typedef typename Site_2::Point_2         Point_2;
  typedef Site_2                           result_type;

public:
  result_type operator()(const Point_2& p) const {
    return Site_2(p);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1) const {
    return Site_2(p0, p1);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1,
                         const Point_2& q0, const Point_2& q1) const {
    return Site_2(p0, p1, q0, q1);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1,
                         const Point_2& q0, const Point_2& q1,
                         bool b) const {
    return Site_2(p0, p1, q0, q1, b);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1,
                         const Point_2& q0, const Point_2& q1,
                         const Point_2& r0, const Point_2& r1) const {
    return Site_2(p0, p1, q0, q1, r0, r1);
  }
};


template<class Site>
class Construct_sdg_site_2<Site,Tag_false>
{
public:
  typedef Site                             Site_2;
  typedef typename Site_2::Point_2         Point_2;
  typedef Site_2                           result_type;

public:
  result_type operator()(const Point_2& p) const {
    return Site_2(p);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1) const {
    return Site_2(p0, p1);
  }
};




//-----------------------------------------------------------------------
//                  Segment Delaunay graph Voronoi vertex
//-----------------------------------------------------------------------

template<class K, class M>
class Construct_svd_vertex_2
{
public:
  typedef typename K::Site_2                Site_2;
  typedef Voronoi_vertex_C2<K,M>            Voronoi_vertex_2;
  typedef typename K::Point_2               Point_2;
  typedef Point_2                           result_type;

public:
  Point_2 operator()(const Site_2& s1, const Site_2& s2,
                     const Site_2& s3) const
  {
    Voronoi_vertex_2 v(s1, s2, s3);
    return v.point();
  }
};



//-----------------------------------------------------------------------
//                  Segment Delaunay graph Voronoi circle
//-----------------------------------------------------------------------


template<class Gt, class M>
class Construct_sdg_circle_2
{
public:
  typedef typename Gt::Site_2                 Site_2;
  typedef Voronoi_vertex_C2<Gt,M>             Voronoi_vertex_2;
  typedef typename Gt::Iso_rectangle_2        Iso_rectangle_2;
  typedef Iso_rectangle_2                     result_type;

public:
  Iso_rectangle_2 operator() (const Site_2& s1, const Site_2& s2,
                       const Site_2& s3) const
  {
    Voronoi_vertex_2 v(s1, s2, s3);
    return v.circle();
  }
};



//-----------------------------------------------------------------------
//                    Segment Delaunay graph Voronoi bisector
//-----------------------------------------------------------------------


template<class Gt, class M>
class Construct_sdg_bisector_2
{
public:
  typedef typename Gt::Site_2                 Site_2;
  typedef typename Gt::Point_2                Point_2;
  typedef typename Gt::Direction_2            Direction_2;
  typedef typename Gt::Segment_2              Segment_2;
  typedef typename Gt::Line_2                 Line_2;

  typedef CGAL::Polychainline_2<Gt>           Polychainline;
  typedef Polychainline                       result_type;

  typedef typename Gt::Compare_x_2            Compare_x_2;
  typedef typename Gt::Compare_y_2            Compare_y_2;
  typedef typename Gt::Equal_2                Equal_2;

  typedef typename Gt::FT                     FT;

public:
  result_type
      operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    CGAL_SDG_DEBUG( std::cout << "debug construct bisector line "
              << "p=" << p << " q=" << q << std::endl; );

    if ( p.is_point() && q.is_point() ) {
      Point_2 pp = p.point();
      Point_2 pq = q.point();

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;
      CGAL_assertion_code ( Equal_2 are_same_points );
      CGAL_assertion( !(are_same_points(p, q)) );
      Comparison_result cmpx = compare_x_2(pp, pq);
      Comparison_result cmpy = compare_y_2(pp, pq);
      Comparison_result cmpabsdxy =
          CGAL::compare( CGAL::abs(pp.x()-pq.x()),
                         CGAL::abs(pp.y()-pq.y()) );
      unsigned int npts;
      Point_2 points[2];

      // (final) direction of bisector d=(-cmpy, cmpx)
      // the bisector should leave p to the right and q to the left
      Direction_2 d (
              (cmpy == EQUAL)? 0 :
              (  cmpy  == SMALLER )? +1 : -1,
              (cmpx == EQUAL)? 0 :
              (  cmpx  == SMALLER )? -1 : +1);

      // midpoint m of two points p and q
      Point_2 m = midpoint(pp, pq);

      if ((cmpabsdxy == EQUAL) || (cmpx == EQUAL) || (cmpy == EQUAL)) {
          // bisector is line going through m with direction d;
          // we will store this line as the union of two rays starting
          // at m with directions -d (incoming) and d (outgoing)
          npts = 1;
          points[0] = m;
      } else {
          // bisector consists of two rays and a middle segment;

          npts = 2;

          // compute length of middle segment
          FT half(0.5);
          FT seglenhalf ( half *
                  CGAL::abs(
                      CGAL::abs(pp.x()-pq.x()) -
                      CGAL::abs(pp.y()-pq.y()))   );

          // construct endpoints of middle segment of bisector
          Point_2 p1, p2;
          if (cmpabsdxy == SMALLER) {
              // middle segment is horizontal
              Point_2 p1temp (m.x() - seglenhalf, m.y());
              Point_2 p2temp (m.x() + seglenhalf, m.y());
              p1 = p1temp;
              p2 = p2temp;
          } else { // cmpabsdxy is LARGER
              // middle segment is vertical
              Point_2 p1temp (m.x(), m.y() - seglenhalf);
              Point_2 p2temp (m.x(), m.y() + seglenhalf);
              p1 = p1temp;
              p2 = p2temp;
          }

          // swap endpoints of segment if necessary
          if ( (cmpabsdxy == SMALLER ? cmpy : -cmpx) == LARGER ) {
              std::swap(p1, p2);
          }

          points[0] = p1;
          points[1] = p2;
      }

      Polychainline pcl(-d, points, points+npts, d);

      //CGAL_SDG_DEBUG( std::cout << "debug construct bisector line is "
      //  << pcl << std::endl; );

      return pcl;
    }// end of point - point case

    if ( (p.is_point() && q.is_segment())
        || (p.is_segment() && q.is_point()) ) {

      Point_2 pnt = (p.is_point()) ? p.point() : q.point();
      Segment_2 seg = (p.is_segment()) ? p.segment() : q.segment();

      Site_2 sitep = (p.is_point()) ? p : q;
      Site_2 sites = (p.is_point()) ? q : p;

      Point_2 points[3];
      unsigned int npts;

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;
      Equal_2     are_same_points;

      if (are_same_points(sitep,sites.source_site())
          || are_same_points(sitep,sites.target_site())) {
        npts = 1;
        points[0] = pnt;
        Point_2 pseg = (are_same_points(sitep,sites.source_site())) ?
          seg.target() : seg.source();

        Comparison_result cmpx = compare_x_2(pnt, pseg);
        Comparison_result cmpy = compare_y_2(pnt, pseg);

        Direction_2 d (
                       (cmpy == EQUAL)? 0 :
                       (  cmpy  == SMALLER )? +1 : -1,
                       (cmpx == EQUAL)? 0 :
                       (  cmpx  == SMALLER )? -1 : +1);
        if (q.is_point()) {
          d = -d;
        }

        Polychainline pcl(-d, points, points+npts, d);
        return pcl;
      } else {
        // pnt is not end point of seg
        // seg must not be horizontal or vertical
        CGAL_assertion(
            !( sites.supporting_site().segment().is_horizontal() ||
               sites.supporting_site().segment().is_vertical()     ) );
        Line_2 lseg = sites.supporting_site().segment().supporting_line();
        Point_2 phor, pver;
        phor = Point_2(lseg.x_at_y(pnt.y()), pnt.y());
        pver = Point_2(pnt.x(), lseg.y_at_x(pnt.x()));
        //pfirst and plast are points on the supporting line of seg
        Point_2 pfirst, plast;
        //pcfirst and pclast are points on the bisector
        Point_2 pcfirst, pclast;
        // segment with positive slope will have pfirst as phor
        // segment with negative slope will have pfirst as pver
        pfirst = (compare_x_2(seg.source(),seg.target())
                  == compare_y_2(seg.source(),seg.target()))
                 ? phor : pver;
        plast = (compare_x_2(seg.source(),seg.target())
                 == compare_y_2(seg.source(),seg.target()))
                ? pver : phor;

        FT half = FT(0.5);
        Point_2 pmid_pfirst_pnt = midpoint(pfirst, pnt);
        Point_2 pmid_plast_pnt = midpoint(plast, pnt);
        FT seglenhalffirst ( half *
                            CGAL::abs(
                                      CGAL::abs(pnt.x()-pfirst.x()) -
                                      CGAL::abs(pnt.y()-pfirst.y()))   );
        FT seglenhalflast ( half *
                           CGAL::abs(
                                     CGAL::abs(pnt.x()-plast.x()) -
                                     CGAL::abs(pnt.y()-plast.y()))   );

        if (compare_x_2(seg.source(),seg.target())
            == compare_y_2(seg.source(),seg.target())) {
          //segment with positive slope
          if ( (compare_x_2(seg.source(),seg.target()) == SMALLER
                && lseg.has_on_positive_side(pnt))
              || (compare_x_2(seg.source(),seg.target()) == LARGER
                  && lseg.has_on_negative_side(pnt)) ) {
                //pcfirst is center of square,
                //pfirst = phor, upward direction
                //pclast is center of sqaure, plast = pver, left direction
                pcfirst = Point_2(pmid_pfirst_pnt.x(),
                                  pmid_pfirst_pnt.y()+seglenhalffirst);
                pclast = Point_2(pmid_plast_pnt.x()-seglenhalflast,
                                 pmid_plast_pnt.y());
              } else {
                //pfirst = phor , pcfirst in downward direction
                //plast = pvor , pclast in right direction
                pcfirst = Point_2(pmid_pfirst_pnt.x(),
                                  pmid_pfirst_pnt.y()-seglenhalffirst);
                pclast = Point_2(pmid_plast_pnt.x()+seglenhalflast,
                                 pmid_plast_pnt.y());
              }
        }
        else {
          //segment with negative slope
          if ( (compare_x_2(seg.source(),seg.target()) == SMALLER
                && lseg.has_on_positive_side(pnt))
              || (compare_x_2(seg.source(),seg.target()) == LARGER
                  && lseg.has_on_negative_side(pnt)) ) {
                //pcfirst is center of square,
                //pfirst = pver, right direction
                //pclast is center of sqaure, plast = phor, upward direction
                pcfirst = Point_2(pmid_pfirst_pnt.x()+seglenhalffirst,
                                  pmid_pfirst_pnt.y());
                pclast = Point_2(pmid_plast_pnt.x(),
                                 pmid_plast_pnt.y()+seglenhalflast);
              } else {
                //pfirst = pver , pcfirst in left direction
                //plast = phor , pclast in downward direction
                pcfirst = Point_2(pmid_pfirst_pnt.x()-seglenhalffirst,
                                  pmid_pfirst_pnt.y());
                pclast = Point_2(pmid_plast_pnt.x(),
                                 pmid_plast_pnt.y()-seglenhalflast);
              }
        }//end of pcfirst and pclast

        //compute pmid and then pcmid = mid point of pmid and pnt
        Line_2 lmid(pnt,pcfirst);
        Line_2 lmidp = lmid.perpendicular(pnt);
        CGAL::Object pmidobject = intersection(lmidp, lseg);
        Point_2 pmid;
        if(CGAL::assign(pmid, pmidobject)){
          points[1] = midpoint(pmid, pnt);
        }
        npts = 3;
        points[0]=pcfirst;
        points[2]=pclast;

        CGAL_SDG_DEBUG( std::cout << "debug: point1 = " << points[0]
            << " point2 = " << points[1] << " point3 = "
            << points[2] << std::endl; );

        if (p.is_segment()) {
          std::swap(points[0], points[2]);
        }
        CGAL_SDG_DEBUG( std::cout << "debug: point1 = " << points[0]
            << " point2 = " << points[1] << " point3 = "
            << points[2] << std::endl; );

        Line_2 ld(pnt,pcfirst);

        Direction_2 d(ld.perpendicular(pcfirst));

        Polychainline pcl(d, points, points+npts, d);
        CGAL_SDG_DEBUG( std::cout
            << "debug about to return pcl=" << pcl << std::endl ;);
        return pcl;

      }
    }
    // this part should never be reached
    // philaris: avoid complaining by some compilers
    return Polychainline();
  }
};

//-----------------------------------------------------------------------
//                 Segment Delaunay graph Voronoi bisecting ray
//-----------------------------------------------------------------------

template<class Gt, class M>
class Construct_sdg_bisector_ray_2
{
public:
  typedef typename Gt::Site_2                  Site_2;
  typedef typename Gt::Point_2                 Point_2;
  typedef typename Gt::Direction_2             Direction_2;
  typedef typename Gt::Line_2                  Line_2;
  typedef typename Gt::Segment_2               Segment_2;
  typedef typename Gt::Construct_svd_vertex_2  Construct_svd_vertex_2;
  typedef CGAL::Polychainline_2<Gt>            Polychainline;
  typedef CGAL::Polychainray_2<Gt>             Polychainray;
  typedef Polychainray                         result_type;

  typedef typename Gt::Compare_x_2       Compare_x_2;
  typedef typename Gt::Compare_y_2       Compare_y_2;
  typedef typename Gt::Equal_2           Equal_2;

  typedef typename Gt::FT           FT;

  result_type
      operator()(const Site_2& p, const Site_2& q,
                   const Site_2& r) const
  {

    CGAL_SDG_DEBUG( std::cout << "debug construct bisector ray "
              << "p=" << p << " q=" << q << " r=" << r << std::endl; );

    // compute pqr vertex
    Point_2 v = Construct_svd_vertex_2()(p, q, r);

    CGAL_SDG_DEBUG( std::cout << "debug construct bisector ray "
              << "p=" << p << " q=" << q << " r=" << r
              << " has v(pqr)=" << v << std::endl; );

    //CGAL_SDG_DEBUG( std::cout
    //    << "debug in ray computing bisector" << std::endl; );

    // compute oriented bisector between p and q

    if ( p.is_point() && q.is_point() ) {
      Point_2 pp = p.point();
      Point_2 pq = q.point();

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;
      CGAL_assertion_code ( Equal_2 are_same_points );
      CGAL_assertion( !(are_same_points(p, q)) );
      Comparison_result cmpx = compare_x_2(pp, pq);
      Comparison_result cmpy = compare_y_2(pp, pq);
      Comparison_result cmpabsdxy =
          CGAL::compare( CGAL::abs(pp.x()-pq.x()),
                         CGAL::abs(pp.y()-pq.y()) );
      unsigned int npts;
      Point_2 points[3];

      // (final) direction of bisector d=(-cmpy, cmpx)
      // the bisector should leave p to the right and q to the left
      Direction_2 d (
              (cmpy == EQUAL)? 0 :
              (  cmpy  == SMALLER )? +1 : -1,
              (cmpx == EQUAL)? 0 :
              (  cmpx  == SMALLER )? -1 : +1);

      // bisector ray always starts with v
      points[0] = v;

      if ((cmpabsdxy == EQUAL) || (cmpx == EQUAL) || (cmpy == EQUAL)) {
          // bisector ray consists of a ray starting at v
          npts = 1;
      } else {
          // midpoint m of two points p and q
          Point_2 m = midpoint(pp, pq);

          // compute length of middle segment
          FT half(0.5);
          FT seglenhalf ( half *
                  CGAL::abs(
                      CGAL::abs(pp.x()-pq.x()) -
                      CGAL::abs(pp.y()-pq.y()))   );

          // construct endpoints of middle segment of bisector
          Point_2 p1, p2;
          if (cmpabsdxy == SMALLER) {
              // middle segment is horizontal
              Point_2 p1temp (m.x() - seglenhalf, m.y());
              Point_2 p2temp (m.x() + seglenhalf, m.y());
              p1 = p1temp;
              p2 = p2temp;
          } else { // cmpabsdxy is LARGER
              // middle segment is vertical
              Point_2 p1temp (m.x(), m.y() - seglenhalf);
              Point_2 p2temp (m.x(), m.y() + seglenhalf);
              p1 = p1temp;
              p2 = p2temp;
          }

          // swap endpoints of segment if necessary
          if ( (cmpabsdxy == SMALLER ? cmpy : -cmpx) == LARGER ) {
              std::swap(p1, p2);
          }

          // oriented line from p1 to p2
          Line_2 l(p1, p2);

          if (l.perpendicular(p1).has_on_positive_side(v)) {
              npts = 3;
              points[1] = p1;
              points[2] = p2;
          } else if (l.perpendicular(p2).has_on_positive_side(v)) {
              npts = 2;
              points[1] = p2;
          } else {
              npts = 1;
          }

      }

      Polychainray pcr(points, points+npts, d);

      CGAL_SDG_DEBUG( std::cout
          << "debug construct bisector ray is " << pcr << std::endl; );

      return pcr;
    }

    if ( (p.is_point() && q.is_segment())
        || (p.is_segment() && q.is_point()) ) {

      Point_2 pnt = (p.is_point()) ? p.point() : q.point();
      Segment_2 seg = (p.is_segment()) ? p.segment() : q.segment();

      Site_2 sitep = (p.is_point()) ? p : q;
      Site_2 sites = (p.is_point()) ? q : p;

      Point_2 points[4];
      unsigned int npts;

      Compare_x_2 compare_x_2;
      Compare_y_2 compare_y_2;
      Equal_2     are_same_points;

      if (are_same_points(sitep,sites.source_site())
          || are_same_points(sitep,sites.target_site())) {
        npts = 1;
        points[0] = v;
        Point_2 pseg =
          (are_same_points(sitep,sites.source_site())) ?
            seg.target() : seg.source();

        Comparison_result cmpx = compare_x_2(pnt, pseg);
        Comparison_result cmpy = compare_y_2(pnt, pseg);

        Direction_2 d (
                       (cmpy == EQUAL)? 0 :
                       (  cmpy  == SMALLER )? +1 : -1,
                       (cmpx == EQUAL)? 0 :
                       (  cmpx  == SMALLER )? -1 : +1);

        if (q.is_point()) {
          d = -d;
        }

        Polychainray pcr(points, points+npts, d);

        CGAL_SDG_DEBUG( std::cout
            << "debug construct bisector ray is " << pcr << std::endl; );

        return pcr;

      } else {
        // pnt is not end point of seg
        // seg must not be horizontal or vertical
        CGAL_assertion(
            !( sites.supporting_site().segment().is_horizontal() ||
               sites.supporting_site().segment().is_vertical()     ) );
        // bisector ray always starts with v
        points[0] = v;
        Line_2 lseg = sites.supporting_site().segment().supporting_line();
        Point_2 phor, pver;
        phor = Point_2(lseg.x_at_y(pnt.y()), pnt.y());
        pver = Point_2(pnt.x(), lseg.y_at_x(pnt.x()));
        //pfirst and plast are points on the supporting line of seg
        Point_2 pfirst, plast;
        //pcfirst and pclast are points on the bisector
        Point_2 pcfirst, pclast;
        // segment with positive slope will have pfirst as phor
        // segment with negative slope will have pfirst as pver
        pfirst = (compare_x_2(seg.source(),seg.target())
                  == compare_y_2(seg.source(),seg.target()))
        ? phor : pver;
        plast = (compare_x_2(seg.source(),seg.target())
                 == compare_y_2(seg.source(),seg.target()))
        ? pver : phor;

        FT half = FT(0.5);
        Point_2 pmid_pfirst_pnt = midpoint(pfirst, pnt);
        Point_2 pmid_plast_pnt = midpoint(plast, pnt);
        FT seglenhalffirst ( half *
                            CGAL::abs(
                                      CGAL::abs(pnt.x()-pfirst.x()) -
                                      CGAL::abs(pnt.y()-pfirst.y()))   );
        FT seglenhalflast ( half *
                           CGAL::abs(
                                     CGAL::abs(pnt.x()-plast.x()) -
                                     CGAL::abs(pnt.y()-plast.y()))   );

        if (compare_x_2(seg.source(),seg.target())
            == compare_y_2(seg.source(),seg.target())) {
          //segment with positive slope
          if ( (compare_x_2(seg.source(),seg.target()) == SMALLER
                && lseg.has_on_positive_side(pnt))
              || (compare_x_2(seg.source(),seg.target()) == LARGER
                  && lseg.has_on_negative_side(pnt)) ) {
            //pcfirst is center of square , pfirst = phor, upward direction
            //pclast is center of sqaure, plast = pver, left direction
            pcfirst = Point_2(pmid_pfirst_pnt.x(),
                              pmid_pfirst_pnt.y()+seglenhalffirst);
            pclast = Point_2(pmid_plast_pnt.x()-seglenhalflast,
                              pmid_plast_pnt.y());
          } else {
            //pfirst = phor , pcfirst in downward direction
            //plast = pvor , pclast in right direction
            pcfirst = Point_2(pmid_pfirst_pnt.x(),
                              pmid_pfirst_pnt.y()-seglenhalffirst);
            pclast = Point_2(pmid_plast_pnt.x()+seglenhalflast,
                              pmid_plast_pnt.y());
          }
        }
        else {
          //segment with negative slope
          if ( (compare_x_2(seg.source(),seg.target()) == SMALLER
                && lseg.has_on_positive_side(pnt))
              || (compare_x_2(seg.source(),seg.target()) == LARGER
                  && lseg.has_on_negative_side(pnt)) ) {
            //pcfirst is center of square , pfirst = pver, right direction
            //pclast is center of sqaure, plast = phor, upward direction
            pcfirst = Point_2(pmid_pfirst_pnt.x()+seglenhalffirst,
                              pmid_pfirst_pnt.y());
            pclast = Point_2(pmid_plast_pnt.x(),
                              pmid_plast_pnt.y()+seglenhalflast);
          } else {
            //pfirst = pver , pcfirst in left direction
            //plast = phor , pclast in downward direction
            pcfirst = Point_2(pmid_pfirst_pnt.x()-seglenhalffirst,
                              pmid_pfirst_pnt.y());
            pclast = Point_2(pmid_plast_pnt.x(),
                              pmid_plast_pnt.y()-seglenhalflast);
          }
        }//end of pcfirst and pclast

        //compute pmid and then pcmid = mid point of pmid and pnt
        Line_2 lmid(pnt,pcfirst);
        Line_2 lmidp = lmid.perpendicular(pnt);
        CGAL::Object pmidobject = intersection(lmidp, lseg);
        Point_2 pmid;
        if(CGAL::assign(pmid, pmidobject)){
          points[2] = midpoint(pmid, pnt);
        }

        points[1]=pcfirst;
        points[3]=pclast;

        CGAL_SDG_DEBUG( std::cout << "debug: point1 = " << points[1] <<
        " point2 = " << points[2] << " point3 = "
        << points[3] << std::endl; );

        if (p.is_segment()) {
          std::swap(points[1], points[3]);
        }
        CGAL_SDG_DEBUG( std::cout << "debug: point1 = " << points[1] <<
        " point2 = " << points[2] << " point3 = "
        << points[3] << std::endl; );

        //oriented line from pcfirst to pnt
        Line_2 l(points[1], pnt);

        Direction_2 d(l.perpendicular(pnt));
        if (p.is_point()) {
          d = -d;
        }

        if (l.perpendicular(pnt).has_on_positive_side(v)) {
          if ( (p.is_segment() && l.has_on_positive_side(v))
              ||(p.is_point() && l.has_on_negative_side(v)) ) {
            npts = 4;
          } else {
            npts = 3;
            points[1] = points[2];
            points[2] = points[3];
          }
        }
        else {
          if ( (p.is_segment() && l.has_on_negative_side(v))
              ||(p.is_point() && l.has_on_positive_side(v)) ) {
            npts = 2;
            points[1] = points[3];
          }
          else {
            npts = 1;
          }
        }

        Polychainray pcr(points, points+npts, d);
        CGAL_SDG_DEBUG( std::cout << "debug construct bisector ray is "
            << pcr << std::endl; );
        return pcr;

      }

    }

    // the following code should never be executed
    CGAL_SDG_DEBUG( std::cout
        << "debug construct bisector ray error " << std::endl; );
    CGAL_assertion(false);
    return Polychainray();
  }
};


//-----------------------------------------------------------------------
//              Segment Delaunay graph Voronoi bisecting segment
//-----------------------------------------------------------------------


template<class Gt, class M>
class Construct_sdg_bisector_segment_2
{
public:
  typedef typename Gt::Site_2                   Site_2;
  typedef typename Gt::Point_2                  Point_2;
  typedef typename Gt::Direction_2              Direction_2;
  typedef typename Gt::Line_2                   Line_2;
  typedef typename Gt::Segment_2                Segment_2;
  typedef typename Gt::Construct_svd_vertex_2   Construct_svd_vertex_2;
  typedef CGAL::Polychainsegment_2<Gt>          Polychainsegment;

  typedef CGAL::Object                          Object_2;
  typedef Object_2                              result_type;

  typedef typename Gt::Compare_x_2       Compare_x_2;
  typedef typename Gt::Compare_y_2       Compare_y_2;
  typedef typename Gt::Equal_2           Equal_2;

  typedef typename Gt::FT                FT;

  result_type operator()(const Site_2& p, const Site_2& q,
                         const Site_2& r, const Site_2& s) const
  {
    return CGAL::make_object(pcsbis(p, q, r, s));
  }

  Polychainsegment pcsbis(const Site_2& p, const Site_2& q,
                          const Site_2& r, const Site_2& s) const
  {
    CGAL_SDG_DEBUG( std::cout << "debug construct bisector segment "
        << "p=" << p << " q=" << q << " r=" << r << " s=" << s
        << std::endl; );

    Construct_svd_vertex_2 circumcenter;
    Point_2 vpqr = circumcenter(p, q, r);
    Point_2 vqps = circumcenter(q, p, s);

    CGAL_SDG_DEBUG( std::cout << "debug construct bisector segment "
        << "p=" << p << " q=" << q << " r=" << r << " s=" << s
        << " has v(pqr)=" << vpqr << " v(qps)=" << vqps << std::endl; );

    Compare_x_2 compare_x_2;
    Compare_y_2 compare_y_2;

    Comparison_result cmpx_vpqr_vqps =
      compare_x_2(vpqr, vqps);

    Comparison_result cmpy_vpqr_vqps =
      compare_y_2(vpqr, vqps);

    CGAL_SDG_DEBUG( std::cout << "debug bis segment vpqr="
      << vpqr << " vqps=" << vqps
      << " cmpx=" << cmpx_vpqr_vqps
      << " cmpy=" << cmpy_vpqr_vqps
      << std::endl; );

    if ( (cmpx_vpqr_vqps == EQUAL) &&
         (cmpy_vpqr_vqps == EQUAL)    ) {
      // philaris: here vpqr and vqps are the same point,
      // therefore the bisector segment is (vpqr,vqps),
      // which is trivial
      Point_2 points[2];
      points[0] = vpqr;
      points[1] = vqps;
      Polychainsegment pcs(points, points+2);

      CGAL_SDG_DEBUG( std::cout
          << "debug construct bisector segment is (trivial) "
          << pcs << std::endl; );
      return pcs;
    }

    if ( p.is_point() && q.is_point() ) {
      Point_2 pp = p.point();
      Point_2 pq = q.point();

      CGAL_assertion_code( Equal_2 are_same_points );
      CGAL_assertion( ! (are_same_points(p, q)) );
      Comparison_result cmpx = compare_x_2(pp, pq);
      Comparison_result cmpy = compare_y_2(pp, pq);
      Comparison_result cmpabsdxy =
          CGAL::compare( CGAL::abs(pp.x()-pq.x()),
                         CGAL::abs(pp.y()-pq.y()) );
      unsigned int npts;
      Point_2 points[4];

      // (final) direction of bisector d=(-cmpy, cmpx)
      // the bisector should leave p to the right and q to the left
      Direction_2 d (
              (cmpy == EQUAL)? 0 :
              (  cmpy  == SMALLER )? +1 : -1,
              (cmpx == EQUAL)? 0 :
              (  cmpx  == SMALLER )? -1 : +1);

      // bisector segment always starts with vpqr
      points[0] = vpqr;

      if ((cmpabsdxy == EQUAL) || (cmpx == EQUAL) || (cmpy == EQUAL)) {
        // bisector ray is line segment (vpqr,vqps)
        npts = 2;
      } else {
        // midpoint m of two points p and q
        Point_2 m = midpoint(pp, pq);

        // compute length of middle segment
        FT half(0.5);
        FT seglenhalf ( half *
            CGAL::abs(
              CGAL::abs(pp.x()-pq.x()) -
              CGAL::abs(pp.y()-pq.y()))   );

        // construct endpoints of middle segment of bisector
        Point_2 p1, p2;
        if (cmpabsdxy == SMALLER) {
          // middle segment is horizontal
          Point_2 p1temp (m.x() - seglenhalf, m.y());
          Point_2 p2temp (m.x() + seglenhalf, m.y());
          p1 = p1temp;
          p2 = p2temp;
        } else { // cmpabsdxy is LARGER
          // middle segment is vertical
          Point_2 p1temp (m.x(), m.y() - seglenhalf);
          Point_2 p2temp (m.x(), m.y() + seglenhalf);
          p1 = p1temp;
          p2 = p2temp;
        }

        // swap endpoints of segment if necessary
        if ( (cmpabsdxy == SMALLER ? cmpy : -cmpx) == LARGER ) {
          std::swap(p1, p2);
        }

        // oriented line from p1 to p2
        Line_2 l(p1, p2);

        if (l.perpendicular(p1).has_on_positive_side(vpqr)) {
          npts = 4;
          points[1] = p1;
          points[2] = p2;
        } else if (l.perpendicular(p2).has_on_positive_side(vpqr)) {
          npts = 3;
          points[1] = p2;
        } else {
          npts = 2;
        }

        CGAL_SDG_DEBUG( std::cout << "debug p1=" << p1 << " p2=" << p2
          << " vpqr=" << vpqr << " vqps=" << vqps << std::endl; );

        CGAL_assertion(l.perpendicular(vpqr).
                         has_on_negative_side(vqps));

        while (npts > 2) {
          if (l.perpendicular(points[npts-2]).
                has_on_negative_side(vqps)) {
            break;
          } else {
            npts = npts - 1;
          }
        }
      }

      // bisector segment always ends with vpqr
      points[npts-1] = vqps;

      Polychainsegment pcs(points, points+npts);

      //CGAL_SDG_DEBUG( std::cout << "debug construct bisector segment is "
      //  << pcs << std::endl; );

      return pcs;
    } // end of two points case
    else if ( (p.is_point() && q.is_segment()) ||
              (q.is_point() && p.is_segment()) ) {
      //one site is point and one site is segment
      unsigned int npts;
      Point_2 points[5];
      points[0] = vpqr;
      Equal_2 are_same_points;
      // check if p is point and is an endpoint of q
      // or if q is a point and is an end point of p
      if (   ( p.is_point() &&
              (are_same_points(p, q.source_site()) ||
               are_same_points(p, q.target_site())   ) )
          || (q.is_point() &&
              (are_same_points(q, p.source_site()) ||
               are_same_points(q, p.target_site())   ) ) )
      {
        npts = 2;
      }
      else {
        CGAL_SDG_DEBUG(std::cout
            << "debug bisector segment "
            << "p=" << p << " q=" << q << std::endl; );
        //pnt is the point site and seg is the segment site
        Point_2 pnt = (p.is_point()) ? p.point() : q.point();
        Segment_2 seg = (p.is_segment()) ? p.segment() : q.segment();
        Site_2 siteseg = (p.is_point()) ? q : p;
        // lseg is the suporting line of the segment site
        Line_2 lseg = siteseg.supporting_site().segment().supporting_line();
        // segment site is horizontal
        if (lseg.is_horizontal()) {
        //pver is vertical projection from point site on to segment site
          Point_2 pver;
          pver = Point_2(pnt.x(), lseg.y_at_x(pnt.x()));

          Point_2 m = midpoint(pnt, pver);
          FT half = FT(0.5);
          FT seglenhalf ( half * CGAL::abs(pnt.y()-pver.y()) );

          Comparison_result cmp = compare_y_2(pnt, pver);
          if (cmp == LARGER) {
            points[1] = Point_2(m.x() + seglenhalf, m.y());
            points[2] = Point_2(m.x() - seglenhalf, m.y());
          } else {
            CGAL_assertion(cmp == SMALLER);
            points[1] = Point_2(m.x() - seglenhalf, m.y());
            points[2] = Point_2(m.x() + seglenhalf, m.y());
          }

          if (p.is_segment()) {
            std::swap(points[1], points[2]);
          }
          // oriented line from p1 to p2
          Line_2 l(points[1], points[2]);

          if (l.perpendicular(points[1]).has_on_positive_side(vpqr)) {
            npts = 4;
            //points[1] = p1;
            //points[2] = p2;
          } else if (l.perpendicular(points[2]).
                     has_on_positive_side(vpqr)) {
            npts = 3;
            points[1] = points[2];
          } else {
            npts = 2;
          }

          CGAL_SDG_DEBUG(std::cout
              << "debug construct bisector segment npts="
              << npts << std::endl; );

          /*CGAL_assertion((l.perpendicular(vpqr).
                         has_on_negative_side(vqps)) ||
                         ((compare_x_2(vpqr, vqps) == EQUAL) and
                          (compare_y_2(vpqr, vqps) == EQUAL))
                         );*/

          while (npts > 2) {
            if (l.perpendicular(points[npts-2]).
                has_on_negative_side(vqps)) {
              break;
            } else {
              npts = npts - 1;
            }
          }
          CGAL_SDG_DEBUG( std::cout
              << "debug bisector segment after cutting points, npts="
              << npts << std::endl; );
        }//end of horizontal segment case
        else if (lseg.is_vertical()) {
          //segment site is vertical
          // phor is the projection of pnt on seg
          Point_2 phor;
          phor = Point_2(lseg.x_at_y(pnt.y()), pnt.y());

          Point_2 m = midpoint(pnt, phor);
          FT half = FT(0.5);
          FT seglenhalf ( half * CGAL::abs(pnt.x()-phor.x()) );

          npts = 4;

          Comparison_result cmp = compare_x_2(pnt, phor);
          if (cmp == LARGER) {
            points[1] = Point_2(m.x(), m.y() - seglenhalf);
            points[2] = Point_2(m.x(), m.y() + seglenhalf);
          } else {
            points[1] = Point_2(m.x(), m.y() + seglenhalf);
            points[2] = Point_2(m.x(), m.y() - seglenhalf);
          }

          if (p.is_segment()) {
            std::swap(points[1], points[2]);
          }
          Line_2 l(points[1], points[2]);

          if (l.perpendicular(points[1]).has_on_positive_side(vpqr)) {
            npts = 4;
            //points[1] = p1;
            //points[2] = p2;
          } else if (l.perpendicular(points[2]).
                     has_on_positive_side(vpqr)) {
            npts = 3;
            points[1] = points[2];
          } else {
            npts = 2;
          }

          CGAL_SDG_DEBUG( std::cout
              << "debug bis segment npts=" << npts << " vpqr="
              << vpqr << " vqps=" << vqps << std::endl; );


          // philaris: assertion does not work always
          // because of inexact arithmetic
          //CGAL_assertion(l.perpendicular(vpqr).
          //                has_on_negative_side(vqps));

          while (npts > 2) {
            if (l.perpendicular(points[npts-2]).
                has_on_negative_side(vqps)) {
              break;
            } else {
              npts = npts - 1;
            }
          }
        }// end of the vertical segment case
        else {//the segment is neither horizontal nor vertical
          Point_2 phor,pver;
          phor = Point_2(lseg.x_at_y(pnt.y()), pnt.y());
          pver = Point_2(pnt.x(), lseg.y_at_x(pnt.x()));
          //pfirst and plast are points on the supporting line of seg
          Point_2 pfirst, plast;
          //pcfirst and pclast are points on the bisector
          Point_2 pcfirst, pclast;

          // segment with positive slope will have pfirst as phor
          // segment with negative slope will have pfirst as pver
          pfirst = (compare_x_2(seg.source(),seg.target())
                    == compare_y_2(seg.source(),seg.target()))
                    ? phor : pver;
          plast = (compare_x_2(seg.source(),seg.target())
                    == compare_y_2(seg.source(),seg.target()))
                    ? pver : phor;

          FT half = FT(0.5);
          Point_2 pmid_pfirst_pnt = midpoint(pfirst, pnt);
          Point_2 pmid_plast_pnt = midpoint(plast, pnt);
          FT seglenhalffirst ( half *
                              CGAL::abs(
                                      CGAL::abs(pnt.x()-pfirst.x()) -
                                      CGAL::abs(pnt.y()-pfirst.y()))   );
          FT seglenhalflast ( half *
                                CGAL::abs(
                                      CGAL::abs(pnt.x()-plast.x()) -
                                      CGAL::abs(pnt.y()-plast.y()))   );

          if (compare_x_2(seg.source(),seg.target())
              == compare_y_2(seg.source(),seg.target())) {
            //segment with positive slope
            if ( (compare_x_2(seg.source(),seg.target()) == SMALLER
                  && lseg.has_on_positive_side(pnt))
                || (compare_x_2(seg.source(),seg.target()) == LARGER
                  && lseg.has_on_negative_side(pnt)) ) {
              //pcfirst is center of square, pfirst = phor, upward direction
              //pclast is center of sqaure, plast = pver, left direction
              pcfirst = Point_2(pmid_pfirst_pnt.x(),
                                pmid_pfirst_pnt.y()+seglenhalffirst);
              pclast = Point_2(pmid_plast_pnt.x()-seglenhalflast,
                                pmid_plast_pnt.y());
            } else {
              //pfirst = phor , pcfirst in downward direction
              //plast = pvor , pclast in right direction
              pcfirst = Point_2(pmid_pfirst_pnt.x(),
                                pmid_pfirst_pnt.y()-seglenhalffirst);
              pclast = Point_2(pmid_plast_pnt.x()+seglenhalflast,
                               pmid_plast_pnt.y());
            }
          }
          else {
            //segment with negative slope
            if ( (compare_x_2(seg.source(),seg.target()) == SMALLER
                  && lseg.has_on_positive_side(pnt))
                || (compare_x_2(seg.source(),seg.target()) == LARGER
                  && lseg.has_on_negative_side(pnt)) ) {
              //pcfirst is center of square , pfirst = pver, right direction
              //pclast is center of sqaure, plast = phor, upward direction
              pcfirst = Point_2(pmid_pfirst_pnt.x()+seglenhalffirst,
                                pmid_pfirst_pnt.y());
              pclast = Point_2(pmid_plast_pnt.x(),
                               pmid_plast_pnt.y()+seglenhalflast);
            } else {
              //pfirst = pver , pcfirst in left direction
              //plast = phor , pclast in downward direction
              pcfirst = Point_2(pmid_pfirst_pnt.x()-seglenhalffirst,
                                pmid_pfirst_pnt.y());
              pclast = Point_2(pmid_plast_pnt.x(),
                               pmid_plast_pnt.y()-seglenhalflast);
            }
          }//end of pcfirst and pclast

          //compute pmid and then pcmid = mid point of pmid and pnt
          Line_2 lmid(pnt,pcfirst);
          Line_2 lmidp = lmid.perpendicular(pnt);
          CGAL::Object pmidobject = intersection(lmidp, lseg);
          Point_2 pmid;
          if(CGAL::assign(pmid, pmidobject)){
            points[2] = midpoint(pmid, pnt);
          }

          points[1]=pcfirst;
          points[3]=pclast;

          CGAL_SDG_DEBUG( std::cout << "debug: point1 = " << points[1] <<
          " point2 = " << points[2] << " point3 = "
          << points[3] << std::endl; );

          if (p.is_segment()) {
            std::swap(points[1], points[3]);
          }

          CGAL_SDG_DEBUG( std::cout << "debug: after possible swap: "
                    << "point1 = " << points[1] <<
                    " point2 = " << points[2] << " point3 = "
                    << points[3] << std::endl; );

          //oriented line from pcfirst to pnt
          Line_2 l(points[1], pnt);

          // philaris: probably not needed, but check more
          //CGAL_assertion((l.perpendicular(vpqr).
          //                has_on_negative_side(vqps)) ||
          //               (vpqr == vqps));

          if (l.perpendicular(pnt).has_on_positive_side(vpqr)) {
            if ( (p.is_segment() && l.has_on_positive_side(vpqr))
                ||(p.is_point() && l.has_on_negative_side(vpqr)) ) {
              npts = 5;
            } else {
              npts = 4;
              points[1] = points[2];
              points[2] = points[3];
            }
          }
          else {
            if ( (p.is_segment() && l.has_on_negative_side(vpqr))
                ||(p.is_point() && l.has_on_positive_side(vpqr)) ) {
              npts = 3;
              points[1] = points[3];
            }
            else {
              npts = 2;
            }
          }

          CGAL_SDG_DEBUG( std::cout << "debug first check npts="
              << npts << std::endl; );

          if ((npts == 3) &&
              ((p.is_segment() && ! l.has_on_positive_side(vqps))
               ||(p.is_point() && ! l.has_on_negative_side(vqps)))) {
            npts = 2;
          }
          if (npts == 4) {
            if (! l.perpendicular(pnt).has_on_negative_side(vqps)) {
              npts = 2;
            } else {
              if ((p.is_segment() && ! l.has_on_positive_side(vqps))
                  ||(p.is_point() && ! l.has_on_negative_side(vqps))) {
                npts = 3;
              }
            }
          }
          if (npts == 5) {
            if (l.perpendicular(pnt).has_on_negative_side(vqps)) {
              if ( (p.is_segment() && l.has_on_positive_side(vqps))
                  ||(p.is_point() && l.has_on_negative_side(vqps)) ) {
                npts = 5;
              } else {
                npts = 4;
              }
            }
            else {
              if ( (p.is_segment() && l.has_on_negative_side(vpqr))
                  ||(p.is_point() && l.has_on_positive_side(vpqr)) ) {
                npts = 3;
              }
              else {
                npts = 2;
              }
            }
          }

        }//end of general segment case, seg != hor or ver
      }// end of point-seg case where pnt is not end point of segment
      points[npts-1] = vqps;
      Polychainsegment pcs(points, points+npts);

      CGAL_SDG_DEBUG( std::cout << "debug: construct bisector segment is "
          << pcs << " with npts = " << npts << std::endl; );

      return pcs;
    }//end of point segment case
    else { //  p and q both are segment
      unsigned int npts = 2;
      Point_2 points[2];
      points[0] = vpqr;
      points[npts-1] = vqps;
      Polychainsegment pcs(points, points+npts);

      CGAL_SDG_DEBUG( std::cout << "debug construct bisector segment is "
        << pcs << std::endl; );

      return pcs;
    } // end of segment segment case

  }
};


//-----------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_CONSTRUCTIONS_C2_H
