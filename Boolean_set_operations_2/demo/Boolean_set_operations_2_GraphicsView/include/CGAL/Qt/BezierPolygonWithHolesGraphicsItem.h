// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Polyline_simplification_2/demo/Polyline_simplification_2/include/CGAL/Qt/Polyline_simplification_2_graphics_item.h $
// $Id: Polyline_simplification_2_graphics_item.h 48710 2009-04-07 21:41:12Z fcacciola $
// 
//
// Author(s) : Fernando Cacciola <Fernando.Cacciola @geometryfactory.com>

#ifndef CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
#define CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H

#include <list>

#include <boost/shared_ptr.hpp>

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QPainter>
#include <QBrush>
#include <QPen>

namespace CGAL {

namespace Qt {

template <class Bezier_polygon_with_holes_, class Converter_, class Path_builder_>
class Bezier_to_path
{
  typedef Bezier_polygon_with_holes_ Bezier_polygon_with_holes ;
  typedef Converter_                 Converter ;
  typedef Path_builder_              Path_builder ;
  
  typedef typename Bezier_polygon_with_holes::Hole_const_iterator Hole_const_iterator ;
  typedef typename Bezier_polygon_with_holes::General_polygon_2   Bezier_polygon ;
  typedef typename Bezier_polygon::Curve_const_iterator           Curve_const_iterator ;
  typedef typename Bezier_polygon::X_monotone_curve_2             Bezier_X_monotone_curve ;
  typedef typename Bezier_X_monotone_curve::Point_2               Bezier_point ;
  typedef typename Bezier_X_monotone_curve::Curve_2               Bezier_curve ;
  
  typedef Simple_cartesian<double> Dbl_kernel;
  typedef Dbl_kernel::Line_2       Dbl_Line ;
  typedef Dbl_kernel::Point_2      Dbl_point;
  typedef std::deque<Dbl_point>    Dbl_control_points ;

public :
  
  Bezier_to_path( Converter& aConverter, Path_builder& aPath_builder ) 
    : convert(aConverter)
    , path_builder(aPath_builder) 
    , mParameterApproximationError(0.001)
    , mSubdivisionErrorThreshold(0.1)
  {}
   
  template<class OutputIterator> 
  void recursive_subdivision( Dbl_control_points const& aCtrlPts, double aMin, double aMid, double aMax, double aErrThreshold, OutputIterator aOut )
  {
    Dbl_control_points lLeft;
    Dbl_control_points lRight;
    
    Dbl_point lP = de_Casteljau_2(aCtrlPts.begin(), aCtrlPts.end(), aMid, std::back_inserter(lLeft), std::front_inserter(lRight) ) ;
    
    Dbl_point lFirst = aCtrlPts.front();
    Dbl_point lLast  = aCtrlPts.back ();
    
    double lErr = CGAL::square_distance(lP, Line(lFirst, lLast));    
    
    if ( lErr > aErrThreshold )
    {
      recursive_subdivision(lLeft,aMin,(aMin+aMid)/2,aMid,aOut);
      
      *aOut ++ = lP ;
      
      recursive_subdivision(lLeft,aMid,(aMid+aMax)/2,aMax,aOut);
    }
    else
    {
      *aOut ++ = lP ;
    }
    
  }
   
  void process( Bezier_X_monotone_curve const& aBXMC, bool first_point )
  {
    Bezier_cuve const& lBC = aBXMC.supporting_curve();
    
    double lST = get_approximate_endpoint_parameter(aBXMC.source(), lBC, aBXMC.xid(), mParameterApproximationError ) ;
    double lET = get_approximate_endpoint_parameter(aBXMC.target(), lBC, aBXMC.xid(), mParameterApproximationError ) ;
    
    Dbl_control_points lCtrlPoints ;
    for ( int k = 0 ; k < lBC.number_of_control_points(); ++ k )
    {
      Dbl_point p ( CGAL::to_double ( lBC.control_point(k).x() )
                  , CGAL::to_double ( lBC.control_point(k).y() )
                  ) ;
      
    }
    
    std::vector<Dbl_point> lSample ;
    
    recursive_subdivision(lCtrlPoints, lST, ( lST + lET ) / 2.0 , lET, mSubdivisionErrorThreshold, std::back_inserter(lSample) ) ;
    
    for( typename std::vector<Dbl_point>::const_iterator it = lSample.begin(); it != lSample.end() ; ++ it )
    {
       
      if ( first_point && it == lSample.begin() )
            path_builder.moveTo( convert(*it) ) ;
      else  path_builder.lineTo( convert(*it) ) ;
    }
    
  }
  
  void process( Bezier_polygon const& aBP)
  {
    bool first_point = true ;
    for( Curve_const_iterator cit = aBP.curves_begin(); cit != aBP.curves_end(); first_point = false, ++ cit )
      process(*cit,first_point);    
  }
  
  void process( Bezier_polygon_with_holes const& aBPWH )
  {
    process(aBPWH.outer_boundary());
    
    for( Hole_const_iterator hit = aBPWH.holes_begin(); hit != aBPWH.holes_end(); ++ hit )
      process(*hit);    
  }
  
private:
  
  double get_approximate_endpoint_parameter( Bezier_point const& p, Bezier_curve const& curve, unsigned int xid, double err )
  {
    typedef typename Bezier_point::Originator_iterator Originator_iterator ;
    
    // First try to use the approximate representation of the endpoints.
    Originator_iterator org = p.get_originator (curve, xid);
  
    double t_min = CGAL::to_double(org->point_bound().t_min) ;
    double t_max = CGAL::to_double(org->point_bound().t_max) ;
    
    bool can_refine = !p.is_exact();
    
    do
    {
      if ( t_max - t_min <= err )
        break ;
      
      if ( can_refine )
      {
        can_refine = p.refine();
        
        t_min = CGAL::to_double(org->point_bound().t_min) ;
        t_max = CGAL::to_double(org->point_bound().t_max) ;
      }  
    }
    while ( can_refine ) ;
    
    return ( t_min + t_max) / 2.0 ;
  }
  
private:
  
  Converter&    convert ;
  Path_builder& path_builder ;
  
  double        mParameterApproximationError ;
  double        mSubdivisionErrorThreshold;
} ;

template <typename Bezier_polygon_wiht_holes_>
class BezierPolygonWithHolesGraphicsItem : public GraphicsItem
{
  typedef Bezier_polygon_wiht_holes_ Bezier_polygon_wiht_holes ;
  
  typedef std::list<Bezier_polygon_wiht_holes> Bezier_polygon_wiht_holes_list ;
  
  typedef typename Bezier_polygon_wiht_holes_list::iterator list_iterator ;
  
  typedef typename Bezier_polygon_wiht_holes::General_polygon_2 Bezier_polygon ;
  
  typedef typename Bezier_polygon::Traits_2 Traits;
  
  typedef typename Traits::Point_2 Point;
  
public:

  BezierPolygonWithHolesGraphicsItem( Bezier_polygon_wiht_holes_list* aList );

  void modelChanged();

public:

  QRectF boundingRect() const { return mBounding_rect ; }
  
  void paint(QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget);
  
  const QBrush& brush() const { return mBrush; }
  
  void setBrush(const QBrush& aBrush ) { mBrush = aBrush; }

  const QPen& pen() const{ return mPen; }

  void setPen(const QPen& aPen) { mPen = aPen; }

protected:

  void updateBoundingBox();

  bool isModelEmpty() const { return !mList || mList->size() == 0 || mList->front().outer_boundary().size() == 0 ; }

protected:

  Bezier_polygon_wiht_holes_list* mList;
  QRectF                          mBounding_rect;
  QBrush                          mBrush;
  QPen                            mPen;
};


template <typename Bezier_polygon_wiht_holes>
BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::BezierPolygonWithHolesGraphicsItem(Bezier_polygon_wiht_holes_list* aList)
  : mList(aList)
{
  if( ! isModelEmpty() )
       updateBoundingBox();
  else this->hide();
}

template <typename Bezier_polygon_wiht_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::paint( QPainter*                       aPainter
                                                                         , const QStyleOptionGraphicsItem* aOption
                                                                         , QWidget*                        aWidget
                                                                         )
{
  if ( ! isModelEmpty() )
  {
    Converter<Traits> lConverter;
    
    QPainterPath lPath;
    
    Bezier_to_path<Bezier_polygon_wiht_holes, Converter<Traits>, QPainterPath > to_path(lConverter,lPath);
    
    for ( list_iterator it = mList->begin() ; it != mList->end() ; ++ it )
      to_path.process(*it);
    
    aPainter->setPen  (mPen );
    aPainter->setBrush(mBrush);
    aPainter->drawPath(lPath);
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename Bezier_polygon_wiht_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::updateBoundingBox()
{
  if ( ! isModelEmpty() )
  {
    prepareGeometryChange();
    Converter<Traits> convert;
//    mBounding_rect = convert(mList->front().outer_boundary().bbox());
  }
}


template <typename Bezier_polygon_wiht_holes>
void BezierPolygonWithHolesGraphicsItem<Bezier_polygon_wiht_holes>::modelChanged()
{
  if ( !isModelEmpty() && this->isVisible() )
       this->show();
  else this->hide();
  
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_BEZIER_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
