// Copyright (c) 2009  GeometryFactory Sarl (France).
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
// $URL:$
// $Id:$
// 
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>

#ifndef CGAL_QT_GENERAL_POLYGON_SET_GRAPHICS_ITEM_H
#define CGAL_QT_GENERAL_POLYGON_SET_GRAPHICS_ITEM_H

#include <vector>
#include <list>

#include <boost/optional.hpp>
#include <boost/utility.hpp>

#include <CGAL/function_objects.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>


#include <QPainter>
#include <QBrush>
#include <QPen>

namespace CGAL {

namespace Qt {

  namespace internal
  {
    template<class General_polygon_set> 
    class Polygon_set_traits 
    {
      typedef typename General_polygon_set::Base base ;
      
    public:
      
      typedef typename base::Polygon_with_holes_2    Polygon_with_holes_2 ;
      typedef typename base::Polygon_2               Polygon_2 ;
      typedef typename General_polygon_set::Traits_2 Traits_2;
    } ;
    
    template<class K, class C, class D>
    class Polygon_set_traits< Polygon_set_2<K,C,D> >
    {
      typedef Polygon_set_2<K,C,D> PS ;
      
    public:
      
      typedef typename PS::Polygon_with_holes_2 Polygon_with_holes_2 ;
      typedef typename PS::Polygon_2            Polygon_2 ;
      typedef typename PS::Traits_2             Traits_2;
    } ;
    
    template<class TPoint,class SPoint>
    struct point_cast_helper
    {
      static TPoint convert( SPoint const& p )
      {
        return TPoint( to_double(p.x()), to_double(p.y()) ) ;
      }
      
    };
    
    template <class TPoint, class R, class A, class N, class B>
    struct point_cast_helper<TPoint, CGAL::_Bezier_point_2<R,A,N,B> >
    {
      static TPoint convert( CGAL::_Bezier_point_2<R,A,N,B> const& p )
      {
        std::pair<double,double> xy = p.approximate();
        return TPoint(xy.first, xy.second);
      }
    } ;
    
    template<class TPoint, class SPoint>
    TPoint point_cast ( SPoint const& aP ) { return point_cast_helper<TPoint,SPoint>::convert(aP); }
  }

template <class General_polygon_set_, class Compute_XM_curve_bbox_, class Draw_XM_curve_>
class GeneralPolygonSetGraphicsItem : public GraphicsItem
{
  typedef internal::Polygon_set_traits<General_polygon_set_> GPS_traits ;
    
  typedef Compute_XM_curve_bbox_                    Compute_XM_curve_bbox ;
  typedef Draw_XM_curve_                            Draw_XM_curve ;
  typedef General_polygon_set_                      General_polygon_set ;
  typedef typename GPS_traits::Polygon_with_holes_2 General_polygon_with_holes ;
  typedef typename GPS_traits::Polygon_2            General_polygon ;
  typedef typename GPS_traits::Traits_2             General_traits;
  typedef typename General_traits::Point_2          General_point;
  typedef std::vector<General_polygon_with_holes>   General_polygon_with_holes_vector ;
   
  typedef typename General_polygon_with_holes_vector::const_iterator General_pwh_const_iterator ; 
  typedef typename General_polygon_with_holes::Hole_const_iterator   General_hole_const_itertator ;
  typedef typename General_polygon::Curve_const_iterator             General_curve_const_iterator ;
  typedef typename General_polygon::X_monotone_curve_2               General_X_monotone_curve_2 ;

  typedef Converter< Simple_cartesian<double> > ToQtConverter;
  
public:

  GeneralPolygonSetGraphicsItem( General_polygon_set*         aSet
                               , Compute_XM_curve_bbox const& aBBoxComputer = Compute_XM_curve_bbox()
                               , Draw_XM_curve  const&        aDrawe        = Draw_XM_curve()
                               );

  void modelChanged();

public:

  bool isModelEmpty() const { return !mSet || mSet->is_empty() ; }
  
  QRectF boundingRect() const { return mBounding_rect ; }
  
  void paint(QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget);
  
  const QBrush& brush() const { return mBrush; }
  
  void setBrush(const QBrush& aBrush ) { mBrush = aBrush; }

  const QPen& pen() const{ return mPen; }

  void setPen(const QPen& aPen) { mPen = aPen; }

protected:

  struct ComputeBBox
  {
    ComputeBBox( Compute_XM_curve_bbox const& aComputer ) : mComputer(aComputer) {}
     
    template<class XM_curve>
    void operator() ( XM_curve const& aCurve, int ) 
    {
      Bbox_2 lB = mComputer(aCurve); 
      
      if ( mBBox )
           mBBox = *mBBox + lB;
      else mBBox =          lB;
    }
    
    Compute_XM_curve_bbox   mComputer ;
    boost::optional<Bbox_2> mBBox ;
  } ;
  
  struct DrawCurve
  {
    DrawCurve( QPainterPath* aPath, Draw_XM_curve const& aDrawer ) : mPath(aPath), mDrawer(aDrawer) {}
    
    template<class XM_curve>
    void operator() ( XM_curve const& aCurve, int aIdx ) 
    {
      mDrawer(aCurve,mPath,ToQtConverter(),aIdx);
    }
    
    QPainterPath*        mPath ;
    Draw_XM_curve const& mDrawer ;
  } ;
  
  template<class Visitor> void traverse_polygon           ( General_polygon const& aP, Visitor& aVisitor) ;
  template<class Visitor> void traverse_polygon_with_holes( General_polygon_with_holes const& aPWH, Visitor& aVisitor) ;
  template<class Visitor> void traverse_set               ( Visitor& aVisitor ) ;
  
  void updateBoundingBox();

protected:

  General_polygon_set*  mSet;
  Compute_XM_curve_bbox mBBoxComputer ;
  Draw_XM_curve         mDrawer ;
  QRectF                mBounding_rect;
  QBrush                mBrush;
  QPen                  mPen;
};


template <class G, class B, class D>
GeneralPolygonSetGraphicsItem<G,B,D>::GeneralPolygonSetGraphicsItem( General_polygon_set*         aSet
                                                                   , Compute_XM_curve_bbox const& aBBoxComputer
                                                                   , Draw_XM_curve  const&        aDrawer       
                                                                   )
  : mSet         (aSet)
  , mBBoxComputer(aBBoxComputer)
  , mDrawer      (aDrawer)
{
  updateBoundingBox();
}

template <class G, class B, class D>
void GeneralPolygonSetGraphicsItem<G,B,D>::paint( QPainter* aPainter, const QStyleOptionGraphicsItem* aOption, QWidget* aWidget )
{
  if ( ! isModelEmpty() )
  {
    QPainterPath lPath ;
    
    DrawCurve lVisitor(&lPath,mDrawer);
    
    traverse_set(lVisitor);
    
    aPainter->setPen  (mPen );
    aPainter->setBrush(mBrush);
    aPainter->drawPath(lPath);
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <class G, class B, class D>
void GeneralPolygonSetGraphicsItem<G,B,D>::updateBoundingBox()
{
  if ( ! isModelEmpty() )
  {
    prepareGeometryChange();
    
    ComputeBBox lVisitor(mBBoxComputer);
    
    traverse_set(lVisitor);  
    
    if ( lVisitor.mBBox ) 
    {
      ToQtConverter to_Qt ;
      mBounding_rect = to_Qt(*lVisitor.mBBox);
    }  
  }
}


template <class G, class B, class D>
void GeneralPolygonSetGraphicsItem<G,B,D>::modelChanged()
{
  updateBoundingBox();
  update();
}

template <class G, class B, class D>
template<class Visitor>
void GeneralPolygonSetGraphicsItem<G,B,D>::traverse_polygon( General_polygon const& aP, Visitor& aVisitor)
{
  int c = 0 ;
  for( General_curve_const_iterator cit = aP.curves_begin(); cit != aP.curves_end(); ++ cit )
      aVisitor(*cit,c++);
}

template <class G, class B, class D>
template<class Visitor>
void GeneralPolygonSetGraphicsItem<G,B,D>::traverse_polygon_with_holes( General_polygon_with_holes const& aPWH, Visitor& aVisitor)
{
  traverse_polygon(aPWH.outer_boundary(), aVisitor) ;
  
  for( General_hole_const_itertator hit = aPWH.holes_begin(); hit != aPWH.holes_end(); ++ hit )
    traverse_polygon(*hit, aVisitor);
}

template <class G, class B, class D>
template<class Visitor>
void GeneralPolygonSetGraphicsItem<G,B,D>::traverse_set( Visitor& aVisitor )
{
  General_polygon_with_holes_vector vec ;
  mSet->polygons_with_holes( std::back_inserter(vec) ) ;
  for( General_pwh_const_iterator lit = vec.begin(); lit != vec.end(); ++ lit )
    traverse_polygon_with_holes(*lit,aVisitor);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GENERAL_POLYGON_SET_GRAPHICS_ITEM_H
