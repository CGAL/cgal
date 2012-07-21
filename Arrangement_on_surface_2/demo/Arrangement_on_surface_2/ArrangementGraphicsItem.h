// Copyright (c) 2008  GeometryFactory Sarl (France).
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

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
//#include <CGAL/apply_to_range.h>
//#include <CGAL/Kernel/global_functions.h> // TODO: should be included in PainterOstream.h
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
//#include <CGAL/Arr_segment_traits_2.h>
//#include <CGAL/Arr_polyline_traits_2.h>

#include <QGraphicsScene>
#include <QPainter>
//#include <QStyleOption>

#include "ArrangementPainterOstream.h"
#include <iostream>

class QGraphicsScene;

namespace CGAL {
namespace Qt {

class ArrangementGraphicsItemBase : public GraphicsItem
{
public:
    ArrangementGraphicsItemBase( );

    const QPen& getVerticesPen( ) const;
    const QPen& getEdgesPen( ) const;
    void setVerticesPen( const QPen& pen );
    void setEdgesPen( const QPen& pen );
    bool visibleVertices( ) const;
    void setVisibleVertices( const bool b );
    bool visibleEdges( ) const;
    void setVisibleEdges( const bool b );
    void setScene( QGraphicsScene* scene_ );

protected:
    QRectF getViewportRect( ) const;

    CGAL::Bbox_2 bb;
    bool bb_initialized;

    QPen verticesPen;
    QPen edgesPen;
    bool visible_edges;
    bool visible_vertices;

    QGraphicsScene* scene;

}; // class ArrangementGraphicsItemBase

template < class Arr_, class ArrTraits = typename Arr_::Geometry_traits_2 >
class ArrangementGraphicsItem : public ArrangementGraphicsItemBase
{
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Arrangement::Curve_iterator Curve_iterator;
    typedef typename Arrangement::Edge_iterator Edge_iterator;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;

public:
    ArrangementGraphicsItem( Arrangement* t_ );

public:
    void modelChanged( );
    QRectF boundingRect( ) const;
    virtual void paint( QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget );

protected:
    void cacheCurveBoundingRects( );
    void updateBoundingBox( );

    Arrangement* arr;
    ArrangementPainterOstream< Traits > painterostream;
    CGAL::Qt::Converter< Kernel > convert;
    std::map< Curve_iterator, CGAL::Bbox_2 > curveBboxMap;
}; // class ArrangementGraphicsItem

template < class Arr_, class ArrTraits >
ArrangementGraphicsItem< Arr_, ArrTraits >::
ArrangementGraphicsItem( Arrangement* arr_ ):
    arr( arr_ ),
    painterostream( 0 )
{
    if ( this->arr->number_of_vertices( ) == 0 ) {
        this->hide( );
    }
    this->updateBoundingBox( );
    this->setZValue( 3 );
}

template < class Arr_, class ArrTraits >
QRectF 
ArrangementGraphicsItem< Arr_, ArrTraits >::
boundingRect( ) const
{
    QRectF rect = this->convert( this->bb );
    return rect;
}

template < class Arr_, class ArrTraits >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{

    painter->setPen( this->verticesPen );
    this->painterostream = ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
    this->painterostream.setScene( this->scene );

    for ( Vertex_iterator it = this->arr->vertices_begin( ); it != this->arr->vertices_end( ); ++it )
    {
        Point_2 pt = it->point( );
        this->painterostream << pt;
    }
    painter->setPen( this->edgesPen );
    for ( Edge_iterator it = this->arr->edges_begin( ); it != this->arr->edges_end( ); ++it )
    {
        X_monotone_curve_2 curve = it->curve( );
        this->painterostream << curve;
    }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template < class Arr_, class ArrTraits >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::updateBoundingBox( )
{
    this->prepareGeometryChange( );
    if ( this->arr->number_of_vertices( ) == 0 )
    {
        this->bb = Bbox_2( 0, 0, 0, 0 );
        this->bb_initialized = false;
        return;
    }
    else
    {
        this->bb = this->arr->vertices_begin( )->point( ).bbox( );
        this->bb_initialized = true;
    }

    for ( Curve_iterator it = this->arr->curves_begin( );
        it != this->arr->curves_end( );
        ++it )
    {
        if ( this->curveBboxMap.count( it ) == 0 )
        {
            this->curveBboxMap[ it ] = it->bbox( );
        }
        this->bb = this->bb + this->curveBboxMap[ it ];
    }
}

template < class Arr_, class ArrTraits >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::modelChanged( )
{
    if ( this->arr->is_empty( ) )
    {
        this->hide( );
    }
    else
    {
        this->show( );
    }
    this->updateBoundingBox( );
    this->update( );
}

/**
 * Why is this not being used?
Specialized methods:
    updateBoundingBox
*/
template < class Arr_, class Kernel_ >
class ArrangementGraphicsItem< Arr_, CGAL::Arr_linear_traits_2< Kernel_ > >  : public ArrangementGraphicsItemBase
{
    typedef Arr_ Arrangement;
    typedef ArrangementGraphicsItemBase Superclass;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Arrangement::Curve_iterator Curve_iterator;
    typedef typename Arrangement::Edge_iterator Edge_iterator;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;

public:
    ArrangementGraphicsItem( Arrangement* arr_ ):
        arr( arr_ ),
        painterostream( 0 )
    {
        if ( this->arr->number_of_vertices( ) == 0 ) {
            this->hide( );
        }
        this->updateBoundingBox( );
        this->setZValue( 3 );
    }

    void modelChanged( )
    {
        if ( this->arr->is_empty( ) )
        {
            this->hide( );
        }
        else
        {
            this->show( );
        }
        this->updateBoundingBox( );
        this->update( );
    }


public: // methods
    // @override QGraphicsItem::boundingRect
    QRectF boundingRect( ) const
    {
        QRectF rect = this->convert( this->bb );
        return rect;
    }

    // @override QGraphicsItem::paint
    virtual void paint( QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget )
    {
        this->updateBoundingBox( );
        painter->setPen( this->verticesPen );
        this->painterostream = ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
        this->painterostream.setScene( this->scene );

        for ( Vertex_iterator it = this->arr->vertices_begin( ); it != this->arr->vertices_end( ); ++it )
        {
            Point_2 pt = it->point( );
            this->painterostream << pt;
        }
        painter->setPen( this->edgesPen );
        for ( Edge_iterator it = this->arr->edges_begin( ); it != this->arr->edges_end( ); ++it )
        {
            X_monotone_curve_2 curve = it->curve( );
            this->painterostream << curve;
        }
    }

protected: // methods
    void updateBoundingBox( )
    {
        this->prepareGeometryChange( );
        QRectF clipRect = this->getViewportRect( );
        this->convert = Converter<Kernel>( clipRect );

        if ( ! clipRect.isValid( ) /*|| this->arr->number_of_vertices( ) == 0*/ )
        {
            this->bb = Bbox_2( 0, 0, 0, 0 );
            this->bb_initialized = false;
            return;
        }
        else
        {
            this->bb = this->convert( clipRect ).bbox( );
            this->bb_initialized = true;
        }

        for ( Curve_iterator it = this->arr->curves_begin( );
            it != this->arr->curves_end( );
            ++it )
        {
            if ( it->is_segment( ) )
            {
                this->bb = this->bb + it->segment( ).bbox( );
            }
            else if ( it->is_ray( ) )
            {
                QLineF qclippedRay = this->convert( it->ray( ) );
                Segment_2 clippedRay = this->convert( qclippedRay );
                this->bb = this->bb + clippedRay.bbox( );
            }
            else // ( it->is_line( ) )
            {
                QLineF qclippedLine = this->convert( it->line( ) );
                Segment_2 clippedLine = this->convert( qclippedLine );
                this->bb = this->bb + clippedLine.bbox( );
            }
        }
    }

protected: // fields
    Arrangement* arr;
    ArrangementPainterOstream< Traits > painterostream;
    CGAL::Qt::Converter< Kernel > convert;
}; // class ArrangementGraphicsItem

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
