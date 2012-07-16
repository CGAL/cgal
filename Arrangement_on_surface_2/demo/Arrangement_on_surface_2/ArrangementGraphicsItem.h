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
    CGAL::Bbox_2 bb;
    bool bb_initialized;

    QPen verticesPen;
    QPen edgesPen;
    bool visible_edges;
    bool visible_vertices;

    QGraphicsScene* scene;

}; // class ArrangementGraphicsItemBase

template < class Arr_ >
class ArrangementGraphicsItem : public ArrangementGraphicsItemBase
{
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Arrangement::Edge_iterator Edge_iterator;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;

public:
    ArrangementGraphicsItem( Arrangement* t_ );
    void modelChanged( );

public:
    // QGraphicsItem overrides
    QRectF boundingRect( ) const;
    virtual void paint( QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget );

protected:
    void updateBoundingBox( );

    Arrangement* arr;
    ArrangementPainterOstream< Traits > painterostream;
    CGAL::Qt::Converter< Kernel > convert;
}; // class ArrangementGraphicsItem

template < class Arr_ >
ArrangementGraphicsItem< Arr_ >::
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

template < class Arr_ >
QRectF 
ArrangementGraphicsItem< Arr_ >::
boundingRect( ) const
{
    QRectF rect = this->convert( this->bb );
    return rect;
}

template < class Arr_ >
void 
ArrangementGraphicsItem< Arr_ >::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
    painter->setClipping( true );
    //painter->drawRect( this->boundingRect( ) );

#if 0
    QTransform transform = painter->worldTransform( );
    if ( transform.isScaling( ) )
    {
        std::cout << "Cool, we're scaling." << std::endl;
        std::cout << transform.m11( ) << std::endl;
    }
#endif
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
#if 0
        Point_2 p1 = it->source( )->point( );
        Point_2 p2 = it->target( )->point( );
        Segment_2 edge( p1, p2 );
        this->painterostream << edge;
#endif
        X_monotone_curve_2 curve = it->curve( );
        this->painterostream << curve;
    }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template < class Arr_ >
void 
ArrangementGraphicsItem< Arr_ >::updateBoundingBox( )
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

    for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( );
        ++it )
    {
        this->bb = this->bb + it->point( ).bbox( );
    }
}

template < class Arr_ >
void 
ArrangementGraphicsItem< Arr_ >::modelChanged( )
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

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
