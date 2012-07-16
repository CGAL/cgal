#ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#include <QRectF>

#include <CGAL/Kernel/global_functions.h> // TODO: should be included in PainterOstream.h
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>

#include "Utils.h"

class QPainter;

namespace CGAL {
namespace Qt {

template < class ArrTraits >
class ArrangementPainterOstreamBase
{
public: // typedefs
    typedef ArrTraits Traits;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::Ray_2 Ray_2;
    typedef typename Kernel::Line_2 Line_2;
    typedef typename Kernel::Triangle_2 Triangle_2;
    typedef typename Kernel::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename Kernel::Circle_2 Circle_2;

public: // constructors
    ArrangementPainterOstreamBase( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        painterOstream( p, clippingRectangle ),
        qp( p ),
        convert( clippingRectangle ),
        scene( NULL ),
        clippingRect( QRectF( ) ), // null rectangle
        scale( 1.0 )
    {
        if ( p != 0 )
        {
            this->scale = p->worldTransform( ).m11( );
        }
    }

public: // methods
    template < class T >
    ArrangementPainterOstreamBase& operator<<( const T& t )
    {
        this->painterOstream << t;
        return *this;
    }

    void setScene( QGraphicsScene* scene_ )
    {
        this->scene = scene_;

        // set the clipping rectangle
        if ( scene_ == NULL )
        {
            return;
        }
        this->clippingRect = this->getViewportRect( );
    }

protected: // methods
    QRectF getViewportRect( ) const
    {
        // assumes scene is not null and attached to exactly one view
        QGraphicsView* view = this->scene->views( ).first( );
        QPointF p1 = view->mapToScene( 0, 0 );
        QPointF p2 = view->mapToScene( view->width( ), view->height( ) );
        QRectF clipRect = QRectF( p1, p2 );

        return clipRect;
    }

protected: // fields
    PainterOstream< Kernel > painterOstream;
    QPainter* qp;
    Converter< Kernel > convert;
    QGraphicsScene* scene;
    QRectF clippingRect;
    double scale;

}; // class ArrangementPainterOstreamBase

template < class ArrTraits >
class ArrangementPainterOstream:
    public ArrangementPainterOstreamBase< ArrTraits >
{
public:
    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        ArrangementPainterOstreamBase< ArrTraits >( p, clippingRectangle )
    { }
};

template < class Kernel_ >
class ArrangementPainterOstream< CGAL::Arr_segment_traits_2< Kernel_ > >:
    public ArrangementPainterOstreamBase< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public: // typedefs
    typedef Kernel_ Kernel;
    typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
    typedef ArrangementPainterOstreamBase< Traits > Superclass;
    typedef typename Superclass::Point_2 Point_2;
    typedef typename Superclass::Segment_2 Segment_2;
    typedef typename Superclass::Ray_2 Ray_2;
    typedef typename Superclass::Line_2 Line_2;
    typedef typename Superclass::Triangle_2 Triangle_2;
    typedef typename Superclass::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename Superclass::Circle_2 Circle_2;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public: // constructors
    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        Superclass( p, clippingRectangle )
    { }

public: // methods
    ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
    {
        const Point_2& p1 = curve.source( );
        const Point_2& p2 = curve.target( );
        Segment_2 seg( p1, p2 );

        // skip segments outside our view
        QRectF seg_bb = this->convert( seg.bbox( ) );
        if ( this->clippingRect.isValid( ) &&
            ! this->clippingRect.intersects( seg_bb ) )
        {
            return *this;
        }

        this->painterOstream << seg;
        return *this;
    }

    ArrangementPainterOstream& operator<<( const Point_2& p )
    {
        QPointF qpt = this->convert( p );
        // clip the point if possible
        if ( this->clippingRect.isValid( ) &&
            ! this->clippingRect.contains( qpt ) )
        {
            return *this;
        }

        QPen savePen = this->qp->pen( );
        this->qp->setBrush( QBrush( savePen.color( ) ) );
        double radius = savePen.width( ) / 2.0;
        radius /= this->scale;

        this->qp->drawEllipse( qpt, radius, radius );

        this->qp->setBrush( QBrush( ) );
        this->qp->setPen( savePen );
        return *this;
    }

    template < class T >
    ArrangementPainterOstream& operator<<( const T& p )
    {
        (*(static_cast< Superclass* >(this)) << p);
        return *this;
    }
};

template < class SegmentTraits >
class ArrangementPainterOstream< CGAL::Arr_polyline_traits_2< SegmentTraits > > :
    public ArrangementPainterOstreamBase< CGAL::Arr_polyline_traits_2< SegmentTraits > >
{
public: // typedefs
    typedef ArrangementPainterOstreamBase< CGAL::Arr_polyline_traits_2< SegmentTraits > > Superclass;
    typedef typename Superclass::Traits Traits;
    typedef typename Superclass::Kernel Kernel;
    typedef typename Superclass::Point_2 Point_2;
    typedef typename Superclass::Segment_2 Segment_2;
    typedef typename Superclass::Ray_2 Ray_2;
    typedef typename Superclass::Line_2 Line_2;
    typedef typename Superclass::Triangle_2 Triangle_2;
    typedef typename Superclass::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename Superclass::Circle_2 Circle_2;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public: // constructors
    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        Superclass( p, clippingRectangle )
    { }

public: // methods
    ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
    {
        for ( int i = 0; i < curve.size( ); ++i )
        {
            Segment_2 segment = curve[ i ];
            this->painterOstream << segment;
        }
        // TODO: implement polyline painting
#if 0
        const Point_2& p1 = curve.source( );
        const Point_2& p2 = curve.target( );
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
#endif
        return *this;
    }

    // cloned from segtraits painter
    ArrangementPainterOstream& operator<<( const Point_2& p )
    {
        QPointF qpt = this->convert( p );
        QPen savePen = this->qp->pen( );
        this->qp->setBrush( QBrush( savePen.color( ) ) );
        double radius = savePen.width( ) / 2.0;
        radius /= this->scale;

        this->qp->drawEllipse( qpt, radius, radius );

        this->qp->setBrush( QBrush( ) );
        this->qp->setPen( savePen );
        return *this;
    }

    template < class T >
    ArrangementPainterOstream& operator<<( const T& p )
    {
        (*(static_cast< Superclass* >(this)) << p);
        return *this;
    }
};

template < class RatKernel, class AlgKernel, class NtTraits >
class ArrangementPainterOstream< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >:
    public ArrangementPainterOstreamBase< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >
{
public: // typedefs
    typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
    typedef ArrangementPainterOstreamBase< Traits > Superclass;
    typedef typename Superclass::Point_2 Point_2;
    typedef typename Superclass::Segment_2 Segment_2;
    typedef typename Superclass::Ray_2 Ray_2;
    typedef typename Superclass::Line_2 Line_2;
    typedef typename Superclass::Triangle_2 Triangle_2;
    typedef typename Superclass::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename Superclass::Circle_2 Circle_2;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public: // constructors
    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        Superclass( p, clippingRectangle )
    { }

public: // methods
    ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
    {
        // TODO: clip the curve to the visible viewing area
        
        int n;
        if ( this->scene == NULL )
            n = 100; // TODO: get an adaptive approximation
        else
        {
            QGraphicsView* view = this->scene->views( ).first( );
            CGAL::Bbox_2 bb = curve.bbox( );
            int xmin, xmax;
            xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
            xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
            n = xmax - xmin;
        }
        if ( n == 0 )
        {
            return *this;
        }

        std::pair< double, double >* app_pts = new std::pair< double, double >[ n + 1 ];
        std::pair< double, double >* end_pts = curve.polyline_approximation( n, app_pts );
        std::pair< double, double >* p_curr = app_pts;
        std::pair< double, double >* p_next = p_curr + 1;
        int count = 0;
        do
        {
            QPointF p1( p_curr->first, p_curr->second );
            QPointF p2( p_next->first, p_next->second );
#if 0
            Segment_2 seg( p1, p2 );
            this->painterOstream << seg;
#endif
            this->qp->drawLine( p1, p2 );
            p_curr++;
            p_next++;
            ++count;
        }
        while ( p_next != end_pts );
        //std::cout << count << " approximation points" << std::endl;

        return *this;
    }

    // cloned from segtraits painter
    ArrangementPainterOstream& operator<<( const Point_2& p )
    {
        QPointF qpt = this->convert( p );
        QPen savePen = this->qp->pen( );
        this->qp->setBrush( QBrush( savePen.color( ) ) );
        double radius = savePen.width( ) / 2.0;
        radius /= this->scale;

        this->qp->drawEllipse( qpt, radius, radius );

        this->qp->setBrush( QBrush( ) );
        this->qp->setPen( savePen );
        return *this;
    }

    template < class T >
    ArrangementPainterOstream& operator<<( const T& p )
    {
        (*(static_cast< Superclass* >(this)) << p);
        return *this;
    }
};

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
