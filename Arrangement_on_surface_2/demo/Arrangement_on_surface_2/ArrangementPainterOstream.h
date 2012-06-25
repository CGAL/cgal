#ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#include <CGAL/Kernel/global_functions.h> // TODO: should be included in PainterOstream.h
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <QRectF>
#include <CGAL/Arr_segment_traits_2.h>

#include "Utils.h"

class QPainter;

namespace CGAL {
namespace Qt {

template < class ArrTraits >
class ArrangementPainterOstreamBase
{
public:
    typedef ArrTraits Traits;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::Ray_2 Ray_2;
    typedef typename Kernel::Line_2 Line_2;
    typedef typename Kernel::Triangle_2 Triangle_2;
    typedef typename Kernel::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename Kernel::Circle_2 Circle_2;

    ArrangementPainterOstreamBase( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        painterOstream( p, clippingRectangle ),
        qp( p ),
        convert( clippingRectangle )
    { }

    template < class T >
    ArrangementPainterOstreamBase& operator<<( const T& t )
    {
        this->painterOstream << t;
        return *this;
    }

protected:
    PainterOstream< Kernel > painterOstream;
    QPainter* qp;
    Converter< Kernel > convert;
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
public:
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

    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        Superclass( p, clippingRectangle )
    { }

    ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
    {
        std::cout << "ArrangementPainterOstream::paint curve" << std::endl;
        const Point_2& p1 = curve.source( );
        const Point_2& p2 = curve.target( );
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
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
public:
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

    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        Superclass( p, clippingRectangle )
    { }

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
