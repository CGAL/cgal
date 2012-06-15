#ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#include <CGAL/Kernel/global_functions.h> // TODO: should be included in PainterOstream.h
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <QRectF>
#include <CGAL/Arr_segment_traits_2.h>

class QPainter;

namespace CGAL {
namespace Qt {

template < class ArrTraits >
class ArrangementPainterOstreamBase
{
public:
    typedef ArrTraits Traits;
    typedef typename Traits::Kernel Kernel;
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

    ArrangementPainterOstreamBase& operator<<( const Point_2& p )
    {
        this->painterOstream << p;
        return *this;
    }

    ArrangementPainterOstreamBase& operator<<( const Segment_2& s )
    {
        this->painterOstream << s;
        return *this;
    }


    ArrangementPainterOstreamBase& operator<<( const Ray_2& r )
    {
        this->painterOstream << r;
        return *this;
    }


    ArrangementPainterOstreamBase& operator<<( const Line_2& l )
    {
        this->painterOstream << l;
        return *this;
    }


    ArrangementPainterOstreamBase& operator<<( const Triangle_2& t )
    {
        this->painterOstream << t;
        return *this;
    }

    ArrangementPainterOstreamBase& operator<<( const Iso_rectangle_2& r )
    {
        this->painterOstream << r;
        return *this;
    }

    ArrangementPainterOstreamBase& operator<<( const Bbox_2& bb )
    {
        this->painterOstream << bb;
        return *this;
    }

    ArrangementPainterOstreamBase& operator<<( const Circle_2& c )
    {
        this->painterOstream << c;
        return *this;
    }

protected:
    PainterOstream< Kernel > painterOstream;
    QPainter* qp;
    Converter< Kernel > convert;
};

template < class ArrTraits >
class ArrangementPainterOstream : public ArrangementPainterOstreamBase< ArrTraits >
{
public:
    ArrangementPainterOstream( QPainter* p, QRectF clippingRectangle = QRectF( ) ):
        ArrangementPainterOstreamBase< ArrTraits >( p, clippingRectangle )
    { }
};

template < class Kernel_ >
class ArrangementPainterOstream< CGAL::Arr_segment_traits_2< Kernel_ > > :
    public ArrangementPainterOstreamBase< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
    typedef ArrangementPainterOstreamBase< CGAL::Arr_segment_traits_2< Kernel_ > > Superclass;
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
        const Point_2& p1 = curve.source( );
        const Point_2& p2 = curve.target( );
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
        return *this;
    }
};


} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
