#ifndef CGAL_QT_WIDGET_MIN_ELLIPSE_2_H
#define CGAL_QT_WIDGET_MIN_ELLIPSE_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Min_ellipse_2.h>

namespace CGAL{

template< class Traits_ >
Qt_widget&
operator<<(Qt_widget &ws,
              const CGAL::Min_ellipse_2<Traits_>& min_ellipse)
{
    typedef CGAL::Min_ellipse_2<Traits_>::Point_iterator  Point_iterator;

    Point_iterator  first( min_ellipse.points_begin());
    Point_iterator  last ( min_ellipse.points_end());
    for ( ; first != last; ++first)
        ws << *first;
    return ws;
}

}//end namespace CGAL

#endif
