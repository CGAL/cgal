#ifndef CGAL_CIRCLE_2_CIRCLE_2_INTERSECTION_H
#define CGAL_CIRCLE_2_CIRCLE_2_INTERSECTION_H

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CEP/Visibility_complex/predicates/Circle_2_Circle_2_intersection_ftC2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class Rep_ >
bool do_intersect( const Circle_2<Rep_>& c1, const Circle_2<Rep_>& c2 )
{
    return circle_2_do_intersectC2(c1.center().x(),c1.center().y(),
				   c1.squared_radius(),
				   c2.center().x(),c2.center().y(),
				   c2.squared_radius());
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
