#ifndef CGAL_CIRCLE_2_BITANGENT_2_INTERSECTION_H
#define CGAL_CIRCLE_2_BITANGENT_2_INTERSECTION_H

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CEP/Visibility_complex/Bitangent_2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class R_ , class C_ >
bool do_intersect( const Circle_2<R_>& c1, const Bitangent_2<C_>& c2 )
{
    cerr << "Not implemented Circle_2 - Bitangent_2 intersection !" << endl;
    return false;
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_CIRCLE_2_BITANGENT_2_INTERSECTION_H
