#ifndef BOX_INTERSECTION_D_DEFINITIONS_H
#define BOX_INTERSECTION_D_DEFINITIONS_H

#include <vector>
#include <algorithm> // for pair
#include <iostream>
#include <iterator>
#include <cstdio>
#include <cmath>

#include <CGAL/Box_intersection_d.h>
#include <CGAL/Timer.h>
#include <CGAL/Random.h>

template< class NT, unsigned int DIM, bool CLOSED = true >
struct Definitions {
    typedef NT Number_type;
    typedef CGAL::Box_intersection_d::Box_d< Number_type, DIM >  Box;
    typedef CGAL::Box_intersection_d::Box_traits_d< Box > Box_adapter;
    typedef CGAL::Box_intersection_d::Box_predicate_traits_d<
                                             Box_adapter, CLOSED > Traits;
    typedef std::vector< Box >      Box_container;
    typedef std::pair< Box, Box >   Box_pair;
    typedef std::vector< Box_pair > Result_container;

    static const unsigned int dim = DIM;
};

#endif
