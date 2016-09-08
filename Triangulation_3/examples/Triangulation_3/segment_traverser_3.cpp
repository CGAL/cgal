#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <assert.h>
#include <iostream>

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >                DT;

typedef DT::Cell_handle                                         Cell_handle;

typedef CGAL::Triangulation_segment_cell_iterator_3< DT >       Traverser;

int main()
{
    std::vector< Point_3 > points;
    points.reserve( 6 );
    points.push_back( Point_3(-.9063, .4226, 1.0 ) );
    points.push_back( Point_3( .8192, .5736, 1.0 ) );
    points.push_back( Point_3( .0872,-.9992, 1.0 ) );
    points.push_back( Point_3(-.8192, .5736,-1.0 ) );
    points.push_back( Point_3( .9063, .4226,-1.0 ) );
    points.push_back( Point_3( .0872,-.9962,-1.0 ) );
    
    // Construct the Delaunay triangulation.
    DT dt( points.begin(), points.end() );
    assert( dt.is_valid() );

    // Construct a traverser.
    Traverser st( dt, Point_3(-3,0,0), Point_3(3,0,0) );

    // Count the number of finite cells traversed.
    unsigned int inf = 0, fin = 0;
    for( ; st != st.end(); ++st ) {
        if( dt.is_infinite(st) )
            ++inf;
        else
            ++fin;
    }

    std::cout << "While traversing from " << st.source() << " to " << st.target() << std::endl;
    std::cout << inf << " infinite and " << fin << " finite cells were visited." << std::endl;

    return 0;
}