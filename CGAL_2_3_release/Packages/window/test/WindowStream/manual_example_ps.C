
// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : include/CGAL/test/WindowStream/manual_example_ps.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#if !defined(CGAL_USE_LEDA)
int main() { return 0; }
#else
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Postscript_file_stream.h>

typedef CGAL::Point_2< CGAL::Cartesian<double> >     Point;
typedef CGAL::Segment_2< CGAL::Cartesian<double> >   Segment;

int main()
{
    Point p(0,1), q(2,2);
    Segment s(p,q);

    CGAL::Postscript_file_stream PS(100,100);
    PS.init(0,10,0);
    CGAL::cgalize( PS);
    PS.display();

    PS << CGAL::RED << s << CGAL::BLACK << p << q ;
    return 0;
}
#endif // USE_LEDA
