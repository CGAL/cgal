#ifndef CGAL_SVD_INFINITE_EDGE_INTERIOR_2_H
#define CGAL_SVD_INFINITE_EDGE_INTERIOR_2_H

#include <CGAL/predicates/Svd_basic_predicates_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>


CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------



template<class R, class Method_tag>
class Svd_infinite_edge_interior_2
{
public:
  typedef typename R::Site_2      Site_2;

public:
  bool operator()(const Site_2& q, const Site_2& s, const Site_2& r,
		  const Site_2& t, Sign sgn) const
  {
    if ( t.is_segment() ) {
#if PRED_PRINT
      std::cout << false << std::endl;
      return false;
#endif
    }

#if 0
    if ( q.is_segment() ) {
      // in this case r and s must be endpoints of q
      return ( sgn == NEGATIVE );
    }
#endif

    return ( sgn == NEGATIVE );
  }

};


//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_INFINITE_EDGE_INTERIOR_2_H
