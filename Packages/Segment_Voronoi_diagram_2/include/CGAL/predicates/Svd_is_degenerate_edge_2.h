#ifndef CGAL_SVD_IS_DEGENERATE_EDGE_2_H
#define CGAL_SVD_IS_DEGENERATE_EDGE_2_H

//#include <CGAL/predicates/Svd_basic_predicates_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>


CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------



template<class R, class Method_tag>
class Svd_is_degen_edge_2
{
public:
  typedef typename R::Site_2      Site_2;

public:
  bool operator()(const Site_2& p, const Site_2& q,
		  const Site_2& r, const Site_2& s)
  {
    return false;
  }

};


//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_IS_DEGENERATE_EDGE_2_H
