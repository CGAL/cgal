// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann

#ifndef CGAL__TEST_D_C
#define CGAL__TEST_D_C

#include <CGAL/_test_cls_vector_d.C>
#include <CGAL/_test_fct_vector_d.C>
#include <CGAL/_test_cls_point_d.C>
#include <CGAL/_test_fct_point_vector_d.C>
#include <CGAL/_test_fct_point_d.C>
#include <CGAL/_test_cls_direction_d.C>
#include <CGAL/_test_cls_plane_d.C>
#include <CGAL/_test_cls_line_d.C>
#include <CGAL/_test_cls_segment_d.C>
#include <CGAL/_test_cls_ray_d.C>
#include <CGAL/_test_cls_triangle_d.C>
#include <CGAL/_test_cls_tetrahedron_d.C>
#include <CGAL/_test_cls_simplex_d.C>
#include <CGAL/_test_cls_aff_transformation_d.C>

template <class R>
bool
_test_d(const R& r)
{
 return
    _test_cls_vector_d(r)
 && _test_fct_vector_d(r)
 && _test_cls_point_d(r)
 && _test_fct_point_vector_d(r)
 && _test_cls_direction_d(r)
 && _test_cls_plane_d( r )
 && _test_cls_line_d( r )
 && _test_cls_ray_d( r )
 && _test_cls_segment_d( r )
 && _test_cls_triangle_d( r )
 && _test_cls_tetrahedron_d( r )
 && _test_fct_point_d(r)
 && _test_cls_simplex_d( r )
 && _test_cls_aff_transformation_d( r )
 ;
}

#endif // CGAL__TEST_D_C
