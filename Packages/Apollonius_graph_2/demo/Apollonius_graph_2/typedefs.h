#ifndef _MK___TYPEDEFS_H
#define _MK___TYPEDEFS_H

#include <CGAL/basic.h>

//-------- choosing number type --- start ---------

#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float  exact_type;
//#include <CGAL/leda_real.h>
//typedef leda_real  exact_type;

#include <CGAL/Filtered_exact.h>
typedef double inexact_type;
typedef CGAL::Filtered_exact< inexact_type, exact_type >  number_type;

//-------- choosing number type --- end ---------

//#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian< number_type >        Rep;


#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>


typedef
CGAL::Apollonius_graph_traits_2<Rep,CGAL::Ring_tag>  Gt;

typedef Gt::Point_2                           Point_2;
typedef Rep::Circle_2                         Circle_2;
typedef Gt::Site_2                            Apollonius_site_2;
typedef Gt::Site_2::Weight                    Weight;

//typedef CGAL::Apollonius_graph_2<Gt> AG_2;
typedef CGAL::Apollonius_graph_hierarchy_2<Gt> AG_2;


typedef AG_2::Vertex_handle                 Vertex_handle;



#endif  // _MK___TYPEDEFS_H
