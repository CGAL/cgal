//CGAL
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h> 
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h> 
#include <CGAL/Alpha_shape_face_base_2.h> 
#include <CGAL/Weighted_alpha_shape_euclidean_traits_2.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Kernel/function_objects.h>


typedef double                      Coord_type;
typedef CGAL::Cartesian<Coord_type> Rep;
typedef Rep::Point_2                Point_2;
typedef Rep::Segment_2              Segment;
typedef Rep::Line_2                 Line;
typedef Rep::Triangle_2             Triangle; 

typedef Rep::Less_xy_2              Point_compare;
typedef CGAL::Triangulation_2<Rep>  Triangulation;
typedef std::list<Point_2>          CGALPointlist;

//Weighted alpha_shape
typedef CGAL::Weighted_alpha_shape_euclidean_traits_2<Rep> Gt_w;
typedef CGAL::Alpha_shape_vertex_base_2<Gt_w>              Av;
typedef CGAL::Regular_triangulation_vertex_base_2<Gt_w, Av> Av_w;
typedef CGAL::Regular_triangulation_face_base_2<Gt_w>      Rf_w;
typedef CGAL::Alpha_shape_face_base_2<Gt_w,Rf_w>           Af_w;
typedef CGAL::Triangulation_default_data_structure_2<Gt_w,Av_w,Af_w> 
                                                           Tds_w;
typedef CGAL::Regular_triangulation_2<Gt_w,Tds_w>          Rt_w;
typedef CGAL::Alpha_shape_2<Rt_w>                          Alpha_shape_w;
typedef CGAL::Weighted_point<Point_2, double>              Wpoint;

//Delaunay triangulation
typedef Rep                                                Gt;
typedef CGAL::Alpha_shape_vertex_base_2<Gt>                Vb;
typedef CGAL::Triangulation_face_base_2<Gt>                Df;
typedef CGAL::Alpha_shape_face_base_2<Gt, Df>              Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb>
                                                           Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>             Delaunay;

//Alpha shape and the types
typedef CGAL::Alpha_shape_2<Delaunay>                      Alpha_shape;
typedef Alpha_shape::Face                                  Face;
typedef Alpha_shape::Vertex                                Vertex;
typedef Alpha_shape::Edge                                  Edge;
typedef Alpha_shape::Face_handle                           Face_handle;
typedef Alpha_shape::Vertex_handle                         Vertex_handle;
typedef Alpha_shape::Face_circulator                       Face_circulator;
typedef Alpha_shape::Vertex_circulator                     Vertex_circulator;
typedef Alpha_shape::Locate_type                           Locate_type;
typedef Alpha_shape::Face_iterator                         Face_iterator;
typedef Alpha_shape::Vertex_iterator                       Vertex_iterator;
typedef Alpha_shape::Edge_iterator                         Edge_iterator;
typedef Alpha_shape::Edge_circulator                       Edge_circulator;
typedef Alpha_shape::Alpha_iterator                        Alpha_iterator;


