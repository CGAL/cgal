#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Viewer_stream.h>
#include "text.h"


typedef CGAL::Cartesian<double> rep_t;
typedef CGAL::Point_3<rep_t> point_t;
typedef CGAL::Tetrahedron_3<rep_t> tetra;


typedef CGAL::Triangulation_geom_traits_3<rep_t>  traits_3;
typedef CGAL::Triangulation_vertex_base_3<traits_3>     Vb ;
typedef CGAL::Triangulation_cell_base_3<traits_3>       Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> TDS3 ;
typedef CGAL::Triangulation_3< traits_3 , TDS3> Triangulation_3;
typedef CGAL::Delaunay_triangulation_3<traits_3,TDS3> Delaunay_3;


int main(int argc, char *argv[]) 
{
  CGAL::Viewer_3 W(500);


  Delaunay_3 tr;
  
  tr.insert(point_t(100,100,100));
  tr.insert(point_t(100,100,100));
  tr.insert(point_t(400,100,300));
  tr.insert(point_t(-100,100,-100));
  tr.insert(point_t(100,-300,0));
  tr.insert(point_t(500,-100,-200));
  tr.insert(point_t(-140,400,-200));
  tr.insert(point_t(100,-350,0));
  tr.insert(point_t(500,-300,-250));

  CGAL::Drawable_triangulation_3<Delaunay_3> dtr(tr,CGAL::PURPLE);
  W.add_drawable(&dtr,1);
 W.main_loop();

}

	  
