#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Viewer_stream.h>

typedef CGAL::Cartesian<double> rep_t;
typedef CGAL::Point_3<rep_t> point_t;
typedef CGAL::Tetrahedron_3<rep_t> tetra;

typedef CGAL::Triangle_3<rep_t> triangle;
typedef CGAL::Line_3<rep_t> line;
typedef CGAL::Segment_3<rep_t> segment;


typedef CGAL::Triangulation_geom_traits_3<rep_t>  traits_3;
typedef CGAL::Triangulation_vertex_base_3<traits_3>     Vb ;
typedef CGAL::Triangulation_cell_base_3<traits_3>       Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> TDS3 ;
typedef CGAL::Triangulation_3< traits_3 , TDS3> Triangulation_3;
typedef CGAL::Delaunay_triangulation_3<traits_3,TDS3> Delaunay_3;

CGAL::Color change(CGAL::Color &c)
{
  if (c==CGAL::RED) return CGAL::YELLOW;
  if (c==CGAL::YELLOW) return CGAL::ORANGE;
  if (c==CGAL::ORANGE) return CGAL::RED;
  return CGAL::WHITE;
}



tetra schrink(tetra tet)
{
  point_t p1,p2,p3,p4;
  double sx,sy,sz;
  double bx= (tet[0].x()+ tet[1].x()+ tet[2].x()+tet[3].x())/4;
  double by= (tet[0].y()+ tet[1].y()+ tet[2].y()+tet[3].y())/4;
  double bz= (tet[0].z()+ tet[1].z()+ tet[2].z()+tet[3].z())/4;
  sx= (bx - tet[0].x())/3;
  sy= (by - tet[0].y())/3;
  sz= (bz - tet[0].z())/3;
  p1 = point_t(tet[0].x()+sx,tet[0].y()+sy,tet[0].z()+sz);  
  sx= (bx - tet[1].x())/3;
  sy= (by - tet[1].y())/3;
  sz= (bz - tet[1].z())/3;
  p2 = point_t(tet[1].x()+sx,tet[1].y()+sy,tet[1].z()+sz);
  sx= (bx - tet[2].x())/3;
  sy= (by - tet[2].y())/3;
  sz= (bz - tet[2].z())/3;
  p3 = point_t(tet[2].x()+sx,tet[2].y()+sy,tet[2].z()+sz);
  sx= (bx - tet[3].x())/3;
  sy= (by - tet[3].y())/3;
  sz= (bz - tet[3].z())/3;
  p4 = point_t(tet[3].x()+sx,tet[3].y()+sy,tet[3].z()+sz);
  return tetra(p1,p2,p3,p4);
}

int main(int argc, char *argv[]) 
{

  CGAL::Viewer_3 W(-200, 600, -400, 400, -250, 300);

  Delaunay_3 tr;
  
  tr.insert(point_t(100,100,100));
  //  tr.insert(point_t(100,100,100));
  tr.insert(point_t(400,100,300));
  tr.insert(point_t(-100,100,-100));
  tr.insert(point_t(100,-300,0));
  tr.insert(point_t(500,-100,-200));

  tr.insert(point_t(-140,400,-200));
  tr.insert(point_t(100,-350,0));
  tr.insert(point_t(500,-300,-250));

  CGAL::Color c=CGAL::ORANGE ;
  Delaunay_3::Cell_iterator cit;
  int i=1;
  tetra t;
  for (cit = tr.finite_cells_begin(); cit != tr.cells_end(); cit++) {
    t= schrink(tr.tetrahedron(cit->handle()));
    CGAL::Drawable_tetrahedron_3<Delaunay_3::Tetrahedron>* tet = new CGAL::Drawable_tetrahedron_3<Delaunay_3::Tetrahedron>(t,c,CGAL::FILL);
    W.add_drawable(tet,i);
    i++;
    c=change(c);
  }
  

 W.main_loop();
 return 0;


}

	  
