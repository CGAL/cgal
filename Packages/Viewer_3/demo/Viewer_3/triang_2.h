
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>




typedef CGAL::Cartesian<double> rep_t;


typedef CGAL::Triangulation_euclidean_traits_2< rep_t > Ttraits;
typedef CGAL::Triangulation_vertex_base_2<Ttraits>      Vertex_base ;
typedef CGAL::Triangulation_face_base_2<Ttraits>        Face_base ;
typedef CGAL::Triangulation_default_data_structure_2<Ttraits,Vertex_base,Face_base> TDS ;
typedef CGAL::Triangulation_2< Ttraits , TDS> Triangulation_2;
typedef CGAL::Delaunay_triangulation_2< Ttraits , TDS> Delaunay_2;


CGAL_BEGIN_NAMESPACE


template<class triangulation_2>
class Drawable_voronoi_2: public  Drawable_object
{
private:
  triangulation_2 tr;
 
public:  
  

Drawable_voronoi_2(){type="voronoi_2";}

Drawable_voronoi_2(const triangulation_2 &tet,Color c, Style
			 sty=WIRE, Size s=5, Precision prec=15)
  {
    tr=tet;
    color = c; size=s;
      set_center();lind=0;type="voronoi_2";style=sty;precision=prec;
    }

  void set_center()
    {

      typename triangulation_2::Vertex_iterator vit;
      o_center[0]=0;o_center[1]=0;o_center[2]=0;
      for (vit=tr.finite_vertices_begin() ; vit != tr.vertices_end(); vit++) {
	o_center[0]= o_center[0] + vit->point().x();
	o_center[1]= o_center[1] + vit->point().y();
      }
      o_center[0] = o_center[0]/tr.number_of_vertices();
      o_center[1] = o_center[1]/tr.number_of_vertices();
    }

  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	set_color(color);
        glLineWidth(size);

	typedef typename triangulation_2::Edge_iterator Edge_iterator;
    
	Edge_iterator   it = tr.edges_begin(),
	  beyond = tr.edges_end();
	glBegin( GL_LINES );
	for ( ;it != beyond; ++it) {
	  CGAL::Object o = tr.dual(it);
	 typename  triangulation_2::Line l;
	 typename  triangulation_2::Ray r;
	 typename  triangulation_2::Segment s;
	  if (assign(s,o)) {
	    glVertex2f(s.source().x(),s.source().y() );
	    glVertex2f(s.target().x(),s.target().y());
	  }
	  if (assign(r,o)) {
	    glVertex2f(r.source().x(),r.source().y() );
	    glVertex2f(r.point(1000).x(),r.point(1000).y());
	  }
	  if (assign(l,o)) {
	    glVertex2f(l.point(-1000).x(),l.point(-1000).y());
	    glVertex2f(l.point(1000).x(),l.point(1000).y());
	  }
	}
	glEnd();
      }
      glEndList();
    }


 void add_point(double x, double y ,double z)
    {
      glDeleteLists(lind,1);
      lind=0;
      typedef typename triangulation_2::Point Point;
      tr.insert(Point(x,y));
      set_center();
    }
};

CGAL_END_NAMESPACE
