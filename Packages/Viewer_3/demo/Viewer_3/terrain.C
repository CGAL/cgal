#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Viewer_stream.h>
#include "demo_win.h"
#ifndef V_UTILS
#include <CGAL/v_utils.h>
#endif 

using namespace CGAL;



typedef Homogeneous<CGAL::Gmpz>  Rp;
typedef Point_3<Rp> point3;
typedef Segment_3<Rp> segment;
typedef Tetrahedron_3<Rp> tetra;

typedef Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef Triangulation_vertex_base_2<Gt> Vb;
typedef Triangulation_face_base_2<Gt> Fb;
typedef Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef Delaunay_triangulation_2<Gt, Tds> Delaunay;




std::list<point3> lpt;


std::list<point3> set_z_zero(std::list<point3> lp)
{

  std::list<point3>::iterator it;
  std::list<point3> res;
  for (it=lp.begin(); it!=lp.end(); it++) 
    res.push_back(point3(it->hx(),it->hy(),-200,1));
  return res;
}

 





template<class triangulation_3>
class Drawable_triangulation_3: public  CGAL::Drawable_object
{
private:
  triangulation_3 tr;
 
public:  


Drawable_triangulation_3(){type="Triangulation_3";}


  Drawable_triangulation_3(const triangulation_3 &tet,Color c1, Color
			   c2=BLACK, Style sty=WIRE, Size s=3, Precision prec=10)

    {
      tr=tet;
      color = c1; col2 = c2; size=s;
      set_center();lind=0;type="Triangulation_3";style=sty;precision=prec;
    }

  void set_center()
    {

      typename triangulation_3::Vertex_iterator vit;
      o_center[0]=0;o_center[1]=0;o_center[2]=0;
      for (vit=tr.finite_vertices_begin() ; vit != tr.vertices_end(); vit++) {
      o_center[0]= o_center[0] + to_double(vit->point().x());
      o_center[1]= o_center[1] +  to_double(vit->point().y());
      o_center[2]= o_center[2] +  to_double(vit->point().z());
      }
      o_center[0] = o_center[0]/tr.number_of_vertices();
      o_center[1] = o_center[1]/tr.number_of_vertices();
      o_center[2] = o_center[2]/tr.number_of_vertices();
    }

  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	int no;
	typedef typename triangulation_3::Edge_iterator Edge_iterator;
	typedef typename triangulation_3::Vertex_handle Vertex_handle;
	typedef typename triangulation_3::Face_handle Face_handle;
	Vertex_handle v1, v2;
	Face_handle f;
	Edge_iterator it = tr.finite_edges_begin();
	Edge_iterator   beyond = tr.edges_end();

	glNewList(lind,GL_COMPILE_AND_EXECUTE);
        if (style==WIRE) {
	  set_color(color);
	  for ( ;it != beyond; ++it) {
	    f = (*it).first;
	    no = (*it).second;
	    v1 = f->vertex(f->ccw(no));
	    v2 = f->vertex(f->cw(no));
	    draw_tube(to_double(v1->point().x()),to_double(v1->point().y()),to_double(v1->point().z()),to_double(v2->point().x()),to_double(v2->point().y()),to_double(v2->point().z()),size, precision);
	  }
	}
	else {
	 typename triangulation_3::Face_iterator fit;
	 set_color(color);
	 glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	 std::vector<float> vp1(3), vp2(3), vp3(3);
	 for (fit=tr.finite_faces_begin(); fit!=tr.faces_end();fit++) {
	     vp1[0] =to_double(fit->vertex(0)->point().x());
	     vp1[1] =to_double(fit->vertex(0)->point().y());
	     vp1[2] =to_double(fit->vertex(0)->point().z());
	     vp2[0] =to_double(fit->vertex(1)->point().x());
	     vp2[1] =to_double(fit->vertex(1)->point().y());
	     vp2[2] =to_double(fit->vertex(1)->point().z());
	     vp3[0] =to_double(fit->vertex(2)->point().x());
	     vp3[1] =to_double(fit->vertex(2)->point().y());
	     vp3[2] =to_double(fit->vertex(2)->point().z());
	     draw_triangle(vp1[0],vp1[1],vp1[2]-1,vp2[0],vp2[1],vp2[2]-1,vp3[0],vp3[1],vp3[2]-1);
	   }
     
	 it = tr.finite_edges_begin();
	 set_color(col2);
	 glLineWidth(2);
	 glBegin(GL_LINES);
	 for ( ;it != beyond; ++it) {
	   f = (*it).first;
	   no = (*it).second;
	   v1 = f->vertex(f->ccw(no));
	   v2 = f->vertex(f->cw(no));
	   glVertex3f(to_double(v1->point().x()),to_double(v1->point().y()),to_double(v1->point().z()));
	   glVertex3f(to_double(v2->point().x()),to_double(v2->point().y()),to_double(v2->point().z()));
	 }
	 glEnd();
	}
	glEndList();
      }
    }

};







int main()
{
  
  point3 p0 = point3(50,50,0,1);
  point3 p1 = point3(100,50,0,1);
  point3 p2 = point3(150,50,10,1);
  point3 p3 = point3(200,50,20,1);
  point3 p4 = point3(261,50,30,1);
  point3 p5 = point3(300,50,50,1);
  point3 p6 = point3(341,50,40,1);
  point3 p7 = point3(400,50,30,1);
  point3 p8 = point3(450,50,20,1);
  point3 p9 = point3(500,50,10,1);

  point3 p10 = point3(50,100,0,1);
  point3 p11 = point3(100,100,0,1);
  point3 p12 = point3(155,100,0,1);
  point3 p13 = point3(200,100,20,1);
  point3 p14 = point3(250,100,40,1);
  point3 p15 = point3(320,100,50,1);
  point3 p16 = point3(350,101,70,1);
  point3 p17 = point3(400,100,50,1);
  point3 p18 = point3(420,104,20,1);
  point3 p19 = point3(500,100,0,1);

  point3 p20 = point3(50,150,0,1);
  point3 p21 = point3(100,150,0,1);
  point3 p22 = point3(150,150,0,1);
  point3 p23 = point3(230,150,50,1);
  point3 p24 = point3(250,150,40,1);
  point3 p25 = point3(300,150,20,1);
  point3 p26 = point3(350,150,30,1);
  point3 p27 = point3(400,150,0,1);
  point3 p28 = point3(450,150,10,1);
  point3 p29 = point3(500,150,0,1);

  point3 p30 = point3(50,200,0,1);
  point3 p31 = point3(100,200,0,1);
  point3 p32 = point3(150,200,20,1);
  point3 p33 = point3(230,200,50,1);
  point3 p34 = point3(250,200,40,1);
  point3 p35 = point3(320,204,20,1);
  point3 p36 = point3(350,200,20,1);
  point3 p37 = point3(400,200,50,1);
  point3 p38 = point3(450,200,20,1);
  point3 p39 = point3(500,200,0,1);

  point3 p40 = point3(50,250,0,1);
  point3 p41 = point3(100,250,0,1);
  point3 p42 = point3(154,250,20,1);
  point3 p43 = point3(220,255,40,1);
  point3 p44 = point3(270,250,60,1);
  point3 p45 = point3(306,250,100,1);
  point3 p46 = point3(330,250,50,1);
  point3 p47 = point3(350,258,40,1);
  point3 p48 = point3(450,250,20,1);
  point3 p49 = point3(500,250,0,1);

  point3 p50 = point3(50,300,0,1);
  point3 p51 = point3(100,300,0,1);
  point3 p52 = point3(150,300,0,1);
  point3 p53 = point3(230,300,20,1);
  point3 p54 = point3(250,300,60,1);
  point3 p55 = point3(320,300,80,1);
  point3 p56 = point3(370,300,40,1);
  point3 p57 = point3(400,300,20,1);
  point3 p58 = point3(450,300,0,1);
  point3 p59 = point3(500,300,0,1);

  point3 p60 = point3(70,350,0,1);
  point3 p61 = point3(100,350,0,1);
  point3 p62 = point3(150,350,10,1);
  point3 p63 = point3(210,350,20,1);
  point3 p64 = point3(260,350,40,1);
  point3 p65 = point3(300,350,50,1);
  point3 p66 = point3(350,350,50,1);
  point3 p67 = point3(400,350,40,1);
  point3 p68 = point3(450,350,0,1);
  point3 p69 = point3(500,350,0,1);

  point3 p70 = point3(50,400,0,1);
  point3 p71 = point3(100,400,0,1);
  point3 p72 = point3(150,400,10,1);
  point3 p73 = point3(200,400,0,1);
  point3 p74 = point3(250,400,30,1);
  point3 p75 = point3(300,400,50,1);
  point3 p76 = point3(350,400,30,1);
  point3 p77 = point3(400,400,10,1);
  point3 p78 = point3(450,400,0,1);
  point3 p79 = point3(500,400,0,1);

  point3 p80 = point3(50,450,0,1);
  point3 p81 = point3(100,450,0,1);
  point3 p82 = point3(150,450,0,1);
  point3 p83 = point3(200,450,0,1);
  point3 p84 = point3(250,450,20,1);
  point3 p85 = point3(300,450,20,1);
  point3 p86 = point3(350,450,10,1);
  point3 p87 = point3(400,450,0,1);
  point3 p88 = point3(450,450,0,1);
  point3 p89 = point3(500,450,0,1);

  point3 p90 = point3(50,500,0,1);
  point3 p91 = point3(100,500,0,1);
  point3 p92 = point3(150,500,0,1);
  point3 p93 = point3(200,500,0,1);
  point3 p94 = point3(260,500,10,1);
  point3 p95 = point3(300,500,0,1);
  point3 p96 = point3(370,500,0,1);
  point3 p97 = point3(400,500,0,1);
  point3 p98 = point3(450,500,0,1);
  point3 p99 = point3(500,500,0,1);


  lpt.push_back(p0);lpt.push_back(p1);lpt.push_back(p2);lpt.push_back(p3);
  lpt.push_back(p4);lpt.push_back(p5);lpt.push_back(p6);lpt.push_back(p7);
  lpt.push_back(p8);lpt.push_back(p9);

  lpt.push_back(p10);lpt.push_back(p11);lpt.push_back(p12);lpt.push_back(p13);
  lpt.push_back(p14);lpt.push_back(p15);lpt.push_back(p16);lpt.push_back(p17);
  lpt.push_back(p18);lpt.push_back(p19);
 
  lpt.push_back(p20);lpt.push_back(p21);lpt.push_back(p22);lpt.push_back(p23);
  lpt.push_back(p24);lpt.push_back(p25);lpt.push_back(p26);lpt.push_back(p27);
  lpt.push_back(p28);lpt.push_back(p29);

  lpt.push_back(p30);lpt.push_back(p31);lpt.push_back(p32);lpt.push_back(p33);
  lpt.push_back(p34);lpt.push_back(p35);lpt.push_back(p36);lpt.push_back(p37);
  lpt.push_back(p38);lpt.push_back(p39);

  lpt.push_back(p40);lpt.push_back(p41);lpt.push_back(p42);lpt.push_back(p43);
  lpt.push_back(p44);lpt.push_back(p45);lpt.push_back(p46);lpt.push_back(p47);
  lpt.push_back(p48);lpt.push_back(p49);

  lpt.push_back(p50);lpt.push_back(p51);lpt.push_back(p52);lpt.push_back(p53);
  lpt.push_back(p54);lpt.push_back(p55);lpt.push_back(p56);lpt.push_back(p57);
  lpt.push_back(p58);lpt.push_back(p59);

  lpt.push_back(p60);lpt.push_back(p61);lpt.push_back(p62);lpt.push_back(p63);
  lpt.push_back(p64);lpt.push_back(p65);lpt.push_back(p66);lpt.push_back(p67);
  lpt.push_back(p68);lpt.push_back(p69);

  lpt.push_back(p70);lpt.push_back(p71);lpt.push_back(p72);lpt.push_back(p73);
  lpt.push_back(p74);lpt.push_back(p75);lpt.push_back(p76);lpt.push_back(p77);
  lpt.push_back(p78);lpt.push_back(p79);

  lpt.push_back(p80);lpt.push_back(p81);lpt.push_back(p82);lpt.push_back(p83);
  lpt.push_back(p84);lpt.push_back(p85);lpt.push_back(p86);lpt.push_back(p87);
  lpt.push_back(p88);lpt.push_back(p89);

  lpt.push_back(p90);lpt.push_back(p91);lpt.push_back(p92);lpt.push_back(p93);
  lpt.push_back(p94);lpt.push_back(p95);lpt.push_back(p96);lpt.push_back(p97);
  lpt.push_back(p98);lpt.push_back(p99);

  std::list<point3>::iterator first=lpt.begin(), last=lpt.end();
  Drawable_points_set_3<std::list<point3>::iterator,point3>
    dlp(first,last,RED,FILL,5,15);

  std::list<point3> lplan=set_z_zero(lpt);
  first=lplan.begin();last=lplan.end();
  Drawable_points_set_3<std::list<point3>::iterator,point3>
    dlplan(first,last,GRAY,FILL,5,15);


  Delaunay dt;
  dt.insert(lpt.begin(),lpt.end());
  Drawable_triangulation_3<Delaunay> ddt(dt,ORANGE,BLACK,FILL);
  
  Delaunay dtp;
  dtp.insert(lplan.begin(),lplan.end());
  Drawable_triangulation_3<Delaunay> ddtp(dtp,PURPLE);

  CGAL::Viewer_3 W(500);
  W.init_window_thread(); 
  W.set_custom_panel(demo_panel);
  stop();
  W.add_drawable(&dlp,1);
  W.display();
  stop();
  W.add_drawable(&dlplan,2);
  W.display();
  stop();
  W.add_drawable(&ddtp,3);
  W.display();
  stop();
  W.add_drawable(&ddt,4);
  W.display();
  pthread_join(W.get_window_thread(), NULL);
}
