

#include "readtex.c" 



typedef CGAL::Point_3<CGAL::Cartesian<double> > point_3;


#define TABLE_TEXTURE "psy.rgb"





static GLubyte *Image = NULL;



static int ImgWidth, ImgHeight;
static GLenum ImgFormat;



void draw_tri_tex(double x1, double y1, double z1,double x2, double
		   y2, double z2,double x3, double y3, double z3)

{
  glPolygonMode(GL_FRONT,GL_FILL);
  glBegin(GL_TRIANGLES);

  vector<double> v1(3);
  vector<double> v2(3);
  v1 = CGAL::normal(x1,y1,z1,x2,y2,z2,x3,y3,z3);
       
  v2[0]= -v1[0];  v2[1]=-v1[1] ; v2[2]= -v1[2];
  glNormal3d(v1[0],v1[1],v1[2]);
  glTexCoord2f( 0.0, 0.0 );
  glVertex3f(x1,y1,z1);
  glNormal3d(v2[0],v2[1],v2[2]);
  glTexCoord2f( 1.0, 0.0 ); 
  glVertex3f(x2,y2,z2);
 glNormal3d(v1[0],v1[1],v1[2]);
 glTexCoord2f( 1,1 );
  glVertex3f(x3,y3,z3);
  glEnd();
}


void schrink_point(point_3 &p1, point_3 &p2, point_3 &p3, point_3 &p4)
{
double sx,sy,sz;
  double bx= (p1.x()+ p2.x()+ p3.x()+p4.x())/4;
  double by= (p1.y()+ p2.y()+ p3.y()+p4.y())/4;
  double bz= (p1.z()+ p2.z()+ p3.z()+p4.z())/4;
  sx= (bx - p1.x())/3;
  sy= (by - p1.y())/3;
  sz= (bz - p1.z())/3;
  p1 = point_3(p1.x()+sx,p1.y()+sy,p1.z()+sz);  
  sx= (bx - p2.x())/3;
  sy= (by - p2.y())/3;
  sz= (bz - p2.z())/3;
  p2 = point_3(p2.x()+sx,p2.y()+sy,p2.z()+sz);
  sx= (bx - p3.x())/3;
  sy= (by - p3.y())/3;
  sz= (bz - p3.z())/3;
  p3 = point_3(p3.x()+sx,p3.y()+sy,p3.z()+sz);
  sx= (bx - p4.x())/3;
  sy= (by - p4.y())/3;
  sz= (bz - p4.z())/3;
  p4 = point_3(p4.x()+sx,p4.y()+sy,p4.z()+sz);
}

CGAL_BEGIN_NAMESPACE

template<class tetrahedron>
class Drawable_tetrahedron_tex: public  Drawable_object
{
private:
  tetrahedron tr;

public:

  Drawable_tetrahedron_tex(){type="Tetrahedron";}
  Drawable_tetrahedron_tex(const tetrahedron &tet,Color c, Style
			 sty=WIRE, Size s=5, Precision prec=15)

    {
      tr=tet;
      color = c; size=s;
      set_center();lind=0;type="Tetrahedron";style=sty;precision=prec;
    }

  void set_center()
    {
      o_center[0]=(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x())/4; 
      o_center[1]=(tr[0].y()+tr[1].y()+tr[2].y()+tr[3].y())/4;
      o_center[2]=(tr[0].z()+tr[1].z()+tr[2].z()+tr[3].z())/4;
    }

void draw()
  {
    glEnable( GL_TEXTURE_2D );
  if(lind)
    glCallList(lind);
  else {
   lind = glGenLists(1);
   Image = LoadRGBImage( TABLE_TEXTURE, &ImgWidth, &ImgHeight, &ImgFormat
   		      );

      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, ImgWidth, ImgHeight,
                          GL_RGB, GL_UNSIGNED_BYTE, Image);
   //glTexImage2D(GL_TEXTURE_2D,0,3, ImgWidth, ImgHeight,0,GL_RGB, GL_UNSIGNED_//BYTE, Image);
     glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
   glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
   //   glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
   //   glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST
   //		    );
   // glEnable( GL_TEXTURE_2D );

   glNewList(lind,GL_COMPILE_AND_EXECUTE);
     glLineWidth(size);
     set_color(color);
       draw_tri_tex(tr[0].x(),tr[0].y(),tr[0].z(),tr[1].x(),tr[1].y(),
		     tr[1].z(),tr[2].x(),tr[2].y(),tr[2].z());
       draw_tri_tex(tr[0].x(),tr[0].y(),tr[0].z(),tr[1].x(),tr[1].y(),
		     tr[1].z(),tr[3].x(),tr[3].y(),tr[3].z());
       draw_tri_tex(tr[0].x(),tr[0].y(),tr[0].z(),tr[3].x(),tr[3].y(),
		     tr[3].z(),tr[2].x(),tr[2].y(),tr[2].z());
       draw_tri_tex(tr[1].x(),tr[1].y(),tr[1].z(),tr[2].x(),tr[2].y(),
		     tr[2].z(),tr[3].x(),tr[3].y(),tr[3].z());
     }

   glEndList();
 glDisable( GL_TEXTURE_2D );
  }

};


template<class triangulation_3>
class Drawable_triangulation_3: public  CGAL::Drawable_object
{
private:
  triangulation_3 tr;
 
public:  


Drawable_triangulation_3(){type="Tetrahedron";}


  Drawable_triangulation_3(const triangulation_3 &tet,Color c, Style
			 sty=WIRE, Size s=5, Precision prec=15)

    {
      tr=tet;
      color = c; size=s;
      set_center();lind=0;type="Triangulation_3";style=sty;precision=prec;
      //      Image = LoadRGBImage( TABLE_TEXTURE, &ImgWidth, &ImgHeight, &ImgFormat);
      //      int i =gluBuild2DMipmaps(GL_TEXTURE_2D, 3,ImgWidth,ImgHeight ,
      //                GL_RGB, GL_UNSIGNED_BYTE, Image);
      //	        glTexImage2D(GL_TEXTURE_2D,0,3, ImgWidth, ImgHeight,0,G// L_RGB, GL_UNSIGNED_BYTE, Image);

      //           glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
      //            glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
      //            glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      //            glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }

  void set_center()
    {

      typename triangulation_3::Vertex_iterator vit;
      o_center[0]=0;o_center[1]=0;o_center[2]=0;
      for (vit=tr.finite_vertices_begin() ; vit != tr.vertices_end(); vit++) {
      o_center[0]= o_center[0] + vit->point().x();
      o_center[1]= o_center[1] + vit->point().y();
      o_center[2]= o_center[2] + vit->point().z();
      }
      o_center[0] = o_center[0]/tr.number_of_vertices();
      o_center[1] = o_center[1]/tr.number_of_vertices();
      o_center[2] = o_center[2]/tr.number_of_vertices();
    }

  void draw()
    {
      //                glEnable( GL_TEXTURE_2D );
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	point_3 p0,p1,p2, p3;
	set_color(color);
	typename triangulation_3::Cell_iterator cit;
	for (cit = tr.finite_cells_begin(); cit != tr.cells_end(); cit++) 
	  {
	    p0=(cit->vertex(0))->point();
	    p1=(cit->vertex(1))->point();
	    p2=(cit->vertex(2))->point();
	    p3=(cit->vertex(3))->point();
	    schrink_point(p0,p1,p2,p3);

	    draw_tri_tex(p1.x(),p1.y(),p1.z(),p3.x(),p3.y(),p3.z(),p2.x(),p2.y(),p2.z());

	    draw_tri_tex(p0.x(),p0.y(),p0.z(),p2.x(),p2.y(),p2.z(),p3.x(),p3.y(),p3.z());
	    draw_tri_tex(p0.x(),p0.y(),p0.z(),p3.x(),p3.y(),p3.z(),p1.x(),p1.y(),p1.z());
	    draw_tri_tex(p0.x(),p0.y(),p0.z(),p1.x(),p1.y(),p1.z(),p2.x(),p2.y(),p2.z());
	  }
      }
      glEndList();
      //        glDisable( GL_TEXTURE_2D );       
    }

 void add_point(double x, double y ,double z)
    {
      glDeleteLists(lind,1);
      lind=0;
      typedef typename triangulation_3::Point Point;
      tr.insert(Point(x,y,z));
      set_center();
    }
};

CGAL_END_NAMESPACE
