// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/GL_win.h
// revision      : $Revision$
//
// author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
// #include <config.h> // fltk's config
#include <FL/Fl.H>  // order of includes is important on Windows.
#include <FL/gl.h>
#include <FL/Fl_Gl_Window.H>
#include <CGAL/scene_graph.h>

#include <stdlib.h>
#include <fstream.h>
#include <iostream.h>

#ifdef USE_THREAD
#include <CGAL/Threads.h>
#endif

//POSTSCRIPT
#include "PS_Stream.C"
#include "PS_Stream_3.C"
#include "PS_edge_3.C"
#include "PS_facet_3.C"

#include <math.h>

typedef CGAL::Cartesian<double> D;
typedef CGAL::Bbox_3 PS_BBox3;
typedef CGAL::Direction_3< D > Direction;
typedef CGAL::Point_3< D > Point3;
typedef CGAL::Point_2< D > Point2;
typedef CGAL::Plane_3< D > Plane3;typedef CGAL::Direction_3< D > Direction;
typedef CGAL::Line_3< D > Line3;

CGAL_BEGIN_NAMESPACE
class GL_win : public Fl_Gl_Window {
  
public:
  typedef void (*Mouse_click)(int , int , int , GL_win *);
  
  typedef void (*Mouse_grab)(int , int , int , int , int ,int,int ,
			     GL_win *) ;

private:

  float Tclip;
  float Wclip;
  double clip0[4];
  double clip1[4];
  
  bool rotate;
  float zplan;
  float d_zplan;
  
  int *scale;
  Color bg_color;
  Scene_graph SCG;
  std::vector<double> add_point;
  
  double global_matrix[16];
  double last_rot[16];
  
  bool Projection;
  float deep;
  float angle;

  float proj_near;
  float proj_far;
  
  enum group_mode {all=0, selected, move, see};
  int group;
  enum mode {View=0, Insert, Slice, User};
  
  int choice;
  
  Mouse_click    M_P, M_R;
  
  Mouse_grab    M_G;
  
  float col_diff[4];
  
  bool light;
  float var;
  float shy;
  
  float Xlight;
  float Ylight;
  float Zlight;

  
  // private functions :
  int handle(int);
  
  void draw();
  
  void draw_scene();
  
  void m_grab(int,int , int ) ;
  
  void m_push(int , int , int ) ;
  
  void m_release(int , int , int ) ;
  
  void set_light(float, float, float, float*, float,float);

  std::vector<double> apply_last_rotation(std::vector<double>, double ,double , double);

  std::vector<double> project_on_plan(const std::vector<double> &, double ,double , double);

  void draw_plan(double, double, double);


public:

  std::vector<double> get_real_point(int,std::vector<double>);

  std::vector<double> get_point();

  GL_win(int ,int ,int ,int ,int *sc ,const char*);

  void reshape();

  void mouse_rotate(float,float);

  void reset();

  void reset_group(int);

  void set_angle(float);

  void set_deep(float);

  void set_projection(bool);

  void set_mouse_push_handler(Mouse_click);

  void set_mouse_release_handler(Mouse_click);

  void set_mouse_grab_handler(Mouse_grab);

  void set_mode(int);

  void add_zplan(int);

  void add_drawable(Drawable_object*, int);

  void remove_drawable(int , int);

  void free_drawable(int , int);

  Drawable_object* get_drawable(int, int);

  double* get_transform();

  void redraw();

  void make_visible(int);

  void change_visibility(int , int , bool );

  void add_point_to_object(int,int,std::vector<double>);

  void add_new_group();

  void change_group(int);
  
  void delete_group(int);

  void group_visible(bool , int );

  void group_rotate(float , float , int );

  Scene_graph* get_scene_graph();

  void clean_graph();

  void set_light_diff(float,float,float);

  void set_light_variation(float);

  void set_light_shy(float);

  void set_X_light_pos(float);

  void set_Y_light_pos(float);

  void set_Z_light_pos(float);

  //void draw_slice(double);

  void set_clip_planes(int);

  void set_clip_width(float,int);

  void set_clip_move(float,int);

  void set_bgcolor(Color);

  void switch_light();

  void set_light();

  void set_light(bool);

  void draw_ps();

  int get_group() {return group;}

};



void mouse_push(int x, int y, int but, GL_win *W) {std::cerr << "No action : default mouse function" << std::endl;}

void mouse_grab(int x_o, int y_o, int x_e, int y_e, int dx, int dy,
		int but, GL_win *W) 
{std::cerr << "No action : default mouse function" << std::endl;}



void mouse_release(int x, int y, int but, GL_win *W) 
{std::cerr << "No action : default mouse function" <<std:: endl;}




// -----------------  Private Part ----------------------------------

void GL_win::set_light(float x,float y , float z, float* diff, float
		       v, float s)
{
  GLfloat ambient[] = { 0.4, 0.4, 0.3, 1.0 };
  
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambient);
  GLfloat position[] = { x, y, z, 1.0 };
  glLightf(GL_LIGHT0,GL_CONSTANT_ATTENUATION, v);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);

  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  
  // define material reflective properties
  
  GLfloat ref[] = { 1.0,1.0,0.0};
  glMaterialfv(GL_FRONT, GL_SPECULAR,ref);
  glMaterialf(GL_FRONT, GL_SHININESS, s);
  glEnable(GL_COLOR_MATERIAL);
}

//#### POST SCRIPT : envoie le graph de scene dans le post script.
void GL_win::draw_ps()
{
#ifdef USE_THREAD
  pthread_mutex_lock(&Synchronizer::sgMutex);
#endif

  Scene_graph::iterator it;
  double* t= SCG.get_translation();
  double* rot;
  double rot_inv[16];

  glLoadIdentity();
  glPushMatrix();
  glTranslated(t[0],-t[1],0);
  glTranslatef(SCG.get_center(1),SCG.get_center(2),SCG.get_center(3));
  glMultMatrixd(SCG.get_rotation());
  rot = SCG.get_rotation();
  for (int i = 0; i < 16; i++) rot_inv[i] = rot[i];
  invert(rot_inv);
  glTranslatef(-SCG.get_center(1),-SCG.get_center(2),-SCG.get_center(3));

  Direction dir(rot_inv[8], rot_inv[9], rot_inv[10]);
  Direction light(0,0,1);
 
  int sca = 0;
  for (int i = 0; i < 6; i++) {
    if (abs(scale[i]) > sca) {
      sca = scale[i];
    }
  }

  PS_BBox3 bb3(-sca, -sca, -sca, sca, sca, sca);
  CGAL::PS_Stream_3 ps(bb3,dir,light,300,"Postcript.ps",CGAL::PS_Stream::READABLE_EPS);

  // on parcoure le graph de scene, groupe apres groupe.
  for (it=SCG.begin() ; it!=SCG.end() ; it++) {
    // On s'occupe de la matrice de transformation globale.
    it->group_to_ps(ps);
  }
  glPopMatrix();

  ps.display();

#ifdef USE_THREAD
  pthread_mutex_unlock(&Synchronizer::sgMutex);
#endif


}




void GL_win::draw_scene()
{
  static bool blink=false;
#ifdef USE_THREAD
  pthread_mutex_lock(&Synchronizer::sgMutex);
#endif
  Scene_graph::iterator it;
  double* t= SCG.get_translation();

  for (it=SCG.begin() ; it!=SCG.end() ; it++) {
    glLoadIdentity();
    glPushMatrix();
    glTranslated(t[0],-t[1],0);
    glTranslatef(SCG.get_center(1),SCG.get_center(2),SCG.get_center(3));
    glMultMatrixd(SCG.get_rotation());
    glTranslatef(-SCG.get_center(1),-SCG.get_center(2),-SCG.get_center(3));
    if (!group)
      glGetDoublev(GL_MODELVIEW_MATRIX,global_matrix);
    glMultMatrixd(it->get_group_translation());
    glMultMatrixd(it->get_group_rotation());
    switch (group) {
    case all:
      it->draw_group();
      break;
    case selected:
      it->draw_visible();
      break;
    case see:
      it->draw_invisible();
      if (blink)
	it->draw_visible();
      break;
    default:
      it->draw_group();
      break;
    }
    glPopMatrix();
  }
  blink=!blink;
#ifdef USE_THREAD
  pthread_mutex_unlock(&Synchronizer::sgMutex);
#endif
}




void GL_win::draw_plan(double x, double y, double z) {
  double v1[3]={0,0,1};
  double v2[3]={0,0,-1};
  glColor4f(1,1,1,0.5);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  glBegin(GL_POLYGON);
  glNormal3dv(v1);
  glVertex3f(x-(scale[1]-scale[0])/2,y-(scale[3]-scale[2])/2,z);
  glNormal3dv(v2);
  glVertex3f(x+(scale[1]-scale[0])/2,y-(scale[3]-scale[2])/2,z);
  glNormal3dv(v1);
  glVertex3f(x+(scale[1]-scale[0])/2,y+(scale[3]-scale[2])/2,z);
  glNormal3dv(v2);
  glVertex3f(x-(scale[1]-scale[0])/2,y+(scale[3]-scale[2])/2,z);
  glEnd();
}


std::vector<double> GL_win::apply_last_rotation(std::vector<double> point, double x
						,double y, double z)
{
  std::vector<double> r_point(3);
  r_point=apply_mat(last_rot, point[0]-x, point[1]-y, point[2]-z);
  r_point=translate(r_point,x,y,z);
  return r_point;
}


std::vector<double> GL_win::project_on_plan(const std::vector<double> &point, double x
				 ,double y, double z)
{
  std::vector<double> v1 = apply_mat(SCG.get_rotation(),-x,-y,SCG.get_center(3)+zplan-z);
  std::vector<double> v2 = apply_mat(SCG.get_rotation(),100-x,-y,SCG.get_center(3)+zplan-z);
  std::vector<double> v3 = apply_mat(SCG.get_rotation(),-x,100-y,SCG.get_center(3)+zplan-z);
  std::vector<double> vr(3);
  v1[0] = v1[0]+x;
  v2[0] = v2[0]+x;
  v3[0] = v3[0]+x;

  v1[1] = v1[1]+y;
  v2[1] = v2[1]+y;
  v3[1] = v3[1]+y;

  v1[2] = v1[2]+z;
  v2[2] = v2[2]+z;
  v3[2] = v3[2]+z;

  std::vector<double> plan = compute_plan(v1[0],v1[1],v1[2],
				     v2[0],v2[1],v2[2],
				     v3[0],v3[1],v3[2]);

  if (d_zplan) {
    double den=sqrt(pow(plan[0],2) + pow(plan[1],2) +
		    pow(plan[2],2) );
    vr=translate(point,d_zplan*(plan[0]/den),d_zplan*(plan[1]/den),
		 d_zplan*(plan[2]/den));
    d_zplan=0;
  } 
  else {
    vr=point;
    vr[2]=-intersect_plan(plan,point[0],point[1]);
  }
  return vr;
}





void GL_win::draw()
{

  if (!valid()) 
    reshape();

  glPushMatrix();
  glDrawBuffer(GL_FRONT_AND_BACK);
  glClearColor((float)
	       bg_color.red()/255,(float)bg_color.green()/255,(float)
	       bg_color.blue()/255,bg_color.alpha());

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPopMatrix();

  if (light)
    set_light(Xlight,Ylight,Zlight,col_diff,var,shy);
  double z_c;
  if (choice == Insert) {
    z_c = SCG.get_center(3);
    double* t = SCG.get_translation();
    
    glPushMatrix();
    glTranslated(t[0],-t[1],0);
    if (rotate) 
      add_point=apply_last_rotation(add_point,SCG.get_center(1),SCG.get_center(2),z_c);
    else 
      add_point=project_on_plan(add_point,SCG.get_center(1),SCG.get_center(2),z_c);
    
    
    
    glTranslatef(add_point[0],add_point[1],add_point[2]);
    glColor4f(1.0,1.0,1.0,1);
    GLUquadricObj *q= gluNewQuadric();
    gluQuadricNormals(q, GL_SMOOTH);
    gluSphere(q,(float) (scale[5]-scale[4]/100) + 1,10,8);
    gluDeleteQuadric(q);
    glPopMatrix();
  }

  if (choice == Slice) {
    glClipPlane(GL_CLIP_PLANE0,clip0);
    glClipPlane(GL_CLIP_PLANE1,clip1);
  }

  draw_scene();
  if (choice == Insert) { 
    glLoadIdentity();
    glPushMatrix();
    glMultMatrixd(global_matrix);
    draw_plan(SCG.get_center(1),SCG.get_center(2), SCG.get_center(3)+zplan); 
    glPopMatrix();
  } 
  
  glFlush();
  rotate=false;
}


void GL_win::m_grab(int dx,int dy, int but) 
{
  double mx,my;
  switch(but) {
  case 1:
    if (group==move) {
      Scene_graph::iterator it;
      for (it = SCG.begin() ; it < SCG.end() ; it++) 
	if (it->group_visible()) 
 	  it->set_translation(dx,dy,global_matrix);
      
    }
    else {
      mx = (double) dx*( (scale[1] - scale[0]) /500.0);
      my = (double) dy*( (scale[3] - scale[2]) /500.0);
      SCG.add_translation(mx,my);
    }
    
    break;
  case 2: 
    if (group==move) {
      Scene_graph::iterator it;
      for (it = SCG.begin() ; it < SCG.end() ; it++) 
	if (it->group_visible()) 
	  it->set_rotation(dx,dy,global_matrix);
      //,SCG.get_center(1),SCG.get_center(2), SCG.get_center(3));
      
    }
    else
      mouse_rotate(dx,dy); 
    break;
  case 3:
    if (choice == Insert) {
      mx = (double) dx*( (scale[1] - scale[0]) /500.0);
      my = (double) dy*( (scale[3] - scale[2]) /500.0);
      add_point[0]=add_point[0] + mx; 	
      add_point[1]=add_point[1] - my;
    }
    break;
  }
}

void GL_win::m_push(int x, int y, int but) {}

void GL_win::m_release(int x, int y, int but) {}


int GL_win::handle(int event)
{
  static int button;
  static int x1, y1, dx, dy;
  int x2, y2;

  switch (event) { 
  case FL_PUSH:
    button = Fl::event_button();
    dx = Fl::event_x(); 
    dy = Fl::event_y();
    x1 = Fl::event_x(); 
    y1 = Fl::event_y();

    if (choice == User) 
      M_P(x1,h()-y1,button,this);
    else  
      m_push(x1,y1,button);

    return 1;
  case FL_DRAG:
    x2 = Fl::event_x();
    y2 = Fl::event_y();
    if (choice == User) 
      M_G(x1,y1,x2,y2,x2-dx,y2-dy,button,this);
    else
      m_grab(x2-dx,y2-dy,button);
    dx=x2; dy=y2;
    //   redraw();
    
    return 1;
  case FL_RELEASE:
    if (choice == User) 
      M_R(x1,h()-y1,button,this);
    else  
      m_release(x1,h()-y1,button);
        redraw();
  }

  return 0;
}




// ------------------------- Public Part ------------------------



GL_win::GL_win(int x,int y,int w,int h,int *sc, const char *l=0)
  : Fl_Gl_Window(x,y,w,h,l), add_point(3)
{
  bg_color=BLACK; bg_color.set_alpha(255);
  angle=40;Projection=true; M_P = mouse_push; 
  M_R =mouse_release; M_G = mouse_grab;choice = 0; group=0;rotate=false;
  zplan=0; col_diff[0]=1; col_diff[1]=1; col_diff[2]=1;
  col_diff[3]=0.5; 
  var=1; shy=50; Xlight = 0; Ylight = 0; 
  Zlight =1000; 
  
  scale = sc;
  
  deep = 2*(scale[5] - scale[4]); 
  Tclip=0;Wclip=(scale[5]-scale[4])/5;
  clip0[0]=0; clip0[1]=0 ; clip0[2]=1 ; clip0[3]= -Tclip + Wclip/2;
  clip1[0]=0; clip1[1]=0 ; clip1[2]=-1 ; clip1[3]= Tclip +
					   Wclip/2;
  light=true;
  Fl::gl_visual(FL_RGB8);
}



void GL_win::add_drawable(Drawable_object* obj, int i)     
{
  SCG.add_drawable(obj,i);
}
    
void GL_win::redraw() {draw();swap_buffers();}


void GL_win::reshape()
{
  static bool first;
  
  glViewport(0, 0, w(), h());
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  float fact = (float) ((scale[5]-scale[4])/100) + 0.1;
  
  if (Projection) {
    if (deep>2*(scale[5] - scale[4])) {
      proj_near = (deep -2*(scale[5] - scale[4]))+fact;
      proj_far = (deep +2*(scale[5] - scale[4]));
    } else {
      proj_near = fact;
      proj_far = 4*(scale[5] - scale[4])+deep;
    }
    gluPerspective(angle,(float) w()/h(), proj_near, proj_far);
   
    glTranslatef(-(scale[0]+scale[1])/2,
		 -(scale[2]+scale[3])/2,-deep);
    

  }
  else {
    glOrtho( (3*scale[0]-scale[2])/2, (3*scale[1]-scale[3])/2, 
	     (3*scale[2]-scale[0])/2, (3*scale[3]-scale[1])/2,
	     (3*scale[4]-scale[5])/2, (3*scale[5]-scale[4])/2);
  }
  

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();
  if (!first) {
    glGetDoublev(GL_MODELVIEW_MATRIX, global_matrix);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    first=true;
    set_light();
  }
}


void GL_win::mouse_rotate(float x, float y)
{

  if ( (-5 < x) && (x < 5) )
    x=x/5 ; 
  if ( ((-10 < x) && (x <= -5)) || ((x >= 5) && (x < 10)))
    x=x/4;
  if ( ((-20 < x) && (x <= -10)) || ((x >= 10) && (x < 20)))
    x=x/3;
  if ( (x <= -20) || (x >= 20) )
    x=x/2;

  if ( (-5 < y) && (y < 5) )
    y=y/5 ; 
  if ( ((-10 < y) && (y <= -5)) || ((y >= 5) && (y < 10)))
    y=y/4;
  if ( ((-20 < y) && (y <= -10)) || ((y >= 10) && (y < 20)))
    y=y/3;
  if ( (y <= -20) || (y >= 20) )
    y=y/2;
  rotate=true;
  glLoadIdentity();
  glPushMatrix();
  glRotated(y,1,0,0);
  glRotated(x,0,1,0);
  glGetDoublev(GL_MODELVIEW_MATRIX,last_rot);
  glMultMatrixd(SCG.get_rotation());
  SCG.set_rotation();
  glPopMatrix();

}





 

void GL_win::reset()
{
 glLoadIdentity();
 glPushMatrix();
 glGetDoublev(GL_MODELVIEW_MATRIX, last_rot);
 glGetDoublev(GL_MODELVIEW_MATRIX, global_matrix);
 SCG.set_rotation();
 glPopMatrix();
 SCG.set_translation(0,0);
 Scene_graph::iterator it = SCG.begin();
 for(; it!=SCG.end(); it++)
   it->set_identity();
 
}

void GL_win::set_projection(bool p)
{
 Projection = p;
}

void GL_win::set_angle(float a) 
{
 angle=a;
}

void GL_win::set_deep(float d) 
{
  deep=d;
}


void GL_win::set_mouse_push_handler(Mouse_click hand)
{
  M_P = hand;
}

void GL_win::set_mouse_release_handler(Mouse_click hand)
{
  M_R = hand;
}

void GL_win::set_mouse_grab_handler(Mouse_grab hand)
{
  M_G = hand;
}

void GL_win::set_mode(int i)
{
  choice=i;
  if (choice == Insert) {
    glEnable(GL_BLEND); 
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    zplan=0;
    add_point[0] = SCG.get_center(1);
    add_point[1] = SCG.get_center(2);
    add_point[2] = SCG.get_center(3);
  }
  else 
    glDisable(GL_BLEND);
  if (choice == Slice) {
   glEnable(GL_CLIP_PLANE0);
   glEnable(GL_CLIP_PLANE1);
  }
  else {
   glDisable(GL_CLIP_PLANE0);
   glDisable(GL_CLIP_PLANE1);
  }
  redraw();
}

void GL_win::add_zplan(int z)
{
  d_zplan  = z- zplan;
  if (z == 0) {
    zplan = 0;
    d_zplan  = 0;
  }
  else
    zplan=z;
  redraw();
}

void GL_win::make_visible(int i)
{
  SCG.make_visible(i);
}


void GL_win::change_visibility(int i, int j, bool b)
{
  SCG.change_visibility(i,j,b);
}

std::vector<double> GL_win::get_point()
{
  return add_point;
}

void GL_win::add_point_to_object(int i, int j, std::vector<double> v )
{
 std::vector<double> p;
  p= get_real_point(i,v);
  SCG.add_point_to_object(i, j, p[0], p[1], p[2]);
}


std::vector<double> GL_win::get_real_point(int i, std::vector<double> v )
{
  double* t=SCG.get_translation();
  double* g_mat_t;
  double* g_mat_r;
  double res_mat[16];
  std::vector<double> res;
    g_mat_t=SCG.get_group_translation(i);
    g_mat_r=SCG.get_group_rotation(i);
    if ((g_mat_t != NULL) &&  (g_mat_r != NULL)){
      glLoadMatrixd(global_matrix);
      glPushMatrix();
      glMultMatrixd(g_mat_t);
      glMultMatrixd(g_mat_r);
      glGetDoublev(GL_MODELVIEW_MATRIX,res_mat);
      glPopMatrix();
      invert(res_mat);
    }
    else 
      invert(global_matrix,res_mat);
    if (v==add_point)
      res=apply_mat(res_mat, v[0]+t[0],v[1]-t[1],v[2]);
    else
      res=apply_mat(res_mat, v[0],v[1],v[2]);
      
  //    add_point=apply_mat(global_matrix, SCG.get_center(1),
  //			SCG.get_center(2),SCG.get_center(3)+ zplan);
    
  //    res[0] += t[0];
  //    res[1] -= t[1];   

  return res;
}

Drawable_object* GL_win::get_drawable(int gr, int i)
{
  return(SCG.get_drawable(gr, i));
}

void GL_win::remove_drawable(int g, int i)
{
  SCG.remove_drawable(g,i);
}

void GL_win::free_drawable(int g, int i)
{
  SCG.free_drawable(g,i);
}

void GL_win::add_new_group()
{
  SCG.add_new_group();
}

void GL_win::change_group(int c)
{
 group=c;
}

void GL_win::group_visible(bool b, int g)
{
  SCG.group_visible(b,g);
}

void GL_win::delete_group(int i)
{
  SCG.delete_group(i);
}

Scene_graph* GL_win::get_scene_graph()
{
  return &SCG;
}

void GL_win::clean_graph()
{
  SCG.clean_graph();
}

void GL_win::reset_group(int gr)
{
  int i=1;
  Scene_graph::iterator it = SCG.begin();
  while (i!=gr) {
    it++;
    i++;
  }
  it->set_identity();
}
  
void GL_win::set_light_diff(float r , float g , float b)
{
  col_diff[0]=r; col_diff[1]=g; col_diff[2]=b; 
  set_light();
}

void GL_win::set_light_variation(float v)
{
 var=v;
 set_light();
}

void GL_win::set_light_shy(float s)
{
 shy=s;
 set_light();
}

void GL_win::set_X_light_pos(float x)
{
  Xlight = x;
  set_light();
}

void GL_win::set_Y_light_pos(float y)
{
  Ylight = y;
  set_light();
}

void GL_win::set_Z_light_pos(float z)
{
  Zlight = z;
  set_light();
}

void GL_win::set_clip_planes(int v)
{
  switch(v) {
  case 0:
    clip0[0]=0; clip0[1]=0 ; clip0[2]=1 ; clip0[3]= -Tclip + Wclip/2;
    clip1[0]=0; clip1[1]=0 ; clip1[2]=-1 ; clip1[3]= Tclip + Wclip/2;
    break;
  case 1:
    clip0[0]=0; clip0[1]=1 ; clip0[2]=0 ; clip0[3]= -Tclip + Wclip/2;
    clip1[0]=0; clip1[1]=-1 ; clip1[2]=0 ; clip1[3]= Tclip + Wclip/2;
    break;
  case 2:
    clip0[0]=1; clip0[1]=0 ; clip0[2]=0 ; clip0[3]= -Tclip + Wclip/2;
    clip1[0]=-1; clip1[1]=0 ; clip1[2]=0 ; clip1[3]= Tclip + Wclip/2;
    break;

  }


}

void GL_win::set_clip_width(float v, int c)
{
  Wclip=v;
  set_clip_planes(c);
}

void GL_win::set_clip_move(float m,int c)
{
  Tclip =m ;
 set_clip_planes(c);
}

void GL_win::set_bgcolor(Color c)
{
  c.set_alpha(255);
  bg_color=c;
}

void GL_win::switch_light()
{
  light =!light;
  set_light();
}

void GL_win::set_light(bool b)
{
  if (light!=b) {
    light=b;
    set_light();
  }
    
}

void GL_win::set_light()
{
  if (light)
    set_light(Xlight,Ylight,Zlight,col_diff,var,shy);
  else {
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_COLOR_MATERIAL);
  }
}

double* GL_win::get_transform()
{
  return global_matrix;
}


CGAL_END_NAMESPACE
