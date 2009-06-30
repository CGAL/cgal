#include "Scene.h"

void Scene::init() {
  // undo from QGLViewer internal initializeGL function
  glDisable(GL_COLOR_MATERIAL);

  // camera
  ui->viewer->camera()->setPosition(Vec(0.5,0.5,2.7));
  ui->viewer->camera()->lookAt(Vec(0.5,0.5,0.5));

  // scene inits
  ui->viewer->setSceneCenter(qglviewer::Vec(0.5,0.5,0.5));
  ui->viewer->setSceneRadius(2.0);
  ui->viewer->setBackgroundColor(Qt::white);
  ui->viewer->setForegroundColor(Qt::red);
  pQuadric = gluNewQuadric();

  // OpenGL inits
  glPointSize(10.0);
  glLineWidth(1.0);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  l_triangulation = glGenLists(1);
  l_domain = glGenLists(2);

  // Scene OpenGL state
  gl_draw_domain();
  init_scene(EMPTY);
}

// Draws the triangulation
void Scene::draw() {

  // Draw the triangulation itself that is stored in the list.
  glCallList(l_triangulation);

  if (flying_ball) {
    change_material(materials[FLYING_BALL_COLOR]);
    if (wireframe) {
      glBegin(GL_POINTS);
      glVertex3f(moving_point.x(),moving_point.y(),moving_point.z());
      glEnd();
    } else {
    // draw the moving point
      glPushMatrix();
        glTranslated(moving_point.x(),moving_point.y(),moving_point.z());
        gluSphere(pQuadric, 0.02,15,15);
        glFlush();
      glPopMatrix();
    }
  }

  // draw the domain (unit square / cube)
  if (ddomain) glCallList(l_domain);

  glDisable(GL_LIGHTING);
  if (dlocate) gl_draw_location();
  if (dconflict) gl_draw_conflict();
  if (!wireframe) glEnable(GL_LIGHTING);
  glFlush();
}

void Scene::load_points(const QString& fileName) {
  p3dt.clear();
  std::vector<Point> points;
  std::ifstream ifs(fileName.toAscii().data() );
  std::copy(std::istream_iterator<Point>(ifs), 
  std::istream_iterator<Point>(),
  std::back_inserter(points));
  std::random_shuffle(points.begin(), points.end());
  p3dt.insert(points.begin(), points.end());
  Vertex_iterator vit = p3dt.vertices_begin();

  make_draw_list();

  QString snv;
  int nv = p3dt.number_of_vertices();
  snv.setNum(nv);
  emit message(QString("|V| = ") + snv, 0);

  draw();
}

// update the position of the moving point
void Scene::update_position()
{
  double x = moving_point.x() +0.01023;
  double y = moving_point.y() +0.003123;
  double z = (in_plane ? 0.0 : moving_point.z() +0.02567);
  if(x>1.)x-=1.;
  if(y>1.)y-=1.;
  if(z>1.)z-=1.;
  moving_point = Point(x,y,z);
  // TODO try to find a better possibility
  ui->viewer->update();
}

void Scene::make_draw_list()
{
  // Prepare set of segments
  Segment_set segments_to_draw;
  primitives_from_geom_it(segments_to_draw);
  if (cube_clipping && !two_color_clipping) segment_clipping(segments_to_draw);

  // Create new list
  glNewList(l_triangulation, GL_COMPILE); 

  // Draw vertices
  change_material(materials[VERTEX_COLOR]);
  if (wireframe) glBegin(GL_POINTS);
  for (Point_iterator pit = p3dt.periodic_points_begin(it_type) ;
       pit != p3dt.periodic_points_end(it_type) ; pit++)
    gl_draw_vertex(pit->first, 0.02);
  if (wireframe) glEnd();

  // Draw segments
  if (wireframe) glBegin(GL_LINES);
  if (cube_clipping && two_color_clipping) {
      change_material(materials[CLIPPING_COLOR]);
      segment_2color_clipping(segments_to_draw);
  }
  change_material(materials[EDGE_COLOR]);
  for (Segment_set::iterator it = segments_to_draw.begin() ;
       it != segments_to_draw.end(); it++)
    gl_draw_edge(it->source(),it->target(),0.005);
  if (wireframe) glEnd();

  glFlush();
  glEndList();

  // TODO: viewer should be specialized and then the respective
  // viewer->paintGL should be called.
  ui->viewer->update();
}

// some initialization templates
void Scene::init_scene(Init ID) {
  bool temp_flags[] = {dlocate, dconflict};
  dlocate = false;
  dconflict = false;
  p3dt.clear();
  make_draw_list();
  RandPts rp(0.5);
  Point pt2;
  switch (ID) {
  case GRID:
    p3dt.insert_dummy_points();
    break;
  case SINGLE:
    p3dt.insert(Point(0.3,0.4,0.5));
    break;
  case PLANE:
    for (int i=0 ; i<10 ; i++) {
      pt2 = *rp+Vector(0.5,0.5,0.5);
      rp++;
      p3dt.insert(Point(pt2.x(),pt2.y(),0.0));
    }
    break;
  case RANDOM:
    do {
      p3dt.insert(*rp+Vector(0.5,0.5,0.5));
      rp++;
    } 
    while (p3dt.number_of_vertices()<30);
  default:
    break;
  }
  dlocate = temp_flags[0];
  dconflict = temp_flags[1];
  make_draw_list();
}

// Draw a unit square in the x/y-plane
void Scene::gl_draw_square() {
  glNewList(l_domain, GL_COMPILE);
  change_material(materials[DOMAIN_COLOR]);
  if (!wireframe) {
    glPushMatrix();
    glTranslated(0.5,0.5,0.5);
    glRotated(90,1.0,0.0,0.0);
    glTranslated(-0.5,-0.5,-0.5);
    for (float x=0.0; x<2.0; x+=1.0) {
      for (float y=0.0; y<1.0; y+=1.0) {
	glPushMatrix();        
	glTranslated(x,y,0.0);
	gluCylinder(pQuadric,0.01,0.01,1.0,25,2);
	glPopMatrix();
      }
    }
    glPopMatrix();
    glPushMatrix();
    glTranslated(0.5,0.5,0.5);
    glRotated(90,0.0,1.0,0.0);
    glTranslated(-0.5,-0.5,-0.5);
    for (float x=1.0; x<2.0; x+=1.0) {
      for (float y=0.0; y<2.0; y+=1.0) {
	glPushMatrix();        
	glTranslated(x,y,0.0);
	gluCylinder(pQuadric,0.01,0.01,1.0,25,2);
	glPopMatrix();
      }
    }
    glPopMatrix();
  } else {
    glBegin(GL_LINES);
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(1.0,0.0,0.0);
    glVertex3d(0.0,1.0,0.0);
    glVertex3d(1.0,1.0,0.0);
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,1.0,0.0);
    glVertex3d(1.0,0.0,0.0);
    glVertex3d(1.0,1.0,0.0);
    glEnd();
  }
  glEndList();
}

// Draw a unit cube
void Scene::gl_draw_cube() {
  glNewList(l_domain, GL_COMPILE);
  change_material(materials[DOMAIN_COLOR]);
  if (!wireframe) {
    glPushMatrix();        
    for (float x=0.0; x<2.0; x+=1.0) {
      for (float y=0.0; y<2.0; y+=1.0) {
	glPushMatrix();        
	glTranslated(x,y,0.0);
	gluCylinder(pQuadric,0.01,0.01,1.0,25,2);
	glPopMatrix();
      }
    }
    glPopMatrix();
    glPushMatrix();
    glTranslated(0.5,0.5,0.5);
    glRotated(90,1.0,0.0,0.0);
    glTranslated(-0.5,-0.5,-0.5);
    for (float x=0.0; x<2.0; x+=1.0) {
      for (float y=0.0; y<2.0; y+=1.0) {
	glPushMatrix();        
	glTranslated(x,y,0.0);
	gluCylinder(pQuadric,0.01,0.01,1.0,25,2);
	glPopMatrix();
      }
    }
    glPopMatrix();
    glPushMatrix();
    glTranslated(0.5,0.5,0.5);
    glRotated(90,0.0,1.0,0.0);
    glTranslated(-0.5,-0.5,-0.5);
    for (float x=0.0; x<2.0; x+=1.0) {
      for (float y=0.0; y<2.0; y+=1.0) {
	glPushMatrix();        
	glTranslated(x,y,0.0);
	gluCylinder(pQuadric,0.01,0.01,1.0,25,2);
	glPopMatrix();
      }
    }
    glPopMatrix();
  } else {
    glBegin(GL_LINES);
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(1.0,0.0,0.0);
    glVertex3d(0.0,1.0,0.0);
    glVertex3d(1.0,1.0,0.0);
    glVertex3d(0.0,0.0,1.0);
    glVertex3d(1.0,0.0,1.0);
    glVertex3d(0.0,1.0,1.0);
    glVertex3d(1.0,1.0,1.0);
  
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,1.0,0.0);
    glVertex3d(1.0,0.0,0.0);
    glVertex3d(1.0,1.0,0.0);
    glVertex3d(0.0,0.0,1.0);
    glVertex3d(0.0,1.0,1.0);
    glVertex3d(1.0,0.0,1.0);
    glVertex3d(1.0,1.0,1.0);

    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,0.0,1.0);
    glVertex3d(1.0,0.0,0.0);
    glVertex3d(1.0,0.0,1.0);
    glVertex3d(0.0,1.0,0.0);
    glVertex3d(0.0,1.0,1.0);
    glVertex3d(1.0,1.0,0.0);
    glVertex3d(1.0,1.0,1.0);
    glEnd();
  }
  glEndList();
}

// get the offset that is common to all points of a triangle in the
// triangulation
inline void Scene::get_tri_offsets(const Cell_handle ch, int i,
    Offset &off0, Offset &off1, Offset &off2) const {
  off0 = p3dt.get_offset(ch,(i+1)&3);
  off1 = p3dt.get_offset(ch,(i+2)&3);
  off2 = p3dt.get_offset(ch,(i+3)&3);
  if (it_type == P3DT::UNIQUE || it_type == P3DT::UNIQUE_COVER_DOMAIN) {
    int diff_offx = (std::min)((std::min)(off0.x(),off1.x()),off2.x());
    int diff_offy = (std::min)((std::min)(off0.y(),off1.y()),off2.y());
    int diff_offz = (std::min)((std::min)(off0.z(),off1.z()),off2.z());
    Offset diff_off(diff_offx, diff_offy, diff_offz);
    off0 -= diff_off;
    off1 -= diff_off;
    off2 -= diff_off;
  }
}

// get the offset that is common to all points of a tetrahedron in the
// triangulation
inline void Scene::get_tet_offsets(const Cell_handle ch,
    Offset &off0, Offset &off1, Offset &off2, Offset &off3) const {
  off0 = p3dt.get_offset(ch,0);
  off1 = p3dt.get_offset(ch,1);
  off2 = p3dt.get_offset(ch,2);
  off3 = p3dt.get_offset(ch,3);
  if (it_type == P3DT::UNIQUE || it_type == P3DT::UNIQUE_COVER_DOMAIN) {
    int diff_offx = (std::min)((std::min)(off0.x(),off1.x()),
			     (std::min)(off2.x(),off3.x()));
    int diff_offy = (std::min)((std::min)(off0.y(),off1.y()),
			     (std::min)(off2.y(),off3.y()));
    int diff_offz = (std::min)((std::min)(off0.z(),off1.z()),
			     (std::min)(off2.z(),off3.z()));
    Offset diff_off(diff_offx, diff_offy, diff_offz);
    off0 -= diff_off;
    off1 -= diff_off;
    off2 -= diff_off;
    off3 -= diff_off;
  }
}

// return an integer that encodes the translations which have to be
// applied to the triangle to draw
inline int Scene::get_tri_drawing_offsets(const Cell_handle ch, int i) const {
  Offset off0, off1, off2;
  // if drawing boundary cells multiply is not activated then there is
  // nothing to do.
  switch( it_type ) {
  case P3DT::UNIQUE_COVER_DOMAIN:
    get_tri_offsets(ch,i,off0,off1,off2);
    break;
  case P3DT::STORED_COVER_DOMAIN:
    off0 = p3dt.int_to_off(ch->offset((i+1)&3));
    off1 = p3dt.int_to_off(ch->offset((i+2)&3));
    off2 = p3dt.int_to_off(ch->offset((i+3)&3));
    break;
  default:
    return 0;
  }

  CGAL_assertion(off0.x() == 0 || off0.x() == 1);
  CGAL_assertion(off0.y() == 0 || off0.y() == 1);
  CGAL_assertion(off0.z() == 0 || off0.z() == 1);
  CGAL_assertion(off1.x() == 0 || off1.x() == 1);
  CGAL_assertion(off1.y() == 0 || off1.y() == 1);
  CGAL_assertion(off1.z() == 0 || off1.z() == 1);
  CGAL_assertion(off2.x() == 0 || off2.x() == 1);
  CGAL_assertion(off2.y() == 0 || off2.y() == 1);
  CGAL_assertion(off2.z() == 0 || off2.z() == 1);
    
  int offx = ( ((off0.x() == 0 && off1.x() == 0 && off2.x() == 0)
		|| (off0.x() == 1 && off1.x() == 1 && off2.x() == 1)) ? 0 : 1);
  int offy = ( ((off0.y() == 0 && off1.y() == 0 && off2.y() == 0)
		|| (off0.y() == 1 && off1.y() == 1 && off2.y() == 1)) ? 0 : 1);
  int offz = ( ((off0.z() == 0 && off1.z() == 0 && off2.z() == 0)
		|| (off0.z() == 1 && off1.z() == 1 && off2.z() == 1)) ? 0 : 1);
    
  return( 4*offx + 2*offy + offz );
}

// return an integer that encodes the translations which have to be
// applied to the tetrahedron to draw
inline int Scene::get_tet_drawing_offsets(const Cell_handle ch) const {
  Offset off0, off1, off2, off3;
  // if drawing boundary cells multiply is not activated then there is
  // nothing to do.
  switch( it_type ) {
  case P3DT::UNIQUE_COVER_DOMAIN:
    get_tet_offsets(ch,off0,off1,off2,off3);
    break;
  case P3DT::STORED_COVER_DOMAIN:
    off0 = p3dt.int_to_off(ch->offset(0));
    off1 = p3dt.int_to_off(ch->offset(1));
    off2 = p3dt.int_to_off(ch->offset(2));
    off3 = p3dt.int_to_off(ch->offset(3));
    break;
  default:
    return 0;
  }

  CGAL_assertion(off0.x() == 0 || off0.x() == 1);
  CGAL_assertion(off0.y() == 0 || off0.y() == 1);
  CGAL_assertion(off0.z() == 0 || off0.z() == 1);
  CGAL_assertion(off1.x() == 0 || off1.x() == 1);
  CGAL_assertion(off1.y() == 0 || off1.y() == 1);
  CGAL_assertion(off1.z() == 0 || off1.z() == 1);
  CGAL_assertion(off2.x() == 0 || off2.x() == 1);
  CGAL_assertion(off2.y() == 0 || off2.y() == 1);
  CGAL_assertion(off2.z() == 0 || off2.z() == 1);
  CGAL_assertion(off3.x() == 0 || off3.x() == 1);
  CGAL_assertion(off3.y() == 0 || off3.y() == 1);
  CGAL_assertion(off3.z() == 0 || off3.z() == 1);
    
  int offx = ( ((off0.x() == 0 && off1.x() == 0 
		 && off2.x() == 0 && off3.x() == 0)
		|| (off0.x() == 1 && off1.x() == 1 
		    && off2.x() == 1 && off3.x() == 1)) ? 0 : 1);
  int offy = ( ((off0.y() == 0 && off1.y() == 0 
		 && off2.y() == 0 && off3.y() == 0)
		|| (off0.y() == 1 && off1.y() == 1 
		    && off2.y() == 1 && off3.y() == 1)) ? 0 : 1);
  int offz = ( ((off0.z() == 0 && off1.z() == 0 
		 && off2.z() == 0 && off3.z() == 0)
		|| (off0.z() == 1 && off1.z() == 1 
		    && off2.z() == 1 && off3.z() == 1)) ? 0 : 1);
    
  return( 4*offx + 2*offy + offz );
}

// construct a triangle from a given facet, given vertex offsets and a
// common offset
inline Triangle Scene::construct_triangle(const Cell_handle ch, int i,
    const Offset& off0, const Offset& off1, const Offset& off2, int off) const {
  if (it_type == P3DT::STORED || it_type == P3DT::UNIQUE) {
    CGAL_assertion( off == 0 );
    return p3dt.construct_triangle(
        ch->vertex((i+1)&3)->point(), ch->vertex((i+2)&3)->point(),
	ch->vertex((i+3)&3)->point(), off0, off1, off2);
  }
  Offset diff_off((off>>2)&1,(off>>1)&1,off&1);
  switch (it_type) {
  case P3DT::STORED_COVER_DOMAIN:
    return p3dt.construct_triangle(
        ch->vertex((i+1)&3)->point(), ch->vertex((i+2)&3)->point(),
	ch->vertex((i+3)&3)->point(),
        p3dt.combine_offsets(off0,-diff_off),
        p3dt.combine_offsets(off1,-diff_off),
        p3dt.combine_offsets(off2,-diff_off));
    break;
  case P3DT::UNIQUE_COVER_DOMAIN:
    return p3dt.construct_triangle(
        ch->vertex((i+1)&3)->point(), ch->vertex((i+2)&3)->point(),
	ch->vertex((i+3)&3)->point(),
        off0-diff_off, off1-diff_off, off2-diff_off);
    break;
  default:
    CGAL_assertion(false);
    return Triangle();
  }
}

// construct a triangle from a given cell, given vertex offsets and a
// common offset
inline Tetrahedron Scene::construct_tetrahedron(const Cell_handle ch,
    const Offset& off0, const Offset& off1, const Offset& off2,
    const Offset& off3, int off) const {
  if (it_type == P3DT::STORED || it_type == P3DT::UNIQUE) {
    CGAL_assertion( off == 0 );
    return p3dt.construct_tetrahedron(
        ch->vertex(0)->point(), ch->vertex(1)->point(),
	ch->vertex(2)->point(), ch->vertex(3)->point(), 
	off0, off1, off2, off3);
  }
  Offset diff_off((off>>2)&1,(off>>1)&1,off&1);
  switch (it_type) {
  case P3DT::STORED_COVER_DOMAIN:
    return p3dt.construct_tetrahedron(
	ch->vertex(0)->point(), ch->vertex(1)->point(),
	ch->vertex(2)->point(), ch->vertex(3)->point(), 
        p3dt.combine_offsets(off0,-diff_off),
        p3dt.combine_offsets(off1,-diff_off),
        p3dt.combine_offsets(off2,-diff_off),
	p3dt.combine_offsets(off3,-diff_off));
    break;
  case P3DT::UNIQUE_COVER_DOMAIN:
    return p3dt.construct_tetrahedron(
        ch->vertex(0)->point(), ch->vertex(1)->point(), 
	ch->vertex(2)->point(), ch->vertex(3)->point(),
        off0-diff_off, off1-diff_off, off2-diff_off, off3-diff_off);
    break;
  default:
    CGAL_assertion(false);
    return Tetrahedron();
  }
}

// Draw a point, either es vertex or as sphere
inline void Scene::gl_draw_vertex(const Point& p, const FT& radius) const {
  if (p.z()!=0 && in_plane) return;
  if (wireframe) {
    glVertex3f(p.x(),p.y(),p.z());
  } else {
    glPushMatrix();
    glTranslated(p.x(),p.y(),p.z());
    gluSphere(pQuadric,radius,15,15);
    glFlush();
    glPopMatrix();
  }
}

// Draw an edge, either as line or as cylinder
inline void Scene::gl_draw_edge(Point p1, Point p2, FT radius ) {
  if (in_plane && (p1.z()!=0. || p2.z()!=0.)) return;
  if (wireframe) {
    glVertex3d(p1.x(),p1.y(),p1.z());
    glVertex3d(p2.x(),p2.y(),p2.z());
    return;
  }

  Point p = p1;
  Vector v = p2-p1;

  FT len = (FT)std::sqrt(CGAL_NTS to_double(v*v));

  // normalize
  v = v / len;
  double angle = 0.0;
  if(std::sqrt(CGAL_NTS to_double(v.x()*v.x()+v.y()*v.y())) > 1)
    angle = 90.0f;
  else
    angle = asin(std::sqrt(CGAL_NTS to_double(v.x()*v.x()+v.y()*v.y())))/3.14159*180.0;
  if(v.z() < 0)
    angle = 180.0-angle;
  Vector axis;
  if (v.x() == 0.0 && v.y() == 0.0)
    axis = Vector(1.0,0.0,0.0);
  else
    axis = Vector(v.y(),-v.x(),0.0);

  glPushMatrix();        
  glTranslated(CGAL_NTS to_double(p.x()),
	       CGAL_NTS to_double(p.y()),
	       CGAL_NTS to_double(p.z()));
  glPushMatrix();        
  glRotated(-angle,CGAL_NTS to_double(axis.x()),
	    CGAL_NTS to_double(axis.y()),
	    CGAL_NTS to_double(axis.z()));

  // draw a cylinder
  gluCylinder(pQuadric,CGAL_NTS to_double(radius),
	      CGAL_NTS to_double(radius),
	      CGAL_NTS to_double(len),25,2);
  glPopMatrix();
  glPopMatrix();
  glPopMatrix();
}

// collect primitives (segments, triangles, tetrahedra) from the
// triangulation using the geometric iterators and store them in the
// given segment set
inline void Scene::primitives_from_geom_it(Segment_set& sset) {
  Point p0,p1,p2,p3;
  switch(draw_type) {
  case SEGMENT:
    for ( Segment_iterator sit = p3dt.periodic_segments_begin(it_type) ;
	  sit != p3dt.periodic_segments_end(it_type) ; ++sit ) {
      sset.insert(p3dt.segment(*sit));
    }
    break;
  case TRIANGLE:
    for ( Triangle_iterator tit = p3dt.periodic_triangles_begin(it_type) ;
	  tit != p3dt.periodic_triangles_end(it_type) ; ++tit ) {
      p0 = p3dt.point(tit->at(0));
      p1 = p3dt.point(tit->at(1));
      p2 = p3dt.point(tit->at(2));
      sset.insert(p0 < p1 ? Segment(p0,p1) : Segment(p1,p0));
      sset.insert(p0 < p2 ? Segment(p0,p2) : Segment(p2,p0));
      sset.insert(p1 < p2 ? Segment(p1,p2) : Segment(p2,p1));
    }
    break;
  case TETRAHEDRON:
    for ( Tetrahedron_iterator tit = p3dt.periodic_tetrahedra_begin(it_type) ;
	  tit != p3dt.periodic_tetrahedra_end(it_type) ; ++tit ) {
      p0 = p3dt.point(tit->at(0));
      p1 = p3dt.point(tit->at(1));
      p2 = p3dt.point(tit->at(2));
      p3 = p3dt.point(tit->at(3));
      sset.insert((p0 < p1) ? Segment(p0,p1) : Segment(p1,p0));
      sset.insert((p0 < p2) ? Segment(p0,p2) : Segment(p2,p0));
      sset.insert((p0 < p3) ? Segment(p0,p3) : Segment(p3,p0));
      sset.insert((p1 < p2) ? Segment(p1,p2) : Segment(p2,p1));
      sset.insert((p1 < p3) ? Segment(p1,p3) : Segment(p3,p1));
      sset.insert((p2 < p3) ? Segment(p2,p3) : Segment(p3,p2));
    }
    break;
  }
}

// clip segments from the given segment set that are partially outside
// of the unit cube/square. Eliminate those who are completely outside
inline void Scene::segment_clipping(Segment_set& sset) {
  Segment_clipper clipper;
  Segment_set sset_tmp;
  for (Segment_set::iterator it = sset.begin() ; it != sset.end() ; ++it) {
    Point s = it->source();
    Point t = it->target();
    if (clipper(s,t)) sset_tmp.insert((s<t?Segment(s,t):Segment(t,s)));
  }
  std::swap(sset, sset_tmp);
}

// clip segments from the given segment set that are partially outside
// of the unit cube/square. Draw their outside part in a different
// color.
// TODO: don't eliminate segments that are completely outside but draw
// them in the different color as well
inline void Scene::segment_2color_clipping (Segment_set& sset) {
  Segment_clipper clipper;
  Segment_set sset_tmp, sset_out;
  for (Segment_set::iterator it = sset.begin() ; it != sset.end() ; ++it) {
    Point s = it->source();
    Point t = it->target();
    if (clipper(s,t)) {
      sset_tmp.insert((s<t?Segment(s,t):Segment(t,s)));
      Point p = it->source();
      Point q = it->target();
      if (Segment(p,s).squared_length() > Segment(p,t).squared_length())
	std::swap(s,t);
      if (p!=s) sset_out.insert((p<s?Segment(p,s):Segment(s,p)));
      if (q!=t) sset_out.insert((q<t?Segment(q,t):Segment(t,q)));
    }
  }

  for (Segment_set::iterator sit = sset_out.begin() ;
       sit != sset_out.end() ; ++sit)
    gl_draw_edge(sit->source(), sit->target(), 0.005);

  std::swap(sset, sset_tmp);
}

// Draw the faces of the tetrahedron in which the moving point is currently
// located transparently. It depends on it_type which periodic copies
// of the respective cell will be drawn. In general it will be all
// cells that occur in the draw list
void Scene::gl_draw_location() {
  if (p3dt.number_of_vertices() == 0) return;
  // Do the point location
  Cell_handle ch = p3dt.locate(moving_point);
  std::vector<Projected_triangle> cf;

  // Transparency
  glEnable(GL_BLEND);
  glColor4f(0.f, 0.f , 0.5f, 0.5f);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  if (in_plane) {
    int i=0;
    int count = 0;
    // Figure out whether there is a facet that is completly contained
    // in the z=0 plane
    for (int j=0 ; j<4 ; j++) {
      if (ch->vertex(j)->point().z() != 0.0 ||
	  p3dt.get_offset(ch,j).z() != 0) {
	i=j;
	count++;
      }
    }
    // If so, compute its triangle(s) and insert it in cf
    if (count==1) {
      Offset off0, off1, off2;
      get_tri_offsets(ch, i, off0, off1, off2);
      int diff_off = get_tri_drawing_offsets(ch, i);
      for (int offs=0 ; offs<=diff_off ; offs++) {
	if ((((~offs)|diff_off)&7)!=7) continue;
	Triangle tri_to_draw = construct_triangle(ch,i,off0,off1,off2,offs);
	Point p = tri_to_draw.vertex(0);
	Point q = tri_to_draw.vertex(1);
	Point r = tri_to_draw.vertex(2);
	cf.push_back(Projected_triangle(.0,Triangle(p,q,r)));
      }
    }
  } else {
    double modelMatrix[16];
    double projMatrix[16];
    int viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    // Compute the triangles that are the facets of the cell and
    // insert them in cf
    Offset off0,off1,off2,off3;
    get_tet_offsets(ch, off0, off1, off2, off3);
    int diff_off = get_tet_drawing_offsets(ch);

    for (int offs=0 ; offs<=diff_off ; offs++) {
      if ((((~offs)|diff_off)&7)!=7) continue;
      Tetrahedron tet_to_draw = construct_tetrahedron(
	  ch, off0, off1, off2, off3, offs);

      for(int i=0; i < 4; i++){ 
	Point p = tet_to_draw.vertex((i+1)&3);
	Point q = tet_to_draw.vertex((i+2)&3);
	Point r = tet_to_draw.vertex((i+3)&3);
	Vector c= (Vector(Point(),p)+Vector(Point(),q)+Vector(Point(),r))/3.;
	Point cp = Point(c.x(),c.y(),c.z());
	// project facet center 
	double px,py,pz;
	gluProject(cp.x(),cp.y(),cp.z(),
		   modelMatrix, projMatrix, viewport,
		   &px,&py,&pz);
	cf.push_back(Projected_triangle(pz,Triangle(p,q,r)));
      }
    }

    // Sort cf according to their z coordinates to enable transparency
    std::sort(cf.begin(), cf.end(), Projected_triangle::closer);
  }

  // Draw all triangles from cf
  glBegin(GL_TRIANGLES);
  for (std::vector<Projected_triangle >::iterator cfit = cf.begin() ; cfit != cf.end() ; cfit++) {
    Point p = cfit->t().vertex(0);
    Point q = cfit->t().vertex(1);
    Point r = cfit->t().vertex(2);
    glVertex3d(p.x(), p.y(), p.z());
    glVertex3d(q.x(), q.y(), q.z());
    glVertex3d(r.x(), r.y(), r.z());
  }
  glEnd();

  glDisable(GL_BLEND);
}

// Draw the boundary faces of the current conflict region of the
// moving point transparently. It depends on it_type which periodic
// copies of the respective cell will be drawn. In general it will be
// all cells that occur in the draw list.
void Scene::gl_draw_conflict() {
  if (p3dt.number_of_vertices() == 0) return;
  Cell_handle ch;
  std::vector<Cell_handle> cic;
  std::vector<Facet> boundary_facets;
  // Find the conflict region
  Cell_handle c = p3dt.locate(moving_point);
  p3dt.find_conflicts(moving_point,c,std::back_inserter(boundary_facets),std::back_inserter(cic),CGAL::Emptyset_iterator());

  std::vector<Projected_triangle> bfm;

  // Transparency
  glEnable(GL_BLEND);
  glColor4f(.69f, 0.18f , 0.26f, 0.6f);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  if (in_plane) {
    for (unsigned int k=0 ; k<cic.size(); k++) {
      ch = cic[k];

      int i = 0;
      int count = 0;
      // Figure out whether there is a facet that is completely
      // contained in the z=0 plane
      for (int j=0 ; j<4 ; j++) {
	if (ch->vertex(j)->point().z() != 0.0 ||
	    p3dt.get_offset(ch,j).z() != 0) {
	  i=j;
	  count++;
	}
      }
      // If so, compute its triangle(s) and insert it in bfm
      if (count==1) {
	Offset off0, off1, off2;
	get_tri_offsets(ch, i, off0, off1, off2);
	int diff_off = get_tri_drawing_offsets(ch,i);
	for (int offs = 0 ; offs<=diff_off ; offs++) {
	  if ((((~offs)|diff_off)&7)!=7) continue;
	  Triangle tri_to_draw = construct_triangle(ch,i,off0,off1,off2,offs);
	  Point p = tri_to_draw.vertex(0);
	  Point q = tri_to_draw.vertex(1);
	  Point r = tri_to_draw.vertex(2);
	  bfm.push_back(Projected_triangle(.0,Triangle(p,q,r)));
	}
      }
    }
  } else {
    double modelMatrix[16];
    double projMatrix[16];
    int viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    for (unsigned int i=0 ; i<boundary_facets.size(); i++) {
      ch = boundary_facets[i].first;
      int j=boundary_facets[i].second;

      // Compute the triangle(s) of the facet and insert them in bfm
      Offset off0, off1, off2;
      get_tri_offsets(ch, j, off0, off1, off2);
      int diff_off = get_tri_drawing_offsets(ch, j);
      for (int offs=0 ; offs<=diff_off ; offs++) {
	if ((((~offs)|diff_off)&7)!=7) continue;
	Triangle tri_to_draw = construct_triangle(ch,j,off0,off1,off2,offs);
	Point p = tri_to_draw.vertex(0);
	Point q = tri_to_draw.vertex(1);
	Point r = tri_to_draw.vertex(2);
	Vector c= (Vector(Point(),p)+Vector(Point(),q)+Vector(Point(),r))/3.;
	Point cp = Point(c.x(),c.y(),c.z());
	// project facet center 
	double px,py,pz;
	gluProject(cp.x(),cp.y(),cp.z(),
		   modelMatrix, projMatrix, viewport,
		   &px,&py,&pz);
	bfm.push_back(Projected_triangle(pz,Triangle(p,q,r)));
      }
    }
    
    // Sort bfm according to their z coordinates to enable transparency
    std::sort(bfm.begin(), bfm.end(), Projected_triangle::closer);
  }

  // Draw all triangles from bfm
  glBegin(GL_TRIANGLES);
  for (std::vector<Projected_triangle >::iterator bfmit = bfm.begin() ; 
       bfmit != bfm.end() ; bfmit++) {
    Point p = bfmit->t().vertex(0);
    Point q = bfmit->t().vertex(1);
    Point r = bfmit->t().vertex(2);
    glVertex3d(p.x(), p.y(), p.z());
    glVertex3d(q.x(), q.y(), q.z());
    glVertex3d(r.x(), r.y(), r.z());
  }
  glEnd();
  
  glDisable(GL_BLEND);
}

// provide some color constants for the GLU primitive appearance.
void Scene::change_material(const QString &string) {
    float    ambient[]  = {0.0f,0.0f,0.0f,1.0f};
    float    diffuse[]  = {0.0f,0.0f,0.0f,1.0f};
    float    specular[]  = {0.0f,0.0f,0.0f,1.0f};
    float    emission[]  = {0.3f,0.3f,0.3f,1.0f};
    float shininess[] = {0.0f};

    // Change
    if(string == "Silver")
    {
        // Ambient
        ambient[0] = 0.19225f;
        ambient[1] = 0.19225f;
        ambient[2] = 0.19225f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.50754f;
        diffuse[1] = 0.50754f;
        diffuse[2] = 0.50754f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.508273f;
        specular[1] = 0.508273f;
        specular[2] = 0.508273f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 51.2f;
    }

    else

    // Change
    if(string == "Light silver")
    {
        // Ambient
        ambient[0] = 0.49225f;
        ambient[1] = 0.49225f;
        ambient[2] = 0.49225f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.80754f;
        diffuse[1] = 0.80754f;
        diffuse[2] = 0.80754f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.808273f;
        specular[1] = 0.808273f;
        specular[2] = 0.808273f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 51.2f;
    }

    else

    // Change
    if(string == "Gold")
    {
        // Ambient
        ambient[0] = 0.24725f;
        ambient[1] = 0.1995f;
        ambient[2] = 0.0745f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.75164f;
        diffuse[1] = 0.60648f;
        diffuse[2] = 0.22648f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.928281f;
        specular[1] = 0.855802f;
        specular[2] = 0.666065f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 51.2f;
    }

    else

    // Change
    if(string == "Jade")
    {
        // Ambient
        ambient[0] = 0.135f;
        ambient[1] = 0.2225f;
        ambient[2] = 0.1575f;
        ambient[3] = 0.95f;
        // Diffuse
        diffuse[0] = 0.54f;
        diffuse[1] = 0.89f;
        diffuse[2] = 0.63f;
        diffuse[3] = 0.95f;
        // Specular
        specular[0] = 0.316228f;
        specular[1] = 0.316228f;
        specular[2] = 0.316228f;
        specular[3] = 0.95f;
        // Shininess
        shininess[0] = 12.8f;
    }

    else

    // Change
    if(string == "Light blue")
    {
        // Ambient
        ambient[0] = 0.0f;
        ambient[1] = 0.5f;
        ambient[2] = 0.75f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.0f;
        diffuse[1] = 0.5f;
        diffuse[2] = 1.0f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.75f;
        specular[1] = 0.75f;
        specular[2] = 0.75f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 64.0f;
    }

    else

    // Change
    if(string == "Pink")
    {
        // Ambient
        ambient[0] = 0.75f;
        ambient[1] = 0.0f;
        ambient[2] = 0.5f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 1.0f;
        diffuse[1] = 0.0f;
        diffuse[2] = 0.5f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.75f;
        specular[1] = 0.75f;
        specular[2] = 0.75f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 64.0f;
    }

    else

    // Change
    if(string == "Red")
    {
        // Ambient
        ambient[0] = 0.75f;
        ambient[1] = 0.0f;
        ambient[2] = 0.0f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 1.0f;
        diffuse[1] = 0.0f;
        diffuse[2] = 0.0f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.75f;
        specular[1] = 0.75f;
        specular[2] = 0.75f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 64.0f;
    }

    else

    // Change
    if(string == "Emerald")
    {
        // Ambient
        ambient[0] = 0.0215f;
        ambient[1] = 0.1745f;
        ambient[2] = 0.0215f;
        ambient[3] = 0.55f;
        // Diffuse
        diffuse[0] = 0.07568f;
        diffuse[1] = 0.61424f;
        diffuse[2] = 0.07568f;
        diffuse[3] = 0.55f;
        // Specular
        specular[0] = 0.633f;
        specular[1] = 0.727811f;
        specular[2] = 0.633f;
        specular[3] = 0.55f;
        // Shininess
        shininess[0] = 76.8f;
    }

    else

    // Change
    if(string == "Polished silver")
    {
        // Ambient
        ambient[0] = 0.23125f;
        ambient[1] = 0.23125f;
        ambient[2] = 0.23125f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.2775f;
        diffuse[1] = 0.2775f;
        diffuse[2] = 0.2775f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.773911f;
        specular[1] = 0.773911f;
        specular[2] = 0.773911f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 89.6f;
    }

    else

    // Change
    if(string == "Chrome")
    {
        // Ambient
        ambient[0] = 0.25f;
        ambient[1] = 0.25f;
        ambient[2] = 0.25f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.4f;
        diffuse[1] = 0.4f;
        diffuse[2] = 0.4f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.774597f;
        specular[1] = 0.774597f;
        specular[2] = 0.774597f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 76.8f;
    }

    else

    // Change
    if(string == "Copper")
    {
        // Ambient
        ambient[0] = 0.19125f;
        ambient[1] = 0.0735f;
        ambient[2] = 0.0225f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.7038f;
        diffuse[1] = 0.27048f;
        diffuse[2] = 0.0828f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.256777f;
        specular[1] = 0.137622f;
        specular[2] = 0.086014f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 12.8f;
    }

    else

    // Change
    if(string == "Green")
    {
        // Ambient
        ambient[0] = 0.0225f;
        ambient[1] = 0.19125f;
        ambient[2] = 0.0735f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.0828f;
        diffuse[1] = 0.3038f;
        diffuse[2] = 0.14048f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.086014f;
        specular[1] = 0.306777f;
        specular[2] = 0.117622f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 12.8f;
    }

    else

    // Change
    if(string == "Blue")
    {
        // Ambient
        ambient[0] = 0.0225f;
        ambient[1] = 0.0735f;
        ambient[2] = 0.19125f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.0828f;
        diffuse[1] = 0.14048f;
        diffuse[2] = 0.6038f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.086014f;
        specular[1] = 0.117622f;
        specular[2] = 0.306777f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 12.8f;
    }

    else

    // Change
    if(string == "Polished gold")
    {
        // Ambient
        ambient[0] = 0.24725f;
        ambient[1] = 0.2245f;
        ambient[2] = 0.0645f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.34615f;
        diffuse[1] = 0.3143f;
        diffuse[2] = 0.0903f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.797357f;
        specular[1] = 0.723991f;
        specular[2] = 0.208006f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 83.2f;
    }

    else

    // Change
    if(string == "Pewter")
    {
        // Ambient
        ambient[0] = 0.105882f;
        ambient[1] = 0.058824f;
        ambient[2] = 0.113725f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.427451f;
        diffuse[1] = 0.470588f;
        diffuse[2] = 0.541176f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.333333f;
        specular[1] = 0.333333f;
        specular[2] = 0.521569f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 9.84615f;
    }

    else

    // Change
    if(string == "Obsidian")
    {
        // Ambient
        ambient[0] = 0.05375f;
        ambient[1] = 0.05f;
        ambient[2] = 0.06625f;
        ambient[3] = 0.82f;
        // Diffuse
        diffuse[0] = 0.18275f;
        diffuse[1] = 0.17f;
        diffuse[2] = 0.22525f;
        diffuse[3] = 0.82f;
        // Specular
        specular[0] = 0.332741f;
        specular[1] = 0.328634f;
        specular[2] = 0.346435f;
        specular[3] = 0.82f;
        // Shininess
        shininess[0] = 38.4f;
    }

    else

    // Change
    if(string == "Black plastic")
    {
        // Ambient
        ambient[0] = 0.0f;
        ambient[1] = 0.0f;
        ambient[2] = 0.0f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.01f;
        diffuse[1] = 0.01f;
        diffuse[2] = 0.01f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.5f;
        specular[1] = 0.5f;
        specular[2] = 0.5f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 32.0f;
    }

    else

    // Change
    if(string == "Polished bronze")
    {
        // Ambient
        ambient[0] = 0.25f;
        ambient[1] = 0.148f;
        ambient[2] = 0.006475f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.4f;
        diffuse[1] = 0.2368f;
        diffuse[2] = 0.1036f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.774597f;
        specular[1] = 0.458561f;
        specular[2] = 0.200621f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 76.8f;
    }

    
    else

    // Change
    if(string == "Polished copper")
    {
        // Ambient
        ambient[0] = 0.2295f;
        ambient[1] = 0.08825f;
        ambient[2] = 0.0275f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.5508f;
        diffuse[1] = 0.2118f;
        diffuse[2] = 0.066f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.580594f;
        specular[1] = 0.223257f;
        specular[2] = 0.0695701f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 51.2f;
    }

    else

    // Change
    if(string == "Pearl")
    {
        // Ambient
        ambient[0] = 0.25f;
        ambient[1] = 0.20725f;
        ambient[2] = 0.20725f;
        ambient[3] = 0.922f;
        // Diffuse
        diffuse[0] = 1.0f;
        diffuse[1] = 0.829f;
        diffuse[2] = 0.829f;
        diffuse[3] = 0.922f;
        // Specular
        specular[0] = 0.296648f;
        specular[1] = 0.296648f;
        specular[2] = 0.296648f;
        specular[3] = 0.922f;
        // Shininess
        shininess[0] = 11.264f;
    }

    else

    // Change
    if(string == "Ruby")
    {
        // Ambient
        ambient[0] = 0.1745f;
        ambient[1] = 0.01175f;
        ambient[2] = 0.01175f;
        ambient[3] = 0.55f;
        // Diffuse
        diffuse[0] = 0.61424f;
        diffuse[1] = 0.04136f;
        diffuse[2] = 0.04136f;
        diffuse[3] = 0.55f;
        // Specular
        specular[0] = 0.727811f;
        specular[1] = 0.626959f;
        specular[2] = 0.626959f;
        specular[3] = 0.55f;
        // Shininess
        shininess[0] = 76.8f;
    }

    else

    // Change
    if(string == "Turquoise")
    {
        // Ambient
        ambient[0] = 0.1f;
        ambient[1] = 0.18725f;
        ambient[2] = 0.1745f;
        ambient[3] = 0.8f;
        // Diffuse
        diffuse[0] = 0.396f;
        diffuse[1] = 0.74151f;
        diffuse[2] = 0.69102f;
        diffuse[3] = 0.8f;
        // Specular
        specular[0] = 0.297254f;
        specular[1] = 0.30829f;
        specular[2] = 0.306678f;
        specular[3] = 0.8f;
        // Shininess
        shininess[0] = 12.8f;
    }

    else

    // Change
    if(string == "Brass")
    {
        // Ambient
        ambient[0] = 0.329412f;
        ambient[1] = 0.223529f;
        ambient[2] = 0.027451f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.780392f;
        diffuse[1] = 0.268627f;
        diffuse[2] = 0.113725f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.992157f;
        specular[1] = 0.741176f;
        specular[2] = 0.807843f;
        specular[3] = 1.0f;
        // Shininess
        shininess[0] = 27.8974f;
    }

    // apply
    if (wireframe) {
      glColor4f(diffuse[0],diffuse[1],diffuse[2],diffuse[3]);
    } else {
      glMaterialfv( GL_FRONT, GL_AMBIENT,   ambient);
      glMaterialfv( GL_FRONT, GL_DIFFUSE,   diffuse);
      glMaterialfv( GL_FRONT, GL_SPECULAR,  specular);
      glMaterialfv( GL_FRONT, GL_SHININESS, shininess);
      glMaterialfv( GL_FRONT, GL_EMISSION,  emission); 
    }
}
