// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 *
 * a viewer of spherical arrangements
 */

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Spherical_map.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Sphere_point_location.h>

#include "ArcLines.h"

#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;

typedef CGAL::Sphere_traits<Kernel>                     Traits_2;
typedef Kernel::Direction_3                             Direction_3;
typedef CGAL::SphereTopologicalMap<Kernel>              SphereType;
typedef CGAL::Spherical_map<SphereType, Traits_2>       SphereMap;
typedef SphereMap::X_monotone_curve_2                   X_monotone_curve_2;
typedef CGAL::Sphere_walk_along_line_point_location<SphereMap>
                                                        PointLocation;
typedef SphereMap::Halfedge_iterator                    Halfedge_iterator;
typedef SphereMap::Halfedge_const_iterator              Halfedge_const_iterator;
typedef SphereMap::Halfedge_handle                      Halfedge_handle;
typedef SphereMap::Face_handle                          Face_handle;
typedef SphereMap::Vertex_handle                        Vertex_handle;
typedef SphereMap::CGM                                  CGM;
typedef CGM::Arrangement                                CGM_arrangement;
typedef CGM_arrangement::Halfedge_iterator              CGM_halfedge_iterator;
typedef CGM_arrangement::Halfedge_const_iterator
  CGM_halfedge_const_iterator;
typedef CGM::Projected_normal                           Projected_normal;

// pi value
#define PI 3.1415926535897

// the spherical map arrangement
static SphereMap s_spMap; // the spherical map
static GLint s_sphereDisplayList; // display list for sphere map
static GLint s_axisCircsDisplayList; // display list for circles on axes
static GLint s_cubeDisplayList; // display list for cubical map
static GLint s_pointLocateDisplayList; // display list for PointLocation
static GLint s_pointLocateCubeDisplayList; // display list for PointLocation in cube
// the angular differnce between lines approximating a spherical arc
static float s_arcDelta = 0.1;
static float s_angX = 0,
  s_angY = 0; // the rotation angles around the x and y axis
static bool s_displayAxis=true; // should the axis circles be displayed
static bool s_displaySphere=true; // should display sphere or cube
static bool s_showPointLocation=false; // show point located object
static PointLocation s_pointLocate; // the point location object
float s_scale = 0.8; // the scaling of the cube size

/*
 read a spherical arrangement from file,
 file format:
   first line - number of curves
   other lines - each line represents a curve using 6
     numbers (in rational format), first three numbers states curve's start direction,
   last 3 numbers state's curve's end direction
*/
void readSpericalMap(const char *fileName, SphereMap &spMap) {
  std::ifstream inStream;
  inStream.open(fileName);
  unsigned int numCurves; // will hold number of curves in data file
  inStream >> numCurves;
  unsigned int i;
  for (i=0;i<numCurves;++i) {
    // get one curve
    Direction_3 stDir, enDir;
    inStream >> stDir;
    inStream >> enDir;
    try {
      spMap.insert(X_monotone_curve_2(stDir, enDir));
    } catch (...) {
      // probably an overlapping segments exception uppon insert
      std::cerr << "there was an internal error when inserting the curve" <<
        std::endl << "from (" << stDir << ") to (" << enDir <<
        ") to the cubical gaussian map" << std::endl <<
        "this may have been caused due to intersecting curves" << std::endl;
    }
  }
  inStream.close();
}

/*
 write scene to screen window
*/
void renderScene(void) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (s_displaySphere) {
    // dislpay the sphere
    glCallList(s_sphereDisplayList);
    if (s_showPointLocation) {
      glCallList(s_pointLocateDisplayList);
    }
    if (s_displayAxis) {
      glCallList(s_axisCircsDisplayList);
    }
  } else {
    // display the cube
    glCallList(s_cubeDisplayList);
    if (s_showPointLocation) {
      glCallList(s_pointLocateCubeDisplayList);
    }
  }
  glutSwapBuffers();
}

/*
  draw a spherical arc of a halfedge
*/
void drawArc(const Halfedge_const_iterator &he) {
  float st[3], en[3];
  Direction_3 dirSt, dirEn;
  // translate to float arrays
  dirSt = (he->source())->direction();
  dirEn = (he->target())->direction();
  st[0] = CGAL::to_double(dirSt.dx());
  st[1] = CGAL::to_double(dirSt.dy());
  st[2] = CGAL::to_double(dirSt.dz());
  normalize(st);
  en[0] = CGAL::to_double(dirEn.dx());
  en[1] = CGAL::to_double(dirEn.dy());
  en[2] = CGAL::to_double(dirEn.dz());
  normalize(en);
  // create arc segmentation object for the arc
  ArcLines arcLines(st, en);
  float curAng = 0.f;
  glBegin(GL_LINE_STRIP);
  float curPnt[3];
  // draw lines
  while (curAng < arcLines.getAngle()) {
    arcLines.getPoint(curAng, curPnt);
    curAng += s_arcDelta;
    glVertex3fv(curPnt);
  }
  glVertex3fv(en);
  glEnd();
}

/*
  draw a spherical arc of a halfedge projection on the cubical map
*/
void drawCubeArc(const Halfedge_const_iterator &he) {
  const CGM *cgm = s_spMap.get_cgm();
  CGM::Arr_halfedge_const_iterator hit;
  hit = he->cubeHalfedge();
  // find cubical twin, start at twin of cube halfedge which is twin of last halfedge
  if (((hit->twin())->face())->is_unbounded()) {
    hit = cgm->get_adjacent_halfedge_handle(hit);
  } else {
    hit = hit->twin();
  }

  // walk from last edge twin until a real vertex and print
  unsigned int i = (hit->source())->get_face_id(); // get cubical map id
  unsigned int id1, id2, id3;
  int mapConst;
  // draw last halfedge
  float src[3], dst[3];
  if (i<3) {
    id1 = i % 3;
    id2 = (i+2) % 3;
    id3 = (i+1) % 3;
    mapConst = -1;
  } else {
    id1 = i % 3;
    id2 = (i+1) % 3;
    id3 = (i+2) % 3;
    mapConst = 1;
  }
  src[id1]=mapConst*s_scale; dst[id1]=mapConst*s_scale;
  src[id2] = CGAL::to_double(((hit->source())->point()).x())*s_scale;
  src[id3] = CGAL::to_double(((hit->source())->point()).y())*s_scale;
  dst[id2] = CGAL::to_double(((hit->target())->point()).x())*s_scale;
  dst[id3] = CGAL::to_double(((hit->target())->point()).y())*s_scale;
  glBegin(GL_LINES);
  glVertex3fv(src);
  glVertex3fv(dst);
  glEnd();
  // draw rest of the halfedges for the spherical arc
  while (!((hit->target())->getReal())) {
    hit = hit->next();
    while (!hit->get_is_real()) {
      hit = cgm->get_adjacent_halfedge_handle(hit);
      hit = hit->next();
    }
    unsigned int i = (hit->source())->get_face_id(); // get cubical map id
    unsigned int id1, id2, id3;
    int mapConst;
    float src[3], dst[3];
    if (i<3) {
      id1 = i % 3;
      id2 = (i+2) % 3;
      id3 = (i+1) % 3;
      mapConst = -1;
    } else {
      id1 = i % 3;
      id2 = (i+1) % 3;
      id3 = (i+2) % 3;
      mapConst = 1;
    }
    src[id1]=mapConst*s_scale; dst[id1]=mapConst*s_scale;
    src[id2] = CGAL::to_double(((hit->source())->point()).x())*s_scale;
    src[id3] = CGAL::to_double(((hit->source())->point()).y())*s_scale;
    dst[id2] = CGAL::to_double(((hit->target())->point()).x())*s_scale;
    dst[id3] = CGAL::to_double(((hit->target())->point()).y())*s_scale;
    glBegin(GL_LINES);
    glVertex3fv(src);
    glVertex3fv(dst);
    glEnd();
  }
}

/*
 draw the spherical map
*/
void drawSphereMap(void) {

  Halfedge_iterator hit;
  glColor3f(1.f,1.f,1.f);
  // draw all all halfedges on the sphere
  for (hit = s_spMap.halfedges_begin(); hit != s_spMap.halfedges_end(); ++hit) {
    drawArc(hit);
  }
}

/*
 draw the cubical map
*/
void drawCubicalMap(void) {
  unsigned int i; // map number
  CGM *theCGM = s_spMap.get_cgm();
  unsigned int id1, id2, id3;
  int mapConst; // the map constant coordinate value
  float src[3], dst[3];
  glBegin(GL_LINES);
  for (i=0;i<6;++i) { // for all cubical maps
    CGM_arrangement &curMap = theCGM->get_arrangement(i);
    if (i<3) {
      id1 = i % 3;
      id2 = (i+2) % 3;
      id3 = (i+1) % 3;
      mapConst = -1;
    } else {
      id1 = i % 3;
      id2 = (i+1) % 3;
      id3 = (i+2) % 3;
      mapConst = 1;
    }
    src[id1]=mapConst*s_scale; dst[id1]=mapConst*s_scale;
    CGM_halfedge_iterator hit;
    // iterate all map halfedges
    for (hit = curMap.halfedges_begin(); hit != curMap.halfedges_end(); ++hit) {
      // additional vertex info, mark target vertex face id
      (hit->target())->set_face_id(i);
      // choose color according to halfedge "real" state
      if (hit->get_is_real()) {
        glColor3f(1.f, 1.f, 1.f);
      } else {
        glColor3f(0.f, 0.f, 1.f);
      }
      // draw current segment
      src[id2] = CGAL::to_double(((hit->source())->point()).x())*s_scale;
      src[id3] = CGAL::to_double(((hit->source())->point()).y())*s_scale;
      dst[id2] = CGAL::to_double(((hit->target())->point()).x())*s_scale;
      dst[id3] = CGAL::to_double(((hit->target())->point()).y())*s_scale;
      glVertex3fv(src);
      glVertex3fv(dst);
    }
  }
  glEnd();
}

/*
 draw the axis parallel circles
*/
void drawAxisCircles(void) {
  glColor4f(0.f,0.f,0.9f,0.9f);
  // draw axis circles
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  // the axis parallel circles are transparent
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  GLint circle_points = 40;
  float angle;

  int i;
  glBegin(GL_LINE_LOOP);  // plane xy
    for (i = 0; i < circle_points; i++) {
    angle = 2*PI*i/circle_points;
    glVertex3f(0.98*cos(angle), 0.98*sin(angle), 0.f);
  }
  glEnd();
  glBegin(GL_LINE_LOOP); // plane yz
    for (i = 0; i < circle_points; i++) {
    angle = 2*PI*i/circle_points;
    glVertex3f(0.f, 0.98*cos(angle), 0.98*sin(angle));
  }
  glEnd();
  glBegin(GL_LINE_LOOP); // plane xz
    for (i = 0; i < circle_points; i++) {
    angle = 2*PI*i/circle_points;
    glVertex3f(0.98*cos(angle), 0.f, 0.98*sin(angle));
  }
  glEnd();
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
}

/*
 create display list to draw spherical map
*/
GLuint createMapDL() {
  GLuint sphereMapDL;

  // Create the id for the list
  sphereMapDL = glGenLists(1);

  // start list
  glNewList(sphereMapDL, GL_COMPILE);

  // call the function that contains
  // the rendering commands
  drawSphereMap();
  // endList
  glEndList();
  return(sphereMapDL);
}

/*
 create display list to draw the cubical map
*/
GLuint createCubeMapDL() {
  GLuint cubeMapDL;
  cubeMapDL = glGenLists(1);
  glNewList(cubeMapDL, GL_COMPILE);
  drawCubicalMap();
  glEndList();
  return(cubeMapDL);
}

/*
 create a display list for circles on the axes planes
*/
GLuint createAxisCircsDL() {
  GLuint circlesDL;
  circlesDL = glGenLists(1);
  glNewList(circlesDL, GL_COMPILE);
  drawAxisCircles();
  glEndList();
  return circlesDL;
}

/*
 initialize display data
*/
void initScene() {
  s_pointLocate.init(&s_spMap);
  s_pointLocateDisplayList = glGenLists(1);
  s_pointLocateCubeDisplayList = glGenLists(1);
  glEnable(GL_DEPTH_TEST);
  s_sphereDisplayList = createMapDL();
  s_axisCircsDisplayList = createAxisCircsDL();
  s_cubeDisplayList = createCubeMapDL();
}

/*
 action uppon window resize
*/
void changeSize(int w, int h) {
  if(h == 0)
    h = 1;
  float ratio = 1.0* w / h;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport (0, 0, w, h);
  gluPerspective(45,ratio,1,1000);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0,0.0,4.0f,
      0.0,0.0,0.0,
    0.0f,1.0f,0.0f);
  glRotatef(s_angY, 0.f,1.f,0.f);
  glRotatef(s_angX, 1.f,0.f,0.f);

}

/*
 special input keys actions
*/
void inputKey(int key, int x, int y) {
  switch (key) {
    case GLUT_KEY_LEFT :
      s_angY -= 1.f;
      break;
    case GLUT_KEY_RIGHT :
      s_angY +=1.f;
      break;
    case GLUT_KEY_UP :
      s_angX -=1.f;
      break;
    case GLUT_KEY_DOWN :
      s_angX += 1.f;
      break;
  }
  glLoadIdentity();
  gluLookAt(0.0,0.0,4.0f,
      0.0,0.0,0.0,
    0.0f,1.0f,0.0f);
  glRotatef(s_angY, 0.f,1.f,0.f);
  glRotatef(s_angX, 1.f,0.f,0.f);
  glutPostRedisplay();
}

/*
 process point location, get direction and show it's object
*/
void pointLocate(void) {
  // get direction to locate
  Direction_3 locateDir;
  std::cout << "enter point location direction: ";
  std::cin >> locateDir;
  CGAL::Object location;
  location = s_pointLocate.locate(locateDir);
  // the possible object to locate
  Face_handle foundFace;
  Vertex_handle foundVertex;
  Halfedge_handle foundHalfedge;

  // create point location display list for sphere
  glDeleteLists(s_pointLocateDisplayList, 1);
  s_pointLocateDisplayList = glGenLists(1);
  glLoadIdentity();
  // start list
  glNewList(s_pointLocateDisplayList, GL_COMPILE);

  // diaplay location point on sphere
  float pnt[3];
  pnt[0] = CGAL::to_double(locateDir.dx());
  pnt[1] = CGAL::to_double(locateDir.dy());
  pnt[2] = CGAL::to_double(locateDir.dz());
  normalize(pnt);
  glColor3f(0.f, 1.f, 0.f);
  glPushMatrix();
  GLUquadricObj *qobj=gluNewQuadric();
  glTranslatef(pnt[0],pnt[1],pnt[2]);
  gluSphere(qobj, 0.019, 4, 4);
  gluDeleteQuadric(qobj);
  glPopMatrix();


  glColor3f(1.f,0.f,0.f);
  glLineWidth(3);

  if (CGAL::assign(foundFace, location)) {
    // if found a face, show face
    std::cout << "found a face" << std::endl;
    if (foundFace->does_outer_ccb_exist()) {
      // show outer ccb
      SphereMap::Ccb_halfedge_circulator heif;
      heif = foundFace->outer_ccb();
      SphereMap::Ccb_halfedge_circulator circEn = heif;
      do {
        std::cout << heif->curve() << std::endl;
        drawArc(heif);
        ++heif;
      } while (heif != circEn);
    } // end of outer ccb display
    // display holes
    SphereMap::Holes_iterator holesIt;
    for (holesIt = foundFace->holes_begin(); holesIt != foundFace->holes_end();
        ++holesIt) {
      std::cout << "a hole: " << std::endl;
      SphereMap::Ccb_halfedge_circulator heif;
      heif = *holesIt;
      SphereMap::Ccb_halfedge_circulator circEn=heif;
      do {
        std::cout << "edge: (" << heif->source()->direction() << ") - (" <<
          heif->target()->direction() << ")" << std::endl;
        // draw halfedge
        drawArc(heif);
        ++heif;
      } while (heif !=circEn);
    }
  }
  if (CGAL::assign(foundHalfedge, location)) {
    // located an arc
    drawArc(foundHalfedge);
  }
  if (CGAL::assign(foundVertex, location)) {
    // located a vertex
    float pnt[3];
    pnt[0] = CGAL::to_double((foundVertex->direction()).dx());
    pnt[1] = CGAL::to_double((foundVertex->direction()).dy());
    pnt[2] = CGAL::to_double((foundVertex->direction()).dz());
    normalize(pnt);
    //glPointSize(7);
    glPushMatrix();
    GLUquadricObj *qobj=gluNewQuadric();
    glTranslatef(pnt[0],pnt[1],pnt[2]);
    gluSphere(qobj, 0.02, 4, 4);
    gluDeleteQuadric(qobj);
    glPopMatrix();
  }
  glLineWidth(1);
  glEndList();

  // update display in cube mode
  glDeleteLists(s_pointLocateCubeDisplayList, 1);
  s_pointLocateCubeDisplayList = glGenLists(1);
  glNewList(s_pointLocateCubeDisplayList, GL_COMPILE);
  // start list
  // diaplay location point
  glColor3f(0.f, 1.f, 0.f);
  Projected_normal pn;
  pn.compute_projection(locateDir.vector());
  pnt[0]=CGAL::to_double((pn.get_projected_normal()).x())*s_scale;
  pnt[1]=CGAL::to_double((pn.get_projected_normal()).y())*s_scale;
  pnt[2]=CGAL::to_double((pn.get_projected_normal()).z())*s_scale;
  glPushMatrix();
  qobj=gluNewQuadric();
  glTranslatef(pnt[0],pnt[1],pnt[2]);
  gluSphere(qobj, 0.019, 4, 4);
  gluDeleteQuadric(qobj);
  glPopMatrix();

  glColor3f(1.f, 0.f, 0.f);
  glLineWidth(3);
  if (CGAL::assign(foundVertex, location)) {
    // located a spherical vertex
    Projected_normal pn;
    pn.compute_projection((foundVertex->direction()).vector());
    float pnt[3];
    pnt[0]=CGAL::to_double((pn.get_projected_normal()).x())*s_scale;
    pnt[1]=CGAL::to_double((pn.get_projected_normal()).y())*s_scale;
    pnt[2]=CGAL::to_double((pn.get_projected_normal()).z())*s_scale;
    glPushMatrix();
    GLUquadricObj *qobj=gluNewQuadric();
    glTranslatef(pnt[0],pnt[1],pnt[2]);
    gluSphere(qobj, 0.02, 4, 4);
    gluDeleteQuadric(qobj);
    glPopMatrix();
  }
  if (CGAL::assign(foundHalfedge, location)) {
    // located an arc
    drawCubeArc(foundHalfedge);
  }
  if (CGAL::assign(foundFace, location)) {
    // located a face
    if (foundFace->does_outer_ccb_exist()) {
      // draw outer ccb
      SphereMap::Ccb_halfedge_circulator heif;
      heif = foundFace->outer_ccb();
      SphereMap::Ccb_halfedge_circulator circEn = heif;
      do {
        drawCubeArc(heif);
        ++heif;
      } while (heif != circEn);
    } // end of outer ccb display
    // display holes
    SphereMap::Holes_iterator holesIt;
    for (holesIt = foundFace->holes_begin(); holesIt != foundFace->holes_end();
        ++holesIt) {
      SphereMap::Ccb_halfedge_circulator heif;
      heif = *holesIt;
      SphereMap::Ccb_halfedge_circulator circEn=heif;
      do {
        // draw halfedge
        drawCubeArc(heif);
        ++heif;
      } while (heif !=circEn);
    }
  }
  glLineWidth(1);
  glEndList();

  s_showPointLocation = true;
  glLoadIdentity();
  gluLookAt(0.0,0.0,4.0f,
      0.0,0.0,0.0,
    0.0f,1.0f,0.0f);
  glRotatef(s_angY, 0.f,1.f,0.f);
  glRotatef(s_angX, 1.f,0.f,0.f);
  glutPostRedisplay();
}

/*
 output display keys help
*/
void displayHelpKeys(void) {
  std::cout << "keys functions:" << std::endl <<
    " 'A' - show/hide axis parallel circles" << std::endl <<
    " 'D' - toggle between display of spherical map and cubical map" << std::endl <<
    " 'P' - point locate, requires entering direction in rational format" << std::endl <<
    " 'R' - reset display, set starting position and clear point location data" << std::endl <<
    " 'C' - don't show the point location data" << std::endl <<
    " 'H' - display this help screen" << std::endl <<
    " 'Esc' - close viewer" << std::endl <<
    " 'Up/Down' - rotate around the X axis" << std::endl <<
    " 'Left/Right' - rotate around the Y axis" << std::endl;
}

/*
 handle ordinary keys (not the special up/down/left/right
*/
void processNormalKeys(unsigned char key, int x, int y) {

  switch (key) {
    case 27:
      exit(0);
      break;
    case 'a':
    case 'A':
      s_displayAxis = !s_displayAxis;
      glutPostRedisplay();
      break;
    case 'd':
    case 'D':
      s_displaySphere = !s_displaySphere;
      glutPostRedisplay();
      break;
    case 'p':
    case 'P':
      // point location
      pointLocate();
      break;
    case 'r':
    case 'R':
      glLoadIdentity();
      gluLookAt(0.0,0.0,4.0f,
      0.0,0.0,0.0,
      0.0f,1.0f,0.0f);
      s_angX = 0; s_angY = 0;
      s_showPointLocation = false;
      s_displayAxis = true;
      glutPostRedisplay();
      break;
    case 'c':
    case 'C':
      s_showPointLocation = false;
      glutPostRedisplay();
      break;
    case 'h':
    case 'H':
      displayHelpKeys();
      break;
  }
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " <data_file>" << std::endl;
    exit(1);
  }

  readSpericalMap(argv[1], s_spMap); // read map

  // set up the openGL display window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);

  glutInitWindowPosition(100,100);
  glutInitWindowSize(320,320);
  glutCreateWindow("Spherical map viewer");
  initScene(); // initialize drawing information
  glutKeyboardFunc(processNormalKeys);
  glutSpecialFunc(inputKey);
  glutDisplayFunc(renderScene);
  glutReshapeFunc(changeSize);
  glutMainLoop();

  return 0;
}
