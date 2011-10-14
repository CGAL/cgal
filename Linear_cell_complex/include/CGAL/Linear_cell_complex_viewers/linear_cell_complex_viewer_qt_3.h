// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_MAP_VIEWER_QT_H
#define CGAL_MAP_VIEWER_QT_H

#include <QApplication>
#include <QKeyEvent>

#include <QGLViewer/qglviewer.h>
#include <GL/gl.h>
#include <CGAL/Combinatorial_map_with_points_operations.h>

template<class Map>
CGAL::Bbox_3 bbox(Map& amap)
{
  CGAL::Bbox_3 bb;
  typename Map::Vertex_iterator it = amap.vertices_begin();
  if ( it!=amap.vertices_end() )
    {
      bb = it->point().bbox();
      for( ++it; it != amap.vertices_end(); ++it)
	{
	  bb = bb + it->point().bbox();
	}
    }
  
  return bb;
}

template<class Map>
class SimpleMapViewerQt : public QGLViewer
{
  typedef typename Map::Dart_handle Dart_handle;
  
 public:

  // Constructor/Destructor
  SimpleMapViewerQt(Map& amap) : QGLViewer(), m(amap),
    wireframe(false), flatShading(true), edges(true), vertices(true)
  {
    setWindowTitle("3D map viewer");	
    resize(500, 450);
  }

protected :
  // Draw the facet given by ADart
  void drawFacet(Dart_handle ADart, int AMark)
  {
    ::glBegin(GL_POLYGON);
    ::glColor3f(.7,.7,.7);

    // If Flat shading: 1 normal per polygon
    if (flatShading)
      {
	typename Map::Vector n = CGAL::compute_facet_normal(m,ADart);
	::glNormal3d(n.x(),n.y(),n.z());
      }

    for (typename Map::Dart_iterator_of_beta1 it(m,ADart); it.cont(); ++it)
      {	
	// If Gouraud shading: 1 normal per vertex
	if (!flatShading)
	  {
	    typename Map::Vector n = CGAL::compute_vertex_normal<Map>(m,*it);
	    ::glNormal3d(n.x(),n.y(),n.z());
	  }
	
	typename Map::Point p = it->vertex()->point();
	::glVertex3d( p.x(),p.y(),p.z());
	
	m.set_mark(*it,AMark);
	if ( !it->is_free(3) ) m.set_mark(it->beta(3),AMark);
      }
    ::glEnd();
  }

  /// Draw all the edge of the facet given by ADart
  void drawEdges(Dart_handle ADart)
  {
    glBegin(GL_LINES);
    glColor3f(.2,.2,.6);
    for (typename Map::Dart_iterator_of_beta1 it(m,ADart); it.cont(); ++it)
      {
	typename Map::Point p = it->vertex()->point();
	Dart_handle d2 = it->opposite();
	if ( d2!=NULL )
	  {
	    typename Map::Point p2 = d2->vertex()->point();
	    glVertex3f( p.x(),p.y(),p.z());
	    glVertex3f( p2.x(),p2.y(),p2.z());
	  }	
      }
    glEnd();
  }   
  
  virtual void draw()
  {

    unsigned int facettreated = m.get_new_mark();
    unsigned int vertextreated = -1;

    if ( vertices) vertextreated=m.get_new_mark();

    for(typename Map::Dart_iterator_of_all it(m); it.cont(); ++it)
      {
	if ( !m.is_marked(*it,facettreated) )
	  {
	    drawFacet(*it,facettreated);
	    if ( edges) drawEdges(*it);
	  }

	if (vertices)
	  {
	    if ( !m.is_marked(*it, vertextreated) )
	      {	    
		typename Map::Point p = it->vertex()->point();
		
		glBegin(GL_POINTS);
		glVertex3f( p.x(),p.y(),p.z());
		glEnd();
		
		mark_orbit(m,*it,Map::VERTEX_ORBIT,vertextreated);
	      }
	  }
      }

    assert(m.is_whole_map_marked(facettreated));

    if ( vertices)
      {
	assert(m.is_whole_map_marked(vertextreated));
	m.free_mark(vertextreated);
      }
    
    m.free_mark(facettreated);  
  }
  
  virtual void init()
  {
    // Restore previous viewer state.
    restoreStateFromFile();

    // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
    setShortcut(EXIT_VIEWER, Qt::CTRL+Qt::Key_Q);

    // Add custom key description (see keyPressEvent).
    setKeyDescription(Qt::Key_W, "Toggles wire frame display");
    setKeyDescription(Qt::Key_F, "Toggles flat shading display");
    setKeyDescription(Qt::Key_E, "Toggles edges display");
    setKeyDescription(Qt::Key_V, "Toggles vertices display");

    // Light default parameters
    ::glLineWidth(1.0f);
    ::glPointSize(3.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(1.0f,1.0f);
    ::glClearColor(1.0f,1.0f,1.0f,0.0f);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

    ::glEnable(GL_LIGHTING);
    
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    // ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

    if (flatShading)
      {
	::glShadeModel(GL_FLAT);
	::glDisable(GL_BLEND); 
	::glDisable(GL_LINE_SMOOTH); 
	::glDisable(GL_POLYGON_SMOOTH_HINT); 
	::glBlendFunc(GL_ONE, GL_ZERO); 
	::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
      }
    else
      {
	::glShadeModel(GL_SMOOTH);
	::glEnable(GL_BLEND);
	::glEnable(GL_LINE_SMOOTH);
	::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }


    CGAL::Bbox_3 bb = bbox(m);
    
    this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
						       bb.ymin(),
						       bb.zmin()),
					qglviewer::Vec(bb.xmax(),
						       bb.ymax(),
						       bb.zmax()));
    
    this->showEntireScene();
  }
  
  void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    bool handled = false;
    if ((e->key()==Qt::Key_W) && (modifiers==Qt::NoButton))
      {
	wireframe = !wireframe;
	if (wireframe)
	  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
	  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	handled = true;
	updateGL();
      }
    else if ((e->key()==Qt::Key_F) && (modifiers==Qt::NoButton))
	{
	  flatShading = !flatShading;
	  if (flatShading)
	    {
	      ::glShadeModel(GL_FLAT);
	      ::glDisable(GL_BLEND); 
	      ::glDisable(GL_LINE_SMOOTH); 
	      ::glDisable(GL_POLYGON_SMOOTH_HINT); 
	      ::glBlendFunc(GL_ONE, GL_ZERO); 
	      ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
	    }
	  else
	    {
	      ::glShadeModel(GL_SMOOTH);
	      ::glEnable(GL_BLEND);
	      ::glEnable(GL_LINE_SMOOTH);
	      ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	      ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	    }
	  handled = true;
	  updateGL();
	}
    else if ((e->key()==Qt::Key_E) && (modifiers==Qt::NoButton))
	{
	  edges = !edges;
	  handled = true;
	  updateGL();
	}
    else if ((e->key()==Qt::Key_V) && (modifiers==Qt::NoButton))
	{
	  vertices = !vertices;
	  handled = true;
	  updateGL();
	}
    
    if (!handled)
      QGLViewer::keyPressEvent(e);
  }


  virtual QString helpString() const
  {
    QString text("<h2>M a p   V i e w e r</h2>");
    text += "Use the mouse to move the camera around the object. ";
    text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
    text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
    text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
    text += "Simply press the function key again to restore it. Several keyFrames define a ";
    text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
    text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
    text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
    text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
    text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
    text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
    text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
    text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
    text += "Press <b>Escape</b> to exit the viewer.";
    return text;
  }
 private:
  Map& m;
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;
};

template<class Map>
void display_map(Map& amap)
{
  int argc=1;
  typedef char* s;
  
  const char* argv[2]={"mapviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));
  
  SimpleMapViewerQt<Map> mainwindow(amap);
  mainwindow.show();

  app.exec();
};

#endif // CGAL_MAP_VIEWER_QT_H

