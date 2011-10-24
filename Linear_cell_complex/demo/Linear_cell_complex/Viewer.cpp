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
#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <QGLViewer/vec.h>
#include <CGAL/Linear_cell_complex_operations.h>

#define NB_FILLED_MODE   4
#define FILLED_ALL       0
#define FILLED_NON_FREE3 1
#define FILLED_VOL       2
#define FILLED_VOL_AND_V 3

template<class LCC>
CGAL::Bbox_3 bbox(LCC& alcc)
{
  CGAL::Bbox_3 bb;
  typename LCC::Vertex_attribute_range::iterator it = alcc.vertex_attributes().begin(),
    itend=alcc.vertex_attributes().end();
  if ( it!=itend )
    {
      bb = it->point().bbox();
      for( ++it; it != itend; ++it)
	{
	  bb = bb + it->point().bbox();
	}
    }
  
  return bb;
}

void
Viewer::sceneChanged()
{
  iteratorAllDarts = scene->lcc->darts().begin();
  scene->lcc->unmark_all(markVolume);
  
  CGAL::Bbox_3 bb = bbox(*scene->lcc);
   
  this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
						     bb.ymin(),
						     bb.zmin()),
				      qglviewer::Vec(bb.xmax(),
						     bb.ymax(),
						     bb.zmax()));
    
  this->showEntireScene();
}

// Draw the facet given by ADart
void Viewer::drawFacet(Dart_handle ADart, int AMark)
{  
  LCC &m = *scene->lcc;
  ::glBegin(GL_POLYGON);
#ifdef COLOR_VOLUME
  assert( ADart->attribute<3>()!=NULL );

  //  double r = (double)ADart->attribute<3>()->info().r()/255.0;
  double r = (double)ADart->attribute<3>()->info().r()/255.0;
  double g = (double)ADart->attribute<3>()->info().g()/255.0;
  double b = (double)ADart->attribute<3>()->info().b()/255.0;
  if ( !ADart->is_free(3) )
    {
      r += (double)ADart->beta(3)->attribute<3>()->info().r()/255.0;
      g += (double)ADart->beta(3)->attribute<3>()->info().g()/255.0;
      b += (double)ADart->beta(3)->attribute<3>()->info().b()/255.0;
      r /= 2; g /= 2; b /= 2;
    }

  ::glColor3f(r,g,b);
#else 
  ::glColor3f(.7,.7,.7);
#endif

  // If Flat shading: 1 normal per polygon
  if (flatShading)
    {
      LCC::Vector n = CGAL::compute_normal_of_cell_2(m,ADart);
      n = n/(CGAL::sqrt(n*n));
      ::glNormal3d(n.x(),n.y(),n.z());
    }

  for ( LCC::Dart_of_orbit_range<1>::iterator it(m,ADart); it.cont(); ++it)
    {	
      // If Gouraud shading: 1 normal per vertex
      if (!flatShading)
	{
	  LCC::Vector n = CGAL::compute_normal_of_cell_0<LCC>(m,it);
	  n = n/(CGAL::sqrt(n*n));	  
	  ::glNormal3d(n.x(),n.y(),n.z());
	}
	
      LCC::Point p = m.point(it);
      ::glVertex3d( p.x(),p.y(),p.z());
	
      m.mark(it,AMark);
      if ( !it->is_free(3) ) m.mark(it->beta(3),AMark);
    }
  ::glEnd();
}

/// Draw all the edge of the facet given by ADart
void Viewer::drawEdges(Dart_handle ADart)
{ 
  LCC &m = *scene->lcc;
  glBegin(GL_LINES);
  glColor3f(.2,.2,.6);
  for ( LCC::Dart_of_orbit_range<1>::iterator it(m,ADart); it.cont(); ++it)
    {
      LCC::Point p = m.point(it);
      Dart_handle d2 = it->other_extremity();
      if ( d2!=NULL )
	{
	  LCC::Point p2 = m.point(d2);
	  glVertex3f( p.x(),p.y(),p.z());
	  glVertex3f( p2.x(),p2.y(),p2.z());
	}	
    }
  glEnd();
}
    
void Viewer::draw_one_vol_filled_facets(Dart_handle adart,
				       int amarkvol, int amarkfacet)
{
  LCC &m = *scene->lcc;
  
  for (CGAL::CMap_dart_iterator_basic_of_cell<LCC,3> it(m,adart,amarkvol); it.cont(); ++it)
    {
      if ( !m.is_marked(it,amarkfacet) )
	{
	  drawFacet(it,amarkfacet);
	}
    }  
}

void Viewer::draw_current_vol_filled_facets(Dart_handle adart)
{
  LCC &m = *scene->lcc;
  unsigned int facettreated = m.get_new_mark();
  unsigned int volmark     = m.get_new_mark();

  draw_one_vol_filled_facets(adart,volmark,facettreated); 

  m.negate_mark(volmark);  
  
  for (CGAL::CMap_dart_iterator_basic_of_cell<LCC,3> it(m,adart,volmark); it.cont(); ++it)
    {
      m.unmark(it,facettreated);
      if ( !it->is_free(3) ) m.unmark(it->beta(3),facettreated);
    }

  m.negate_mark(volmark);  
  
  assert(m.is_whole_map_unmarked(volmark));
  assert(m.is_whole_map_unmarked(facettreated));
  
  m.free_mark(volmark);
  m.free_mark(facettreated);  
}

void Viewer::draw_current_vol_and_neighboors_filled_facets(Dart_handle adart)
{
  LCC &m = *scene->lcc;
  unsigned int facettreated = m.get_new_mark();
  unsigned int volmark     = m.get_new_mark();
  
  draw_one_vol_filled_facets(adart,volmark,facettreated);

  CGAL::CMap_dart_iterator_of_cell<LCC,3> it(m,adart);
  for (; it.cont(); ++it)
    {
      if ( !it->is_free(3) && !m.is_marked(it->beta(3),volmark) )
	{
	  draw_one_vol_filled_facets(it->beta(3),volmark,facettreated);
	}
    }

  m.negate_mark(volmark);
  
  for (it.rewind(); it.cont(); ++it)
    {
      m.mark(it,volmark);
	    
      if ( m.is_marked(it,facettreated))
	CGAL::unmark_cell<LCC,2>(m,it,facettreated);
      
      if ( !it->is_free(3) && !m.is_marked(it->beta(3),volmark) )
	{
	  CGAL::CMap_dart_iterator_basic_of_cell<LCC,3> it2(m,it->beta(3),volmark);
	  for (; it2.cont(); ++it2)
	    {
	      if ( m.is_marked(it2,facettreated))
		CGAL::unmark_cell<LCC,2>(m,it2,facettreated);
	    }
	}
    }

  m.negate_mark(volmark);

  assert(m.is_whole_map_unmarked(volmark));
  assert(m.is_whole_map_unmarked(facettreated));

  m.free_mark(volmark);
  m.free_mark(facettreated);  
}

void Viewer::draw()
{
  LCC &m = *scene->lcc;

  if ( m.is_empty() ) return;

  unsigned int facettreated = m.get_new_mark();
  unsigned int vertextreated = -1;

  if ( vertices) vertextreated=m.get_new_mark();

  for(LCC::Dart_range::iterator it=m.darts().begin(); it!=m.darts().end(); ++it)
    {
      if ( !m.is_marked(it,facettreated) )
	{
	  if ( modeFilledFacet==FILLED_ALL ||
	      modeFilledFacet==FILLED_NON_FREE3 && !it->is_free(3) )
	    drawFacet(it,facettreated);
	  else
	    CGAL::mark_cell<LCC,2>(m,it,facettreated);

	  if ( edges) drawEdges(it);
	}

      if (vertices)
	{
	  if ( !m.is_marked(it, vertextreated) )
	    {	    
	      LCC::Point p = m.point(it);
		
	      glBegin(GL_POINTS);
	      glColor3f(.6,.2,.8);
	      glVertex3f( p.x(),p.y(),p.z());
	      glEnd();
		
	      CGAL::mark_cell<LCC,0>(m,it,vertextreated);
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

  if ( modeFilledFacet==FILLED_VOL) 
    draw_current_vol_filled_facets(iteratorAllDarts);
  else if ( modeFilledFacet==FILLED_VOL_AND_V)
    draw_current_vol_and_neighboors_filled_facets(iteratorAllDarts);
}

/*
  void
  Viewer::draw()
  {

  // define material
  float	ambient[]  =   { 0.25f,
  0.20725f,
  0.20725f,
  0.922f };
  float	diffuse[]  =   { 1.0f,
  0.829f,
  0.829f,
  0.922f };

  float	specular[]  = {  0.296648f,
  0.296648f,
  0.296648f,
  0.522f };

  float	emission[]  = {  0.3f,
  0.3f,
  0.3f,
  1.0f };
  float shininess[] = {  11.264f };

  // apply material
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);

  // anti-aliasing (if the OpenGL driver permits that)
  ::glEnable(GL_LINE_SMOOTH);

  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 
  // draw surface mesh
  bool m_view_surface = true;
  bool draw_triangles_edges = true;

  if(m_view_surface)
  {
  ::glEnable(GL_LIGHTING);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  ::glColor3f(0.2f, 0.2f, 1.f);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(3.0f,-3.0f);
  gl_draw_surface();

  if(draw_triangles_edges)
  {
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1.);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  ::glColor3ub(0,0,0);
  ::glDisable(GL_POLYGON_OFFSET_FILL);
  gl_draw_surface();
  }
  }

  }


  void 
  Viewer::gl_draw_surface()
  {
  ::glColor3f(1.0f, 0.0f, 0.0f);
  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_POINT_SMOOTH);
  ::glPointSize(5);
  ::glBegin(GL_POINTS);

  for(std::list<Point_3>::iterator it = scene->points.begin();
  it != scene->points.end();
  ++it){
  ::glVertex3d(it->x(), it->y(), it->z());
  }

  ::glEnd();
  ::glDisable(GL_POINT_SMOOTH);

  ::glEnable(GL_LIGHTING);
  ::glBegin(GL_TRIANGLES);

  ::glColor3f(0.2f, 1.0f, 0.2f);

  std::list<Facett> facetts;
  scene->alpha_shape.get_alpha_shape_facetts(std::back_inserter(facetts), Alpha_shape_3::REGULAR);
  
  for(std::list<Facett>::iterator fit = facetts.begin();
  fit != facetts.end();
  ++fit) {
  const Cell_handle& ch = fit->first;
  const int index = fit->second;
    
  //const Vector_3& n = ch->normal(index); // must be unit vector
    
  const Point_3& a = ch->vertex((index+1)&3)->point();
  const Point_3& b = ch->vertex((index+2)&3)->point();
  const Point_3& c = ch->vertex((index+3)&3)->point();
   
  Vector_3 v = CGAL::unit_normal(a,b,c);


  ::glNormal3d(v.x(),v.y(),v.z());
  ::glVertex3d(a.x(),a.y(),a.z());
  ::glVertex3d(b.x(),b.y(),b.z());
  ::glVertex3d(c.x(),c.y(),c.z());
  }

  
  ::glEnd();

  }
*/


void Viewer::init()
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
  setKeyDescription(Qt::Key_Z, "Next mode filled facet");
  setKeyDescription(Qt::Key_R, "Select next volume, used for filled facet");

  // Light default parameters
  ::glLineWidth(1.4f);
  ::glPointSize(4.f);
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
}

void Viewer::keyPressEvent(QKeyEvent *e)
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
  else if ((e->key()==Qt::Key_Z) && (modifiers==Qt::NoButton))
    {
      modeFilledFacet = (modeFilledFacet+1)%NB_FILLED_MODE;
      handled = true;
      updateGL();
    }
  else if ((e->key()==Qt::Key_R) && (modifiers==Qt::NoButton))
    {      
      CGAL::mark_cell<LCC,3>(*scene->lcc, iteratorAllDarts, markVolume);

      while ( iteratorAllDarts!=scene->lcc->darts().end() &&
	      scene->lcc->is_marked(iteratorAllDarts,markVolume) )
	{
	  ++iteratorAllDarts;
	}
      
      if ( iteratorAllDarts==scene->lcc->darts().end() )
	{
	  scene->lcc->negate_mark(markVolume);
	  assert( scene->lcc->is_whole_map_unmarked(markVolume) );
	  iteratorAllDarts=scene->lcc->darts().begin();
	}

      handled = true;
      updateGL();
    }
    
  if (!handled)
    QGLViewer::keyPressEvent(e);
}

QString Viewer::helpString() const
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
  text += "In filled facet, there are four modes: all facets are filled; only facets between two volumes are filles; only the facets of current volume are filled; only the facets of current volume and all its adjacent volumes are filled.";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}

#include "Viewer.moc"
