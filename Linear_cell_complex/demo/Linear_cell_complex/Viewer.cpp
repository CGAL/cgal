// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//                 Kumar Snehasish <kumar.snehasish@gmail.com>
//
#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <QGLViewer/vec.h>
#include <CGAL/Linear_cell_complex_operations.h>

void
Viewer::sceneChanged()
{
  this->initDraw();

  this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
						     bb.ymin(),
						     bb.zmin()),
				      qglviewer::Vec(bb.xmax(),
						     bb.ymax(),
						     bb.zmax()));
  this->showEntireScene();
}

void Viewer::drawAllFaces(bool flat)
{
  LCC &lcc = *scene->lcc;

  for (LCC::Attribute_range<3>::type::iterator
         it=lcc.attributes<3>().begin(),
         itend=lcc.attributes<3>().end(); it!=itend; ++it )
  {
    if ( it->info().is_visible() )
    {
      for(LCC::One_dart_per_incident_cell_range<2,3>::iterator
            dartIter=lcc.one_dart_per_incident_cell<2,3>
            (lcc.dart_of_attribute<3>(it)).begin(); dartIter.cont(); ++dartIter)
      {
        // We draw the polygon
        ::glBegin(GL_POLYGON);

        //  double r = (double)dartIter->attribute<3>()->info().r()/255.0;
        double r = (double)lcc.info<3>(dartIter).color().r()/255.0;
        double g = (double)lcc.info<3>(dartIter).color().g()/255.0;
        double b = (double)lcc.info<3>(dartIter).color().b()/255.0;
        if ( !lcc.is_free(dartIter, 3) )
        {
          r += (double)lcc.info<3>(lcc.beta(dartIter,3)).color().r()/255.0;
          g += (double)lcc.info<3>(lcc.beta(dartIter,3)).color().g()/255.0;
          b += (double)lcc.info<3>(lcc.beta(dartIter,3)).color().b()/255.0;
          r /= 2; g /= 2; b /= 2;
        }

        ::glColor3f(r,g,b);

        if(flat)
        {
          LCC::Vector normal = CGAL::compute_normal_of_cell_2(lcc,dartIter);
          normal = normal/(CGAL::sqrt(normal*normal));
          ::glNormal3d(normal.x(), normal.y(), normal.z());
        }

        for (LCC::Dart_of_orbit_range<1>::const_iterator
               orbitIter = lcc.darts_of_orbit<1>(dartIter).begin();
             orbitIter.cont(); ++orbitIter)
        {
          if(!flat)
          {
            // If Gouraud shading: 1 normal per vertex
            LCC::Vector normal = CGAL::compute_normal_of_cell_0(lcc,orbitIter);
            normal = normal/(CGAL::sqrt(normal*normal));
            ::glNormal3d(normal.x(), normal.y(), normal.z());
          }

          const LCC::Point& p = lcc.point(orbitIter);
          ::glVertex3d(p.x(),p.y(),p.z());
        }

        ::glEnd();
      }
    }
  }
}

void Viewer::drawAllEdges()
{
  LCC &lcc = *scene->lcc;

  if ( lcc.is_empty() ) return;

  //    ::glDepthRange(0.0, 1.0-0.005);
  ::glBegin(GL_LINES);
  //::glColor3f(0.0f, 0.0f, 0.0f);
  ::glColor3f(.2f,.2f,.6f);

  for (LCC::Attribute_range<3>::type::iterator
         it=lcc.attributes<3>().begin(),
         itend=lcc.attributes<3>().end(); it!=itend; ++it )
  {
    if ( it->info().is_visible() )
    {
      for(LCC::One_dart_per_incident_cell_range<1,3>::iterator
            dartIter=lcc.one_dart_per_incident_cell<1,3>
            (lcc.dart_of_attribute<3>(it)).begin(); dartIter.cont(); ++dartIter)
      {
        const LCC::Point& p =  lcc.point(dartIter);
        Dart_handle d2 = lcc.other_extremity(dartIter);
        if ( d2!=NULL )
        {
          const LCC::Point& p2 = lcc.point(d2);
          glVertex3f( p.x(),p.y(),p.z());
          glVertex3f( p2.x(),p2.y(),p2.z());
        }
      }
    }
  }

  ::glEnd();
}

void Viewer::drawAllVertices()
{
  LCC &lcc = *scene->lcc;

  if ( lcc.is_empty() )
  {
    bb = LCC::Point(CGAL::ORIGIN).bbox();
    bb = bb + LCC::Point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    return;
  }

  //    ::glDepthRange(0.0, 1.0-0.005);
  ::glPointSize(7.0);
  ::glBegin(GL_POINTS);
  ::glColor3f(0.2f, 0.2f, 0.7f);

  bool empty = true;
  for (LCC::Attribute_range<3>::type::iterator
         it=lcc.attributes<3>().begin(),
         itend=lcc.attributes<3>().end(); it!=itend; ++it )
  {
    if ( it->info().is_visible() )
    {
      for(LCC::One_dart_per_incident_cell_range<0,3>::iterator
            dartIter=lcc.one_dart_per_incident_cell<0,3>
            (lcc.dart_of_attribute<3>(it)).begin();
          dartIter.cont(); ++dartIter)
      {
        const LCC::Point& p =  lcc.point(dartIter);
        ::glVertex3f(p.x(), p.y(), p.z());

        if ( empty )
        {
          bb = p.bbox();
          empty = false;
        }
        else
          bb = bb + p.bbox();
      }
    }
  }

  ::glEnd();

  if ( lcc.is_empty() )
  {
    bb = LCC::Point(CGAL::ORIGIN).bbox();
    bb = bb + LCC::Point(1,1,1).bbox(); // To avoid a warning from Qglviewer
  }
}

void Viewer::initDraw()
{
    //Compile drawFacet
    //    std::cout << "Compile Display Lists : Faces, " << std::flush;
    m_dlFaces = ::glGenLists(1);
    ::glNewList(m_dlFaces, GL_COMPILE);
    drawAllFaces(false);
    ::glEndList();

    //Compile drawFacet with flat shading
    //    std::cout << "Faces (flat shading), " << std::flush;
    m_dlFacesFlat = ::glGenLists(1);
    ::glNewList(m_dlFacesFlat, GL_COMPILE);
    drawAllFaces(true);
    ::glEndList();

    //Compile drawEdge
    //    std::cout << "edges, " << std::flush;
    m_dlEdges = ::glGenLists(1);
    ::glNewList(m_dlEdges, GL_COMPILE);
    drawAllEdges();
    ::glEndList();

    //Compile drawvertices
    //    std::cout << "vertices" << std::flush;
    m_dlVertices = ::glGenLists(1);
    ::glNewList(m_dlVertices, GL_COMPILE);
    drawAllVertices();
    ::glEndList();

    //    std::cout << ". DONE." << std::endl;
    m_displayListCreated = true;
}

void Viewer::draw()
{
  if(!m_displayListCreated) initDraw();

  if ( !wireframe )
  {
    if(flatShading) ::glCallList(m_dlFacesFlat);
    else ::glCallList(m_dlFaces);
  }

  if(edges) ::glCallList(m_dlEdges);

  if(vertices) ::glCallList(m_dlVertices);
}

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

  if (!handled)
    QGLViewer::keyPressEvent(e);
}

QString Viewer::helpString() const
{
  QString text("<h2>L C C   V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with "
    "the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the "
    "camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys "
    "(<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several "
    "keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and "
    "restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the "
    "world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to "
    "save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete "
    "shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button "
    "double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and "
    "the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed "
    "defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for "
    "details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}

#include "Viewer.moc"
