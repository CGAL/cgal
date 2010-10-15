// This file is part of GlSplat, a simple splatting C++ library
//
// Copyright (C) 2008-2009 Gael Guennebaud <g.gael@free.fr>
//
// GlSplat is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// GlSplat is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with GlSplat. If not, see <http://www.gnu.org/licenses/>.

#include "GlSplat.h"
#include <QGLViewer/qglviewer.h>

class Viewer : public QGLViewer
{
protected :
  virtual void draw();
  virtual void init();
  virtual QString helpString() const;

  virtual void drawpoints();
  GlSplat::SplatRenderer mRenderer;
  std::vector<float> mNormals;
  std::vector<float> mPositions;
  std::vector<float> mRadii;
  std::vector<unsigned int> mColors;
  int mNbPoints;
};

void Viewer::draw()
{
  mRenderer.beginVisibilityPass();
  drawpoints();
  mRenderer.beginAttributePass();
  drawpoints();
  mRenderer.finalize();
}

void Viewer::drawpoints()
{
  glBegin(GL_POINTS);
  for (int i=0; i<mNbPoints; ++i)
  {
    glMultiTexCoord1f(GL_TEXTURE2, mRadii[i]);
    glNormal3fv(&mNormals[i*3]);
    glColor4ubv((const unsigned char *)&mColors[i]);
    glVertex3fv(&mPositions[i*3]);
  }
  glEnd();
}

void Viewer::init()
{
  // Restore previous viewer state.
  restoreStateFromFile();

  mNbPoints = 1000;
  mPositions.resize(3*mNbPoints);
  mNormals.resize(3*mNbPoints);
  mRadii.resize(mNbPoints);
  mColors.resize(mNbPoints);

  for (int i=0; i<mNbPoints; ++i)
  {
    mPositions[i*3 + 0] = drand48() - 0.5;
    mPositions[i*3 + 1] = drand48() - 0.5;
    mPositions[i*3 + 2] = drand48() - 0.5;
    float l = 1./sqrtf(mPositions[i*3+0]*mPositions[i*3+0] + mPositions[i*3+1]*mPositions[i*3+1] + mPositions[i*3+2]*mPositions[i*3+2]);
    unsigned char rgba[4];
    for (int k=0; k<3; ++k)
    {
      mNormals[i*3+k] = (mPositions[i*3+k] *= l);
      rgba[k] = (unsigned char)((mNormals[i*3+k]*0.4+0.6)*255);
    }
    rgba[3] = 0xff;
    mRadii[i] = 0.1;
    mColors[i] = *(unsigned int*)rgba;
  }

  mRenderer.init();

  // Opens help window
//   help();
}

QString Viewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
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


#include <qapplication.h>

int main(int argc, char** argv)
{
  // Read command lines arguments.
  QApplication application(argc,argv);

  // Instantiate the viewer.
  Viewer viewer;

  viewer.setWindowTitle("simpleViewer");

  // Make the viewer window visible on screen.
  viewer.show();

  // Run main loop.
  return application.exec();
}

