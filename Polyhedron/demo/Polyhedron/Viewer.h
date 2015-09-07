#ifndef VIEWER_H
#define VIEWER_H

#include "Viewer_config.h"
#include <CGAL_demo/Viewer_interface.h>

#include <QGLViewer/qglviewer.h>
#include <QPoint>

// forward declarations
class QWidget;
class Scene_draw_interface;
class QMouseEvent;
class QKeyEvent;

class Viewer_impl;

class VIEWER_EXPORT Viewer : public Viewer_interface {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  ~Viewer();

  // overload several QGLViewer virtual functions
  //! Deprecated and does nothing.
  void draw();
  //! Deprecated and does nothing.
  void fastDraw();
  //! Initializes the OpenGL functions and sets the backGround color.
  void initializeGL();
  //! Deprecated and does nothing.
  void drawWithNames();
  /*! Uses the parameter pixel's coordinates to get the corresponding point
   * in the World frame. If this point is found, emits selectedPoint, selected,
   * and selectionRay signals.
   */
  void postSelection(const QPoint&);
  //! Sets the picking matrix to allow the picking.
  void beginSelection(const QPoint &point);
  //! Sets the pick matrix to Identity once the picking is done.
  void endSelection(const QPoint &point);
  //! Sets the scene for the viewer.
  void setScene(Scene_draw_interface* scene);
  //! @returns the antialiasing state.
  bool antiAliasing() const;
  //! @returns the fastDrawing state.
  bool inFastDrawing() const;

public Q_SLOTS:
  //! Sets the antialiasing to true or false.
  void setAntiAliasing(bool b);
  //! If b is true, facets will be ligted from both internal and external sides.
  //! If b is false, only the side that is exposed to the light source will be lighted.
  void setTwoSides(bool b);
  //! Make the camera turn around.
  void turnCameraBy180Degres();
  //! @returns a QString containing the position and orientation of the camera.
  QString dumpCameraCoordinates();
  //!Moves the camera to the new coordinates (position and orientation) through an animation.
  bool moveCameraToCoordinates(QString, 
                               float animation_duration = 0.5f);

protected:
  //!Defines the behaviour for the mouse press events
  void mousePressEvent(QMouseEvent*);
  //!Defines the behaviour for the key press events
  void keyPressEvent(QKeyEvent*);
  /*! \brief Encapsulates the pickMatrix.
  * Source code of gluPickMatrix slightly modified : instead of multiplying the current matrix by this value,
  * sets the viewer's pickMatrix_ so that the drawing area is only around the cursor. This is because since CGAL 4.7,
  * the drawing system changed to use shaders, and these need this value. pickMatrix_ is passed to the shaders in
  * Scene_item::attrib_buffers(Viewer_interface* viewer, int program_name).*/
  void pickMatrix(GLdouble x, GLdouble y, GLdouble width, GLdouble height,
                  GLint viewport[4]);

protected:
  Viewer_impl* d;
}; // end class Viewer

#endif // VIEWER_H
