//! \file Viewer.h

#ifndef VIEWER_H
#define VIEWER_H

#include <CGAL/Three/Viewer_config.h>
#include <CGAL/Three/Scene_interface.h>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>
#include <CGAL/Three/Viewer_interface.h>

#include <QGLViewer/qglviewer.h>
#include <QPoint>
#include <QFont>
#include <QOpenGLFramebufferObject>
#include <CGAL/Three/TextRenderer.h>
// forward declarations
class QWidget;
namespace CGAL{
namespace Three{
class Scene_draw_interface;
}
}
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;

class Viewer_impl;
//!The viewer class. Deals with all the openGL rendering and the mouse/keyboard events.
class VIEWER_EXPORT Viewer : public CGAL::Three::Viewer_interface {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  ~Viewer();
  bool testDisplayId(double, double, double);
  // overload several QGLViewer virtual functions
  //! Draws the scene.
  void draw();
  //!This step happens after draw(). It is here that all the useful information is displayed, like the axis system or the informative text.
  void drawVisualHints();
  //! Deprecated. Does the same as draw().
  void fastDraw();
  //! Initializes the OpenGL functions and sets the backGround color.
  void initializeGL();
  //! Draws the scene "with names" to allow picking.
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
  void setScene(CGAL::Three::Scene_draw_interface* scene);
  //! @returns the antialiasing state.
  bool antiAliasing() const;
  //! @returns the fastDrawing state.
  bool inFastDrawing() const;
  //! Implementation of `Viewer_interface::inDrawWithNames()`
  bool inDrawWithNames() const;
  //! Implementation of `Viewer_interface::attribBuffers()`
  void attribBuffers(int program_name) const;
  //! Implementation of `Viewer_interface::getShaderProgram()`
  QOpenGLShaderProgram* getShaderProgram(int name) const;
  QPainter* getPainter();
  void saveSnapshot(bool , bool overwrite = false);

Q_SIGNALS:
  void sendMessage(QString);
public Q_SLOTS:
  //! Sets the antialiasing to true or false.
  void setAntiAliasing(bool b);
  //! If b is true, facets will be ligted from both internal and external sides.
  //! If b is false, only the side that is exposed to the light source will be lighted.
  void setTwoSides(bool b);
  //! If b is true, some items are displayed in a simplified version when moving the camera.
  //! If b is false, items display is never altered, even when moving.
  void setFastDrawing(bool b);
  //! Makes the camera turn around.
  void turnCameraBy180Degres();
  //! @returns a QString containing the position and orientation of the camera.
  QString dumpCameraCoordinates();
  //!Moves the camera to the new coordinates (position and orientation) through an animation.
  bool moveCameraToCoordinates(QString, 
                               float animation_duration = 0.5f);
  //!Makes the Viewer display a message
  void printMessage(QString message, int ms_delay );
  void displayMessage(const QString &_message, int delay);
  void displayMessage(const QString &_message){displayMessage(_message, 2000);}
  void hideMessage();

protected:
  void postDraw();
  void paintEvent(QPaintEvent *);
  void paintGL();
  //! Holds useful data to draw the axis system
  struct AxisData
  {
      std::vector<float> *vertices;
      std::vector<float> *normals;
      std::vector<float> *colors;
  };

  //! The buffers used to draw the axis system
  QOpenGLBuffer buffers[4];
  //! The VAO used to draw the axis system
  QOpenGLVertexArrayObject vao[2];
  //! The rendering program used to draw the axis system
  QOpenGLShaderProgram rendering_program;
  //! The rendering program used to draw the distance
  QOpenGLShaderProgram rendering_program_dist;
  QList<TextItem*>  distance_text;
  //! Holds the vertices data for the axis system
  std::vector<float> v_Axis;
  //! Holds the normals data for the axis system
  std::vector<float> n_Axis;
  //! Holds the color data for the axis system
  std::vector<float> c_Axis;
  //! Decides if the axis system must be drawn or not
  bool axis_are_displayed;
  //! Decides if the text is displayed in the drawVisualHints function.
  bool has_text;
  //! Decides if the distance between APoint and BPoint must be drawn;
  bool distance_is_displayed;
  //!Defines the behaviour for the mouse press events
  void mousePressEvent(QMouseEvent*);
  void wheelEvent(QWheelEvent *);
  //!Defines the behaviour for the key press events
  void keyPressEvent(QKeyEvent*);
  //!Deal with context menu events
  void contextMenuEvent(QContextMenuEvent*);
  //!Defines the behaviour for the key release events
  void keyReleaseEvent(QKeyEvent *);
  //!Clears the distance display
  void clearDistancedisplay();
  /*!
   * \brief makeArrow creates an arrow and stores it in a struct of vectors.
   * \param R the radius of the arrow.
   * \param prec the precision of the quadric. The lower this value is, the higher precision you get.
   * It can be any int between 1 and 360.
   * \param from the starting point of the arrow.
   * \param to the destination point of the arrow (the pointed extremity).
   * \param color the RGB color of the arrow.
   * \param data the struct of std::vector that will contain the results.
   */

  void makeArrow(double R, int prec, qglviewer::Vec from, qglviewer::Vec to, qglviewer::Vec color, AxisData &data);
  void resizeGL(int w, int h);
  bool i_is_pressed;

  //!Draws the distance between two selected points.
  void showDistance(QPoint);
  qglviewer::Vec APoint;
  qglviewer::Vec BPoint;
  bool is_d_pressed;

protected:
  friend class Viewer_impl;
  Viewer_impl* d;
  double prev_radius;

private:
  // F P S    d i s p l a y
  QTime fpsTime;
  unsigned int fpsCounter;
  QString fpsString;
  float f_p_s;
  // M e s s a g e s
  QString message;
  bool _displayMessage;
  QTimer messageTimer;
}; // end class Viewer


#endif // VIEWER_H
