//! \file Viewer.h

#ifndef VIEWER_H
#define VIEWER_H

#include <CGAL/Three/Viewer_config.h>
#include <CGAL/Three/Scene_interface.h>
#include <QOpenGLBuffer>
#include <QOpenGLDebugMessage>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>
#include <CGAL/Three/Viewer_interface.h>
#include <QPoint>
#include <QFont>
#include <QOpenGLFramebufferObject>
#include <CGAL/Three/TextRenderer.h>
// forward declarations
class QWidget;
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;
class Viewer_impl;
namespace CGAL{
namespace Three{
class Scene_draw_interface;
}
}
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;

class Viewer_impl;
//! The viewer class. Deals with all the openGL rendering and the mouse/keyboard events.
//! It should not be needed in the plugin.
class VIEWER_EXPORT Viewer : public CGAL::Three::Viewer_interface {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  Viewer(QWidget * parent, Viewer *sharedWidget, bool antialiasing = false);
  ~Viewer();
  bool testDisplayId(double, double, double)override;
  void updateIds(CGAL::Three::Scene_item *)override;
  //! overload several CGAL::QGLViewer virtual functions
  //! Draws the scene.
  void draw()override;
  //!This step happens after draw(). It is here that all the useful information is displayed, like the axis system or the informative text.
  void drawVisualHints()override;
  //! Deprecated. Does the same as draw().
  void fastDraw()override;
  bool isExtensionFound()override;
  void initializeGL() override;
  //! Initializes the OpenGL functions and sets the backGround color.
  void init()override;
  //! Draws the scene "with names" to allow picking.
  void drawWithNames()override;
  /*! Uses the parameter pixel's coordinates to get the corresponding point
   * in the World frame. If this point is found, emits selectedPoint, selected,
   * and selectionRay signals.
   */
  void postSelection(const QPoint&)override;
  //! Sets the picking matrix to allow the picking.
  void beginSelection(const QPoint &point)override;
  //! Sets the pick matrix to Identity once the picking is done.
  void endSelection(const QPoint &point)override;
  //! Sets the scene for the viewer.
  void setScene(CGAL::Three::Scene_draw_interface* scene)override;
  //! @returns the antialiasing state.
  bool antiAliasing() const override;
  //! @returns the fastDrawing state.
  bool inFastDrawing() const override;
  //! Implementation of `Viewer_interface::inDrawWithNames()`
  bool inDrawWithNames() const override;
  //! Implementation of `Viewer_interface::attribBuffers()`
  void attribBuffers(int program_name) const override;
  //! Implementation of `Viewer_interface::getShaderProgram()`
  QOpenGLShaderProgram* getShaderProgram(int name) const override;
  //!Declares a program names `name`, using `v_shader` as vertex shader and `f_shader` as fragment shader.
  QOpenGLShaderProgram* declare_program(int name,
                                        const char* v_shader,
                                        const char* f_shader)const;
  QPainter* getPainter()override;


  TextRenderer* textRenderer() override;
  void enableClippingBox(QVector4D box[6]) override;
  void disableClippingBox() override;
  void set2DSelectionMode(bool) override;
  void setStaticImage(QImage image) override;
  const QImage& staticImage() const override;
  //!Set total number of depth peeling passes.
   void setTotalPass(int);
   void resetFov();
   const QVector3D& scaler() const override;
Q_SIGNALS:
  void sendMessage(QString);
  void doneInitGL(CGAL::Three::Viewer_interface*);
  void socketClosed();
public Q_SLOTS:
  //! Sets the antialiasing to true or false.
  void setAntiAliasing(bool b) override;
  //! If b is true, facets will be ligted from both internal and external sides.
  //! If b is false, only the side that is exposed to the light source will be lighted.
  void setTwoSides(bool b) override;
  void setBackFrontShading(bool b) override;
  void SetOrthoProjection( bool b) override;
  //! If b is true, some items are displayed in a simplified version when moving the camera.
  //! If b is false, items display is never altered, even when moving.
  void setFastDrawing(bool b) override;
  //! Makes the camera turn around.
  void turnCameraBy180Degres() override;
  //! @returns a QString containing the position and orientation of the camera.
  QString dumpCameraCoordinates() override;
  //!Moves the camera to the new coordinates (position and orientation) through an animation.
  bool moveCameraToCoordinates(QString,
                               float animation_duration = 0.5f) override;
  //!Makes the Viewer display a message
  void printMessage(QString message, int ms_delay );
  void displayMessage(const QString &_message, int delay);
  void displayMessage(const QString &_message){displayMessage(_message, 2000);}
  void hideMessage();
  void setBindingSelect() override
  {
    setMouseBinding(::Qt::ShiftModifier, ::Qt::LeftButton, CGAL::qglviewer::SELECT);
  }
  virtual void setNoBinding() override
  {
    setMouseBinding(::Qt::ShiftModifier, ::Qt::LeftButton, CGAL::qglviewer::NO_CLICK_ACTION);
  }

  void setLighting();
  void setBackFrontColors();

  void messageLogged(QOpenGLDebugMessage);
#ifdef CGAL_USE_WEBSOCKETS
  void setShareCam(bool, QString);
  void onSocketConnected();
  void onTextMessageSocketReceived(QString message);
#endif
  void scaleScene();
  void showEntireScene()override;
protected:
  void paintEvent(QPaintEvent *)override;
  void paintGL()override;

  //!Defines the behavior for the mouse press events
  void mousePressEvent(QMouseEvent*)override;
  void mouseDoubleClickEvent(QMouseEvent*)override;
  void wheelEvent(QWheelEvent *)override;
  //!Defines the behavior for the key press events
  void keyPressEvent(QKeyEvent*)override;
  //!Deal with context menu events
  void contextMenuEvent(QContextMenuEvent*)override;
  //!Defines the behavior for the key release events
  void keyReleaseEvent(QKeyEvent *)override;

protected:
  friend class Viewer_impl;
  Viewer_impl* d;
  double prev_radius;
  void doBindings();

public:
  QOpenGLFunctions_4_3_Core* openGL_4_3_functions() override;
  void setCurrentPass(int pass) override;
   void setDepthWriting(bool writing_depth) override;
   void setDepthPeelingFbo(QOpenGLFramebufferObject *fbo) override;
   int currentPass()const override;
   bool isDepthWriting()const override;
   QOpenGLFramebufferObject* depthPeelingFbo()override;
   float total_pass()override;
   const GLfloat& getGlPointSize()const override;
   void setGlPointSize(const GLfloat& p) override;
   void makeCurrent() override;
   QVector4D* clipBox() const override;
   bool isClipping() const override;
}; // end class Viewer


#endif // VIEWER_H
