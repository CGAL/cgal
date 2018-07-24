/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.
 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 version 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0

#ifndef QGLVIEWER_QGLVIEWER_H
#define QGLVIEWER_QGLVIEWER_H
#include <CGAL/export/Qt.h>
#include <CGAL/number_type_config.h>
#include <CGAL/Qt/viewer_actions.h>
#include <CGAL/Qt/vec.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Qt/keyFrameInterpolator.h>
#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/manipulatedCameraFrame.h>

#include <QClipboard>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMap>
#include <QVector>
#include <QTime>
#include <QTimer>
#include <QGLContext>
#include <QOpenGLWidget>
#include <QMouseEvent>


class QTabWidget;
class QImage;
class QOpenGLFramebufferObject;

namespace CGAL{
/*! \brief A versatile 3D OpenGL viewer based on QOpenGLWidget.
\class QGLViewer qglviewer.h QGLViewer/qglviewer.h

It features many classical viewer functionalities, such as a camera trackball,
manipulated objects and much <a
href="../features.html">more</a>. Its main goal is to ease the development of
new 3D applications.

New users should read the <a href="../introduction.html">introduction page</a>
to get familiar with important notions such as sceneRadius(), sceneCenter() and
the world coordinate system. Try the numerous simple <a
href="../examples/index.html">examples</a> to discover the possibilities and
understand how it works.

<h3>Usage</h3>

To use a QGLViewer, derive you viewer class from the QGLViewer and overload its
draw() virtual method. See the <a
href="../examples/simpleViewer.html">simpleViewer example</a> for an
illustration.

An other option is to connect your drawing methods to the signals emitted by the
QGLViewer (Qt's callback mechanism). See the <a
href="../examples/callback.html">callback example</a> for a complete
implementation.

\nosubgrouping */
class CGAL_QT_EXPORT QGLViewer : public QOpenGLWidget, public QOpenGLFunctions {
  Q_OBJECT

public:  
  explicit QGLViewer(QGLContext* context, QWidget *parent = 0,
                     ::Qt::WindowFlags flags = 0);
  explicit QGLViewer(QWidget *parent = 0,
                     ::Qt::WindowFlags flags = 0);

  virtual ~QGLViewer();

  /*! @name Display of visual hints */
  //@{
public:
  /*! Returns \c true if the world axis is drawn by the viewer.

  Set by setAxisIsDrawn() or toggleAxisIsDrawn(). Default value is \c false. */
  bool axisIsDrawn() const { return axisIsDrawn_; }
  /*! Returns \c true if a XY grid is drawn by the viewer.

  Set by setGridIsDrawn() or toggleGridIsDrawn(). Default value is \c false. */
  bool gridIsDrawn() const { return gridIsDrawn_; }
  /*! Returns \c true if the viewer displays the current frame rate (Frames Per
  Second).

  Use QApplication::setFont() to define the display font (see drawText()).

  Set by setFPSIsDisplayed() or toggleFPSIsDisplayed(). Use currentFPS() to get
  the current FPS. Default value is \c false. */
  bool FPSIsDisplayed() const { return FPSIsDisplayed_; }
  /*! Returns \c true if text display (see drawText()) is enabled.

  Set by setTextIsEnabled() or toggleTextIsEnabled(). This feature conveniently
  removes all the possibly displayed text, cleaning display. Default value is \c
  true. */
  bool textIsEnabled() const { return textIsEnabled_; }

  /*! Returns \c true if the camera() is being edited in the viewer.

  Set by setCameraIsEdited() or toggleCameraIsEdited(). Default value is \p
  false.

  The current implementation is limited: the defined camera() paths (see
  qglviewer::Camera::keyFrameInterpolator()) are simply displayed using
  qglviewer::Camera::drawAllPaths(). Actual camera and path edition will be
  implemented in the future. */
  bool cameraIsEdited() const { return cameraIsEdited_; }
  
public Q_SLOTS:
  /*! Sets the state of axisIsDrawn(). Emits the axisIsDrawnChanged() signal.
   * See also toggleAxisIsDrawn(). */
  void setAxisIsDrawn(bool draw = true) {
    if(!draw)
      axis_size = 0;
    axisIsDrawn_ = draw;
    Q_EMIT axisIsDrawnChanged(draw);
    update();
  }
  /*! Sets the state of gridIsDrawn(). Emits the gridIsDrawnChanged() signal.
   * See also toggleGridIsDrawn(). */
  void setGridIsDrawn(bool draw = true) {
    if(!draw)
    {
      grid_size=0;
      g_axis_size=0;
    }
    gridIsDrawn_ = draw;
    Q_EMIT gridIsDrawnChanged(draw);
    update();
  }
  /*! Sets the state of FPSIsDisplayed(). Emits the FPSIsDisplayedChanged()
   * signal. See also toggleFPSIsDisplayed(). */
  void setFPSIsDisplayed(bool display = true) {
    FPSIsDisplayed_ = display;
    Q_EMIT FPSIsDisplayedChanged(display);
    update();
  }
  /*! Sets the state of textIsEnabled(). Emits the textIsEnabledChanged()
   * signal. See also toggleTextIsEnabled(). */
  void setTextIsEnabled(bool enable = true) {
    textIsEnabled_ = enable;
    Q_EMIT textIsEnabledChanged(enable);
    update();
  }
  void setCameraIsEdited(bool edit = true);

  /*! Toggles the state of axisIsDrawn(). See also setAxisIsDrawn(). */
  void toggleAxisIsDrawn() { setAxisIsDrawn(!axisIsDrawn()); }
  /*! Toggles the state of gridIsDrawn(). See also setGridIsDrawn(). */
  void toggleGridIsDrawn() { setGridIsDrawn(!gridIsDrawn()); }
  /*! Toggles the state of FPSIsDisplayed(). See also setFPSIsDisplayed(). */
  void toggleFPSIsDisplayed() { setFPSIsDisplayed(!FPSIsDisplayed()); }
  /*! Toggles the state of textIsEnabled(). See also setTextIsEnabled(). */
  void toggleTextIsEnabled() { setTextIsEnabled(!textIsEnabled()); }
  /*! Toggles the state of cameraIsEdited(). See also setCameraIsEdited(). */
  void toggleCameraIsEdited() { setCameraIsEdited(!cameraIsEdited()); }
  //@}

  /*! @name Viewer's colors */
  //@{
public:
  /*! Returns the background color of the viewer.

  This method is provided for convenience since the background color is an
  OpenGL state variable set with \c glClearColor(). However, this internal
  representation has the advantage that it is saved (resp. restored) with
  saveStateToFile() (resp. restoreStateFromFile()).

  Use setBackgroundColor() to define and activate a background color.

  \attention Each QColor component is an integer ranging from 0 to 255. This
  differs from the qreal values used by \c glClearColor() which are in the
  0.0-1.0 range. Default value is (51, 51, 51) (dark gray). You may have to
  change foregroundColor() accordingly.

  \attention This method does not return the current OpenGL clear color as \c
  glGet() does. Instead, it returns the QGLViewer internal variable. If you
  directly use \c glClearColor() or \c qglClearColor() instead of
  setBackgroundColor(), the two results will differ. */
  QColor backgroundColor() const { return backgroundColor_; }

  /*! Returns the foreground color used by the viewer.

  This color is used when FPSIsDisplayed(), gridIsDrawn(), to display the camera
  paths when the cameraIsEdited().

  \attention Each QColor component is an integer in the range 0-255. This
  differs from the qreal values used by \c glColor3f() which are in the range
  0-1. Default value is (180, 180, 180) (light gray).

  Use \c qglColor(foregroundColor()) to set the current OpenGL color to the
  foregroundColor().

  See also backgroundColor(). */
  QColor foregroundColor() const { return foregroundColor_; }
public Q_SLOTS:
  /*! Sets the backgroundColor() of the viewer and calls \c qglClearColor(). See
     also setForegroundColor(). */
  void setBackgroundColor(const QColor &color) {
    backgroundColor_ = color;
    glClearColor(GLclampf(color.redF()),
                 GLclampf(color.greenF()),
                 GLclampf(color.blueF()),
                 GLclampf(color.alphaF()));
  }
  /*! Sets the foregroundColor() of the viewer, used to draw visual hints. See
   * also setBackgroundColor(). */
  void setForegroundColor(const QColor &color) { foregroundColor_ = color; }
  //@}

  /*! @name Scene dimensions */
  //@{
public:
  /*! Returns the scene radius.

  The entire displayed scene should be included in a sphere of radius
  sceneRadius(), centered on sceneCenter().

  This approximate value is used by the camera() to set
  qglviewer::Camera::zNear() and qglviewer::Camera::zFar(). It is also used to
  showEntireScene() or to scale the world axis display..

  Default value is 1.0. This method is equivalent to camera()->sceneRadius().
  See setSceneRadius(). */
  qreal sceneRadius() const;
  /*! Returns the scene center, defined in world coordinates.

  See sceneRadius() for details.

  Default value is (0,0,0). Simply a wrapper for camera()->sceneCenter(). Set
  using setSceneCenter().

  Do not mismatch this value (that only depends on the scene) with the
  qglviewer::Camera::pivotPoint(). */
  qglviewer::Vec sceneCenter() const;

public Q_SLOTS:
  /*! Sets the sceneRadius().

    The camera() qglviewer::Camera::flySpeed() is set to 1% of this value by
    this method. Simple wrapper around camera()->setSceneRadius(). */
  virtual void setSceneRadius(qreal radius);

  /*! Sets the sceneCenter(), defined in world coordinates.

    \attention The qglviewer::Camera::pivotPoint() is set to the sceneCenter()
    value by this method. */
  virtual void setSceneCenter(const qglviewer::Vec &center);

  /*! Convenient way to call setSceneCenter() and setSceneRadius() from a (world
    axis aligned) bounding box of the scene. Takes the offset into account.

    This is equivalent to:
    \code
    setSceneCenter((min+max) / 2.0);
    setSceneRadius((max-min).norm() / 2.0);
    \endcode */
  void setSceneBoundingBox(const qglviewer::Vec &min,
                           const qglviewer::Vec &max);

  /*! Moves the camera so that the entire scene is visible.

    Simple wrapper around qglviewer::Camera::showEntireScene(). */
  void showEntireScene() ;
  //@}

  /*! @name Associated objects */
  //@{
public:
  /*! Returns the associated qglviewer::Camera, never \c NULL. */
  qglviewer::Camera *camera() const { return camera_; }

  /*! Returns the viewer's qglviewer::ManipulatedFrame.

  This qglviewer::ManipulatedFrame can be moved with the mouse when the
  associated mouse bindings are used (default is when pressing the \c Control
  key with any mouse button). Use setMouseBinding() to define new bindings.

  See the <a href="../examples/manipulatedFrame.html">manipulatedFrame
  example</a> for a complete implementation.

  Default value is \c NULL, meaning that no qglviewer::ManipulatedFrame is set.
*/
  qglviewer::ManipulatedFrame *manipulatedFrame() const {
    return manipulatedFrame_;
  }

public Q_SLOTS:
  void setCamera(qglviewer::Camera *const camera);
  void setManipulatedFrame(qglviewer::ManipulatedFrame *frame);
  //@}

  /*! @name Mouse grabbers */
  //@{
public:
  /*! Returns the current qglviewer::MouseGrabber, or \c NULL if no
  qglviewer::MouseGrabber currently grabs mouse events.

  When qglviewer::MouseGrabber::grabsMouse(), the different mouse events are
  sent to the mouseGrabber() instead of their usual targets (camera() or
  manipulatedFrame()).

  See the qglviewer::MouseGrabber documentation for details on MouseGrabber's
  mode of operation.

  In order to use MouseGrabbers, you need to enable mouse tracking (so that
  mouseMoveEvent() is called even when no mouse button is pressed). Add this
  line in init() or in your viewer constructor: \code setMouseTracking(true);
  \endcode
  Note that mouse tracking is disabled by default. Use
  QWidget::hasMouseTracking() to retrieve current state. */
  qglviewer::MouseGrabber *mouseGrabber() const { return mouseGrabber_; }

  void
  setMouseGrabberIsEnabled(const qglviewer::MouseGrabber *const mouseGrabber,
                           bool enabled = true);
  /*! Returns \c true if \p mouseGrabber is enabled.

  Default value is \c true for all MouseGrabbers. When set to \c false using
  setMouseGrabberIsEnabled(), the specified \p mouseGrabber will never become
  the mouseGrabber() of this QGLViewer. This is useful when you use several
  viewers: some MouseGrabbers may only have a meaning for some specific viewers
  and should not be selectable in others.

  You can also use qglviewer::MouseGrabber::removeFromMouseGrabberPool() to
  completely disable a MouseGrabber in all the QGLViewers. */
  bool
  mouseGrabberIsEnabled(const qglviewer::MouseGrabber *const mouseGrabber) {
    return !disabledMouseGrabbers_.contains(
        reinterpret_cast<size_t>(mouseGrabber));
  }
public Q_SLOTS:
  void setMouseGrabber(qglviewer::MouseGrabber *mouseGrabber);
  //@}

  /*! @name State of the viewer */
  //@{
public:
  /*! Returns the aspect ratio of the viewer's widget (width() / height()). */
  qreal aspectRatio() const { return width() / static_cast<qreal>(height()); }
  /*! Returns the current averaged viewer frame rate.

  This value is computed and averaged over 20 successive frames. It only changes
  every 20 draw() (previously computed value is otherwise returned).

  This method is useful for true real-time applications that may adapt their
  computational load accordingly in order to maintain a given frequency.

  This value is meaningful only when draw() is regularly called, either using a
  \c QTimer, when animationIsStarted() or when the camera is manipulated with
  the mouse.  */
  qreal currentFPS() { return f_p_s_; }
  /*!
   * Returns the string used to display the current fps.
   */
  QString fpsString() { return fpsString_; }
  /*! Returns \c true if the viewer is in fullScreen mode.

  Default value is \c false. Set by setFullScreen() or toggleFullScreen().

  Note that if the QGLViewer is embedded in an other QWidget, it returns \c true
  when the top level widget is in full screen mode. */
  bool isFullScreen() const { return fullScreen_; }
  
  /*! Returns the recommended size for the QGLViewer. Default value is 600x400
   * pixels. */
  virtual QSize sizeHint() const { return QSize(600, 400); }
  /*!
   * Sets the offset of the scene. The offset is the difference between the origin 
   * of the world and the origin of the scene. It is relevant when the whole scene is translated
   * of a big number, because there is a useless loss of precision when drawing.
   * 
   * The offset must be added to the drawn coordinates, and substracted from the computation 
   * \attention  the result of pointUnderPixel is the real item translated by the offset.
   * 
   */
  void setOffset(qglviewer::Vec offset);
  
  /*!
   * returns the offset of the scene.
   * \see `setOffset()`
   */
  qglviewer::Vec offset()const;

public Q_SLOTS:
  void setFullScreen(bool fullScreen = true);
  /*! Toggles the state of isFullScreen(). See also setFullScreen(). */
  void toggleFullScreen() { setFullScreen(!isFullScreen()); }
  void toggleCameraMode();

protected:
  bool cameraIsInRotateMode() const;
  //@}

  /*! @name Display methods */
  //@{
public:
  void drawArrow(double r, double R, int prec,
                        qglviewer::Vec from, qglviewer::Vec to, qglviewer::Vec color, std::vector<float> &data);
  void drawAxis(qreal l = 1.0);
  void drawGrid(qreal size= 1.0, int nbSubdivisions = 10);

  virtual void startScreenCoordinatesSystem(bool upward = false) const;
  virtual void stopScreenCoordinatesSystem() const;

  void drawText(int x, int y, const QString &text, const QFont &fnt = QFont());
  void displayMessage(const QString &message, int delay = 2000);

protected:
  virtual void drawLight(GLenum light, qreal scale = 1.0) const;

protected:
  void displayFPS();
  

//@}

#ifdef DOXYGEN
  /*! @name Useful inherited methods */
  //@{
public:
  /*! Returns viewer's widget width (in pixels). See QOpenGLWidget
   * documentation. */
  int width() const;
  /*! Returns viewer's widget height (in pixels). See QOpenGLWidget
   * documentation. */
  int height() const;
  /*! Updates the display. Do not call draw() directly, use this method instead.
   * See QOpenGLWidget documentation. */
  virtual void update();
  /*! Returns \c true if the widget has a valid GL rendering context. See
  QOpenGLWidget documentation. */
  bool isValid() const;
  /*! Makes this widget's rendering context the current OpenGL rendering
  context. Useful with several viewers. See QOpenGLWidget documentation. */
  virtual void makeCurrent();
  /*! Returns \c true if mouseMoveEvent() is called even when no mouse button is
  pressed.

  You need to setMouseTracking() to \c true in order to use MouseGrabber (see
  mouseGrabber()). See details in the QWidget documentation. */
  bool hasMouseTracking() const;
public Q_SLOTS:
  /*! Resizes the widget to size \p width by \p height pixels. See also width()
   * and height(). */
  virtual void resize(int width, int height);
  /*! Sets the hasMouseTracking() value. */
  virtual void setMouseTracking(bool enable);
#endif

  /*! @name Buffer to texture */
  //@{
public:
  GLuint bufferTextureId() const;
  /*! Returns the texture coordinate corresponding to the u extremum of the
  bufferTexture.

  The bufferTexture is created by copyBufferToTexture(). The texture size has
  powers of two dimensions and the buffer image hence only fills a part of it.
  This value corresponds to the u coordinate of the extremum right side of the
  buffer image.

  Use (0,0) to (bufferTextureMaxU(), bufferTextureMaxV()) texture coordinates to
  map the entire texture on a quad. */
  qreal bufferTextureMaxU() const { return bufferTextureMaxU_; }
  /*! Same as bufferTextureMaxU(), but for the v texture coordinate. */
  qreal bufferTextureMaxV() const { return bufferTextureMaxV_; }
#if (QT_VERSION >= QT_VERSION_CHECK(5, 4, 0))
  // These methods are part of the QGLWidget public API.
  // As of version 2.7.0, the use of QOpenGLWidget instead means that they have
  // to be provided for backward compatibility.
  void renderText(int x, int y, const QString &str,
                  const QFont &font = QFont());
  void renderText(double x, double y, double z, const QString &str,
                  const QFont &font = QFont());
#endif

public Q_SLOTS:
  void copyBufferToTexture(GLint, GLenum = GL_NONE);
  //@}

  /*! @name Animation */
  //@{
public:
  /*! Return \c true when the animation loop is started.

  During animation, an infinite loop calls animate() and draw() and then waits
  for animationPeriod() milliseconds before calling animate() and draw() again.
  And again.

  Use startAnimation(), stopAnimation() or toggleAnimation() to change this
  value.

  See the <a href="../examples/animation.html">animation example</a> for
  illustration. */
  bool animationIsStarted() const { return animationStarted_; }
  /*! The animation loop period, in milliseconds.

  When animationIsStarted(), this is delay waited after draw() to call animate()
  and draw() again. Default value is 40 milliseconds (25 Hz).

  This value will define the currentFPS() when animationIsStarted() (provided
  that your animate() and draw() methods are fast enough).

  If you want to know the maximum possible frame rate of your machine on a given
  scene, setAnimationPeriod() to \c 0, and startAnimation() (keyboard shortcut
  is \c Enter). The display will then be updated as often as possible, and the
  frame rate will be meaningful.

  \note This value is taken into account only the next time you call
  startAnimation(). If animationIsStarted(), you should stopAnimation() first.
*/
  int animationPeriod() const { return animationPeriod_; }

public Q_SLOTS:
  /*! Sets the animationPeriod(), in milliseconds. */
  void setAnimationPeriod(int period) { animationPeriod_ = period; }
  virtual void startAnimation();
  virtual void stopAnimation();
  /*! Scene animation method.

    When animationIsStarted(), this method is in charge of the scene update
    before each draw(). Overload it to define how your scene evolves over time.
    The time should either be regularly incremented in this method (frame-rate
    independent animation) or computed from actual time (for instance using
    QTime::elapsed()) for real-time animations.

        Note that KeyFrameInterpolator (which regularly updates a Frame) does
    not use this method to animate a Frame, but rather rely on a QTimer
    signal-slot mechanism.

    See the <a href="../examples/animation.html">animation example</a> for an
    illustration. */
  virtual void animate() { Q_EMIT animateNeeded(); }
  /*! Calls startAnimation() or stopAnimation(), depending on
   * animationIsStarted(). */
  void toggleAnimation() {
    if (animationIsStarted())
      stopAnimation();
    else
      startAnimation();
  }
  //@}
public:
  /*!
   * Prompt a configuration dialog and takes a snapshot.
   */
  void saveSnapshot();

public:
Q_SIGNALS:
  /*! Signal emitted by the default init() method.

  Connect this signal to the methods that need to be called to initialize your
  viewer or overload init(). */
  void viewerInitialized();

  /*! Signal emitted by the default draw() method.

  Connect this signal to your main drawing method or overload draw(). See the <a
  href="../examples/callback.html">callback example</a> for an illustration. */
  void drawNeeded();

  /*! Signal emitted at the end of the QGLViewer::paintGL() method, when frame
  is drawn.

  Can be used to notify an image grabbing process that the image is ready.  */
  void drawFinished(bool automatic);

  /*! Signal emitted by the default animate() method.

  Connect this signal to your scene animation method or overload animate(). */
  void animateNeeded();

  /*! Signal emitted by the default QGLViewer::help() method.

  Connect this signal to your own help method or overload help(). */
  void helpRequired();

  /*! This signal is emitted whenever axisIsDrawn() changes value. */
  void axisIsDrawnChanged(bool drawn);
  /*! This signal is emitted whenever gridIsDrawn() changes value. */
  void gridIsDrawnChanged(bool drawn);
  /*! This signal is emitted whenever FPSIsDisplayed() changes value. */
  void FPSIsDisplayedChanged(bool displayed);
  /*! This signal is emitted whenever textIsEnabled() changes value. */
  void textIsEnabledChanged(bool enabled);
  /*! This signal is emitted whenever cameraIsEdited() changes value.. */
  void cameraIsEditedChanged(bool edited);
  /*! Signal emitted by select().

  Connect this signal to your selection method or overload select(), or more
  probably simply drawWithNames(). */
  void pointSelected(const QMouseEvent *e);

  /*! Signal emitted by setMouseGrabber() when the mouseGrabber() is changed.

  \p mouseGrabber is a pointer to the new MouseGrabber. Note that this signal is
  emitted with a \c NULL parameter each time a MouseGrabber stops grabbing
  mouse. */
  void mouseGrabberChanged(qglviewer::MouseGrabber *mouseGrabber);

  /*! @name Help window */
  //@{
public:
  /*! Returns the QString displayed in the help() window main tab.

  Overload this method to define your own help string, which should shortly
  describe your application and explain how it works. Rich-text (HTML) tags can
  be used (see QStyleSheet() documentation for available tags): \code QString
  myViewer::helpString() const
  {
  QString text("<h2>M y V i e w e r</h2>");
  text += "Displays a <b>Scene</b> using OpenGL. Move the camera using the
  mouse."; return text;
  }
  \endcode

  See also mouseString() and keyboardString(). */
  virtual QString helpString() const { return tr("No help available."); }

  virtual QString mouseString() const;
  virtual QString keyboardString() const;

public Q_SLOTS:
  virtual void help();
  virtual void aboutQGLViewer();

protected:
  /*! Returns a pointer to the help widget.

  Use this only if you want to directly modify the help widget. Otherwise use
  helpString(), setKeyDescription() and setMouseBindingDescription() to
  customize the text displayed in the help window tabs. */
  QTabWidget *helpWidget() { return helpWidget_; }
  //@}

  /*! @name Drawing methods */
  //@{
protected:
  virtual void resizeGL(int width, int height);
  virtual void initializeGL();

  /*! Initializes the viewer OpenGL context.

  This method is called before the first drawing and should be overloaded to
  initialize some of the OpenGL flags. The default implementation is empty. See
  initializeGL().

  Typical usage include camera() initialization (showEntireScene()), previous
  viewer state restoration (restoreStateFromFile()), OpenGL state modification
  and display list creation.

  Note that initializeGL() modifies the standard OpenGL context. These values
  can be restored back in this method.

  \attention You should not call updateGL() (or any method that calls it) in
  this method, as it will result in an infinite loop. The different QGLViewer
  set methods (setAxisIsDrawn(), setFPSIsDisplayed()...) are protected against
  this problem and can safely be called.

  \note All the OpenGL specific initializations must be done in this method: the
  OpenGL context is not yet available in your viewer constructor. */
  virtual void init() { Q_EMIT viewerInitialized(); }

  virtual void paintGL();
  virtual void preDraw();

  /*! The core method of the viewer, that draws the scene.

  If you build a class that inherits from QGLViewer, this is the method you want
  to overload. See the <a href="../examples/simpleViewer.html">simpleViewer
  example</a> for an illustration.

  The camera modelView matrix set in preDraw() converts from the world to the
  camera coordinate systems. Vertices given in draw() can then be considered as
  being given in the world coordinate system. The camera is moved in this world
  using the mouse. This representation is much more intuitive than the default
  camera-centric OpenGL standard.

  \attention The \c GL_PROJECTION matrix should not be modified by this method,
  to correctly display visual hints (axis, grid, FPS...) in postDraw(). Use
  push/pop or call camera()->loadProjectionMatrix() at the end of draw() if you
  need to change the projection matrix (unlikely). On the other hand, the \c
  GL_MODELVIEW matrix can be modified and left in a arbitrary state. */
  virtual void draw() {}
  virtual void fastDraw();
  virtual void postDraw();
  //@}

  /*! @name Mouse, keyboard and event handlers */
  //@{
protected:
  virtual void mousePressEvent(QMouseEvent *);
  virtual void mouseMoveEvent(QMouseEvent *);
  virtual void mouseReleaseEvent(QMouseEvent *);
  virtual void mouseDoubleClickEvent(QMouseEvent *);
  virtual void wheelEvent(QWheelEvent *);
  virtual void keyPressEvent(QKeyEvent *);
  virtual void keyReleaseEvent(QKeyEvent *);
  virtual void timerEvent(QTimerEvent *);
  virtual void closeEvent(QCloseEvent *);
  //@}

  /*! @name Object selection */
  //@{
public:
  /*! Returns the name (an integer value) of the entity that was last selected
  by select(). This value is set by endSelection(). See the select()
  documentation for details.

  As a convention, this method returns -1 if the selectBuffer() was empty,
  meaning that no object was selected.

  Return value is -1 before the first call to select(). This value is modified
  using setSelectedName(). */
  int selectedName() const { return selectedObjectId_; }
  /*! Returns the selectBuffer() size.

  See the select() documentation for details. Use setSelectBufferSize() to
  change this value.

  Default value is 4000 (i.e. 1000 objects in selection region, since each
  object pushes 4 values). This size should be over estimated to prevent a
  buffer overflow when many objects are drawn under the mouse cursor. */
  int selectBufferSize() const { return selectBufferSize_; }

  /*! Returns the width (in pixels) of a selection frustum, centered on the
  mouse cursor, that is used to select objects.

  The height of the selection frustum is defined by selectRegionHeight().

  The objects that will be drawn in this region by drawWithNames() will be
  recorded in the selectBuffer(). endSelection() then analyzes this buffer and
  setSelectedName() to the name of the closest object. See the gluPickMatrix()
  documentation for details.

  The default value is 3, which is adapted to standard applications. A smaller
  value results in a more precise selection but the user has to be careful for
  small feature selection.

  See the <a href="../examples/multiSelect.html">multiSelect example</a> for an
  illustration. */
  int selectRegionWidth() const { return selectRegionWidth_; }
  /*! See the selectRegionWidth() documentation. Default value is 3 pixels. */
  int selectRegionHeight() const { return selectRegionHeight_; }

  /*! Returns a pointer to an array of \c GLuint.

  This buffer is used by the \c GL_SELECT mode in select() to perform object
  selection. The buffer size can be modified using setSelectBufferSize(). If you
  overload endSelection(), you will analyze the content of this buffer. See the
  \c glSelectBuffer() man page for details. */
  GLuint *selectBuffer() { return selectBuffer_; }

public Q_SLOTS:
  virtual void select(const QMouseEvent *event);
  virtual void select(const QPoint &point);

  void setSelectBufferSize(int size);
  /*! Sets the selectRegionWidth(). */
  void setSelectRegionWidth(int width) { selectRegionWidth_ = width; }
  /*! Sets the selectRegionHeight(). */
  void setSelectRegionHeight(int height) { selectRegionHeight_ = height; }
  /*! Set the selectedName() value.

    Used in endSelection() during a selection. You should only call this method
    if you overload the endSelection() method. */
  void setSelectedName(int id) { selectedObjectId_ = id; }

protected:
  virtual void beginSelection(const QPoint &point);
  /*! This method is called by select() and should draw selectable entities.

  Default implementation is empty. Overload and draw the different elements of
your scene you want to be able to select. The default select() implementation
relies on the \c GL_SELECT, and requires that each selectable element is drawn
within a \c glPushName() - \c glPopName() block. A typical usage would be (see
the <a href="../examples/select.html">select example</a>): \code void
Viewer::drawWithNames() { for (int i=0; i<nbObjects; ++i) { glPushName(i);
    object(i)->draw();
    glPopName();
   }
}
\endcode

  The resulting selected name is computed by endSelection(), which
setSelectedName() to the integer id pushed by this method (a value of -1 means
no selection). Use selectedName() to update your selection, probably in the
postSelection() method.

  \attention If your selected objects are points, do not use \c
glBegin(GL_POINTS); and \c glVertex3fv() in the above \c draw() method (not
compatible with raster mode): use \c glRasterPos3fv() instead. */
  virtual void drawWithNames() {}
  virtual void endSelection(const QPoint &point);
  /*! This method is called at the end of the select() procedure. It should
  finalize the selection process and update the data
  structure/interface/computation/display... according to the newly selected
  entity.

  The default implementation is empty. Overload this method if needed, and use
  selectedName() to retrieve the selected entity name (returns -1 if no object
  was selected). See the <a href="../examples/select.html">select example</a>
  for an illustration. */
  virtual void postSelection(const QPoint &point) { Q_UNUSED(point); }
  //@}

  /*! @name Keyboard customization */
  //@{
public:
  unsigned int shortcut(qglviewer::KeyboardAction action) const;

  ::Qt::Key pathKey(unsigned int index) const;
  ::Qt::KeyboardModifiers addKeyFrameKeyboardModifiers() const;
  ::Qt::KeyboardModifiers playPathKeyboardModifiers() const;

public Q_SLOTS:
  void setShortcut(qglviewer::KeyboardAction action, unsigned int key);

  void setKeyDescription(unsigned int key, QString description);
  void clearShortcuts();

// Key Frames shortcut keys

  virtual void setPathKey(int key, unsigned int index = 0);
  virtual void setPlayPathKeyboardModifiers(::Qt::KeyboardModifiers modifiers);
  virtual void setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers modifiers);
  //@}

public:
  /*! @name Mouse customization */
  //@{

  qglviewer::MouseAction mouseAction(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                          ::Qt::MouseButton button) const;
  int mouseHandler(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                   ::Qt::MouseButton button) const;

  void getMouseActionBinding(qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                             bool withConstraint, ::Qt::Key &key,
                             ::Qt::KeyboardModifiers &modifiers,
                             ::Qt::MouseButton &button) const;

  qglviewer::ClickAction clickAction(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                          ::Qt::MouseButton button, bool doubleClick = false,
                          ::Qt::MouseButtons buttonsBefore = ::Qt::NoButton) const;

  void getClickActionBinding(qglviewer::ClickAction action, ::Qt::Key &key,
                             ::Qt::KeyboardModifiers &modifiers,
                             ::Qt::MouseButton &button, bool &doubleClick,
                             ::Qt::MouseButtons &buttonsBefore) const;

  qglviewer::MouseAction wheelAction(::Qt::Key key, ::Qt::KeyboardModifiers modifiers) const;
  int wheelHandler(::Qt::Key key, ::Qt::KeyboardModifiers modifiers) const;

  void getWheelActionBinding(qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                             bool withConstraint, ::Qt::Key &key,
                             ::Qt::KeyboardModifiers &modifiers) const;

public Q_SLOTS:

  void setMouseBinding(::Qt::KeyboardModifiers modifiers, ::Qt::MouseButton buttons,
                       qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                       bool withConstraint = true);
  void setMouseBinding(::Qt::KeyboardModifiers modifiers, ::Qt::MouseButton button,
                       qglviewer::ClickAction action, bool doubleClick = false,
                       ::Qt::MouseButtons buttonsBefore = ::Qt::NoButton);
  void setWheelBinding(::Qt::KeyboardModifiers modifiers, qglviewer::MouseHandler handler,
                       qglviewer::MouseAction action, bool withConstraint = true);
  void
  setMouseBindingDescription(::Qt::KeyboardModifiers modifiers,
                             ::Qt::MouseButton button, QString description,
                             bool doubleClick = false,
                             ::Qt::MouseButtons buttonsBefore = ::Qt::NoButton);

  void setMouseBinding(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                       ::Qt::MouseButton buttons, qglviewer::MouseHandler handler,
                       qglviewer::MouseAction action, bool withConstraint = true);
  void setMouseBinding(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                       ::Qt::MouseButton button, qglviewer::ClickAction action,
                       bool doubleClick = false,
                       ::Qt::MouseButtons buttonsBefore = ::Qt::NoButton);
  void setWheelBinding(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                       qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                       bool withConstraint = true);
  void
  setMouseBindingDescription(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                             ::Qt::MouseButton button, QString description,
                             bool doubleClick = false,
                             ::Qt::MouseButtons buttonsBefore = ::Qt::NoButton);

  void clearMouseBindings();

protected:
  static QString mouseActionString(qglviewer::MouseAction ma);
  static QString clickActionString(qglviewer::ClickAction ca);
  //@}

  /*! @name State persistence */
  //@{
public:
  QString stateFileName() const;
  virtual QDomElement domElement(const QString &name,
                                 QDomDocument &document) const;
Q_SIGNALS:
  void needNewContext();

public Q_SLOTS:
  virtual void initFromDOMElement(const QDomElement &element);
  virtual void saveStateToFile(); // cannot be const because of QMessageBox
  virtual bool restoreStateFromFile();

  /*! Defines the stateFileName() used by saveStateToFile() and
    restoreStateFromFile().

    The file name can have an optional prefix directory (no prefix meaning
    current directory). If the directory does not exist, it will be created by
    saveStateToFile().

    \code
    // Name depends on the displayed 3D model. Saved in current directory.
    setStateFileName(3DModelName() + ".xml");

    // Files are stored in a dedicated directory under user's home directory.
    setStateFileName(QDir::homeDirPath + "/.config/myApp.xml");
    \endcode */
  void setStateFileName(const QString &name) { stateFileName_ = name; }


protected:
  static void saveStateToFileForAllViewers();
  //@}

  /*! @name QGLViewer pool */
  //@{
public:
  /*! Returns a \c QList that contains pointers to all the created QGLViewers.
    Note that this list may contain \c NULL pointers if the associated viewer
  has been deleted.

  Can be useful to apply a method or to connect a signal to all the viewers:
    \code
  foreach (QGLViewer* viewer, QGLViewer::QGLViewerPool())
    connect(myObject, SIGNAL(IHaveChangedSignal()), viewer, SLOT(update()));
  \endcode
*/
  static QList<QGLViewer *> &QGLViewerPool();

  /*! Returns the index of the QGLViewer \p viewer in the QGLViewerPool(). This
  index in unique and can be used to identify the different created QGLViewers
  (see stateFileName() for an application example).

  When a QGLViewer is deleted, the QGLViewers' indexes are preserved and NULL is
  set for that index. When a QGLViewer is created, it is placed in the first
  available position in that list. Returns -1 if the QGLViewer could not be
  found (which should not be possible). */
  static int QGLViewerIndex(const QGLViewer *const viewer) {
    return QGLViewer::QGLViewerPool().indexOf(const_cast<QGLViewer *>(viewer));
  }
//@}

#ifndef DOXYGEN
  /*! @name Visual hints */
  //@{
public:
  virtual void setVisualHintsMask(int mask, int delay = 2000);
  virtual void drawVisualHints();
  QOpenGLFramebufferObject* getStoredFrameBuffer();

public Q_SLOTS:
  virtual void resetVisualHints();
//@}
#endif

private Q_SLOTS:
  // Patch for a Qt bug with fullScreen on startup
  void delayedFullScreen() {
    move(prevPos_);
    setFullScreen();
  }
  void hideMessage();

private:
  // Copy constructor and operator= are declared private and undefined
  // Prevents everyone from trying to use them
  QGLViewer(const QGLViewer &v);
  QGLViewer &operator=(const QGLViewer &v);

protected:
  // Set parameters to their default values. Called by the constructors.
  void defaultConstructor();

  void handleKeyboardAction(qglviewer::KeyboardAction id);

  // C a m e r a
  qglviewer::Camera *camera_;
  bool cameraIsEdited_;
  qreal previousCameraZClippingCoefficient_;
  unsigned int previousPathId_; // double key press recognition
  void connectAllCameraKFIInterpolatedSignals(bool connection = true);

  // C o l o r s
  QColor backgroundColor_, foregroundColor_;

  // D i s p l a y    f l a g s
  bool axisIsDrawn_;    // world axis
  bool gridIsDrawn_;    // world XY grid
  bool FPSIsDisplayed_; // Frame Per Seconds
  bool textIsEnabled_;  // drawText() actually draws text or not
  bool fullScreen_;     // full screen mode
  QPoint prevPos_;      // Previous window position, used for full screen mode

  // A n i m a t i o n
  bool animationStarted_; // animation mode started
  int animationPeriod_;   // period in msecs
  int animationTimerId_;

  // F P S    d i s p l a y
  QTime fpsTime_;
  unsigned int fpsCounter_;
  QString fpsString_;
  qreal f_p_s_;

  // M e s s a g e s
  QString message_;
  bool displayMessage_;
  QTimer messageTimer_;

  // M a n i p u l a t e d    f r a m e
  qglviewer::ManipulatedFrame *manipulatedFrame_;
  bool manipulatedFrameIsACamera_;

  // M o u s e   G r a b b e r
  qglviewer::MouseGrabber *mouseGrabber_;
  bool mouseGrabberIsAManipulatedFrame_;
  bool mouseGrabberIsAManipulatedCameraFrame_;
  QMap<size_t, bool> disabledMouseGrabbers_;

  // S e l e c t i o n
  int selectRegionWidth_, selectRegionHeight_;
  int selectBufferSize_;
  GLuint *selectBuffer_;
  int selectedObjectId_;

  // V i s u a l   h i n t s
  int visualHint_;

  // S h o r t c u t   k e y s
  void setDefaultShortcuts();
  QString cameraPathKeysString() const;
  QMap<qglviewer::KeyboardAction, QString> keyboardActionDescription_;
  QMap<qglviewer::KeyboardAction, unsigned int> keyboardBinding_;
  QMap<unsigned int, QString> keyDescription_;

  // K e y   F r a m e s   s h o r t c u t s
  QMap< ::Qt::Key, unsigned int> pathIndex_;
  ::Qt::KeyboardModifiers addKeyFrameKeyboardModifiers_,
      playPathKeyboardModifiers_;

  // B u f f e r   T e x t u r e
  GLuint bufferTextureId_;
  qreal bufferTextureMaxU_, bufferTextureMaxV_;
  int bufferTextureWidth_, bufferTextureHeight_;
  unsigned int previousBufferTextureFormat_;
  int previousBufferTextureInternalFormat_;

#ifndef DOXYGEN
  // M o u s e   a c t i o n s
  struct MouseActionPrivate {
    qglviewer::MouseHandler handler;
    qglviewer::MouseAction action;
    bool withConstraint;
  };

  // M o u s e   b i n d i n g s
  struct MouseBindingPrivate {
    const ::Qt::KeyboardModifiers modifiers;
    const ::Qt::MouseButton button;
    const ::Qt::Key key;

    MouseBindingPrivate(::Qt::KeyboardModifiers m, ::Qt::MouseButton b, ::Qt::Key k)
        : modifiers(m), button(b), key(k) {}

    // This sort order is used in mouseString() to display sorted mouse bindings
    bool operator<(const MouseBindingPrivate &mbp) const {
      if (key != mbp.key)
        return key < mbp.key;
      if (modifiers != mbp.modifiers)
        return modifiers < mbp.modifiers;
      return button < mbp.button;
    }
  };

  // W h e e l   b i n d i n g s
  struct WheelBindingPrivate {
    const ::Qt::KeyboardModifiers modifiers;
    const ::Qt::Key key;

    WheelBindingPrivate(::Qt::KeyboardModifiers m, ::Qt::Key k)
        : modifiers(m), key(k) {}

    // This sort order is used in mouseString() to display sorted wheel bindings
    bool operator<(const WheelBindingPrivate &wbp) const {
      if (key != wbp.key)
        return key < wbp.key;
      return modifiers < wbp.modifiers;
    }
  };

  // C l i c k   b i n d i n g s
  struct ClickBindingPrivate {
    const ::Qt::KeyboardModifiers modifiers;
    const ::Qt::MouseButton button;
    const bool doubleClick;
    const ::Qt::MouseButtons
        buttonsBefore; // only defined when doubleClick is true
    const ::Qt::Key key;

    ClickBindingPrivate(::Qt::KeyboardModifiers m, ::Qt::MouseButton b, bool dc,
                        ::Qt::MouseButtons bb, ::Qt::Key k)
        : modifiers(m), button(b), doubleClick(dc), buttonsBefore(bb), key(k) {}

    // This sort order is used in mouseString() to display sorted mouse bindings
    bool operator<(const ClickBindingPrivate &cbp) const {
      if (key != cbp.key)
        return key < cbp.key;
      if (buttonsBefore != cbp.buttonsBefore)
        return buttonsBefore < cbp.buttonsBefore;
      if (modifiers != cbp.modifiers)
        return modifiers < cbp.modifiers;
      if (button != cbp.button)
        return button < cbp.button;
      return doubleClick != cbp.doubleClick;
    }
  };
#endif
  static QString formatClickActionPrivate(ClickBindingPrivate cbp);
  static bool isValidShortcutKey(int key);

  QMap<ClickBindingPrivate, QString> mouseDescription_;

  void setDefaultMouseBindings();
  void performClickAction(qglviewer::ClickAction ca, const QMouseEvent *const e);
  QMap<MouseBindingPrivate, MouseActionPrivate> mouseBinding_;
  QMap<WheelBindingPrivate, MouseActionPrivate> wheelBinding_;
  QMap<ClickBindingPrivate, qglviewer::ClickAction> clickBinding_;
  ::Qt::Key currentlyPressedKey_;
  
  // S t a t e   F i l e
  QString stateFileName_;

  // H e l p   w i n d o w
  QTabWidget *helpWidget_;
  
  //internal drawing buffers
  enum VBO
  {
    Grid = 0,
    Grid_axis,
    Axis,
    Pivot_point,
    VBO_size
  };
  enum VAO
  {
    GRID = 0,
    GRID_AXIS,
    AXIS,
    PIVOT_POINT,
    VAO_size
  };
  QOpenGLShaderProgram rendering_program;
  QOpenGLShaderProgram rendering_program_light;
  QOpenGLVertexArrayObject vaos[VAO_size];
  QVector<QOpenGLBuffer> vbos;
  std::size_t grid_size;
  std::size_t g_axis_size;
  std::size_t axis_size;
  QOpenGLFramebufferObject* stored_fbo;
  //S n a p s h o t
  QImage* takeSnapshot(qglviewer::SnapShotBackground  background_color, 
                       QSize finalSize, double oversampling, bool expand);
  
  //Internal Projection Matrix
  
  // O f f s e t
  qglviewer::Vec _offset;
  //C o n t e x t
  bool is_ogl_4_3;
public:
  //! Is used to know if the openGL context is 4.3 or ES 2.0.
  //! @returns `true` if the context is 4.3.
  //! @returns `false` if the context is ES 2.0.  
  bool isOpenGL_4_3()const {return is_ogl_4_3; }
  
};

} //end CGAL


#ifdef CGAL_HEADER_ONLY
#include <CGAL/Qt/qglviewer_impl_list.h>
#endif // CGAL_HEADER_ONLY

#endif // QGLVIEWER_QGLVIEWER_H
