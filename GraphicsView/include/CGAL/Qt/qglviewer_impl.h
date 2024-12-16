/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-only

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <cmath>

#include <CGAL/Qt/qglviewer.h>
#include <CGAL/Qt/manipulatedCameraFrame.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Qt/keyFrameInterpolator.h>
#include <CGAL/Qt/image_interface.h>

#include <QApplication>
#include <QDateTime>
#include <QDir>
#include <QFileInfo>
#include <QOpenGLContext>
#include <QImage>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPainter>
#include <QPushButton>
#include <QTabWidget>
#include <QTextEdit>
#include <QTextStream>
#include <QTimer>
#include <QUrl>
#include <QtAlgorithms>
#include <QColorDialog>
#include <QOpenGLFramebufferObject>
#include <QFileDialog>
#include <QElapsedTimer>

namespace CGAL{
// Static private variable
CGAL_INLINE_FUNCTION
QList<CGAL::QGLViewer *> &CGAL::QGLViewer::QGLViewerPool() {
  static QList<CGAL::QGLViewer *> QGLViewerPool_;
  void* p = qApp->property("qglviewer pool").value<void*>();
  if(p == 0) {
    p = (void*)(&QGLViewerPool_);
    qApp->setProperty("qglviewer pool", QVariant::fromValue(p));
  }
  return *static_cast<QList<CGAL::QGLViewer *> * >(p);
}

/*! \mainpage

libCGAL::QGLViewer is a free C++ library based on Qt that enables the quick creation
of OpenGL 3D viewers. It features a powerful camera trackball and simple
applications simply require an implementation of the <code>draw()</code> method.
This makes it a tool of choice for OpenGL beginners and assignments. It provides
mouse manipulated frames, interpolated
keyFrames, object selection, and much more. It is fully
customizable and easy to extend to create complex applications, with a possible
Qt GUI.

libCGAL::QGLViewer is <i>not</i> a 3D viewer that can be used directly to view 3D
scenes in various formats. It is more likely to be the starting point for the
coding of such a viewer.

libCGAL::QGLViewer is based on the Qt toolkit and hence compiles on any architecture
(Unix-Linux, Mac, Windows, ...). Full reference documentation and many examples
are provided.

See the project main page for details on the project and installation steps. */

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::defaultConstructor() {
  setFocusPolicy(::Qt::StrongFocus);

  CGAL::QGLViewer::QGLViewerPool().append(this);
  camera_ = new qglviewer::Camera(this);
  setCamera(camera_);

  setDefaultShortcuts();
  setDefaultMouseBindings();

  fpsTime_.start();
  fpsCounter_ = 0;
  f_p_s_ = 0.0;
  fpsString_ = tr("%1Hz", "Frames per seconds, in Hertz").arg("?");
  visualHint_ = 0;
  previousPathId_ = 0;
  // prevPos_ is not initialized since pos() is not meaningful here.
  // It will be set when setFullScreen(false) is called after
  // setFullScreen(true)

  manipulatedFrame_ = nullptr;
  manipulatedFrameIsACamera_ = false;
  mouseGrabberIsAManipulatedFrame_ = false;
  mouseGrabberIsAManipulatedCameraFrame_ = false;
  displayMessage_ = false;
  connect(&messageTimer_, SIGNAL(timeout()), SLOT(hideMessage()));
  messageTimer_.setSingleShot(true);
  helpWidget_ = nullptr;
  setMouseGrabber(nullptr);

  setSceneRadius(1.0);
  showEntireScene();

  setAxisIsDrawn(false);
  setGridIsDrawn(false);
  setFPSIsDisplayed(false);
  setCameraIsEdited(false);
  setTextIsEnabled(true);
  // Make sure move() is not called, which would call initializeGL()
  fullScreen_ = false;
  setFullScreen(false);

  animationTimerId_ = 0;
  stopAnimation();
  setAnimationPeriod(40); // 25Hz

  selectBuffer_ = nullptr;
  setSelectBufferSize(4 * 1000);
  setSelectRegionWidth(3);
  setSelectRegionHeight(3);
  setSelectedName(-1);

  bufferTextureId_ = 0;
  bufferTextureMaxU_ = 0.0;
  bufferTextureMaxV_ = 0.0;
  bufferTextureWidth_ = 0;
  bufferTextureHeight_ = 0;
  previousBufferTextureFormat_ = 0;
  previousBufferTextureInternalFormat_ = 0;
  currentlyPressedKey_ = ::Qt::Key(0);

  setAttribute(::Qt::WA_NoSystemBackground);
  axisIsDrawn_ = true;

  _offset = CGAL::qglviewer::Vec(0,0,0);
  stored_fbo = nullptr;
  is_sharing = false;
  is_linked = false;
  shared_context = nullptr;
  _first_tick  = true;
}

CGAL_INLINE_FUNCTION
CGAL::QGLViewer::QGLViewer(QWidget *parent,
                     ::Qt::WindowFlags flags)
    : QOpenGLWidget(parent, flags) {
  defaultConstructor();
}
CGAL_INLINE_FUNCTION
CGAL::QGLViewer::QGLViewer(QOpenGLContext* context, QWidget *parent,
                   ::Qt::WindowFlags flags)
  : QOpenGLWidget(parent, flags) {
  defaultConstructor();
  shared_context = context;
  is_sharing = true;
}

/*! Virtual destructor.

The viewer is replaced by \c nullptr in the QGLViewerPool() (in order to preserve
other viewer's indexes) and allocated memory is released. The camera() is
deleted and should be copied before if it is shared by an other viewer. */
CGAL_INLINE_FUNCTION
CGAL::QGLViewer::~QGLViewer() {
  CGAL::QGLViewer::QGLViewerPool().removeAll(this);

  camera()->deleteLater();
  delete[] selectBuffer_;
  if (helpWidget()) {
    // Needed for Qt 4 which has no main widget.
    helpWidget()->close();
    delete helpWidget_;
  }
  disconnect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &CGAL::QGLViewer::contextIsDestroyed);
}


/*! Initializes the CGAL::QGLViewer OpenGL context and then calls user-defined init().

This method is automatically called once, before the first call to paintGL().

Overload init() instead of this method to modify viewer specific OpenGL state.

If a 4.3 context could not be set, an ES 2.0 context will be used instead.
 \see `isOpenGL_4_3()`
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::initializeGL() {
  if(!is_sharing)
  {
    QSurfaceFormat format = context()->format();
    context()->format().setOption(QSurfaceFormat::DebugContext);
    if ( !context()->isValid()
         || format.majorVersion() != 4
         || QCoreApplication::arguments().contains(QStringLiteral("--old")))

    {
      is_ogl_4_3 = false;
    }
    else
    {
      is_ogl_4_3 = true;
    }

    QSurfaceFormat cur_f = QOpenGLContext::currentContext()->format();
    const char* rt =(cur_f.renderableType() == QSurfaceFormat::OpenGLES) ? "GLES" : "GL";
    qDebug().noquote() <<tr("Using OpenGL context %1.%2 %3").arg(cur_f.majorVersion()).arg(cur_f.minorVersion()).arg(rt);
  }
  else
  {
    context()->setFormat(shared_context->format());
    context()->setShareContext(shared_context);
    context()->create();
    makeCurrent();
  }
  connect(context(), &QOpenGLContext::aboutToBeDestroyed,
          this, &CGAL::QGLViewer::contextIsDestroyed);
  QOpenGLFunctions::initializeOpenGLFunctions();
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  // Default colors
  setForegroundColor(QColor(180, 180, 180));
  setBackgroundColor(QColor(51, 51, 51));

  // Clear the buffer where we're going to draw
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // Calls user defined method. Default emits a signal.
  init();

  //OpenGL buffers and programs initialization
  for(int i=0; i<VAO_size; ++i)
  {
    vaos[i].create();
  }
  vbos.resize(VBO_size);
  for(int i=0; i<VBO_size; ++i)
  {
    vbos[i].create();
  }
  if(!is_linked)
  {
    //Vertex source code
    const char vertex_source[] =
    {
      "#version 150\n"
      "in vec4 vertex;\n"
      "uniform mat4 mvp_matrix;\n"
      "void main(void)\n"
      "{\n"
      "  gl_Position = mvp_matrix * vertex;\n"
      "}"
    };
    const char vertex_source_comp[] =
    {
      "attribute highp vec4 vertex;\n"
      "uniform highp mat4 mvp_matrix;\n"
      "void main(void)\n"
      "{\n"
      "  gl_Position = mvp_matrix * vertex;\n"
      "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
      "#version 150\n"
      "uniform vec4 color;\n"
      "out vec4 out_color;\n"
      "void main(void)\n"
      "{\n"
      "  out_color = color;\n"
      "}"
    };
    const char fragment_source_comp[] =
    {
      "uniform highp vec4 color;\n"
      "void main(void)\n"
      "{\n"
      "  gl_FragColor = color;\n"
      "}"
    };

    //It is said in the doc that a QOpenGLShader is
    // only destroyed with the QOpenGLShaderProgram
    //it has been linked with.

    QOpenGLShader vertex_shader(QOpenGLShader::Vertex);
    QOpenGLShader fragment_shader(QOpenGLShader::Fragment);
    if(is_ogl_4_3)
    {
      if(!vertex_shader.compileSourceCode(vertex_source))
      {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
      }

      if(!fragment_shader.compileSourceCode(fragment_source))
      {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
      }
    }
    else
    {
      if(!vertex_shader.compileSourceCode(vertex_source_comp))
      {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
      }

      if(!fragment_shader.compileSourceCode(fragment_source_comp))
      {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
      }
    }

    if(!rendering_program.addShader(&vertex_shader))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(&fragment_shader))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
      qDebug() << rendering_program.log();
    }
  }

  is_linked = true;

  // Give time to glInit to finish and then call setFullScreen().
  if (isFullScreen())
    QTimer::singleShot(100, this, SLOT(delayedFullScreen()));
}

/*! Main paint method, inherited from \c QOpenGLWidget.

Calls the following methods, in that order:
\arg preDraw() : places the
camera in the world coordinate system. \arg draw() (or fastDraw() when the
camera is manipulated) : main drawing method. Should be overloaded. \arg
postDraw() : display of visual hints (world axis, FPS...) */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::paintGL() {
  // Clears screen, set model view matrix...
  preDraw();
  // Used defined method. Default calls draw()
  if (camera()->frame()->isManipulated())
    fastDraw();
  else
    draw();
  // Add visual hints: axis, camera, grid...
  postDraw();
  Q_EMIT drawFinished(true);
}

/*! Sets OpenGL state before draw().

Default behavior clears screen and sets the projection and modelView matrices:
\code
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

camera()->loadProjectionMatrix();
camera()->loadModelViewMatrix();
\endcode

Emits the drawNeeded() signal once this is done (see the <a
href="../examples/callback.html">callback example</a>). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::preDraw() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // GL_PROJECTION matrix
  camera()->loadProjectionMatrix();
  // GL_MODELVIEW matrix
  camera()->loadModelViewMatrix();

  Q_EMIT drawNeeded();
}

/*! Called after draw() to draw viewer visual hints.

Default implementation displays axis, grid, FPS... when the respective flags are
sets.

See the <a href="../examples/multiSelect.html">multiSelect</a> and <a
href="../examples/contribs.html#thumbnail">thumbnail</a> examples for an
overloading illustration.

The GLContext (color, LIGHTING, BLEND...) is \e not modified by this method, so
that in draw(), the user can rely on the OpenGL context he defined. Respect this
convention (by pushing/popping the different attributes) if you overload this
method. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::postDraw() {
  // Pivot point, line when camera rolls, zoom region
  if (gridIsDrawn()) {
    drawGrid(camera()->sceneRadius());
  }
  if (axisIsDrawn()) {
    drawAxis(1.0);
  }


  drawVisualHints();
  // FPS computation
  const unsigned int maxCounter = 20;
  if (++fpsCounter_ == maxCounter) {
    f_p_s_ = 1000.0 * maxCounter / fpsTime_.restart();
    fpsString_ = tr("%1Hz", "Frames per seconds, in Hertz")
                     .arg(f_p_s_, 0, 'f', ((f_p_s_ < 10.0) ? 1 : 0));
    fpsCounter_ = 0;
  }


  glDisable(GL_DEPTH_TEST);

  if (FPSIsDisplayed())
    displayFPS();
  if (displayMessage_)
    drawText(10, height() - 10, message_);

  //zoom region
  if(camera()->frame()->action_ == qglviewer::ZOOM_ON_REGION)
  {
    QPainter painter(this);
    painter.setPen(QColor(120,120,120));
    painter.drawRect(QRect(camera()->frame()->pressPos_, mapFromGlobal(QCursor::pos())));
    painter.end();
  }

  //zoom_fov indicator
  if(camera()->frame()->action_ == qglviewer::ZOOM_FOV)
  {
    QPainter painter(this);
    QPoint bot(width()-30,height()/2-0.33*height()),
           top(width()-30, height()/2+0.33*height());
    int fov_height = (top.y()-bot.y())*camera()->fieldOfView()*4.0/CGAL_PI + bot.y();


    painter.setPen(QColor(120,120,120));
    painter.drawLine(bot, top);
    painter.fillRect(QRect(QPoint(width()-40, fov_height+10),
                           QPoint(width()-20, fov_height-10)),
                     QColor(120,120,120));
    painter.end();
    camera()->frame()->action_= qglviewer::NO_MOUSE_ACTION;
  }
  glEnable(GL_DEPTH_TEST);
}


/*! Draws a simplified version of the scene to guarantee interactive camera
displacements.

This method is called instead of draw() when the CGAL::qglviewer::Camera::frame() is
CGAL::qglviewer::ManipulatedCameraFrame::isManipulated(). Default implementation
simply calls draw().

Overload this method if your scene is too complex to allow for interactive
camera manipulation. See the <a href="../examples/fastDraw.html">fastDraw
example</a> for an illustration. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::fastDraw() { draw(); }

/*! Starts (\p edit = \c true, default) or stops (\p edit=\c false) the edition
of the camera().

Current implementation is limited to paths display. Get current state using
cameraIsEdited().

\attention This method sets the CGAL::qglviewer::Camera::zClippingCoefficient() to 5.0
when \p edit is \c true, so that the Camera paths (see
CGAL::qglviewer::Camera::keyFrameInterpolator()) are not clipped. It restores the
previous value when \p edit is \c false. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setCameraIsEdited(bool edit) {
  cameraIsEdited_ = edit;
  if (edit) {
    previousCameraZClippingCoefficient_ = camera()->zClippingCoefficient();

    camera()->setZClippingCoefficient(5.0);
  } else
    camera()->setZClippingCoefficient(previousCameraZClippingCoefficient_);

  Q_EMIT cameraIsEditedChanged(edit);

  update();
}

// Key bindings. 0 means not defined
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setDefaultShortcuts() {
  // D e f a u l t   a c c e l e r a t o r s
  setShortcut(qglviewer::DRAW_AXIS, ::Qt::Key_A);
  setShortcut(qglviewer::DRAW_GRID, ::Qt::Key_G);
  setShortcut(qglviewer::DISPLAY_FPS, ::Qt::Key_F);
  setShortcut(qglviewer::ENABLE_TEXT, ::Qt::SHIFT, ::Qt::Key_Question);
  setShortcut(qglviewer::EXIT_VIEWER, ::Qt::Key_Escape);
  setShortcut(qglviewer::CAMERA_MODE, ::Qt::Key_Space);
  setShortcut(qglviewer::FULL_SCREEN, ::Qt::ALT, ::Qt::Key_Return);
  setShortcut(qglviewer::ANIMATION, ::Qt::Key_Return);
  setShortcut(qglviewer::HELP, ::Qt::Key_H);
  setShortcut(qglviewer::MOVE_CAMERA_LEFT, ::Qt::Key_Left);
  setShortcut(qglviewer::MOVE_CAMERA_RIGHT, ::Qt::Key_Right);
  setShortcut(qglviewer::MOVE_CAMERA_UP, ::Qt::Key_Up);
  setShortcut(qglviewer::MOVE_CAMERA_DOWN, ::Qt::Key_Down);
  setShortcut(qglviewer::INCREASE_FLYSPEED, ::Qt::Key_Plus);
  setShortcut(qglviewer::DECREASE_FLYSPEED, ::Qt::Key_Minus);

  keyboardActionDescription_[qglviewer::DISPLAY_FPS] =
      tr("Toggles the display of the FPS", "DISPLAY_FPS action description");
  keyboardActionDescription_[qglviewer::FULL_SCREEN] =
      tr("Toggles full screen display", "FULL_SCREEN action description");
  keyboardActionDescription_[qglviewer::DRAW_AXIS] = tr(
      "Toggles the display of the world axis", "DRAW_AXIS action description");
  keyboardActionDescription_[qglviewer::DRAW_GRID] =
      tr("Toggles the display of the XY grid", "DRAW_GRID action description");
  keyboardActionDescription_[qglviewer::CAMERA_MODE] = tr(
      "Changes camera mode (observe or fly)", "CAMERA_MODE action description");
  keyboardActionDescription_[qglviewer::HELP] =
      tr("Opens this help window", "HELP action description");
  keyboardActionDescription_[qglviewer::ANIMATION] =
      tr("Starts/stops the animation", "ANIMATION action description");
  keyboardActionDescription_[qglviewer::ENABLE_TEXT] =
      tr("Toggles the display of the text", "ENABLE_TEXT action description");
  keyboardActionDescription_[qglviewer::EXIT_VIEWER] =
      tr("Exits program", "EXIT_VIEWER action description");
  keyboardActionDescription_[qglviewer::MOVE_CAMERA_LEFT] =
      tr("Moves camera left", "MOVE_CAMERA_LEFT action description");
  keyboardActionDescription_[qglviewer::MOVE_CAMERA_RIGHT] =
      tr("Moves camera right", "MOVE_CAMERA_RIGHT action description");
  keyboardActionDescription_[qglviewer::MOVE_CAMERA_UP] =
      tr("Moves camera up", "MOVE_CAMERA_UP action description");
  keyboardActionDescription_[qglviewer::MOVE_CAMERA_DOWN] =
      tr("Moves camera down", "MOVE_CAMERA_DOWN action description");
  keyboardActionDescription_[qglviewer::INCREASE_FLYSPEED] =
      tr("Increases fly speed", "INCREASE_FLYSPEED action description");
  keyboardActionDescription_[qglviewer::DECREASE_FLYSPEED] =
      tr("Decreases fly speed", "DECREASE_FLYSPEED action description");

  // K e y f r a m e s   s h o r t c u t   k e y s
  setPathKey(::Qt::Key_F1, 1);
  setPathKey(::Qt::Key_F2, 2);
  setPathKey(::Qt::Key_F3, 3);
  setPathKey(::Qt::Key_F4, 4);
  setPathKey(::Qt::Key_F5, 5);
  setPathKey(::Qt::Key_F6, 6);
  setPathKey(::Qt::Key_F7, 7);
  setPathKey(::Qt::Key_F8, 8);
  setPathKey(::Qt::Key_F9, 9);
  setPathKey(::Qt::Key_F10, 10);
  setPathKey(::Qt::Key_F11, 11);
  setPathKey(::Qt::Key_F12, 12);

  setAddKeyFrameKeyboardModifiers(::Qt::AltModifier);
  setPlayPathKeyboardModifiers(::Qt::NoModifier);
}

// M o u s e   b e h a v i o r
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setDefaultMouseBindings() {
  const ::Qt::KeyboardModifiers cameraKeyboardModifiers = ::Qt::NoModifier;
  const ::Qt::KeyboardModifiers frameKeyboardModifiers = ::Qt::ControlModifier;

  //#CONNECTION# toggleCameraMode()
  for (int handler = 0; handler < 2; ++handler) {
    qglviewer::MouseHandler mh = (qglviewer::MouseHandler)(handler);
    ::Qt::KeyboardModifiers modifiers =
        (mh == qglviewer::FRAME) ? frameKeyboardModifiers : cameraKeyboardModifiers;

    setMouseBinding(modifiers, ::Qt::LeftButton, mh, qglviewer::ROTATE);
    setMouseBinding(modifiers, ::Qt::MiddleButton, mh, qglviewer::ZOOM);
    setMouseBinding(modifiers, ::Qt::RightButton, mh, qglviewer::TRANSLATE);

    setMouseBinding(::Qt::Key_R, modifiers, ::Qt::LeftButton, mh, qglviewer::SCREEN_ROTATE);

    setWheelBinding(modifiers, mh, qglviewer::ZOOM);
  }
  setWheelBinding(::Qt::Key_Z, ::Qt::NoModifier, qglviewer::CAMERA, qglviewer::ZOOM_FOV);

  // Z o o m   o n   r e g i o n
  setMouseBinding(::Qt::ShiftModifier, ::Qt::MiddleButton, qglviewer::CAMERA, qglviewer::ZOOM_ON_REGION);

  // S e l e c t
  setMouseBinding(::Qt::ShiftModifier, ::Qt::LeftButton, qglviewer::SELECT);

  setMouseBinding(::Qt::ShiftModifier, ::Qt::RightButton, qglviewer::RAP_FROM_PIXEL);
  // D o u b l e   c l i c k
  setMouseBinding(::Qt::NoModifier, ::Qt::LeftButton, qglviewer::ALIGN_CAMERA, true);
  setMouseBinding(::Qt::NoModifier, ::Qt::MiddleButton, qglviewer::SHOW_ENTIRE_SCENE, true);
  setMouseBinding(::Qt::NoModifier, ::Qt::RightButton, qglviewer::CENTER_SCENE, true);

  setMouseBinding(frameKeyboardModifiers, ::Qt::LeftButton, qglviewer::ALIGN_FRAME, true);
  // middle double click makes no sense for manipulated frame
  setMouseBinding(frameKeyboardModifiers, ::Qt::RightButton, qglviewer::CENTER_FRAME, true);

  // A c t i o n s   w i t h   k e y   m o d i f i e r s
  setMouseBinding(::Qt::Key_Z, ::Qt::NoModifier, ::Qt::LeftButton, qglviewer::ZOOM_ON_PIXEL);
  setMouseBinding(::Qt::Key_Z, ::Qt::NoModifier, ::Qt::RightButton, qglviewer::ZOOM_TO_FIT);


#ifdef Q_OS_MAC
  // Specific Mac bindings for touchpads. Two fingers emulate a wheelEvent which
  // zooms. There is no right button available : make Option key + left emulate
  // the right button. A Control+Left indeed emulates a right click (OS X system
  // configuration), but it does no seem to support dragging. Done at the end to
  // override previous settings.
  const ::Qt::KeyboardModifiers macKeyboardModifiers = ::Qt::AltModifier;

  setMouseBinding(macKeyboardModifiers, ::Qt::LeftButton, qglviewer::CAMERA, qglviewer::TRANSLATE);
  setMouseBinding(macKeyboardModifiers, ::Qt::LeftButton, qglviewer::CENTER_SCENE, true);
  setMouseBinding(frameKeyboardModifiers | macKeyboardModifiers, ::Qt::LeftButton,
                  qglviewer::CENTER_FRAME, true);
  setMouseBinding(frameKeyboardModifiers | macKeyboardModifiers, ::Qt::LeftButton,
                  qglviewer::FRAME, qglviewer::TRANSLATE);
#endif
}

/*! Associates a new CGAL::qglviewer::Camera to the viewer.

You should only use this method when you derive a new class from
CGAL::qglviewer::Camera and want to use one of its instances instead of the original
class.

It you simply want to save and restore Camera positions, use
CGAL::qglviewer::Camera::addKeyFrameToPath() and CGAL::qglviewer::Camera::playPath()
instead.

This method silently ignores \c nullptr \p camera pointers. The calling method is
responsible for deleting the previous camera pointer in order to prevent memory
leaks if needed.

The sceneRadius() and sceneCenter() of \p camera are set to the \e current
CGAL::QGLViewer values.

All the \p camera CGAL::qglviewer::Camera::keyFrameInterpolator()
CGAL::qglviewer::KeyFrameInterpolator::interpolated() signals are connected to the
viewer update() slot. The connections with the previous viewer's camera are
removed. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setCamera(qglviewer::Camera *const camera) {
  if (!camera)
    return;

  camera->setSceneRadius(sceneRadius());
  camera->setSceneCenter(sceneCenter());
  camera->setScreenWidthAndHeight(width(), height(), devicePixelRatio());

  // Disconnect current camera from this viewer.
  disconnect(this->camera()->frame(), SIGNAL(manipulated()), this,
             SLOT(update()));
  disconnect(this->camera()->frame(), SIGNAL(spun()), this, SLOT(update()));

  // Connect camera frame to this viewer.
  connect(camera->frame(), SIGNAL(manipulated()), SLOT(update()));
  connect(camera->frame(), SIGNAL(spun()), SLOT(update()));

  connectAllCameraKFIInterpolatedSignals(false);
  camera_ = camera;
  connectAllCameraKFIInterpolatedSignals();

  previousCameraZClippingCoefficient_ = this->camera()->zClippingCoefficient();
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::connectAllCameraKFIInterpolatedSignals(bool connection) {
  for (QMap<unsigned int, qglviewer::KeyFrameInterpolator *>::ConstIterator
           it = camera()->kfi_.begin(),
           end = camera()->kfi_.end();
       it != end; ++it) {
    if (connection)
      connect(camera()->keyFrameInterpolator(it.key()), SIGNAL(interpolated()),
              SLOT(update()));
    else
      disconnect(camera()->keyFrameInterpolator(it.key()),
                 SIGNAL(interpolated()), this, SLOT(update()));
  }

  if (connection)
    connect(camera()->interpolationKfi_, SIGNAL(interpolated()),
            SLOT(update()));
  else
    disconnect(camera()->interpolationKfi_, SIGNAL(interpolated()), this,
               SLOT(update()));
}

/*! Draws a representation of \p light.

Called in draw(), this method is useful to debug or display your light setup.
Light drawing depends on the type of light (point, spot, directional).

The method retrieves the light setup using \c glGetLightfv. Position and define
your lights before calling this method.

Light is drawn using its diffuse color. Disabled lights are not displayed.

Drawing size is proportional to sceneRadius(). Use \p scale to rescale it.

See the <a href="../examples/drawLight.html">drawLight example</a> for an
illustration.

\attention You need to enable \c GL_COLOR_MATERIAL before calling this method.
\c glColor is set to the light diffuse color. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::drawLight(GLenum, qreal ) const {
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::renderText(int x, int y, const QString &str,
                           const QFont &font) {
  QColor fontColor = QColor(0, 0,
                            0, 255);

  // Render text
  QPainter painter(this);
  painter.setPen(fontColor);
  painter.setFont(font);
  painter.drawText(x, y, str);
  painter.end();
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::renderText(double x, double y, double z, const QString &str,
                           const QFont &font) {
  using CGAL::qglviewer::Vec;
  const Vec proj = camera_->projectedCoordinatesOf(Vec(x, y, z));
  renderText(int(proj.x), int(proj.y), str, font);
}

/*! Draws \p text at position \p x, \p y (expressed in screen coordinates
pixels, origin in the upper left corner of the widget).

The default QApplication::font() is used to render the text when no \p fnt is
specified. Use QApplication::setFont() to define this default font.

You should disable \c GL_LIGHTING and \c GL_DEPTH_TEST before this method so
that colors are properly rendered.

This method can be used in conjunction with the
CGAL::qglviewer::Camera::projectedCoordinatesOf() method to display a text attached to
an object. In your draw() method use: \code CGAL::qglviewer::Vec screenPos =
camera()->projectedCoordinatesOf(myFrame.position());
drawText((int)screenPos[0], (int)screenPos[1], "My Object");
\endcode
See the <a href="../examples/screenCoordSystem.html">screenCoordSystem
example</a> for an illustration.

Text is displayed only when textIsEnabled() (default). This mechanism allows the
user to conveniently remove all the displayed text with a single keyboard
shortcut.

See also displayMessage() to drawText() for only a short amount of time.

Use renderText(x,y,z, text) instead if you want to draw a text located
 at a specific 3D position instead of 2D screen coordinates (fixed size text,
facing the camera).

The \c GL_MODELVIEW and \c GL_PROJECTION matrices are not modified by this
method.
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::drawText(int x, int y, const QString &text, const QFont &fnt) {
  if (!textIsEnabled())
    return;

  renderText(x, y, text, fnt);
}

/*! Briefly displays a message in the lower left corner of the widget.
Convenient to provide feedback to the user.

\p message is displayed during \p delay milliseconds (default is 2 seconds)
using drawText().

This method should not be called in draw(). If you want to display a text in
each draw(), use drawText() instead.

If this method is called when a message is already displayed, the new message
replaces the old one. Use setTextIsEnabled() (default shortcut is '?') to enable
or disable text (and hence messages) display. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::displayMessage(const QString &message, int delay) {
  message_ = message;
  displayMessage_ = true;
  // Was set to single shot in defaultConstructor.
  messageTimer_.start(delay);
  if (textIsEnabled())
    update();
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::hideMessage() {
  displayMessage_ = false;
  if (textIsEnabled())
    update();
}

/*! Displays the averaged currentFPS() frame rate in the upper left corner of
the widget.

update() should be called in a loop in order to have a meaningful value (this is
the case when you continuously move the camera using the mouse or when
animationIsStarted()). setAnimationPeriod(0) to make this loop as fast as
possible in order to reach and measure the maximum available frame rate.

When FPSIsDisplayed() is \c true (default is \c false), this method is called by
postDraw() to display the currentFPS(). Use QApplication::setFont() to define
the font (see drawText()). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::displayFPS() {
  drawText(10,
           int(1.5 * ((QApplication::font().pixelSize() > 0)
                          ? QApplication::font().pixelSize()
                          : QApplication::font().pointSize())),
           fpsString_);
}

/*! Modify the projection matrix so that drawing can be done directly with 2D
screen coordinates.

Once called, the \p x and \p y coordinates passed to \c glVertex are expressed
in pixels screen coordinates. The origin (0,0) is in the upper left corner of
the widget by default. This follows the Qt standards, so that you can directly
use the \c pos() provided by for instance \c QMouseEvent. Set \p upward to \c
true to place the origin in the \e lower left corner, thus following the OpenGL
and mathematical standards. It is always possible to switch between the two
representations using \c newY = height() - \c y.

You need to call stopScreenCoordinatesSystem() at the end of the drawing block
to restore the previous camera matrix.

In practice, this method should be used in draw(). It sets an appropriate
orthographic projection matrix and then sets \c glMatrixMode to \c GL_MODELVIEW.

See the <a href="../examples/screenCoordSystem.html">screenCoordSystem</a>, <a
href="../examples/multiSelect.html">multiSelect</a> and <a
href="../examples/contribs.html#backgroundImage">backgroundImage</a> examples
for an illustration.

You may want to disable \c GL_LIGHTING, to enable \c GL_LINE_SMOOTH or \c
GL_BLEND to draw when this method is used.

If you want to link 2D drawings to 3D objects, use
CGAL::qglviewer::Camera::projectedCoordinatesOf() to compute the 2D projection on
screen of a 3D point (see the <a
href="../examples/screenCoordSystem.html">screenCoordSystem</a> example). See
also drawText().

In this mode, you should use z values that are in the [0.0, 1.0[ range (0.0
corresponding to the near clipping plane and 1.0 being just beyond the far
clipping plane). This interval matches the values that can be read from the
z-buffer. Note that if you use the convenient \c glVertex2i() to provide
coordinates, the implicit 0.0 z coordinate will make your drawings appear \e on
\e top of the rest of the scene. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::startScreenCoordinatesSystem(bool ) const {
}

/*! Stops the pixel coordinate drawing block started by
startScreenCoordinatesSystem().

The \c GL_MODELVIEW and \c GL_PROJECTION matrices modified in
startScreenCoordinatesSystem() are restored. \c glMatrixMode is set to \c
GL_MODELVIEW. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::stopScreenCoordinatesSystem() const {
}

/*! Overloading of the \c QObject method.

If animationIsStarted(), calls animate() and draw(). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::timerEvent(QTimerEvent *) {
  if (animationIsStarted()) {
    animate();
    update();
  }
}

/*! Starts the animation loop. See animationIsStarted(). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::startAnimation() {
  animationTimerId_ = startTimer(animationPeriod());
  animationStarted_ = true;
}

/*! Stops animation. See animationIsStarted(). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::stopAnimation() {
  animationStarted_ = false;
  if (animationTimerId_ != 0)
    killTimer(animationTimerId_);
}


/*! Simple wrapper method: calls \c select(event->pos()).

Emits \c pointSelected(e) which is useful only if you rely on the Qt signal-slot
mechanism and you did not overload CGAL::QGLViewer. If you choose to derive your own
viewer class, simply overload select() (or probably simply drawWithNames(), see
the <a href="../examples/select.html">select example</a>) to implement your
selection mechanism.

This method is called when you use the CGAL::QGLViewer::SELECT mouse binding(s)
(default is Shift + left button). Use setMouseBinding() to change this. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::select(const QMouseEvent *event) {
  // For those who don't derive but rather rely on the signal-slot mechanism.
  Q_EMIT pointSelected(event);
  select(event->pos());
}

/*! This method performs a selection in the scene from pixel coordinates.

It is called when the user clicks on the CGAL::QGLViewer::SELECT
CGAL::QGLViewer::ClickAction binded button(s) (default is Shift + LeftButton).

This template method successively calls four other methods:
\code
beginSelection(point);
drawWithNames();
endSelection(point);
postSelection(point);
\endcode

The default implementation of these methods is as follows (see the methods'
documentation for more details):

\arg beginSelection() sets the \c GL_SELECT mode with the appropriate picking
matrices. A rectangular frustum (of size defined by selectRegionWidth() and
selectRegionHeight()) centered on \p point is created.

\arg drawWithNames() is empty and should be overloaded. It draws each selectable
object of the scene, enclosed by calls to \c glPushName() / \c glPopName() to
tag the object with an integer id.

\arg endSelection() then restores \c GL_RENDER mode and analyzes the
selectBuffer() to set in selectedName() the id of the object that was drawn in
the region. If several object are in the region, the closest one in the depth
buffer is chosen. If no object has been drawn under cursor, selectedName() is
set to -1.

\arg postSelection() is empty and can be overloaded for possible
signal/display/interface update.

See the \c glSelectBuffer() man page for details on this \c GL_SELECT mechanism.

This default implementation is quite limited: only the closer object is
selected, and only one level of names can be pushed. However, this reveals
sufficient in many cases and you usually only have to overload drawWithNames()
to implement a simple object selection process. See the <a
href="../examples/select.html">select example</a> for an illustration.

If you need a more complex selection process (such as a point, edge or triangle
selection, which is easier with a 2 or 3 levels selectBuffer() heap, and which
requires a finer depth sorting to privilege point over edge and edges over
triangles), overload the endSelection() method. Use setSelectRegionWidth(),
setSelectRegionHeight() and setSelectBufferSize() to tune the select buffer
configuration. See the <a href="../examples/multiSelect.html">multiSelect
example</a> for an illustration.

\p point is the center pixel (origin in the upper left corner) of the selection
region. Use CGAL::qglviewer::Camera::convertClickToLine() to transform these
coordinates in a 3D ray if you want to perform an analytical intersection.

\attention \c GL_SELECT mode seems to report wrong results when used in
conjunction with backface culling. If you encounter problems try to \c
glDisable(GL_CULL_FACE). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::select(const QPoint &point) {
  beginSelection(point);
  drawWithNames();
  endSelection(point);
  postSelection(point);
}

/*! This method should prepare the selection. It is called by select() before
drawWithNames().
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::beginSelection(const QPoint &point)
{
  makeCurrent();
  glEnable(GL_SCISSOR_TEST);
  glScissor(point.x() * devicePixelRatio(),
            (camera()->screenHeight() - point.y()) * devicePixelRatio() - 1, 1,
            1);
}

/*! This method is called by select() after scene elements were drawn by
drawWithNames().
 It clears the OpenGL state set by beginSelection*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::endSelection(const QPoint &point) {
  Q_UNUSED(point);
  glDisable(GL_SCISSOR_TEST);
}

/*! Sets the selectBufferSize().

The previous selectBuffer() is deleted and a new one is created. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setSelectBufferSize(int size) {
  if (selectBuffer_)
    delete[] selectBuffer_;
  selectBufferSize_ = size;
  selectBuffer_ = new GLuint[selectBufferSize()];
}

static QString mouseButtonsString(::Qt::MouseButtons b) {
  QString result("");
  bool addAmpersand = false;
  if (b & ::Qt::LeftButton) {
    result += CGAL::QGLViewer::tr("Left", "left mouse button");
    addAmpersand = true;
  }
  if (b & ::Qt::MiddleButton) {
    if (addAmpersand)
      result += " & ";
    result += CGAL::QGLViewer::tr("Middle", "middle mouse button");
    addAmpersand = true;
  }
  if (b & ::Qt::RightButton) {
    if (addAmpersand)
      result += " & ";
    result += CGAL::QGLViewer::tr("Right", "right mouse button");
  }
  return result;
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::performClickAction(qglviewer::ClickAction ca, const QMouseEvent *const e) {
  // the following call is needed to update the pixel ratio
  camera()->setScreenWidthAndHeight(this->width(), this->height(), this->devicePixelRatio());

  // Note: action that need it should call update().
  switch (ca) {
  // # CONNECTION setMouseBinding prevents adding NO_CLICK_ACTION in
  // clickBinding_ This case should hence not be possible. Prevents unused case
  // warning.
  case qglviewer::NO_CLICK_ACTION:
    break;
  case qglviewer::ZOOM_ON_PIXEL:
    camera()->interpolateToZoomOnPixel(e->pos());
    break;
  case qglviewer::ZOOM_TO_FIT:
    camera()->interpolateToFitScene();
    break;
  case qglviewer::SELECT:
    select(e);
    update();
    break;
  case qglviewer::RAP_FROM_PIXEL:
    if (!camera()->setPivotPointFromPixel(e->pos()))
      camera()->setPivotPoint(sceneCenter());
    setVisualHintsMask(1);
    update();
    break;
  case qglviewer::RAP_IS_CENTER:
    camera()->setPivotPoint(sceneCenter());
    setVisualHintsMask(1);
    update();
    break;
  case qglviewer::CENTER_FRAME:
    if (manipulatedFrame())
      manipulatedFrame()->projectOnLine(camera()->position(),
                                        camera()->viewDirection());
    break;
  case qglviewer::CENTER_SCENE:
    camera()->centerScene();
    break;
  case qglviewer::SHOW_ENTIRE_SCENE:
    camera()->showEntireScene();
    break;
  case qglviewer::ALIGN_FRAME:
    if (manipulatedFrame())
      manipulatedFrame()->alignWithFrame(camera()->frame());
    break;
  case qglviewer::ALIGN_CAMERA:
    qglviewer::Frame *frame = new qglviewer::Frame();
    frame->setTranslation(camera()->pivotPoint());
    camera()->frame()->alignWithFrame(frame, true);
    delete frame;
    break;
  }
}

/*! Overloading of the \c QWidget method.

When the user clicks on the mouse:
\arg if a mouseGrabber() is defined, CGAL::qglviewer::MouseGrabber::mousePressEvent()
is called, \arg otherwise, the camera() or the manipulatedFrame() interprets the
mouse displacements, depending on mouse bindings.

Mouse bindings customization can be achieved using setMouseBinding() and
setWheelBinding(). See the <a href="../mouse.html">mouse page</a> for a complete
description of mouse bindings.

See the mouseMoveEvent() documentation for an example of more complex mouse
behavior customization using overloading.

\note When the mouseGrabber() is a manipulatedFrame(), the modifier keys are not
taken into account. This allows for a direct manipulation of the
manipulatedFrame() when the mouse hovers, which is probably what is expected. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::mousePressEvent(QMouseEvent *e) {
  //#CONNECTION# mouseDoubleClickEvent has the same structure
  //#CONNECTION# mouseString() concatenates bindings description in inverse
  // order.
  makeCurrent();
  ClickBindingPrivate cbp(e->modifiers(), e->button(), false,
                          (::Qt::MouseButtons)(e->buttons() & ~(e->button())),
                          currentlyPressedKey_);

  if (clickBinding_.contains(cbp)) {
    performClickAction(clickBinding_[cbp], e);
  } else if (mouseGrabber()) {
    if (mouseGrabberIsAManipulatedFrame_) {
      for (QMap<MouseBindingPrivate, MouseActionPrivate>::ConstIterator
               it = mouseBinding_.begin(),
               end = mouseBinding_.end();
           it != end; ++it)
        if ((it.value().handler == qglviewer::FRAME) && (it.key().button == e->button())) {
          qglviewer::ManipulatedFrame *mf =
              dynamic_cast<qglviewer::ManipulatedFrame *>(mouseGrabber());
          if (mouseGrabberIsAManipulatedCameraFrame_) {
            mf->qglviewer::ManipulatedFrame::startAction(it.value().action,
                                              it.value().withConstraint);
            mf->qglviewer::ManipulatedFrame::mousePressEvent(e, camera());
          } else {
            mf->startAction(it.value().action, it.value().withConstraint);
            mf->mousePressEvent(e, camera());
          }
          break;
        }
    } else
      mouseGrabber()->mousePressEvent(e, camera());
    update();
  } else {
    //#CONNECTION# wheelEvent has the same structure
    const MouseBindingPrivate mbp(e->modifiers(), e->button(),
                                  currentlyPressedKey_);

    if (mouseBinding_.contains(mbp)) {
      MouseActionPrivate map = mouseBinding_[mbp];
      switch (map.handler) {
      case qglviewer::CAMERA:
        camera()->frame()->startAction(map.action, map.withConstraint);
        camera()->frame()->mousePressEvent(e, camera());
        break;
      case qglviewer::FRAME:
        if (manipulatedFrame()) {
          if (manipulatedFrameIsACamera_) {
            manipulatedFrame()->ManipulatedFrame::startAction(
                map.action, map.withConstraint);
            manipulatedFrame()->ManipulatedFrame::mousePressEvent(e, camera());
          } else {
            manipulatedFrame()->startAction(map.action, map.withConstraint);
            manipulatedFrame()->mousePressEvent(e, camera());
          }
        }
        break;
      }
      if (map.action == qglviewer::SCREEN_ROTATE)
        // Display visual hint line
        update();
    } else
      e->ignore();
  }
}

/*! Overloading of the \c QWidget method.

Mouse move event is sent to the mouseGrabber() (if any) or to the camera() or
the manipulatedFrame(), depending on mouse bindings (see setMouseBinding()).

If you want to define your own mouse behavior, do something like this:
\code
CGAL_INLINE_FUNCTION
void Viewer::mousePressEvent(QMouseEvent* e)
{

if ((e->button() == myButton) && (e->modifiers() == myModifiers))
  myMouseBehavior = true;
else
  CGAL::QGLViewer::mousePressEvent(e);
}

CGAL_INLINE_FUNCTION
void Viewer::mouseMoveEvent(QMouseEvent *e)
{
if (myMouseBehavior)
  // Use e->position().x() and e->position().y() as you want...
else
  CGAL::QGLViewer::mouseMoveEvent(e);
}

CGAL_INLINE_FUNCTION
void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
if (myMouseBehavior)
  myMouseBehavior = false;
else
  CGAL::QGLViewer::mouseReleaseEvent(e);
}
\endcode */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::mouseMoveEvent(QMouseEvent *e) {
  if (mouseGrabber()) {
    mouseGrabber()->checkIfGrabsMouse(e->position().x(), e->position().y(), camera());
    if (mouseGrabber()->grabsMouse())
      if (mouseGrabberIsAManipulatedCameraFrame_)
        (dynamic_cast<qglviewer::ManipulatedFrame *>(mouseGrabber()))
            ->qglviewer::ManipulatedFrame::mouseMoveEvent(e, camera());
      else
        mouseGrabber()->mouseMoveEvent(e, camera());
    else
      setMouseGrabber(nullptr);
    update();
  }

  if (!mouseGrabber()) {
    //#CONNECTION# mouseReleaseEvent has the same structure
    if (camera()->frame()->isManipulated()) {
      camera()->frame()->mouseMoveEvent(e, camera());
      // #CONNECTION# manipulatedCameraFrame::mouseMoveEvent specific if at the
      // beginning
      if (camera()->frame()->action_ == qglviewer::ZOOM_ON_REGION)
        update();
    } else // !
        if ((manipulatedFrame()) && (manipulatedFrame()->isManipulated()))
      if (manipulatedFrameIsACamera_)
        manipulatedFrame()->ManipulatedFrame::mouseMoveEvent(e, camera());
      else
        manipulatedFrame()->mouseMoveEvent(e, camera());
    else if (hasMouseTracking()) {
      Q_FOREACH (qglviewer::MouseGrabber *mg, qglviewer::MouseGrabber::MouseGrabberPool()) {
        mg->checkIfGrabsMouse(e->position().x(), e->position().y(), camera());
        if (mg->grabsMouse()) {
          setMouseGrabber(mg);
          // Check that MouseGrabber is not disabled
          if (mouseGrabber() == mg) {
            update();
            break;
          }
        }
      }
    }
  }
}

/*! Overloading of the \c QWidget method.

Calls the mouseGrabber(), camera() or manipulatedFrame \c mouseReleaseEvent
method.

See the mouseMoveEvent() documentation for an example of mouse behavior
customization. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::mouseReleaseEvent(QMouseEvent *e) {
  if (mouseGrabber()) {
    if (mouseGrabberIsAManipulatedCameraFrame_)
      (dynamic_cast<qglviewer::ManipulatedFrame *>(mouseGrabber()))
          ->qglviewer::ManipulatedFrame::mouseReleaseEvent(e, camera());
    else
      mouseGrabber()->mouseReleaseEvent(e, camera());
    mouseGrabber()->checkIfGrabsMouse(e->position().x(), e->position().y(), camera());
    if (!(mouseGrabber()->grabsMouse()))
      setMouseGrabber(nullptr);
    // update();
  } else
      //#CONNECTION# mouseMoveEvent has the same structure
      if (camera()->frame()->isManipulated()) {
    camera()->frame()->mouseReleaseEvent(e, camera());
  } else if ((manipulatedFrame()) && (manipulatedFrame()->isManipulated())) {
    if (manipulatedFrameIsACamera_)
      manipulatedFrame()->ManipulatedFrame::mouseReleaseEvent(e, camera());
    else
      manipulatedFrame()->mouseReleaseEvent(e, camera());
  } else
    e->ignore();

  // Not absolutely needed (see above commented code for the optimal version),
  // but may reveal useful for specific applications.
  update();
}

/*! Overloading of the \c QWidget method.

If defined, the wheel event is sent to the mouseGrabber(). It is otherwise sent
according to wheel bindings (see setWheelBinding()). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::wheelEvent(QWheelEvent *e) {
  makeCurrent();
  if (mouseGrabber()) {
    if (mouseGrabberIsAManipulatedFrame_) {
      for (QMap<WheelBindingPrivate, MouseActionPrivate>::ConstIterator
               it = wheelBinding_.begin(),
               end = wheelBinding_.end();
           it != end; ++it)
        if (it.value().handler == qglviewer::FRAME) {
          qglviewer::ManipulatedFrame *mf =
              dynamic_cast<qglviewer::ManipulatedFrame *>(mouseGrabber());
          if (mouseGrabberIsAManipulatedCameraFrame_) {
            mf->qglviewer::ManipulatedFrame::startAction(it.value().action,
                                              it.value().withConstraint);
            mf->qglviewer::ManipulatedFrame::wheelEvent(e, camera());
          } else {
            mf->startAction(it.value().action, it.value().withConstraint);
            mf->wheelEvent(e, camera());
          }
          break;
        }
    } else
      mouseGrabber()->wheelEvent(e, camera());
    update();
  } else {
    //#CONNECTION# mousePressEvent has the same structure
    WheelBindingPrivate wbp(e->modifiers(), currentlyPressedKey_);

    if (wheelBinding_.contains(wbp)) {
      MouseActionPrivate map = wheelBinding_[wbp];
      switch (map.handler) {
      case qglviewer::CAMERA:
        if(currentlyPressedKey_ == ::Qt::Key_Z && _first_tick)
        {
          _first_tick = false;
          makeCurrent();
          //orient camera to the cursor.
          bool found = false;
          qglviewer::Vec point;
          point = camera()->pointUnderPixel(mapFromGlobal(QCursor::pos()), found);
          if(found)
            camera()->lookAt(point);
        }
        camera()->frame()->startAction(map.action, map.withConstraint);
        camera()->frame()->wheelEvent(e, camera());
        break;
      case qglviewer::FRAME:
        if (manipulatedFrame()) {
          if (manipulatedFrameIsACamera_) {
            manipulatedFrame()->ManipulatedFrame::startAction(
                map.action, map.withConstraint);
            manipulatedFrame()->ManipulatedFrame::wheelEvent(e, camera());
          } else {
            manipulatedFrame()->startAction(map.action, map.withConstraint);
            manipulatedFrame()->wheelEvent(e, camera());
          }
        }
        break;
      }
    } else
      e->ignore();
  }
}

/*! Overloading of the \c QWidget method.

The behavior of the mouse double click depends on the mouse binding. See
setMouseBinding() and the <a href="../mouse.html">mouse page</a>. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::mouseDoubleClickEvent(QMouseEvent *e) {
  //#CONNECTION# mousePressEvent has the same structure
  ClickBindingPrivate cbp(e->modifiers(), e->button(), true,
                          (::Qt::MouseButtons)(e->buttons() & ~(e->button())),
                          currentlyPressedKey_);
  if (clickBinding_.contains(cbp))
    performClickAction(clickBinding_[cbp], e);
  else if (mouseGrabber())
    mouseGrabber()->mouseDoubleClickEvent(e, camera());
  else
    e->ignore();
}


/*! Sets the isFullScreen() state.

If the CGAL::QGLViewer is embedded in an other QWidget (see
CGAL_INLINE_FUNCTION
QWidget::topLevelWidget()), this widget is displayed in full screen instead. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setFullScreen(bool fullScreen) {
  if (fullScreen_ == fullScreen)
    return;

  fullScreen_ = fullScreen;

  QWidget *tlw = topLevelWidget();

  if (isFullScreen()) {
    prevPos_ = topLevelWidget()->pos();
    tlw->showFullScreen();
    tlw->move(0, 0);
  } else {
    tlw->showNormal();
    tlw->move(prevPos_);
  }
}

/*! Directly defines the mouseGrabber().

You should not call this method directly as it bypasses the
CGAL::qglviewer::MouseGrabber::checkIfGrabsMouse() test performed by mouseMoveEvent().

If the MouseGrabber is disabled (see mouseGrabberIsEnabled()), this method
silently does nothing. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseGrabber(qglviewer::MouseGrabber *mouseGrabber) {
  if (!mouseGrabberIsEnabled(mouseGrabber))
    return;

  mouseGrabber_ = mouseGrabber;

  mouseGrabberIsAManipulatedFrame_ =
      (dynamic_cast<qglviewer::ManipulatedFrame *>(mouseGrabber) != nullptr);
  mouseGrabberIsAManipulatedCameraFrame_ =
      ((dynamic_cast<qglviewer::ManipulatedCameraFrame *>(mouseGrabber) != nullptr) &&
       (mouseGrabber != camera()->frame()));
  Q_EMIT mouseGrabberChanged(mouseGrabber);
}

/*! Sets the mouseGrabberIsEnabled() state. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseGrabberIsEnabled(
    const CGAL::qglviewer::MouseGrabber *const mouseGrabber, bool enabled) {
  if (enabled)
    disabledMouseGrabbers_.remove(reinterpret_cast<size_t>(mouseGrabber));
  else
    disabledMouseGrabbers_[reinterpret_cast<size_t>(mouseGrabber)];
}

CGAL_INLINE_FUNCTION
QString CGAL::QGLViewer::mouseActionString(qglviewer::MouseAction ma) {
  switch (ma) {
  case CGAL::qglviewer::NO_MOUSE_ACTION:
    return QString();
  case CGAL::qglviewer::ROTATE:
    return CGAL::QGLViewer::tr("Rotates", "ROTATE mouse action");
  case CGAL::qglviewer::ZOOM:
    return CGAL::QGLViewer::tr("Zooms", "ZOOM mouse action");
  case CGAL::qglviewer::TRANSLATE:
    return CGAL::QGLViewer::tr("Translates", "TRANSLATE mouse action");
  case CGAL::qglviewer::MOVE_FORWARD:
    return CGAL::QGLViewer::tr("Moves forward", "MOVE_FORWARD mouse action");
  case CGAL::qglviewer::LOOK_AROUND:
    return CGAL::QGLViewer::tr("Looks around", "LOOK_AROUND mouse action");
  case CGAL::qglviewer::MOVE_BACKWARD:
    return CGAL::QGLViewer::tr("Moves backward", "MOVE_BACKWARD mouse action");
  case CGAL::qglviewer::SCREEN_ROTATE:
    return CGAL::QGLViewer::tr("Rotates in screen plane",
                         "SCREEN_ROTATE mouse action");
  case CGAL::qglviewer::ROLL:
    return CGAL::QGLViewer::tr("Rolls", "ROLL mouse action");
  case CGAL::qglviewer::DRIVE:
    return CGAL::QGLViewer::tr("Drives", "DRIVE mouse action");
  case CGAL::qglviewer::SCREEN_TRANSLATE:
    return CGAL::QGLViewer::tr("Horizontally/Vertically translates",
                         "SCREEN_TRANSLATE mouse action");
  case CGAL::qglviewer::ZOOM_ON_REGION:
    return CGAL::QGLViewer::tr("Zooms on region for", "ZOOM_ON_REGION mouse action");
  case CGAL::qglviewer::ZOOM_FOV:
    return CGAL::QGLViewer::tr("Changes the FOV to emulate an optical zoom for ", "ZOOM_FOV mouse action");
  }
  return QString();
}

CGAL_INLINE_FUNCTION
QString CGAL::QGLViewer::clickActionString(CGAL::qglviewer::ClickAction ca) {
  switch (ca) {
  case CGAL::qglviewer::NO_CLICK_ACTION:
    return QString();
  case CGAL::qglviewer::ZOOM_ON_PIXEL:
    return CGAL::QGLViewer::tr("Zooms on pixel", "ZOOM_ON_PIXEL click action");
  case CGAL::qglviewer::ZOOM_TO_FIT:
    return CGAL::QGLViewer::tr("Zooms to fit scene", "ZOOM_TO_FIT click action");
  case CGAL::qglviewer::SELECT:
    return CGAL::QGLViewer::tr("Selects", "SELECT click action");
  case CGAL::qglviewer::RAP_FROM_PIXEL:
    return CGAL::QGLViewer::tr("Sets pivot point", "RAP_FROM_PIXEL click action");
  case CGAL::qglviewer::RAP_IS_CENTER:
    return CGAL::QGLViewer::tr("Resets pivot point", "RAP_IS_CENTER click action");
  case CGAL::qglviewer::CENTER_FRAME:
    return CGAL::QGLViewer::tr("Centers manipulated frame",
                         "CENTER_FRAME click action");
  case CGAL::qglviewer::CENTER_SCENE:
    return CGAL::QGLViewer::tr("Centers scene", "CENTER_SCENE click action");
  case CGAL::qglviewer::SHOW_ENTIRE_SCENE:
    return CGAL::QGLViewer::tr("Shows entire scene",
                         "SHOW_ENTIRE_SCENE click action");
  case CGAL::qglviewer::ALIGN_FRAME:
    return CGAL::QGLViewer::tr("Aligns manipulated frame",
                         "ALIGN_FRAME click action");
  case CGAL::qglviewer::ALIGN_CAMERA:
    return CGAL::QGLViewer::tr("Aligns camera", "ALIGN_CAMERA click action");
  }
  return QString();
}

static QString keyString(QKeyCombination key) {
  return QKeySequence(key).toString(QKeySequence::NativeText);
}

CGAL_INLINE_FUNCTION
QString CGAL::QGLViewer::formatClickActionPrivate(ClickBindingPrivate cbp) {
  bool buttonsBefore = cbp.buttonsBefore != ::Qt::NoButton;
  QString keyModifierString = keyString(QKeyCombination(cbp.modifiers, cbp.key));
  if (!keyModifierString.isEmpty()) {
#ifdef Q_OS_MAC
    // modifiers never has a '+' sign. Add one space to clearly separate
    // modifiers (and possible key) from button
    keyModifierString += " ";
#else
    // modifiers might be of the form : 'S' or 'Ctrl+S' or 'Ctrl+'. For
    // consistency, add an other '+' if needed, no spaces
    if (!keyModifierString.endsWith('+'))
      keyModifierString += "+";
#endif
  }

  return tr("%1%2%3%4%5%6", "Modifier / button or wheel / double click / with "
                            "/ button / pressed")
      .arg(keyModifierString)
      .arg(mouseButtonsString(cbp.button) +
           (cbp.button == ::Qt::NoButton ? tr("Wheel", "Mouse wheel") : ""))
      .arg(cbp.doubleClick ? tr(" double click", "Suffix after mouse button")
                           : "")
      .arg(buttonsBefore ? tr(" with ", "As in : Left button with Ctrl pressed")
                         : "")
      .arg(buttonsBefore ? mouseButtonsString(cbp.buttonsBefore) : "")
      .arg(buttonsBefore
               ? tr(" pressed", "As in : Left button with Ctrl pressed")
               : "");
}

CGAL_INLINE_FUNCTION
bool CGAL::QGLViewer::isValidShortcutKey(int key) {
  return (key >= ::Qt::Key_Any && key < ::Qt::Key_Escape) ||
         (key >= ::Qt::Key_F1 && key <= ::Qt::Key_F35);
}


/*! Defines a custom mouse binding description, displayed in the help() window's
 Mouse tab.

 Same as calling setMouseBindingDescription(::Qt::Key, ::Qt::KeyboardModifiers,
 ::Qt::MouseButton, QString, bool, ::Qt::MouseButtons), with a key value of
 ::Qt::Key(0) (i.e. binding description when no regular key needs to be pressed).
 */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseBindingDescription(::Qt::KeyboardModifiers modifiers,
                                           ::Qt::MouseButton button,
                                           QString description,
                                           bool doubleClick,
                                           ::Qt::MouseButtons buttonsBefore) {
  setMouseBindingDescription(::Qt::Key(0), modifiers, button, description,
                             doubleClick, buttonsBefore);
}

/*! Defines a custom mouse binding description, displayed in the help() window's
Mouse tab.

\p modifiers is a combination of ::Qt::KeyboardModifiers (\c ::Qt::ControlModifier,
\c ::Qt::AltModifier, \c ::Qt::ShiftModifier, \c ::Qt::MetaModifier). Possibly
combined using the \c "|" operator.

\p button is one of the ::Qt::MouseButtons (\c ::Qt::LeftButton, \c ::Qt::MiddleButton,
\c ::Qt::RightButton...).

\p doubleClick indicates whether or not the user has to double click this button
to perform the described action. \p buttonsBefore lists the buttons that need to
be pressed before the double click.

Set an empty \p description to \e remove a mouse binding description.

\code
// The R key combined with the Left mouse button rotates the camera in the
screen plane. setMouseBindingDescription(::Qt::Key_R, ::Qt::NoModifier,
::Qt::LeftButton, "Rotates camera in screen plane");

// A left button double click toggles full screen
setMouseBindingDescription(::Qt::NoModifier, ::Qt::LeftButton, "Toggles full screen
mode", true);

// Removes the description of Ctrl+Right button
setMouseBindingDescription(::Qt::ControlModifier, ::Qt::RightButton, "");
\endcode

Overload mouseMoveEvent() and friends to implement your custom mouse behavior
(see the mouseMoveEvent() documentation for an example). See the <a
href="../examples/keyboardAndMouse.html">keyboardAndMouse example</a> for an
illustration.

Use setMouseBinding() and setWheelBinding() to change the standard mouse action
bindings. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseBindingDescription(
    ::Qt::Key key, ::Qt::KeyboardModifiers modifiers, ::Qt::MouseButton button,
    QString description, bool doubleClick, ::Qt::MouseButtons buttonsBefore) {
  ClickBindingPrivate cbp(modifiers, button, doubleClick, buttonsBefore, key);

  if (description.isEmpty())
    mouseDescription_.remove(cbp);
  else
    mouseDescription_[cbp] = description;
}

static QString tableLine(const QString &left, const QString &right) {
  static bool even = false;
  const QString tdtd("</b></td><td>");
  const QString tdtr("</td></tr>\n");

  QString res("<tr bgcolor=\"");

  if (even)
    res += "#eeeeff\">";
  else
    res += "#ffffff\">";
  res += "<td><b>" + left + tdtd + right + tdtr;
  even = !even;

  return res;
}

/*! Returns a QString that describes the application mouse bindings, displayed
in the help() window \c Mouse tab.

Result is a table that describes custom application mouse binding descriptions
defined using setMouseBindingDescription() as well as standard mouse bindings
(defined using setMouseBinding() and setWheelBinding()). See the <a
href="../mouse.html">mouse page</a> for details on mouse bindings.

See also helpString() and keyboardString(). */
CGAL_INLINE_FUNCTION
QString CGAL::QGLViewer::mouseString() const {
  QString text(
      "<center><table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">\n");
  const QString trtd("<tr><td>");
  const QString tdtr("</td></tr>\n");
  const QString tdtd("</td><td>");

  text += QString("<tr bgcolor=\"#aaaacc\"><th align=\"center\">%1</th><th "
                  "align=\"center\">%2</th></tr>\n")
              .arg(tr("Button(s)",
                      "Buttons column header in help window mouse tab"))
              .arg(tr("Description",
                      "Description column header in help window mouse tab"));

  QMap<ClickBindingPrivate, QString> mouseBinding;

  // User-defined mouse bindings come first.
  for (QMap<ClickBindingPrivate, QString>::ConstIterator
           itm = mouseDescription_.begin(),
           endm = mouseDescription_.end();
       itm != endm; ++itm)
    mouseBinding[itm.key()] = itm.value();

  for (QMap<ClickBindingPrivate, QString>::ConstIterator
           it = mouseBinding.begin(),
           end = mouseBinding.end();
       it != end; ++it) {
    // Should not be needed (see setMouseBindingDescription())
    if (it.value().isNull())
      continue;

    text += tableLine(formatClickActionPrivate(it.key()), it.value());
  }

  // Optional separator line
  if (!mouseBinding.isEmpty()) {
    mouseBinding.clear();
    text += QString("<tr bgcolor=\"#aaaacc\"><td colspan=2>%1</td></tr>\n")
                .arg(tr("Standard mouse bindings", "In help window mouse tab"));
  }

  // Then concatenates the descriptions of wheelBinding_, mouseBinding_ and
  // clickBinding_. The order is significant and corresponds to the priorities
  // set in mousePressEvent() (reverse priority order, last one overwrites
  // previous) #CONNECTION# mousePressEvent() order
  for (QMap<MouseBindingPrivate, MouseActionPrivate>::ConstIterator
           itmb = mouseBinding_.begin(),
           endmb = mouseBinding_.end();
       itmb != endmb; ++itmb) {
    ClickBindingPrivate cbp(itmb.key().modifiers, itmb.key().button, false,
                            ::Qt::NoButton, itmb.key().key);

    QString text = mouseActionString(itmb.value().action);

    if (!text.isNull()) {
      switch (itmb.value().handler) {
      case qglviewer::CAMERA:
        text += " " + tr("camera", "Suffix after action");
        break;
      case qglviewer::FRAME:
        text += " " + tr("manipulated frame", "Suffix after action");
        break;
      }
      if (!(itmb.value().withConstraint))
        text += "*";
    }
    mouseBinding[cbp] = text;
  }

  for (QMap<WheelBindingPrivate, MouseActionPrivate>::ConstIterator
           itw = wheelBinding_.begin(),
           endw = wheelBinding_.end();
       itw != endw; ++itw) {
    ClickBindingPrivate cbp(itw.key().modifiers, ::Qt::NoButton, false,
                            ::Qt::NoButton, itw.key().key);

    QString text = mouseActionString(itw.value().action);

    if (!text.isNull()) {
      switch (itw.value().handler) {
      case qglviewer::CAMERA:
        text += " " + tr("camera", "Suffix after action");
        break;
      case qglviewer::FRAME:
        text += " " + tr("manipulated frame", "Suffix after action");
        break;
      }
      if (!(itw.value().withConstraint))
        text += "*";
    }

    mouseBinding[cbp] = text;
  }

  for (QMap<ClickBindingPrivate, qglviewer::ClickAction>::ConstIterator
           itcb = clickBinding_.begin(),
           endcb = clickBinding_.end();
       itcb != endcb; ++itcb)
    mouseBinding[itcb.key()] = clickActionString(itcb.value());

  for (QMap<ClickBindingPrivate, QString>::ConstIterator
           it2 = mouseBinding.begin(),
           end2 = mouseBinding.end();
       it2 != end2; ++it2) {
    if (it2.value().isNull())
      continue;

    text += tableLine(formatClickActionPrivate(it2.key()), it2.value());
  }

  text += "</table></center>";

  return text;
}

/*! Defines a custom keyboard shortcut description, that will be displayed in
the help() window \c Keyboard tab.

The \p key definition is given as an \c QKeyCombination using Qt enumerated values. Set an
empty \p description to remove a shortcut description: \code
setKeyDescription(::Qt::Key_W, "Toggles wireframe display");
setKeyDescription(::Qt::CTRL+::Qt::Key_L, "Loads a new scene");
// Removes a description
setKeyDescription(::Qt::CTRL+::Qt::Key_C, "");
\endcode

See the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse example</a>
for illustration and the <a href="../keyboard.html">keyboard page</a> for
details. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setKeyDescription(QKeyCombination key, QString description) {
  if (description.isEmpty())
    keyDescription_.remove(key);
  else
    keyDescription_[key] = description;
}

CGAL_INLINE_FUNCTION
QString CGAL::QGLViewer::cameraPathKeysString() const {
  if (pathIndex_.isEmpty())
    return QString();

  QVector< ::Qt::Key> keys;
  keys.reserve(pathIndex_.count());
  for (QMap< ::Qt::Key, unsigned int>::ConstIterator i = pathIndex_.begin(),
                                                  endi = pathIndex_.end();
       i != endi; ++i)
    keys.push_back(i.key());
  std::sort(keys.begin(), keys.end());

  QVector< ::Qt::Key>::const_iterator it = keys.begin(), end = keys.end();
  QString res = keyString(*it);

  const int maxDisplayedKeys = 6;
  int nbDisplayedKeys = 0;
  ::Qt::Key previousKey = (*it);
  int state = 0;
  ++it;
  while ((it != end) && (nbDisplayedKeys < maxDisplayedKeys - 1)) {
    switch (state) {
    case 0:
      if ((*it) == previousKey + 1)
        state++;
      else {
        res += ", " + keyString(*it);
        nbDisplayedKeys++;
      }
      break;
    case 1:
      if ((*it) == previousKey + 1)
        state++;
      else {
        res += ", " + keyString(previousKey);
        res += ", " + keyString(*it);
        nbDisplayedKeys += 2;
        state = 0;
      }
      break;
    default:
      if ((*it) != previousKey + 1) {
        res += ".." + keyString(previousKey);
        res += ", " + keyString(*it);
        nbDisplayedKeys += 2;
        state = 0;
      }
      break;
    }
    previousKey = *it;
    ++it;
  }

  if (state == 1)
    res += ", " + keyString(previousKey);
  if (state == 2)
    res += ".." + keyString(previousKey);
  if (it != end)
    res += "...";

  return res;
}

/*! Returns a QString that describes the application keyboard shortcut bindings,
and that will be displayed in the help() window \c Keyboard tab.

Default value is a table that describes the custom shortcuts defined using
setKeyDescription() as well as the \e standard CGAL::QGLViewer::KeyboardAction
shortcuts (defined using setShortcut()). See the <a
href="../keyboard.html">keyboard page</a> for details on key customization.

See also helpString() and mouseString(). */
CGAL_INLINE_FUNCTION
QString CGAL::QGLViewer::keyboardString() const {
  QString text(
      "<center><table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">\n");
  text += QString("<tr bgcolor=\"#aaaacc\"><th align=\"center\">%1</th><th "
                  "align=\"center\">%2</th></tr>\n")
              .arg(CGAL::QGLViewer::tr("Key(s)",
                                 "Keys column header in help window mouse tab"))
              .arg(CGAL::QGLViewer::tr(
                  "Description",
                  "Description column header in help window mouse tab"));

  QHash<QKeyCombination, QString> keyDescription;

  // 1 - User defined key descriptions
  for (auto kd = keyDescription_.begin(), kdend = keyDescription_.end();
       kd != kdend; ++kd)
    keyDescription[kd.key()] = kd.value();

  // Add to text in sorted order
  for (auto kb = keyDescription.begin(), endb = keyDescription.end();
       kb != endb; ++kb)
    text += tableLine(keyString(kb.key()), kb.value());

  // 2 - Optional separator line
  if (!keyDescription.isEmpty()) {
    keyDescription.clear();
    text += QString("<tr bgcolor=\"#aaaacc\"><td colspan=2>%1</td></tr>\n")
                .arg(CGAL::QGLViewer::tr("Standard viewer keys",
                                   "In help window keys tab"));
  }

  // 3 - KeyboardAction bindings description
  for (auto it = keyboardBinding_.begin(), end = keyboardBinding_.end();
       it != end; ++it)
    if ((it.value() != QKeyCombination{}) &&
        ((!cameraIsInRotateMode()) ||
         ((it.key() != qglviewer::INCREASE_FLYSPEED) && (it.key() != qglviewer::DECREASE_FLYSPEED))))
      keyDescription[it.value()] = keyboardActionDescription_[it.key()];

  // Add to text in sorted order
  for (auto kb2 = keyDescription.begin(), endb2 = keyDescription.end();
       kb2 != endb2; ++kb2)
    text += tableLine(keyString(kb2.key()), kb2.value());

  // 4 - Camera paths keys description
  const QString cpks = cameraPathKeysString();
  if (!cpks.isNull()) {
    text += "<tr bgcolor=\"#ccccff\"><td colspan=2>\n";
    text += CGAL::QGLViewer::tr("Camera paths are controlled using the %1 keys "
                          "(noted <i>Fx</i> below):",
                          "Help window key tab camera keys")
                .arg(cpks) +
            "</td></tr>\n";
    text += tableLine(
        keyString(QKeyCombination(playPathKeyboardModifiers())) + "<i>" +
            CGAL::QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>",
        CGAL::QGLViewer::tr("Plays path (or resets saved position)"));
    text += tableLine(
        keyString(QKeyCombination(addKeyFrameKeyboardModifiers())) + "<i>" +
            CGAL::QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>",
        CGAL::QGLViewer::tr("Adds a key frame to path (or defines a position)"));
    text += tableLine(
        keyString(QKeyCombination(addKeyFrameKeyboardModifiers())) + "<i>" +
            CGAL::QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>+<i>" +
            CGAL::QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>",
        CGAL::QGLViewer::tr("Deletes path (or saved position)"));
  }
  text += "</table></center>";

  return text;
}

/*! Displays the help window "About" tab. See help() for details. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::aboutQGLViewer() {
  help();
  helpWidget()->setCurrentIndex(3);
}

/*! Opens a modal help window that includes four tabs, respectively filled with
helpString(), keyboardString(), mouseString() and about libCGAL::QGLViewer.

Rich html-like text can be used (see the QStyleSheet documentation). This method
is called when the user presses the CGAL::QGLViewer::HELP key (default is 'H').

You can use helpWidget() to access to the help widget (to add/remove tabs,
change layout...).

The helpRequired() signal is emitted. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::help() {
  Q_EMIT helpRequired();

  bool resize = false;
  int width = 600;
  int height = 400;

  static QString label[] = {tr("&Help", "Help window tab title"),
                            tr("&Keyboard", "Help window tab title"),
                            tr("&Mouse", "Help window tab title"),
                            tr("&About", "Help window about title")};

  if (!helpWidget()) {
    // Qt4 requires a nullptr parent...
    helpWidget_ = new QTabWidget(nullptr);
    helpWidget()->setWindowTitle(tr("Help", "Help window title"));

    resize = true;
    for (int i = 0; i < 4; ++i) {
      QTextEdit *tab = new QTextEdit(nullptr);
      tab->setReadOnly(true);

      helpWidget()->insertTab(i, tab, label[i]);
      if (i == 3) {
#include "resources/qglviewer-icon.xpm"
        QPixmap pixmap(qglviewer_icon);
        tab->document()->addResource(QTextDocument::ImageResource,
                                     QUrl("mydata://qglviewer-icon.xpm"),
                                     QVariant(pixmap));
      }
    }
  }

  for (int i = 0; i < 4; ++i) {
    QString text;
    switch (i) {
    case 0:
      text = helpString();
      break;
    case 1:
      text = keyboardString();
      break;
    case 2:
      text = mouseString();
      break;
    case 3:
      text = QString("<center><br><img src=\"mydata://qglviewer-icon.xpm\">") +
             tr("<h1>libQGLViewer</h1>"
                "<h3>Forked from version 2.7.0</h3><br>"
                "A versatile 3D viewer based on OpenGL and Qt<br>"
                "Copyright 2002-%2 Gilles Debunne<br>"
                "<code>%3</code>")
                 .arg("2014")
                 .arg("http://www.libqglviewer.com") +
             QString("</center>");
      break;
    default:
      break;
    }

    QTextEdit *textEdit = (QTextEdit *)(helpWidget()->widget(i));
    textEdit->setHtml(text);
    textEdit->setText(text);

    if (resize && (textEdit->height() > height))
      height = textEdit->height();
  }

  if (resize)
    helpWidget()->resize(width, height + 40); // 40 pixels is ~ tabs' height
  helpWidget()->show();
  helpWidget()->raise();
}

/*! Overloading of the \c QWidget method.

Default keyboard shortcuts are defined using setShortcut(). Overload this method
to implement a specific keyboard binding. Call the original method if you do not
catch the event to preserve the viewer default key bindings: \code void
CGAL_INLINE_FUNCTION
Viewer::keyPressEvent(QKeyEvent *e)
{
  // Defines the Alt+R shortcut.
  if ((e->key() == ::Qt::Key_R) && (e->modifiers() == ::Qt::AltModifier))
  {
  myResetFunction();
  update(); // Refresh display
  }
  else
  CGAL::QGLViewer::keyPressEvent(e);
}

// With Qt 2 or 3, you would retrieve modifiers keys using :
// const ::Qt::ButtonState modifiers = (::Qt::ButtonState)(e->state() &
::Qt::KeyButtonMask); \endcode When you define a new keyboard shortcut, use
setKeyDescription() to provide a short description which is displayed in the
help() window Keyboard tab. See the <a
href="../examples/keyboardAndMouse.html">keyboardAndMouse</a> example for an
illustration.

CGAL_INLINE_FUNCTION
See also QOpenGLWidget::keyReleaseEvent(). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::keyPressEvent(QKeyEvent *e) {
  if (e->key() == 0) {
    e->ignore();
    return;
  }

  const ::Qt::Key key = ::Qt::Key(e->key());

  if(key == ::Qt::Key_Z && ! e->isAutoRepeat())
  {
    _first_tick = true;
  }
  const ::Qt::KeyboardModifiers modifiers = e->modifiers();
  auto it = keyboardBinding_.begin(), end = keyboardBinding_.end();
  const QKeyCombination target{modifiers, key};
  while ((it != end) && (it.value() != target))
    ++it;

  if (it != end)
    handleKeyboardAction(it.key());
  else if (pathIndex_.contains(::Qt::Key(key))) {
    // Camera paths
    unsigned int index = pathIndex_[::Qt::Key(key)];

    // not safe, but try to double press on two viewers at the same time !
    static QElapsedTimer doublePress;

    if (modifiers == playPathKeyboardModifiers()) {
      qint64 elapsed = doublePress.restart();
      if ((elapsed < 250) && (index == previousPathId_))
        camera()->resetPath(index);
      else {
        // Stop previous interpolation before starting a new one.
        if (index != previousPathId_) {
          qglviewer::KeyFrameInterpolator *previous =
              camera()->keyFrameInterpolator(previousPathId_);
          if ((previous) && (previous->interpolationIsStarted()))
            previous->resetInterpolation();
        }
        camera()->playPath(index);
      }
      previousPathId_ = index;
    } else if (modifiers == addKeyFrameKeyboardModifiers()) {
      qint64 elapsed = doublePress.restart();
      if ((elapsed < 250) && (index == previousPathId_)) {
        if (camera()->keyFrameInterpolator(index)) {
          disconnect(camera()->keyFrameInterpolator(index),
                     SIGNAL(interpolated()), this, SLOT(update()));
          if (camera()->keyFrameInterpolator(index)->numberOfKeyFrames() > 1)
            displayMessage(
                tr("Path %1 deleted", "Feedback message").arg(index));
          else
            displayMessage(
                tr("Position %1 deleted", "Feedback message").arg(index));
          camera()->deletePath(index);
        }
      } else {
        bool nullBefore = (camera()->keyFrameInterpolator(index) == nullptr);
        camera()->addKeyFrameToPath(index);
        if (nullBefore)
          connect(camera()->keyFrameInterpolator(index), SIGNAL(interpolated()),
                  SLOT(update()));
        int nbKF = camera()->keyFrameInterpolator(index)->numberOfKeyFrames();
        if (nbKF > 1)
          displayMessage(tr("Path %1, position %2 added", "Feedback message")
                             .arg(index)
                             .arg(nbKF));
        else
          displayMessage(
              tr("Position %1 saved", "Feedback message").arg(index));
      }
      previousPathId_ = index;
    }
    update();
  } else {
    if (isValidShortcutKey(key))
      currentlyPressedKey_ = key;
    e->ignore();
  }
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::keyReleaseEvent(QKeyEvent *e) {
  if (isValidShortcutKey(e->key()))
    currentlyPressedKey_ = ::Qt::Key(0);
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::handleKeyboardAction(qglviewer::KeyboardAction id) {
  switch (id) {
  case qglviewer::DRAW_AXIS:
    toggleAxisIsDrawn();
    break;
  case qglviewer::DRAW_GRID:
    toggleGridIsDrawn();
    break;
  case qglviewer::DISPLAY_FPS:
    toggleFPSIsDisplayed();
    break;
  case qglviewer::ENABLE_TEXT:
    toggleTextIsEnabled();
    break;
  case qglviewer::EXIT_VIEWER:
    qApp->closeAllWindows();
    break;
  case qglviewer::FULL_SCREEN:
    toggleFullScreen();
    break;
  case qglviewer::ANIMATION:
    toggleAnimation();
    break;
  case qglviewer::HELP:
    help();
    break;
  case qglviewer::EDIT_CAMERA:
    toggleCameraIsEdited();
    break;
  case qglviewer::CAMERA_MODE:
    toggleCameraMode();
    displayMessage(cameraIsInRotateMode()
                       ? tr("Camera in observer mode", "Feedback message")
                       : tr("Camera in fly mode", "Feedback message"));
    break;

  case qglviewer::MOVE_CAMERA_LEFT:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        qglviewer::Vec(-10.0 * camera()->flySpeed(), 0.0, 0.0)));
    update();
    break;
  case qglviewer::MOVE_CAMERA_RIGHT:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        qglviewer::Vec(10.0 * camera()->flySpeed(), 0.0, 0.0)));
    update();
    break;
  case qglviewer::MOVE_CAMERA_UP:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        qglviewer::Vec(0.0, 10.0 * camera()->flySpeed(), 0.0)));
    update();
    break;
  case qglviewer::MOVE_CAMERA_DOWN:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        qglviewer::Vec(0.0, -10.0 * camera()->flySpeed(), 0.0)));
    update();
    break;

  case qglviewer::INCREASE_FLYSPEED:
    camera()->setFlySpeed(camera()->flySpeed() * 1.5);
    break;
  case qglviewer::DECREASE_FLYSPEED:
    camera()->setFlySpeed(camera()->flySpeed() / 1.5);
    break;
  }
}

/*! Callback method used when the widget size is modified.

If you overload this method, first call the inherited method. Also called when
the widget is created, before its first display. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::resizeGL(int width, int height) {
  QOpenGLWidget::resizeGL(width, height);
  glViewport(0, 0, GLint(width), GLint(height));
  camera()->setScreenWidthAndHeight(this->width(), this->height(), this->devicePixelRatio());
}

//////////////////////////////////////////////////////////////////////////
//              K e y b o a r d   s h o r t c u t s                     //
//////////////////////////////////////////////////////////////////////////

/*! Defines the shortcut() that triggers a given CGAL::QGLViewer::KeyboardAction.

Here are some examples:
\code
// Press 'Q' to exit application
setShortcut(EXIT_VIEWER, ::Qt::Key_Q);

// Alt+M toggles camera mode
setShortcut(CAMERA_MODE, ::Qt::ALT | ::Qt::Key_M);

// The DISPLAY_FPS action is disabled
setShortcut(DISPLAY_FPS, 0);
\endcode

Only one shortcut can be assigned to a given CGAL::QGLViewer::KeyboardAction (new
bindings replace previous ones). If several KeyboardAction are bound to the
same shortcut, only one of them is active. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setShortcut(qglviewer::KeyboardAction action, QKeyCombination key) {
  keyboardBinding_[action] = key;
}

/*! Returns the keyboard shortcut associated to a given
CGAL::QGLViewer::KeyboardAction.

Result is an \c unsigned \c int defined using Qt enumerated values, as in \c
::Qt::Key_Q or \c ::Qt::CTRL + ::Qt::Key_X. Use ::Qt::MODIFIER_MASK to separate the key
from the state keys. Returns \c 0 if the KeyboardAction is disabled (not
binded). Set using setShortcut().

If you want to define keyboard shortcuts for custom actions (say, open a scene
file), overload keyPressEvent() and then setKeyDescription().

These shortcuts and their descriptions are automatically included in the help()
window \c Keyboard tab.

See the <a href="../keyboard.html">keyboard page</a> for details and default
values and the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse</a>
example for a practical illustration. */
CGAL_INLINE_FUNCTION
QKeyCombination CGAL::QGLViewer::shortcut(qglviewer::KeyboardAction action) const {
  if (keyboardBinding_.contains(action))
    return keyboardBinding_[action];
  else
    return {};
}


///////     Key Frames associated keys       ///////

/*! Returns the keyboard key associated to camera Key Frame path \p index.

Default values are F1..F12 for indexes 1..12.

addKeyFrameKeyboardModifiers() (resp. playPathKeyboardModifiers()) define the
state key(s) that must be pressed with this key to add a KeyFrame to (resp. to
play) the associated Key Frame path. If you quickly press twice the pathKey(),
the path is reset (resp. deleted).

Use camera()->keyFrameInterpolator( \p index ) to retrieve the
KeyFrameInterpolator that defines the path.

If several keys are binded to a given \p index (see setPathKey()), one of them
is returned. Returns \c 0 if no key is associated with this index.

See also the <a href="../keyboard.html">keyboard page</a>. */
CGAL_INLINE_FUNCTION
::Qt::Key CGAL::QGLViewer::pathKey(unsigned int index) const {
  for (QMap< ::Qt::Key, unsigned int>::ConstIterator it = pathIndex_.begin(),
                                                  end = pathIndex_.end();
       it != end; ++it)
    if (it.value() == index)
      return it.key();
  return ::Qt::Key(0);
}

/*! Sets the pathKey() associated with the camera Key Frame path \p index.

Several keys can be binded to the same \p index. Use a negated \p key value to
delete the binding (the \p index value is then ignored): \code
// Press 'space' to play/pause/add/delete camera path of index 0.
setPathKey(::Qt::Key_Space, 0);

// Remove this binding
setPathKey(-::Qt::Key_Space);
\endcode */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setPathKey(int key, unsigned int index) {
  ::Qt::Key k = ::Qt::Key(abs(key));
  if (key < 0)
    pathIndex_.remove(k);
  else
    pathIndex_[k] = index;
}

/*! Sets the playPathKeyboardModifiers() value. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setPlayPathKeyboardModifiers(::Qt::KeyboardModifiers modifiers) {
  playPathKeyboardModifiers_ = modifiers;
}

/*! Sets the addKeyFrameKeyboardModifiers() value. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setAddKeyFrameKeyboardModifiers(
    ::Qt::KeyboardModifiers modifiers) {
  addKeyFrameKeyboardModifiers_ = modifiers;
}

/*! Returns the keyboard modifiers that must be pressed with a pathKey() to add
the current camera position to a KeyFrame path.

It can be \c ::Qt::NoModifier, \c ::Qt::ControlModifier, \c ::Qt::ShiftModifier, \c
::Qt::AltModifier, \c ::Qt::MetaModifier or a combination of these (using the
bitwise '|' operator).

Default value is ::Qt::AltModifier. Defined using
setAddKeyFrameKeyboardModifiers().

See also playPathKeyboardModifiers(). */
CGAL_INLINE_FUNCTION
::Qt::KeyboardModifiers CGAL::QGLViewer::addKeyFrameKeyboardModifiers() const {
  return addKeyFrameKeyboardModifiers_;
}

/*! Returns the keyboard modifiers that must be pressed with a pathKey() to play
a camera KeyFrame path.

It can be \c ::Qt::NoModifier, \c ::Qt::ControlModifier, \c ::Qt::ShiftModifier, \c
::Qt::AltModifier, \c ::Qt::MetaModifier or a combination of these (using the
bitwise '|' operator).

Default value is ::Qt::NoModifier. Defined using setPlayPathKeyboardModifiers().

See also addKeyFrameKeyboardModifiers(). */
CGAL_INLINE_FUNCTION
::Qt::KeyboardModifiers CGAL::QGLViewer::playPathKeyboardModifiers() const {
  return playPathKeyboardModifiers_;
}


////////////////////////////////////////////////////////////////////////////////
//              M o u s e   b e h a v i o r   s t a t e   k e y s             //
////////////////////////////////////////////////////////////////////////////////

/*! Defines a MouseAction binding.

  Same as calling setMouseBinding(::Qt::Key, ::Qt::KeyboardModifiers,
  ::Qt::MouseButton, MouseHandler, MouseAction, bool), with a key value of
  ::Qt::Key(0) (i.e. no regular extra key needs to be pressed to perform this
  action). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseBinding(::Qt::KeyboardModifiers modifiers,
                                ::Qt::MouseButton button, qglviewer::MouseHandler handler,
                                qglviewer::MouseAction action, bool withConstraint) {
  setMouseBinding(::Qt::Key(0), modifiers, button, handler, action,
                  withConstraint);
}

/*! Associates a MouseAction to any mouse \p button, while keyboard \p modifiers
and \p key are pressed. The receiver of the mouse events is a MouseHandler
(CGAL::QGLViewer::CAMERA or CGAL::QGLViewer::FRAME).

The parameters should read: when the mouse \p button is pressed, while the
keyboard \p modifiers and \p key are down, activate \p action on \p handler. Use
::Qt::NoModifier to indicate that no modifier key is needed, and a \p key value of
0 if no regular key has to be pressed (or simply use
setMouseBinding(::Qt::KeyboardModifiers, ::Qt::MouseButton, MouseHandler,
MouseAction, bool)).

Use the '|' operator to combine modifiers:
\code
// The R key combined with the Left mouse button rotates the camera in the
screen plane. setMouseBinding(::Qt::Key_R, ::Qt::NoModifier, ::Qt::LeftButton, CAMERA,
SCREEN_ROTATE);

// Alt + Shift and Left button rotates the manipulatedFrame().
setMouseBinding(::Qt::AltModifier | ::Qt::ShiftModifier, ::Qt::LeftButton, FRAME,
ROTATE); \endcode

If \p withConstraint is \c true (default), the possible
CGAL::qglviewer::Frame::constraint() of the associated Frame will be enforced during
motion.

The list of all possible MouseAction, some binding examples and default bindings
are provided in the <a href="../mouse.html">mouse page</a>.

See the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse</a> example
for an illustration.

If no mouse button is specified, the binding is ignored. If an action was
previously associated with this keyboard and button combination, it is silently
overwritten (call mouseAction() before to check).

To remove a specific mouse binding, use \p NO_MOUSE_ACTION as the \p action.

See also setMouseBinding(::Qt::KeyboardModifiers, ::Qt::MouseButtons, ClickAction,
bool, int), setWheelBinding() and clearMouseBindings(). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseBinding(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                                ::Qt::MouseButton button, qglviewer::MouseHandler handler,
                                qglviewer::MouseAction action, bool withConstraint) {
  if ((handler == qglviewer::FRAME) &&
      ((action == qglviewer::MOVE_FORWARD) || (action == qglviewer::MOVE_BACKWARD) ||
       (action == qglviewer::ROLL) || (action == qglviewer::LOOK_AROUND) ||
       (action == qglviewer::ZOOM_ON_REGION))) {
    qWarning("Cannot bind %s to FRAME",
             mouseActionString(action).toLatin1().constData());
    return;
  }

  if (button == ::Qt::NoButton) {
    qWarning("No mouse button specified in setMouseBinding");
    return;
  }

  MouseActionPrivate map;
  map.handler = handler;
  map.action = action;
  map.withConstraint = withConstraint;

  MouseBindingPrivate mbp(modifiers, button, key);
  if (action == qglviewer::NO_MOUSE_ACTION)
    mouseBinding_.remove(mbp);
  else
    mouseBinding_.insert(mbp, map);

  ClickBindingPrivate cbp(modifiers, button, false, ::Qt::NoButton, key);
  clickBinding_.remove(cbp);
}


/*! Defines a ClickAction binding.

 Same as calling setMouseBinding(::Qt::Key, ::Qt::KeyboardModifiers,
 ::Qt::MouseButton, ClickAction, bool, ::Qt::MouseButtons), with a key value of
 ::Qt::Key(0) (i.e. no regular key needs to be pressed to activate this action).
 */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseBinding(::Qt::KeyboardModifiers modifiers,
                                ::Qt::MouseButton button, qglviewer::ClickAction action,
                                bool doubleClick,
                                ::Qt::MouseButtons buttonsBefore) {
  setMouseBinding(::Qt::Key(0), modifiers, button, action, doubleClick,
                  buttonsBefore);
}

/*! Associates a ClickAction to a button and keyboard key and modifier(s)
combination.

The parameters should read: when \p button is pressed, while the \p modifiers
and \p key keys are down, and possibly as a \p doubleClick, then perform \p
action. Use ::Qt::NoModifier to indicate that no modifier key is needed, and a \p
key value of 0 if no regular key has to be pressed (or simply use
setMouseBinding(::Qt::KeyboardModifiers, ::Qt::MouseButton, ClickAction, bool,
::Qt::MouseButtons)).

If \p buttonsBefore is specified (valid only when \p doubleClick is \c true),
then this (or these) other mouse button(s) has (have) to be pressed \e before
the double click occurs in order to execute \p action.

The list of all possible ClickAction, some binding examples and default bindings
are listed in the <a href="../mouse.html">mouse page</a>. See also the
setMouseBinding() documentation.

See the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse example</a>
for an illustration.

The binding is ignored if ::Qt::NoButton is specified as \p buttons.

See also setMouseBinding(::Qt::KeyboardModifiers, ::Qt::MouseButtons, MouseHandler,
MouseAction, bool), setWheelBinding() and clearMouseBindings().
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setMouseBinding(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                                ::Qt::MouseButton button, qglviewer::ClickAction action,
                                bool doubleClick,
                                ::Qt::MouseButtons buttonsBefore) {
  if ((buttonsBefore != ::Qt::NoButton) && !doubleClick) {
    qWarning("Buttons before is only meaningful when doubleClick is true in "
             "setMouseBinding().");
    return;
  }

  if (button == ::Qt::NoButton) {
    qWarning("No mouse button specified in setMouseBinding");
    return;
  }

  ClickBindingPrivate cbp(modifiers, button, doubleClick, buttonsBefore, key);

  // #CONNECTION performClickAction comment on NO_CLICK_ACTION
  if (action == qglviewer::NO_CLICK_ACTION)
    clickBinding_.remove(cbp);
  else
    clickBinding_.insert(cbp, action);

  if ((!doubleClick) && (buttonsBefore == ::Qt::NoButton)) {
    MouseBindingPrivate mbp(modifiers, button, key);
    mouseBinding_.remove(mbp);
  }
}

/*! Defines a mouse wheel binding.

 Same as calling setWheelBinding(::Qt::Key, ::Qt::KeyboardModifiers, MouseHandler,
 MouseAction, bool), with a key value of ::Qt::Key(0) (i.e. no regular key needs
 to be pressed to activate this action). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setWheelBinding(::Qt::KeyboardModifiers modifiers,
                                qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                                bool withConstraint) {
  setWheelBinding(::Qt::Key(0), modifiers, handler, action, withConstraint);
}

/*! Associates a MouseAction and a MouseHandler to a mouse wheel event.

This method is very similar to setMouseBinding(), but specific to the wheel.

In the current implementation only CGAL::QGLViewer::ZOOM can be associated with
CGAL::QGLViewer::FRAME, while CGAL::QGLViewer::CAMERA can receive CGAL::QGLViewer::ZOOM and
CGAL::QGLViewer::MOVE_FORWARD.

The difference between CGAL::QGLViewer::ZOOM and CGAL::QGLViewer::MOVE_FORWARD is that
CGAL::QGLViewer::ZOOM speed depends on the distance to the object, while
CGAL::QGLViewer::MOVE_FORWARD moves at a constant speed defined by
CGAL::qglviewer::Camera::flySpeed(). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setWheelBinding(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                                qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                                bool withConstraint) {
  //#CONNECTION# ManipulatedFrame::wheelEvent and
  // ManipulatedCameraFrame::wheelEvent switches
  if ((action != qglviewer::ZOOM) && (action != qglviewer::ZOOM_FOV) &&
      (action != qglviewer::MOVE_FORWARD) && (action != qglviewer::MOVE_BACKWARD)
      && (action != qglviewer::NO_MOUSE_ACTION)) {
    qWarning("Cannot bind %s to wheel",
             mouseActionString(action).toLatin1().constData());
    return;
  }

  if ((handler == qglviewer::FRAME) && (action != qglviewer::ZOOM)
      && (action != qglviewer::NO_MOUSE_ACTION)) {
    qWarning("Cannot bind %s to FRAME wheel",
             mouseActionString(action).toLatin1().constData());
    return;
  }

  MouseActionPrivate map;
  map.handler = handler;
  map.action = action;
  map.withConstraint = withConstraint;

  WheelBindingPrivate wbp(modifiers, key);
  if (action == qglviewer::NO_MOUSE_ACTION)
    wheelBinding_.remove(wbp);
  else
    wheelBinding_[wbp] = map;
}

/*! Clears all the default mouse bindings.

After this call, you will have to use setMouseBinding() and setWheelBinding() to
restore the mouse bindings you are interested in.
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::clearMouseBindings() {
  mouseBinding_.clear();
  clickBinding_.clear();
  wheelBinding_.clear();
}

/*! Clears all the default keyboard shortcuts.

After this call, you will have to use setShortcut() to define your own keyboard
shortcuts.
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::clearShortcuts() {
  keyboardBinding_.clear();
  pathIndex_.clear();
}

/*! Returns the MouseAction the will be triggered when the mouse \p button is
pressed, while the keyboard \p modifiers and \p key are pressed.

Returns CGAL::QGLViewer::NO_MOUSE_ACTION if no action is associated with this
combination. Use 0 for \p key to indicate that no regular key needs to be
pressed.

For instance, to know which motion corresponds to Alt+LeftButton, do:
\code
MouseAction ma = mouseAction(0, ::Qt::AltModifier, ::Qt::LeftButton);
if (ma != CGAL::QGLViewer::NO_MOUSE_ACTION) ...
\endcode

Use mouseHandler() to know which object (CGAL::QGLViewer::CAMERA or CGAL::QGLViewer::FRAME)
will execute this action. */
CGAL_INLINE_FUNCTION
qglviewer::MouseAction CGAL::QGLViewer::mouseAction(::Qt::Key key,
                                              ::Qt::KeyboardModifiers modifiers,
                                              ::Qt::MouseButton button) const {
  MouseBindingPrivate mbp(modifiers, button, key);
  if (mouseBinding_.contains(mbp))
    return mouseBinding_[mbp].action;
  else
    return qglviewer::NO_MOUSE_ACTION;
}


/*! Returns the MouseHandler which will be activated when the mouse \p button is
pressed, while the \p modifiers and \p key are pressed.

If no action is associated with this combination, returns \c -1. Use 0 for \p
key and ::Qt::NoModifier for \p modifiers to represent the lack of a key press.

For instance, to know which handler receives the Alt+LeftButton, do:
\code
int mh = mouseHandler(0, ::Qt::AltModifier, ::Qt::LeftButton);
if (mh == CGAL::QGLViewer::CAMERA) ...
\endcode

Use mouseAction() to know which action (see the MouseAction enum) will be
performed on this handler. */
CGAL_INLINE_FUNCTION
int CGAL::QGLViewer::mouseHandler(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                            ::Qt::MouseButton button) const {
  MouseBindingPrivate mbp(modifiers, button, key);
  if (mouseBinding_.contains(mbp))
    return mouseBinding_[mbp].handler;
  else
    return -1;
}



/*! Returns the keyboard state that triggers \p action on \p handler \p
withConstraint using the mouse wheel.

If such a binding exists, results are stored in the \p key and \p modifiers
parameters. If the MouseAction \p action is not bound, \p key is set to the
illegal -1 value. If several keyboard states trigger the MouseAction, one of
them is returned.

See also setMouseBinding(), getClickActionBinding() and getMouseActionBinding().
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::getWheelActionBinding(qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                                      bool withConstraint, ::Qt::Key &key,
                                      ::Qt::KeyboardModifiers &modifiers) const {
  for (QMap<WheelBindingPrivate, MouseActionPrivate>::ConstIterator
           it = wheelBinding_.begin(),
           end = wheelBinding_.end();
       it != end; ++it)
    if ((it.value().handler == handler) && (it.value().action == action) &&
        (it.value().withConstraint == withConstraint)) {
      key = it.key().key;
      modifiers = it.key().modifiers;
      return;
    }

  key = ::Qt::Key_unknown;
  modifiers = ::Qt::NoModifier;
}

/*! Returns the mouse and keyboard state that triggers \p action on \p handler
\p withConstraint.

If such a binding exists, results are stored in the \p key, \p modifiers and \p
button parameters. If the MouseAction \p action is not bound, \p button is set
to \c ::Qt::NoButton. If several mouse and keyboard states trigger the
MouseAction, one of them is returned.

See also setMouseBinding(), getClickActionBinding() and getWheelActionBinding().
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::getMouseActionBinding(qglviewer::MouseHandler handler, qglviewer::MouseAction action,
                                      bool withConstraint, ::Qt::Key &key,
                                      ::Qt::KeyboardModifiers &modifiers,
                                      ::Qt::MouseButton &button) const {
  for (QMap<MouseBindingPrivate, MouseActionPrivate>::ConstIterator
           it = mouseBinding_.begin(),
           end = mouseBinding_.end();
       it != end; ++it) {
    if ((it.value().handler == handler) && (it.value().action == action) &&
        (it.value().withConstraint == withConstraint)) {
      key = it.key().key;
      modifiers = it.key().modifiers;
      button = it.key().button;
      return;
    }
  }

  key = ::Qt::Key(0);
  modifiers = ::Qt::NoModifier;
  button = ::Qt::NoButton;
}

/*! Returns the MouseAction (if any) that is performed when using the wheel,
when the \p modifiers and \p key keyboard keys are pressed.

Returns NO_MOUSE_ACTION if no such binding has been defined using
setWheelBinding().

Same as mouseAction(), but for the wheel action. See also wheelHandler().
*/
qglviewer::MouseAction
CGAL_INLINE_FUNCTION
CGAL::QGLViewer::wheelAction(::Qt::Key key, ::Qt::KeyboardModifiers modifiers) const {
  WheelBindingPrivate wbp(modifiers, key);
  if (wheelBinding_.contains(wbp))
    return wheelBinding_[wbp].action;
  else
    return qglviewer::NO_MOUSE_ACTION;
}

/*! Returns the MouseHandler (if any) that receives wheel events when the \p
  modifiers and \p key keyboard keys are pressed.

  Returns -1 if no no such binding has been defined using setWheelBinding(). See
  also wheelAction().
*/
CGAL_INLINE_FUNCTION
int CGAL::QGLViewer::wheelHandler(::Qt::Key key,
                            ::Qt::KeyboardModifiers modifiers) const {
  WheelBindingPrivate wbp(modifiers, key);
  if (wheelBinding_.contains(wbp))
    return wheelBinding_[wbp].handler;
  else
    return -1;
}

/*! Same as mouseAction(), but for the ClickAction set using setMouseBinding().

Returns NO_CLICK_ACTION if no click action is associated with this keyboard and
mouse buttons combination. */
CGAL_INLINE_FUNCTION
CGAL::qglviewer::ClickAction
CGAL::QGLViewer::clickAction(::Qt::Key key, ::Qt::KeyboardModifiers modifiers,
                       ::Qt::MouseButton button, bool doubleClick,
                       ::Qt::MouseButtons buttonsBefore) const {
  ClickBindingPrivate cbp(modifiers, button, doubleClick, buttonsBefore, key);
  if (clickBinding_.contains(cbp))
    return clickBinding_[cbp];
  else
    return qglviewer::NO_CLICK_ACTION;
}

/*! Returns the mouse and keyboard state that triggers \p action.

If such a binding exists, results are stored in the \p key, \p modifiers, \p
button, \p doubleClick and \p buttonsBefore parameters. If the ClickAction \p
action is not bound, \p button is set to \c ::Qt::NoButton. If several mouse
buttons trigger in the ClickAction, one of them is returned.

See also setMouseBinding(), getMouseActionBinding() and getWheelActionBinding().
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::getClickActionBinding(qglviewer::ClickAction action, ::Qt::Key &key,
                                      ::Qt::KeyboardModifiers &modifiers,
                                      ::Qt::MouseButton &button,
                                      bool &doubleClick,
                                      ::Qt::MouseButtons &buttonsBefore) const {
  for (QMap<ClickBindingPrivate, qglviewer::ClickAction>::ConstIterator
           it = clickBinding_.begin(),
           end = clickBinding_.end();
       it != end; ++it)
    if (it.value() == action) {
      modifiers = it.key().modifiers;
      button = it.key().button;
      doubleClick = it.key().doubleClick;
      buttonsBefore = it.key().buttonsBefore;
      key = it.key().key;
      return;
    }

  modifiers = ::Qt::NoModifier;
  button = ::Qt::NoButton;
  doubleClick = false;
  buttonsBefore = ::Qt::NoButton;
  key = ::Qt::Key(0);
}

/*! This function should be used in conjunction with toggleCameraMode(). It
returns \c true when at least one mouse button is binded to the \c ROTATE
mouseAction. This is crude way of determining which "mode" the camera is in. */
CGAL_INLINE_FUNCTION
bool CGAL::QGLViewer::cameraIsInRotateMode() const {
  //#CONNECTION# used in toggleCameraMode() and keyboardString()
  ::Qt::Key key;
  ::Qt::KeyboardModifiers modifiers;
  ::Qt::MouseButton button;
  getMouseActionBinding(qglviewer::CAMERA, qglviewer::ROTATE, true /*constraint*/, key, modifiers,
                        button);
  return button != ::Qt::NoButton;
}

/*! Swaps between two predefined camera mouse bindings.

The first mode makes the camera observe the scene while revolving around the
CGAL::qglviewer::Camera::pivotPoint(). The second mode is designed for walkthrough
applications and simulates a flying camera.

Practically, the three mouse buttons are respectively binded to:
\arg In rotate mode: CGAL::QGLViewer::ROTATE, CGAL::QGLViewer::ZOOM, CGAL::QGLViewer::TRANSLATE.
\arg In fly mode: CGAL::QGLViewer::MOVE_FORWARD, CGAL::QGLViewer::LOOK_AROUND,
CGAL::QGLViewer::MOVE_BACKWARD.

The current mode is determined by checking if a mouse button is binded to
CGAL::QGLViewer::ROTATE for the CGAL::QGLViewer::CAMERA. The state key that was previously
used to move the camera is preserved. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::toggleCameraMode() {
  ::Qt::Key key;
  ::Qt::KeyboardModifiers modifiers;
  ::Qt::MouseButton button;
  getMouseActionBinding(qglviewer::CAMERA, qglviewer::ROTATE, true /*constraint*/, key, modifiers,
                        button);
  bool rotateMode = button != ::Qt::NoButton;

  if (!rotateMode) {
    getMouseActionBinding(qglviewer::CAMERA, qglviewer::MOVE_FORWARD, true /*constraint*/, key,
                          modifiers, button);
  }

  //#CONNECTION# setDefaultMouseBindings()
  if (rotateMode) {
    camera()->frame()->updateSceneUpVector();
    camera()->frame()->stopSpinning();

    setMouseBinding(modifiers, ::Qt::LeftButton, qglviewer::CAMERA, qglviewer::MOVE_FORWARD);
    setMouseBinding(modifiers, ::Qt::MiddleButton, qglviewer::CAMERA, qglviewer::LOOK_AROUND);
    setMouseBinding(modifiers, ::Qt::RightButton, qglviewer::CAMERA, qglviewer::MOVE_BACKWARD);

    setMouseBinding(::Qt::Key_R, modifiers, ::Qt::LeftButton, qglviewer::CAMERA, qglviewer::ROLL);

    setMouseBinding(::Qt::NoModifier, ::Qt::LeftButton, qglviewer::NO_CLICK_ACTION, true);
    setMouseBinding(::Qt::NoModifier, ::Qt::MiddleButton, qglviewer::NO_CLICK_ACTION, true);
    setMouseBinding(::Qt::NoModifier, ::Qt::RightButton, qglviewer::NO_CLICK_ACTION, true);

    setWheelBinding(modifiers, qglviewer::CAMERA, qglviewer::MOVE_FORWARD);
  } else {
    // Should stop flyTimer. But unlikely and not easy.
    setMouseBinding(modifiers, ::Qt::LeftButton, qglviewer::CAMERA, qglviewer::ROTATE);
    setMouseBinding(modifiers, ::Qt::MiddleButton, qglviewer::CAMERA, qglviewer::ZOOM);
    setMouseBinding(modifiers, ::Qt::RightButton, qglviewer::CAMERA, qglviewer::TRANSLATE);

    setMouseBinding(::Qt::Key_R, modifiers, ::Qt::LeftButton, qglviewer::CAMERA,
                    qglviewer::SCREEN_ROTATE);

    setMouseBinding(::Qt::NoModifier, ::Qt::LeftButton, qglviewer::ALIGN_CAMERA, true);
    setMouseBinding(::Qt::NoModifier, ::Qt::MiddleButton, qglviewer::SHOW_ENTIRE_SCENE, true);
    setMouseBinding(::Qt::NoModifier, ::Qt::RightButton, qglviewer::CENTER_SCENE, true);

    setWheelBinding(modifiers, qglviewer::CAMERA, qglviewer::ZOOM);
  }
}

////////////////////////////////////////////////////////////////////////////////
//              M a n i p u l a t e d   f r a m e s                           //
////////////////////////////////////////////////////////////////////////////////

/*! Sets the viewer's manipulatedFrame().

Several objects can be manipulated simultaneously, as is done the <a
href="../examples/multiSelect.html">multiSelect example</a>.

Defining the \e own viewer's camera()->frame() as the manipulatedFrame() is
possible and will result in a classical camera manipulation. See the <a
href="../examples/luxo.html">luxo example</a> for an illustration.

Note that a CGAL::qglviewer::ManipulatedCameraFrame can be set as the
manipulatedFrame(): it is possible to manipulate the camera of a first viewer in
a second viewer. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setManipulatedFrame(qglviewer::ManipulatedFrame *frame) {
  if (manipulatedFrame()) {
    manipulatedFrame()->stopSpinning();

    if (manipulatedFrame() != camera()->frame()) {
      disconnect(manipulatedFrame(), SIGNAL(manipulated()), this,
                 SLOT(update()));
      disconnect(manipulatedFrame(), SIGNAL(spun()), this, SLOT(update()));
    }
  }

  manipulatedFrame_ = frame;

  manipulatedFrameIsACamera_ =
      ((manipulatedFrame() != camera()->frame()) &&
       (dynamic_cast<qglviewer::ManipulatedCameraFrame *>(manipulatedFrame()) != nullptr));

  if (manipulatedFrame()) {
    // Prevent multiple connections, that would result in useless display
    // updates
    if (manipulatedFrame() != camera()->frame()) {
      connect(manipulatedFrame(), SIGNAL(manipulated()), SLOT(update()));
      connect(manipulatedFrame(), SIGNAL(spun()), SLOT(update()));
    }
  }
}

#ifndef DOXYGEN
////////////////////////////////////////////////////////////////////////////////
//                          V i s u a l   H i n t s                           //
////////////////////////////////////////////////////////////////////////////////
/*! Draws viewer related visual hints.

CGAL_INLINE_FUNCTION
Displays the new CGAL::qglviewer::Camera::pivotPoint() when it is changed. See the <a
href="../mouse.html">mouse page</a> for details. Also draws a line between
CGAL::qglviewer::Camera::pivotPoint() and mouse cursor when the camera is rotated
around the camera Z axis.

See also setVisualHintsMask() and resetVisualHints(). The hint color is
foregroundColor().

\note These methods may become more interesting one day. The current design is
too limited and should be improved when other visual hints must be drawn.

Limitation : One needs to have access to visualHint_ to overload this method.

Removed from the documentation for this reason. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::drawVisualHints() {

  // G r i d
  GLdouble mat[16];
  camera()->getModelViewProjectionMatrix(mat);

  QMatrix4x4 mvpMatrix;
  for(int i=0; i < 16; i++)
  {
    mvpMatrix.data()[i] = float(mat[i]);
  }

  QMatrix4x4 mvMatrix;
  for(int i=0; i < 16; i++)
  {
    mvMatrix.data()[i] = float(camera()->orientation().inverse().matrix()[i]);
  }

  rendering_program.bind();
  vaos[GRID].bind();
  rendering_program.setUniformValue("mvp_matrix", mvpMatrix);
  rendering_program.setUniformValue("color", QColor(::Qt::lightGray));
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(grid_size));
  vaos[GRID].release();

  vaos[GRID_AXIS_X].bind();
  rendering_program.setUniformValue("color", QColor(::Qt::red));
//  glEnable(GL_POLYGON_OFFSET_FILL);
//  glPolygonOffset(0.0f,-0.0f);
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(grid_axis_size/3));
  vaos[GRID_AXIS_X].release();

  vaos[GRID_AXIS_Y].bind();
  rendering_program.setUniformValue("color", QColor(::Qt::green));
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(grid_axis_size/3));
  vaos[GRID_AXIS_Y].release();

  //A x i s
  CGAL::qglviewer::Camera::Type camera_type = camera()->type();

  camera()->setType(CGAL::qglviewer::Camera::ORTHOGRAPHIC);
  for(int i=0; i < 16; i++)
  {
    mvMatrix.data()[i] = float(camera()->orientation().inverse().matrix()[i]);
  }
  mvpMatrix.setToIdentity();
  mvpMatrix.ortho(-1,1,-1,1,-1,1);
  mvpMatrix = mvpMatrix*mvMatrix;
  rendering_program.setUniformValue("mvp_matrix", mvpMatrix);

  camera()->setType(camera_type);

  // The viewport and the scissor are changed to fit the upper right corner.
  // Original values are saved.
  int viewport[4];
  int scissor[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetIntegerv(GL_SCISSOR_BOX, scissor);

  // Axis viewport size, in pixels
  int size = 100;
  glViewport(width()*devicePixelRatio()-size, height()*devicePixelRatio()-size, size, size);
  glScissor (width()*devicePixelRatio()-size, height()*devicePixelRatio()-size, size, size);

  vaos[AXIS_X].bind();
  rendering_program.setUniformValue("color", QColor(::Qt::red));
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(axis_size / 3));
  vaos[AXIS_X].release();

  vaos[AXIS_Y].bind();
  rendering_program.setUniformValue("color", QColor(::Qt::green));
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(axis_size / 3));
  vaos[AXIS_Y].release();

  vaos[AXIS_Z].bind();
  rendering_program.setUniformValue("color", QColor(::Qt::blue));
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(axis_size / 3));
  vaos[AXIS_Z].release();

  // The viewport and the scissor are restored.
  glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
  glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);

  //P i v o t - P o i n t
  if (visualHint_ & 1)
  {
    std::vector<float> vertices;
    for(int i=0; i< 4; ++i)
    {
      float x = float(std::pow(-1, i)*(1-i/2));
      float y = float(std::pow(-1, i)*(i/2));
      vertices.push_back(x);
      vertices.push_back(y);
      vertices.push_back(0);
    }

    rendering_program.bind();
    vaos[PIVOT_POINT].bind();
    vbos[Pivot_point].bind();

    vbos[Pivot_point].allocate(vertices.data(),static_cast<int>(vertices.size()*sizeof(float)));
    rendering_program.enableAttributeArray("vertex");
    rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
    vbos[Pivot_point].release();

    mvpMatrix.setToIdentity();
    mvpMatrix.ortho(-1,1,-1,1,-1,1);
    rendering_program.setUniformValue("mvp_matrix", mvpMatrix);

    const auto point_2d = camera()->projectedCoordinatesOf(camera()->pivotPoint());
    size = 30 * devicePixelRatio();

    glViewport(GLint(point_2d.x*devicePixelRatio()-size/2),
               GLint((height() - point_2d.y)*devicePixelRatio()-size/2), size, size);
    glScissor (GLint(point_2d.x*devicePixelRatio()-size/2),
               GLint((height() - point_2d.y)*devicePixelRatio()-size/2), size, size);

    rendering_program.setUniformValue("color", QColor(::Qt::black));
    glDisable(GL_DEPTH_TEST);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(4));
    rendering_program.setUniformValue("color", QColor(::Qt::white));
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(4));
    glEnable(GL_DEPTH_TEST);

    // The viewport and the scissor are restored.
    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);

    vaos[PIVOT_POINT].release();
    rendering_program.release();
  }

}

/*! Defines the mask that will be used to drawVisualHints(). The only available
mask is currently 1, corresponding to the display of the
CGAL::qglviewer::Camera::pivotPoint(). resetVisualHints() is automatically called
after \p delay milliseconds (default is 2 seconds). */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setVisualHintsMask(int mask, int delay) {
  visualHint_ = visualHint_ | mask;
  QTimer::singleShot(delay, this, SLOT(resetVisualHints()));
}

/*! Reset the mask used by drawVisualHints(). Called by setVisualHintsMask()
 * after 2 seconds to reset the display. */
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::resetVisualHints() { visualHint_ = 0; }
#endif

////////////////////////////////////////////////////////////////////////////////
//       A x i s   a n d   G r i d   d i s p l a y   l i s t s                //
////////////////////////////////////////////////////////////////////////////////

/*! Draws a 3D arrow between the 3D point \p from and the 3D point \p to.
\p data is filled with the three components of a point, which makes it filled like this:
[P1.x-P1.y-P1.z|P2.x-P2.y-P2.z|...]
*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::drawArrow(double r, double R, int prec,
                                const CGAL::qglviewer::Vec& from, const CGAL::qglviewer::Vec& to,
                                std::vector<float> &data) {
  using std::cos;
  using std::sin;
  using std::acos;

  CGAL::qglviewer::Vec temp = to-from;
  QVector3D dir = QVector3D(float(temp.x), float(temp.y), float(temp.z));

  QMatrix4x4 mat;
  mat.setToIdentity();
  mat.translate(float(from.x), float(from.y), float(from.z));
  mat.scale(dir.length());

  dir.normalize();
  float angle = 0.f;
  if(std::sqrt((dir.x()*dir.x()+dir.y()*dir.y())) > 1)
      angle = 90.f;
  else
      angle = float(acos(dir.y()/std::sqrt(dir.lengthSquared()))*180.0/CGAL_PI);

  QVector3D axis;
  axis = QVector3D(dir.z(), 0, -dir.x());
  mat.rotate(angle, axis);

  //Head
  const float Rf = static_cast<float>(R);
  for(int d = 0; d<360; d+= 360/prec)
  {
      float D = float(d * CGAL_PI / 180.);
//      float a = float(std::atan(Rf / 0.33));

      //point A1
      QVector4D p(0.f, 1.f, 0.f, 1.f);
      QVector4D pR = mat*p;
//      QVector4D n(Rf*sin(D), sin(a), Rf*cos(D), 1.f);
//      QVector4D nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point B1
      p = QVector4D(Rf*sin(D), 0.66f, Rf* cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), sin(a), cos(D), 1.);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point C1
      D = float((d+360/prec)*CGAL_PI/180.0);
      p = QVector4D(Rf* sin(D), 0.66f, Rf* cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), sin(a), cos(D), 1.f);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
  }

  //cylinder
  //body of the cylinder
  const float rf = static_cast<float>(r);
  for(int d = 0; d<360; d+= 360/prec)
  {
      //point A1
      float D = float(d*CGAL_PI/180.);
      QVector4D p(rf*sin(D), 0.66f, rf*cos(D), 1.f);
      QVector4D pR = mat*p;
//      QVector4D n(sin(D), 0.f, cos(D), 1.f);
//      QVector4D nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point B1
      p = QVector4D(rf * sin(D), 0.f, rf*cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), 0.f, cos(D), 1.f);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point C1
      D = float((d + 360./prec) * CGAL_PI/180.0);
      p = QVector4D(rf * sin(D), 0.f, rf*cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), 0.f, cos(D), 1.f);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point A2
      D = float((d + 360./prec) * CGAL_PI/180.);
      p = QVector4D(rf * sin(D), 0.f, rf*cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), 0.f, cos(D), 1.f);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point B2
      p = QVector4D(rf * sin(D), 0.66f, rf*cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), 0.f, cos(D), 1.f);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());

      //point C2
      D = float(d * CGAL_PI/180.);
      p = QVector4D(rf * sin(D), 0.66f, rf*cos(D), 1.f);
      pR = mat*p;
//      n = QVector4D(sin(D), 0.f, cos(D), 1.f);
//      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
  }
}

/*! Draws an XYZ axis, with a given size (default is 1.0).

The axis orientation matches the current modelView matrix state:
three arrows (red, green and blue) of length \p length are drawn along the
positive X, Y and Z directions in the top right corner of the screen.
X arrow is red, Y arrow is green and Z arrow is blue.*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::drawAxis(qreal length) {
  std::vector<float> axis_x_data, axis_y_data, axis_z_data;

  drawArrow(0.06,0.12,10, CGAL::qglviewer::Vec(0,0,0),CGAL::qglviewer::Vec(length,0,0), axis_x_data);
  drawArrow(0.06,0.12,10, CGAL::qglviewer::Vec(0,0,0),CGAL::qglviewer::Vec(0,length,0), axis_y_data);
  drawArrow(0.06,0.12,10, CGAL::qglviewer::Vec(0,0,0),CGAL::qglviewer::Vec(0,0,length), axis_z_data);

  rendering_program.bind();

  vaos[AXIS_X].bind();
  vbos[Axis_x].bind();
  vbos[Axis_x].allocate(axis_x_data.data(), static_cast<int>(axis_x_data.size() * sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Axis_x].release();
  vaos[AXIS_X].release();

  vaos[AXIS_Y].bind();
  vbos[Axis_y].bind();
  vbos[Axis_y].allocate(axis_y_data.data(), static_cast<int>(axis_y_data.size() * sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Axis_y].release();
  vaos[AXIS_Y].release();

  vaos[AXIS_Z].bind();
  vbos[Axis_z].bind();
  vbos[Axis_z].allocate(axis_z_data.data(), static_cast<int>(axis_z_data.size() * sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Axis_z].release();
  vaos[AXIS_Z].release();

  rendering_program.release();

  axis_size = axis_x_data.size(); // == axis_y_data.size() == axis_z_data.size()
}

/*! Draws a grid in the XY plane, centered on (0,0,0) (defined in the current
coordinate system).

\p size (OpenGL units) and \p nbSubdivisions define its geometry.*/
CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::drawGrid(qreal size, int nbSubdivisions) {

  //The Grid
  std::vector<float> grid_data;
  for (int i=0; i<=nbSubdivisions; ++i)
  {
    const float pos = float(size*(2.0*i/nbSubdivisions-1.0));
    grid_data.push_back(pos);
    grid_data.push_back(float(-size));
    grid_data.push_back(0.f);

    grid_data.push_back(pos);
    grid_data.push_back(float(+size));
    grid_data.push_back(0.f);

    grid_data.push_back(float(-size));
    grid_data.push_back(pos);
    grid_data.push_back(0.f);

    grid_data.push_back(float(size));
    grid_data.push_back(pos);
    grid_data.push_back(0.f);
  }

  rendering_program.bind();
  vaos[GRID].bind();
  vbos[Grid].bind();
  vbos[Grid].allocate(grid_data.data(),static_cast<int>(grid_data.size()*sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Grid].release();
  vaos[GRID].release();
  rendering_program.release();
  grid_size = grid_data.size();

  //The Axis
  std::vector<float> grid_axis_x_data, grid_axis_y_data;
  grid_axis_x_data.reserve(270); // default subdivision parameter yields this amount
  grid_axis_y_data.reserve(270);

  //axis_data is filled by drawArrow always this way : V.x V.y V.z
  drawArrow(0.005,0.02,10, CGAL::qglviewer::Vec(0,0,0),CGAL::qglviewer::Vec(size,0,0), grid_axis_x_data);
  drawArrow(0.005,0.02,10, CGAL::qglviewer::Vec(0,0,0),CGAL::qglviewer::Vec(0,size,0), grid_axis_y_data);

  rendering_program.bind();

  // X
  vaos[GRID_AXIS_X].bind();
  vbos[Grid_axis_x].bind();
  vbos[Grid_axis_x].allocate(grid_axis_x_data.data(), static_cast<int>(grid_axis_x_data.size() * sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Grid_axis_x].release();
  vaos[GRID_AXIS_X].release();

  // Y
  vaos[GRID_AXIS_Y].bind();
  vbos[Grid_axis_y].bind();
  vbos[Grid_axis_y].allocate(grid_axis_y_data.data(), static_cast<int>(grid_axis_y_data.size() * sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Grid_axis_y].release();
  vaos[GRID_AXIS_Y].release();

  rendering_program.release();

  grid_axis_size = grid_axis_x_data.size(); // == grid_axis_y_data.size()
}

////////////////////////////////////////////////////////////////////////////////
//       S t a t i c    m e t h o d s   :  Q G L V i e w e r   P o o l        //
////////////////////////////////////////////////////////////////////////////////


CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::copyBufferToTexture(GLint , GLenum ) {
}

/*! Returns the texture id of the texture created by copyBufferToTexture().

Use glBindTexture() to use this texture. Note that this is already done by
copyBufferToTexture().

Returns \c 0 is copyBufferToTexture() was never called or if the texture was
deleted using glDeleteTextures() since then. */
CGAL_INLINE_FUNCTION
GLuint CGAL::QGLViewer::bufferTextureId() const {
    return 0;
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setOffset(CGAL::qglviewer::Vec offset)
{
  this->_offset = offset;
}

CGAL_INLINE_FUNCTION
const CGAL::qglviewer::Vec& CGAL::QGLViewer::offset() const
{
  return _offset;
}

CGAL_INLINE_FUNCTION
qreal CGAL::QGLViewer::sceneRadius() const { return camera()->sceneRadius(); }

CGAL_INLINE_FUNCTION
CGAL::qglviewer::Vec CGAL::QGLViewer::sceneCenter() const { return camera()->sceneCenter(); }

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setSceneRadius(qreal radius) {
  camera()->setSceneRadius(radius);
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setSceneCenter(const CGAL::qglviewer::Vec &center) {
  camera()->setSceneCenter(center);
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setSceneBoundingBox(const CGAL::qglviewer::Vec &min,
                                          const CGAL::qglviewer::Vec &max) {
  camera()->setSceneBoundingBox(min + offset(), max + offset());
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::showEntireScene() {
  camera()->showEntireScene();
  update();
}


CGAL_INLINE_FUNCTION
QImage* CGAL::QGLViewer::takeSnapshot(CGAL::qglviewer::SnapShotBackground background_color,
                                      QSize finalSize,
                                      double oversampling,
                                      bool expand)
{
  makeCurrent();
  qreal aspectRatio = width() / static_cast<qreal>(height());
  GLfloat alpha = 1.0f;
  QColor previousBGColor = backgroundColor();
  switch(background_color)
  {
  case CGAL::qglviewer::CURRENT_BACKGROUND:
    break;
  case CGAL::qglviewer::TRANSPARENT_BACKGROUND:
    setBackgroundColor(QColor(::Qt::transparent));
    alpha = 0.0f;
    break;
  case CGAL::qglviewer::CHOOSE_BACKGROUND:
    QColor c =  QColorDialog::getColor();
    if(c.isValid()) {
      setBackgroundColor(c);
    }
    else
      return nullptr;
    break;
  }


  QSize subSize(int(width()/oversampling), int(height()/oversampling));
  QSize size=QSize(width(), height());

  qreal newAspectRatio = finalSize.width() / static_cast<qreal>(finalSize.height());

  qreal zNear = camera()->zNear();
  qreal zFar = camera()->zFar();

  qreal xMin, yMin;

  if(camera()->type()==CGAL::qglviewer::Camera::PERSPECTIVE)
  {
    if ((expand && (newAspectRatio>aspectRatio)) || (!expand && (newAspectRatio<aspectRatio)))
    {
      yMin = zNear * tan(camera()->fieldOfView() / 2.0);
      xMin = newAspectRatio * yMin;
    }
    else
    {
      xMin = zNear * tan(camera()->fieldOfView() / 2.0) * aspectRatio;
      yMin = xMin / newAspectRatio;
    }
  }
  else
  {
    double xy[6];
    camera()->getFrustum(xy);
    if ((expand && (newAspectRatio>aspectRatio)) || (!expand && (newAspectRatio<aspectRatio)))
    {
      yMin = -xy[3];
      xMin = newAspectRatio * yMin;
      ;
    }
    else
    {
      xMin = -xy[0] * aspectRatio;
      yMin = xMin / newAspectRatio;
    }
  }
  QImage *image = new QImage(finalSize.width(), finalSize.height(), QImage::Format_ARGB32);

  if (image->isNull())
  {
    QMessageBox::warning(this, "Image saving error",
                         "Unable to create resulting image",
                         QMessageBox::Ok, QMessageBox::NoButton);
    setBackgroundColor(previousBGColor);
    return nullptr;
  }

  qreal scaleX = subSize.width() / static_cast<qreal>(finalSize.width());
  qreal scaleY = subSize.height() / static_cast<qreal>(finalSize.height());

  qreal deltaX = 2.0 * xMin * scaleX;
  qreal deltaY = 2.0 * yMin * scaleY;

  int nbX = finalSize.width() / subSize.width();
  int nbY = finalSize.height() / subSize.height();

  // Extra subimage on the right/bottom border(s) if needed
  if (nbX * subSize.width() < finalSize.width())
    nbX++;
  if (nbY * subSize.height() < finalSize.height())
    nbY++;
  GLdouble frustum[6];
  camera()->getFrustum(frustum);
  QOpenGLFramebufferObject fbo(size,QOpenGLFramebufferObject::CombinedDepthStencil, GL_TEXTURE_2D, GL_RGBA32F);
  stored_fbo = &fbo;
  for (int i=0; i<nbX; i++)
    for (int j=0; j<nbY; j++)
    {
      fbo.bind();
      glClearColor(GLfloat(backgroundColor().redF()),
                   GLfloat(backgroundColor().greenF()),
                   GLfloat(backgroundColor().blueF()),
                   alpha);
      double frustum[6];
      frustum[0]= -xMin + i*deltaX;
      frustum[1]= -xMin + (i+1)*deltaX ;
      frustum[2]= yMin - j*deltaY;
      frustum[3]= yMin - (j+1)*deltaY;
      frustum[4]= zNear;
      frustum[5]= zFar;
      camera()->setFrustum(frustum);
      preDraw();
      draw();
      fbo.release();

      QImage snapshot = fbo.toImage();
      QImage subImage = snapshot.scaled(subSize, ::Qt::IgnoreAspectRatio, ::Qt::SmoothTransformation);
      // Copy subImage in image
      for (int ii=0; ii<subSize.width(); ii++)
      {
        int fi = i*subSize.width() + ii;
        if (fi == image->width())
          break;
        for (int jj=0; jj<subSize.height(); jj++)
        {
          int fj = j*subSize.height() + jj;
          if (fj == image->height())
            break;
          image->setPixel(fi, fj, subImage.pixel(ii,jj));
        }
      }
    }
  if(background_color !=0)
    setBackgroundColor(previousBGColor);
  camera()->setFrustum(frustum);
  stored_fbo = nullptr;
  return image;
}

CGAL_INLINE_FUNCTION
QOpenGLFramebufferObject* CGAL::QGLViewer::getStoredFrameBuffer() const
{
  return stored_fbo;
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::setStoredFrameBuffer(QOpenGLFramebufferObject *fbo)
{
  stored_fbo = fbo;
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::saveSnapshot()
{
  qreal aspectRatio = width() / static_cast<qreal>(height());
  static ImageInterface* imageInterface = nullptr;

  if (!imageInterface)
    imageInterface = new ImageInterface(this, aspectRatio);

  imageInterface->imgWidth->setValue(width());
  imageInterface->imgHeight->setValue(height());

  if (imageInterface->exec() == QDialog::Rejected)
    return;
  QSize finalSize(imageInterface->imgWidth->value(), imageInterface->imgHeight->value());
  bool expand = imageInterface->expandFrustum->isChecked();
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save Snapshot"), "", tr("Image Files (*.png *.jpg *.bmp)"));
  if(fileName.isEmpty())
  {
    return;
  }
  QImage* image= takeSnapshot(qglviewer::SnapShotBackground(imageInterface->color_comboBox->currentIndex()),
        finalSize, imageInterface->oversampling->value(), expand);
  if(image)
  {
    image->save(fileName);
    delete image;
  }
}

CGAL_INLINE_FUNCTION
void CGAL::QGLViewer::saveSnapshot(const QString& fileName,
                                   const qreal finalWidth, const qreal finalHeight,
                                   const bool expand,
                                   const double oversampling,
                                   qglviewer::SnapShotBackground background_color)
{
  if(fileName.isEmpty())
    return;

  QSize finalSize(finalWidth, finalHeight);

  QImage* image = takeSnapshot(qglviewer::SnapShotBackground(background_color),
                               finalSize, oversampling, expand);
  if(image)
  {
    image->save(fileName);
    delete image;
  }
}

} // namespace CGAL

CGAL_INLINE_FUNCTION
bool CGAL::QGLViewer::isSharing() const
{
  return is_sharing;
}
