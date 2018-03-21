/****************************************************************************

 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of the QGLViewer library version 2.7.0.

 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.
 In addition, as a special exception, Gilles Debunne gives you certain 
 additional rights, described in the file GPL_EXCEPTION in this package.

 libQGLViewer uses dual licensing. Commercial/proprietary software must
 purchase a libQGLViewer Commercial License.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/

#include <CGAL/Qt/qglviewer.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Qt/domUtils.h>
#include <CGAL/Qt/keyFrameInterpolator.h>
#include <CGAL/Qt/manipulatedCameraFrame.h>

#include <QApplication>
#include <QDateTime>
#include <QDir>
#include <QFileInfo>
#include <QGLContext>
#include <QImage>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPushButton>
#include <QTabWidget>
#include <QTextEdit>
#include <QTextStream>
#include <QTimer>
#include <QUrl>
#include <QtAlgorithms>

using namespace std;
using namespace qglviewer;

// Static private variable
QList<QGLViewer *> QGLViewer::QGLViewerPool_;

/*! \mainpage

libQGLViewer is a free C++ library based on Qt that enables the quick creation
of OpenGL 3D viewers. It features a powerful camera trackball and simple
applications simply require an implementation of the <code>draw()</code> method.
This makes it a tool of choice for OpenGL beginners and assignments. It provides
screenshot saving, mouse manipulated frames, stereo display, interpolated
keyFrames, object selection, and much more. It is fully
customizable and easy to extend to create complex applications, with a possible
Qt GUI.

libQGLViewer is <i>not</i> a 3D viewer that can be used directly to view 3D
scenes in various formats. It is more likely to be the starting point for the
coding of such a viewer.

libQGLViewer is based on the Qt toolkit and hence compiles on any architecture
(Unix-Linux, Mac, Windows, ...). Full reference documentation and many examples
are provided.

See the project main page for details on the project and installation steps. */

void QGLViewer::defaultConstructor() {
  int poolIndex = QGLViewer::QGLViewerPool_.indexOf(NULL);
  setFocusPolicy(Qt::StrongFocus);

  if (poolIndex >= 0)
    QGLViewer::QGLViewerPool_.replace(poolIndex, this);
  else
    QGLViewer::QGLViewerPool_.append(this);
  camera_ = new Camera();
  setCamera(camera());

  setDefaultShortcuts();
  setDefaultMouseBindings();

  setSnapshotFileName(tr("snapshot", "Default snapshot file name"));
  initializeSnapshotFormats();
  setSnapshotCounter(0);
  setSnapshotQuality(95);

  fpsTime_.start();
  fpsCounter_ = 0;
  f_p_s_ = 0.0;
  fpsString_ = tr("%1Hz", "Frames per seconds, in Hertz").arg("?");
  visualHint_ = 0;
  previousPathId_ = 0;
  // prevPos_ is not initialized since pos() is not meaningful here.
  // It will be set when setFullScreen(false) is called after
  // setFullScreen(true)

  // #CONNECTION# default values in initFromDOMElement()
  manipulatedFrame_ = NULL;
  manipulatedFrameIsACamera_ = false;
  mouseGrabberIsAManipulatedFrame_ = false;
  mouseGrabberIsAManipulatedCameraFrame_ = false;
  displayMessage_ = false;
  connect(&messageTimer_, SIGNAL(timeout()), SLOT(hideMessage()));
  messageTimer_.setSingleShot(true);
  helpWidget_ = NULL;
  setMouseGrabber(NULL);

  setSceneRadius(1.0);
  showEntireScene();
  setStateFileName(".qglviewer.xml");

  // #CONNECTION# default values in initFromDOMElement()
  setAxisIsDrawn(false);
  setGridIsDrawn(false);
  setFPSIsDisplayed(false);
  setCameraIsEdited(false);
  setTextIsEnabled(true);
  setStereoDisplay(false);
  // Make sure move() is not called, which would call initializeGL()
  fullScreen_ = false;
  setFullScreen(false);

  animationTimerId_ = 0;
  stopAnimation();
  setAnimationPeriod(40); // 25Hz

  selectBuffer_ = NULL;
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
  currentlyPressedKey_ = Qt::Key(0);

  setAttribute(Qt::WA_NoSystemBackground);
  axisIsDrawn_ = true;

  tileRegion_ = NULL;
}

#ifndef DOXYGEN
/*! These contructors are deprecated since version 2.7.0, since they are not
 * supported by QOpenGlWidget */

/*! Constructor. See \c QGLWidget documentation for details.

All viewer parameters (display flags, scene parameters, associated objects...)
are set to their default values. See the associated documentation.

If the \p shareWidget parameter points to a valid \c QGLWidget, the QGLViewer
will share the OpenGL context with \p shareWidget (see isSharing()). */
QGLViewer::QGLViewer(QWidget *parent,
                     Qt::WindowFlags flags)
    : QOpenGLWidget(parent, flags) {
  defaultConstructor();
}

QGLViewer::QGLViewer(QGLContext*,
                     QWidget *parent,
                     Qt::WindowFlags flags)
  : QOpenGLWidget(parent, flags) {
  defaultConstructor();
}
#endif // DOXYGEN

/*! Virtual destructor.

The viewer is replaced by \c NULL in the QGLViewerPool() (in order to preserve
other viewer's indexes) and allocated memory is released. The camera() is
deleted and should be copied before if it is shared by an other viewer. */
QGLViewer::~QGLViewer() {
  // See closeEvent comment. Destructor is called (and not closeEvent) only when
  // the widget is embedded. Hence we saveToFile here. It is however a bad idea
  // if virtual domElement() has been overloaded ! if (parent())
  // saveStateToFileForAllViewers();

  QGLViewer::QGLViewerPool_.replace(QGLViewer::QGLViewerPool_.indexOf(this),
                                    NULL);

  delete camera();
  delete[] selectBuffer_;
  if (helpWidget()) {
    // Needed for Qt 4 which has no main widget.
    helpWidget()->close();
    delete helpWidget_;
  }
}

static QString QGLViewerVersionString() {
  return QString::number((QGLVIEWER_VERSION & 0xff0000) >> 16) + "." +
         QString::number((QGLVIEWER_VERSION & 0x00ff00) >> 8) + "." +
         QString::number(QGLVIEWER_VERSION & 0x0000ff);
}

static Qt::KeyboardModifiers keyboardModifiersFromState(unsigned int state) {
  // Convertion of keyboard modifiers and mouse buttons as an int is no longer
  // supported : emulate
  return Qt::KeyboardModifiers(int(state & 0xFF000000));
}

static Qt::MouseButton mouseButtonFromState(unsigned int state) {
  // Convertion of keyboard modifiers and mouse buttons as an int is no longer
  // supported : emulate
  return Qt::MouseButton(state & 0xFFFF);
}

/*! Initializes the QGLViewer OpenGL context and then calls user-defined init().

This method is automatically called once, before the first call to paintGL().

Overload init() instead of this method to modify viewer specific OpenGL state or
to create display lists.

To make beginners' life easier and to simplify the examples, this method
slightly modifies the standard OpenGL state: \code glEnable(GL_LIGHT0);
glEnable(GL_LIGHTING);
glEnable(GL_DEPTH_TEST);
glEnable(GL_COLOR_MATERIAL);
\endcode

If you port an existing application to QGLViewer and your display changes, you
probably want to disable these flags in init() to get back to a standard OpenGL
state. */
void QGLViewer::initializeGL() {
  QOpenGLFunctions_2_1::initializeOpenGLFunctions();
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);

  // Default colors
  setForegroundColor(QColor(180, 180, 180));
  setBackgroundColor(QColor(51, 51, 51));

  // Clear the buffer where we're going to draw
  if (format().stereo()) {
    glDrawBuffer(GL_BACK_RIGHT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawBuffer(GL_BACK_LEFT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  } else
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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
  //program without light
  {
    //Vertex source code
    const char v_s[] =
    {
      "#version 120 \n"
      "attribute highp vec4 vertex;\n"
      "uniform highp mat4 mvp_matrix;\n"
      "void main(void)\n"
      "{\n"
      "   gl_Position = mvp_matrix * vertex; \n"
      "} \n"
      "\n"
    };
    //Fragment source code
    const char f_s[] =
    {
      "#version 120 \n"
      "uniform highp vec4 color; \n"
      "void main(void) { \n"
      "gl_FragColor = color; \n"
      "} \n"
      "\n"
    };
    
    //It is said in the doc that a QOpenGLShader is 
    // only destroyed with the QOpenGLShaderProgram 
    //it has been linked with.
    QOpenGLShader vertex_shader(QOpenGLShader::Vertex);
    if(!vertex_shader.compileSourceCode(v_s))
    {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }
    
    QOpenGLShader fragment_shader(QOpenGLShader::Fragment);
    if(!fragment_shader.compileSourceCode(f_s))
    {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
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
  //program with light
  {
    //Vertex source code
    const char vertex_source[] =
    {
      "#version 120 \n"
      "attribute highp vec4 vertex;\n"
      "attribute highp vec3 normal;\n"
      "attribute highp vec4 colors;\n"
      "uniform highp mat4 mvp_matrix;\n"
      "uniform highp mat4 mv_matrix; \n"
      "varying highp vec4 fP; \n"
      "varying highp vec3 fN; \n"
      "varying highp vec4 color; \n"
      "void main(void)\n"
      "{\n"
      "   color = vec4(colors.xyz, 1.0f); \n"
      "   fP = mv_matrix * vertex; \n"
      "   fN = mat3(mv_matrix)* normal; \n"
      "   gl_Position = vec4(mvp_matrix * vertex); \n"
      "} \n"
      "\n"
    };
    //Fragment source code
    const char fragment_source[] =
    {
      "#version 120 \n"
      "varying highp vec4 color; \n"
      "varying highp vec4 fP; \n"
      "varying highp vec3 fN; \n"  
      "void main(void) { \n"
      "   highp vec4 light_pos = vec4(0.0f, 0.0f, 1.0f, 1.0f);  \n"
      "   highp vec4 light_diff = vec4(1.0f, 1.0f, 1.0f, 1.0f); \n"
      "   highp vec4 light_spec = vec4(0.0f, 0.0f, 0.0f, 1.0f); \n"
      "   highp vec4 light_amb = vec4(0.4f, 0.4f, 0.4f, 0.4f);  \n"
      "   highp float spec_power = 51.8f ; \n"
      "   vec3 L = light_pos.xyz - fP.xyz; \n"
      "   vec3 V = -fP.xyz; \n"
      "   vec3 N; \n"
      "   if(fN == vec3(0.0,0.0,0.0)) \n"
      "       N = vec3(0.0,0.0,0.0); \n"
      "   else \n"
      "       N = normalize(fN); \n"
      "   L = normalize(L); \n"
      "   V = normalize(V); \n"
      "   vec3 R = reflect(-L, N); \n"
      "   vec4 diffuse = max(abs(dot(N,L)),0.0) * light_diff*color; \n"
      "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"
      
      "gl_FragColor = color*light_amb + diffuse + specular; \n"
      "gl_FragColor = vec4(gl_FragColor.xyz, 1.0f); \n"
      "} \n"
      "\n"
    };
    
    //It is said in the doc that a QOpenGLShader is 
    // only destroyed with the QOpenGLShaderProgram 
    //it has been linked with.
    QOpenGLShader vertex_shader(QOpenGLShader::Vertex);
    if(!vertex_shader.compileSourceCode(vertex_source))
    {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }
    
    QOpenGLShader fragment_shader(QOpenGLShader::Fragment);
    if(!fragment_shader.compileSourceCode(fragment_source))
    {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }
    
    if(!rendering_program_light.addShader(&vertex_shader))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_light.addShader(&fragment_shader))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_light.link())
    {
      qDebug() << rendering_program_light.log();
    }
  }
  
  // Calls user defined method. Default emits a signal.
  init();

  // Give time to glInit to finish and then call setFullScreen().
  if (isFullScreen())
    QTimer::singleShot(100, this, SLOT(delayedFullScreen()));
}

/*! Main paint method, inherited from \c QOpenGLWidget.

Calls the following methods, in that order:
\arg preDraw() (or preDrawStereo() if viewer displaysInStereo()) : places the
camera in the world coordinate system. \arg draw() (or fastDraw() when the
camera is manipulated) : main drawing method. Should be overloaded. \arg
postDraw() : display of visual hints (world axis, FPS...) */
void QGLViewer::paintGL() {
  if (displaysInStereo()) {
    for (int view = 1; view >= 0; --view) {
      // Clears screen, set model view matrix with shifted matrix for ith buffer
      preDrawStereo(view);
      // Used defined method. Default is empty
      if (camera()->frame()->isManipulated())
        fastDraw();
      else
        draw();
      postDraw();
    }
  } else {
    // Clears screen, set model view matrix...
    preDraw();
    // Used defined method. Default calls draw()
    if (camera()->frame()->isManipulated())
      fastDraw();
    else
      draw();
    // Add visual hints: axis, camera, grid...
    postDraw();
  }
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
void QGLViewer::preDraw() {
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
void QGLViewer::postDraw() {
  // Pivot point, line when camera rolls, zoom region

  if (gridIsDrawn()) {
    glLineWidth(1.0);
    drawGrid(camera()->sceneRadius());
  }
  if (axisIsDrawn()) {
    glLineWidth(2.0);
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

  // Restore foregroundColor
  float color[4];
  color[0] = foregroundColor().red() / 255.0f;
  color[1] = foregroundColor().green() / 255.0f;
  color[2] = foregroundColor().blue() / 255.0f;
  color[3] = 1.0f;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  if (FPSIsDisplayed())
    displayFPS();
  if (displayMessage_)
    drawText(10, height() - 10, message_);

  // Restore GL state
  glPopAttrib();
  glPopMatrix();
}

/*! Called before draw() (instead of preDraw()) when viewer displaysInStereo().

Same as preDraw() except that the glDrawBuffer() is set to \c GL_BACK_LEFT or \c
GL_BACK_RIGHT depending on \p leftBuffer, and it uses
qglviewer::Camera::loadProjectionMatrixStereo() and
qglviewer::Camera::loadModelViewMatrixStereo() instead. */
void QGLViewer::preDrawStereo(bool leftBuffer) {
  // Set buffer to draw in
  // Seems that SGI and Crystal Eyes are not synchronized correctly !
  // That's why we don't draw in the appropriate buffer...
  if (!leftBuffer)
    glDrawBuffer(GL_BACK_LEFT);
  else
    glDrawBuffer(GL_BACK_RIGHT);

  // Clear the buffer where we're going to draw
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // GL_PROJECTION matrix
  camera()->loadProjectionMatrixStereo(leftBuffer);
  // GL_MODELVIEW matrix
  camera()->loadModelViewMatrixStereo(leftBuffer);

  Q_EMIT drawNeeded();
}

/*! Draws a simplified version of the scene to guarantee interactive camera
displacements.

This method is called instead of draw() when the qglviewer::Camera::frame() is
qglviewer::ManipulatedCameraFrame::isManipulated(). Default implementation
simply calls draw().

Overload this method if your scene is too complex to allow for interactive
camera manipulation. See the <a href="../examples/fastDraw.html">fastDraw
example</a> for an illustration. */
void QGLViewer::fastDraw() { draw(); }

/*! Starts (\p edit = \c true, default) or stops (\p edit=\c false) the edition
of the camera().

Current implementation is limited to paths display. Get current state using
cameraIsEdited().

\attention This method sets the qglviewer::Camera::zClippingCoefficient() to 5.0
when \p edit is \c true, so that the Camera paths (see
qglviewer::Camera::keyFrameInterpolator()) are not clipped. It restores the
previous value when \p edit is \c false. */
void QGLViewer::setCameraIsEdited(bool edit) {
  cameraIsEdited_ = edit;
  if (edit) {
    previousCameraZClippingCoefficient_ = camera()->zClippingCoefficient();
    // #CONNECTION# 5.0 also used in domElement() and in initFromDOMElement().
    camera()->setZClippingCoefficient(5.0);
  } else
    camera()->setZClippingCoefficient(previousCameraZClippingCoefficient_);

  Q_EMIT cameraIsEditedChanged(edit);

  update();
}

// Key bindings. 0 means not defined
void QGLViewer::setDefaultShortcuts() {
  // D e f a u l t   a c c e l e r a t o r s
  setShortcut(DRAW_AXIS, Qt::Key_A);
  setShortcut(DRAW_GRID, Qt::Key_G);
  setShortcut(DISPLAY_FPS, Qt::Key_F);
  setShortcut(ENABLE_TEXT, Qt::SHIFT + Qt::Key_Question);
  setShortcut(EXIT_VIEWER, Qt::Key_Escape);
  setShortcut(SAVE_SCREENSHOT, Qt::CTRL + Qt::Key_S);
  setShortcut(CAMERA_MODE, Qt::Key_Space);
  setShortcut(FULL_SCREEN, Qt::ALT + Qt::Key_Return);
  setShortcut(STEREO, Qt::Key_S);
  setShortcut(ANIMATION, Qt::Key_Return);
  setShortcut(HELP, Qt::Key_H);
  setShortcut(EDIT_CAMERA, Qt::Key_C);
  setShortcut(MOVE_CAMERA_LEFT, Qt::Key_Left);
  setShortcut(MOVE_CAMERA_RIGHT, Qt::Key_Right);
  setShortcut(MOVE_CAMERA_UP, Qt::Key_Up);
  setShortcut(MOVE_CAMERA_DOWN, Qt::Key_Down);
  setShortcut(INCREASE_FLYSPEED, Qt::Key_Plus);
  setShortcut(DECREASE_FLYSPEED, Qt::Key_Minus);
  setShortcut(SNAPSHOT_TO_CLIPBOARD, Qt::CTRL + Qt::Key_C);

  keyboardActionDescription_[DISPLAY_FPS] =
      tr("Toggles the display of the FPS", "DISPLAY_FPS action description");
  keyboardActionDescription_[SAVE_SCREENSHOT] =
      tr("Saves a screenshot", "SAVE_SCREENSHOT action description");
  keyboardActionDescription_[FULL_SCREEN] =
      tr("Toggles full screen display", "FULL_SCREEN action description");
  keyboardActionDescription_[DRAW_AXIS] = tr(
      "Toggles the display of the world axis", "DRAW_AXIS action description");
  keyboardActionDescription_[DRAW_GRID] =
      tr("Toggles the display of the XY grid", "DRAW_GRID action description");
  keyboardActionDescription_[CAMERA_MODE] = tr(
      "Changes camera mode (observe or fly)", "CAMERA_MODE action description");
  keyboardActionDescription_[STEREO] =
      tr("Toggles stereo display", "STEREO action description");
  keyboardActionDescription_[HELP] =
      tr("Opens this help window", "HELP action description");
  keyboardActionDescription_[ANIMATION] =
      tr("Starts/stops the animation", "ANIMATION action description");
  keyboardActionDescription_[EDIT_CAMERA] =
      tr("Toggles camera paths display",
         "EDIT_CAMERA action description"); // TODO change
  keyboardActionDescription_[ENABLE_TEXT] =
      tr("Toggles the display of the text", "ENABLE_TEXT action description");
  keyboardActionDescription_[EXIT_VIEWER] =
      tr("Exits program", "EXIT_VIEWER action description");
  keyboardActionDescription_[MOVE_CAMERA_LEFT] =
      tr("Moves camera left", "MOVE_CAMERA_LEFT action description");
  keyboardActionDescription_[MOVE_CAMERA_RIGHT] =
      tr("Moves camera right", "MOVE_CAMERA_RIGHT action description");
  keyboardActionDescription_[MOVE_CAMERA_UP] =
      tr("Moves camera up", "MOVE_CAMERA_UP action description");
  keyboardActionDescription_[MOVE_CAMERA_DOWN] =
      tr("Moves camera down", "MOVE_CAMERA_DOWN action description");
  keyboardActionDescription_[INCREASE_FLYSPEED] =
      tr("Increases fly speed", "INCREASE_FLYSPEED action description");
  keyboardActionDescription_[DECREASE_FLYSPEED] =
      tr("Decreases fly speed", "DECREASE_FLYSPEED action description");
  keyboardActionDescription_[SNAPSHOT_TO_CLIPBOARD] =
      tr("Copies a snapshot to clipboard",
         "SNAPSHOT_TO_CLIPBOARD action description");

  // K e y f r a m e s   s h o r t c u t   k e y s
  setPathKey(Qt::Key_F1, 1);
  setPathKey(Qt::Key_F2, 2);
  setPathKey(Qt::Key_F3, 3);
  setPathKey(Qt::Key_F4, 4);
  setPathKey(Qt::Key_F5, 5);
  setPathKey(Qt::Key_F6, 6);
  setPathKey(Qt::Key_F7, 7);
  setPathKey(Qt::Key_F8, 8);
  setPathKey(Qt::Key_F9, 9);
  setPathKey(Qt::Key_F10, 10);
  setPathKey(Qt::Key_F11, 11);
  setPathKey(Qt::Key_F12, 12);

  setAddKeyFrameKeyboardModifiers(Qt::AltModifier);
  setPlayPathKeyboardModifiers(Qt::NoModifier);
}

// M o u s e   b e h a v i o r
void QGLViewer::setDefaultMouseBindings() {
  const Qt::KeyboardModifiers cameraKeyboardModifiers = Qt::NoModifier;
  const Qt::KeyboardModifiers frameKeyboardModifiers = Qt::ControlModifier;

  //#CONNECTION# toggleCameraMode()
  for (int handler = 0; handler < 2; ++handler) {
    MouseHandler mh = (MouseHandler)(handler);
    Qt::KeyboardModifiers modifiers =
        (mh == FRAME) ? frameKeyboardModifiers : cameraKeyboardModifiers;

    setMouseBinding(modifiers, Qt::LeftButton, mh, ROTATE);
    setMouseBinding(modifiers, Qt::MidButton, mh, ZOOM);
    setMouseBinding(modifiers, Qt::RightButton, mh, TRANSLATE);

    setMouseBinding(Qt::Key_R, modifiers, Qt::LeftButton, mh, SCREEN_ROTATE);

    setWheelBinding(modifiers, mh, ZOOM);
  }

  // Z o o m   o n   r e g i o n
  setMouseBinding(Qt::ShiftModifier, Qt::MidButton, CAMERA, ZOOM_ON_REGION);

  // S e l e c t
  setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, SELECT);

  setMouseBinding(Qt::ShiftModifier, Qt::RightButton, RAP_FROM_PIXEL);
  // D o u b l e   c l i c k
  setMouseBinding(Qt::NoModifier, Qt::LeftButton, ALIGN_CAMERA, true);
  setMouseBinding(Qt::NoModifier, Qt::MidButton, SHOW_ENTIRE_SCENE, true);
  setMouseBinding(Qt::NoModifier, Qt::RightButton, CENTER_SCENE, true);

  setMouseBinding(frameKeyboardModifiers, Qt::LeftButton, ALIGN_FRAME, true);
  // middle double click makes no sense for manipulated frame
  setMouseBinding(frameKeyboardModifiers, Qt::RightButton, CENTER_FRAME, true);

  // A c t i o n s   w i t h   k e y   m o d i f i e r s
  setMouseBinding(Qt::Key_Z, Qt::NoModifier, Qt::LeftButton, ZOOM_ON_PIXEL);
  setMouseBinding(Qt::Key_Z, Qt::NoModifier, Qt::RightButton, ZOOM_TO_FIT);

#ifdef Q_OS_MAC
  // Specific Mac bindings for touchpads. Two fingers emulate a wheelEvent which
  // zooms. There is no right button available : make Option key + left emulate
  // the right button. A Control+Left indeed emulates a right click (OS X system
  // configuration), but it does no seem to support dragging. Done at the end to
  // override previous settings.
  const Qt::KeyboardModifiers macKeyboardModifiers = Qt::AltModifier;

  setMouseBinding(macKeyboardModifiers, Qt::LeftButton, CAMERA, TRANSLATE);
  setMouseBinding(macKeyboardModifiers, Qt::LeftButton, CENTER_SCENE, true);
  setMouseBinding(frameKeyboardModifiers | macKeyboardModifiers, Qt::LeftButton,
                  CENTER_FRAME, true);
  setMouseBinding(frameKeyboardModifiers | macKeyboardModifiers, Qt::LeftButton,
                  FRAME, TRANSLATE);
#endif
}

/*! Associates a new qglviewer::Camera to the viewer.

You should only use this method when you derive a new class from
qglviewer::Camera and want to use one of its instances instead of the original
class.

It you simply want to save and restore Camera positions, use
qglviewer::Camera::addKeyFrameToPath() and qglviewer::Camera::playPath()
instead.

This method silently ignores \c NULL \p camera pointers. The calling method is
responsible for deleting the previous camera pointer in order to prevent memory
leaks if needed.

The sceneRadius() and sceneCenter() of \p camera are set to the \e current
QGLViewer values.

All the \p camera qglviewer::Camera::keyFrameInterpolator()
qglviewer::KeyFrameInterpolator::interpolated() signals are connected to the
viewer update() slot. The connections with the previous viewer's camera are
removed. */
void QGLViewer::setCamera(Camera *const camera) {
  if (!camera)
    return;

  camera->setSceneRadius(sceneRadius());
  camera->setSceneCenter(sceneCenter());
  camera->setScreenWidthAndHeight(width(), height());

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

void QGLViewer::connectAllCameraKFIInterpolatedSignals(bool connection) {
  for (QMap<unsigned int, KeyFrameInterpolator *>::ConstIterator
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
void QGLViewer::drawLight(GLenum light, qreal scale) const {
}

#if (QT_VERSION >= QT_VERSION_CHECK(5, 4, 0))
void QGLViewer::renderText(int x, int y, const QString &str,
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

void QGLViewer::renderText(double x, double y, double z, const QString &str,
                           const QFont &font) {
  const Vec proj = camera_->projectedCoordinatesOf(Vec(x, y, z));
  renderText(proj.x, proj.y, str, font);
}
#endif

/*! Draws \p text at position \p x, \p y (expressed in screen coordinates
pixels, origin in the upper left corner of the widget).

The default QApplication::font() is used to render the text when no \p fnt is
specified. Use QApplication::setFont() to define this default font.

You should disable \c GL_LIGHTING and \c GL_DEPTH_TEST before this method so
that colors are properly rendered.

This method can be used in conjunction with the
qglviewer::Camera::projectedCoordinatesOf() method to display a text attached to
an object. In your draw() method use: \code qglviewer::Vec screenPos =
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
void QGLViewer::drawText(int x, int y, const QString &text, const QFont &fnt) {
  if (!textIsEnabled())
    return;

  if (tileRegion_ != NULL) {
    renderText(int((x - tileRegion_->xMin) * width() /
                   (tileRegion_->xMax - tileRegion_->xMin)),
               int((y - tileRegion_->yMin) * height() /
                   (tileRegion_->yMax - tileRegion_->yMin)),
               text, scaledFont(fnt));
  } else
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
void QGLViewer::displayMessage(const QString &message, int delay) {
  message_ = message;
  displayMessage_ = true;
  // Was set to single shot in defaultConstructor.
  messageTimer_.start(delay);
  if (textIsEnabled())
    update();
}

void QGLViewer::hideMessage() {
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
void QGLViewer::displayFPS() {
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
qglviewer::Camera::projectedCoordinatesOf() to compute the 2D projection on
screen of a 3D point (see the <a
href="../examples/screenCoordSystem.html">screenCoordSystem</a> example). See
also drawText().

In this mode, you should use z values that are in the [0.0, 1.0[ range (0.0
corresponding to the near clipping plane and 1.0 being just beyond the far
clipping plane). This interval matches the values that can be read from the
z-buffer. Note that if you use the convenient \c glVertex2i() to provide
coordinates, the implicit 0.0 z coordinate will make your drawings appear \e on
\e top of the rest of the scene. */
void QGLViewer::startScreenCoordinatesSystem(bool ) const {
}

/*! Stops the pixel coordinate drawing block started by
startScreenCoordinatesSystem().

The \c GL_MODELVIEW and \c GL_PROJECTION matrices modified in
startScreenCoordinatesSystem() are restored. \c glMatrixMode is set to \c
GL_MODELVIEW. */
void QGLViewer::stopScreenCoordinatesSystem() const {
}

/*! Overloading of the \c QObject method.

If animationIsStarted(), calls animate() and draw(). */
void QGLViewer::timerEvent(QTimerEvent *) {
  if (animationIsStarted()) {
    animate();
    update();
  }
}

/*! Starts the animation loop. See animationIsStarted(). */
void QGLViewer::startAnimation() {
  animationTimerId_ = startTimer(animationPeriod());
  animationStarted_ = true;
}

/*! Stops animation. See animationIsStarted(). */
void QGLViewer::stopAnimation() {
  animationStarted_ = false;
  if (animationTimerId_ != 0)
    killTimer(animationTimerId_);
}

/*! Overloading of the \c QWidget method.

Saves the viewer state using saveStateToFile() and then calls
QOpenGLWidget::closeEvent(). */
void QGLViewer::closeEvent(QCloseEvent *e) {
  // When the user clicks on the window close (x) button:
  // - If the viewer is a top level window, closeEvent is called and then saves
  // to file. - Otherwise, nothing happen s:( When the user press the
  // EXIT_VIEWER keyboard shortcut: - If the viewer is a top level window,
  // saveStateToFile() is also called - Otherwise, closeEvent is NOT called and
  // keyPressEvent does the job.

  /* After tests:
  E : Embedded widget
  N : Widget created with new
  C : closeEvent called
  D : destructor called

  E	N	C	D
  y	y
  y	n		y
  n	y	y
  n	n	y	y

  closeEvent is called iif the widget is NOT embedded.

  Destructor is called iif the widget is created on the stack
  or if widget (resp. parent if embedded) is created with WDestructiveClose
  flag.

  closeEvent always before destructor.

  Close using qApp->closeAllWindows or (x) is identical.
  */

  // #CONNECTION# Also done for EXIT_VIEWER in keyPressEvent().
  saveStateToFile();
  QOpenGLWidget::closeEvent(e);
}

/*! Simple wrapper method: calls \c select(event->pos()).

Emits \c pointSelected(e) which is useful only if you rely on the Qt signal-slot
mechanism and you did not overload QGLViewer. If you choose to derive your own
viewer class, simply overload select() (or probably simply drawWithNames(), see
the <a href="../examples/select.html">select example</a>) to implement your
selection mechanism.

This method is called when you use the QGLViewer::SELECT mouse binding(s)
(default is Shift + left button). Use setMouseBinding() to change this. */
void QGLViewer::select(const QMouseEvent *event) {
  // For those who don't derive but rather rely on the signal-slot mechanism.
  Q_EMIT pointSelected(event);
  select(event->pos());
}

/*! This method performs a selection in the scene from pixel coordinates.

It is called when the user clicks on the QGLViewer::SELECT
QGLViewer::ClickAction binded button(s) (default is Shift + LeftButton).

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
region. Use qglviewer::Camera::convertClickToLine() to transform these
coordinates in a 3D ray if you want to perform an analytical intersection.

\attention \c GL_SELECT mode seems to report wrong results when used in
conjunction with backface culling. If you encounter problems try to \c
glDisable(GL_CULL_FACE). */
void QGLViewer::select(const QPoint &point) {
  beginSelection(point);
  drawWithNames();
  endSelection(point);
  postSelection(point);
}

/*! This method should prepare the selection. It is called by select() before
drawWithNames().
The default function is empty.
*/
void QGLViewer::beginSelection(const QPoint &) 
{
  
}

/*! This method is called by select() after scene elements were drawn by
drawWithNames(). It should analyze the selection result to determine which
object is actually selected.

The default implementation relies on \c GL_SELECT mode (see beginSelection()).
It assumes that names were pushed and popped in drawWithNames(), and analyzes
the selectBuffer() to find the name that corresponds to the closer (z min)
object. It then setSelectedName() to this value, or to -1 if the selectBuffer()
is empty (no object drawn in selection region). Use selectedName() (probably in
the postSelection() method) to retrieve this value and update your data
structure accordingly.

This default implementation, although sufficient for many cases is however
limited and you may have to overload this method. This will be the case if
drawWithNames() uses several push levels in the name heap. A more precise depth
selection, for instance privileging points over edges and triangles to avoid z
precision problems, will also require an overloading. A typical implementation
will look like:
\code
glFlush();

// Get the number of objects that were seen through the pick matrix frustum.
// Resets GL_RENDER mode.
GLint nbHits = glRenderMode(GL_RENDER);

if (nbHits <= 0)
setSelectedName(-1);
else
{
// Interpret results: each object created values in the selectBuffer().
// See the glSelectBuffer() man page for details on the buffer structure.
// The following code depends on your selectBuffer() structure.
for (int i=0; i<nbHits; ++i)
if ((selectBuffer())[i*4+1] < zMin)
setSelectedName((selectBuffer())[i*4+3])
}
\endcode

See the <a href="../examples/multiSelect.html">multiSelect example</a> for
a multi-object selection implementation of this method. */
void QGLViewer::endSelection(const QPoint &point) {
  Q_UNUSED(point);

  // Flush GL buffers
  glFlush();

  // Get the number of objects that were seen through the pick matrix frustum.
  // Reset GL_RENDER mode.
  GLint nbHits = glRenderMode(GL_RENDER);

  if (nbHits <= 0)
    setSelectedName(-1);
  else {
    // Interpret results: each object created 4 values in the selectBuffer().
    // selectBuffer[4*i+1] is the object minimum depth value, while
    // selectBuffer[4*i+3] is the id pushed on the stack. Of all the objects
    // that were projected in the pick region, we select the closest one (zMin
    // comparison). This code needs to be modified if you use several stack
    // levels. See glSelectBuffer() man page.
    GLuint zMin = (selectBuffer())[1];
    setSelectedName(int((selectBuffer())[3]));
    for (int i = 1; i < nbHits; ++i)
      if ((selectBuffer())[4 * i + 1] < zMin) {
        zMin = (selectBuffer())[4 * i + 1];
        setSelectedName(int((selectBuffer())[4 * i + 3]));
      }
  }
}

/*! Sets the selectBufferSize().

The previous selectBuffer() is deleted and a new one is created. */
void QGLViewer::setSelectBufferSize(int size) {
  if (selectBuffer_)
    delete[] selectBuffer_;
  selectBufferSize_ = size;
  selectBuffer_ = new GLuint[selectBufferSize()];
}

static QString mouseButtonsString(Qt::MouseButtons b) {
  QString result("");
  bool addAmpersand = false;
  if (b & Qt::LeftButton) {
    result += QGLViewer::tr("Left", "left mouse button");
    addAmpersand = true;
  }
  if (b & Qt::MidButton) {
    if (addAmpersand)
      result += " & ";
    result += QGLViewer::tr("Middle", "middle mouse button");
    addAmpersand = true;
  }
  if (b & Qt::RightButton) {
    if (addAmpersand)
      result += " & ";
    result += QGLViewer::tr("Right", "right mouse button");
  }
  return result;
}

void QGLViewer::performClickAction(ClickAction ca, const QMouseEvent *const e) {
  // Note: action that need it should call update().
  switch (ca) {
  // # CONNECTION setMouseBinding prevents adding NO_CLICK_ACTION in
  // clickBinding_ This case should hence not be possible. Prevents unused case
  // warning.
  case NO_CLICK_ACTION:
    break;
  case ZOOM_ON_PIXEL:
    camera()->interpolateToZoomOnPixel(e->pos());
    break;
  case ZOOM_TO_FIT:
    camera()->interpolateToFitScene();
    break;
  case SELECT:
    select(e);
    update();
    break;
  case RAP_FROM_PIXEL:
    if (!camera()->setPivotPointFromPixel(e->pos()))
      camera()->setPivotPoint(sceneCenter());
    setVisualHintsMask(1);
    update();
    break;
  case RAP_IS_CENTER:
    camera()->setPivotPoint(sceneCenter());
    setVisualHintsMask(1);
    update();
    break;
  case CENTER_FRAME:
    if (manipulatedFrame())
      manipulatedFrame()->projectOnLine(camera()->position(),
                                        camera()->viewDirection());
    break;
  case CENTER_SCENE:
    camera()->centerScene();
    break;
  case SHOW_ENTIRE_SCENE:
    camera()->showEntireScene();
    break;
  case ALIGN_FRAME:
    if (manipulatedFrame())
      manipulatedFrame()->alignWithFrame(camera()->frame());
    break;
  case ALIGN_CAMERA:
    Frame *frame = new Frame();
    frame->setTranslation(camera()->pivotPoint());
    camera()->frame()->alignWithFrame(frame, true);
    delete frame;
    break;
  }
}

/*! Overloading of the \c QWidget method.

When the user clicks on the mouse:
\arg if a mouseGrabber() is defined, qglviewer::MouseGrabber::mousePressEvent()
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
void QGLViewer::mousePressEvent(QMouseEvent *e) {
  //#CONNECTION# mouseDoubleClickEvent has the same structure
  //#CONNECTION# mouseString() concatenates bindings description in inverse
  // order.
  ClickBindingPrivate cbp(e->modifiers(), e->button(), false,
                          (Qt::MouseButtons)(e->buttons() & ~(e->button())),
                          currentlyPressedKey_);

  if (clickBinding_.contains(cbp)) {
    performClickAction(clickBinding_[cbp], e);
  } else if (mouseGrabber()) {
    if (mouseGrabberIsAManipulatedFrame_) {
      for (QMap<MouseBindingPrivate, MouseActionPrivate>::ConstIterator
               it = mouseBinding_.begin(),
               end = mouseBinding_.end();
           it != end; ++it)
        if ((it.value().handler == FRAME) && (it.key().button == e->button())) {
          ManipulatedFrame *mf =
              dynamic_cast<ManipulatedFrame *>(mouseGrabber());
          if (mouseGrabberIsAManipulatedCameraFrame_) {
            mf->ManipulatedFrame::startAction(it.value().action,
                                              it.value().withConstraint);
            mf->ManipulatedFrame::mousePressEvent(e, camera());
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
      case CAMERA:
        camera()->frame()->startAction(map.action, map.withConstraint);
        camera()->frame()->mousePressEvent(e, camera());
        break;
      case FRAME:
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
      if (map.action == SCREEN_ROTATE)
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
void Viewer::mousePressEvent(QMouseEvent* e)
{

if ((e->button() == myButton) && (e->modifiers() == myModifiers))
  myMouseBehavior = true;
else
  QGLViewer::mousePressEvent(e);
}

void Viewer::mouseMoveEvent(QMouseEvent *e)
{
if (myMouseBehavior)
  // Use e->x() and e->y() as you want...
else
  QGLViewer::mouseMoveEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
if (myMouseBehavior)
  myMouseBehavior = false;
else
  QGLViewer::mouseReleaseEvent(e);
}
\endcode */
void QGLViewer::mouseMoveEvent(QMouseEvent *e) {
  if (mouseGrabber()) {
    mouseGrabber()->checkIfGrabsMouse(e->x(), e->y(), camera());
    if (mouseGrabber()->grabsMouse())
      if (mouseGrabberIsAManipulatedCameraFrame_)
        (dynamic_cast<ManipulatedFrame *>(mouseGrabber()))
            ->ManipulatedFrame::mouseMoveEvent(e, camera());
      else
        mouseGrabber()->mouseMoveEvent(e, camera());
    else
      setMouseGrabber(NULL);
    update();
  }

  if (!mouseGrabber()) {
    //#CONNECTION# mouseReleaseEvent has the same structure
    if (camera()->frame()->isManipulated()) {
      camera()->frame()->mouseMoveEvent(e, camera());
      // #CONNECTION# manipulatedCameraFrame::mouseMoveEvent specific if at the
      // beginning
      if (camera()->frame()->action_ == ZOOM_ON_REGION)
        update();
    } else // !
        if ((manipulatedFrame()) && (manipulatedFrame()->isManipulated()))
      if (manipulatedFrameIsACamera_)
        manipulatedFrame()->ManipulatedFrame::mouseMoveEvent(e, camera());
      else
        manipulatedFrame()->mouseMoveEvent(e, camera());
    else if (hasMouseTracking()) {
      Q_FOREACH (MouseGrabber *mg, MouseGrabber::MouseGrabberPool()) {
        mg->checkIfGrabsMouse(e->x(), e->y(), camera());
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
void QGLViewer::mouseReleaseEvent(QMouseEvent *e) {
  if (mouseGrabber()) {
    if (mouseGrabberIsAManipulatedCameraFrame_)
      (dynamic_cast<ManipulatedFrame *>(mouseGrabber()))
          ->ManipulatedFrame::mouseReleaseEvent(e, camera());
    else
      mouseGrabber()->mouseReleaseEvent(e, camera());
    mouseGrabber()->checkIfGrabsMouse(e->x(), e->y(), camera());
    if (!(mouseGrabber()->grabsMouse()))
      setMouseGrabber(NULL);
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
void QGLViewer::wheelEvent(QWheelEvent *e) {
  if (mouseGrabber()) {
    if (mouseGrabberIsAManipulatedFrame_) {
      for (QMap<WheelBindingPrivate, MouseActionPrivate>::ConstIterator
               it = wheelBinding_.begin(),
               end = wheelBinding_.end();
           it != end; ++it)
        if (it.value().handler == FRAME) {
          ManipulatedFrame *mf =
              dynamic_cast<ManipulatedFrame *>(mouseGrabber());
          if (mouseGrabberIsAManipulatedCameraFrame_) {
            mf->ManipulatedFrame::startAction(it.value().action,
                                              it.value().withConstraint);
            mf->ManipulatedFrame::wheelEvent(e, camera());
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
      case CAMERA:
        camera()->frame()->startAction(map.action, map.withConstraint);
        camera()->frame()->wheelEvent(e, camera());
        break;
      case FRAME:
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
void QGLViewer::mouseDoubleClickEvent(QMouseEvent *e) {
  //#CONNECTION# mousePressEvent has the same structure
  ClickBindingPrivate cbp(e->modifiers(), e->button(), true,
                          (Qt::MouseButtons)(e->buttons() & ~(e->button())),
                          currentlyPressedKey_);
  if (clickBinding_.contains(cbp))
    performClickAction(clickBinding_[cbp], e);
  else if (mouseGrabber())
    mouseGrabber()->mouseDoubleClickEvent(e, camera());
  else
    e->ignore();
}

/*! Sets the state of displaysInStereo(). See also toggleStereoDisplay().

First checks that the display is able to handle stereovision using
QOpenGLWidget::format(). Opens a warning message box in case of failure. Emits
the stereoChanged() signal otherwise. */
void QGLViewer::setStereoDisplay(bool stereo) {
  if (format().stereo()) {
    stereo_ = stereo;
    if (!displaysInStereo()) {
      glDrawBuffer(GL_BACK_LEFT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glDrawBuffer(GL_BACK_RIGHT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

    Q_EMIT stereoChanged(stereo_);

    update();
  } else if (stereo)
    QMessageBox::warning(this,
                         tr("Stereo not supported", "Message box window title"),
                         tr("Stereo is not supported on this display."));
  else
    stereo_ = false;
}

/*! Sets the isFullScreen() state.

If the QGLViewer is embedded in an other QWidget (see
QWidget::topLevelWidget()), this widget is displayed in full screen instead. */
void QGLViewer::setFullScreen(bool fullScreen) {
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
qglviewer::MouseGrabber::checkIfGrabsMouse() test performed by mouseMoveEvent().

If the MouseGrabber is disabled (see mouseGrabberIsEnabled()), this method
silently does nothing. */
void QGLViewer::setMouseGrabber(MouseGrabber *mouseGrabber) {
  if (!mouseGrabberIsEnabled(mouseGrabber))
    return;

  mouseGrabber_ = mouseGrabber;

  mouseGrabberIsAManipulatedFrame_ =
      (dynamic_cast<ManipulatedFrame *>(mouseGrabber) != NULL);
  mouseGrabberIsAManipulatedCameraFrame_ =
      ((dynamic_cast<ManipulatedCameraFrame *>(mouseGrabber) != NULL) &&
       (mouseGrabber != camera()->frame()));
  Q_EMIT mouseGrabberChanged(mouseGrabber);
}

/*! Sets the mouseGrabberIsEnabled() state. */
void QGLViewer::setMouseGrabberIsEnabled(
    const qglviewer::MouseGrabber *const mouseGrabber, bool enabled) {
  if (enabled)
    disabledMouseGrabbers_.remove(reinterpret_cast<size_t>(mouseGrabber));
  else
    disabledMouseGrabbers_[reinterpret_cast<size_t>(mouseGrabber)];
}

QString QGLViewer::mouseActionString(QGLViewer::MouseAction ma) {
  switch (ma) {
  case QGLViewer::NO_MOUSE_ACTION:
    return QString::null;
  case QGLViewer::ROTATE:
    return QGLViewer::tr("Rotates", "ROTATE mouse action");
  case QGLViewer::ZOOM:
    return QGLViewer::tr("Zooms", "ZOOM mouse action");
  case QGLViewer::TRANSLATE:
    return QGLViewer::tr("Translates", "TRANSLATE mouse action");
  case QGLViewer::MOVE_FORWARD:
    return QGLViewer::tr("Moves forward", "MOVE_FORWARD mouse action");
  case QGLViewer::LOOK_AROUND:
    return QGLViewer::tr("Looks around", "LOOK_AROUND mouse action");
  case QGLViewer::MOVE_BACKWARD:
    return QGLViewer::tr("Moves backward", "MOVE_BACKWARD mouse action");
  case QGLViewer::SCREEN_ROTATE:
    return QGLViewer::tr("Rotates in screen plane",
                         "SCREEN_ROTATE mouse action");
  case QGLViewer::ROLL:
    return QGLViewer::tr("Rolls", "ROLL mouse action");
  case QGLViewer::DRIVE:
    return QGLViewer::tr("Drives", "DRIVE mouse action");
  case QGLViewer::SCREEN_TRANSLATE:
    return QGLViewer::tr("Horizontally/Vertically translates",
                         "SCREEN_TRANSLATE mouse action");
  case QGLViewer::ZOOM_ON_REGION:
    return QGLViewer::tr("Zooms on region for", "ZOOM_ON_REGION mouse action");
  }
  return QString::null;
}

QString QGLViewer::clickActionString(QGLViewer::ClickAction ca) {
  switch (ca) {
  case QGLViewer::NO_CLICK_ACTION:
    return QString::null;
  case QGLViewer::ZOOM_ON_PIXEL:
    return QGLViewer::tr("Zooms on pixel", "ZOOM_ON_PIXEL click action");
  case QGLViewer::ZOOM_TO_FIT:
    return QGLViewer::tr("Zooms to fit scene", "ZOOM_TO_FIT click action");
  case QGLViewer::SELECT:
    return QGLViewer::tr("Selects", "SELECT click action");
  case QGLViewer::RAP_FROM_PIXEL:
    return QGLViewer::tr("Sets pivot point", "RAP_FROM_PIXEL click action");
  case QGLViewer::RAP_IS_CENTER:
    return QGLViewer::tr("Resets pivot point", "RAP_IS_CENTER click action");
  case QGLViewer::CENTER_FRAME:
    return QGLViewer::tr("Centers manipulated frame",
                         "CENTER_FRAME click action");
  case QGLViewer::CENTER_SCENE:
    return QGLViewer::tr("Centers scene", "CENTER_SCENE click action");
  case QGLViewer::SHOW_ENTIRE_SCENE:
    return QGLViewer::tr("Shows entire scene",
                         "SHOW_ENTIRE_SCENE click action");
  case QGLViewer::ALIGN_FRAME:
    return QGLViewer::tr("Aligns manipulated frame",
                         "ALIGN_FRAME click action");
  case QGLViewer::ALIGN_CAMERA:
    return QGLViewer::tr("Aligns camera", "ALIGN_CAMERA click action");
  }
  return QString::null;
}

static QString keyString(unsigned int key) {
#if QT_VERSION >= 0x040100
  return QKeySequence(int(key)).toString(QKeySequence::NativeText);
#else
  return QString(QKeySequence(key));
#endif
}

QString QGLViewer::formatClickActionPrivate(ClickBindingPrivate cbp) {
  bool buttonsBefore = cbp.buttonsBefore != Qt::NoButton;
  QString keyModifierString = keyString(cbp.modifiers + cbp.key);
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
           (cbp.button == Qt::NoButton ? tr("Wheel", "Mouse wheel") : ""))
      .arg(cbp.doubleClick ? tr(" double click", "Suffix after mouse button")
                           : "")
      .arg(buttonsBefore ? tr(" with ", "As in : Left button with Ctrl pressed")
                         : "")
      .arg(buttonsBefore ? mouseButtonsString(cbp.buttonsBefore) : "")
      .arg(buttonsBefore
               ? tr(" pressed", "As in : Left button with Ctrl pressed")
               : "");
}

bool QGLViewer::isValidShortcutKey(int key) {
  return (key >= Qt::Key_Any && key < Qt::Key_Escape) ||
         (key >= Qt::Key_F1 && key <= Qt::Key_F35);
}

#ifndef DOXYGEN
/*! This method is deprecated since version 2.5.0

 Use setMouseBindingDescription(Qt::KeyboardModifiers, Qt::MouseButtons,
 QString, bool, Qt::MouseButtons) instead.
*/
void QGLViewer::setMouseBindingDescription(unsigned int state,
                                           QString description,
                                           bool doubleClick,
                                           Qt::MouseButtons buttonsBefore) {
  qWarning("setMouseBindingDescription(int state,...) is deprecated. Use the "
           "modifier/button equivalent");
  setMouseBindingDescription(keyboardModifiersFromState(state),
                             mouseButtonFromState(state), description,
                             doubleClick, buttonsBefore);
}
#endif

/*! Defines a custom mouse binding description, displayed in the help() window's
 Mouse tab.

 Same as calling setMouseBindingDescription(Qt::Key, Qt::KeyboardModifiers,
 Qt::MouseButton, QString, bool, Qt::MouseButtons), with a key value of
 Qt::Key(0) (i.e. binding description when no regular key needs to be pressed).
 */
void QGLViewer::setMouseBindingDescription(Qt::KeyboardModifiers modifiers,
                                           Qt::MouseButton button,
                                           QString description,
                                           bool doubleClick,
                                           Qt::MouseButtons buttonsBefore) {
  setMouseBindingDescription(Qt::Key(0), modifiers, button, description,
                             doubleClick, buttonsBefore);
}

/*! Defines a custom mouse binding description, displayed in the help() window's
Mouse tab.

\p modifiers is a combination of Qt::KeyboardModifiers (\c Qt::ControlModifier,
\c Qt::AltModifier, \c Qt::ShiftModifier, \c Qt::MetaModifier). Possibly
combined using the \c "|" operator.

\p button is one of the Qt::MouseButtons (\c Qt::LeftButton, \c Qt::MidButton,
\c Qt::RightButton...).

\p doubleClick indicates whether or not the user has to double click this button
to perform the described action. \p buttonsBefore lists the buttons that need to
be pressed before the double click.

Set an empty \p description to \e remove a mouse binding description.

\code
// The R key combined with the Left mouse button rotates the camera in the
screen plane. setMouseBindingDescription(Qt::Key_R, Qt::NoModifier,
Qt::LeftButton, "Rotates camera in screen plane");

// A left button double click toggles full screen
setMouseBindingDescription(Qt::NoModifier, Qt::LeftButton, "Toggles full screen
mode", true);

// Removes the description of Ctrl+Right button
setMouseBindingDescription(Qt::ControlModifier, Qt::RightButton, "");
\endcode

Overload mouseMoveEvent() and friends to implement your custom mouse behavior
(see the mouseMoveEvent() documentation for an example). See the <a
href="../examples/keyboardAndMouse.html">keyboardAndMouse example</a> for an
illustration.

Use setMouseBinding() and setWheelBinding() to change the standard mouse action
bindings. */
void QGLViewer::setMouseBindingDescription(
    Qt::Key key, Qt::KeyboardModifiers modifiers, Qt::MouseButton button,
    QString description, bool doubleClick, Qt::MouseButtons buttonsBefore) {
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
QString QGLViewer::mouseString() const {
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
                            Qt::NoButton, itmb.key().key);

    QString text = mouseActionString(itmb.value().action);

    if (!text.isNull()) {
      switch (itmb.value().handler) {
      case CAMERA:
        text += " " + tr("camera", "Suffix after action");
        break;
      case FRAME:
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
    ClickBindingPrivate cbp(itw.key().modifiers, Qt::NoButton, false,
                            Qt::NoButton, itw.key().key);

    QString text = mouseActionString(itw.value().action);

    if (!text.isNull()) {
      switch (itw.value().handler) {
      case CAMERA:
        text += " " + tr("camera", "Suffix after action");
        break;
      case FRAME:
        text += " " + tr("manipulated frame", "Suffix after action");
        break;
      }
      if (!(itw.value().withConstraint))
        text += "*";
    }

    mouseBinding[cbp] = text;
  }

  for (QMap<ClickBindingPrivate, ClickAction>::ConstIterator
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

The \p key definition is given as an \c int using Qt enumerated values. Set an
empty \p description to remove a shortcut description: \code
setKeyDescription(Qt::Key_W, "Toggles wireframe display");
setKeyDescription(Qt::CTRL+Qt::Key_L, "Loads a new scene");
// Removes a description
setKeyDescription(Qt::CTRL+Qt::Key_C, "");
\endcode

See the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse example</a>
for illustration and the <a href="../keyboard.html">keyboard page</a> for
details. */
void QGLViewer::setKeyDescription(unsigned int key, QString description) {
  if (description.isEmpty())
    keyDescription_.remove(key);
  else
    keyDescription_[key] = description;
}

QString QGLViewer::cameraPathKeysString() const {
  if (pathIndex_.isEmpty())
    return QString::null;

  QVector<Qt::Key> keys;
  keys.reserve(pathIndex_.count());
  for (QMap<Qt::Key, unsigned int>::ConstIterator i = pathIndex_.begin(),
                                                  endi = pathIndex_.end();
       i != endi; ++i)
    keys.push_back(i.key());
  qSort(keys);

  QVector<Qt::Key>::const_iterator it = keys.begin(), end = keys.end();
  QString res = keyString(*it);

  const int maxDisplayedKeys = 6;
  int nbDisplayedKeys = 0;
  Qt::Key previousKey = (*it);
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
setKeyDescription() as well as the \e standard QGLViewer::KeyboardAction
shortcuts (defined using setShortcut()). See the <a
href="../keyboard.html">keyboard page</a> for details on key customization.

See also helpString() and mouseString(). */
QString QGLViewer::keyboardString() const {
  QString text(
      "<center><table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">\n");
  text += QString("<tr bgcolor=\"#aaaacc\"><th align=\"center\">%1</th><th "
                  "align=\"center\">%2</th></tr>\n")
              .arg(QGLViewer::tr("Key(s)",
                                 "Keys column header in help window mouse tab"))
              .arg(QGLViewer::tr(
                  "Description",
                  "Description column header in help window mouse tab"));

  QMap<unsigned int, QString> keyDescription;

  // 1 - User defined key descriptions
  for (QMap<unsigned int, QString>::ConstIterator kd = keyDescription_.begin(),
                                                  kdend = keyDescription_.end();
       kd != kdend; ++kd)
    keyDescription[kd.key()] = kd.value();

  // Add to text in sorted order
  for (QMap<unsigned int, QString>::ConstIterator kb = keyDescription.begin(),
                                                  endb = keyDescription.end();
       kb != endb; ++kb)
    text += tableLine(keyString(kb.key()), kb.value());

  // 2 - Optional separator line
  if (!keyDescription.isEmpty()) {
    keyDescription.clear();
    text += QString("<tr bgcolor=\"#aaaacc\"><td colspan=2>%1</td></tr>\n")
                .arg(QGLViewer::tr("Standard viewer keys",
                                   "In help window keys tab"));
  }

  // 3 - KeyboardAction bindings description
  for (QMap<KeyboardAction, unsigned int>::ConstIterator
           it = keyboardBinding_.begin(),
           end = keyboardBinding_.end();
       it != end; ++it)
    if ((it.value() != 0) &&
        ((!cameraIsInRotateMode()) ||
         ((it.key() != INCREASE_FLYSPEED) && (it.key() != DECREASE_FLYSPEED))))
      keyDescription[it.value()] = keyboardActionDescription_[it.key()];

  // Add to text in sorted order
  for (QMap<unsigned int, QString>::ConstIterator kb2 = keyDescription.begin(),
                                                  endb2 = keyDescription.end();
       kb2 != endb2; ++kb2)
    text += tableLine(keyString(kb2.key()), kb2.value());

  // 4 - Camera paths keys description
  const QString cpks = cameraPathKeysString();
  if (!cpks.isNull()) {
    text += "<tr bgcolor=\"#ccccff\"><td colspan=2>\n";
    text += QGLViewer::tr("Camera paths are controlled using the %1 keys "
                          "(noted <i>Fx</i> below):",
                          "Help window key tab camera keys")
                .arg(cpks) +
            "</td></tr>\n";
    text += tableLine(
        keyString(playPathKeyboardModifiers()) + "<i>" +
            QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>",
        QGLViewer::tr("Plays path (or resets saved position)"));
    text += tableLine(
        keyString(addKeyFrameKeyboardModifiers()) + "<i>" +
            QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>",
        QGLViewer::tr("Adds a key frame to path (or defines a position)"));
    text += tableLine(
        keyString(addKeyFrameKeyboardModifiers()) + "<i>" +
            QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>+<i>" +
            QGLViewer::tr("Fx", "Generic function key (F1..F12)") + "</i>",
        QGLViewer::tr("Deletes path (or saved position)"));
  }
  text += "</table></center>";

  return text;
}

/*! Displays the help window "About" tab. See help() for details. */
void QGLViewer::aboutQGLViewer() {
  help();
  helpWidget()->setCurrentIndex(3);
}

/*! Opens a modal help window that includes four tabs, respectively filled with
helpString(), keyboardString(), mouseString() and about libQGLViewer.

Rich html-like text can be used (see the QStyleSheet documentation). This method
is called when the user presses the QGLViewer::HELP key (default is 'H').

You can use helpWidget() to access to the help widget (to add/remove tabs,
change layout...).

The helpRequired() signal is emitted. */
void QGLViewer::help() {
  Q_EMIT helpRequired();

  bool resize = false;
  int width = 600;
  int height = 400;

  static QString label[] = {tr("&Help", "Help window tab title"),
                            tr("&Keyboard", "Help window tab title"),
                            tr("&Mouse", "Help window tab title"),
                            tr("&About", "Help window about title")};

  if (!helpWidget()) {
    // Qt4 requires a NULL parent...
    helpWidget_ = new QTabWidget(NULL);
    helpWidget()->setWindowTitle(tr("Help", "Help window title"));

    resize = true;
    for (int i = 0; i < 4; ++i) {
      QTextEdit *tab = new QTextEdit(NULL);
      tab->setReadOnly(true);

      helpWidget()->insertTab(i, tab, label[i]);
      if (i == 3) {
#include "qglviewer-icon.xpm"
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
                "<h3>Version %1</h3><br>"
                "A versatile 3D viewer based on OpenGL and Qt<br>"
                "Copyright 2002-%2 Gilles Debunne<br>"
                "<code>%3</code>")
                 .arg(QGLViewerVersionString())
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
Viewer::keyPressEvent(QKeyEvent *e)
{
  // Defines the Alt+R shortcut.
  if ((e->key() == Qt::Key_R) && (e->modifiers() == Qt::AltModifier))
  {
  myResetFunction();
  update(); // Refresh display
  }
  else
  QGLViewer::keyPressEvent(e);
}

// With Qt 2 or 3, you would retrieve modifiers keys using :
// const Qt::ButtonState modifiers = (Qt::ButtonState)(e->state() &
Qt::KeyButtonMask); \endcode When you define a new keyboard shortcut, use
setKeyDescription() to provide a short description which is displayed in the
help() window Keyboard tab. See the <a
href="../examples/keyboardAndMouse.html">keyboardAndMouse</a> example for an
illustration.

See also QOpenGLWidget::keyReleaseEvent(). */
void QGLViewer::keyPressEvent(QKeyEvent *e) {
  if (e->key() == 0) {
    e->ignore();
    return;
  }

  const Qt::Key key = Qt::Key(e->key());

  const Qt::KeyboardModifiers modifiers = e->modifiers();

  QMap<KeyboardAction, unsigned int>::ConstIterator it = keyboardBinding_
                                                             .begin(),
                                                    end =
                                                        keyboardBinding_.end();
  const unsigned int target = key | modifiers;
  while ((it != end) && (it.value() != target))
    ++it;

  if (it != end)
    handleKeyboardAction(it.key());
  else if (pathIndex_.contains(Qt::Key(key))) {
    // Camera paths
    unsigned int index = pathIndex_[Qt::Key(key)];

    // not safe, but try to double press on two viewers at the same time !
    static QTime doublePress;

    if (modifiers == playPathKeyboardModifiers()) {
      int elapsed = doublePress.restart();
      if ((elapsed < 250) && (index == previousPathId_))
        camera()->resetPath(index);
      else {
        // Stop previous interpolation before starting a new one.
        if (index != previousPathId_) {
          KeyFrameInterpolator *previous =
              camera()->keyFrameInterpolator(previousPathId_);
          if ((previous) && (previous->interpolationIsStarted()))
            previous->resetInterpolation();
        }
        camera()->playPath(index);
      }
      previousPathId_ = index;
    } else if (modifiers == addKeyFrameKeyboardModifiers()) {
      int elapsed = doublePress.restart();
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
        bool nullBefore = (camera()->keyFrameInterpolator(index) == NULL);
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

void QGLViewer::keyReleaseEvent(QKeyEvent *e) {
  if (isValidShortcutKey(e->key()))
    currentlyPressedKey_ = Qt::Key(0);
}

void QGLViewer::handleKeyboardAction(KeyboardAction id) {
  switch (id) {
  case DRAW_AXIS:
    toggleAxisIsDrawn();
    break;
  case DRAW_GRID:
    toggleGridIsDrawn();
    break;
  case DISPLAY_FPS:
    toggleFPSIsDisplayed();
    break;
  case ENABLE_TEXT:
    toggleTextIsEnabled();
    break;
  case EXIT_VIEWER:
    saveStateToFileForAllViewers();
    qApp->closeAllWindows();
    break;
  case SAVE_SCREENSHOT:
    saveSnapshot(false, false);
    break;
  case FULL_SCREEN:
    toggleFullScreen();
    break;
  case STEREO:
    toggleStereoDisplay();
    break;
  case ANIMATION:
    toggleAnimation();
    break;
  case HELP:
    help();
    break;
  case EDIT_CAMERA:
    toggleCameraIsEdited();
    break;
  case SNAPSHOT_TO_CLIPBOARD:
    snapshotToClipboard();
    break;
  case CAMERA_MODE:
    toggleCameraMode();
    displayMessage(cameraIsInRotateMode()
                       ? tr("Camera in observer mode", "Feedback message")
                       : tr("Camera in fly mode", "Feedback message"));
    break;

  case MOVE_CAMERA_LEFT:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        Vec(-10.0 * camera()->flySpeed(), 0.0, 0.0)));
    update();
    break;
  case MOVE_CAMERA_RIGHT:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        Vec(10.0 * camera()->flySpeed(), 0.0, 0.0)));
    update();
    break;
  case MOVE_CAMERA_UP:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        Vec(0.0, 10.0 * camera()->flySpeed(), 0.0)));
    update();
    break;
  case MOVE_CAMERA_DOWN:
    camera()->frame()->translate(camera()->frame()->inverseTransformOf(
        Vec(0.0, -10.0 * camera()->flySpeed(), 0.0)));
    update();
    break;

  case INCREASE_FLYSPEED:
    camera()->setFlySpeed(camera()->flySpeed() * 1.5);
    break;
  case DECREASE_FLYSPEED:
    camera()->setFlySpeed(camera()->flySpeed() / 1.5);
    break;
  }
}

/*! Callback method used when the widget size is modified.

If you overload this method, first call the inherited method. Also called when
the widget is created, before its first display. */
void QGLViewer::resizeGL(int width, int height) {
  QOpenGLWidget::resizeGL(width, height);
  glViewport(0, 0, GLint(width), GLint(height));
  camera()->setScreenWidthAndHeight(this->width(), this->height());
}

//////////////////////////////////////////////////////////////////////////
//              K e y b o a r d   s h o r t c u t s                     //
//////////////////////////////////////////////////////////////////////////

/*! Defines the shortcut() that triggers a given QGLViewer::KeyboardAction.

Here are some examples:
\code
// Press 'Q' to exit application
setShortcut(EXIT_VIEWER, Qt::Key_Q);

// Alt+M toggles camera mode
setShortcut(CAMERA_MODE, Qt::ALT + Qt::Key_M);

// The DISPLAY_FPS action is disabled
setShortcut(DISPLAY_FPS, 0);
\endcode

Only one shortcut can be assigned to a given QGLViewer::KeyboardAction (new
bindings replace previous ones). If several KeyboardAction are binded to the
same shortcut, only one of them is active. */
void QGLViewer::setShortcut(KeyboardAction action, unsigned int key) {
  keyboardBinding_[action] = key;
}

/*! Returns the keyboard shortcut associated to a given
QGLViewer::KeyboardAction.

Result is an \c unsigned \c int defined using Qt enumerated values, as in \c
Qt::Key_Q or \c Qt::CTRL + Qt::Key_X. Use Qt::MODIFIER_MASK to separate the key
from the state keys. Returns \c 0 if the KeyboardAction is disabled (not
binded). Set using setShortcut().

If you want to define keyboard shortcuts for custom actions (say, open a scene
file), overload keyPressEvent() and then setKeyDescription().

These shortcuts and their descriptions are automatically included in the help()
window \c Keyboard tab.

See the <a href="../keyboard.html">keyboard page</a> for details and default
values and the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse</a>
example for a practical illustration. */
unsigned int QGLViewer::shortcut(KeyboardAction action) const {
  if (keyboardBinding_.contains(action))
    return keyboardBinding_[action];
  else
    return 0;
}

#ifndef DOXYGEN
void QGLViewer::setKeyboardAccelerator(KeyboardAction action,
                                       unsigned int key) {
  qWarning("setKeyboardAccelerator is deprecated. Use setShortcut instead.");
  setShortcut(action, key);
}

unsigned int QGLViewer::keyboardAccelerator(KeyboardAction action) const {
  qWarning("keyboardAccelerator is deprecated. Use shortcut instead.");
  return shortcut(action);
}
#endif

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
Qt::Key QGLViewer::pathKey(unsigned int index) const {
  for (QMap<Qt::Key, unsigned int>::ConstIterator it = pathIndex_.begin(),
                                                  end = pathIndex_.end();
       it != end; ++it)
    if (it.value() == index)
      return it.key();
  return Qt::Key(0);
}

/*! Sets the pathKey() associated with the camera Key Frame path \p index.

Several keys can be binded to the same \p index. Use a negated \p key value to
delete the binding (the \p index value is then ignored): \code
// Press 'space' to play/pause/add/delete camera path of index 0.
setPathKey(Qt::Key_Space, 0);

// Remove this binding
setPathKey(-Qt::Key_Space);
\endcode */
void QGLViewer::setPathKey(int key, unsigned int index) {
  Qt::Key k = Qt::Key(abs(key));
  if (key < 0)
    pathIndex_.remove(k);
  else
    pathIndex_[k] = index;
}

/*! Sets the playPathKeyboardModifiers() value. */
void QGLViewer::setPlayPathKeyboardModifiers(Qt::KeyboardModifiers modifiers) {
  playPathKeyboardModifiers_ = modifiers;
}

/*! Sets the addKeyFrameKeyboardModifiers() value. */
void QGLViewer::setAddKeyFrameKeyboardModifiers(
    Qt::KeyboardModifiers modifiers) {
  addKeyFrameKeyboardModifiers_ = modifiers;
}

/*! Returns the keyboard modifiers that must be pressed with a pathKey() to add
the current camera position to a KeyFrame path.

It can be \c Qt::NoModifier, \c Qt::ControlModifier, \c Qt::ShiftModifier, \c
Qt::AltModifier, \c Qt::MetaModifier or a combination of these (using the
bitwise '|' operator).

Default value is Qt::AltModifier. Defined using
setAddKeyFrameKeyboardModifiers().

See also playPathKeyboardModifiers(). */
Qt::KeyboardModifiers QGLViewer::addKeyFrameKeyboardModifiers() const {
  return addKeyFrameKeyboardModifiers_;
}

/*! Returns the keyboard modifiers that must be pressed with a pathKey() to play
a camera KeyFrame path.

It can be \c Qt::NoModifier, \c Qt::ControlModifier, \c Qt::ShiftModifier, \c
Qt::AltModifier, \c Qt::MetaModifier or a combination of these (using the
bitwise '|' operator).

Default value is Qt::NoModifier. Defined using setPlayPathKeyboardModifiers().

See also addKeyFrameKeyboardModifiers(). */
Qt::KeyboardModifiers QGLViewer::playPathKeyboardModifiers() const {
  return playPathKeyboardModifiers_;
}

#ifndef DOXYGEN
// Deprecated methods
Qt::KeyboardModifiers QGLViewer::addKeyFrameStateKey() const {
  qWarning("addKeyFrameStateKey has been renamed addKeyFrameKeyboardModifiers");
  return addKeyFrameKeyboardModifiers();
}

Qt::KeyboardModifiers QGLViewer::playPathStateKey() const {
  qWarning("playPathStateKey has been renamed playPathKeyboardModifiers");
  return playPathKeyboardModifiers();
}

void QGLViewer::setAddKeyFrameStateKey(unsigned int buttonState) {
  qWarning("setAddKeyFrameStateKey has been renamed "
           "setAddKeyFrameKeyboardModifiers");
  setAddKeyFrameKeyboardModifiers(keyboardModifiersFromState(buttonState));
}

void QGLViewer::setPlayPathStateKey(unsigned int buttonState) {
  qWarning("setPlayPathStateKey has been renamed setPlayPathKeyboardModifiers");
  setPlayPathKeyboardModifiers(keyboardModifiersFromState(buttonState));
}

Qt::Key QGLViewer::keyFrameKey(unsigned int index) const {
  qWarning("keyFrameKey has been renamed pathKey.");
  return pathKey(index);
}

Qt::KeyboardModifiers QGLViewer::playKeyFramePathStateKey() const {
  qWarning(
      "playKeyFramePathStateKey has been renamed playPathKeyboardModifiers.");
  return playPathKeyboardModifiers();
}

void QGLViewer::setKeyFrameKey(unsigned int index, int key) {
  qWarning("setKeyFrameKey is deprecated, use setPathKey instead, with swapped "
           "parameters.");
  setPathKey(key, index);
}

void QGLViewer::setPlayKeyFramePathStateKey(unsigned int buttonState) {
  qWarning("setPlayKeyFramePathStateKey has been renamed "
           "setPlayPathKeyboardModifiers.");
  setPlayPathKeyboardModifiers(keyboardModifiersFromState(buttonState));
}
#endif

////////////////////////////////////////////////////////////////////////////////
//              M o u s e   b e h a v i o r   s t a t e   k e y s             //
////////////////////////////////////////////////////////////////////////////////
#ifndef DOXYGEN
/*! This method has been deprecated since version 2.5.0

Associates keyboard modifiers to MouseHandler \p handler.

The \p modifiers parameter is \c Qt::AltModifier, \c Qt::ShiftModifier, \c
Qt::ControlModifier, \c Qt::MetaModifier or a combination of these using the '|'
bitwise operator.

\e All the \p handler's associated bindings will then need the specified \p
modifiers key(s) to be activated.

With this code,
\code
setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::AltModifier);
setHandlerKeyboardModifiers(QGLViewer::FRAME, Qt::NoModifier);
\endcode
you will have to press the \c Alt key while pressing mouse buttons in order to
move the camera(), while no key will be needed to move the manipulatedFrame().

This method has a very basic implementation: every action binded to \p handler
has its keyboard modifier replaced by \p modifiers. If \p handler had some
actions binded to different modifiers, these settings will be lost. You should
hence consider using setMouseBinding() for finer tuning.

The default binding associates \c Qt::ControlModifier to all the
QGLViewer::FRAME actions and \c Qt::NoModifier to all QGLViewer::CAMERA actions.
See <a href="../mouse.html">mouse page</a> for details.

\attention This method calls setMouseBinding(), which ensures that only one
action is binded to a given modifiers. If you want to \e swap the
QGLViewer::CAMERA and QGLViewer::FRAME keyboard modifiers, you have to use a
temporary dummy modifier (as if you were swapping two variables) or else the
first call will overwrite the previous settings: \code
// Associate FRAME with Alt (temporary value)
setHandlerKeyboardModifiers(QGLViewer::FRAME, Qt::AltModifier);
// Control is associated with CAMERA
setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::ControlModifier);
// And finally, FRAME can be associated with NoModifier
setHandlerKeyboardModifiers(QGLViewer::FRAME, Qt::NoModifier);
\endcode */
void QGLViewer::setHandlerKeyboardModifiers(MouseHandler handler,
                                            Qt::KeyboardModifiers modifiers) {
  qWarning("setHandlerKeyboardModifiers is deprecated, call setMouseBinding() "
           "instead");

  QMap<MouseBindingPrivate, MouseActionPrivate> newMouseBinding;
  QMap<WheelBindingPrivate, MouseActionPrivate> newWheelBinding;
  QMap<ClickBindingPrivate, ClickAction> newClickBinding_;

  QMap<MouseBindingPrivate, MouseActionPrivate>::Iterator mit;
  QMap<WheelBindingPrivate, MouseActionPrivate>::Iterator wit;

  // First copy unchanged bindings.
  for (mit = mouseBinding_.begin(); mit != mouseBinding_.end(); ++mit)
    if ((mit.value().handler != handler) ||
        (mit.value().action == ZOOM_ON_REGION))
      newMouseBinding[mit.key()] = mit.value();

  for (wit = wheelBinding_.begin(); wit != wheelBinding_.end(); ++wit)
    if (wit.value().handler != handler)
      newWheelBinding[wit.key()] = wit.value();

  // Then, add modified bindings, that can overwrite the previous ones.
  for (mit = mouseBinding_.begin(); mit != mouseBinding_.end(); ++mit)
    if ((mit.value().handler == handler) &&
        (mit.value().action != ZOOM_ON_REGION)) {
      MouseBindingPrivate mbp(modifiers, mit.key().button, mit.key().key);
      newMouseBinding[mbp] = mit.value();
    }

  for (wit = wheelBinding_.begin(); wit != wheelBinding_.end(); ++wit)
    if (wit.value().handler == handler) {
      WheelBindingPrivate wbp(modifiers, wit.key().key);
      newWheelBinding[wbp] = wit.value();
    }

  // Same for button bindings
  for (QMap<ClickBindingPrivate, ClickAction>::ConstIterator
           cb = clickBinding_.begin(),
           end = clickBinding_.end();
       cb != end; ++cb)
    if (((handler == CAMERA) &&
         ((cb.value() == CENTER_SCENE) || (cb.value() == ALIGN_CAMERA))) ||
        ((handler == FRAME) &&
         ((cb.value() == CENTER_FRAME) || (cb.value() == ALIGN_FRAME)))) {
      ClickBindingPrivate cbp(modifiers, cb.key().button, cb.key().doubleClick,
                              cb.key().buttonsBefore, cb.key().key);
      newClickBinding_[cbp] = cb.value();
    } else
      newClickBinding_[cb.key()] = cb.value();

  mouseBinding_ = newMouseBinding;
  wheelBinding_ = newWheelBinding;
  clickBinding_ = newClickBinding_;
}

void QGLViewer::setHandlerStateKey(MouseHandler handler,
                                   unsigned int buttonState) {
  qWarning("setHandlerStateKey has been renamed setHandlerKeyboardModifiers");
  setHandlerKeyboardModifiers(handler, keyboardModifiersFromState(buttonState));
}

void QGLViewer::setMouseStateKey(MouseHandler handler,
                                 unsigned int buttonState) {
  qWarning("setMouseStateKey has been renamed setHandlerKeyboardModifiers.");
  setHandlerKeyboardModifiers(handler, keyboardModifiersFromState(buttonState));
}

/*! This method is deprecated since version 2.5.0

 Use setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButtons, MouseHandler,
 MouseAction, bool) instead.
*/
void QGLViewer::setMouseBinding(unsigned int state, MouseHandler handler,
                                MouseAction action, bool withConstraint) {
  qWarning("setMouseBinding(int state, MouseHandler...) is deprecated. Use the "
           "modifier/button equivalent");
  setMouseBinding(keyboardModifiersFromState(state),
                  mouseButtonFromState(state), handler, action, withConstraint);
}
#endif

/*! Defines a MouseAction binding.

  Same as calling setMouseBinding(Qt::Key, Qt::KeyboardModifiers,
  Qt::MouseButton, MouseHandler, MouseAction, bool), with a key value of
  Qt::Key(0) (i.e. no regular extra key needs to be pressed to perform this
  action). */
void QGLViewer::setMouseBinding(Qt::KeyboardModifiers modifiers,
                                Qt::MouseButton button, MouseHandler handler,
                                MouseAction action, bool withConstraint) {
  setMouseBinding(Qt::Key(0), modifiers, button, handler, action,
                  withConstraint);
}

/*! Associates a MouseAction to any mouse \p button, while keyboard \p modifiers
and \p key are pressed. The receiver of the mouse events is a MouseHandler
(QGLViewer::CAMERA or QGLViewer::FRAME).

The parameters should read: when the mouse \p button is pressed, while the
keyboard \p modifiers and \p key are down, activate \p action on \p handler. Use
Qt::NoModifier to indicate that no modifier key is needed, and a \p key value of
0 if no regular key has to be pressed (or simply use
setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButton, MouseHandler,
MouseAction, bool)).

Use the '|' operator to combine modifiers:
\code
// The R key combined with the Left mouse button rotates the camera in the
screen plane. setMouseBinding(Qt::Key_R, Qt::NoModifier, Qt::LeftButton, CAMERA,
SCREEN_ROTATE);

// Alt + Shift and Left button rotates the manipulatedFrame().
setMouseBinding(Qt::AltModifier | Qt::ShiftModifier, Qt::LeftButton, FRAME,
ROTATE); \endcode

If \p withConstraint is \c true (default), the possible
qglviewer::Frame::constraint() of the associated Frame will be enforced during
motion.

The list of all possible MouseAction, some binding examples and default bindings
are provided in the <a href="../mouse.html">mouse page</a>.

See the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse</a> example
for an illustration.

If no mouse button is specified, the binding is ignored. If an action was
previously associated with this keyboard and button combination, it is silently
overwritten (call mouseAction() before to check).

To remove a specific mouse binding, use \p NO_MOUSE_ACTION as the \p action.

See also setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButtons, ClickAction,
bool, int), setWheelBinding() and clearMouseBindings(). */
void QGLViewer::setMouseBinding(Qt::Key key, Qt::KeyboardModifiers modifiers,
                                Qt::MouseButton button, MouseHandler handler,
                                MouseAction action, bool withConstraint) {
  if ((handler == FRAME) &&
      ((action == MOVE_FORWARD) || (action == MOVE_BACKWARD) ||
       (action == ROLL) || (action == LOOK_AROUND) ||
       (action == ZOOM_ON_REGION))) {
    qWarning("Cannot bind %s to FRAME",
             mouseActionString(action).toLatin1().constData());
    return;
  }

  if (button == Qt::NoButton) {
    qWarning("No mouse button specified in setMouseBinding");
    return;
  }

  MouseActionPrivate map;
  map.handler = handler;
  map.action = action;
  map.withConstraint = withConstraint;

  MouseBindingPrivate mbp(modifiers, button, key);
  if (action == NO_MOUSE_ACTION)
    mouseBinding_.remove(mbp);
  else
    mouseBinding_.insert(mbp, map);

  ClickBindingPrivate cbp(modifiers, button, false, Qt::NoButton, key);
  clickBinding_.remove(cbp);
}

#ifndef DOXYGEN
/*! This method is deprecated since version 2.5.0

 Use setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButtons, MouseHandler,
 MouseAction, bool) instead.
*/
void QGLViewer::setMouseBinding(unsigned int state, ClickAction action,
                                bool doubleClick,
                                Qt::MouseButtons buttonsBefore) {
  qWarning("setMouseBinding(int state, ClickAction...) is deprecated. Use the "
           "modifier/button equivalent");
  setMouseBinding(keyboardModifiersFromState(state),
                  mouseButtonFromState(state), action, doubleClick,
                  buttonsBefore);
}
#endif

/*! Defines a ClickAction binding.

 Same as calling setMouseBinding(Qt::Key, Qt::KeyboardModifiers,
 Qt::MouseButton, ClickAction, bool, Qt::MouseButtons), with a key value of
 Qt::Key(0) (i.e. no regular key needs to be pressed to activate this action).
 */
void QGLViewer::setMouseBinding(Qt::KeyboardModifiers modifiers,
                                Qt::MouseButton button, ClickAction action,
                                bool doubleClick,
                                Qt::MouseButtons buttonsBefore) {
  setMouseBinding(Qt::Key(0), modifiers, button, action, doubleClick,
                  buttonsBefore);
}

/*! Associates a ClickAction to a button and keyboard key and modifier(s)
combination.

The parameters should read: when \p button is pressed, while the \p modifiers
and \p key keys are down, and possibly as a \p doubleClick, then perform \p
action. Use Qt::NoModifier to indicate that no modifier key is needed, and a \p
key value of 0 if no regular key has to be pressed (or simply use
setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButton, ClickAction, bool,
Qt::MouseButtons)).

If \p buttonsBefore is specified (valid only when \p doubleClick is \c true),
then this (or these) other mouse button(s) has (have) to be pressed \e before
the double click occurs in order to execute \p action.

The list of all possible ClickAction, some binding examples and default bindings
are listed in the <a href="../mouse.html">mouse page</a>. See also the
setMouseBinding() documentation.

See the <a href="../examples/keyboardAndMouse.html">keyboardAndMouse example</a>
for an illustration.

The binding is ignored if Qt::NoButton is specified as \p buttons.

See also setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButtons, MouseHandler,
MouseAction, bool), setWheelBinding() and clearMouseBindings().
*/
void QGLViewer::setMouseBinding(Qt::Key key, Qt::KeyboardModifiers modifiers,
                                Qt::MouseButton button, ClickAction action,
                                bool doubleClick,
                                Qt::MouseButtons buttonsBefore) {
  if ((buttonsBefore != Qt::NoButton) && !doubleClick) {
    qWarning("Buttons before is only meaningful when doubleClick is true in "
             "setMouseBinding().");
    return;
  }

  if (button == Qt::NoButton) {
    qWarning("No mouse button specified in setMouseBinding");
    return;
  }

  ClickBindingPrivate cbp(modifiers, button, doubleClick, buttonsBefore, key);

  // #CONNECTION performClickAction comment on NO_CLICK_ACTION
  if (action == NO_CLICK_ACTION)
    clickBinding_.remove(cbp);
  else
    clickBinding_.insert(cbp, action);

  if ((!doubleClick) && (buttonsBefore == Qt::NoButton)) {
    MouseBindingPrivate mbp(modifiers, button, key);
    mouseBinding_.remove(mbp);
  }
}

/*! Defines a mouse wheel binding.

 Same as calling setWheelBinding(Qt::Key, Qt::KeyboardModifiers, MouseHandler,
 MouseAction, bool), with a key value of Qt::Key(0) (i.e. no regular key needs
 to be pressed to activate this action). */
void QGLViewer::setWheelBinding(Qt::KeyboardModifiers modifiers,
                                MouseHandler handler, MouseAction action,
                                bool withConstraint) {
  setWheelBinding(Qt::Key(0), modifiers, handler, action, withConstraint);
}

/*! Associates a MouseAction and a MouseHandler to a mouse wheel event.

This method is very similar to setMouseBinding(), but specific to the wheel.

In the current implementation only QGLViewer::ZOOM can be associated with
QGLViewer::FRAME, while QGLViewer::CAMERA can receive QGLViewer::ZOOM and
QGLViewer::MOVE_FORWARD.

The difference between QGLViewer::ZOOM and QGLViewer::MOVE_FORWARD is that
QGLViewer::ZOOM speed depends on the distance to the object, while
QGLViewer::MOVE_FORWARD moves at a constant speed defined by
qglviewer::Camera::flySpeed(). */
void QGLViewer::setWheelBinding(Qt::Key key, Qt::KeyboardModifiers modifiers,
                                MouseHandler handler, MouseAction action,
                                bool withConstraint) {
  //#CONNECTION# ManipulatedFrame::wheelEvent and
  // ManipulatedCameraFrame::wheelEvent switches
  if ((action != ZOOM) && (action != MOVE_FORWARD) &&
      (action != MOVE_BACKWARD) && (action != NO_MOUSE_ACTION)) {
    qWarning("Cannot bind %s to wheel",
             mouseActionString(action).toLatin1().constData());
    return;
  }

  if ((handler == FRAME) && (action != ZOOM) && (action != NO_MOUSE_ACTION)) {
    qWarning("Cannot bind %s to FRAME wheel",
             mouseActionString(action).toLatin1().constData());
    return;
  }

  MouseActionPrivate map;
  map.handler = handler;
  map.action = action;
  map.withConstraint = withConstraint;

  WheelBindingPrivate wbp(modifiers, key);
  if (action == NO_MOUSE_ACTION)
    wheelBinding_.remove(wbp);
  else
    wheelBinding_[wbp] = map;
}

/*! Clears all the default mouse bindings.

After this call, you will have to use setMouseBinding() and setWheelBinding() to
restore the mouse bindings you are interested in.
*/
void QGLViewer::clearMouseBindings() {
  mouseBinding_.clear();
  clickBinding_.clear();
  wheelBinding_.clear();
}

/*! Clears all the default keyboard shortcuts.

After this call, you will have to use setShortcut() to define your own keyboard
shortcuts.
*/
void QGLViewer::clearShortcuts() {
  keyboardBinding_.clear();
  pathIndex_.clear();
}

/*! This method is deprecated since version 2.5.0

 Use mouseAction(Qt::Key, Qt::KeyboardModifiers, Qt::MouseButtons) instead.
*/
QGLViewer::MouseAction QGLViewer::mouseAction(unsigned int state) const {
  qWarning("mouseAction(int state,...) is deprecated. Use the modifier/button "
           "equivalent");
  return mouseAction(Qt::Key(0), keyboardModifiersFromState(state),
                     mouseButtonFromState(state));
}

/*! Returns the MouseAction the will be triggered when the mouse \p button is
pressed, while the keyboard \p modifiers and \p key are pressed.

Returns QGLViewer::NO_MOUSE_ACTION if no action is associated with this
combination. Use 0 for \p key to indicate that no regular key needs to be
pressed.

For instance, to know which motion corresponds to Alt+LeftButton, do:
\code
QGLViewer::MouseAction ma = mouseAction(0, Qt::AltModifier, Qt::LeftButton);
if (ma != QGLViewer::NO_MOUSE_ACTION) ...
\endcode

Use mouseHandler() to know which object (QGLViewer::CAMERA or QGLViewer::FRAME)
will execute this action. */
QGLViewer::MouseAction QGLViewer::mouseAction(Qt::Key key,
                                              Qt::KeyboardModifiers modifiers,
                                              Qt::MouseButton button) const {
  MouseBindingPrivate mbp(modifiers, button, key);
  if (mouseBinding_.contains(mbp))
    return mouseBinding_[mbp].action;
  else
    return NO_MOUSE_ACTION;
}

/*! This method is deprecated since version 2.5.0

 Use mouseHanler(Qt::Key, Qt::KeyboardModifiers, Qt::MouseButtons) instead.
*/
int QGLViewer::mouseHandler(unsigned int state) const {
  qWarning("mouseHandler(int state,...) is deprecated. Use the modifier/button "
           "equivalent");
  return mouseHandler(Qt::Key(0), keyboardModifiersFromState(state),
                      mouseButtonFromState(state));
}

/*! Returns the MouseHandler which will be activated when the mouse \p button is
pressed, while the \p modifiers and \p key are pressed.

If no action is associated with this combination, returns \c -1. Use 0 for \p
key and Qt::NoModifier for \p modifiers to represent the lack of a key press.

For instance, to know which handler receives the Alt+LeftButton, do:
\code
int mh = mouseHandler(0, Qt::AltModifier, Qt::LeftButton);
if (mh == QGLViewer::CAMERA) ...
\endcode

Use mouseAction() to know which action (see the MouseAction enum) will be
performed on this handler. */
int QGLViewer::mouseHandler(Qt::Key key, Qt::KeyboardModifiers modifiers,
                            Qt::MouseButton button) const {
  MouseBindingPrivate mbp(modifiers, button, key);
  if (mouseBinding_.contains(mbp))
    return mouseBinding_[mbp].handler;
  else
    return -1;
}

#ifndef DOXYGEN
/*! This method is deprecated since version 2.5.0

 Use mouseButtons() and keyboardModifiers() instead.
*/
int QGLViewer::mouseButtonState(MouseHandler handler, MouseAction action,
                                bool withConstraint) const {
  qWarning("mouseButtonState() is deprecated. Use mouseButtons() and "
           "keyboardModifiers() instead");
  for (QMap<MouseBindingPrivate, MouseActionPrivate>::ConstIterator
           it = mouseBinding_.begin(),
           end = mouseBinding_.end();
       it != end; ++it)
    if ((it.value().handler == handler) && (it.value().action == action) &&
        (it.value().withConstraint == withConstraint))
      return (int)it.key().modifiers | (int)it.key().button;

  return Qt::NoButton;
}
#endif

/*! Returns the keyboard state that triggers \p action on \p handler \p
withConstraint using the mouse wheel.

If such a binding exists, results are stored in the \p key and \p modifiers
parameters. If the MouseAction \p action is not bound, \p key is set to the
illegal -1 value. If several keyboard states trigger the MouseAction, one of
them is returned.

See also setMouseBinding(), getClickActionBinding() and getMouseActionBinding().
*/
void QGLViewer::getWheelActionBinding(MouseHandler handler, MouseAction action,
                                      bool withConstraint, Qt::Key &key,
                                      Qt::KeyboardModifiers &modifiers) const {
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

  key = Qt::Key(-1);
  modifiers = Qt::NoModifier;
}

/*! Returns the mouse and keyboard state that triggers \p action on \p handler
\p withConstraint.

If such a binding exists, results are stored in the \p key, \p modifiers and \p
button parameters. If the MouseAction \p action is not bound, \p button is set
to \c Qt::NoButton. If several mouse and keyboard states trigger the
MouseAction, one of them is returned.

See also setMouseBinding(), getClickActionBinding() and getWheelActionBinding().
*/
void QGLViewer::getMouseActionBinding(MouseHandler handler, MouseAction action,
                                      bool withConstraint, Qt::Key &key,
                                      Qt::KeyboardModifiers &modifiers,
                                      Qt::MouseButton &button) const {
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

  key = Qt::Key(0);
  modifiers = Qt::NoModifier;
  button = Qt::NoButton;
}

/*! Returns the MouseAction (if any) that is performed when using the wheel,
when the \p modifiers and \p key keyboard keys are pressed.

Returns NO_MOUSE_ACTION if no such binding has been defined using
setWheelBinding().

Same as mouseAction(), but for the wheel action. See also wheelHandler().
*/
QGLViewer::MouseAction
QGLViewer::wheelAction(Qt::Key key, Qt::KeyboardModifiers modifiers) const {
  WheelBindingPrivate wbp(modifiers, key);
  if (wheelBinding_.contains(wbp))
    return wheelBinding_[wbp].action;
  else
    return NO_MOUSE_ACTION;
}

/*! Returns the MouseHandler (if any) that receives wheel events when the \p
  modifiers and \p key keyboard keys are pressed.

  Returns -1 if no no such binding has been defined using setWheelBinding(). See
  also wheelAction().
*/
int QGLViewer::wheelHandler(Qt::Key key,
                            Qt::KeyboardModifiers modifiers) const {
  WheelBindingPrivate wbp(modifiers, key);
  if (wheelBinding_.contains(wbp))
    return wheelBinding_[wbp].handler;
  else
    return -1;
}

/*! Same as mouseAction(), but for the ClickAction set using setMouseBinding().

Returns NO_CLICK_ACTION if no click action is associated with this keyboard and
mouse buttons combination. */
QGLViewer::ClickAction
QGLViewer::clickAction(Qt::Key key, Qt::KeyboardModifiers modifiers,
                       Qt::MouseButton button, bool doubleClick,
                       Qt::MouseButtons buttonsBefore) const {
  ClickBindingPrivate cbp(modifiers, button, doubleClick, buttonsBefore, key);
  if (clickBinding_.contains(cbp))
    return clickBinding_[cbp];
  else
    return NO_CLICK_ACTION;
}

#ifndef DOXYGEN
/*! This method is deprecated since version 2.5.0

  Use wheelAction(Qt::Key key, Qt::KeyboardModifiers modifiers) instead. */
QGLViewer::MouseAction
QGLViewer::wheelAction(Qt::KeyboardModifiers modifiers) const {
  qWarning("wheelAction() is deprecated. Use the new wheelAction() method with "
           "a key parameter instead");
  return wheelAction(Qt::Key(0), modifiers);
}

/*! This method is deprecated since version 2.5.0

  Use wheelHandler(Qt::Key key, Qt::KeyboardModifiers modifiers) instead. */
int QGLViewer::wheelHandler(Qt::KeyboardModifiers modifiers) const {
  qWarning("wheelHandler() is deprecated. Use the new wheelHandler() method "
           "with a key parameter instead");
  return wheelHandler(Qt::Key(0), modifiers);
}

/*! This method is deprecated since version 2.5.0

  Use wheelAction() and wheelHandler() instead. */
unsigned int QGLViewer::wheelButtonState(MouseHandler handler,
                                         MouseAction action,
                                         bool withConstraint) const {
  qWarning("wheelButtonState() is deprecated. Use the wheelAction() and "
           "wheelHandler() instead");
  for (QMap<WheelBindingPrivate, MouseActionPrivate>::ConstIterator
           it = wheelBinding_.begin(),
           end = wheelBinding_.end();
       it != end; ++it)
    if ((it.value().handler == handler) && (it.value().action == action) &&
        (it.value().withConstraint == withConstraint))
      return it.key().key + it.key().modifiers;

  return -1;
}

/*! This method is deprecated since version 2.5.0

 Use clickAction(Qt::KeyboardModifiers, Qt::MouseButtons, bool,
 Qt::MouseButtons) instead.
*/
QGLViewer::ClickAction
QGLViewer::clickAction(unsigned int state, bool doubleClick,
                       Qt::MouseButtons buttonsBefore) const {
  qWarning("clickAction(int state,...) is deprecated. Use the modifier/button "
           "equivalent");
  return clickAction(Qt::Key(0), keyboardModifiersFromState(state),
                     mouseButtonFromState(state), doubleClick, buttonsBefore);
}

/*! This method is deprecated since version 2.5.0

 Use getClickActionState(ClickAction, Qt::Key, Qt::KeyboardModifiers,
 Qt::MouseButton, bool, Qt::MouseButtons) instead.
*/
void QGLViewer::getClickButtonState(ClickAction action, unsigned int &state,
                                    bool &doubleClick,
                                    Qt::MouseButtons &buttonsBefore) const {
  qWarning("getClickButtonState(int state,...) is deprecated. Use the "
           "modifier/button equivalent");
  Qt::KeyboardModifiers modifiers;
  Qt::MouseButton button;
  Qt::Key key;
  getClickActionBinding(action, key, modifiers, button, doubleClick,
                        buttonsBefore);
  state = (unsigned int)modifiers | (unsigned int)button | (unsigned int)key;
}
#endif

/*! Returns the mouse and keyboard state that triggers \p action.

If such a binding exists, results are stored in the \p key, \p modifiers, \p
button, \p doubleClick and \p buttonsBefore parameters. If the ClickAction \p
action is not bound, \p button is set to \c Qt::NoButton. If several mouse
buttons trigger in the ClickAction, one of them is returned.

See also setMouseBinding(), getMouseActionBinding() and getWheelActionBinding().
*/
void QGLViewer::getClickActionBinding(ClickAction action, Qt::Key &key,
                                      Qt::KeyboardModifiers &modifiers,
                                      Qt::MouseButton &button,
                                      bool &doubleClick,
                                      Qt::MouseButtons &buttonsBefore) const {
  for (QMap<ClickBindingPrivate, ClickAction>::ConstIterator
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

  modifiers = Qt::NoModifier;
  button = Qt::NoButton;
  doubleClick = false;
  buttonsBefore = Qt::NoButton;
  key = Qt::Key(0);
}

/*! This function should be used in conjunction with toggleCameraMode(). It
returns \c true when at least one mouse button is binded to the \c ROTATE
mouseAction. This is crude way of determining which "mode" the camera is in. */
bool QGLViewer::cameraIsInRotateMode() const {
  //#CONNECTION# used in toggleCameraMode() and keyboardString()
  Qt::Key key;
  Qt::KeyboardModifiers modifiers;
  Qt::MouseButton button;
  getMouseActionBinding(CAMERA, ROTATE, true /*constraint*/, key, modifiers,
                        button);
  return button != Qt::NoButton;
}

/*! Swaps between two predefined camera mouse bindings.

The first mode makes the camera observe the scene while revolving around the
qglviewer::Camera::pivotPoint(). The second mode is designed for walkthrough
applications and simulates a flying camera.

Practically, the three mouse buttons are respectively binded to:
\arg In rotate mode: QGLViewer::ROTATE, QGLViewer::ZOOM, QGLViewer::TRANSLATE.
\arg In fly mode: QGLViewer::MOVE_FORWARD, QGLViewer::LOOK_AROUND,
QGLViewer::MOVE_BACKWARD.

The current mode is determined by checking if a mouse button is binded to
QGLViewer::ROTATE for the QGLViewer::CAMERA. The state key that was previously
used to move the camera is preserved. */
void QGLViewer::toggleCameraMode() {
  Qt::Key key;
  Qt::KeyboardModifiers modifiers;
  Qt::MouseButton button;
  getMouseActionBinding(CAMERA, ROTATE, true /*constraint*/, key, modifiers,
                        button);
  bool rotateMode = button != Qt::NoButton;

  if (!rotateMode) {
    getMouseActionBinding(CAMERA, MOVE_FORWARD, true /*constraint*/, key,
                          modifiers, button);
  }

  //#CONNECTION# setDefaultMouseBindings()
  if (rotateMode) {
    camera()->frame()->updateSceneUpVector();
    camera()->frame()->stopSpinning();

    setMouseBinding(modifiers, Qt::LeftButton, CAMERA, MOVE_FORWARD);
    setMouseBinding(modifiers, Qt::MidButton, CAMERA, LOOK_AROUND);
    setMouseBinding(modifiers, Qt::RightButton, CAMERA, MOVE_BACKWARD);

    setMouseBinding(Qt::Key_R, modifiers, Qt::LeftButton, CAMERA, ROLL);

    setMouseBinding(Qt::NoModifier, Qt::LeftButton, NO_CLICK_ACTION, true);
    setMouseBinding(Qt::NoModifier, Qt::MidButton, NO_CLICK_ACTION, true);
    setMouseBinding(Qt::NoModifier, Qt::RightButton, NO_CLICK_ACTION, true);

    setWheelBinding(modifiers, CAMERA, MOVE_FORWARD);
  } else {
    // Should stop flyTimer. But unlikely and not easy.
    setMouseBinding(modifiers, Qt::LeftButton, CAMERA, ROTATE);
    setMouseBinding(modifiers, Qt::MidButton, CAMERA, ZOOM);
    setMouseBinding(modifiers, Qt::RightButton, CAMERA, TRANSLATE);

    setMouseBinding(Qt::Key_R, modifiers, Qt::LeftButton, CAMERA,
                    SCREEN_ROTATE);

    setMouseBinding(Qt::NoModifier, Qt::LeftButton, ALIGN_CAMERA, true);
    setMouseBinding(Qt::NoModifier, Qt::MidButton, SHOW_ENTIRE_SCENE, true);
    setMouseBinding(Qt::NoModifier, Qt::RightButton, CENTER_SCENE, true);

    setWheelBinding(modifiers, CAMERA, ZOOM);
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

Note that a qglviewer::ManipulatedCameraFrame can be set as the
manipulatedFrame(): it is possible to manipulate the camera of a first viewer in
a second viewer. */
void QGLViewer::setManipulatedFrame(ManipulatedFrame *frame) {
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
       (dynamic_cast<ManipulatedCameraFrame *>(manipulatedFrame()) != NULL));

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

Displays the new qglviewer::Camera::pivotPoint() when it is changed. See the <a
href="../mouse.html">mouse page</a> for details. Also draws a line between
qglviewer::Camera::pivotPoint() and mouse cursor when the camera is rotated
around the camera Z axis.

See also setVisualHintsMask() and resetVisualHints(). The hint color is
foregroundColor().

\note These methods may become more interesting one day. The current design is
too limited and should be improved when other visual hints must be drawn.

Limitation : One needs to have access to visualHint_ to overload this method.

Removed from the documentation for this reason. */
void QGLViewer::drawVisualHints() {
  rendering_program.bind();
  vaos[GRID].bind();
  QMatrix4x4 mvpMatrix;
  double mat[16];
  camera()->getModelViewProjectionMatrix(mat);
  for(int i=0; i < 16; i++)
  {
      mvpMatrix.data()[i] = (float)mat[i];
  }
  QMatrix4x4 mvMatrix;
  for(int i=0; i < 16; i++)
  {
    mvMatrix.data()[i] = camera()->orientation().inverse().matrix()[i];
  }
  rendering_program.setUniformValue("mvp_matrix", mvpMatrix);
  rendering_program.setUniformValue("color", QColor(Qt::lightGray));
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(grid_size));
  vaos[GRID].release();
  rendering_program.release();
  
  rendering_program_light.bind();
  vaos[GRID_AXIS].bind();
  
  rendering_program_light.setUniformValue("mvp_matrix", mvpMatrix);
  rendering_program_light.setUniformValue("mv_matrix", mvMatrix);
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(g_axis_size/9));
  vaos[GRID_AXIS].release();
  
  qglviewer::Camera::Type camera_type = camera()->type();
  camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
  for(int i=0; i < 16; i++)
  {
    mvMatrix.data()[i] = camera()->orientation().inverse().matrix()[i];
  }
  mvpMatrix.setToIdentity();
  mvpMatrix.ortho(-1,1,-1,1,-1,1);
  mvpMatrix = mvpMatrix*mvMatrix;
  rendering_program_light.setUniformValue("mvp_matrix", mvpMatrix);
  rendering_program_light.setUniformValue("mv_matrix", mvMatrix);
  camera()->setType(camera_type);
  vaos[AXIS].bind();
  int viewport[4];
  int scissor[4];
  
  // The viewport and the scissor are changed to fit the upper right
  // corner. Original values are saved.
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetIntegerv(GL_SCISSOR_BOX, scissor);
  
  // Axis viewport size, in pixels
  const int size = 100;
  glViewport(width()*devicePixelRatio()-size, height()*devicePixelRatio()-size, size, size);
  glScissor (width()*devicePixelRatio()-size, height()*devicePixelRatio()-size, size, size);
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(axis_size / 9));
  // The viewport and the scissor are restored.
  glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
  glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
  vaos[AXIS].release();
  rendering_program_light.release();
}

/*! Defines the mask that will be used to drawVisualHints(). The only available
mask is currently 1, corresponding to the display of the
qglviewer::Camera::pivotPoint(). resetVisualHints() is automatically called
after \p delay milliseconds (default is 2 seconds). */
void QGLViewer::setVisualHintsMask(int mask, int delay) {
  visualHint_ = visualHint_ | mask;
  QTimer::singleShot(delay, this, SLOT(resetVisualHints()));
}

/*! Reset the mask used by drawVisualHints(). Called by setVisualHintsMask()
 * after 2 seconds to reset the display. */
void QGLViewer::resetVisualHints() { visualHint_ = 0; }
#endif

////////////////////////////////////////////////////////////////////////////////
//       A x i s   a n d   G r i d   d i s p l a y   l i s t s                //
////////////////////////////////////////////////////////////////////////////////

/*! Draws a 3D arrow between the 3D point \p from and the 3D point \p to. 
\p data is filled with the three components of a point, then its normal, and then its color, which makes it filled like this:
[P1.x-P1.y-P1.z-N1.x-N1.y-N1.z-C1.r-C1.g-C1.b|P2.x-P2.y-P2.z-N2.x-N2.y-N2.z-C2.r-C2.g-C2.b|...]
*/
void QGLViewer::drawArrow(double r,double R, int prec, qglviewer::Vec from,
                          qglviewer::Vec to, qglviewer::Vec color, 
                          std::vector<float> &data) {
  qglviewer::Vec temp = to-from;
  QVector3D dir = QVector3D(temp.x, temp.y, temp.z);
  QMatrix4x4 mat;
  mat.setToIdentity();
  mat.translate(from.x, from.y, from.z);
  mat.scale(dir.length());
  dir.normalize();
  float angle = 0.0;
  if(std::sqrt((dir.x()*dir.x()+dir.y()*dir.y())) > 1)
      angle = 90.0f;
  else
      angle =acos(dir.y()/std::sqrt(dir.x()*dir.x()+dir.y()*dir.y()+dir.z()*dir.z()))*180.0/M_PI;

  QVector3D axis;
  axis = QVector3D(dir.z(), 0, -dir.x());
  mat.rotate(angle, axis);

  //Head
  const float Rf = static_cast<float>(R);
  for(int d = 0; d<360; d+= 360/prec)
  {
      float D = (float) (d * M_PI / 180.);
      float a = (float) std::atan(Rf / 0.33);
      QVector4D p(0., 1., 0, 1.);
      QVector4D n(Rf*sin(D), sin(a), Rf*cos(D), 1.);
      QVector4D pR = mat*p;
      QVector4D nR = mat*n;

      //point A1
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back((float)color.x);
      data.push_back((float)color.y);
      data.push_back((float)color.z);

      //point B1
      p = QVector4D(Rf*sin(D), 0.66f, Rf* cos(D), 1.f);
      n = QVector4D(sin(D), sin(a), cos(D), 1.);
      pR = mat*p;
      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back((float)color.x);
      data.push_back((float)color.y);
      data.push_back((float)color.z);
      //point C1
      D = (d+360/prec)*M_PI/180.0;
      p = QVector4D(Rf* sin(D), 0.66f, Rf* cos(D), 1.f);
      n = QVector4D(sin(D), sin(a), cos(D), 1.0);
      pR = mat*p;
      nR = mat*n;

      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back((float)color.x);
      data.push_back((float)color.y);
      data.push_back((float)color.z);

  }

  //cylinder
  //body of the cylinder
  const float rf = static_cast<float>(r);
  for(int d = 0; d<360; d+= 360/prec)
  {
      //point A1
      double D = d*M_PI/180.0;
      QVector4D p(rf*sin(D), 0.66f, rf*cos(D), 1.f);
      QVector4D n(sin(D), 0.f, cos(D), 1.f);
      QVector4D pR = mat*p;
      QVector4D nR = mat*n;

      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back(color.x);
      data.push_back(color.y);
      data.push_back(color.z);
      //point B1
      p = QVector4D(rf * sin(D),0,rf*cos(D), 1.0);
      n = QVector4D(sin(D), 0, cos(D), 1.0);
      pR = mat*p;
      nR = mat*n;


      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back(color.x);
      data.push_back(color.y);
      data.push_back(color.z);
        //point C1
      D = (d+360/prec)*M_PI/180.0;
      p = QVector4D(rf * sin(D),0,rf*cos(D), 1.0);
      n = QVector4D(sin(D), 0, cos(D), 1.0);
      pR = mat*p;
      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back(color.x);
      data.push_back(color.y);
      data.push_back(color.z);
      //point A2
      D = (d+360/prec)*M_PI/180.0;

      p = QVector4D(rf * sin(D),0,rf*cos(D), 1.0);
      n = QVector4D(sin(D), 0, cos(D), 1.0);
      pR = mat*p;
      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back((float)color.x);
      data.push_back((float)color.y);
      data.push_back((float)color.z);
      //point B2
      p = QVector4D(rf * sin(D), 0.66f, rf*cos(D), 1.f);
      n = QVector4D(sin(D), 0, cos(D), 1.0);
      pR = mat*p;
      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back((float)color.x);
      data.push_back((float)color.y);
      data.push_back((float)color.z);
      //point C2
      D = d*M_PI/180.0;
      p = QVector4D(rf * sin(D), 0.66f, rf*cos(D), 1.f);
      n = QVector4D(sin(D), 0.f, cos(D), 1.f);
      pR = mat*p;
      nR = mat*n;
      data.push_back(pR.x());
      data.push_back(pR.y());
      data.push_back(pR.z());
      data.push_back(nR.x());
      data.push_back(nR.y());
      data.push_back(nR.z());
      data.push_back(color.x);
      data.push_back(color.y);
      data.push_back(color.z);

  }
}

/*! Draws an XYZ axis, with a given size (default is 1.0).

The axis orientation matches the current modelView matrix state:
three arrows (red, green and blue) of length \p length are drawn along the
positive X, Y and Z directions in the top right corner of the screen. 
X arrow is red, Y arrow is green and Z arrow is blue.*/
void QGLViewer::drawAxis(qreal length) {
  std::vector<float> data;
  data.resize(0);
  drawArrow(0.06,0.12,10, qglviewer::Vec(0,0,0),qglviewer::Vec(length,0,0),qglviewer::Vec(1,0,0), data);
  drawArrow(0.06,0.12,10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,length,0),qglviewer::Vec(0,1,0), data);
  drawArrow(0.06,0.12,10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,0,length),qglviewer::Vec(0,0,1), data);
  rendering_program_light.bind();
  vaos[AXIS].bind();
  vbos[Axis].bind();
  vbos[Axis].allocate(data.data(), static_cast<int>(data.size()) * sizeof(float));
  rendering_program_light.enableAttributeArray("vertex");
  rendering_program_light.setAttributeBuffer("vertex",GL_FLOAT,0,3,
                                             static_cast<int>(9*sizeof(float)));
  
  rendering_program_light.enableAttributeArray("normal");
  rendering_program_light.setAttributeBuffer("normal",GL_FLOAT,3*sizeof(float),3,
                                             static_cast<int>(9*sizeof(float)));
  
  rendering_program_light.enableAttributeArray("colors");
  rendering_program_light.setAttributeBuffer("colors",GL_FLOAT,6*sizeof(float),3,
                                             static_cast<int>(9*sizeof(float)));
  vbos[Axis].release();
  vaos[AXIS].release();
  axis_size = data.size();
  rendering_program_light.release();
}

/*! Draws a grid in the XY plane, centered on (0,0,0) (defined in the current
coordinate system).

\p size (OpenGL units) and \p nbSubdivisions define its geometry.*/
void QGLViewer::drawGrid(qreal size, int nbSubdivisions) {
  
  //The Grid
  std::vector<float> v_Grid;
  for (int i=0; i<=nbSubdivisions; ++i)
  {
          const float pos = size*(2.0*i/nbSubdivisions-1.0);
          v_Grid.push_back(pos);
          v_Grid.push_back(-size);
          v_Grid.push_back(0.0);

          v_Grid.push_back(pos);
          v_Grid.push_back(+size);
          v_Grid.push_back(0.0);

          v_Grid.push_back(-size);
          v_Grid.push_back(pos);
          v_Grid.push_back(0.0);

          v_Grid.push_back( size);
          v_Grid.push_back( pos);
          v_Grid.push_back( 0.0);
  }
  rendering_program.bind();
  vaos[GRID].bind();
  vbos[Grid].bind();
  vbos[Grid].allocate(v_Grid.data(),static_cast<int>(v_Grid.size()*sizeof(float)));
  rendering_program.enableAttributeArray("vertex");
  rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  vbos[Grid].release();
  vaos[GRID].release();
  rendering_program.release();
  grid_size = v_Grid.size();
  
  //The Axis
  std::vector<float> d_axis;
  d_axis.resize(0);
  //d_axis is filled by drawArrow always this way : V.x V.y V.z N.x N.y N.z C.r C.g C.b, so it is possible
  //to use a single buffer with offset and stride
  drawArrow(0.005*size,0.02*size,10, qglviewer::Vec(0,0,0),qglviewer::Vec(size,0,0),qglviewer::Vec(1,0,0), d_axis);
  drawArrow(0.005*size,0.02*size,10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,size,0),qglviewer::Vec(0,1,0), d_axis);

  rendering_program_light.bind();
  vaos[GRID_AXIS].bind();
  vbos[Grid_axis].bind();
  vbos[Grid_axis].allocate(d_axis.data(), static_cast<int>(d_axis.size()) * sizeof(float));
  rendering_program_light.enableAttributeArray("vertex");
  rendering_program_light.setAttributeBuffer("vertex",GL_FLOAT,0,3,
                                             static_cast<int>(9*sizeof(float)));
  
  rendering_program_light.enableAttributeArray("normal");
  rendering_program_light.setAttributeBuffer("normal",GL_FLOAT,3*sizeof(float),3,
                                             static_cast<int>(9*sizeof(float)));
  
  rendering_program_light.enableAttributeArray("colors");
  rendering_program_light.setAttributeBuffer("colors",GL_FLOAT,6*sizeof(float),3,
                                             static_cast<int>(9*sizeof(float)));
  vbos[Grid_axis].release();

  vaos[GRID_AXIS].release();
  rendering_program_light.release();

  g_axis_size = d_axis.size();
}

////////////////////////////////////////////////////////////////////////////////
//       S t a t i c    m e t h o d s   :  Q G L V i e w e r   P o o l        //
////////////////////////////////////////////////////////////////////////////////

/*! saveStateToFile() is called on all the QGLViewers using the QGLViewerPool().
 */
void QGLViewer::saveStateToFileForAllViewers() {
  Q_FOREACH (QGLViewer *viewer, QGLViewer::QGLViewerPool()) {
    if (viewer)
      viewer->saveStateToFile();
  }
}

//////////////////////////////////////////////////////////////////////////
//       S a v e   s t a t e   b e t w e e n    s e s s i o n s         //
//////////////////////////////////////////////////////////////////////////

/*! Returns the state file name. Default value is \c .qglviewer.xml.

This is the name of the XML file where saveStateToFile() saves the viewer state
(camera state, widget geometry, display flags... see domElement()) on exit. Use
restoreStateFromFile() to restore this state later (usually in your init()
method).

Setting this value to \c QString::null will disable the automatic state file
saving that normally occurs on exit.

If more than one viewer are created by the application, this function will
return a numbered file name (as in ".qglviewer1.xml", ".qglviewer2.xml"... using
QGLViewer::QGLViewerIndex()) for extra viewers. Each viewer will then read back
its own information in restoreStateFromFile(), provided that the viewers are
created in the same order, which is usually the case. */
QString QGLViewer::stateFileName() const {
  QString name = stateFileName_;

  if (!name.isEmpty() && QGLViewer::QGLViewerIndex(this) > 0) {
    QFileInfo fi(name);
    if (fi.suffix().isEmpty())
      name += QString::number(QGLViewer::QGLViewerIndex(this));
    else
      name = fi.absolutePath() + '/' + fi.completeBaseName() +
             QString::number(QGLViewer::QGLViewerIndex(this)) + "." +
             fi.suffix();
  }

  return name;
}

/*! Saves in stateFileName() an XML representation of the QGLViewer state,
obtained from domElement().

Use restoreStateFromFile() to restore this viewer state.

This method is automatically called when a viewer is closed (using Escape or
using the window's upper right \c x close button). setStateFileName() to \c
QString::null to prevent this. */
void QGLViewer::saveStateToFile() {
  QString name = stateFileName();

  if (name.isEmpty())
    return;

  QFileInfo fileInfo(name);

  if (fileInfo.isDir()) {
    QMessageBox::warning(
        this, tr("Save to file error", "Message box window title"),
        tr("State file name (%1) references a directory instead of a file.")
            .arg(name));
    return;
  }

  const QString dirName = fileInfo.absolutePath();
  if (!QFileInfo(dirName).exists()) {
    QDir dir;
    if (!(dir.mkdir(dirName))) {
      QMessageBox::warning(this,
                           tr("Save to file error", "Message box window title"),
                           tr("Unable to create directory %1").arg(dirName));
      return;
    }
  }

  // Write the DOM tree to file
  QFile f(name);
  if (f.open(QIODevice::WriteOnly)) {
    QTextStream out(&f);
    QDomDocument doc("QGLVIEWER");
    doc.appendChild(domElement("QGLViewer", doc));
    doc.save(out, 2);
    f.flush();
    f.close();
  } else
    QMessageBox::warning(
        this, tr("Save to file error", "Message box window title"),
        tr("Unable to save to file %1").arg(name) + ":\n" + f.errorString());
}

/*! Restores the QGLViewer state from the stateFileName() file using
initFromDOMElement().

States are saved using saveStateToFile(), which is automatically called on
viewer exit.

Returns \c true when the restoration is successful. Possible problems are an non
existing or unreadable stateFileName() file, an empty stateFileName() or an XML
syntax error.

A manipulatedFrame() should be defined \e before calling this method, so that
its state can be restored. Initialization code put \e after this function will
override saved values: \code void Viewer::init()
{
// Default initialization goes here (including the declaration of a possible
manipulatedFrame).

if (!restoreStateFromFile())
showEntireScene(); // Previous state cannot be restored: fit camera to scene.

// Specific initialization that overrides file savings goes here.
}
\endcode */
bool QGLViewer::restoreStateFromFile() {
  QString name = stateFileName();

  if (name.isEmpty())
    return false;

  QFileInfo fileInfo(name);

  if (!fileInfo.isFile())
    // No warning since it would be displayed at first start.
    return false;

  if (!fileInfo.isReadable()) {
    QMessageBox::warning(
        this, tr("Problem in state restoration", "Message box window title"),
        tr("File %1 is not readable.").arg(name));
    return false;
  }

  // Read the DOM tree form file
  QFile f(name);
  if (f.open(QIODevice::ReadOnly)) {
    QDomDocument doc;
    doc.setContent(&f);
    f.close();
    QDomElement main = doc.documentElement();
    initFromDOMElement(main);
  } else {
    QMessageBox::warning(
        this, tr("Open file error", "Message box window title"),
        tr("Unable to open file %1").arg(name) + ":\n" + f.errorString());
    return false;
  }

  return true;
}

/*! Returns an XML \c QDomElement that represents the QGLViewer.

Used by saveStateToFile(). restoreStateFromFile() uses initFromDOMElement() to
restore the QGLViewer state from the resulting \c QDomElement.

\p name is the name of the QDomElement tag. \p doc is the \c QDomDocument
factory used to create QDomElement.

The created QDomElement contains state values (axisIsDrawn(), FPSIsDisplayed(),
isFullScreen()...), viewer geometry, as well as camera() (see
qglviewer::Camera::domElement()) and manipulatedFrame() (if defined, see
qglviewer::ManipulatedFrame::domElement()) states.

Overload this method to add your own attributes to the state file:
\code
QDomElement Viewer::domElement(const QString& name, QDomDocument& document)
const
{
// Creates a custom node for a light
QDomElement de = document.createElement("Light");
de.setAttribute("state", (lightIsOn()?"on":"off"));
// Note the include of the ManipulatedFrame domElement method.
de.appendChild(lightManipulatedFrame()->domElement("LightFrame", document));

// Get default state domElement and append custom node
QDomElement res = QGLViewer::domElement(name, document);
res.appendChild(de);
return res;
}
\endcode
See initFromDOMElement() for the associated restoration code.

\attention For the manipulatedFrame(), qglviewer::Frame::constraint() and
qglviewer::Frame::referenceFrame() are not saved. See
qglviewer::Frame::domElement(). */
QDomElement QGLViewer::domElement(const QString &name,
                                  QDomDocument &document) const {
  QDomElement de = document.createElement(name);
  de.setAttribute("version", QGLViewerVersionString());

  QDomElement stateNode = document.createElement("State");
  // hasMouseTracking() is not saved
  stateNode.appendChild(DomUtils::QColorDomElement(
      foregroundColor(), "foregroundColor", document));
  stateNode.appendChild(DomUtils::QColorDomElement(
      backgroundColor(), "backgroundColor", document));
  DomUtils::setBoolAttribute(stateNode, "stereo", displaysInStereo());
  // Revolve or fly camera mode is not saved
  de.appendChild(stateNode);

  QDomElement displayNode = document.createElement("Display");
  DomUtils::setBoolAttribute(displayNode, "axisIsDrawn", axisIsDrawn());
  DomUtils::setBoolAttribute(displayNode, "gridIsDrawn", gridIsDrawn());
  DomUtils::setBoolAttribute(displayNode, "FPSIsDisplayed", FPSIsDisplayed());
  DomUtils::setBoolAttribute(displayNode, "cameraIsEdited", cameraIsEdited());
  // textIsEnabled() is not saved
  de.appendChild(displayNode);

  QDomElement geometryNode = document.createElement("Geometry");
  DomUtils::setBoolAttribute(geometryNode, "fullScreen", isFullScreen());
  if (isFullScreen()) {
    geometryNode.setAttribute("prevPosX", QString::number(prevPos_.x()));
    geometryNode.setAttribute("prevPosY", QString::number(prevPos_.y()));
  } else {
    QWidget *tlw = topLevelWidget();
    geometryNode.setAttribute("width", QString::number(tlw->width()));
    geometryNode.setAttribute("height", QString::number(tlw->height()));
    geometryNode.setAttribute("posX", QString::number(tlw->pos().x()));
    geometryNode.setAttribute("posY", QString::number(tlw->pos().y()));
  }
  de.appendChild(geometryNode);

  // Restore original Camera zClippingCoefficient before saving.
  if (cameraIsEdited())
    camera()->setZClippingCoefficient(previousCameraZClippingCoefficient_);
  de.appendChild(camera()->domElement("Camera", document));
  if (cameraIsEdited())
    // #CONNECTION# 5.0 from setCameraIsEdited()
    camera()->setZClippingCoefficient(5.0);

  if (manipulatedFrame())
    de.appendChild(
        manipulatedFrame()->domElement("ManipulatedFrame", document));

  return de;
}

/*! Restores the QGLViewer state from a \c QDomElement created by domElement().

Used by restoreStateFromFile() to restore the QGLViewer state from a file.

Overload this method to retrieve custom attributes from the QGLViewer state
file. This code corresponds to the one given in the domElement() documentation:
\code
void Viewer::initFromDOMElement(const QDomElement& element)
{
// Restore standard state
QGLViewer::initFromDOMElement(element);

QDomElement child=element.firstChild().toElement();
while (!child.isNull())
{
if (child.tagName() == "Light")
{
if (child.hasAttribute("state"))
setLightOn(child.attribute("state").toLower() == "on");

// Assumes there is only one child. Otherwise you need to parse child's children
recursively. QDomElement lf = child.firstChild().toElement(); if (!lf.isNull()
&& lf.tagName() == "LightFrame")
lightManipulatedFrame()->initFromDomElement(lf);
}
child = child.nextSibling().toElement();
}
}
\endcode

See also qglviewer::Camera::initFromDOMElement(),
qglviewer::ManipulatedFrame::initFromDOMElement().

\note The manipulatedFrame() \e pointer is not modified by this method. If
defined, its state is simply set from the \p element values. */
void QGLViewer::initFromDOMElement(const QDomElement &element) {
  const QString version = element.attribute("version");
  // if (version != QGLViewerVersionString())
  if (version[0] != '2')
    // Patches for previous versions should go here when the state file syntax
    // is modified.
    qWarning("State file created using QGLViewer version %s may not be "
             "correctly read.",
             version.toLatin1().constData());

  QDomElement child = element.firstChild().toElement();
  bool tmpCameraIsEdited = cameraIsEdited();
  while (!child.isNull()) {
    if (child.tagName() == "State") {
      // #CONNECTION# default values from defaultConstructor()
      // setMouseTracking(DomUtils::boolFromDom(child, "mouseTracking", false));
      setStereoDisplay(DomUtils::boolFromDom(child, "stereo", false));
      // if ((child.attribute("cameraMode", "revolve") == "fly") &&
      // (cameraIsInRevolveMode())) 	toggleCameraMode();

      QDomElement ch = child.firstChild().toElement();
      while (!ch.isNull()) {
        if (ch.tagName() == "foregroundColor")
          setForegroundColor(DomUtils::QColorFromDom(ch));
        if (ch.tagName() == "backgroundColor")
          setBackgroundColor(DomUtils::QColorFromDom(ch));
        ch = ch.nextSibling().toElement();
      }
    }

    if (child.tagName() == "Display") {
      // #CONNECTION# default values from defaultConstructor()
      setAxisIsDrawn(DomUtils::boolFromDom(child, "axisIsDrawn", false));
      setGridIsDrawn(DomUtils::boolFromDom(child, "gridIsDrawn", false));
      setFPSIsDisplayed(DomUtils::boolFromDom(child, "FPSIsDisplayed", false));
      // See comment below.
      tmpCameraIsEdited = DomUtils::boolFromDom(child, "cameraIsEdited", false);
      // setTextIsEnabled(DomUtils::boolFromDom(child, "textIsEnabled", true));
    }

    if (child.tagName() == "Geometry") {
      setFullScreen(DomUtils::boolFromDom(child, "fullScreen", false));

      if (isFullScreen()) {
        prevPos_.setX(DomUtils::intFromDom(child, "prevPosX", 0));
        prevPos_.setY(DomUtils::intFromDom(child, "prevPosY", 0));
      } else {
        int width = DomUtils::intFromDom(child, "width", 600);
        int height = DomUtils::intFromDom(child, "height", 400);
        topLevelWidget()->resize(width, height);
        camera()->setScreenWidthAndHeight(this->width(), this->height());

        QPoint pos;
        pos.setX(DomUtils::intFromDom(child, "posX", 0));
        pos.setY(DomUtils::intFromDom(child, "posY", 0));
        topLevelWidget()->move(pos);
      }
    }

    if (child.tagName() == "Camera") {
      connectAllCameraKFIInterpolatedSignals(false);
      camera()->initFromDOMElement(child);
      connectAllCameraKFIInterpolatedSignals();
    }

    if ((child.tagName() == "ManipulatedFrame") && (manipulatedFrame()))
      manipulatedFrame()->initFromDOMElement(child);

    child = child.nextSibling().toElement();
  }

  // The Camera always stores its "real" zClippingCoef in domElement(). If it is
  // edited, its "real" coef must be saved and the coef set to 5.0, as is done
  // in setCameraIsEdited(). BUT : Camera and Display are read in an arbitrary
  // order. We must initialize Camera's "real" coef BEFORE calling
  // setCameraIsEdited. Hence this temp cameraIsEdited and delayed call
  cameraIsEdited_ = tmpCameraIsEdited;
  if (cameraIsEdited_) {
    previousCameraZClippingCoefficient_ = camera()->zClippingCoefficient();
    // #CONNECTION# 5.0 from setCameraIsEdited.
    camera()->setZClippingCoefficient(5.0);
  }
}

#ifndef DOXYGEN
/*! This method is deprecated since version 1.3.9-5. Use saveStateToFile() and
setStateFileName() instead. */
void QGLViewer::saveToFile(const QString &fileName) {
  if (!fileName.isEmpty())
    setStateFileName(fileName);

  qWarning("saveToFile() is deprecated, use saveStateToFile() instead.");
  saveStateToFile();
}

/*! This function is deprecated since version 1.3.9-5. Use
restoreStateFromFile() and setStateFileName() instead. */
bool QGLViewer::restoreFromFile(const QString &fileName) {
  if (!fileName.isEmpty())
    setStateFileName(fileName);

  qWarning(
      "restoreFromFile() is deprecated, use restoreStateFromFile() instead.");
  return restoreStateFromFile();
}
#endif


void QGLViewer::copyBufferToTexture(GLint , GLenum ) {
}

/*! Returns the texture id of the texture created by copyBufferToTexture().

Use glBindTexture() to use this texture. Note that this is already done by
copyBufferToTexture().

Returns \c 0 is copyBufferToTexture() was never called or if the texure was
deleted using glDeleteTextures() since then. */
GLuint QGLViewer::bufferTextureId() const {
    return 0;
}
