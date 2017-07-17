#include "Viewer.h"
#include <CGAL/gl.h>
#include <CGAL/Three/Scene_draw_interface.h>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QGLViewer/manipulatedCameraFrame.h>
#include <QDebug>
#include <QOpenGLShader>
#include <QFileDialog>
#include <QOpenGLShaderProgram>
#include <QOpenGLFramebufferObject>
#include <QMessageBox>
#include <QColorDialog>
#include <QInputDialog>
#include <cmath>
#include <QApplication>

#if defined(_WIN32)
#include <QMimeData>
#include <QByteArray>
#include <QBuffer>
#endif

class Viewer_impl {
public:
  CGAL::Three::Scene_draw_interface* scene;
  Viewer *viewer;
  bool antialiasing;
  bool twosides;
  bool macro_mode;
  bool inFastDrawing;
  bool inDrawWithNames;
  QPainter *painter;
  // F P S    d i s p l a y
  QTime fpsTime;
  unsigned int fpsCounter;
  QString fpsString;
  float f_p_s;
  // M e s s a g e s
  QString message;
  bool _displayMessage;
  QTimer messageTimer;

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
  bool i_is_pressed;
  bool z_is_pressed;
  //!Draws the distance between two selected points.
  void showDistance(QPoint);
  qglviewer::Vec APoint;
  qglviewer::Vec BPoint;
  qglviewer::Vec offset;
  bool is_d_pressed;
  bool extension_is_found;

  TextRenderer *textRenderer;
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
  //!Clears the distance display
  void clearDistancedisplay();
  void draw_aux(bool with_names, Viewer*);
  //! Contains all the programs for the item rendering.
  mutable std::vector<QOpenGLShaderProgram*> shader_programs;
  QMatrix4x4 projectionMatrix;
  void setFrustum(double l, double r, double t, double b, double n, double f);
  QImage* takeSnapshot(Viewer* viewer, int quality, int background_color, QSize size, double oversampling, bool expand);
  void sendSnapshotToClipboard(Viewer*);
};
Viewer::Viewer(QWidget* parent, bool antialiasing)
  : CGAL::Three::Viewer_interface(parent)
{
  d = new Viewer_impl;
  d->scene = 0;
  d->antialiasing = antialiasing;
  d->twosides = false;
  this->setProperty("draw_two_sides", false);
  d->macro_mode = false;
  d->inFastDrawing = true;
  d->inDrawWithNames = false;
  d->shader_programs.resize(NB_OF_PROGRAMS);
  d->offset = qglviewer::Vec(0,0,0);
  d->textRenderer = new TextRenderer();
  connect( d->textRenderer, SIGNAL(sendMessage(QString,int)),
           this, SLOT(printMessage(QString,int)) );
  connect(&d->messageTimer, SIGNAL(timeout()), SLOT(hideMessage()));
  setShortcut(EXIT_VIEWER, 0);
  setShortcut(DRAW_AXIS, 0);
  setKeyDescription(Qt::Key_T,
                    tr("Turn the camera by 180 degrees"));
  setKeyDescription(Qt::Key_M,
                    tr("Toggle macro mode: useful to view details very near from the camera, "
                       "but decrease the z-buffer precision"));
  setKeyDescription(Qt::Key_A,
                      tr("Toggle the axis system visibility."));
  setKeyDescription(Qt::Key_I + Qt::CTRL,
                      tr("Toggle the primitive IDs visibility of the selected Item."));
  setKeyDescription(Qt::Key_D,
                      tr("Disable the distance between two points  visibility."));

#if QGLVIEWER_VERSION >= 0x020501
  //modify mouse bindings that have been updated
  setMouseBinding(Qt::Key(0), Qt::NoModifier, Qt::LeftButton, RAP_FROM_PIXEL, true, Qt::RightButton);
  setMouseBindingDescription(Qt::ShiftModifier, Qt::RightButton,
                             tr("Select and pop context menu"));
  setMouseBinding(Qt::Key_R, Qt::NoModifier, Qt::LeftButton, RAP_FROM_PIXEL);
  //use the new API for these
  setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, SELECT);

  setMouseBindingDescription(Qt::Key(0), Qt::ShiftModifier, Qt::LeftButton,
                             tr("Selects and display context "
                                "menu of the selected item"));
  setMouseBindingDescription(Qt::Key_I, Qt::NoModifier, Qt::LeftButton,
                             tr("Show/hide the primitive ID."));
  setMouseBindingDescription(Qt::Key_D, Qt::NoModifier, Qt::LeftButton,
                             tr("Selects a point. When the second point is selected,  "
                                "displays the two points and the distance between them."));
  setMouseBindingDescription(Qt::Key_O, Qt::NoModifier, Qt::LeftButton,
                             tr("Move the camera orthogonally to the picked facet of a Scene_polyhedron_item or "
                                "to the current selection of a Scene_points_with_normal_item."));
#else
  setMouseBinding(Qt::SHIFT + Qt::LeftButton, SELECT);
  setMouseBindingDescription(Qt::SHIFT + Qt::RightButton,
                             tr("Selects and display context "
                                "menu of the selected item"));

#endif // QGLVIEWER_VERSION >= 2.5.0
  prev_radius = sceneRadius();
  d->axis_are_displayed = true;
  d->has_text = false;
  d->i_is_pressed = false;
  d->z_is_pressed = false;
  d->fpsTime.start();
  d->fpsCounter=0;
  d->f_p_s=0.0;
  d->fpsString=tr("%1Hz", "Frames per seconds, in Hertz").arg("?");
  d->distance_is_displayed = false;
  d->is_d_pressed = false;
  d->viewer = this;
}

Viewer::~Viewer()
{
  delete d;
}

void Viewer::setScene(CGAL::Three::Scene_draw_interface* scene)
{
  d->scene = scene;
}

bool Viewer::antiAliasing() const
{
  return d->antialiasing; 
}

void Viewer::setAntiAliasing(bool b)
{
  d->antialiasing = b;
  update();
}

void Viewer::setTwoSides(bool b)
{
  this->setProperty("draw_two_sides", b);
  d->twosides = b;
  update();
}


void Viewer::setFastDrawing(bool b)
{
  d->inFastDrawing = b;
  update();
}

bool Viewer::inFastDrawing() const
{
  return (d->inFastDrawing
          && (camera()->frame()->isSpinning()
              || camera()->frame()->isManipulated()));
}

void Viewer::draw()
{ 
  makeCurrent();
  glEnable(GL_DEPTH_TEST);
  d->draw_aux(false, this);
}

void Viewer::fastDraw()
{
  d->draw_aux(false, this);
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  initializeOpenGLFunctions();
  glDrawArraysInstanced = (PFNGLDRAWARRAYSINSTANCEDARBPROC)this->context()->getProcAddress("glDrawArraysInstancedARB");
  if(!glDrawArraysInstanced)
  {
      qDebug()<<"glDrawArraysInstancedARB : extension not found. Spheres will be displayed as points.";
      d->extension_is_found = false;
  }
  else
      d->extension_is_found = true;

  glVertexAttribDivisor = (PFNGLVERTEXATTRIBDIVISORARBPROC)this->context()->getProcAddress("glVertexAttribDivisorARB");
  if(!glDrawArraysInstanced)
  {
      qDebug()<<"glVertexAttribDivisorARB : extension not found. Spheres will be displayed as points.";
      d->extension_is_found = false;
  }
  else
      d->extension_is_found = true;


  setBackgroundColor(::Qt::white);
  d->vao[0].create();
  for(int i=0; i<3; i++)
    d->buffers[i].create();

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
      "   color = colors; \n"
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
      "uniform highp vec4 light_pos;  \n"
      "uniform highp vec4 light_diff; \n"
      "uniform highp vec4 light_spec; \n"
      "uniform highp vec4 light_amb;  \n"
      "uniform highp float spec_power ; \n"

      "void main(void) { \n"

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
      "} \n"
      "\n"
      };

      QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader->compileSourceCode(vertex_source))
      {
          std::cerr<<"Compiling vertex source FAILED"<<std::endl;
      }

      QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader->compileSourceCode(fragment_source))
      {
          std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
      }

      if(!d->rendering_program.addShader(vertex_shader))
      {
          std::cerr<<"adding vertex shader FAILED"<<std::endl;
      }
      if(!d->rendering_program.addShader(fragment_shader))
      {
          std::cerr<<"adding fragment shader FAILED"<<std::endl;
      }
      d->rendering_program.bindAttributeLocation("colors", 1);
      if(!d->rendering_program.link())
      {
          //std::cerr<<"linking Program FAILED"<<std::endl;
          qDebug() << d->rendering_program.log();
      }
  //setting the program used for the distance
     {
         d->vao[1].create();
         d->buffers[3].create();
         //Vertex source code
         const char vertex_source_dist[] =
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
         const char fragment_source_dist[] =
         {
             "#version 120 \n"
             "void main(void) { \n"
             "gl_FragColor = vec4(0.0,0.0,0.0,1.0); \n"
             "} \n"
             "\n"
         };
         vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
         if(!vertex_shader->compileSourceCode(vertex_source_dist))
         {
             std::cerr<<"Compiling vertex source FAILED"<<std::endl;
         }

         fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
         if(!fragment_shader->compileSourceCode(fragment_source_dist))
         {
             std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
         }

         if(!d->rendering_program_dist.addShader(vertex_shader))
         {
             std::cerr<<"adding vertex shader FAILED"<<std::endl;
         }
         if(!d->rendering_program_dist.addShader(fragment_shader))
         {
             std::cerr<<"adding fragment shader FAILED"<<std::endl;
         }
         if(!d->rendering_program_dist.link())
         {
             qDebug() << d->rendering_program_dist.log();
         }
     }

      Viewer_impl::AxisData data;
      d->v_Axis.resize(0);
      d->n_Axis.resize(0);
      d->c_Axis.resize(0);
      data.vertices = &d->v_Axis;
      data.normals =  &d->n_Axis;
      data.colors =   &d->c_Axis;
      GLdouble l = 1.0;
      d->makeArrow(0.06,10, qglviewer::Vec(0,0,0),qglviewer::Vec(l,0,0),qglviewer::Vec(1,0,0), data);
      d->makeArrow(0.06,10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,l,0),qglviewer::Vec(0,1,0), data);
      d->makeArrow(0.06,10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,0,l),qglviewer::Vec(0,0,1), data);

      d->rendering_program.bind();
      d->vao[0].bind();
      d->buffers[0].bind();
      d->buffers[0].allocate(d->v_Axis.data(), static_cast<int>(d->v_Axis.size()) * sizeof(float));
      d->rendering_program.enableAttributeArray("vertex");
      d->rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
      d->buffers[0].release();

      d->buffers[1].bind();
      d->buffers[1].allocate(d->n_Axis.data(), static_cast<int>(d->n_Axis.size() * sizeof(float)));
      d->rendering_program.enableAttributeArray("normal");
      d->rendering_program.setAttributeBuffer("normal",GL_FLOAT,0,3);
      d->buffers[1].release();

      d->buffers[2].bind();
      d->buffers[2].allocate(d->c_Axis.data(), static_cast<int>(d->c_Axis.size() * sizeof(float)));
      d->rendering_program.enableAttributeArray("colors");
      d->rendering_program.setAttributeBuffer("colors",GL_FLOAT,0,3);
      d->buffers[2].release();
      d->vao[0].release();

      d->rendering_program.release();

  d->painter = new QPainter(this);
}

#include <QMouseEvent>

void Viewer::mousePressEvent(QMouseEvent* event)
{
  if(event->button() == Qt::RightButton &&
     event->modifiers().testFlag(Qt::ShiftModifier)) 
  {
    select(event->pos());
    requestContextMenu(event->globalPos());
    event->accept();
  }
  else if(!event->modifiers()
          && event->button() == Qt::LeftButton
          && d->i_is_pressed)
  {
      d->scene->printPrimitiveId(event->pos(), this);
  }
  else if(!event->modifiers()
          && event->button() == Qt::LeftButton
          && d->z_is_pressed)
  {
      d->scene->zoomToPosition(event->pos(), this);
  }
  else if(!event->modifiers()
          && event->button() == Qt::LeftButton
          && d->is_d_pressed)
  {
      d->showDistance(event->pos());
      event->accept();
  }
  else {
    QGLViewer::mousePressEvent(event);
  }
}

#include <QContextMenuEvent>
void Viewer::contextMenuEvent(QContextMenuEvent* event)
{
  if(event->reason() != QContextMenuEvent::Mouse) {
    requestContextMenu(event->globalPos());
    event->accept();
  }
  else {
    QGLViewer::contextMenuEvent(event);
  }
}

void Viewer::keyPressEvent(QKeyEvent* e)
{
  if(!e->modifiers()) {
    if(e->key() == Qt::Key_T) {
      turnCameraBy180Degres();
      return;
    }
    else if(e->key() == Qt::Key_M) {
      d->macro_mode = ! d->macro_mode;

      if(d->macro_mode) {
          camera()->setZNearCoefficient(0.0005f);
      } else {
        camera()->setZNearCoefficient(0.005f);
      }
      this->displayMessage(tr("Macro mode: %1").
                           arg(d->macro_mode ? tr("on") : tr("off")));



      return;
    }
    else if(e->key() == Qt::Key_A) {
          d->axis_are_displayed = !d->axis_are_displayed;
          update();
        }
    else if(e->key() == Qt::Key_I) {
          d->i_is_pressed = true;
        }
    else if(e->key() == Qt::Key_O) {
          d->z_is_pressed = true;
        }
    else if(e->key() == Qt::Key_D) {
        if(e->isAutoRepeat())
        {
            return;
        }
        if(!d->is_d_pressed)
        {
            d->clearDistancedisplay();
        }
        d->is_d_pressed = true;
        update();
        return;
    }
  }
  else if(e->key() == Qt::Key_I && e->modifiers() & Qt::ControlModifier){
    d->scene->printPrimitiveIds(this);
    update();
    return;
  }

  else if(e->key() == Qt::Key_C && e->modifiers() & Qt::ControlModifier){
    d->sendSnapshotToClipboard(this);
    return;
  }

  else if(e->key() == Qt::Key_S && e->modifiers() & Qt::ControlModifier){
    this->saveSnapshot(true,true);
    return;
  }

  //forward the event to the scene (item handling of the event)
  if (! d->scene->keyPressEvent(e) )
    QGLViewer::keyPressEvent(e);
}

void Viewer::keyReleaseEvent(QKeyEvent *e)
{
  if(e->key() == Qt::Key_I) {
    d->i_is_pressed = false;
  }
  else if(e->key() == Qt::Key_O) {
    d->z_is_pressed = false;
  }
  else if(!e->modifiers() && e->key() == Qt::Key_D)
  {
    if(e->isAutoRepeat())
    {
      return;
    }
    d->is_d_pressed = false;
  }
  QGLViewer::keyReleaseEvent(e);
}

void Viewer::turnCameraBy180Degres() {
  qglviewer::Camera* camera = this->camera();
  using qglviewer::ManipulatedCameraFrame;

  ManipulatedCameraFrame frame_from(*camera->frame());
  camera->setViewDirection(-camera->viewDirection());
  ManipulatedCameraFrame frame_to(*camera->frame());

  camera->setOrientation(frame_from.orientation());
  camera->interpolateTo(frame_to, 0.5f);
}

void Viewer_impl::draw_aux(bool with_names, Viewer* viewer)
{
  if(scene == 0)
    return;
  viewer->glLineWidth(1.0f);
  viewer->glPointSize(2.f);
  viewer->glEnable(GL_POLYGON_OFFSET_FILL);
  viewer->glPolygonOffset(1.0f,1.0f);
  viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  viewer->glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

  if(twosides)
    viewer->glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    viewer->glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  if(!with_names && antialiasing)
  {
    viewer->glEnable(GL_BLEND);
    viewer->glEnable(GL_LINE_SMOOTH);
    viewer->glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    viewer->glDisable(GL_BLEND);
    viewer->glDisable(GL_LINE_SMOOTH);
    viewer->glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    viewer->glBlendFunc(GL_ONE, GL_ZERO);
  }
  inDrawWithNames = with_names;
  if(with_names)
    scene->drawWithNames(viewer);
  else
    scene->draw(viewer);
  viewer->glDisable(GL_POLYGON_OFFSET_FILL);
  viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
}

bool Viewer::inDrawWithNames() const {
  return d->inDrawWithNames;
}

void Viewer::drawWithNames()
{
  QGLViewer::draw();
  d->draw_aux(true, this);
}

void Viewer::postSelection(const QPoint& pixel)
{
  Q_EMIT selected(this->selectedName());
  bool found = false;
  qglviewer::Vec point = camera()->pointUnderPixel(pixel, found) - d->offset;
  if(found) {
    Q_EMIT selectedPoint(point.x,
                       point.y,
                       point.z);
    const qglviewer::Vec orig = camera()->position() - d->offset;
    const qglviewer::Vec dir = point - orig;
    Q_EMIT selectionRay(orig.x, orig.y, orig.z,
                      dir.x, dir.y, dir.z);
  }
}
bool CGAL::Three::Viewer_interface::readFrame(QString s, qglviewer::Frame& frame)
{
  QStringList list = s.split(" ", QString::SkipEmptyParts);
  if(list.size() != 7)
    return false;
  float vec[3];
  for(int i = 0; i < 3; ++i)
  {
    bool ok;
    vec[i] = list[i].toFloat(&ok);
    if(!ok) return false;
  }
  double orient[4];
  for(int i = 0; i < 4; ++i)
  {
    bool ok;
    orient[i] = list[i + 3].toDouble(&ok);
    if(!ok) return false;
  }
  frame.setPosition(qglviewer::Vec(vec[0],
                                   vec[1],
                                   vec[2]));
  frame.setOrientation(orient[0],
                       orient[1],
                       orient[2],
                       orient[3]);
  return true;
}

QString CGAL::Three::Viewer_interface::dumpFrame(const qglviewer::Frame& frame) {
  const qglviewer::Vec pos = frame.position();
  const qglviewer::Quaternion q = frame.orientation();

  return QString("%1 %2 %3 %4 %5 %6 %7")
    .arg(pos[0])
    .arg(pos[1])
    .arg(pos[2])
    .arg(q[0])
    .arg(q[1])
    .arg(q[2])
    .arg(q[3]);
}

bool Viewer::moveCameraToCoordinates(QString s, float animation_duration) {
  qglviewer::Frame new_frame;
  if(readFrame(s, new_frame)) {
    camera()->interpolateTo(new_frame, animation_duration); 
    return true;
  }
  else
    return false;
}

QString Viewer::dumpCameraCoordinates()
{
  if(camera()->frame()) {
    return dumpFrame(*camera()->frame());
  } else {
    return QString();
  }
}

void Viewer::attribBuffers(int program_name) const {
    GLint is_both_sides = 0;
    //ModelViewMatrix used for the transformation of the camera.
    QMatrix4x4 mvp_mat;
    // ModelView Matrix used for the lighting system
    QMatrix4x4 mv_mat;
    // transformation of the manipulated frame
    QMatrix4x4 f_mat;

    f_mat.setToIdentity();
    //fills the MVP and MV matrices.
    GLdouble d_mat[16];

    this->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mv_mat.data()[i] = GLfloat(d_mat[i]);
    mvp_mat = d->projectionMatrix*mv_mat;

    const_cast<Viewer*>(this)->glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE,
                                             &is_both_sides);

    QVector4D position(0.0f,0.0f,1.0f, 1.0f );
    QVector4D ambient(0.4f, 0.4f, 0.4f, 0.4f);
    // Diffuse
    QVector4D diffuse(1.0f, 1.0f, 1.0f, 1.0f);
    // Specular
    QVector4D specular(0.0f, 0.0f, 0.0f, 1.0f);
    QOpenGLShaderProgram* program = getShaderProgram(program_name);
    program->bind();
    program->setUniformValue("mvp_matrix", mvp_mat);
    switch(program_name)
    {
    case PROGRAM_WITH_LIGHT:
    case PROGRAM_C3T3:
    case PROGRAM_PLANE_TWO_FACES:
    case PROGRAM_INSTANCED:
    case PROGRAM_WITH_TEXTURE:
    case PROGRAM_CUTPLANE_SPHERES:
    case PROGRAM_SPHERES:
    case PROGRAM_C3T3_TETS:
    case PROGRAM_FLAT:
        program->setUniformValue("light_pos", position);
        program->setUniformValue("light_diff",diffuse);
        program->setUniformValue("light_spec", specular);
        program->setUniformValue("light_amb", ambient);
        program->setUniformValue("spec_power", 51.8f);
        program->setUniformValue("is_two_side", is_both_sides);
        break;
    }
    switch(program_name)
    {
    case PROGRAM_WITH_LIGHT:
    case PROGRAM_C3T3:
    case PROGRAM_PLANE_TWO_FACES:
    case PROGRAM_INSTANCED:
    case PROGRAM_CUTPLANE_SPHERES:
    case PROGRAM_SPHERES:
    case PROGRAM_C3T3_TETS:
    case PROGRAM_FLAT:
      program->setUniformValue("mv_matrix", mv_mat);
      break;
    case PROGRAM_WITHOUT_LIGHT:
      program->setUniformValue("f_matrix",f_mat);
      break;
    case PROGRAM_WITH_TEXTURE:
      program->setUniformValue("mv_matrix", mv_mat);
      program->setUniformValue("s_texture",0);
      program->setUniformValue("f_matrix",f_mat);
      break;
    case PROGRAM_WITH_TEXTURED_EDGES:
        program->setUniformValue("s_texture",0);
        break;
    case PROGRAM_NO_SELECTION:
        program->setUniformValue("f_matrix",f_mat);
        break;
    }
    program->release();
}

void Viewer::beginSelection(const QPoint &point)
{
    makeCurrent();
    glEnable(GL_SCISSOR_TEST);
    glScissor(point.x(), camera()->screenHeight()-1-point.y(), 1, 1);
    d->scene->setPickedPixel(point);
}
void Viewer::endSelection(const QPoint&)
{
    glDisable(GL_SCISSOR_TEST);
    //redraw the true scene for the glReadPixel in postSelection();
    d->draw_aux(false, this);
}

void Viewer_impl::makeArrow(double R, int prec, qglviewer::Vec from, qglviewer::Vec to, qglviewer::Vec color, Viewer_impl::AxisData &data)
{
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
        QVector4D n(Rf*2.*sin(D), sin(a), Rf*2.*cos(D), 1.);
        QVector4D pR = mat*p;
        QVector4D nR = mat*n;

        //point A1
        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back((float)color.x);
        data.colors->push_back((float)color.y);
        data.colors->push_back((float)color.z);

        //point B1
        p = QVector4D(Rf*2.*sin(D), 0.66f, Rf*2.* cos(D), 1.f);
        n = QVector4D(sin(D), sin(a), cos(D), 1.);
        pR = mat*p;
        nR = mat*n;
        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back((float)color.x);
        data.colors->push_back((float)color.y);
        data.colors->push_back((float)color.z);
        //point C1
        D = (d+360/prec)*M_PI/180.0;
        p = QVector4D(Rf*2.* sin(D), 0.66f, Rf *2.* cos(D), 1.f);
        n = QVector4D(sin(D), sin(a), cos(D), 1.0);
        pR = mat*p;
        nR = mat*n;

        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back((float)color.x);
        data.colors->push_back((float)color.y);
        data.colors->push_back((float)color.z);

    }

    //cylinder
    //body of the cylinder
    for(int d = 0; d<360; d+= 360/prec)
    {
        //point A1
        double D = d*M_PI/180.0;
        QVector4D p(Rf*sin(D), 0.66f, Rf*cos(D), 1.f);
        QVector4D n(sin(D), 0.f, cos(D), 1.f);
        QVector4D pR = mat*p;
        QVector4D nR = mat*n;

        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back(color.x);
        data.colors->push_back(color.y);
        data.colors->push_back(color.z);
        //point B1
        p = QVector4D(Rf * sin(D),0,Rf*cos(D), 1.0);
        n = QVector4D(sin(D), 0, cos(D), 1.0);
        pR = mat*p;
        nR = mat*n;


        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back(color.x);
        data.colors->push_back(color.y);
        data.colors->push_back(color.z);
          //point C1
        D = (d+360/prec)*M_PI/180.0;
        p = QVector4D(Rf * sin(D),0,Rf*cos(D), 1.0);
        n = QVector4D(sin(D), 0, cos(D), 1.0);
        pR = mat*p;
        nR = mat*n;
        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back(color.x);
        data.colors->push_back(color.y);
        data.colors->push_back(color.z);
        //point A2
        D = (d+360/prec)*M_PI/180.0;

        p = QVector4D(Rf * sin(D),0,Rf*cos(D), 1.0);
        n = QVector4D(sin(D), 0, cos(D), 1.0);
        pR = mat*p;
        nR = mat*n;
        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back((float)color.x);
        data.colors->push_back((float)color.y);
        data.colors->push_back((float)color.z);
        //point B2
        p = QVector4D(Rf * sin(D), 0.66f, Rf*cos(D), 1.f);
        n = QVector4D(sin(D), 0, cos(D), 1.0);
        pR = mat*p;
        nR = mat*n;
        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back((float)color.x);
        data.colors->push_back((float)color.y);
        data.colors->push_back((float)color.z);
        //point C2
        D = d*M_PI/180.0;
        p = QVector4D(Rf * sin(D), 0.66f, Rf*cos(D), 1.f);
        n = QVector4D(sin(D), 0.f, cos(D), 1.f);
        pR = mat*p;
        nR = mat*n;
        data.vertices->push_back(pR.x());
        data.vertices->push_back(pR.y());
        data.vertices->push_back(pR.z());
        data.normals->push_back(nR.x());
        data.normals->push_back(nR.y());
        data.normals->push_back(nR.z());
        data.colors->push_back(color.x);
        data.colors->push_back(color.y);
        data.colors->push_back(color.z);

    }
}

void Viewer::drawVisualHints()
{

    QGLViewer::drawVisualHints();
    if(d->axis_are_displayed)
    {
      d->rendering_program.bind();
      qglviewer::Camera::Type camera_type = camera()->type();
      camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
        QMatrix4x4 mvpMatrix;
        QMatrix4x4 mvMatrix;
        for(int i=0; i < 16; i++)
        {
            mvMatrix.data()[i] = camera()->orientation().inverse().matrix()[i];
        }
        mvpMatrix.ortho(-1,1,-1,1,-1,1);
        mvpMatrix = mvpMatrix*mvMatrix;
        camera()->setType(camera_type);
        QVector4D	position(0.0f,0.0f,1.0f,1.0f );
        // define material
        QVector4D	ambient;
        QVector4D	diffuse;
        QVector4D	specular;
        GLfloat      shininess ;
        // Ambient
        ambient[0] = 0.29225f;
        ambient[1] = 0.29225f;
        ambient[2] = 0.29225f;
        ambient[3] = 1.0f;
        // Diffuse
        diffuse[0] = 0.50754f;
        diffuse[1] = 0.50754f;
        diffuse[2] = 0.50754f;
        diffuse[3] = 1.0f;
        // Specular
        specular[0] = 0.0f;
        specular[1] = 0.0f;
        specular[2] = 0.0f;
        specular[3] = 0.0f;
        // Shininess
        shininess = 51.2f;

        d->rendering_program.setUniformValue("light_pos", position);
        d->rendering_program.setUniformValue("mvp_matrix", mvpMatrix);
        d->rendering_program.setUniformValue("mv_matrix", mvMatrix);
        d->rendering_program.setUniformValue("light_diff", diffuse);
        d->rendering_program.setUniformValue("light_spec", specular);
        d->rendering_program.setUniformValue("light_amb", ambient);
        d->rendering_program.setUniformValue("spec_power", shininess);

        d->vao[0].bind();
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
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->v_Axis.size() / 3));
        // The viewport and the scissor are restored.
        glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
        glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
        d->vao[0].release();
        d->rendering_program.release();
    }

    if(d->distance_is_displayed)
    {
        glDisable(GL_DEPTH_TEST);

        glLineWidth(3.0f);
        glPointSize(6.0f);
        //draws the distance
        QMatrix4x4 mvpMatrix;
        double mat[16];
        //camera()->frame()->rotation().getMatrix(mat);
        camera()->getModelViewProjectionMatrix(mat);
        //nullifies the translation
        for(int i=0; i < 16; i++)
        {
            mvpMatrix.data()[i] = (float)mat[i];
        }
        d->rendering_program_dist.bind();
        d->rendering_program_dist.setUniformValue("mvp_matrix", mvpMatrix);
        d->vao[1].bind();
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(2));
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(2));
        d->vao[1].release();
        d->rendering_program_dist.release();
        glEnable(GL_DEPTH_TEST);
        glPointSize(1.0f);
        glLineWidth(1.0f);

    }
    if (!d->painter->isActive())
      d->painter->begin(this);
    //So that the text is drawn in front of everything
    d->painter->beginNativePainting();
    glDisable(GL_DEPTH_TEST);
    d->painter->endNativePainting();
    //prints FPS
    TextItem *fps_text = new TextItem(20, int(1.5*((QApplication::font().pixelSize()>0)?QApplication::font().pixelSize():QApplication::font().pointSize())),0,d->fpsString,false, QFont(), Qt::gray);
    if(FPSIsDisplayed())
    {
      d->textRenderer->addText(fps_text);
    }
    //Prints the displayMessage
    QFont font = QFont();
    QFontMetrics fm(font);
    TextItem *message_text = new TextItem(10 + fm.width(d->message)/2, height()-20, 0, d->message, false, QFont(), Qt::gray );
    if (d->_displayMessage)
    {
      d->textRenderer->addText(message_text);
    }
    d->textRenderer->draw(this);
    if(FPSIsDisplayed())
      d->textRenderer->removeText(fps_text);
    if (d->_displayMessage)
      d->textRenderer->removeText(message_text);
}

void Viewer::resizeGL(int w, int h)
{
    QGLViewer::resizeGL(w,h);
}
QOpenGLShaderProgram* Viewer::declare_program(int name,
                                      const char* v_shader,
                                      const char* f_shader) const
{
  // workaround constness issues in Qt
  Viewer* viewer = const_cast<Viewer*>(this);

  if(d->shader_programs[name])
  {
    return d->shader_programs[name];
  }

  else
  {

    QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
    if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,v_shader))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,f_shader))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    program->bindAttributeLocation("colors", 1);
    program->link();
    d->shader_programs[name] = program;
    return program;
  }
}
QOpenGLShaderProgram* Viewer::getShaderProgram(int name) const
{
    switch(name)
    {
    case PROGRAM_C3T3:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3.v" , ":/cgal/Polyhedron_3/resources/shader_c3t3.f");
        break;
    case PROGRAM_C3T3_EDGES:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3_edges.v" , ":/cgal/Polyhedron_3/resources/shader_c3t3_edges.f");
        break;
    case PROGRAM_WITH_LIGHT:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_light.v" , ":/cgal/Polyhedron_3/resources/shader_with_light.f");
        break;
    case PROGRAM_WITHOUT_LIGHT:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_without_light.v" , ":/cgal/Polyhedron_3/resources/shader_without_light.f");
       break;
    case PROGRAM_NO_SELECTION:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_without_light.v" , ":/cgal/Polyhedron_3/resources/shader_no_light_no_selection.f");
        break;
    case PROGRAM_WITH_TEXTURE:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_texture.v" , ":/cgal/Polyhedron_3/resources/shader_with_texture.f");
      break;
    case PROGRAM_PLANE_TWO_FACES:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_without_light.v" , ":/cgal/Polyhedron_3/resources/shader_plane_two_faces.f");
        break;
    case PROGRAM_WITH_TEXTURED_EDGES:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_textured_edges.v" , ":/cgal/Polyhedron_3/resources/shader_with_textured_edges.f");
        break;
    case PROGRAM_INSTANCED:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_instanced.v" , ":/cgal/Polyhedron_3/resources/shader_with_light.f");
        break;
    case PROGRAM_INSTANCED_WIRE:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_instanced.v" , ":/cgal/Polyhedron_3/resources/shader_without_light.f");
        break;
    case PROGRAM_CUTPLANE_SPHERES:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3_spheres.v" , ":/cgal/Polyhedron_3/resources/shader_c3t3.f");
     break;
    case PROGRAM_C3T3_TETS:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3_tets.v" , ":/cgal/Polyhedron_3/resources/shader_with_light.f");
     break;
    case PROGRAM_SPHERES:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_spheres.v" , ":/cgal/Polyhedron_3/resources/shader_with_light.f");
      break;
    case PROGRAM_FLAT:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_light.v", ":/cgal/Polyhedron_3/resources/shader_flat.f");
      break;

    default:
        std::cerr<<"ERROR : Program not found."<<std::endl;
        return 0;
    }
}
void Viewer::wheelEvent(QWheelEvent* e)
{
    if(e->modifiers().testFlag(Qt::ShiftModifier))
    {
        double delta = e->delta();
        if(delta>0)
        {
            camera()->setZNearCoefficient(camera()->zNearCoefficient() * 1.01);
        }
        else
            camera()->setZNearCoefficient(camera()->zNearCoefficient() / 1.01);
        update();
    }
    else
        QGLViewer::wheelEvent(e);
}

bool Viewer::testDisplayId(double x, double y, double z)
{
    return d->scene->testDisplayId(x,y,z,this);
}

QPainter* Viewer::getPainter(){return d->painter;}

void Viewer::paintEvent(QPaintEvent *)
{
    paintGL();
}

void Viewer::paintGL()
{
  if (!d->painter->isActive())
    d->painter->begin(this);
  d->painter->beginNativePainting();
  glClearColor(backgroundColor().redF(), backgroundColor().greenF(), backgroundColor().blueF(), 1.0);

  //set the default frustum
  GLdouble d_mat[16];
  this->camera()->getProjectionMatrix(d_mat);
  //Convert the GLdoubles matrices in GLfloats
  for (int i=0; i<16; ++i)
      d->projectionMatrix.data()[i] = GLfloat(d_mat[i]);

  preDraw();
  draw();
  postDraw();
  d->painter->endNativePainting();
  d->painter->end();
}

void Viewer::postDraw()
{

#ifdef GL_RESCALE_NORMAL  // OpenGL 1.2 Only...
  glEnable(GL_RESCALE_NORMAL);
#endif

  if (cameraIsEdited())
    camera()->drawAllPaths();

  // Pivot point, line when camera rolls, zoom region
  drawVisualHints();

  if (gridIsDrawn()) { glLineWidth(1.0); drawGrid(camera()->sceneRadius()); }
  if (axisIsDrawn()) { glLineWidth(2.0); drawAxis(camera()->sceneRadius()); }

  // FPS computation
  const unsigned int maxCounter = 20;
  if (++d->fpsCounter == maxCounter)
  {
    d->f_p_s = 1000.0 * maxCounter / d->fpsTime.restart();
    d->fpsString = tr("%1Hz", "Frames per seconds, in Hertz").arg(d->f_p_s, 0, 'f', ((d->f_p_s < 10.0)?1:0));
    d->fpsCounter = 0;
  }

}
void Viewer::displayMessage(const QString &_message, int delay)
{
          d->message = _message;
          d->_displayMessage = true;
          // Was set to single shot in defaultConstructor.
          d->messageTimer.start(delay);
          if (textIsEnabled())
                  update();
}
void Viewer::hideMessage()
{
        d->_displayMessage = false;
        if (textIsEnabled())
                update();
}
void Viewer::printMessage(QString _message, int ms_delay)
{
  displayMessage(_message, ms_delay);
}

void Viewer_impl::showDistance(QPoint pixel)
{
    static bool isAset = false;
    bool found;
    qglviewer::Vec point;
    point = viewer->camera()->pointUnderPixel(pixel, found);
    if(!isAset && found)
    {
        //set APoint
        APoint = point;
        isAset = true;
        clearDistancedisplay();
    }
    else if (found)
    {
        //set BPoint
        BPoint = point;
        isAset = false;

        // fills the buffers
        std::vector<float> v;
        v.resize(6);
        v[0] = APoint.x; v[1] = APoint.y; v[2] = APoint.z;
        v[3] = BPoint.x; v[4] = BPoint.y; v[5] = BPoint.z;
        rendering_program_dist.bind();
        vao[1].bind();
        buffers[3].bind();
        buffers[3].allocate(v.data(),6*sizeof(float));
        rendering_program_dist.enableAttributeArray("vertex");
        rendering_program_dist.setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[3].release();
        vao[1].release();
        rendering_program_dist.release();
        distance_is_displayed = true;
        double dist = std::sqrt((BPoint.x-APoint.x)*(BPoint.x-APoint.x) + (BPoint.y-APoint.y)*(BPoint.y-APoint.y) + (BPoint.z-APoint.z)*(BPoint.z-APoint.z));
        QFont font;
        font.setBold(true);
        TextItem *ACoord = new TextItem(APoint.x, APoint.y, APoint.z,QString("A(%1,%2,%3)").arg(APoint.x-offset.x).arg(APoint.y-offset.y).arg(APoint.z-offset.z), true, font, Qt::red, true);
        distance_text.append(ACoord);
        TextItem *BCoord = new TextItem(BPoint.x, BPoint.y, BPoint.z,QString("B(%1,%2,%3)").arg(BPoint.x-offset.x).arg(BPoint.y-offset.y).arg(BPoint.z-offset.z), true, font, Qt::red, true);
        distance_text.append(BCoord);
        qglviewer::Vec centerPoint = 0.5*(BPoint+APoint);
        TextItem *centerCoord = new TextItem(centerPoint.x, centerPoint.y, centerPoint.z,QString(" distance: %1").arg(dist), true, font, Qt::red, true);

        distance_text.append(centerCoord);
        Q_FOREACH(TextItem* ti, distance_text)
          textRenderer->addText(ti);
        Q_EMIT(viewer->sendMessage(QString("First point : A(%1,%2,%3), second point : B(%4,%5,%6), distance between them : %7")
                  .arg(APoint.x-offset.x)
                  .arg(APoint.y-offset.y)
                  .arg(APoint.z-offset.z)
                  .arg(BPoint.x-offset.x)
                  .arg(BPoint.y-offset.y)
                  .arg(BPoint.z-offset.z)
                  .arg(dist)));
    }

}

void Viewer_impl::clearDistancedisplay()
{
  distance_is_displayed = false;
  Q_FOREACH(TextItem* ti, distance_text)
  {
    textRenderer->removeText(ti);
    delete ti;
  }
  distance_text.clear();
}

void Viewer_impl::setFrustum(double l, double r, double t, double b, double n, double f)
{
  double A = 2*n/(r-l);
  double B = (r+l)/(r-l);
  double C = 2*n/(t-b);
  double D = (t+b)/(t-b);
  float E = -(f+n)/(f-n);
  float F = -2*(f*n)/(f-n);
  projectionMatrix.setRow(0, QVector4D(A,0,B,0));
  projectionMatrix.setRow(1, QVector4D(0,C,D,0));
  projectionMatrix.setRow(2, QVector4D(0,0,E,F));
  projectionMatrix.setRow(3, QVector4D(0,0,-1,0));

}

#include "ui_ImageInterface.h"
class ImageInterface: public QDialog, public Ui::ImageInterface
{
  Q_OBJECT
  qreal ratio;
  QWidget *currentlyFocused;
public:
  ImageInterface(QWidget *parent, qreal ratio)
    : QDialog(parent), ratio(ratio)
  {
    currentlyFocused = NULL;
    setupUi(this);
    connect(imgHeight, SIGNAL(valueChanged(int)),
            this, SLOT(imgHeightValueChanged(int)));

    connect(imgWidth, SIGNAL(valueChanged(int)),
            this, SLOT(imgWidthValueChanged(int)));

    connect(qApp, SIGNAL(focusChanged(QWidget*, QWidget*)),
            this, SLOT(onFocusChanged(QWidget*, QWidget*)));
  }
private Q_SLOTS:
  void imgHeightValueChanged(int i)
  {
    if(currentlyFocused == imgHeight
       && ratioCheckBox->isChecked())
    {imgWidth->setValue(i*ratio);}
  }

  void imgWidthValueChanged(int i)
  {
    if(currentlyFocused == imgWidth
       && ratioCheckBox->isChecked())
    {imgHeight->setValue(i/ratio);}
  }

  void onFocusChanged(QWidget*, QWidget* now)
  {
    currentlyFocused = now;
  }
};

void Viewer::saveSnapshot(bool, bool)
{
  qreal aspectRatio = width() / static_cast<qreal>(height());
  static ImageInterface* imageInterface = NULL;

  if (!imageInterface)
    imageInterface = new ImageInterface(this, aspectRatio);

  imageInterface->imgWidth->setValue(width());
  imageInterface->imgHeight->setValue(height());
  imageInterface->imgQuality->setValue(snapshotQuality());

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
  QImage* image= d->takeSnapshot(this, imageInterface->imgQuality->value(), imageInterface->color_comboBox->currentIndex(),
        finalSize, imageInterface->oversampling->value(), expand);
  if(image)
  {
    image->save(fileName);
    delete image;
  }

}
//copy a snapshot with transparent background with arbitrary quality values.
QImage* Viewer_impl::takeSnapshot(Viewer *viewer, int quality, int background_color, QSize finalSize, double oversampling, bool expand)
{
  qreal aspectRatio = viewer->width() / static_cast<qreal>(viewer->height());
  viewer->setSnapshotQuality(quality);
  GLfloat alpha = 1.0f;
  QColor previousBGColor = viewer->backgroundColor();
  switch(background_color)
  {
  case 0:
    break;
  case 1:
    viewer->setBackgroundColor(QColor(Qt::transparent));
    alpha = 0.0f;
    break;
  case 2:
    QColor c =  QColorDialog::getColor();
    if(c.isValid()) {
      viewer->setBackgroundColor(c);
    }
    else
      return NULL;
    break;
  }


  QSize subSize(int(viewer->width()/oversampling), int(viewer->height()/oversampling));
  QSize size=QSize(viewer->width(), viewer->height());


  qreal newAspectRatio = finalSize.width() / static_cast<qreal>(finalSize.height());

  qreal zNear = viewer->camera()->zNear();
  qreal zFar = viewer->camera()->zFar();

  qreal xMin, yMin;

  if ((expand && (newAspectRatio>aspectRatio)) || (!expand && (newAspectRatio<aspectRatio)))
  {
    yMin = zNear * tan(viewer->camera()->fieldOfView() / 2.0);
    xMin = newAspectRatio * yMin;
  }
  else
  {
    xMin = zNear * tan(viewer->camera()->fieldOfView() / 2.0) * aspectRatio;
    yMin = xMin / newAspectRatio;
  }

  QImage *image = new QImage(finalSize.width(), finalSize.height(), QImage::Format_ARGB32);

  if (image->isNull())
  {
    QMessageBox::warning(viewer, "Image saving error",
                         "Unable to create resulting image",
                         QMessageBox::Ok, QMessageBox::NoButton);
    viewer->setBackgroundColor(previousBGColor);
    return NULL;
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

  QOpenGLFramebufferObject* fbo = new QOpenGLFramebufferObject(size, QOpenGLFramebufferObject::CombinedDepthStencil);
  viewer->makeCurrent();
  int count=0;
  for (int i=0; i<nbX; i++)
    for (int j=0; j<nbY; j++)
    {
      setFrustum(-xMin + i*deltaX, -xMin + (i+1)*deltaX, yMin - j*deltaY, yMin - (j+1)*deltaY, zNear, zFar);
      fbo->bind();
      viewer->glClearColor(viewer->backgroundColor().redF(), viewer->backgroundColor().greenF(), viewer->backgroundColor().blueF(), alpha);
      viewer->preDraw();
      viewer->draw();
      viewer->postDraw();
      fbo->release();

      QImage snapshot = fbo->toImage();
      QImage subImage = snapshot.scaled(subSize, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
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
      count++;
    }
  if(background_color !=0)
    viewer->setBackgroundColor(previousBGColor);
  return image;
}
void Viewer_impl::sendSnapshotToClipboard(Viewer *viewer)
{
  QImage * snap = takeSnapshot(viewer, 95, 1, 2*viewer->size(), 1, true);
  if(snap)
  {
#if defined(_WIN32)
    QApplication::clipboard()->setImage(*snap);
    QMimeData *mimeData = new QMimeData();
    QByteArray ba;
    QBuffer buffer(&ba);
    buffer.open(QIODevice::WriteOnly);
    snap->save(&buffer, "PNG"); // writes image into ba in PNG format
    buffer.close();
    mimeData->setData("PNG", ba);
    //According to the doc, the ownership of mime_data is transferred to
    //clipboard, so this is not a memory leak.
    QApplication::clipboard()->setMimeData(mimeData);
#else
    QApplication::clipboard()->setImage(*snap);
#endif
    delete snap;
  }
}
void Viewer::SetOrthoProjection(bool b)
{
  if(b)
    camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
  else
    camera()->setType(qglviewer::Camera::PERSPECTIVE);
  update();
}

void Viewer::setOffset(qglviewer::Vec offset){ d->offset = offset; }
qglviewer::Vec Viewer::offset()const { return d->offset; }
void Viewer::setSceneBoundingBox(const qglviewer::Vec &min, const qglviewer::Vec &max)
{
  QGLViewer::setSceneBoundingBox(min+d->offset, max+d->offset);
}

void Viewer::updateIds(CGAL::Three::Scene_item * item)
{
  //all ids are computed when they are displayed the first time.
  //Calling printPrimitiveIds twice hides and show the ids again, so they are re-computed.

  d->scene->updatePrimitiveIds(this, item);
  d->scene->updatePrimitiveIds(this, item);
}

TextRenderer* Viewer::textRenderer()
{
  return d->textRenderer;
}

bool Viewer::isExtensionFound()
{
  return d->extension_is_found;
}
 #include "Viewer.moc"
