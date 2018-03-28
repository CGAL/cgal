#include "Viewer.h"
#include <CGAL/Three/Scene_draw_interface.h>
#include <QMouseEvent>
#include <QKeyEvent>
#include <CGAL/Qt/manipulatedCameraFrame.h>
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
  bool clipping;
  QVector4D clipbox[6];
  QPainter *painter;
  // M e s s a g e s
  QString message;
  bool _displayMessage;
  QTimer messageTimer;
  QOpenGLFunctions_4_3_Compatibility* _recentFunctions;
  bool is_2d_selection_mode;

  //! The buffers used to draw the axis system
  QOpenGLBuffer buffer;
  //! The VAO used to draw the axis system
  QOpenGLVertexArrayObject vao;
  //! The rendering program used to draw the distance
  QOpenGLShaderProgram rendering_program_dist;
  QList<TextItem*>  distance_text;
  //! Decides if the text is displayed in the drawVisualHints function.
  bool has_text;
  //! Decides if the distance between APoint and BPoint must be drawn;
  bool distance_is_displayed;
  bool i_is_pressed;
  bool initialized;
  bool z_is_pressed;
  QImage static_image;
  //!Draws the distance between two selected points.
  void showDistance(QPoint);
  qglviewer::Vec APoint;
  qglviewer::Vec BPoint;
  bool is_d_pressed;
  bool extension_is_found;
  int quality;

  TextRenderer *textRenderer;
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
  d->initialized = false;
  d->antialiasing = antialiasing;
  d->twosides = false;
  this->setProperty("draw_two_sides", false);
  d->macro_mode = false;
  d->inFastDrawing = true;
  d->inDrawWithNames = false;
  d->clipping = false;
  d->shader_programs.resize(NB_OF_PROGRAMS);
  d->textRenderer = new TextRenderer();
  d->is_2d_selection_mode = false;
  connect( d->textRenderer, SIGNAL(sendMessage(QString,int)),
           this, SLOT(printMessage(QString,int)) );
  connect(&d->messageTimer, SIGNAL(timeout()), SLOT(hideMessage()));
  setShortcut(EXIT_VIEWER, 0);
  setKeyDescription(Qt::Key_T,
                    tr("Turn the camera by 180 degrees"));
  setKeyDescription(Qt::Key_M,
                    tr("Toggle macro mode: useful to view details very near from the camera, "
                       "but decrease the z-buffer precision"));
  setKeyDescription(Qt::Key_I + Qt::CTRL,
                      tr("Toggle the primitive IDs visibility of the selected Item."));
  setKeyDescription(Qt::Key_D,
                      tr("Disable the distance between two points  visibility."));
  setKeyDescription(Qt::Key_F5,
                    tr("Reload selected items if possible."));

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
  d->has_text = false;
  d->i_is_pressed = false;
  d->z_is_pressed = false;
  d->distance_is_displayed = false;
  d->is_d_pressed = false;
  d->viewer = this;
  setTextIsEnabled(true);
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
  glEnable(GL_DEPTH_TEST);
  d->draw_aux(false, this);
}

void Viewer::fastDraw()
{
  d->draw_aux(false, this);
}

void Viewer::init()
{
  if(!isOpenGL_4_3())
  {
    std::cerr<<"The openGL context initialization failed "
    "and the default context (2.1) will be used" <<std::endl;
  }
  else
  {
    d->_recentFunctions = new QOpenGLFunctions_4_3_Compatibility();
    d->_recentFunctions->initializeOpenGLFunctions();
  }
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
  d->vao.create();
  d->buffer.create();

  QOpenGLShader *vertex_shader, *fragment_shader;
  
  //setting the program used for the distance
     {
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
  d->painter = new QPainter();
  d->initialized = true;
}

#include <QMouseEvent>

void Viewer::mousePressEvent(QMouseEvent* event)
{
  makeCurrent();
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
    makeCurrent();
    QGLViewer::mousePressEvent(event);
  }
}
void Viewer::mouseDoubleClickEvent(QMouseEvent* event)
{
  makeCurrent(); 
  QGLViewer::mouseDoubleClickEvent(event);
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
    else if(e->key() == Qt::Key_C) {
      QVector4D box[6];
      for(int i=0; i<6; ++i)
        box[i] = QVector4D(1,0,0,0);
          enableClippingBox(box);
        }
  }
  else if(e->key() == Qt::Key_I && e->modifiers() & Qt::ControlModifier){
    d->scene->printAllIds(this);
    update();
    return;
  }

  else if(e->key() == Qt::Key_C && e->modifiers() & Qt::ControlModifier){
    d->sendSnapshotToClipboard(this);
    return;
  }

  else if(e->key() == Qt::Key_S && e->modifiers() & Qt::ControlModifier){
    this->saveSnapshot();
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
  qglviewer::Vec point = camera()->pointUnderPixel(pixel, found) - offset();
  if(found) {
    Q_EMIT selectedPoint(point.x,
                       point.y,
                       point.z);
    const qglviewer::Vec orig = camera()->position() - offset();
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
    program->setUniformValue("is_clipbox_on", d->clipping);
    if(d->clipping)
    {
      QMatrix4x4 clipbox1;
      QMatrix4x4 clipbox2;
      for(int i=0;i<12;++i)
      {
        clipbox1.data()[i]=d->clipbox[i/4][i%4];
        clipbox2.data()[i]=d->clipbox[(i+12)/4][(i+12)%4];
      }
      program->setUniformValue("clipbox1", clipbox1);
      program->setUniformValue("clipbox2", clipbox2);
    }
    switch(program_name)
    {
    case PROGRAM_WITH_LIGHT:
    case PROGRAM_C3T3:
    case PROGRAM_PLANE_TWO_FACES:
    case PROGRAM_INSTANCED:
    case PROGRAM_WITH_TEXTURE:
    case PROGRAM_CUTPLANE_SPHERES:
    case PROGRAM_SPHERES:
    case PROGRAM_OLD_FLAT:
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
    case PROGRAM_OLD_FLAT:
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
  QGLViewer::beginSelection(point);
  d->scene->setPickedPixel(point);
}
void Viewer::endSelection(const QPoint& point)
{
  QGLViewer::endSelection(point);
    //redraw the true scene for the glReadPixel in postSelection();
    d->draw_aux(false, this);
}

void Viewer::drawVisualHints()
{

    QGLViewer::drawVisualHints();

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
        d->vao.bind();
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(2));
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(2));
        d->vao.release();
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
    //Prints the displayMessage
    QFont font = QFont();
    QFontMetrics fm(font);
    TextItem *message_text = new TextItem(10 + fm.width(d->message)/2, height()-20, 0, d->message, false, QFont(), Qt::gray );
    if (d->_displayMessage)
    {
      d->textRenderer->addText(message_text);
    }
    d->textRenderer->draw(this);
    
    if (d->_displayMessage)
      d->textRenderer->removeText(message_text);
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
    if(strcmp(f_shader,":/cgal/Polyhedron_3/resources/shader_flat.f" ) == 0)
    {
      if(!program->addShaderFromSourceFile(QOpenGLShader::Geometry,":/cgal/Polyhedron_3/resources/shader_flat.g" ))
      {
        std::cerr<<"adding geometry shader FAILED"<<std::endl;
      }
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
    case PROGRAM_SPHERES:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_spheres.v" , ":/cgal/Polyhedron_3/resources/shader_with_light.f");
      break;
    case PROGRAM_FLAT:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_flat.v", ":/cgal/Polyhedron_3/resources/shader_flat.f");
    case PROGRAM_OLD_FLAT:
      return declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_light.v", ":/cgal/Polyhedron_3/resources/shader_old_flat.f");
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
  if(!d->initialized)
    initializeGL();
  paintGL();
}

void Viewer::paintGL()
{
  makeCurrent();
  if (!d->painter->isActive())
    d->painter->begin(this);
  if(d->is_2d_selection_mode)
  {
    d->painter->drawImage(QPoint(0,0), d->static_image);
  }
  else
  {
    d->painter->beginNativePainting();
    glClearColor(backgroundColor().redF(), backgroundColor().greenF(), backgroundColor().blueF(), 1.0);
    glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
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
  }
  d->painter->end();
  doneCurrent();
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
        vao.bind();
        buffer.bind();
        buffer.allocate(v.data(),6*sizeof(float));
        rendering_program_dist.enableAttributeArray("vertex");
        rendering_program_dist.setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffer.release();
        vao.release();
        rendering_program_dist.release();
        distance_is_displayed = true;
        double dist = std::sqrt((BPoint.x-APoint.x)*(BPoint.x-APoint.x) + (BPoint.y-APoint.y)*(BPoint.y-APoint.y) + (BPoint.z-APoint.z)*(BPoint.z-APoint.z));
        QFont font;
        font.setBold(true);
        TextItem *ACoord = new TextItem(APoint.x, APoint.y, APoint.z,QString("A(%1,%2,%3)").arg(APoint.x-viewer->offset().x).arg(APoint.y-viewer->offset().y).arg(APoint.z-viewer->offset().z), true, font, Qt::red, true);
        distance_text.append(ACoord);
        TextItem *BCoord = new TextItem(BPoint.x, BPoint.y, BPoint.z,QString("B(%1,%2,%3)").arg(BPoint.x-viewer->offset().x).arg(BPoint.y-viewer->offset().y).arg(BPoint.z-viewer->offset().z), true, font, Qt::red, true);
        distance_text.append(BCoord);
        qglviewer::Vec centerPoint = 0.5*(BPoint+APoint);
        TextItem *centerCoord = new TextItem(centerPoint.x, centerPoint.y, centerPoint.z,QString(" distance: %1").arg(dist), true, font, Qt::red, true);

        distance_text.append(centerCoord);
        Q_FOREACH(TextItem* ti, distance_text)
          textRenderer->addText(ti);
        Q_EMIT(viewer->sendMessage(QString("First point : A(%1,%2,%3), second point : B(%4,%5,%6), distance between them : %7")
                  .arg(APoint.x-viewer->offset().x)
                  .arg(APoint.y-viewer->offset().y)
                  .arg(APoint.z-viewer->offset().z)
                  .arg(BPoint.x-viewer->offset().x)
                  .arg(BPoint.y-viewer->offset().y)
                  .arg(BPoint.z-viewer->offset().z)
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

void Viewer::saveSnapshot()
{
  qreal aspectRatio = width() / static_cast<qreal>(height());
  static ImageInterface* imageInterface = NULL;

  if (!imageInterface)
    imageInterface = new ImageInterface(this, aspectRatio);

  imageInterface->imgWidth->setValue(width());
  imageInterface->imgHeight->setValue(height());
  imageInterface->imgQuality->setValue(d->quality);

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
  viewer->makeCurrent();
  this->quality = quality;
  qreal aspectRatio = viewer->width() / static_cast<qreal>(viewer->height());
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

  QOpenGLFramebufferObject fbo(size, QOpenGLFramebufferObject::CombinedDepthStencil);
  for (int i=0; i<nbX; i++)
    for (int j=0; j<nbY; j++)
    {
      setFrustum(-xMin + i*deltaX, -xMin + (i+1)*deltaX, yMin - j*deltaY, yMin - (j+1)*deltaY, zNear, zFar);
      fbo.bind();
      viewer->glClearColor(viewer->backgroundColor().redF(), viewer->backgroundColor().greenF(), viewer->backgroundColor().blueF(), alpha);
      viewer->preDraw();
      viewer->draw();
      fbo.release();

      QImage snapshot = fbo.toImage();
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

void Viewer::disableClippingBox()
{
  d->clipping = false;
}

void Viewer::enableClippingBox(QVector4D box[6])
{
  d->clipping = true;
  for(int i=0; i<6; ++i)
    d->clipbox[i] = box[i];
}

QOpenGLFunctions_4_3_Compatibility* Viewer::openGL_4_3_functions() { return d->_recentFunctions; }

void Viewer::set2DSelectionMode(bool b) { d->is_2d_selection_mode = b; }

void Viewer::setStaticImage(QImage image) { d->static_image = image; }

const QImage& Viewer:: staticImage() const { return d->static_image; }



#include "Viewer.moc"
