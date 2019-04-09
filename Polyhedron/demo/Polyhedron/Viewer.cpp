#include "Viewer.h"
#include <CGAL/Three/Scene_draw_interface.h>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QSettings>
#include <QDebug>
#include <QSettings>
#include <QOpenGLShader>
#include <QFileDialog>
#include <QOpenGLShaderProgram>
#include <QOpenGLFramebufferObject>
#include <QMessageBox>
#include <QColorDialog>
#include <QInputDialog>
#include <cmath>
#include <QApplication>
#include <QOpenGLDebugLogger>
#include <QStyleFactory>

#include "ui_LightingDialog.h"

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
  bool projection_is_ortho;
  GLfloat gl_point_size;
  QVector4D clipbox[6];
  QPainter *painter;

  // L i g h t i n g
  QVector4D position;
  QVector4D ambient;
  QVector4D diffuse;
  QVector4D specular;
  float spec_power;
  
  // M e s s a g e s
  QString message;
  bool _displayMessage;
  QTimer messageTimer;
  QOpenGLFunctions_4_3_Core* _recentFunctions;
  bool is_2d_selection_mode;

  // D e p t h  P e e l i n g
  // \param pass the current pass in the Depth Peeling (transparency) algorithm.
  // -1 means that no depth peeling is applied.
  // \param writing_depth means that the color of the faces will be drawn in a grayscale
  // according to the depth of the fragment in the shader. It is used by the transparency.
  // \param fbo contains the texture used by the Depth Peeling algorithm.
  // Should be NULL if pass <= 0;
  int current_pass;
  bool writing_depth;
  int total_pass;
  int current_total_pass;
  QOpenGLFramebufferObject* dp_fbo;
  QOpenGLDebugLogger *logger;


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
  CGAL::qglviewer::Vec APoint;
  CGAL::qglviewer::Vec BPoint;
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
  void sendSnapshotToClipboard(Viewer*);
};

class LightingDialog :
    public QDialog,
    public Ui::LightingDialog
{
  Q_OBJECT
public: 
  QColor ambient, diffuse, specular;
  LightingDialog(Viewer_impl* d)
  {
    setupUi(this);
    position_lineEdit->setText(QString("%1,%2,%3")
                               .arg(d->position.x())
                               .arg(d->position.y())
                               .arg(d->position.z()));    
    QPalette palette;
    ambient=QColor(255*d->ambient.x(), 
                   255*d->ambient.y(), 
                   255*d->ambient.z());
    palette.setColor(QPalette::Button,ambient);
    ambientButton->setPalette(palette);
    ambientButton->setStyle(QStyleFactory::create("Fusion"));
    
    diffuse=QColor(255*d->diffuse.x(), 
                   255*d->diffuse.y(), 
                   255*d->diffuse.z());
    palette.setColor(QPalette::Button,diffuse);
    diffuseButton->setPalette(palette);
    diffuseButton->setStyle(QStyleFactory::create("Fusion"));
    
    specular=QColor(255*d->specular.x(), 
                    255*d->specular.y(), 
                    255*d->specular.z());
    palette.setColor(QPalette::Button,specular);
    specularButton->setPalette(palette);
    specularButton->setStyle(QStyleFactory::create("Fusion"));
    spec_powrSlider->setValue(static_cast<int>(d->spec_power));
    
    connect(&ambient_dial, &QColorDialog::currentColorChanged, this, &LightingDialog::ambient_changed );
    connect(&diffuse_dial, &QColorDialog::currentColorChanged, this, &LightingDialog::diffuse_changed );
    connect(&spec_dial, &QColorDialog::currentColorChanged, this,    &LightingDialog::specular_changed);
    
    connect(ambientButton, &QPushButton::clicked, 
            [this](){
      ambient_dial.setCurrentColor(ambient);
      ambient_dial.exec();
      ambient = ambient_dial.selectedColor();
      QPalette palette;
      palette.setColor(QPalette::Button, ambient);
      ambientButton->setPalette(palette);
    });
    connect(diffuseButton, &QPushButton::clicked, 
            [this](){
      diffuse_dial.setCurrentColor(diffuse);
      diffuse_dial.exec();
      diffuse = diffuse_dial.selectedColor();
      QPalette palette;
      palette.setColor(QPalette::Button, diffuse);
      diffuseButton->setPalette(palette);
    });
    connect(specularButton, &QPushButton::clicked, 
            [this](){
      spec_dial.setCurrentColor(specular);
      spec_dial.exec();
      specular = spec_dial.selectedColor();
      QPalette palette;
      palette.setColor(QPalette::Button, specular);
      specularButton->setPalette(palette);
    });
    
    //D e f a u l t - S e t t i n g s
    connect(buttonBox->button(QDialogButtonBox::StandardButton::RestoreDefaults), &QPushButton::clicked,
            [this](){
      position_lineEdit->setText(QString("0,0,1"));
      ambient=QColor(77,77,77);
      diffuse=QColor(204,204,204);
      specular=QColor(0,0,0);
      spec_powrSlider->setValue(51);
      QPalette palette;
      palette.setColor(QPalette::Button, ambient);
      ambientButton->setPalette(palette);
      palette.setColor(QPalette::Button, diffuse);
      diffuseButton->setPalette(palette);
      palette.setColor(QPalette::Button, specular);
      specularButton->setPalette(palette);
    });
  }
private Q_SLOTS:
  void diffuse_changed()
  {
    diffuse = diffuse_dial.currentColor();
    s_diffuse_changed(); 
  }
  void ambient_changed()
  {
    ambient = ambient_dial.currentColor();
    s_ambient_changed(); 
  }
  void specular_changed()
  {
    specular = spec_dial.currentColor();
    s_specular_changed();
  }
Q_SIGNALS:
  void s_diffuse_changed();
  void s_ambient_changed();
  void s_specular_changed();
private:
  QColorDialog diffuse_dial;
  QColorDialog ambient_dial;
  QColorDialog spec_dial;
};

Viewer::Viewer(QWidget* parent, bool antialiasing)
  : CGAL::Three::Viewer_interface(parent)
{
  d = new Viewer_impl;
  QSettings viewer_settings;
  // enable anti-aliasing
  QString cam_pos = viewer_settings.value("cam_pos", QString("0.0,0.0,1.0")).toString();
  d->position = QVector4D(cam_pos.split(",").at(0).toFloat(),
                          cam_pos.split(",").at(1).toFloat(),
                          cam_pos.split(",").at(2).toFloat(),
                          1.0f);
  
  QString ambient = viewer_settings.value("ambient", QString("0.4,0.4,0.4")).toString();
  d->ambient = QVector4D(ambient.split(",").at(0).toFloat(),
                         ambient.split(",").at(1).toFloat(),
                         ambient.split(",").at(2).toFloat(),
                         1.0f);
  
  QString diffuse = viewer_settings.value("diffuse", QString("1.0,1.0,1.0")).toString();
  d->diffuse = QVector4D(diffuse.split(",").at(0).toFloat(),
                         diffuse.split(",").at(1).toFloat(),
                         diffuse.split(",").at(2).toFloat(),
                         1.0f);
  
  QString specular = viewer_settings.value("specular", QString("0.0,0.0,0.0")).toString();
  d->specular = QVector4D(specular.split(",").at(0).toFloat(),
                          specular.split(",").at(1).toFloat(),
                          specular.split(",").at(2).toFloat(),
                          1.0f);
  
  d->spec_power = viewer_settings.value("spec_power", 51.8).toFloat();
  d->scene = 0;
  d->projection_is_ortho = false;
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
  setShortcut(CGAL::qglviewer::EXIT_VIEWER, 0);
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

  //modify mouse bindings that have been updated
  setMouseBinding(Qt::Key(0), Qt::NoModifier, Qt::LeftButton, CGAL::qglviewer::RAP_FROM_PIXEL, true, Qt::RightButton);
  setMouseBindingDescription(Qt::ShiftModifier, Qt::RightButton,
                             tr("Select and pop context menu"));
  setMouseBinding(Qt::Key_R, Qt::NoModifier, Qt::LeftButton, CGAL::qglviewer::RAP_FROM_PIXEL);
  //use the new API for these
  setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, CGAL::qglviewer::SELECT);

  setMouseBindingDescription(Qt::Key(0), Qt::ShiftModifier, Qt::LeftButton,
                             tr("Selects and display context "
                                "menu of the selected item"));
  setMouseBindingDescription(Qt::Key_I, Qt::NoModifier, Qt::LeftButton,
                             tr("Show/hide the primitive ID."));
  setMouseBindingDescription(Qt::Key_D, Qt::NoModifier, Qt::LeftButton,
                             tr("Selects a point. When the second point is selected,  "
                                "displays the two points and the distance between them."));
  setMouseBindingDescription(Qt::Key_O, Qt::NoModifier, Qt::LeftButton,
                             tr("Move the camera orthogonally to the picked facet of a Scene_surface_mesh_item or "
                                "to the current selection of a Scene_points_with_normal_item."));
  setKeyDescription(Qt::Key_F5,
                    tr("Reloads the selected item if possible."));
  setKeyDescription(Qt::Key_F11,
                    tr("Toggle the viewer's fullscreen mode."));

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
    QSettings viewer_settings;
    viewer_settings.setValue("cam_pos",
                             QString("%1,%2,%3")
                             .arg(d->position.x())
                             .arg(d->position.y())
                             .arg(d->position.z()));
    viewer_settings.setValue("ambient",
                             QString("%1,%2,%3")
                             .arg(d->ambient.x())
                             .arg(d->ambient.y())
                             .arg(d->ambient.z()));
    viewer_settings.setValue("diffuse",
                             QString("%1,%2,%3")
                             .arg(d->diffuse.x())
                             .arg(d->diffuse.y())
                             .arg(d->diffuse.z()));
    viewer_settings.setValue("specular",
                             QString("%1,%2,%3")
                             .arg(d->specular.x())
                             .arg(d->specular.y())
                             .arg(d->specular.z()));
    viewer_settings.setValue("spec_power",
                             d->spec_power);
    if(d->_recentFunctions)
      delete d->_recentFunctions;
    if(d->painter)
      delete d->painter;
    if(d->textRenderer)
      d->textRenderer->deleteLater();
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
    "and the default context (2.0 ES) will be used. \n"
    " This means, among other things, that no widelines can be displayed,"
    " which makes selected edges harder to see." <<std::endl;
  }
  else
  {
    d->_recentFunctions = new QOpenGLFunctions_4_3_Core();
    d->_recentFunctions->initializeOpenGLFunctions();
  }
  d->logger = new QOpenGLDebugLogger(this);
  if(!d->logger->initialize())
    qDebug()<<"logger could not init.";
  else{
    connect(d->logger, SIGNAL(messageLogged(QOpenGLDebugMessage)), this, SLOT(messageLogged(QOpenGLDebugMessage)));
    d->logger->startLogging();
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
  QSettings settings;
  QString colorname = settings.value("background_color", "#ffffff").toString();
  QColor bc(colorname);
  setBackgroundColor(bc);
  d->vao.create();
  d->buffer.create();
  
  //setting the program used for the distance
     {
         //Vertex source code
         const char vertex_source_dist[] =
         {
             "#version 150  \n"
             "in vec4 vertex;\n"
             "uniform mat4 mvp_matrix;\n"
             "uniform float point_size;\n"
             "void main(void)\n"
             "{\n"
             "   gl_PointSize = point_size; \n"
             "   gl_Position = mvp_matrix * vertex; \n"
             "} \n"
             "\n"
         };
         const char vertex_source_comp_dist[] =
         {
             "attribute highp vec4 vertex;\n"
             "uniform highp mat4 mvp_matrix;\n"
             "uniform highp float point_size;\n"
             "void main(void)\n"
             "{\n"
             "   gl_PointSize = point_size; \n"
             "   gl_Position = mvp_matrix * vertex; \n"
             "} \n"
             "\n"
         };
         //Fragment source code
         const char fragment_source_dist[] =
         {
             "#version 150  \n"
             "out vec4 out_color; \n"
             "void main(void) { \n"
             "out_color = vec4(0.0,0.0,0.0,1.0); \n"
             "} \n"
             "\n"
         };
         const char fragment_source_comp_dist[] =
         {
             "void main(void) { \n"
             "gl_FragColor = vec4(0.0,0.0,0.0,1.0); \n"
             "} \n"
             "\n"
         };
         QOpenGLShader vertex_shader(QOpenGLShader::Vertex);
         QOpenGLShader fragment_shader(QOpenGLShader::Fragment);
         if(isOpenGL_4_3())
         {
           if(!vertex_shader.compileSourceCode(vertex_source_dist))
           {
             std::cerr<<"Compiling vertex source FAILED"<<std::endl;
           }
           
           if(!fragment_shader.compileSourceCode(fragment_source_dist))
           {
             std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
           }
         }
         else
         {
           if(!vertex_shader.compileSourceCode(vertex_source_comp_dist))
           {
             std::cerr<<"Compiling vertex source FAILED"<<std::endl;
           }
           
           if(!fragment_shader.compileSourceCode(fragment_source_comp_dist))
           {
             std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
           }
         }
         if(!d->rendering_program_dist.addShader(&vertex_shader))
         {
             std::cerr<<"adding vertex shader FAILED"<<std::endl;
         }
         if(!d->rendering_program_dist.addShader(&fragment_shader))
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
    CGAL::QGLViewer::mousePressEvent(event);
  }
}
void Viewer::mouseDoubleClickEvent(QMouseEvent* event)
{
  makeCurrent(); 
  CGAL::QGLViewer::mouseDoubleClickEvent(event);
}

#include <QContextMenuEvent>
void Viewer::contextMenuEvent(QContextMenuEvent* event)
{
  if(event->reason() != QContextMenuEvent::Mouse) {
    requestContextMenu(event->globalPos());
    event->accept();
  }
  else {
    CGAL::QGLViewer::contextMenuEvent(event);
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
    CGAL::QGLViewer::keyPressEvent(e);
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
  CGAL::QGLViewer::keyReleaseEvent(e);
}

void Viewer::turnCameraBy180Degres() {
  CGAL::qglviewer::Camera* camera = this->camera();
  using CGAL::qglviewer::ManipulatedCameraFrame;

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
  current_total_pass = viewer->inFastDrawing() ? total_pass/2 : total_pass;
  viewer->setGlPointSize(2.f);
  viewer->glEnable(GL_POLYGON_OFFSET_FILL);
  viewer->glPolygonOffset(1.0f,1.0f);

  if(!with_names && antialiasing)
  {
    viewer->glEnable(GL_BLEND);
    viewer->glEnable(GL_LINE_SMOOTH);
    viewer->glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    viewer->glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    viewer->glDisable(GL_BLEND);
    viewer->glDisable(GL_LINE_SMOOTH);
    viewer->glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    //viewer->glBlendFunc(GL_ONE, GL_ZERO);
    viewer->glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  }
  inDrawWithNames = with_names;
  if(with_names)
    scene->drawWithNames(viewer);
  else
    scene->draw(viewer);
  viewer->glDisable(GL_POLYGON_OFFSET_FILL);
}

bool Viewer::inDrawWithNames() const {
  return d->inDrawWithNames;
}

void Viewer::drawWithNames()
{
  CGAL::QGLViewer::draw();
  d->draw_aux(true, this);
}

void Viewer::postSelection(const QPoint& pixel)
{
  Q_EMIT selected(this->selectedName());
  CGAL::qglviewer::Vec point;
  bool found = true;
  if(property("picked_point").isValid()) {
    if(!property("picked_point").toList().isEmpty())
    {
      QList<QVariant> picked_point = property("picked_point").toList();
      point = CGAL::qglviewer::Vec (picked_point[0].toDouble(),
          picked_point[1].toDouble(),
          picked_point[2].toDouble());
    }
    else{
      found = false;
    }
  }
  else{
    point = camera()->pointUnderPixel(pixel, found) - offset();
  }
  if(found) {
    Q_EMIT selectedPoint(point.x,
                         point.y,
                         point.z);
    CGAL::qglviewer::Vec dir;
    CGAL::qglviewer::Vec orig;
    if(d->projection_is_ortho)
    {
      dir = camera()->viewDirection();
      orig = point;
    }
    else{
      orig = camera()->position() - offset();
      dir = point - orig;
    }
    this->setProperty("performing_selection", true);
    Q_EMIT selectionRay(orig.x, orig.y, orig.z,
                        dir.x, dir.y, dir.z);
    this->setProperty("performing_selection", false);
  }
}
bool CGAL::Three::Viewer_interface::readFrame(QString s, CGAL::qglviewer::Frame& frame)
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
  frame.setPosition(CGAL::qglviewer::Vec(vec[0],
                                   vec[1],
                                   vec[2]));
  frame.setOrientation(orient[0],
                       orient[1],
                       orient[2],
                       orient[3]);
  return true;
}

QString CGAL::Three::Viewer_interface::dumpFrame(const CGAL::qglviewer::Frame& frame) {
  const CGAL::qglviewer::Vec pos = frame.position();
  const CGAL::qglviewer::Quaternion q = frame.orientation();

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
  CGAL::qglviewer::Frame new_frame;
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
    this->camera()->getModelViewProjectionMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mvp_mat.data()[i] = GLfloat(d_mat[i]);
   
    QOpenGLShaderProgram* program = getShaderProgram(program_name);
    program->bind();
    program->setUniformValue("point_size", getGlPointSize());
    program->setUniformValue("mvp_matrix", mvp_mat);
    QMatrix4x4 id_mat;
    id_mat.setToIdentity();
    program->setUniformValue("f_matrix", id_mat);
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
    QVector4D light_pos(d->position.x(),
                        d->position.y(),
                        d->position.z(),
                        1.0f);
    switch(program_name)
    {
    case PROGRAM_WITH_LIGHT:
    case PROGRAM_SPHERES:
    case PROGRAM_CUTPLANE_SPHERES:
    case PROGRAM_HEAT_INTENSITY:
      
      program->setUniformValue("alpha", 1.0f); //overriden in item draw() if necessary
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
    case PROGRAM_HEAT_INTENSITY:
        program->setUniformValue("light_pos", light_pos);
        program->setUniformValue("light_diff",d->diffuse);
        program->setUniformValue("light_spec", d->specular);
        program->setUniformValue("light_amb", d->ambient);
        program->setUniformValue("spec_power", d->spec_power);
        program->setUniformValue("is_two_side", d->twosides);
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
    case PROGRAM_HEAT_INTENSITY:
      program->setUniformValue("mv_matrix", mv_mat);
      break;
    case PROGRAM_WITHOUT_LIGHT:
    case PROGRAM_SOLID_WIREFRAME:
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
  CGAL::QGLViewer::beginSelection(point);
  d->scene->setPickedPixel(point);
}
void Viewer::endSelection(const QPoint& point)
{
  CGAL::QGLViewer::endSelection(point);
    //redraw the true scene for the glReadPixel in postSelection();
    d->draw_aux(false, this);
}

void Viewer::drawVisualHints()
{

    CGAL::QGLViewer::drawVisualHints();

    if(d->distance_is_displayed)
    {
        glDisable(GL_DEPTH_TEST);
        QMatrix4x4 mvpMatrix;
        double mat[16];
        camera()->getModelViewProjectionMatrix(mat);
        for(int i=0; i < 16; i++)
        {
          mvpMatrix.data()[i] = (float)mat[i];
        }
        if(!isOpenGL_4_3())
        {
          //draws the distance
          //nullifies the translation
          d->rendering_program_dist.bind();
          d->rendering_program_dist.setUniformValue("mvp_matrix", mvpMatrix);
          d->rendering_program_dist.setUniformValue("point_size", GLfloat(6.0f));
          d->vao.bind();
          glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(2));
          glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(2));
          d->vao.release();
          d->rendering_program_dist.release();
          glEnable(GL_DEPTH_TEST);
        }
        else
        {          
          QOpenGLShaderProgram* program = getShaderProgram(PROGRAM_SOLID_WIREFRAME);
          program->bind();
          QVector2D vp(width(), height());
          program->setUniformValue("viewport", vp);
          program->setUniformValue("near",(GLfloat)camera()->zNear());
          program->setUniformValue("far",(GLfloat)camera()->zFar());
          program->setUniformValue("width", GLfloat(3.0f));
          program->setAttributeValue("colors", QColor(Qt::black));
          program->setUniformValue("mvp_matrix", mvpMatrix);
          QMatrix4x4 f_mat;
          f_mat.setToIdentity();
          program->setUniformValue("f_matrix", f_mat);
          d->vao.bind();
          glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(2));
          d->vao.release();
          program->release();
          
          program = getShaderProgram(PROGRAM_NO_SELECTION);
          program->bind();
          program->setAttributeValue("colors", QColor(Qt::black));
          program->setAttributeValue("point_size", 6.0f);
          program->setUniformValue("mvp_matrix", mvpMatrix);
          program->setUniformValue("f_matrix", f_mat);
          d->vao.bind();
          glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(2));
          d->vao.release();
          program->release();
        }

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
    TextItem *message_text = new TextItem(float(10 + fm.width(d->message)/2),
                                          float(height()-20),
                                          0, d->message, false,
                                          QFont(), Qt::gray );
    if (d->_displayMessage)
    {
      d->textRenderer->addText(message_text);
    }
    d->textRenderer->draw(this);
    
    if (d->_displayMessage)
      d->textRenderer->removeText(message_text);
    delete message_text;
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
    if(isOpenGL_4_3())
    {
      if(strcmp(f_shader,":/cgal/Polyhedron_3/resources/shader_flat.frag" ) == 0)
      {
        if(!program->addShaderFromSourceFile(QOpenGLShader::Geometry,":/cgal/Polyhedron_3/resources/shader_flat.geom" ))
        {
          std::cerr<<"adding geometry shader FAILED"<<std::endl;
        }
      }
      if(strcmp(f_shader,":/cgal/Polyhedron_3/resources/solid_wireframe_shader.frag" ) == 0)
      {
        if(!program->addShaderFromSourceFile(QOpenGLShader::Geometry,":/cgal/Polyhedron_3/resources/solid_wireframe_shader.geom" ))
        {
          std::cerr<<"adding geometry shader FAILED"<<std::endl;
        }
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
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3.vert" , ":/cgal/Polyhedron_3/resources/shader_c3t3.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_c3t3.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_c3t3.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasCutPlane", true);
    program->setProperty("hasTransparency", true);
    return program;
  }
  case PROGRAM_C3T3_EDGES:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3_edges.vert" , ":/cgal/Polyhedron_3/resources/shader_c3t3_edges.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_c3t3_edges.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_c3t3_edges.frag");
    program->setProperty("hasCutPlane", true);
    return program;
  }
  case PROGRAM_WITH_LIGHT:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_light.vert" , ":/cgal/Polyhedron_3/resources/shader_with_light.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_light.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_light.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasTransparency", true);
    program->setProperty("hasFMatrix", true);
    return program;
  }
  case PROGRAM_HEAT_INTENSITY:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/heat_intensity_shader.vert" , ":/cgal/Polyhedron_3/resources/heat_intensity_shader.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/heat_intensity_shader.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/heat_intensity_shader.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasTransparency", true);
    program->setProperty("hasDistanceValues", true);
    return program;
  }
  case PROGRAM_WITHOUT_LIGHT:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_without_light.vert" , ":/cgal/Polyhedron_3/resources/shader_without_light.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_without_light.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_without_light.frag");
    program->setProperty("hasFMatrix", true);
    return program;
  }
  case PROGRAM_NO_SELECTION:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_without_light.vert" , ":/cgal/Polyhedron_3/resources/shader_no_light_no_selection.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_without_light.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_no_light_no_selection.frag");
    program->setProperty("hasFMatrix", true);
    return program;
  }
  case PROGRAM_WITH_TEXTURE:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_texture.vert" , ":/cgal/Polyhedron_3/resources/shader_with_texture.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_texture.vert" ,
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_texture.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasFMatrix", true);
    program->setProperty("hasTexture", true);
    return program;
  }
  case PROGRAM_PLANE_TWO_FACES:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ?declare_program(name, ":/cgal/Polyhedron_3/resources/shader_without_light.vert" , ":/cgal/Polyhedron_3/resources/shader_plane_two_faces.frag")
       : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_without_light.vert" ,
                         ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_plane_two_faces.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    return program;
  }
  case PROGRAM_WITH_TEXTURED_EDGES:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_textured_edges.vert" , ":/cgal/Polyhedron_3/resources/shader_with_textured_edges.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_textured_edges.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_textured_edges.frag");
    program->setProperty("hasFMatrix", true);
    program->setProperty("hasTexture", true);
    return program;
  }
  case PROGRAM_INSTANCED:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_instanced.vert" , ":/cgal/Polyhedron_3/resources/shader_with_light.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_instanced.vert" ,
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_light.frag");
    
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("isInstanced", true);
    return program;
  }
  case PROGRAM_INSTANCED_WIRE:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_instanced.vert" , ":/cgal/Polyhedron_3/resources/shader_without_light.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_instanced.vert" ,
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_without_light.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasBarycenter", true);
    program->setProperty("isInstanced", true);
    return program;
  }
  case PROGRAM_CUTPLANE_SPHERES:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_c3t3_spheres.vert" , ":/cgal/Polyhedron_3/resources/shader_c3t3.frag")
        : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_c3t3_spheres.vert" , 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_c3t3.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasBarycenter", true);
    program->setProperty("hasRadius", true);
    program->setProperty("isInstanced", true);
    return program;
  }
  case PROGRAM_SPHERES:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ?declare_program(name, ":/cgal/Polyhedron_3/resources/shader_spheres.vert" , ":/cgal/Polyhedron_3/resources/shader_with_light.frag")
       : declare_program(name, ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_spheres.vert" ,
                         ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_light.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    program->setProperty("hasBarycenter", true);
    program->setProperty("hasRadius", true);
    program->setProperty("hasTransparency", true);
    program->setProperty("isInstanced", true);
    return program;
  }
  case PROGRAM_FLAT:
  {
    if(!isOpenGL_4_3())
    {
      std::cerr<<"An OpenGL context of version 4.3 is required for the program ("<<name<<")."<<std::endl;
      return 0;
    }
    QOpenGLShaderProgram* program = declare_program(name, ":/cgal/Polyhedron_3/resources/shader_flat.vert", ":/cgal/Polyhedron_3/resources/shader_flat.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    return program;
  }
  case PROGRAM_OLD_FLAT:
  {
    QOpenGLShaderProgram* program = isOpenGL_4_3() 
        ? declare_program(name, ":/cgal/Polyhedron_3/resources/shader_with_light.vert", ":/cgal/Polyhedron_3/resources/shader_old_flat.frag")
        : declare_program(name, 
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_with_light.vert",
                          ":/cgal/Polyhedron_3/resources/compatibility_shaders/shader_old_flat.frag");
    program->setProperty("hasLight", true);
    program->setProperty("hasNormals", true);
    return program;
  }
  case PROGRAM_SOLID_WIREFRAME:
    if(!isOpenGL_4_3())
    {
      std::cerr<<"An OpenGL context of version 4.3 is required for the program ("<<name<<")."<<std::endl;
      return 0;
    }
    return declare_program(name,
                           ":/cgal/Polyhedron_3/resources/solid_wireframe_shader.vert", 
                           ":/cgal/Polyhedron_3/resources/solid_wireframe_shader.frag");
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
        CGAL::QGLViewer::wheelEvent(e);
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
    glClearColor(GLfloat(backgroundColor().redF()),
                 GLfloat(backgroundColor().greenF()),
                 GLfloat(backgroundColor().blueF()),
                 1.f);
    glClearDepthf(1.0f);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //set the default frustum
    if(d->projection_is_ortho)
      camera()->setType(CGAL::qglviewer::Camera::ORTHOGRAPHIC);
    else
      camera()->setType(CGAL::qglviewer::Camera::PERSPECTIVE);
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
    CGAL::qglviewer::Vec point;
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
        v[0] = float(APoint.x); v[1] = float(APoint.y); v[2] = float(APoint.z);
        v[3] = float(BPoint.x); v[4] = float(BPoint.y); v[5] = float(BPoint.z);
       
        vao.bind();
        buffer.bind();
        buffer.allocate(v.data(),6*sizeof(float));
        rendering_program_dist.enableAttributeArray("vertex");
        rendering_program_dist.setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffer.release();
        vao.release();
        
        distance_is_displayed = true;
        double dist = std::sqrt((BPoint.x-APoint.x)*(BPoint.x-APoint.x) + (BPoint.y-APoint.y)*(BPoint.y-APoint.y) + (BPoint.z-APoint.z)*(BPoint.z-APoint.z));
        QFont font;
        font.setBold(true);
        TextItem *ACoord = new TextItem(float(APoint.x),
                                        float(APoint.y),
                                        float(APoint.z),
                                        QString("A(%1,%2,%3)").arg(APoint.x-viewer->offset().x).arg(APoint.y-viewer->offset().y).arg(APoint.z-viewer->offset().z), true, font, Qt::red, true);
        distance_text.append(ACoord);
        TextItem *BCoord = new TextItem(float(BPoint.x),
                                        float(BPoint.y),
                                        float(BPoint.z),
                                        QString("B(%1,%2,%3)").arg(BPoint.x-viewer->offset().x).arg(BPoint.y-viewer->offset().y).arg(BPoint.z-viewer->offset().z), true, font, Qt::red, true);
        distance_text.append(BCoord);
        CGAL::qglviewer::Vec centerPoint = 0.5*(BPoint+APoint);
        TextItem *centerCoord = new TextItem(float(centerPoint.x),
                                             float(centerPoint.y),
                                             float(centerPoint.z),
                                             QString(" distance: %1").arg(dist), true, font, Qt::red, true);

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

void Viewer_impl::sendSnapshotToClipboard(Viewer *viewer)
{
  QImage * snap = viewer->takeSnapshot(CGAL::qglviewer::TRANSPARENT_BACKGROUND, 2*viewer->size(), 1, true);
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
  d->projection_is_ortho = b;
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

QOpenGLFunctions_4_3_Core *Viewer::openGL_4_3_functions() { return d->_recentFunctions; }

void Viewer::set2DSelectionMode(bool b) { d->is_2d_selection_mode = b; }

void Viewer::setStaticImage(QImage image) { d->static_image = image; }

const QImage& Viewer:: staticImage() const { return d->static_image; }


void Viewer::setCurrentPass(int pass) { d->current_pass = pass; }

void Viewer::setDepthWriting(bool writing_depth) { d->writing_depth = writing_depth; }

void Viewer::setDepthPeelingFbo(QOpenGLFramebufferObject* fbo) { d->dp_fbo = fbo; }

int Viewer::currentPass()const{ return d->current_pass; }
bool Viewer::isDepthWriting()const{ return d->writing_depth; }
QOpenGLFramebufferObject *Viewer::depthPeelingFbo(){ return d->dp_fbo; }
float Viewer::total_pass()
{
  return d->current_total_pass * 1.0f;
}
void Viewer::setTotalPass(int p)
{
  d->total_pass = p;
  update();
}

void Viewer::messageLogged(QOpenGLDebugMessage msg)
{
  QString error;

  // Format based on severity
  switch (msg.severity())
  {
  case QOpenGLDebugMessage::NotificationSeverity:
    return;
    break;
  case QOpenGLDebugMessage::HighSeverity:
    error += "GL ERROR :";
    break;
  case QOpenGLDebugMessage::MediumSeverity:
    error += "GL WARNING :";
    break;
  case QOpenGLDebugMessage::LowSeverity:
    error += "GL NOTE :";
    break;
  default:
    break;
  }

  error += " (";

  // Format based on source
#define CASE(c) case QOpenGLDebugMessage::c: error += #c; break
  switch (msg.source())
  {
  CASE(APISource);
  CASE(WindowSystemSource);
  CASE(ShaderCompilerSource);
  CASE(ThirdPartySource);
  CASE(ApplicationSource);
  CASE(OtherSource);
  CASE(InvalidSource);
  default:
    break;
  }
#undef CASE

  error += " : ";

  // Format based on type
#define CASE(c) case QOpenGLDebugMessage::c: error += #c; break
  switch (msg.type())
  {
  CASE(ErrorType);
  CASE(DeprecatedBehaviorType);
  CASE(UndefinedBehaviorType);
  CASE(PortabilityType);
  CASE(PerformanceType);
  CASE(OtherType);
  CASE(MarkerType);
  CASE(GroupPushType);
  CASE(GroupPopType);
  default:
    break;
  }
#undef CASE

  error += ")";
  qDebug() << qPrintable(error) << "\n" << qPrintable(msg.message()) << "\n";
}

void Viewer::setLighting()
{
  
  //save current settings;
  float prev_spec = d->spec_power;
  QVector4D prev_pos = d->position;
  QVector4D prev_ambient = d->ambient;
  QVector4D prev_diffuse = d->diffuse;
  QVector4D prev_spec_color = d->specular;
  //open dialog
  LightingDialog* dialog = new LightingDialog(d);
  //set specular
  connect(dialog->spec_powrSlider, &QSlider::valueChanged,
          [this, dialog]()
  {
    d->spec_power = dialog->spec_powrSlider->value();
    update();
  });
  //set position
  connect(dialog->position_lineEdit, &QLineEdit::editingFinished,
          [this, dialog]()
  {
    QStringList list = dialog->position_lineEdit->text().split(QRegExp(","), QString::SkipEmptyParts);
    if (list.isEmpty()) return;
    if (list.size()!=3){
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Input should consists of 3 floats.");
      msgBox->exec();
      return;
    }
    double coords[3];
    for(int j=0; j<3; ++j)
    {
      bool ok;
      coords[j] = list.at(j).toFloat(&ok);
      if(!ok)
      {
          QMessageBox *msgBox = new QMessageBox;
          msgBox->setWindowTitle("Error");
          msgBox->setText("ERROR : Coordinates are invalid.");
          msgBox->exec();
          return;
      }
    }
    d->position = QVector4D(coords[0], coords[1], coords[2], 1.0f);
    update();
  });
  

  //set ambient
  connect(dialog, &LightingDialog::s_ambient_changed,
          [this, dialog](){
    d->ambient=QVector4D(dialog->ambient.redF(), 
                         dialog->ambient.greenF(), 
                         dialog->ambient.blueF(), 
                         1.0f);
    update();
  });
          
  //set diffuse
  connect(dialog, &LightingDialog::s_diffuse_changed,
          [this, dialog](){
    d->diffuse=QVector4D(dialog->diffuse.redF(), 
                         dialog->diffuse.greenF(), 
                         dialog->diffuse.blueF(), 
                         1.0f);
    update();
  });
  //set specular
  connect(dialog, &LightingDialog::s_specular_changed,
          [this, dialog](){
    d->specular=QVector4D(dialog->specular.redF() , 
                         dialog->specular.greenF(), 
                         dialog->specular.blueF() , 
                         1.0f);
    update();
    
  });
  
  //reset default
  connect(dialog->buttonBox->button(QDialogButtonBox::StandardButton::RestoreDefaults), &QPushButton::clicked,
          [this](){
    d->position = QVector4D(0,0,1,1);
    d->ambient=QVector4D(77.0/255,77.0/255,77.0/255, 1.0);
    d->diffuse=QVector4D(204.0/255,204.0/255,204.0/255,1.0);
    d->specular=QVector4D(0,0,0,1.0);
    d->spec_power = 51;    
    update();
    
  });
  if(!dialog->exec())
  {
    //restore previous settings
    d->spec_power = prev_spec;
    d->position = prev_pos;
    d->ambient = prev_ambient;
    d->diffuse = prev_diffuse;
    d->specular = prev_spec_color;
    return;
  }
}

void Viewer::setGlPointSize(const GLfloat &p) { d->gl_point_size = p; }

const GLfloat& Viewer::getGlPointSize() const { return d->gl_point_size; }

QVector4D* Viewer::clipBox() const
{
  return d->clipbox;
}

bool Viewer::isClipping() const
{
  return d->clipping;
}

#include "Viewer.moc"

