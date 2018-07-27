#ifndef CGAL_VOLUME_PLANE_H
#define CGAL_VOLUME_PLANE_H

#include <CGAL/Three/Scene_item.h>

#include <vector>
#include <cassert>

#include <QDebug>
#include <CGAL/Qt/qglviewer.h>

#include "Volume_plane_interface.h"
#include "create_sphere.h"
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QMouseEvent>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Qt/constraint.h>
#include <CGAL/Qt/debug.h>

using namespace CGAL::Three;

#if !defined(NDEBUG)
inline
void printGlError(CGAL::Three::Viewer_interface* viewer, unsigned int line) {
  CGAL::Qt::opengl_check_errors(viewer, line);
}
#else
inline
void printGlError(CGAL::Three::Viewer_interface*, unsigned int) {
}
#endif

template<int Dim>
class Length_constraint : public CGAL::qglviewer::WorldConstraint {
public:
  Length_constraint(double max_) : max_(max_) { }

  void constrainTranslation(CGAL::qglviewer::Vec& t, CGAL::qglviewer::Frame* const frame) {
    WorldConstraint::constrainTranslation(t, frame);
    CGAL::qglviewer::Vec pos = frame->position();
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    double start = pos[Dim] - offset[Dim];
    double end = t[Dim];
    start += end;
    if(start > max_ || (start < 0)) {
      t[Dim] = 0.0;
    }
  }

private:
  double max_;
};

struct x_tag {};
struct y_tag {};
struct z_tag {};

template<typename Tag>
class Volume_plane : public Volume_plane_interface, public Tag {
private : 
  float tx, ty, tz;
public:
 Volume_plane(float tx, float ty, float tz);

 void invalidateOpenGLBuffers()
 {
     const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
     mFrame_->setPosition(offset.x, offset.y, offset.z);
 }
 bool eventFilter(QObject *, QEvent *);
 void compute_bbox()const
 {
   compute_bbox(*this);
 }

 void compute_bbox(x_tag)const
 {
   _bbox = Bbox(tx,ty,tz,
               tx, ty+(adim_ - 1) * yscale_, tz+(bdim_ - 1) * zscale_);
 }

 void compute_bbox(y_tag)const
 {
   _bbox = Bbox(tx,ty,tz,
               tx+(adim_ - 1) * xscale_, ty, tz+(bdim_ - 1) * zscale_);
 }

 void compute_bbox(z_tag)const
 {
   _bbox = Bbox(tx,ty,tz,
               tx+(adim_ - 1) * xscale_,ty+ (bdim_ - 1) * yscale_, tz);
 }

 void setData(unsigned int adim, unsigned int bdim, unsigned int cdim,
                 float xscale, float yscale, float zscale, std::vector<float>& colors);

  virtual ~Volume_plane();

  Volume_plane* clone() const { return NULL; }

  virtual RenderingMode renderingMode() const { return Flat; }
  bool supportsRenderingMode(RenderingMode m) const { return m == Flat; }
  
  QString toolTip() const { return "Plane through a volume"; }
  QString name() const { return name(*this); }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  bool manipulatable() const { return true; }

  unsigned int cube() {return currentCube; }

  void draw(Viewer_interface* viewer)const;

  unsigned int aDim() const { return adim_; }
  unsigned int bDim() const { return bdim_; }
  unsigned int cDim() const { return cdim_; }

  CGAL::qglviewer::Vec translationVector() const { return translationVector(*this); }

  unsigned int getCurrentCube() const { return currentCube; }

  // uses a public init function to enable construction in
  // threads without gl-context
  void init(Viewer_interface* viewer);

private:
  bool is_grabbing;



  static const char* vertexShader_source;

  static const char* fragmentShader_source;

  static const char* vertexShader_bordures_source;

  static const char* fragmentShader_bordures_source;
  
  static const char* vertexShader_source_comp;

  static const char* fragmentShader_source_comp;

  CGAL::qglviewer::Vec translationVector(x_tag) const {
    return CGAL::qglviewer::Vec(xscale_, 0.0, 0.0);
  }
  CGAL::qglviewer::Vec translationVector(y_tag) const {
    return CGAL::qglviewer::Vec(0.0, yscale_, 0.0);
  }
  CGAL::qglviewer::Vec translationVector(z_tag) const {
    return CGAL::qglviewer::Vec(0.0, 0.0, zscale_);
  }

  void initShaders();

  void buildVertex(std::vector<float>& out, unsigned int i, unsigned int j) {
    buildVertex(out, i, j, *this);
  }

  void buildVertex(std::vector<float>& out, unsigned int i, unsigned int j, x_tag) {
    out.push_back(0.0f);
    out.push_back(i * yscale_);
    out.push_back(j * zscale_);
  }

  void buildVertex(std::vector<float>& out, unsigned int i, unsigned int j, y_tag) {
    out.push_back(i * xscale_);
    out.push_back(0.0f);
    out.push_back(j * zscale_);
  }

  void buildVertex(std::vector<float>& out, unsigned int i, unsigned int j, z_tag) {
    out.push_back(i * xscale_);
    out.push_back(j * yscale_);
    out.push_back(0.0f);
  }

  unsigned int adim_, bdim_, cdim_;
  double xscale_, yscale_, zscale_;
  mutable int currentCube;
  mutable QOpenGLBuffer vVBO;
  mutable QOpenGLBuffer cbuffer;
  mutable QOpenGLBuffer rectBuffer;
  mutable std::vector<float> v_rec;
  mutable std::vector<float> v_spheres;
  mutable std::vector<float> n_spheres;
  mutable std::vector<float> c_spheres;
  mutable QOpenGLShaderProgram* program_bordures;
  mutable QOpenGLShaderProgram program;
  mutable QOpenGLShaderProgram* spheres_program;
  mutable std::vector< std::pair<QOpenGLBuffer, unsigned int> > ebos;
  mutable double sphere_radius;
  std::vector< float > colors_;

  QString name(x_tag) const { return tr("X Slice for %1").arg(name_); }
  QString name(y_tag) const { return tr("Y Slice for %2").arg(name_); }
  QString name(z_tag) const { return tr("Z Slice for %2").arg(name_); }

  double compute_maxDim() const
  {
    double ax((adim_ - 1) * xscale_), ay((adim_ - 1) * yscale_), az((adim_ - 1) * zscale_),
           bx((bdim_ - 1) * xscale_), by((bdim_ - 1) * yscale_), bz((bdim_ - 1) * zscale_),
           cx((cdim_ - 1) * xscale_), cy((cdim_ - 1) * yscale_), cz((cdim_ - 1) * zscale_);

    double max_a = (std::max)((std::max)(ax, ay), az);
    double max_b = (std::max)((std::max)(bx, by), bz);
    double max_c = (std::max)((std::max)(cx, cy), cz);
    return (std::max)((std::max)(max_a, max_b), max_c);

  }
  void drawRectangle(x_tag, bool is_loop) const {

    
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back(0.0f); v_rec.push_back((adim_ - 1) * yscale_); v_rec.push_back(0.0f);
      if(!is_loop){
        v_rec.push_back(0.0f); v_rec.push_back((adim_ - 1) * yscale_); v_rec.push_back(0.0f);
      }
      v_rec.push_back(0.0f); v_rec.push_back((adim_ - 1) * yscale_); v_rec.push_back((bdim_ - 1) * zscale_);
      if(!is_loop){
        v_rec.push_back(0.0f); v_rec.push_back((adim_ - 1) * yscale_); v_rec.push_back((bdim_ - 1) * zscale_);
      }
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * zscale_);
      if(!is_loop)
      {
        v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * zscale_);
        v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      }


  }

  void drawRectangle(y_tag, bool is_loop) const {
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back((adim_ - 1) * xscale_);v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      if(!is_loop){
        v_rec.push_back((adim_ - 1) * xscale_);v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      }
      v_rec.push_back((adim_ - 1) * xscale_);v_rec.push_back(0.0f);v_rec.push_back( (bdim_ - 1) * zscale_);
      if(!is_loop){
        v_rec.push_back((adim_ - 1) * xscale_);v_rec.push_back(0.0f);v_rec.push_back( (bdim_ - 1) * zscale_);
      }
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * zscale_);
      if(!is_loop){
        v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * zscale_);
        v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      }
  }

  void drawRectangle(z_tag, bool is_loop) const {
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back((adim_ - 1) * xscale_); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      if(!is_loop){
        v_rec.push_back((adim_ - 1) * xscale_); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      }
      v_rec.push_back((adim_ - 1) * xscale_); v_rec.push_back((bdim_ - 1) * yscale_); v_rec.push_back(0.0f);
      if(!is_loop){
        v_rec.push_back((adim_ - 1) * xscale_); v_rec.push_back((bdim_ - 1) * yscale_); v_rec.push_back(0.0f);
      }
      v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * yscale_); v_rec.push_back(0.0f);
      if(!is_loop){
        v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * yscale_); v_rec.push_back(0.0f);
        v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      }
  }

  void drawSpheres(x_tag) const
  {
      double max_dim = compute_maxDim();
      sphere_radius = max_dim / 40.0f;
      create_flat_sphere(1.0f, v_spheres, n_spheres, 18);

      c_spheres.push_back(0.0f); c_spheres.push_back((adim_ - 1) * yscale_/2.0f + max_dim/15.0f); c_spheres.push_back(0.0f);
      c_spheres.push_back(0.0f); c_spheres.push_back((adim_ - 1) * yscale_ ); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f + max_dim/15.0f);
      c_spheres.push_back(0.0f); c_spheres.push_back((adim_ - 1) * yscale_/2.0f + max_dim/15.0f); c_spheres.push_back((bdim_ - 1 ) * zscale_);
      c_spheres.push_back(0.0f); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f + max_dim/15.0f);

  }

  void drawSpheres(y_tag) const
  {
      double max_dim = compute_maxDim();
      sphere_radius = max_dim / 40.0f;
      create_flat_sphere(1.0f, v_spheres, n_spheres,18);

      c_spheres.push_back((adim_ - 1) * xscale_/2.0f); c_spheres.push_back(0.0f); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f);
      c_spheres.push_back((adim_ - 1) * xscale_/2.0f); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_);
      c_spheres.push_back(0.0f); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f);
  }

  void drawSpheres(z_tag) const
  {
      double max_dim = compute_maxDim();
      sphere_radius = max_dim / 40.0f;
      create_flat_sphere(1.0f, v_spheres, n_spheres,18);

      c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * yscale_/2.0f - max_dim/15.0f); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_/2.0f-max_dim/15.0f); c_spheres.push_back((bdim_ - 1) * yscale_); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_); c_spheres.push_back((bdim_ - 1) * yscale_/2.0f-max_dim/15.0f); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_/2.0f-max_dim/15.0f); c_spheres.push_back(0.0f); c_spheres.push_back(0.0f);
  }

  CGAL::qglviewer::Constraint* setConstraint(x_tag) {
    CGAL::qglviewer::AxisPlaneConstraint* c = new Length_constraint<0>(cdim_ * xscale_);
    c->setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(CGAL::qglviewer::AxisPlaneConstraint::AXIS, CGAL::qglviewer::Vec(1.0f, 0.0f, 0.0f));
    return c;
  }
  
  CGAL::qglviewer::Constraint* setConstraint(y_tag) {
    CGAL::qglviewer::AxisPlaneConstraint* c = new Length_constraint<1>(cdim_ * yscale_);
    c->setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(CGAL::qglviewer::AxisPlaneConstraint::AXIS, CGAL::qglviewer::Vec(0.0f, 1.0f, 0.0f));
    return c;
  }

  CGAL::qglviewer::Constraint* setConstraint(z_tag) {
    CGAL::qglviewer::AxisPlaneConstraint* c = new Length_constraint<2>(cdim_ * zscale_);
    c->setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(CGAL::qglviewer::AxisPlaneConstraint::AXIS, CGAL::qglviewer::Vec(0.0f, 0.0f, 1.0f));
    return c;
  }

  void updateCurrentCube() const {
      currentCube = getTranslation();

  }

  GLdouble getTranslation() const { return getTranslation(*this); }
  GLdouble getTranslation(x_tag) const {
      const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
      return (mFrame_->matrix()[12]  - offset.x)/ xscale_;
  }
  GLdouble getTranslation(y_tag) const {
      const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
      return (mFrame_->matrix()[13]  - offset.y)/ yscale_;
  }
  GLdouble getTranslation(z_tag) const {
      const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
      return (mFrame_->matrix()[14] - offset.z)/ zscale_ ;
  }

  void initializeBuffers(CGAL::Three::Viewer_interface* viewer) const
  {
      spheres_program = getShaderProgram(PROGRAM_SPHERES, viewer);

      spheres_program->bind();
      vaos[0]->bind();
      buffers[0].bind();
      buffers[0].allocate(v_spheres.data(),
                                 static_cast<int>(v_spheres.size()*sizeof(float)));
      spheres_program->enableAttributeArray("vertex");
      spheres_program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
      buffers[0].release();

      buffers[1].bind();
      buffers[1].allocate(n_spheres.data(),
                                static_cast<int>(n_spheres.size()*sizeof(float)));
      spheres_program->enableAttributeArray("normals");
      spheres_program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
      buffers[1].release();

      buffers[2].bind();
      buffers[2].allocate(c_spheres.data(),
                               static_cast<int>(c_spheres.size()*sizeof(float)));
      spheres_program->enableAttributeArray("center");
      spheres_program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
      buffers[2].release();

      viewer->glVertexAttribDivisor(spheres_program->attributeLocation("center"), 1);
      vaos[0]->release();
      spheres_program->release();
      are_buffers_filled = true;
  }
};

template<typename T>
const char* Volume_plane<T>::vertexShader_source =
      "#version 150  \n"
      "in vec4 vertex; \n"
      "in float color; \n"
      "uniform mat4 mvp_matrix; \n"
      "uniform mat4 f_matrix; \n"
      "out vec4 fullColor; \n"
      "void main() \n"
      "{ gl_Position = mvp_matrix * f_matrix * vertex; \n"
      " fullColor = vec4(color, color, color, 1.0); } \n";
template<typename T>
const char* Volume_plane<T>::vertexShader_source_comp =
      "attribute highp vec4 vertex; \n"
      "attribute highp float color; \n"
      "uniform highp mat4 mvp_matrix; \n"
      "uniform highp mat4 f_matrix; \n"
      "varying highp vec4 fullColor; \n"
      "void main() \n"
      "{ gl_Position = mvp_matrix * f_matrix * vertex; \n"
      " fullColor = vec4(color, color, color, 1.0); } \n";

template<typename T>
const char* Volume_plane<T>::fragmentShader_source =
      "#version 150 \n"
      "in vec4 fullColor; \n"
      "out vec4 out_color; \n"
      "void main() { out_color = fullColor; } \n";
template<typename T>
const char* Volume_plane<T>::fragmentShader_source_comp =
      "varying highp vec4 fullColor; \n"
      "void main() { gl_FragColor = fullColor; } \n";

template<typename T>
const char* Volume_plane<T>::vertexShader_bordures_source =
      "attribute highp vec4 vertex; \n"
      "attribute highp vec4 colors; \n"
      "uniform highp mat4 mvp_matrix; \n"
      "uniform highp mat4 f_matrix; \n"
      "varying highp vec4 fullColor; \n"
      "void main() \n"
      "{ gl_Position = mvp_matrix * f_matrix * vertex; \n"
      " fullColor = colors; } \n";


template<typename T>
const char* Volume_plane<T>::fragmentShader_bordures_source =
      "varying highp vec4 fullColor; \n"
      "void main() { gl_FragColor = fullColor; } \n";


template<typename T>
Volume_plane<T>::Volume_plane(float tx, float ty, float tz)
  : Volume_plane_interface(new CGAL::qglviewer::ManipulatedFrame),
    tx(tx), ty(ty), tz(tz)    
 {
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    mFrame_->setPosition(offset.x, offset.y, offset.z);
    sphere_Slider = new QSlider(Qt::Horizontal);
    sphere_Slider->setValue(40);
    sphere_Slider->setMinimum(1);
    sphere_Slider->setMaximum(100);
    program_bordures = NULL;
 }
template<typename T>
void Volume_plane<T>::setData(unsigned int adim, unsigned int bdim, unsigned int cdim, float xscale, float yscale, float zscale, std::vector<float> &colors)
{
 adim_ = adim;
 bdim_ = bdim;
 cdim_= cdim;
 xscale_ = xscale;
 yscale_ =yscale;
 zscale_ = zscale;
 currentCube = 0;
 colors_.swap(colors);
 mFrame_->setConstraint(setConstraint(*this));
 v_spheres.resize(0);
 n_spheres.resize(0);
 c_spheres.resize(0);
 drawSpheres(*this);
}
template<typename T>
Volume_plane<T>::~Volume_plane() {
  for(std::vector< std::pair< QOpenGLBuffer, unsigned int> >::iterator it = ebos.begin();
      it != ebos.end(); ++it) { 
      it->first.destroy();
  }
  program.release();
  delete sphere_Slider;
}

template<typename T>
void Volume_plane<T>::draw(Viewer_interface *viewer) const {
  updateCurrentCube();


  viewer->glDepthRangef(0.00001f,0.99999f);
  GLdouble mat[16];
  viewer->camera()->getModelViewProjectionMatrix(mat);
  QMatrix4x4 mvp;
  QMatrix4x4 f;
  for(int i=0; i<16; i++)
  {
      mvp.data()[i] = (float)mat[i];
      f.data()[i] = (float)mFrame_->matrix()[i];
  }
  mvp.translate(QVector3D(tx, ty, tz));

 if(!program_bordures)
  {
    if(viewer->isOpenGL_4_3())
      program_bordures = viewer->getShaderProgram(PROGRAM_SOLID_WIREFRAME);
    else
    {
      program_bordures = new QOpenGLShaderProgram(viewer);
      program_bordures->addShaderFromSourceCode(QOpenGLShader::Vertex,vertexShader_bordures_source);
      program_bordures->addShaderFromSourceCode(QOpenGLShader::Fragment,fragmentShader_bordures_source);
      
      program_bordures->link();
    }
  }
  program_bordures->bind();
  program_bordures->setUniformValue("mvp_matrix", mvp);
  program_bordures->setUniformValue("f_matrix", f);
  program_bordures->release();
  GLint renderMode;
  viewer->glGetIntegerv(GL_RENDER_MODE, &renderMode);
  printGlError(viewer, __LINE__);

  v_rec.resize(0);
  drawRectangle(*this, !viewer->isOpenGL_4_3());

  program_bordures->bind();
  vaos[2]->bind();
  rectBuffer.create();
  rectBuffer.bind();
  rectBuffer.allocate(v_rec.data(), static_cast<int>(v_rec.size()*sizeof(float)));
  program_bordures->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  program_bordures->enableAttributeArray("vertex");
  program_bordures->disableAttributeArray("colors");
  program_bordures->setAttributeValue("colors",this->color());
  if(!viewer->isOpenGL_4_3())
  {
    viewer->glDrawArrays(GL_LINE_LOOP, 0, static_cast<GLsizei>(v_rec.size()/3));
  }
  else
  {
    QVector2D vp(viewer->width(), viewer->height());
    program_bordures->setUniformValue("viewport", vp);
    program_bordures->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
    program_bordures->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
    program_bordures->setUniformValue("width", 4.0f);
    viewer->glDepthRangef(0.00005f, 0.99995f);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_rec.size()/3));
    viewer->glDepthRangef(0.0,1.0);
  }
  rectBuffer.release();
  vaos[1]->release();
  program_bordures->release();

  program.bind();
  vaos[1]->bind();
  int mvpLoc = program.uniformLocation("mvp_matrix");
  int fLoc = program.uniformLocation("f_matrix");
  program.setUniformValue(mvpLoc, mvp);
  program.setUniformValue(fLoc, f);
  vVBO.bind();
  int vloc = program.attributeLocation("vertex");
  program.enableAttributeArray(vloc);
  program.setAttributeBuffer(vloc, GL_FLOAT, 0, 3);
  vVBO.release();
  cbuffer.bind();
  int colorLoc = program.attributeLocation("color");
  program.enableAttributeArray(colorLoc);
  program.setAttributeBuffer(colorLoc, GL_FLOAT, (currentCube*int(sizeof(float))) * (bdim_) * (adim_), 1, 0 );
  cbuffer.release();



  printGlError(viewer, __LINE__);

 for(unsigned int i = 0; i < ebos.size(); ++i)
  {
      ebos[i].first.bind();
      viewer->glDrawElements(GL_TRIANGLES, ebos[i].second, GL_UNSIGNED_INT, 0);
      ebos[i].first.release();
  }

  cbuffer.release();
  printGlError(viewer, __LINE__);
  vaos[1]->release();
  program.release();

  printGlError(viewer, __LINE__);

  if(!are_buffers_filled)
      initializeBuffers(viewer);
  mvp = mvp*f;
  //hide spheres if only 1 plane.
  if(aDim() <= 1 ||
     bDim() <= 1 ||
     cDim() <=1)
    return;
  vaos[0]->bind();
  spheres_program = getShaderProgram(PROGRAM_SPHERES, viewer);
  attribBuffers(viewer, PROGRAM_SPHERES);
  spheres_program->bind();
  double max_dim = compute_maxDim();
  sphere_radius = max_dim/20.0f * sphere_Slider->value()/100.0f;
  spheres_program->setAttributeValue("radius", sphere_radius);
  spheres_program->disableAttributeArray("colors");
  spheres_program->setAttributeValue("colors", this->color());
  spheres_program->setUniformValue("mvp_matrix", mvp);
  viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                static_cast<GLsizei>(v_spheres.size()/3),
                                static_cast<GLsizei>(4));
  spheres_program->release();
  vaos[0]->release();
}

template<typename T>
void Volume_plane<T>::init(Viewer_interface* viewer) {
  viewer->makeCurrent();
  is_grabbing = false;
  initShaders();

  // for each vertex
  std::vector< float > vertices;
  vertices.reserve(bdim_ * adim_ * 3);
  for(unsigned int i = 0; i < adim_; ++i) 
  {
    for(unsigned int j = 0; j < bdim_; ++j)
    {
      buildVertex(vertices, i, j);
    }
  }
    
  assert(vertices.size() == (3 * adim_ * bdim_));

  vVBO.create();
  vaos[1]->bind();
  vVBO.bind();
  vVBO.allocate(vertices.data(),static_cast<int>(sizeof(float) * vertices.size()));
  vVBO.release();
  vaos[1]->release();
  printGlError(viewer, __LINE__);

  // for each patch
  std::vector<unsigned int> indices;
  for(unsigned int j = 0; j < adim_ - 1; ++j) {
    for(unsigned int k = 0; k < bdim_ - 1; ++k) {
        //0
        indices.push_back( j * bdim_ + k );
        assert(indices.back() < (vertices.size() / 3));

        //1
        indices.push_back( j * bdim_ + (k + 1) );
        assert(indices.back() < (vertices.size() / 3));

        //3
        indices.push_back( (j+1) * bdim_ + (k+1) );
        assert(indices.back() < (vertices.size() / 3));

        //0
        indices.push_back( j * bdim_ + k );
        assert(indices.back() < (vertices.size() / 3));

        //3
        indices.push_back( (j+1) * bdim_ + (k+1) );
        assert(indices.back() < (vertices.size() / 3));

        //2
        indices.push_back( (j+1) * bdim_ + (k) );
        assert(indices.back() < (vertices.size() / 3));

    }
  }

  assert((indices.size() / 6) == (adim_ - 1) * (bdim_ - 1));
  //slice must be multiple of 3.
  const unsigned int slice = 63399;
  for(unsigned int i = 0; i < indices.size(); i+=slice)
  {
    QOpenGLBuffer ebo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
    unsigned int left_over =
      (i + slice) > indices.size()  ?
      (unsigned int)(std::distance(indices.begin() + i, indices.end())) :
      slice;
    ebo.create();
    ebo.bind();
    ebo.allocate(&indices[i],static_cast<int>(sizeof(unsigned int) * left_over));
    ebo.release();
    ebos.push_back(std::make_pair(ebo, left_over));
  }

  cbuffer.create();
  cbuffer.bind();
  cbuffer.allocate(colors_.data(),static_cast<int>(colors_.size()*sizeof(float)));
  cbuffer.release();

  printGlError(viewer, __LINE__);
}

template<typename T>
void Volume_plane<T>::initShaders() {
  if(QOpenGLContext::currentContext()->format().majorVersion() >= 3)
  {
    program.addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShader_source  );
    program.addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShader_source);
  }
  else
  {
    program.addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShader_source_comp  );
    program.addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShader_source_comp);
  }
  program.link();

}


//intercept events for picking
template<typename T>
bool Volume_plane<T>::eventFilter(QObject *, QEvent *event)
{
    CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
    if(event->type() == QEvent::MouseButtonPress)
    {
        QMouseEvent* e = static_cast<QMouseEvent*>(event);
        if(e->modifiers() == Qt::NoModifier)
        {
            //pick
            bool found = false;
            viewer->makeCurrent();
            CGAL::qglviewer::Vec pos = viewer->camera()->pointUnderPixel(e->pos(), found);
            if(!found)
                return false;
            found = false;
            //check the found point is on a sphere
            pos = manipulatedFrame()->inverse().inverseCoordinatesOf(pos);
            float sq_radius = sphere_radius * sphere_radius;
            CGAL::qglviewer::Vec center;
            //is the picked point on the sphere ?
            for (int i=0; i<4; ++i)
            {
                center = CGAL::qglviewer::Vec(c_spheres[i*3]+tx, c_spheres[i*3+1]+ty, c_spheres[i*3+2]+tz);
                if(
                        (center.x - pos.x) * (center.x - pos.x) +
                        (center.y - pos.y) * (center.y - pos.y) +
                        (center.z - pos.z) * (center.z - pos.z)
                        <= sq_radius)
                {
                    found = true;
                    break;
                }
            }
            if(!found)
                return false;
            is_grabbing = true;
            emitSelection();
            viewer->setManipulatedFrame(manipulatedFrame());
            viewer->setMouseBinding(
                        Qt::NoModifier,
                        Qt::LeftButton,
                        CGAL::qglviewer::FRAME,
                        CGAL::qglviewer::TRANSLATE);
        }
    }
    else if(event->type() == QEvent::MouseButtonRelease && is_grabbing)
    {
        is_grabbing = false;
        viewer->setMouseBinding(
                    Qt::NoModifier,
                    Qt::LeftButton,
                    CGAL::qglviewer::CAMERA,
                    CGAL::qglviewer::ROTATE);
    }
    return false;
}

#endif /* CGAL_VOLUME_PLANE_H */
