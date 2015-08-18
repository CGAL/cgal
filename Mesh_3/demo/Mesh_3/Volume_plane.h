#ifndef CGAL_VOLUME_PLANE_H
#define CGAL_VOLUME_PLANE_H

#include <CGAL_demo/Scene_item.h>

#include <vector>
#include <cassert>

#include <QDebug>
#include <QGLViewer/qglviewer.h>

#include "Volume_plane_interface.h"
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#if !defined(NDEBUG)
inline
void printGlError(unsigned int line) {
  GLenum error = glGetError();
  if(error != GL_NO_ERROR)
    std::cerr << gluErrorString(error) << "@" << line << std::endl;
}
#else
inline
void printGlError(unsigned int) {
}
#endif

template<int Dim>
class Length_constraint : public qglviewer::WorldConstraint {
public:
  Length_constraint(double max_) : max_(max_) { }

  void constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const frame) {
    WorldConstraint::constrainTranslation(t, frame);
    qglviewer::Vec pos = frame->position();
    double start = pos[Dim];
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
public:
  Volume_plane(unsigned int adim, unsigned int bdim, unsigned int cdim, 
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

  void draw(Viewer* viewer)const;

  unsigned int aDim() const { return adim_; }
  unsigned int bDim() const { return bdim_; }
  unsigned int cDim() const { return cdim_; }

  qglviewer::Vec translationVector() const { return translationVector(*this); }

  unsigned int getCurrentCube() const { return currentCube; }

  // uses a public init function to make enable construction in
  // threads without gl-context
  void init();

private:
  static const char* vertexShader_source;

  static const char* fragmentShader_source;

  static const char* vertexShader_bordures_source;

  static const char* fragmentShader_bordures_source;


  qglviewer::Vec translationVector(x_tag) const {
    return qglviewer::Vec(xscale_, 0.0, 0.0);
  }
  qglviewer::Vec translationVector(y_tag) const {
    return qglviewer::Vec(0.0, yscale_, 0.0);
  }
  qglviewer::Vec translationVector(z_tag) const {
    return qglviewer::Vec(0.0, 0.0, zscale_);
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

  const unsigned int adim_, bdim_, cdim_;
  const double xscale_, yscale_, zscale_;
  mutable int currentCube;

  mutable QOpenGLBuffer vVBO;
  mutable QOpenGLBuffer cbuffer;
  mutable QOpenGLBuffer rectBuffer;
  mutable std::vector<float> v_rec;
  mutable QOpenGLShaderProgram program_bordures;
  mutable QOpenGLShaderProgram program;
  mutable std::vector< std::pair<QOpenGLBuffer, unsigned int> > ebos;
  std::vector< float > colors_;

  QString name(x_tag) const { return tr("X Slice for %1").arg(name_); }
  QString name(y_tag) const { return tr("Y Slice for %2").arg(name_); }
  QString name(z_tag) const { return tr("Z Slice for %2").arg(name_); }

  void drawRectangle(x_tag) const {


      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back(0.0f); v_rec.push_back((adim_ - 1) * yscale_); v_rec.push_back(0.0f);
      v_rec.push_back(0.0f); v_rec.push_back((adim_ - 1) * yscale_); v_rec.push_back((bdim_ - 1) * zscale_);
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * zscale_);

  }

  void drawRectangle(y_tag) const {
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back((adim_ - 1) * xscale_);v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back((adim_ - 1) * xscale_);v_rec.push_back(0.0f);v_rec.push_back( (bdim_ - 1) * zscale_);
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * zscale_);
  }

  void drawRectangle(z_tag) const {
      v_rec.push_back(0.0f); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back((adim_ - 1) * xscale_); v_rec.push_back(0.0f); v_rec.push_back(0.0f);
      v_rec.push_back((adim_ - 1) * xscale_); v_rec.push_back((bdim_ - 1) * yscale_); v_rec.push_back(0.0f);
      v_rec.push_back(0.0f); v_rec.push_back((bdim_ - 1) * yscale_); v_rec.push_back(0.0f);
  }

  qglviewer::Constraint* setConstraint(x_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<0>(cdim_ * xscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(1.0f, 0.0f, 0.0f));
    return c;
  }
  
  qglviewer::Constraint* setConstraint(y_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<1>(cdim_ * yscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(0.0f, 1.0f, 0.0f));
    return c;
  }

  qglviewer::Constraint* setConstraint(z_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<2>(cdim_ * zscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(0.0f, 0.0f, 1.0f));
    return c;
  }

  void updateCurrentCube() const { currentCube = getTranslation(); }

  GLdouble getTranslation() const { return getTranslation(*this); }
  GLdouble getTranslation(x_tag) const { return mFrame_->matrix()[12] / xscale_; }
  GLdouble getTranslation(y_tag) const { return mFrame_->matrix()[13] / yscale_; }
  GLdouble getTranslation(z_tag) const { return mFrame_->matrix()[14] / zscale_; }

};

template<typename T>
const char* Volume_plane<T>::vertexShader_source =
      "#version 330 \n"
      "in vec4 vertex; \n"
      "in float color; \n"
      "uniform highp mat4 mvp_matrix; \n"
      "uniform highp mat4 f_matrix; \n"
      "out vec4 fullColor; \n"
      "void main() \n"
      "{ gl_Position = mvp_matrix * f_matrix * vertex; \n"
      " fullColor = vec4(color, color, color, 1.0f); } \n";

template<typename T>
const char* Volume_plane<T>::fragmentShader_source =
      "#version 330\n"
      "in vec4 fullColor; \n"
      "void main() { gl_FragColor = fullColor; } \n";

template<typename T>
const char* Volume_plane<T>::vertexShader_bordures_source =
      "#version 330 \n"
      "in vec4 vertex; \n"
      "uniform vec4 color; \n"
      "uniform highp mat4 mvp_matrix; \n"
      "uniform highp mat4 f_matrix; \n"
      "out vec4 fullColor; \n"
      "void main() \n"
      "{ gl_Position = mvp_matrix * f_matrix * vertex; \n"
      " fullColor = color; } \n";

template<typename T>
const char* Volume_plane<T>::fragmentShader_bordures_source =
      "#version 330\n"
      "in vec4 fullColor; \n"
      "void main() { gl_FragColor = fullColor; } \n";



template<typename T>
Volume_plane<T>::Volume_plane(unsigned int adim, unsigned int bdim, unsigned int cdim, 
                              float xscale, float yscale, float zscale, std::vector<float>& colors)
  : Volume_plane_interface(new qglviewer::ManipulatedFrame), adim_(adim), bdim_(bdim), cdim_(cdim), xscale_(xscale), 
    yscale_(yscale), zscale_(zscale), currentCube(0)
{
  colors_.swap(colors);
  mFrame_->setConstraint(setConstraint(*this));
}

template<typename T>
Volume_plane<T>::~Volume_plane() {
  for(std::vector< std::pair< QOpenGLBuffer, unsigned int> >::iterator it = ebos.begin();
      it != ebos.end(); ++it) { 
      it->first.destroy();
  }
  program.release();
}

template<typename T>
void Volume_plane<T>::draw(Viewer* viewer) const {
  updateCurrentCube();



  GLdouble mat[16];
  viewer->camera()->getModelViewProjectionMatrix(mat);
  QMatrix4x4 mvp;
  QMatrix4x4 f;
  for(int i=0; i<16; i++)
  {
      mvp.data()[i] = (float)mat[i];
      f.data()[i] = (float)mFrame_->matrix()[i];
  }


  program_bordures.bind();
  program_bordures.setUniformValue("mvp_matrix", mvp);
  program_bordures.setUniformValue("f_matrix", f);
  program_bordures.release();
  GLint renderMode;
  glGetIntegerv(GL_RENDER_MODE, &renderMode);
  printGlError(__LINE__);


  glLineWidth(4.0f);
  v_rec.resize(0);
  drawRectangle(*this);

  program_bordures.bind();
  rectBuffer.create();
  rectBuffer.bind();
  rectBuffer.allocate(v_rec.data(), static_cast<int>(v_rec.size()*sizeof(float)));
  program_bordures.setAttributeBuffer("vertex",GL_FLOAT,0,3);
  program_bordures.enableAttributeArray("vertex");
  float current_color[4];
  glGetFloatv(GL_CURRENT_COLOR, current_color);
  QColor color;
  color.setRgbF(current_color[0], current_color[1], current_color[2]);
  program_bordures.setUniformValue("color",color);
  glDrawArrays(GL_LINE_LOOP, 0, static_cast<GLsizei>(v_rec.size()/3));
  rectBuffer.release();
  program_bordures.release();
  glLineWidth(1.0f);

  program.bind();
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
  program.setAttributeBuffer(colorLoc, GL_FLOAT, (currentCube*sizeof(float)) * (bdim_) * (adim_), 1, 0 );
  cbuffer.release();


  printGlError(__LINE__);

 for(unsigned int i = 0; i < ebos.size(); ++i)
  {
      ebos[i].first.bind();
      glDrawElements(GL_TRIANGLES, ebos[i].second, GL_UNSIGNED_INT, 0);
      ebos[i].first.release();
  }

  cbuffer.release();
  printGlError(__LINE__);
  program.release();

  printGlError(__LINE__);
}

template<typename T>
void Volume_plane<T>::init() {
  
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

  int maxi, maxv;
  glGetIntegerv(GL_MAX_ELEMENTS_INDICES, &maxi);
  glGetIntegerv(GL_MAX_ELEMENTS_VERTICES, &maxv);
  assert((vertices.size( ) / 3) < (unsigned int)maxi);
  vVBO.create();
  vVBO.bind();
  vVBO.allocate(vertices.data(),static_cast<int>(sizeof(float) * vertices.size()));
  vVBO.release();
  printGlError(__LINE__);

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
    unsigned int left_over = (i + slice) > indices.size()  ? std::distance(indices.begin() + i, indices.end()) : slice;
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

  printGlError(__LINE__);
}

template<typename T>
void Volume_plane<T>::initShaders() {
  QOpenGLShader *vertex = new QOpenGLShader(QOpenGLShader::Vertex);

  vertex->compileSourceCode(vertexShader_source);
  QOpenGLShader *fragment= new QOpenGLShader(QOpenGLShader::Fragment);
  fragment->compileSourceCode(fragmentShader_source);
  program.addShader(vertex);
  program.addShader(fragment);
  program.link();


  QOpenGLShader *vertex_bordures = new QOpenGLShader(QOpenGLShader::Vertex);

  vertex_bordures->compileSourceCode(vertexShader_bordures_source);
  QOpenGLShader *fragment_bordures= new QOpenGLShader(QOpenGLShader::Fragment);
  fragment_bordures->compileSourceCode(fragmentShader_bordures_source);
  program_bordures.addShader(vertex_bordures);
  program_bordures.addShader(fragment_bordures);
  program_bordures.link();
}


#endif /* CGAL_VOLUME_PLANE_H */
