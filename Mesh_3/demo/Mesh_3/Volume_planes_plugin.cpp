#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include <CGAL_demo/Scene_interface.h>
#include <CGAL_demo/Scene_item.h>

#include <GL/glew.h>
#include <CGAL/gl.h>
#include <CGAL/Image_3.h>

#include "Scene_segmented_image_item.h"

#include <QAction>
#include <QMenu>
#include <QList>
#include <QMainWindow>
#include <QString>
#include <QDebug>

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <map>
#include <string>

#include <QGLViewer/qglviewer.h>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/unordered_map.hpp>

namespace {
#if !defined(NDEBUG)
  void printGlError(unsigned int line) {
    GLenum error = glGetError();
    if(error != GL_NO_ERROR)
      qDebug() << gluErrorString(error) << "@" << line;
  }
#else
  void printGlError(unsigned int) {
  }
#endif

template<typename T>
struct Clamp_to_one_zero_range {
  std::pair<T, T> min_max;

  T operator()(const T& inVal) {
    T inValNorm = inVal - min_max.first;
    T aUpperNorm = min_max.second - min_max.first;
    T bValNorm = inValNorm / aUpperNorm;
    return bValNorm;
  }
};

const char* vertexShader = 
                     "#version 330 \n"
                     "layout(location = 0) in vec3 position; \n"
                     "layout(location = 1) in float color; \n"
                     "uniform mat4 modelviewMatrix; \n"
                     "uniform mat4 projectionMatrix; \n"
                     "out vec4 fullColor; \n"
                     "void main() \n"
                     "{ gl_Position = (projectionMatrix * modelviewMatrix) * vec4(position.x, position.y, position.z, 1.0); \n"
                     "  fullColor = vec4(color, color, color, 1.0f); } \n";

const char* fragmentShader = 
                     "#version 330 \n"
                     "smooth in vec4 fullColor; \n"
                     "out vec4 outputColor; \n"
                     "void main() { outputColor = fullColor; } \n";

struct XSel { double& operator()(qglviewer::Vec& v) const { return v.x; } };
struct YSel { double& operator()(qglviewer::Vec& v) const { return v.y; } };
struct ZSel { double& operator()(qglviewer::Vec& v) const { return v.z; } };

template<typename DimSelector>
class Length_constraint : public qglviewer::WorldConstraint {
public:
  Length_constraint(double max) : max(max) { }

  void constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* frame) {
    WorldConstraint::constrainTranslation(t, frame);
    qglviewer::Vec pos = frame->position();
    double start = DimSelector()(pos);
    double end = DimSelector()(t);
    start += end;

    if(start > max || (start < 0)) {
        DimSelector()(t) = 0.0;
    }
  }
private:
  double max;
};

struct x_tag {};
struct y_tag {};
struct z_tag {};

} // unnamed namespace

template<typename Tag>
class Volume_plane : public Scene_item, public Tag {
public:
  Volume_plane(unsigned int adim, unsigned int bdim, unsigned int cdim, 
               float xscale, float yscale, float zscale)
    : adim_(adim), bdim_(bdim), cdim_(cdim), xscale_(xscale), 
      yscale_(yscale), zscale_(zscale), currentCube(0),
      mFrame_(new qglviewer::ManipulatedFrame) { 

    mFrame_->setConstraint(setConstraint(*this));

    initShaders();

    // for each vertex
    vertices.reserve(bdim_ * adim_ * 3);
    for(unsigned int i = 0; i < adim_; ++i) 
    {
      for(unsigned int j = 0; j < bdim_; ++j)
      {
        buildVertex(vertices, i, j);
      }
    }
    
    std::cout << "vertex size: " << vertices.size() << std::endl;
    assert(vertices.size() == (3 * adim_ * bdim_));

    int maxi, maxv;
    glGetIntegerv(GL_MAX_ELEMENTS_INDICES, &maxi);
    glGetIntegerv(GL_MAX_ELEMENTS_VERTICES, &maxv);
    assert((vertices.size( ) / 3) < (unsigned int)maxi);

    // glGenBuffers(1, &vVBO);
    // glBindBuffer(GL_ARRAY_BUFFER, vVBO);

    // glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), &vertices[0], GL_STATIC_DRAW);
    // glBindBuffer(GL_ARRAY_BUFFER, 0);
    printGlError(__LINE__);

    // for each patch
    std::vector<unsigned int> indices;
    for(unsigned int j = 0; j < adim_ - 1; ++j) {
      for(unsigned int k = 0; k < bdim_ - 1; ++k) {
        indices.push_back( j * bdim_ + k );
        assert(indices.back() < (vertices.size() / 3));
    
        indices.push_back( j * bdim_ + (k + 1) );
        assert(indices.back() < (vertices.size() / 3));
    
        indices.push_back( (j+1) * bdim_ + (k+1) );
        assert(indices.back() < (vertices.size() / 3));
    
        indices.push_back( (j+1) * bdim_ + (k) );
        assert(indices.back() < (vertices.size() / 3));
      }
    }

    assert((indices.size() / 4) == (adim_ - 1) * (bdim_ - 1));

    const unsigned int slice = 64000;
    for(unsigned int i = 0; i < indices.size(); i+=slice)
    {
      GLuint ebo;
      unsigned int left_over = (i + slice) > indices.size()  ? std::distance(indices.begin() + i, indices.end()) : slice;
      glGenBuffers(1, &ebo);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * left_over, &indices[i], GL_STATIC_DRAW);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
      ebos.push_back(std::make_pair(ebo, left_over));
    }

    printGlError(__LINE__);
  }

  virtual RenderingMode renderingMode() const { return Flat; }

  Volume_plane* clone() const {
    return NULL;
  }

  bool supportsRenderingMode(RenderingMode m) const { return m == Flat; }
  
  QString toolTip() const { return "Plane through a volume"; }
  QString name() const { return name(*this); }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  bool manipulatable() const { return true; }
  ManipulatedFrame* manipulatedFrame() { return mFrame_; }

  void draw() const {
    updateCurrentCube();

    glDisable(GL_LIGHTING);
    glPushMatrix();
    glMultMatrixd(mFrame_->matrix());

    drawRectangle(*this);

    glUseProgram(program);
    printGlError(__LINE__);

    float matrices[32];
    glGetFloatv(GL_MODELVIEW_MATRIX, &matrices[0]);
    glGetFloatv(GL_PROJECTION_MATRIX, &matrices[16]);
    glUniformMatrix4fv(glGetUniformLocation(program, "modelviewMatrix"), 1, GL_FALSE, &matrices[0]);
    glUniformMatrix4fv(glGetUniformLocation(program, "projectionMatrix"), 1, GL_FALSE, &matrices[16]);

    // glBindBuffer(GL_ARRAY_BUFFER, vVBO);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, &vertices[0]);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, &(colors[ currentCube * adim_ * bdim_ ]));
    glEnableVertexAttribArray(1);

    printGlError(__LINE__);

    for(unsigned int i = 0; i < ebos.size(); ++i)
    {
       glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebos[i].first);
       glDrawElements(GL_QUADS, ebos[i].second, GL_UNSIGNED_INT, 0);
    }

    printGlError(__LINE__);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glPopMatrix();

    glUseProgram(0);
    // glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glEnable(GL_LIGHTING);

    printGlError(__LINE__);
  }

private:
  void initShaders() {
    program = glCreateProgram();
    GLint status;
    char buffer[1024];
    int len;

    GLuint vertex = glCreateShader(GL_VERTEX_SHADER);
    GLint* lengths = { 0 };
    glShaderSource(vertex, 1, &vertexShader, lengths);
    glCompileShader(vertex);
    glGetShaderiv(vertex, GL_COMPILE_STATUS, &status);

    glGetShaderInfoLog(vertex, 1024, &len, buffer);
    std::cout << "VertexShader: " << buffer << std::endl;
    assert(status == GL_TRUE);

    printGlError(__LINE__);

    GLuint fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, &fragmentShader, lengths);
    glCompileShader(fragment);
    glGetShaderiv(fragment, GL_COMPILE_STATUS, &status);

    glGetShaderInfoLog(fragment, 1024, &len, buffer);
    std::cout << "FragmentShader: " << buffer << std::endl;

    assert(status == GL_TRUE);

    printGlError(__LINE__);

    glAttachShader(program, vertex);
    glAttachShader(program, fragment);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    assert(status == GL_TRUE);

    printGlError(__LINE__);
  }

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
  ManipulatedFrame* mFrame_;
  std::vector<float> vertices;

  // currently disabled as qglviewer doesn't like selecting with a vbo
  // enabled

  // GLuint vVBO;
  GLuint program;
  std::vector< std::pair<GLuint, unsigned int> > ebos;

  QString name(x_tag) const { return "X Slice"; }
  QString name(y_tag) const { return "Y Slice"; }
  QString name(z_tag) const { return "Z Slice"; }

  void drawRectangle(x_tag) const {
    glLineWidth(3.0f);
    glBegin(GL_LINE_LOOP);
    {
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, (adim_ - 1) * yscale_, 0.0f);
      glVertex3f(0.0f, (adim_ - 1) * yscale_, (bdim_ - 1) * zscale_);
      glVertex3f(0.0f, 0.0f, (bdim_ - 1) * zscale_);
    }
    glEnd();
    glLineWidth(1.0f);
  }

  void drawRectangle(y_tag) const {
    glLineWidth(3.0f);
    glBegin(GL_LINE_LOOP);
    {
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f((adim_ - 1) * xscale_, 0.0f, 0.0f);
      glVertex3f((adim_ - 1) * xscale_, 0.0f, (bdim_ - 1) * zscale_);
      glVertex3f(0.0f, 0.0f, (bdim_ - 1) * zscale_);
    }
    glEnd();
    glLineWidth(1.0f);
  }

  void drawRectangle(z_tag) const {
    glLineWidth(3.0f);
    glBegin(GL_LINE_LOOP);
    {
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f((adim_ - 1) * xscale_, 0.0f, 0.0f);
      glVertex3f((adim_ - 1) * xscale_, (bdim_ - 1) * yscale_, 0.0f);
      glVertex3f(0.0f, (bdim_ - 1) * yscale_, 0.0f);
    }
    glEnd();
    glLineWidth(1.0f);
  }

  qglviewer::Constraint* setConstraint(x_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<XSel>(cdim_ * xscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(1.0f, 0.0f, 0.0f));
    return c;
  }
  qglviewer::Constraint* setConstraint(y_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<YSel>(cdim_ * yscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(0.0f, 1.0f, 0.0f));
    return c;
  }
  qglviewer::Constraint* setConstraint(z_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<ZSel>(cdim_ * zscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(0.0f, 0.0f, 1.0f));
    return c;
  }

  void updateCurrentCube() const {
    // to which cube does this translation take us?
    currentCube = getTranslation(*this);
  }

  GLdouble getTranslation(x_tag) const {
    return mFrame_->matrix()[12] / xscale_;
  }
  GLdouble getTranslation(y_tag) const {
    return mFrame_->matrix()[13] / yscale_;
  }
  GLdouble getTranslation(z_tag) const {
    return mFrame_->matrix()[14] / zscale_;
  }
public:
  std::vector< float > colors;
};

class Volume_plane_plugin :
  public QObject,
  public Plugin_interface 
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface);
public:
  Volume_plane_plugin() : planeSwitch(NULL), sc(NULL), mw(NULL) {
  }
  
  virtual ~Volume_plane_plugin() { }

  virtual void init(QMainWindow* mw, Scene_interface* sc) {
    assert(mw != NULL);
    assert(sc != NULL);
    this->sc = sc;
    this->mw = mw;

    QList<QMenu*> menus = mw->findChildren<QMenu*>();

    planeSwitch = new QAction(mw);
    planeSwitch->setText("Add Volume Planes");

    QMenu* viewMenu = mw->findChild<QMenu*>("menuView");
    assert(viewMenu != NULL);
    viewMenu->addAction(planeSwitch);

    connect(planeSwitch, SIGNAL(triggered()), this, SLOT(addVolumePlanes()));
  }

  virtual QList<QAction*> actions() const {
    return QList<QAction*>() << planeSwitch;
  }

public slots:
  void addVolumePlanes() {
    Scene_segmented_image_item* seg_img = NULL;
    
    for(unsigned int i = 0; i < sc->numberOfEntries(); ++i) {
      if((seg_img = qobject_cast<Scene_segmented_image_item*>(sc->item(i))) != NULL)
        break;
    }

    if(!(seg_img == NULL)) {
      const CGAL::Image_3* img = seg_img->image();

      switch(img->image()->wordKind) {
      case WK_FIXED: 
      {
        typedef unsigned char Word; 
        this->addImageToSc<Word>(img);
        break;
      }
      case WK_FLOAT:
      {
        typedef double Word; 
        this->addImageToSc<Word>(img);
        break;
      }
      default: 
      { 
        std::cout << "Unknown word type" << std::endl;
      }
      }

    } else {
      throw std::runtime_error("No suitable data for Volume planes found");
    }
  }
 
private:
  QAction* planeSwitch;
  Scene_interface* sc;
  QMainWindow* mw;

  template<typename Word> 
  void addImageToSc(const CGAL::Image_3* img) {
    const Word* begin = (const Word*)img->data();
    const Word* end = (const Word*)img->data() + img->size();
        
    std::vector< float > shaped;

    Clamp_to_one_zero_range<float> clamper = { std::make_pair(*std::max_element(begin, end),
                                                              *std::min_element(begin, end)) };

    // mangle the array into a shape suitable for the different axes
    // x
    shaped.reserve(img->size() * 3);
    for(unsigned int i = 0; i < img->xdim(); ++i) {
      for(unsigned int j = 0; j < img->ydim(); ++j) {
        for(unsigned int k = 0; k < img->zdim(); ++k) {
          float x = CGAL::IMAGEIO::static_evaluate<Word>(img->image(), i, j, k);
          x = clamper(x);
          shaped.push_back(x);
        }
      }
    }

    Volume_plane<x_tag>* x = new Volume_plane<x_tag>(img->ydim(), img->zdim(), img->xdim(), 
                                                     img->vx(), img->vy(), img->vz());
    (x->colors).swap(shaped);
    sc->addItem(x);

    // y
    shaped.reserve(img->size() * 3);
    for(unsigned int i = 0; i < img->ydim(); ++i) {
      for(unsigned int j = 0; j < img->xdim(); ++j) {
        for(unsigned int k = 0; k < img->zdim(); ++k) {
          float x = CGAL::IMAGEIO::static_evaluate<Word>(img->image(), j, i, k);
          x = clamper(x);
          shaped.push_back(x);
        }
      }
    }
    
    Volume_plane<y_tag>* y = new Volume_plane<y_tag>(img->xdim(), img->zdim(), img->ydim(), 
                                                     img->vx(), img->vy(), img->vz());
    
    (y->colors).swap(shaped);
    sc->addItem(y);

    shaped.reserve(img->size() * 3);
    for(unsigned int i = 0; i < img->zdim(); ++i) {
      for(unsigned int j = 0; j < img->xdim(); ++j) {
        for(unsigned int k = 0; k < img->ydim(); ++k) {
          float x = CGAL::IMAGEIO::static_evaluate<Word>(img->image(), j, k, i);
          x = clamper(x);
          shaped.push_back(x);
        }
      }
    }
    
    Volume_plane<z_tag>* z = new Volume_plane<z_tag>(img->xdim(), img->ydim(), img->zdim(), 
                                                     img->vx(), img->vy(), img->vz());
    
    (z->colors).swap(shaped);
    sc->addItem(z);
  }
};


Q_EXPORT_PLUGIN2(Volume_plane_plugin, Volume_plane_plugin)

#include "Volume_planes_plugin.moc"
