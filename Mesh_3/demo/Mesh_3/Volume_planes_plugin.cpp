#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include <CGAL_demo/Scene_interface.h>
#include <CGAL_demo/Scene_item.h>
#include "Scene_segmented_image_item.h"

#include <CGAL/Image_3.h>


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
  void printGlError(unsigned int line) {
    GLenum eenum = glGetError();
    if(eenum != GL_NO_ERROR) {
      QString s = (const char*)gluErrorString(eenum);
      qDebug() << s << "@" << line;
    }
  }
}

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

template<typename Tag>
class Volume_plane : public Scene_item, public Tag {
public:
  Volume_plane(unsigned int adim, unsigned int bdim, unsigned int cdim, 
               float xscale, float yscale, float zscale)
    : adim_(adim), bdim_(bdim), cdim_(cdim), xscale_(xscale), 
      yscale_(yscale), zscale_(zscale), currentCube(0),
      mFrame_(new qglviewer::ManipulatedFrame){ 
    mFrame_->setConstraint(setConstraint(*this));
    
    // for each patch
    vertices.reserve(bdim_ * adim_ * 3);
    for(unsigned int j = 0; j < bdim_ - 1; ++j)
    {
      for(unsigned int i = 0; i < adim_ - 1; ++i) {
        buildPatch(i, j, *this);
       }
    }
  }

  virtual RenderingMode renderingMode() const { 
    return Flat;
  }

  Volume_plane* clone() const {
    return NULL;
  }

  bool supportsRenderingMode(RenderingMode m) const {
    if(m == Flat)
      return true;
    else
      return false;
  }
  
  QString toolTip() const {
    return "Delicious planes";
  }

  QString name() const {
    return name(*this);
  }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  bool manipulatable() const { return true; }
  ManipulatedFrame* manipulatedFrame() { 
    return mFrame_;
  }

  void draw() const {
    printGlError(__LINE__);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(mFrame_->matrix());

    updateCurrentCube();

    assert( vertices.size() == colors.size() );

    glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);
    glEnableClientState(GL_VERTEX_ARRAY);

    enableColorPointer();
    glEnableClientState(GL_COLOR_ARRAY);

    glNormalPointer(GL_FLOAT, 0, &normals[0]);
    glEnableClientState(GL_NORMAL_ARRAY);

    printGlError(__LINE__);

    glDrawArrays(GL_QUADS, 0, vertices.size() / 3);

    printGlError(__LINE__);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
        
    printGlError(__LINE__);

    glPopMatrix();
    printGlError(__LINE__);
  }

private:
  const unsigned int adim_, bdim_, cdim_;
  const double xscale_, yscale_, zscale_;
protected:
  mutable int currentCube;
private:
  ManipulatedFrame* mFrame_;
  std::vector< float > vertices;
  std::vector< float > normals;

  std::pair<float, float> min_max, zero_one;
  
  QString name(x_tag) const { return "X Slice"; }
  QString name(y_tag) const { return "Y Slice"; }
  QString name(z_tag) const { return "Z Slice"; }

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
    int tmp = getTranslation(*this);
    // to which cube does this translation take us?
    if(tmp != currentCube) {
      currentCube = tmp;
      buildColors(adim_, bdim_);
    }
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

  void buildPatch(int i, int j, x_tag) {
    for(int k = 0; k < 4; ++k)
    {
      normals.push_back(1.0f);
      normals.push_back(0.0f);
      normals.push_back(0.0f);
    }
    
    vertices.push_back(0.0f);
    vertices.push_back(i * yscale_);
    vertices.push_back(j* zscale_);
        
    vertices.push_back(0.0f);
    vertices.push_back((i + 1) * yscale_);
    vertices.push_back(j * zscale_);

    vertices.push_back(0.0f);
    vertices.push_back((i + 1) * yscale_);
    vertices.push_back((j + 1) * zscale_);

    vertices.push_back(0.0f);
    vertices.push_back(i * yscale_);
    vertices.push_back((j + 1) * zscale_);
  }

  void buildPatch(int i, int j, y_tag) {
    for(int k = 0; k < 4; ++k)
    {
      normals.push_back(0.0f);
      normals.push_back(1.0f);
      normals.push_back(0.0f);
    }
    
    vertices.push_back(i * xscale_);
    vertices.push_back(0.0f);
    vertices.push_back(j * zscale_);

    vertices.push_back((i + 1) * xscale_);
    vertices.push_back( 0.0f);
    vertices.push_back( j * zscale_);

    vertices.push_back((i + 1) * xscale_);
    vertices.push_back( 0.0f);
    vertices.push_back( (j + 1) * zscale_);

    vertices.push_back(i * xscale_);
    vertices.push_back( 0.0f);
    vertices.push_back( (j + 1) * zscale_);
  }

  void buildPatch(int i, int j, z_tag) {
    for(int k = 0; k < 4; ++k)
    {
      normals.push_back(0.0f);
      normals.push_back(0.0f);
      normals.push_back(1.0f);
    }

    vertices.push_back(i * xscale_);
    vertices.push_back(j * yscale_);
    vertices.push_back(0.0f);

    vertices.push_back((i + 1) * xscale_);
    vertices.push_back(j * yscale_);
    vertices.push_back(0.0f);

    vertices.push_back((i + 1) * xscale_);
    vertices.push_back((j + 1) * yscale_);
    vertices.push_back(0.0f);

    vertices.push_back(i * xscale_);
    vertices.push_back((j + 1) * yscale_);
    vertices.push_back(0.0f);
  }

protected:
  virtual void buildColors(unsigned int, unsigned int) const = 0;
  virtual void enableColorPointer() const = 0;
};

template<typename>
struct Type_to_gl_enum;

template<>
struct Type_to_gl_enum<unsigned char> {
  const static GLenum gl = GL_UNSIGNED_BYTE;
};

template<>
struct Type_to_gl_enum<float> {
  const static GLenum gl = GL_FLOAT;
};

template<>
struct Type_to_gl_enum<double> {
  const static GLenum gl = GL_DOUBLE;
};


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

template<typename Word>
struct ImageWrap {
  boost::shared_ptr< std::vector<Word> > img_;
  unsigned int xdim, ydim;

  Word image_data(unsigned int i, unsigned int j, unsigned int k) const {
    return (*img_)[(k * ydim + j) * xdim + i];
  }
};

template<typename Word, typename Tag>
class Volume_plane_impl : public Volume_plane<Tag>  {
public:
  Volume_plane_impl(unsigned int adim, unsigned int bdim, unsigned int cdim,
                    float xscale, float yscale, float zscale, ImageWrap<Word> img) 
    : Volume_plane<Tag>(adim, bdim, cdim, xscale, yscale, zscale), wrap(img)
    {
      buildColors(adim, bdim);
    }
protected:
  void buildColors(unsigned int x, unsigned y) const {
    colors.clear();
      for(unsigned int i = 0; i < (y - 1); ++i)
      {
        for(unsigned int j = 0; j < (x - 1); ++j)
        {
          Word v = image_data(j, i, Volume_plane<Tag>::currentCube, *this);
          for(int k = 0; k < 4; ++k)
          {
            colors.push_back(v);
            colors.push_back(0.0);
            colors.push_back(0.0);
        // v = image_data(j + 1, i, Volume_plane<Tag>::currentCube, *this);
        //   colors.push_back(v);
        //   colors.push_back(v);
        //   colors.push_back(v);
        // v = image_data(j + 1, i + 1, Volume_plane<Tag>::currentCube, *this);
        //   colors.push_back(v);
        //   colors.push_back(v);
        //   colors.push_back(v);
        // v = image_data(j, i + 1, Volume_plane<Tag>::currentCube, *this);
        //   colors.push_back(v);
        //   colors.push_back(v);
        //   colors.push_back(v);
          }
        }
      }
  }

  void enableColorPointer() const {
    glColorPointer(3, Type_to_gl_enum<Word>::gl, 0, &colors[0]);
  }

private:
  ImageWrap<Word> wrap;
  mutable std::vector< Word > colors;

  Word image_data(unsigned int i, unsigned int j, unsigned int k, x_tag) const {
    return wrap.image_data(k, i, j);
  }
  Word image_data(unsigned int i, unsigned int j, unsigned int k, y_tag) const {
    return wrap.image_data(i, k, j);
  }
  Word image_data(unsigned int i, unsigned int j, unsigned int k, z_tag) const {
    return wrap.image_data(i, j, k);
  }
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
        this->addImageToSc<Word>(img, false);
        break;
      }
      case WK_FLOAT:
      {
        typedef double Word; 
        this->addImageToSc<Word>(img, true);
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
  void addImageToSc(const CGAL::Image_3* img, bool clamp) {
    const Word* begin = (const Word*)img->data();
    const Word* end = (const Word*)img->data() + img->size();
    Clamp_to_one_zero_range<Word> clamper = { std::make_pair(*std::max_element(begin, end),
                                                           *std::min_element(begin, end)) };
        
    boost::shared_ptr< std::vector< Word> > clamped(new std::vector<Word>);
    clamped->reserve(img->size());


    if(clamp)
      std::transform(begin, end, std::back_inserter(*clamped), clamper);
    else
      clamped->assign(begin, end);

    assert(clamped->size() == img->size());

    ImageWrap<Word> wrap = { clamped, img->xdim(), img->ydim() };

    sc->addItem(new Volume_plane_impl<Word, x_tag>(img->ydim(), img->zdim(), img->xdim(), 
                                                            img->vx(), img->vy(), img->vz(), wrap));
    sc->addItem(new Volume_plane_impl<Word, y_tag>(img->xdim(), img->zdim(), img->ydim(), 
                                                            img->vx(), img->vy(), img->vz(), wrap));
    sc->addItem(new Volume_plane_impl<Word, z_tag>(img->xdim(), img->ydim(), img->zdim(), 
                                                            img->vx(), img->vy(), img->vz(), wrap));
  }
};


Q_EXPORT_PLUGIN2(Volume_plane_plugin, Volume_plane_plugin)

#include "Volume_planes_plugin.moc"
