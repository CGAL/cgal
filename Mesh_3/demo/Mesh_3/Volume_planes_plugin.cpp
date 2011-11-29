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
#include <boost/scoped_ptr.hpp>

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
  Volume_plane(double adim, double bdim, float xscale, float yscale, float zscale, const CGAL::Image_3* img)
    : adim_(adim), bdim_(bdim), xscale_(xscale), 
      yscale_(yscale), zscale_(zscale), img_(img), currentCube(0),
      mFrame_(new qglviewer::ManipulatedFrame){ 

    buildColors();
    mFrame_->setConstraint(setConstraint(*this));
  }

  Volume_plane* clone() const {
    return NULL;
  }

  bool supportsRenderingMode(RenderingMode) const {
    return true;
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
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(mFrame_->matrix());
    printGlError(__LINE__);

    updateCurrentCube();

    glBegin(GL_QUADS);
    {
      for(int j = 0; j < bdim_ - 1; ++j)
      {
        for(int i = 0; i < adim_ - 1; ++i)
        {
          patch(i, j, *this);
        }
      }
    }
    glEnd();
    printGlError(__LINE__);


    glPopMatrix();
    printGlError(__LINE__);
  }

private:

  const double adim_, bdim_, xscale_, yscale_, zscale_;
  const CGAL::Image_3* img_;
  mutable int currentCube;
  ManipulatedFrame* mFrame_;
  mutable std::vector< std::vector<unsigned char> > colors;
  
  QString name(x_tag) const {return "X Slice"; }
  QString name(y_tag) const {return "Y Slice"; }
  QString name(z_tag) const {return "Z Slice"; }

  void buildColors() const {
    colors.clear();
    for(int i = 0; i < adim_; ++i)
    {
      colors.push_back(std::vector<unsigned char>());
      for(int j = 0; j < bdim_; ++j)
      {
        colors.back().push_back(image_data(i,j, *this));
      }
    }
  }

  qglviewer::Constraint* setConstraint(x_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<XSel>(img_->xdim() * xscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(1.0f, 0.0f, 0.0f));
    return c;
  }
  qglviewer::Constraint* setConstraint(y_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<YSel>(img_->ydim() * yscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(0.0f, 1.0f, 0.0f));
    return c;
  }
  qglviewer::Constraint* setConstraint(z_tag) {
    qglviewer::AxisPlaneConstraint* c = new Length_constraint<ZSel>(img_->zdim() * zscale_);
    c->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    c->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS, qglviewer::Vec(0.0f, 0.0f, 1.0f));
    return c;
  }

  void updateCurrentCube() const {
    int tmp = getTranslation(*this);
    // to which cube does this translation take us?
    if(tmp != currentCube) {
      currentCube = tmp;
      buildColors();
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

  unsigned char image_data(unsigned int a, unsigned int b, x_tag) const {
    return CGAL::IMAGEIO::static_evaluate<unsigned char>(img_->image(),currentCube, a, b);
  }
  unsigned char image_data(unsigned int a, unsigned int b, y_tag) const {
    return CGAL::IMAGEIO::static_evaluate<unsigned char>(img_->image(),a, currentCube, b);
  }
  unsigned char image_data(unsigned int a, unsigned int b, z_tag) const {
    // if(a < adim_ || b < bdim_) { std::cout << a << b << currentCube << std::endl;
    //   }
    return CGAL::IMAGEIO::static_evaluate<unsigned char>(img_->image(),a ,b, currentCube);
  }

  void patch(int i, int j, x_tag) const {
    glNormal3f(1.0f, 0.0f, 0.0f);
    glColor3ub(colors[i][j], colors[i][j], colors[i][j]);
    glVertex3f(0.0f, i * yscale_, j * zscale_);

    glNormal3f(1.0f, 0.0f, 0.0f);
    // glColor3ub(colors[i + 1][j], colors[i + 1][j], colors[i + 1][j]);
    glVertex3f(0.0f, (i + 1) * yscale_, j * zscale_);

    glNormal3f(1.0f, 0.0f, 0.0f);
    // glColor3ub(colors[i + 1][j + j], colors[i + 1][j + j], colors[i + 1][j + j]);
    glVertex3f(0.0f, (i + 1) * yscale_, (j + 1) * zscale_);

    glNormal3f(1.0f, 0.0f, 0.0f);
    // glColor3ub(colors[i][j + j], colors[i][j + j], colors[i][j + j]);
    glVertex3f(0.0f, i * yscale_, (j + 1) * zscale_);
  }

  void patch(int i, int j, y_tag) const {
    glNormal3f(0.0f, 1.0f, 0.0f);
    glColor3ub(colors[i][j], colors[i][j], colors[i][j]);
    glVertex3f(i * xscale_, 0.0f, j * zscale_);

    glNormal3f(0.0f, 1.0f, 0.0f);
    // glColor3ub(colors[i + 1][j], colors[i + 1][j], colors[i + 1][j]);
    glVertex3f((i + 1) * xscale_, 0.0f, j * zscale_);

    glNormal3f(0.0f, 1.0f, 0.0f);
    // glColor3ub(colors[i + 1][j + j], colors[i + 1][j + j], colors[i + 1][j + j]);
    glVertex3f((i + 1) * xscale_, 0.0f, (j + 1) * zscale_);

    glNormal3f(0.0f, 1.0f, 0.0f);
    // glColor3ub(colors[i][j + j], colors[i][j + j], colors[i][j + j]);
    glVertex3f(i * xscale_, 0.0f, (j + 1) * zscale_);
  }

  void patch(int i, int j, z_tag) const {
    glNormal3f(0.0f, 0.0f, 1.0f);
    glColor3ub(colors[i][j], colors[i][j], colors[i][j]);
    glVertex3f(i * xscale_, j * yscale_, 0.0f);

    glNormal3f(0.0f, 0.0f, 1.0f);
    // glColor3ub(colors[i + 1][j], colors[i + 1][j], colors[i + 1][j]);
    glVertex3f((i + 1) * xscale_, j * yscale_, 0.0f);

    glNormal3f(0.0f, 0.0f, 1.0f);
    // glColor3ub(colors[i + 1][j + j], colors[i + 1][j + j], colors[i + 1][j + j]);
    glVertex3f((i + 1) * xscale_, (j + 1) * yscale_, 0.0f);

    glNormal3f(0.0f, 0.0f, 1.0f);
    // glColor3ub(colors[i][j + j], colors[i][j + j], colors[i][j + j]);
    glVertex3f(i * xscale_, (j + 1) * yscale_, 0.0f);
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
      CGAL::Image_3 img = *seg_img->image();
      std::cout << img.xdim() << " " << img.ydim() << " " << img.zdim() << std::endl;


      // x
      sc->addItem(new Volume_plane<x_tag>(img.ydim(), img.zdim(), img.vx(), img.vy(), img.vz(), seg_img->image()));
      sc->addItem(new Volume_plane<y_tag>(img.xdim(), img.zdim(), img.vx(), img.vy(), img.vz(), seg_img->image()));
      sc->addItem(new Volume_plane<z_tag>(img.xdim(), img.ydim(), img.vx(), img.vy(), img.vz(), seg_img->image()));
    } else {
      throw std::runtime_error("No suitable data for Volume planes found");
    }
  }
 
private:
  QAction* planeSwitch;
  Scene_interface* sc;
  QMainWindow* mw;
};


Q_EXPORT_PLUGIN2(Volume_plane_plugin, Volume_plane_plugin)

#include "Volume_planes_plugin.moc"
