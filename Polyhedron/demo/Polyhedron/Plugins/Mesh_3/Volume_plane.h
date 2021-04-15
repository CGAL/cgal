#ifndef CGAL_VOLUME_PLANE_H
#define CGAL_VOLUME_PLANE_H

#include <CGAL/Three/Scene_item_rendering_helper.h>

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
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Qt/constraint.h>
#include <CGAL/Qt/debug.h>

using namespace CGAL::Three;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Viewer_interface Vi;

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
   setBbox(Bbox(tx,ty,tz,
                tx, ty+(adim_ - 1) * yscale_, tz+(bdim_ - 1) * zscale_));
 }

 void compute_bbox(y_tag)const
 {
   setBbox(Bbox(tx,ty,tz,
                tx+(adim_ - 1) * xscale_, ty, tz+(bdim_ - 1) * zscale_));
 }

 void compute_bbox(z_tag)const
 {
   setBbox(Bbox(tx,ty,tz,
                tx+(adim_ - 1) * xscale_,ty+ (bdim_ - 1) * yscale_, tz));
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
  void init(Viewer_interface* viewer)
  {
    if(!isInit(viewer))
      initGL(viewer);
    viewer->makeCurrent();
    getTriangleContainer(1)->reset_vbos(ALL);
    drawSpheres(*this);
    getTriangleContainer(1)->allocate(
          Tc::Flat_vertices,
          v_spheres.data(),
          static_cast<int>(v_spheres.size()*sizeof(float)));
    getTriangleContainer(1)->allocate(
          Tc::Flat_normals,
          n_spheres.data(),
          static_cast<int>(n_spheres.size()*sizeof(float)));
    getTriangleContainer(1)->allocate(
          Tc::Facet_centers,
          c_spheres.data(),
          static_cast<int>(c_spheres.size()*sizeof(float)));

    nb_vertices = v_spheres.size();
    nb_centers = c_spheres.size();
    v_rec.resize(0);
    getEdgeContainer(0)->reset_vbos(ALL);
    drawRectangle(*this, !viewer->isOpenGL_4_3());
    getEdgeContainer(0)->allocate(Ec::Vertices,
                                  v_rec.data(),
                                  static_cast<int>(v_rec.size()*sizeof(float)));
    nb_edges = v_rec.size();
    is_grabbing = false;

    // for each vertex
    getTriangleContainer(0)->reset_vbos(ALL);
    vertices.resize(0);
    indices.resize(0);
    vertices.reserve(bdim_ * adim_ * 3);
    for(unsigned int i = 0; i < adim_; ++i)
    {
      for(unsigned int j = 0; j < bdim_; ++j)
      {
        buildVertex(vertices, i, j);
      }
    }

    assert(vertices.size() == (3 * adim_ * bdim_));

    // for each patch
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
    getTriangleContainer(0)->allocate(
          Tc::Vertex_indices,
          indices.data(),
          static_cast<int>(sizeof(unsigned int) * indices.size()));
    getTriangleContainer(0)->allocate(
          Tc::Smooth_vertices,
          vertices.data(),
          static_cast<int>(sizeof(float) * vertices.size()));
    getTriangleContainer(0)->allocate(
          Tc::VColors,
          colors_.data(),
          static_cast<int>(colors_.size()*sizeof(float)));
    nb_idx = indices.size();
  }

private:
  bool is_grabbing;

  CGAL::qglviewer::Vec translationVector(x_tag) const {
    return CGAL::qglviewer::Vec(xscale_, 0.0, 0.0);
  }
  CGAL::qglviewer::Vec translationVector(y_tag) const {
    return CGAL::qglviewer::Vec(0.0, yscale_, 0.0);
  }
  CGAL::qglviewer::Vec translationVector(z_tag) const {
    return CGAL::qglviewer::Vec(0.0, 0.0, zscale_);
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

  unsigned int adim_, bdim_, cdim_;
  double xscale_, yscale_, zscale_;
  mutable int currentCube;
  mutable std::vector<float> v_rec;
  mutable std::vector<float> v_spheres;
  mutable std::vector<float> n_spheres;
  mutable std::vector<float> c_spheres;
  mutable std::size_t nb_vertices;
  mutable std::size_t nb_centers;
  mutable std::size_t nb_edges;
  mutable std::size_t nb_idx;
  mutable double sphere_radius;
  mutable std::vector< float > colors_;
  mutable std::vector< float > vertices;
  mutable std::vector<unsigned int> indices;


  QString name(x_tag) const { return tr("X Slice for %1").arg(name_); }
  QString name(y_tag) const { return tr("Y Slice for %2").arg(name_); }
  QString name(z_tag) const { return tr("Z Slice for %2").arg(name_); }

  //according to the tag, a,b,c dim change but not the scale. We look for the max dimension of the whole image.
  //A high scale factor will often go with a low dimesion, to compensate it. So we don't want a max being the
  //higher scale * the higher dim, hence the tag specialisation.
//TODO: set the scale factors according to the dimensipon to avoid doing that.
  double compute_maxDim(x_tag) const
  {
    double max_a((adim_ - 1) * yscale_),
        max_b((bdim_ - 1) * zscale_),
        max_c((cdim_ - 1) * xscale_);

    return (std::max)((std::max)(max_a, max_b), max_c);
  }

  double compute_maxDim(y_tag) const
  {
    double max_a((adim_ - 1) * xscale_),
        max_b((bdim_ - 1) * zscale_),
        max_c((cdim_ - 1) * yscale_);

    return (std::max)((std::max)(max_a, max_b), max_c);
  }

  double compute_maxDim(z_tag) const
  {
    double max_a((adim_ - 1) * xscale_),
        max_b((bdim_ - 1) * yscale_),
        max_c((cdim_ - 1) * zscale_);

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
      double max_dim = compute_maxDim(x_tag());
      sphere_radius = max_dim / 40.0f;
      create_flat_sphere(1.0f, v_spheres, n_spheres, 18);

      c_spheres.push_back(0.0f); c_spheres.push_back((adim_ - 1) * yscale_/2.0f + 1.1*sphere_radius); c_spheres.push_back(0.0f);

      c_spheres.push_back(0.0f); c_spheres.push_back((adim_ - 1) * yscale_ ); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f + 1.1*sphere_radius);

      c_spheres.push_back(0.0f); c_spheres.push_back((adim_ - 1) * yscale_/2.0f + 1.1*sphere_radius); c_spheres.push_back((bdim_ - 1 ) * zscale_);

      c_spheres.push_back(0.0f); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f + 1.1*sphere_radius);

  }

  void drawSpheres(y_tag) const
  {
      double max_dim = compute_maxDim(y_tag());
      sphere_radius = max_dim / 40.0f;
      create_flat_sphere(1.0f, v_spheres, n_spheres,18);

      c_spheres.push_back((adim_ - 1) * xscale_/2.0f); c_spheres.push_back(0.0f); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f);
      c_spheres.push_back((adim_ - 1) * xscale_/2.0f); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_);
      c_spheres.push_back(0.0f); c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * zscale_/2.0f);
  }

  void drawSpheres(z_tag) const
  {
      double max_dim = compute_maxDim(z_tag());
      sphere_radius = max_dim / 40.0f;
      create_flat_sphere(1.0f, v_spheres, n_spheres,18);

      c_spheres.push_back(0.0f); c_spheres.push_back((bdim_ - 1) * yscale_/2.0f - 1.1*sphere_radius); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_/2.0f-1.1*sphere_radius); c_spheres.push_back((bdim_ - 1) * yscale_); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_); c_spheres.push_back((bdim_ - 1) * yscale_/2.0f-1.1*sphere_radius); c_spheres.push_back(0.0f);
      c_spheres.push_back((adim_ - 1) * xscale_/2.0f-1.1*sphere_radius); c_spheres.push_back(0.0f); c_spheres.push_back(0.0f);
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
    getTriangleContainer(0)->getVbo(Tc::VColors)->offset
        = (currentCube*int(sizeof(float))) * (bdim_) * (adim_) * 3;
      getTriangleContainer(1)->initializeBuffers(viewer);
      getTriangleContainer(1)->setFlatDataSize(nb_vertices);
      getTriangleContainer(1)->setCenterSize(nb_centers);
      getTriangleContainer(0)->initializeBuffers(viewer);
      getTriangleContainer(0)->setIdxSize(nb_idx);
      getEdgeContainer(0)->initializeBuffers(viewer);
      getEdgeContainer(0)->setFlatDataSize(nb_edges);
      v_spheres.clear();
      v_spheres.shrink_to_fit();
      n_spheres.clear();
      n_spheres.shrink_to_fit();
      //not the centers, they are used for picking
      v_rec.clear();
      v_rec.shrink_to_fit();
      vertices.clear();
      vertices.shrink_to_fit();
      indices.clear();
      indices.shrink_to_fit();
  }
};

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
    setTriangleContainer(1, new Tc(Vi::PROGRAM_SPHERES, false));
    setTriangleContainer(0, new Tc(Vi::PROGRAM_NO_SELECTION, true));
    setEdgeContainer(0, new Ec(Three::mainViewer()->isOpenGL_4_3()
                               ? PROGRAM_SOLID_WIREFRAME
                               : PROGRAM_NO_SELECTION, false));

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
}
template<typename T>
Volume_plane<T>::~Volume_plane() {
  delete sphere_Slider;
}

template<typename T>
void Volume_plane<T>::draw(Viewer_interface *viewer) const {
  double cur_cube = currentCube;
  updateCurrentCube();
  if(cur_cube != currentCube)
  {
    Q_FOREACH(CGAL::QGLViewer*v, CGAL::QGLViewer::QGLViewerPool()){
      v->setProperty("need_update", true);
    }
  }
  if(viewer->property("need_update").toBool())
  {
    setBuffersInit(viewer, false);
    viewer->setProperty("need_update", false);
  }
  if(!isInit(viewer))
    initGL(viewer);
  {
    setBuffersInit(viewer, false);
  }
  if ( ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }


  viewer->glDepthRangef(0.00001f,0.99999f);
  GLdouble mat[16];
  viewer->camera()->getModelViewProjectionMatrix(mat);
  QMatrix4x4 f;
  for(int i=0; i<16; i++)
  {
    f.data()[i] = (float)mFrame_->matrix()[i];
  }

  f.translate(QVector3D(tx, ty, tz));
  getEdgeContainer(0)->setFrameMatrix(f);
  printGlError(viewer, __LINE__);
  QVector2D vp(viewer->width(), viewer->height());
  getEdgeContainer(0)->setViewport(vp);
  getEdgeContainer(0)->setWidth(4.0f);
  getEdgeContainer(0)->setColor(QColor(Qt::black));
  viewer->glDepthRangef(0.00005f, 0.99995f);
  getEdgeContainer(0)->draw(viewer, true);
  viewer->glDepthRangef(0.0,1.0);


  getTriangleContainer(0)->setFrameMatrix(f);
  getTriangleContainer(0)->draw(viewer, false);
  printGlError(viewer, __LINE__);

  //hide spheres if only 1 plane.
  if(aDim() <= 1 ||
     bDim() <= 1 ||
     cDim() <=1)
    return;
  double max_dim = compute_maxDim(*this);
  sphere_radius = max_dim/20.0f * sphere_Slider->value()/100.0f;
  getTriangleContainer(1)->getVbo(Tc::Radius)->bind();
  getTriangleContainer(1)->getVao(viewer)->program->setAttributeValue("radius", sphere_radius);
  getTriangleContainer(1)->getVbo(Tc::Radius)->release();
  getTriangleContainer(1)->setColor(this->color());
  getTriangleContainer(1)->setFrameMatrix(f);
  getTriangleContainer(1)->draw(viewer, true);
}



//intercept events for picking
template<typename T>
bool Volume_plane<T>::eventFilter(QObject *sender, QEvent *event)
{
    CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(sender);
    if(!viewer)
      return false;
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
        redraw();
    }
    return false;
}

#endif /* CGAL_VOLUME_PLANE_H */
