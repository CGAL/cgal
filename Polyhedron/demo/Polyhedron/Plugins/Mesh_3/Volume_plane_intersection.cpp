#include "config.h"

#include "Volume_plane_intersection.h"
#include "Volume_plane_interface.h"

#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>
#include <QApplication>

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Edge_container Ec;
struct Volume_plane_intersection_priv
{
  Volume_plane_intersection_priv(float x, float y, float z,
                                 float tx, float ty, float tz,
                                 Volume_plane_intersection* parent)
    : a(NULL), b(NULL), c(NULL), x(x), y(y), z(z),
      tx(tx), ty(ty), tz(tz), item(parent)
  {
    item->are_buffers_filled = false;
  }
  Volume_plane_interface *a, *b, *c;
  float x, y, z, tx, ty, tz;
  mutable std::vector<float> a_vertex;
  mutable std::vector<float> b_vertex;
  mutable std::vector<float> c_vertex;
  mutable std::size_t a_size;
  mutable std::size_t b_size;
  mutable std::size_t c_size;
  mutable QOpenGLShaderProgram *program;
  void computeElements();
  void initializeBuffers(Viewer_interface*)const;
  void attribBuffers(Viewer_interface*) const;
  Volume_plane_intersection* item;
};


Volume_plane_intersection::~Volume_plane_intersection()
{
  delete d;
}

void Volume_plane_intersection_priv::computeElements()
{
   QApplication::setOverrideCursor(Qt::WaitCursor);
   a_vertex.resize(0);
   b_vertex.resize(0);
   c_vertex.resize(0);
   const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();

   a_vertex.push_back(0.0-offset.x+tx);
   a_vertex.push_back(0.0-offset.y+ty);
   a_vertex.push_back(0.0-offset.z+tz);
   a_vertex.push_back(x  -offset.x+tx);
   a_vertex.push_back(0.0-offset.y+ty);
   a_vertex.push_back(0.0-offset.z+tz);

   b_vertex.push_back(0.0-offset.x+tx);
   b_vertex.push_back(0.0-offset.y+ty);
   b_vertex.push_back(0.0-offset.z+tz);
   b_vertex.push_back(0.0-offset.x+tx);
   b_vertex.push_back(y  -offset.y+ty);
   b_vertex.push_back(0.0-offset.z+tz);

   c_vertex.push_back(0.0-offset.x+tx);
   c_vertex.push_back(0.0-offset.y+ty);
   c_vertex.push_back(0.0-offset.z+tz);
   c_vertex.push_back(0.0-offset.x+tx);
   c_vertex.push_back(0.0-offset.y+ty);
   c_vertex.push_back(z  -offset.z+tz);

   item->_bbox =  Scene_item::Bbox( tx, ty, tz,
                  tz+x, ty+y, tz+z);

   item->getEdgeContainer(0)->allocate(Ec::Vertices,
                                       a_vertex.data(),
                                       static_cast<int>(a_vertex.size()*sizeof(float)));
   item->getEdgeContainer(1)->allocate(Ec::Vertices,
                                       b_vertex.data(),
                                       static_cast<int>(b_vertex.size()*sizeof(float)));
   item->getEdgeContainer(2)->allocate(Ec::Vertices,
                                       c_vertex.data(),
                                       static_cast<int>(c_vertex.size()*sizeof(float)));
   a_size = a_vertex.size();
   b_size = b_vertex.size();
   c_size = c_vertex.size();
   QApplication::restoreOverrideCursor();
}

void Volume_plane_intersection_priv::initializeBuffers(Viewer_interface* viewer)const
{

  item->getEdgeContainer(0)->initializeBuffers(viewer);
  item->getEdgeContainer(0)->setFlatDataSize(a_size);
  item->getEdgeContainer(1)->initializeBuffers(viewer);
  item->getEdgeContainer(1)->setFlatDataSize(b_size);
  item->getEdgeContainer(2)->initializeBuffers(viewer);
  item->getEdgeContainer(2)->setFlatDataSize(c_size);

  a_vertex.clear();
  a_vertex.shrink_to_fit();
  b_vertex.clear();
  b_vertex.shrink_to_fit();
  c_vertex.clear();
  c_vertex.shrink_to_fit();
}

void Volume_plane_intersection::draw(Viewer_interface* viewer) const {
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  viewer->glDepthRangef(0.00001f, 0.99999f);

  QVector2D vp(viewer->width(), viewer->height());
  if(d->b && d->c) {
    GLdouble mat[16];
    d->b->manipulatedFrame()->getMatrix(mat);
    QMatrix4x4 b_mat, c_mat;
    for(int i=0; i<16; i++)
    {
       b_mat.data()[i] = (float)mat[i];
    }
    d->c->manipulatedFrame()->getMatrix(mat);
    for(int i=0; i<16; i++)
    {
       c_mat.data()[i] = (float)mat[i];
    }
    getEdgeContainer(0)->setFrameMatrix(b_mat*c_mat);
    getEdgeContainer(0)->setColor(this->color());
    if(viewer->isOpenGL_4_3())
    {
      getEdgeContainer(0)->setViewport(vp);
      getEdgeContainer(0)->setWidth(4.0f);
    }
    getEdgeContainer(0)->draw(viewer, true);
  }

  if(d->a && d->c) {
      GLdouble mat[16];
      d->a->manipulatedFrame()->getMatrix(mat);
      QMatrix4x4 a_mat, c_mat;
      for(int i=0; i<16; i++)
      {
         a_mat.data()[i] = (float)mat[i];
      }
      d->c->manipulatedFrame()->getMatrix(mat);
      for(int i=0; i<16; i++)
      {
         c_mat.data()[i] = (float)mat[i];
      }
      getEdgeContainer(1)->setFrameMatrix(a_mat*c_mat);
      getEdgeContainer(1)->setColor(this->color());
      if(viewer->isOpenGL_4_3())
      {
        getEdgeContainer(1)->setViewport(vp);
        getEdgeContainer(1)->setWidth(4.0f);
      }
      getEdgeContainer(1)->draw(viewer, true);
  }

  if(d->a && d->b) {
      GLdouble mat[16];
      d->a->manipulatedFrame()->getMatrix(mat);
      QMatrix4x4 a_mat, b_mat;
      for(int i=0; i<16; i++)
      {
         a_mat.data()[i] = (float)mat[i];
      }
      d->b->manipulatedFrame()->getMatrix(mat);
      for(int i=0; i<16; i++)
      {
         b_mat.data()[i] = (float)mat[i];
      }
      getEdgeContainer(2)->setFrameMatrix(a_mat*b_mat);
      getEdgeContainer(2)->setColor(this->color());
      if(viewer->isOpenGL_4_3())
      {
        getEdgeContainer(2)->setViewport(vp);
        getEdgeContainer(2)->setWidth(4.0f);
      }
      getEdgeContainer(2)->draw(viewer, true);
  }
  if(!viewer->isOpenGL_4_3())
    viewer->glDepthRangef(0.0f,1.0f);
}

Volume_plane_intersection::Volume_plane_intersection(float x, float y, float z,
                                                     float tx, float ty, float tz)
  :d(new Volume_plane_intersection_priv(x,y,z,tx,ty,tz,this))
{
  setColor(QColor(255, 128, 0));
  setName("Volume plane intersection");
  setEdgeContainer(2, new Ec(
                     Three::mainViewer()->isOpenGL_4_3()
                     ? Vi::PROGRAM_SOLID_WIREFRAME
                     : Vi::PROGRAM_NO_SELECTION,
                     false));
  setEdgeContainer(1, new Ec(
                     Three::mainViewer()->isOpenGL_4_3()
                     ? Vi::PROGRAM_SOLID_WIREFRAME
                     : Vi::PROGRAM_NO_SELECTION,
                     false));
  setEdgeContainer(0, new Ec(
                     Three::mainViewer()->isOpenGL_4_3()
                     ? Vi::PROGRAM_SOLID_WIREFRAME
                     : Vi::PROGRAM_NO_SELECTION,
                     false));
}

void Volume_plane_intersection::setX(Volume_plane_interface* x) { d->a = x; }
void Volume_plane_intersection::setY(Volume_plane_interface* x) { d->b = x; }
void Volume_plane_intersection::setZ(Volume_plane_interface* x) { d->c = x; }
void Volume_plane_intersection::planeRemoved(Volume_plane_interface* i) {
  if(d->a == i) {
    d->a = NULL;
  } else if(d->b == i) {
    d->b = NULL;
  } else if(d->c == i) {
    d->c = NULL;
  }
}

void Volume_plane_intersection::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);

  for(int i=0; i<3; ++i)
    getEdgeContainer(i)->reset_vbos(ALL);
}

void Volume_plane_intersection::initializeBuffers(Viewer_interface *v) const
{
  d->initializeBuffers(v);
}

void Volume_plane_intersection::computeElements() const
{
  d->computeElements();
  setBuffersFilled(true);
}
