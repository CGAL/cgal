#include "config.h"

#include "Volume_plane_intersection.h"
#include "Volume_plane_interface.h"

#include <QApplication>
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

  enum VAOs{
    AArray = 0,
    BArray,
    CArray,
    NumberOfVaos
  };
  enum VBOs{
    AVertex = 0,
    BVertex,
    CVertex,
    NumberOfVbos
  };
  std::vector<float> a_vertex;
  std::vector<float> b_vertex;
  std::vector<float> c_vertex;
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
   QApplication::restoreOverrideCursor();
}

void Volume_plane_intersection_priv::initializeBuffers(Viewer_interface* viewer)const
{
  if(!viewer->isOpenGL_4_3())
    program = item->getShaderProgram(Volume_plane_intersection::PROGRAM_NO_SELECTION, viewer);
  else
    program = item->getShaderProgram(Volume_plane_intersection::PROGRAM_SOLID_WIREFRAME, viewer);
  program->bind();
  item->vaos[AArray]->bind();
  item->buffers[AVertex].bind();
  item->buffers[AVertex].allocate(a_vertex.data(), static_cast<int>(a_vertex.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  item->buffers[AVertex].release();
  item->vaos[AArray]->release();

  item->vaos[BArray]->bind();
  item->buffers[BVertex].bind();
  item->buffers[BVertex].allocate(b_vertex.data(), static_cast<int>(b_vertex.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  item->buffers[BVertex].release();
  item->vaos[BArray]->release();

  item->vaos[CArray]->bind();
  item->buffers[CVertex].bind();
  item->buffers[CVertex].allocate(c_vertex.data(), static_cast<int>(c_vertex.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  item->buffers[CVertex].release();
  item->vaos[CArray]->release();

  program->release();
  item->are_buffers_filled = true;

}

void Volume_plane_intersection::draw(Viewer_interface* viewer) const {
  if(!are_buffers_filled)
  {
    d->computeElements();
    d->initializeBuffers(viewer);
  }
  if(!viewer->isOpenGL_4_3())
  {
    attribBuffers(viewer, PROGRAM_NO_SELECTION);
  }
  else
  {
    attribBuffers(viewer, PROGRAM_SOLID_WIREFRAME);
  }     
  viewer->glDepthRangef(0.00001f, 0.99999f);
  
  QVector2D vp(viewer->width(), viewer->height());
  if(d->b && d->c) {

    vaos[Volume_plane_intersection_priv::AArray]->bind();
    d->program->bind();
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
    d->program->setUniformValue("f_matrix", b_mat*c_mat);
    d->program->setAttributeValue("colors", this->color());
    if(viewer->isOpenGL_4_3())
    {
      d->program->setUniformValue("viewport", vp);
      d->program->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
      d->program->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
      d->program->setUniformValue("width", 4.0f);
    }
    viewer->glDrawArrays(GL_LINES, 0, 2);
    d->program->release();
    vaos[Volume_plane_intersection_priv::AArray]->release();
  }

  if(d->a && d->c) {
      vaos[Volume_plane_intersection_priv::BArray]->bind();
      d->program->bind();
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
      d->program->setUniformValue("f_matrix", a_mat*c_mat);
      d->program->setAttributeValue("colors", this->color());
      if(viewer->isOpenGL_4_3())
      {
        d->program->setUniformValue("viewport", vp);
        d->program->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
        d->program->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
        d->program->setUniformValue("width", 4.0f);
      }
      viewer->glDrawArrays(GL_LINES, 0, 2);
      d->program->release();
      vaos[Volume_plane_intersection_priv::BArray]->release();
  }

  if(d->a && d->b) {
      vaos[Volume_plane_intersection_priv::CArray]->bind();
      d->program->bind();
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
      d->program->setUniformValue("f_matrix", a_mat*b_mat);
      d->program->setAttributeValue("colors", this->color());
      if(viewer->isOpenGL_4_3())
      {
        d->program->setUniformValue("viewport", vp);
        d->program->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
        d->program->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
        d->program->setUniformValue("width", 4.0f);
      }
      viewer->glDrawArrays(GL_LINES, 0, 2);
      d->program->release();
      vaos[Volume_plane_intersection_priv::CArray]->release();
  }
  if(!viewer->isOpenGL_4_3())
  viewer->glDepthRangef(0.0f,1.0f);
}

Volume_plane_intersection::Volume_plane_intersection(float x, float y, float z,
                                                     float tx, float ty, float tz)
  :Scene_item(Volume_plane_intersection_priv::NumberOfVbos, Volume_plane_intersection_priv::NumberOfVaos),
    d(new Volume_plane_intersection_priv(x,y,z,tx,ty,tz,this))
{
  setColor(QColor(255, 128, 0));
  setName("Volume plane intersection");
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
 are_buffers_filled = false;
}
