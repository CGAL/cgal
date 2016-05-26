#include "Scene_plane_item.h"
using namespace CGAL::Three;


struct Scene_plane_item_priv
{
  Scene_plane_item_priv(const CGAL::Three::Scene_interface* scene_interface, Scene_plane_item* parent)
      :scene(scene_interface),
      manipulable(false),
      can_clone(true),
      frame(new Scene_item::ManipulatedFrame())
  {
    item = parent;
  }

  ~Scene_plane_item_priv() {
    frame = 0;
    delete frame;
  }

  double scene_diag() const {
    const Scene_item::Bbox& bbox = scene->bbox();
    const double& xdelta = bbox.xmax()-bbox.xmin();
    const double& ydelta = bbox.ymax()-bbox.ymin();
    const double& zdelta = bbox.zmax()-bbox.zmin();
    const double diag = std::sqrt(xdelta*xdelta +
                            ydelta*ydelta +
                            zdelta*zdelta);
    return diag * 0.7;
  }

  void initializeBuffers(Viewer_interface*)const;
  void compute_normals_and_vertices(void);

  enum VAOs {
      Facets = 0,
      Edges,
      NbOfVaos
  };
  enum VBOs {
      Facets_vertices = 0,
      Edges_vertices,
      NbOfVbos
  };

  const CGAL::Three::Scene_interface* scene;
  bool manipulable;
  bool can_clone;
  qglviewer::ManipulatedFrame* frame;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_quad;
  mutable GLint sampler_location;
  mutable bool smooth_shading;
  mutable QOpenGLShaderProgram *program;
  mutable bool are_buffers_filled;
  Scene_plane_item* item;

};


Scene_plane_item::Scene_plane_item(const CGAL::Three::Scene_interface* scene_interface)
    :CGAL::Three::Scene_item(Scene_plane_item_priv::NbOfVbos,Scene_plane_item_priv::NbOfVaos)
{
  setNormal(0., 0., 1.);
  d = new Scene_plane_item_priv(scene_interface, this);
  invalidateOpenGLBuffers();
}

Scene_plane_item::~Scene_plane_item() {
  delete d;
}

void Scene_plane_item_priv::initializeBuffers(Viewer_interface *viewer) const
{
    program = item->getShaderProgram(Scene_plane_item::PROGRAM_WITHOUT_LIGHT, viewer);
    program->bind();
    item->vaos[Facets]->bind();

    item->buffers[Facets_vertices].bind();
    item->buffers[Facets_vertices].allocate(positions_quad.data(),
                        static_cast<int>(positions_quad.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    item->buffers[Facets_vertices].release();
    item->vaos[Facets]->release();


    item->vaos[Edges]->bind();
    item->buffers[Edges_vertices].bind();
    item->buffers[Edges_vertices].allocate(positions_lines.data(),
                        static_cast<int>(positions_lines.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    item->buffers[Edges_vertices].release();
    item->vaos[Edges]->release();

    program->release();
    are_buffers_filled = true;

}

void Scene_plane_item_priv::compute_normals_and_vertices(void)
{
    positions_quad.resize(0);
    positions_lines.resize(0);

    const double diag = scene_diag();
    //The quad
    {

    positions_quad.push_back(-diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(-diag);
    positions_quad.push_back(diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);

    positions_quad.push_back(-diag);
    positions_quad.push_back(diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(diag);
    positions_quad.push_back(diag);
    positions_quad.push_back(0.0);

}
    //The grid
    float x = (2*diag)/10.0;
    float y = (2*diag)/10.0;
    {
        for(int u = 0; u < 11; u++)
        {

            positions_lines.push_back(-diag + x* u);
            positions_lines.push_back(-diag);
            positions_lines.push_back(0.0);

            positions_lines.push_back(-diag + x* u);
            positions_lines.push_back(diag);
            positions_lines.push_back(0.0);
        }
        for(int v=0; v<11; v++)
        {

            positions_lines.push_back(-diag);
            positions_lines.push_back(-diag + v * y);
            positions_lines.push_back(0.0);

            positions_lines.push_back(diag);
            positions_lines.push_back(-diag + v * y);
            positions_lines.push_back(0.0);
        }

    }
}

void Scene_plane_item::draw(Viewer_interface* viewer)const
{
    if(!d->are_buffers_filled)
        d->initializeBuffers(viewer);
    vaos[Scene_plane_item_priv::Facets]->bind();
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    QMatrix4x4 f_matrix;
    for(int i=0; i<16; i++)
        f_matrix.data()[i] = (float)d->frame->matrix()[i];
    d->program->bind();
    d->program->setUniformValue("f_matrix", f_matrix);
    d->program->setAttributeValue("colors",this->color());
    d->program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->positions_quad.size()/3));
    d->program->release();
    vaos[Scene_plane_item_priv::Facets]->release();

}

void Scene_plane_item::drawEdges(CGAL::Three::Viewer_interface* viewer)const
{
    if(!d->are_buffers_filled)
        d->initializeBuffers(viewer);
    vaos[Scene_plane_item_priv::Edges]->bind();
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    QMatrix4x4 f_matrix;
    for(int i=0; i<16; i++)
        f_matrix.data()[i] = (float)d->frame->matrix()[i];
    d->program->bind();
    d->program->setUniformValue("f_matrix", f_matrix);
    d->program->setAttributeValue("colors",QVector3D(0,0,0));
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->positions_lines.size()/3));
    d->program->release();
    vaos[Scene_plane_item_priv::Edges]->release();
}

void Scene_plane_item::flipPlane()
{
  qglviewer::Quaternion q;
  qglviewer::Vec axis(0,1,0);
  if(d->frame->orientation().axis() == axis)
    q.setAxisAngle(qglviewer::Vec(1,0,0), M_PI);
  else
    q.setAxisAngle(axis, M_PI);
  d->frame->rotate(q.normalized());
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();

}

bool Scene_plane_item::manipulatable() const {
  return d->manipulable;
}
Scene_item::ManipulatedFrame* Scene_plane_item::manipulatedFrame() {
  return d->frame;
}

Scene_plane_item* Scene_plane_item::clone() const {
  if(d->can_clone)
  {
    Scene_plane_item* item = new Scene_plane_item(d->scene);
    item->d->manipulable = d->manipulable;
    item->d->can_clone = true;
    item->d->frame = new ManipulatedFrame;
    item->d->frame->setPosition(d->frame->position());
    item->d->frame->setOrientation(d->frame->orientation());
    return item;
  }
  else
    return 0;
}

QString Scene_plane_item::toolTip() const {
  const qglviewer::Vec& pos = d->frame->position();
  const qglviewer::Vec& n = d->frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return
    tr("<p><b>%1</b> (mode: %2, color: %3)<br />")
    .arg(this->name())
    .arg(this->renderingModeName())
    .arg(this->color().name())

    +
    tr("<i>Plane</i></p>"
       "<p>Equation: %1*x + %2*y + %3*z + %4 = 0<br />"
       "Normal vector: (%1, %2, %3)<br />"
       "Point: (%5, %6, %7)</p>")
    .arg(n[0]).arg(n[1]).arg(n[2])
    .arg( - pos * n)
    .arg(pos[0]).arg(pos[1]).arg(pos[2])
    +
    tr("<p>Can clone: %1<br />"
       "Manipulatable: %2</p>")
    .arg(d->can_clone?tr("true"):tr("false"))
    .arg(d->manipulable?tr("true"):tr("false"));
}

Plane_3 Scene_plane_item::plane() const {
  const qglviewer::Vec& pos = d->frame->position();
  const qglviewer::Vec& n =
    d->frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return Plane_3(n[0], n[1],  n[2], - n * pos);
}

void Scene_plane_item::invalidateOpenGLBuffers()
{
    d->compute_normals_and_vertices();
    d->are_buffers_filled = false;
    compute_bbox();
}

void Scene_plane_item::setPosition(float x, float y, float z) {
  d->frame->setPosition(x, y, z);
}

void Scene_plane_item::setPosition(double x, double y, double z) {
  d->frame->setPosition((float)x, (float)y, (float)z);
}

void Scene_plane_item::setNormal(float x, float y, float z) {
  QVector3D normal(x,y,z);
  if(normal == QVector3D(0,0,0))
    return;
  QVector3D origin(0,0,1);
  qglviewer::Quaternion q;
  if(origin == normal)
  {
    return;
  }
   if(origin == -normal)
  {
    q.setAxisAngle(qglviewer::Vec(0,1,0),M_PI);
    d->frame->setOrientation(q);
    return;
  }

  QVector3D cp = QVector3D::crossProduct(origin, normal);
  cp.normalize();
  q.setAxisAngle(qglviewer::Vec(cp.x(),cp.y(), cp.z()),acos(QVector3D::dotProduct(origin, normal)/(normal.length()*origin.length())));

  d->frame->setOrientation(q.normalized());
}

void Scene_plane_item::setNormal(double x, double y, double z) {
  setNormal((float)x, (float)y, (float)z);
}

void Scene_plane_item::setClonable(bool b) {
  d->can_clone = b;
}

void Scene_plane_item::setManipulatable(bool b) {
  d->manipulable = b;
}
