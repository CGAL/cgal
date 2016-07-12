#include "Scene_plane_item.h"
#include <QApplication>
using namespace CGAL::Three;



Scene_plane_item::Scene_plane_item(const CGAL::Three::Scene_interface* scene_interface)
      :CGAL::Three::Scene_item(NbOfVbos,NbOfVaos),
      scene(scene_interface),
      manipulable(false),
      can_clone(true),
      frame(new ManipulatedFrame())
  {
    setNormal(0., 0., 1.);
    //Generates an integer which will be used as ID for each buffer
    invalidateOpenGLBuffers();
  }
Scene_plane_item::~Scene_plane_item() {
  delete frame;
}

void Scene_plane_item::initializeBuffers(Viewer_interface *viewer) const
{
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
    program->bind();
    vaos[Facets]->bind();

    buffers[Facets_vertices].bind();
    buffers[Facets_vertices].allocate(positions_quad.data(),
                        static_cast<int>(positions_quad.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[Facets_vertices].release();
    vaos[Facets]->release();


    vaos[Edges]->bind();
    buffers[Edges_vertices].bind();
    buffers[Edges_vertices].allocate(positions_lines.data(),
                        static_cast<int>(positions_lines.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[Edges_vertices].release();
    vaos[Edges]->release();

    program->release();
    are_buffers_filled = true;

}

void Scene_plane_item::compute_normals_and_vertices(void) const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
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

    positions_quad.push_back(diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(-diag);
    positions_quad.push_back(diag);
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
    QApplication::restoreOverrideCursor();
}

void Scene_plane_item::draw(Viewer_interface* viewer)const
{
    if(!are_buffers_filled)
        initializeBuffers(viewer);
    vaos[Facets]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    QMatrix4x4 f_matrix;
    for(int i=0; i<16; i++)
        f_matrix.data()[i] = (float)frame->matrix()[i];
    program->bind();
    program->setUniformValue("f_matrix", f_matrix);
    program->setAttributeValue("colors",this->color());
    program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(positions_quad.size()/3));
    program->release();
    vaos[Facets]->release();

}

void Scene_plane_item::drawEdges(CGAL::Three::Viewer_interface* viewer)const
{
    if(!are_buffers_filled)
        initializeBuffers(viewer);
    vaos[Edges]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    QMatrix4x4 f_matrix;
    for(int i=0; i<16; i++)
        f_matrix.data()[i] = (float)frame->matrix()[i];
    program->bind();
    program->setUniformValue("f_matrix", f_matrix);
    program->setAttributeValue("colors",QVector3D(0,0,0));
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
    program->release();
    vaos[Edges]->release();
}

void Scene_plane_item::flipPlane()
{
  qglviewer::Quaternion q;
  qglviewer::Vec axis(0,1,0);
  if(frame->orientation().axis() == axis)
    q.setAxisAngle(qglviewer::Vec(1,0,0), M_PI);
  else
    q.setAxisAngle(axis, M_PI);
  frame->rotate(q.normalized());
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();

}

bool Scene_plane_item::manipulatable() const {
  return manipulable;
}
Scene_item::ManipulatedFrame* Scene_plane_item::manipulatedFrame() {
  return frame;
}

Scene_plane_item* Scene_plane_item::clone() const {
  if(can_clone)
  {
    Scene_plane_item* item = new Scene_plane_item(scene);
    item->manipulable = manipulable;
    item->can_clone = true;
    item->frame = new ManipulatedFrame;
    item->frame->setPosition(frame->position());
    item->frame->setOrientation(frame->orientation());
    return item;
  }
  else
    return 0;
}

QString Scene_plane_item::toolTip() const {
  const qglviewer::Vec& pos = frame->position();
  const qglviewer::Vec& n = frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
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
    .arg(can_clone?tr("true"):tr("false"))
    .arg(manipulable?tr("true"):tr("false"));
}

Plane_3 Scene_plane_item::plane() const {
  const qglviewer::Vec& pos = frame->position();
  const qglviewer::Vec& n =
    frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return Plane_3(n[0], n[1],  n[2], - n * pos);
}

void Scene_plane_item::invalidateOpenGLBuffers()
{
    compute_normals_and_vertices();
    are_buffers_filled = false;
    compute_bbox();
}

void Scene_plane_item::setPosition(float x, float y, float z) {
  frame->setPosition(x, y, z);
}

void Scene_plane_item::setPosition(double x, double y, double z) {
  frame->setPosition((float)x, (float)y, (float)z);
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
    frame->setOrientation(q);
    return;
  }

  QVector3D cp = QVector3D::crossProduct(origin, normal);
  cp.normalize();
  q.setAxisAngle(qglviewer::Vec(cp.x(),cp.y(), cp.z()),acos(QVector3D::dotProduct(origin, normal)/(normal.length()*origin.length())));

  frame->setOrientation(q.normalized());
}

void Scene_plane_item::setNormal(double x, double y, double z) {
  setNormal((float)x, (float)y, (float)z);
}

void Scene_plane_item::setClonable(bool b) {
  can_clone = b;
}

void Scene_plane_item::setManipulatable(bool b) {
  manipulable = b;
}
