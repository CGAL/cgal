#include "Scene_spheres_item.h"
#include <QApplication>

struct Scene_spheres_item_priv
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef std::pair<CGAL::Sphere_3<Kernel>*, CGAL::Color> Sphere_pair ;

  Scene_spheres_item_priv(bool planed, Scene_spheres_item* parent)
    :precision(36)
    ,has_plane(planed)

  {
    item = parent;
    create_flat_and_wire_sphere(1.0f,vertices,normals, edges);
  }

  ~Scene_spheres_item_priv() {
    Q_FOREACH(Sphere_pair sphere, spheres)
      delete sphere.first;
  }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const;
  enum Vbos
  {
    Vertices = 0,
    Edge_vertices,
    Normals,
    Center,
    Radius,
    Color,
    Edge_color,
    NbOfVbos
  };
  enum Vaos
  {
    Facets = 0,
    Edges,
    NbOfVaos
  };


  int precision;
  mutable CGAL::Plane_3<Kernel> plane;
  bool has_plane;

  QList<Sphere_pair> spheres;
  mutable std::vector<float> vertices;
  mutable std::vector<float> normals;
  mutable std::vector<float> edges;
  mutable std::vector<float> colors;
  mutable std::vector<float> edges_colors;
  mutable std::vector<float> centers;
  mutable std::vector<float> radius;
  mutable QOpenGLShaderProgram *program;
  mutable int nb_centers;
  Scene_spheres_item* item;

};
Scene_spheres_item::Scene_spheres_item(Scene_group_item* parent, bool planed)
  :CGAL::Three::Scene_item(Scene_spheres_item_priv::NbOfVbos,Scene_spheres_item_priv::NbOfVaos)

{
  setParent(parent);
  d = new Scene_spheres_item_priv(planed, this);
}

Scene_spheres_item::~Scene_spheres_item()
{
  delete d;
}
void Scene_spheres_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  d->colors.clear();
  d->edges_colors.clear();
  d->centers.clear();
  d->radius.clear();
  Q_FOREACH(Sphere_pair sp, d->spheres)
  {
    d->colors.push_back((float)sp.second.red()/255);
    d->colors.push_back((float)sp.second.green()/255);
    d->colors.push_back((float)sp.second.blue()/255);

    d->edges_colors.push_back((float)sp.second.red()/255);
    d->edges_colors.push_back((float)sp.second.green()/255);
    d->edges_colors.push_back((float)sp.second.blue()/255);

    d->centers.push_back(sp.first->center().x());
    d->centers.push_back(sp.first->center().y());
    d->centers.push_back(sp.first->center().z());

    d->radius.push_back(sp.first->squared_radius());

  }
  QApplication::restoreOverrideCursor();
}

void Scene_spheres_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  if(has_plane)
  {
    program = item->getShaderProgram(Scene_spheres_item::PROGRAM_CUTPLANE_SPHERES, viewer);
    item->attribBuffers(viewer, Scene_spheres_item::PROGRAM_CUTPLANE_SPHERES);
  }
  else
  {
    program = item->getShaderProgram(Scene_spheres_item::PROGRAM_SPHERES, viewer);
    item->attribBuffers(viewer, Scene_spheres_item::PROGRAM_SPHERES);
  }

  program->bind();
  item->vaos[Facets]->bind();
  item->buffers[Vertices].bind();
  item->buffers[Vertices].allocate(vertices.data(),
                             static_cast<int>(vertices.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
  item->buffers[Vertices].release();

  item->buffers[Normals].bind();
  item->buffers[Normals].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(float)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
  item->buffers[Normals].release();

  item->buffers[Color].bind();
  item->buffers[Color].allocate(colors.data(),
                          static_cast<int>(colors.size()*sizeof(float)));
  program->enableAttributeArray("colors");
  program->setAttributeBuffer("colors", GL_FLOAT, 0, 3);
  item->buffers[Color].release();

  item->buffers[Radius].bind();
  item->buffers[Radius].allocate(radius.data(),
                           static_cast<int>(radius.size()*sizeof(float)));
  program->enableAttributeArray("radius");
  program->setAttributeBuffer("radius", GL_FLOAT, 0, 1);
  item->buffers[Radius].release();

  item->buffers[Center].bind();
  item->buffers[Center].allocate(centers.data(),
                           static_cast<int>(centers.size()*sizeof(float)));
  program->enableAttributeArray("center");
  program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
  item->buffers[Center].release();

  viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("radius"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
  item->vaos[Facets]->release();


  item->vaos[Edges]->bind();
  item->buffers[Edge_vertices].bind();
  item->buffers[Edge_vertices].allocate(edges.data(),
                                  static_cast<int>(edges.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
  item->buffers[Edge_vertices].release();

  item->buffers[Normals].bind();
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
  item->buffers[Normals].release();

  item->buffers[Edge_color].bind();
  item->buffers[Edge_color].allocate(edges_colors.data(),
                               static_cast<int>(edges_colors.size()*sizeof(float)));
  program->enableAttributeArray("colors");
  program->setAttributeBuffer("colors", GL_FLOAT, 0, 3);
  item->buffers[Edge_color].release();

  item->buffers[Radius].bind();
  program->enableAttributeArray("radius");
  program->setAttributeBuffer("radius", GL_FLOAT, 0, 1);
  item->buffers[Radius].release();

  item->buffers[Center].bind();
  program->enableAttributeArray("center");
  program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
  item->buffers[Center].release();

  viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("radius"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
  item->vaos[Edges]->release();

  program->release();

  nb_centers = static_cast<int>(centers.size());
  centers.clear();
  centers.swap(centers);
  colors.clear();
  colors.swap(colors);
  radius.clear();
  radius.swap(radius);
  edges_colors.clear();
  edges_colors.swap(edges_colors);

  item->are_buffers_filled = true;
}

void Scene_spheres_item::draw(Viewer_interface *viewer) const
{
  if (!are_buffers_filled)
  {
    computeElements();
    d->initializeBuffers(viewer);
  }
  vaos[Scene_spheres_item_priv::Facets]->bind();
  if(d->has_plane)
  {
    d->program = getShaderProgram(PROGRAM_CUTPLANE_SPHERES, viewer);
    attribBuffers(viewer, PROGRAM_CUTPLANE_SPHERES);
    d->program->bind();
    QVector4D cp(d->plane.a(),d->plane.b(),d->plane.c(),d->plane.d());
    d->program->setUniformValue("cutplane", cp);

  }
  else
  {
    d->program = getShaderProgram(PROGRAM_SPHERES, viewer);
    attribBuffers(viewer, PROGRAM_SPHERES);
    d->program->bind();
  }
  viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                static_cast<GLsizei>(d->vertices.size()/3),
                                static_cast<GLsizei>(d->nb_centers));
  d->program->release();
  vaos[Scene_spheres_item_priv::Facets]->release();
}
void Scene_spheres_item::drawEdges(Viewer_interface *viewer) const
{
  if (!are_buffers_filled)
  {
    computeElements();
    d->initializeBuffers(viewer);
  }
  vaos[Scene_spheres_item_priv::Edges]->bind();
  if(d->has_plane)
  {
    d->program = getShaderProgram(PROGRAM_CUTPLANE_SPHERES, viewer);
    attribBuffers(viewer, PROGRAM_CUTPLANE_SPHERES);
    d->program->bind();
    QVector4D cp(d->plane.a(),d->plane.b(),d->plane.c(),d->plane.d());
    d->program->setUniformValue("cutplane", cp);
  }
  else
  {
    d->program = getShaderProgram(PROGRAM_SPHERES, viewer);
    attribBuffers(viewer, PROGRAM_SPHERES);
    d->program->bind();
  }
  viewer->glDrawArraysInstanced(GL_LINES, 0,
                                static_cast<GLsizei>(d->edges.size()/3),
                                static_cast<GLsizei>(d->nb_centers));
  d->program->release();
  vaos[Scene_spheres_item_priv::Edges]->release();
}
void Scene_spheres_item::add_sphere(CGAL::Sphere_3<Kernel> *sphere, CGAL::Color color)
{
  Scene_spheres_item::Sphere_pair pair_(sphere, color);
  d->spheres.append(pair_);
}

void Scene_spheres_item::remove_sphere(CGAL::Sphere_3<Kernel> *sphere)
{
  Q_FOREACH(Sphere_pair pair_, d->spheres)
    if(pair_.first == sphere)
    {
      d->spheres.removeAll(pair_);
      break;
    }
}

void Scene_spheres_item::clear_spheres()
{
  d->spheres.clear();
}
void Scene_spheres_item::setPrecision(int prec) { d->precision = prec; }
void Scene_spheres_item::setPlane(Kernel::Plane_3 p_plane) { d->plane = p_plane; }
void Scene_spheres_item::invalidateOpenGLBuffers(){are_buffers_filled = false;}
