#include "Scene_spheres_item.h"

#include <QApplication>
#include <CGAL/Three/Buffer_objects.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>

typedef CGAL::Three::Viewer_interface VI;
struct Scene_spheres_item_priv
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

  Scene_spheres_item_priv(bool planed, Scene_spheres_item* parent)
    :precision(36)
    ,has_plane(planed)

  {
    item = parent;
    isinit=false;
    create_flat_and_wire_sphere(1.0f,vertices,normals, edges);
    colors.clear();
    edges_colors.clear();
    centers.clear();
    radius.clear();
    if(has_plane)
    {
      face_container = new Triangle_container(VI::PROGRAM_CUTPLANE_SPHERES, false);
      edge_container = new Edge_container(VI::PROGRAM_CUTPLANE_SPHERES, false);
    }
    else
    {
      face_container = new Triangle_container(VI::PROGRAM_SPHERES, false);
      edge_container = new Edge_container(VI::PROGRAM_SPHERES, false);
    }
  }

  ~Scene_spheres_item_priv()
  {
    delete face_container;
    delete edge_container;
  }

  int precision;
  mutable CGAL::Plane_3<Kernel> plane;
  bool has_plane;
  bool isinit;

  mutable std::vector<float> vertices;
  mutable std::vector<float> normals;
  mutable std::vector<float> edges;
  mutable std::vector<float> colors;
  mutable std::vector<float> edges_colors;
  mutable std::vector<float> centers;
  mutable std::vector<float> radius;
  mutable QMap<CGAL::Three::Viewer_interface*, bool> buffers_init;
  mutable QOpenGLShaderProgram *program;
  Triangle_container* face_container;
  Edge_container* edge_container;
  mutable int nb_centers;
  Scene_spheres_item* item;
  QString tooltip;

};
Scene_spheres_item::Scene_spheres_item(Scene_group_item* parent, bool planed)
  :CGAL::Three::Scene_item()

{
  setParent(parent);
  d = new Scene_spheres_item_priv(planed, this);
}

Scene_spheres_item::~Scene_spheres_item()
{
  delete d;
}

void Scene_spheres_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  d->face_container->initializeBuffers(viewer);
  d->edge_container->initializeBuffers(viewer);
  d->face_container->flat_size = static_cast<int>(d->vertices.size());
  d->edge_container->flat_size = static_cast<int>(d->edges.size());
  d->face_container->center_size = d->centers.size();
  d->edge_container->center_size = d->centers.size();

//  centers.clear();
//  centers.swap(centers);
//  colors.clear();
//  colors.swap(colors);
//  radius.clear();
//  radius.swap(radius);
//  edges_colors.clear();
//  edges_colors.swap(edges_colors);


}

void Scene_spheres_item::draw(Viewer_interface *viewer) const
{
  if(!d->isinit)
    initGL();
  if(!is_locked && are_buffers_filled &&
     ! d->buffers_init[viewer])
  {
    initializeBuffers(viewer);
    d->buffers_init[viewer] = true;
  }
  if(d->has_plane)
  {
    QVector4D cp(d->plane.a(),d->plane.b(),d->plane.c(),d->plane.d());
    d->face_container->plane = cp;
    d->face_container->draw(*this, viewer,false);
  }
  else
  {
    d->face_container->draw(*this, viewer,false);
  }
}
void Scene_spheres_item::drawEdges(Viewer_interface *viewer) const
{
  if(!d->isinit)
    initGL();
  if(!is_locked && are_buffers_filled &&
     ! d->buffers_init[viewer])
  {
    initializeBuffers(viewer);
    d->buffers_init[viewer] = true;
  }
  if(d->has_plane)
  {
    QVector4D cp(d->plane.a(),d->plane.b(),d->plane.c(),d->plane.d());
    d->edge_container->plane = cp;
    d->edge_container->draw(*this, viewer,false);
  }
  else
  {
    d->edge_container->draw(*this, viewer,false);
  }
}
void Scene_spheres_item::add_sphere(const CGAL::Sphere_3<Kernel>& sphere, CGAL::Color color)
{
    d->colors.push_back((float)color.red()/255);
    d->colors.push_back((float)color.green()/255);
    d->colors.push_back((float)color.blue()/255);

    d->edges_colors.push_back((float)color.red()/255);
    d->edges_colors.push_back((float)color.green()/255);
    d->edges_colors.push_back((float)color.blue()/255);

    d->centers.push_back(sphere.center().x());
    d->centers.push_back(sphere.center().y());
    d->centers.push_back(sphere.center().z());

    d->radius.push_back(CGAL::sqrt(sphere.squared_radius()));
}

void Scene_spheres_item::clear_spheres()
{
  d->colors.clear();
  d->edges_colors.clear();
  d->centers.clear();
  d->radius.clear();
}
void Scene_spheres_item::setPrecision(int prec) { d->precision = prec; }
void Scene_spheres_item::setPlane(Kernel::Plane_3 p_plane) { d->plane = p_plane; }
void Scene_spheres_item::invalidateOpenGLBuffers(){are_buffers_filled = false;}

QString
Scene_spheres_item::toolTip() const {
    return d->tooltip;
}

void Scene_spheres_item::setToolTip(QString s)
{
  d->tooltip = s;
}

void Scene_spheres_item::setColor(QColor c)
{
  CGAL::Three::Scene_item::setColor(c);
  this->on_color_changed();
}

void Scene_spheres_item::computeElements() const
{}
