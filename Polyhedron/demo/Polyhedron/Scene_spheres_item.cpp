#include "Scene_spheres_item.h"

#include <QApplication>
#include <CGAL/Three/Buffer_objects.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>

typedef CGAL::Three::Viewer_interface VI;
typedef Triangle_container Tri;
typedef Edge_container Ed;
struct Scene_spheres_item_priv
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

  Scene_spheres_item_priv(bool planed, Scene_spheres_item* parent)
    :precision(36)
    ,has_plane(planed)

  {
    item = parent;
    create_flat_and_wire_sphere(1.0f,vertices,normals, edges);
    colors.clear();
    edges_colors.clear();
    centers.clear();
    radius.clear();
    if(has_plane)
    {
      item->setTriangleContainer(0, new Triangle_container(VI::PROGRAM_CUTPLANE_SPHERES, false));
      item->setEdgeContainer(0, new Edge_container(VI::PROGRAM_CUTPLANE_SPHERES, false));
    }
    else
    {
      item->setTriangleContainer(0, new Triangle_container(VI::PROGRAM_SPHERES, false));
      item->setEdgeContainer(0, new Edge_container(VI::PROGRAM_SPHERES, false));
    }
  }


  int precision;
  mutable CGAL::Plane_3<Kernel> plane;
  bool has_plane;

  mutable std::vector<float> vertices;
  mutable std::vector<float> normals;
  mutable std::vector<float> edges;
  mutable std::vector<float> colors;
  mutable std::vector<float> edges_colors;
  mutable std::vector<float> centers;
  mutable std::vector<float> radius;
  mutable int nb_centers;
  Scene_spheres_item* item;
  QString tooltip;

};
Scene_spheres_item::Scene_spheres_item(Scene_group_item* parent, bool planed)
  :CGAL::Three::Scene_item_rendering_helper()

{
  setParent(parent);
  d = new Scene_spheres_item_priv(planed, this);
}

Scene_spheres_item::~Scene_spheres_item()
{
  delete d;
}

void Scene_spheres_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer)
{
  getTriangleContainer(0)->initializeBuffers(viewer);
  getEdgeContainer(0)->initializeBuffers(viewer);
  d->centers.clear();
  d->centers.shrink_to_fit();
  d->colors.clear();
  d->colors.shrink_to_fit();
  d->radius.clear();
  d->radius.shrink_to_fit();
  d->edges_colors.clear();
  d->edges_colors.shrink_to_fit();


}

void Scene_spheres_item::draw(Viewer_interface *viewer, int pass, bool depth_writing, QOpenGLFramebufferObject *fbo)
{
  if(!isWriting() && !isInit())
    initGL();
  if (!isWriting() && getBuffersFilled() &&
     ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  float near(viewer->camera()->zNear()), far(viewer->camera()->zFar());
  getTriangleContainer(0)->setComparing(pass > 0);
  getTriangleContainer(0)->setWidth(viewer->width()*1.0f);
  getTriangleContainer(0)->setHeight(viewer->height()*1.0f);
  getTriangleContainer(0)->setNear(near);
  getTriangleContainer(0)->setFar(far);
  getTriangleContainer(0)->setDepthWriting(depth_writing);
  getTriangleContainer(0)->setSelected(isSelected());
  getTriangleContainer(0)->setAlpha(alpha());
  if(d->has_plane)
  {
    QVector4D cp(d->plane.a(),d->plane.b(),d->plane.c(),d->plane.d());
    getTriangleContainer(0)->setPlane(cp);
    getTriangleContainer(0)->draw(viewer,false, fbo);
  }
  else
  {
    getTriangleContainer(0)->draw(viewer,false, fbo);
  }
}

void Scene_spheres_item::drawEdges(Viewer_interface *viewer)
{
  if(!isWriting() && !isInit())
    initGL();
  if (!isWriting() && getBuffersFilled() &&
      ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  getTriangleContainer(0)->setSelected(isSelected());
  if(d->has_plane)
  {
    QVector4D cp(d->plane.a(),d->plane.b(),d->plane.c(),d->plane.d());
    getEdgeContainer(0)->setPlane(cp);
    getEdgeContainer(0)->draw(viewer, false);
  }
  else
  {
    getEdgeContainer(0)->draw(viewer, false);
  }
}

void Scene_spheres_item::add_sphere(const CGAL::Sphere_3<Kernel>& sphere, CGAL::Color color)
{
    d->colors.push_back((float)color.red()/255);
    d->colors.push_back((float)color.green()/255);
    d->colors.push_back((float)color.blue()/255);
    d->colors.push_back(1.0f);
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
void Scene_spheres_item::invalidate(Gl_data_names){
  if(isWriting() || isReading())
    return;
  setBuffersFilled(false);
  Q_FOREACH(QGLViewer* v, QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(viewer == NULL)
      continue;
    setBuffersInit(viewer, false);
    viewer->update();
  }
  getTriangleContainer(0)->reset_vbos(GEOMETRY|NORMALS|COLORS);
  getEdgeContainer(0)->reset_vbos(GEOMETRY|NORMALS|COLORS);
  if(!isInit())
    initGL();
}

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

void Scene_spheres_item::computeElements(Gl_data_names)
{
  getTriangleContainer(0)->allocate(
        Tri::Flat_vertices,
        d->vertices.data(),
        static_cast<int>(d->vertices.size()*sizeof(float)));
  getTriangleContainer(0)->allocate(
        Tri::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));
  getTriangleContainer(0)->allocate(
        Tri::FColors,
        d->colors.data(),
        static_cast<int>(d->colors.size()*sizeof(float)));
  getTriangleContainer(0)->allocate(
        Tri::Facet_barycenters,
        d->centers.data(),
        static_cast<int>(d->centers.size()*sizeof(float)));
  getTriangleContainer(0)->allocate(
        Tri::Radius,
        d->radius.data(),
        static_cast<int>(d->radius.size()*sizeof(float)));

  getEdgeContainer(0)->allocate(
        Ed::Vertices,
        d->edges.data(),
        static_cast<int>(d->edges.size()*sizeof(float)));
  getEdgeContainer(0)->allocate(
        Ed::Normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));
  getEdgeContainer(0)->allocate(
        Ed::Colors,
        d->edges_colors.data(),
        static_cast<int>(d->edges_colors.size()*sizeof(float)));
  getEdgeContainer(0)->allocate(
        Ed::Barycenters,
        d->centers.data(),
        static_cast<int>(d->centers.size()*sizeof(float)));
  getEdgeContainer(0)->allocate(
        Ed::Radius,
        d->radius.data(),
        static_cast<int>(d->radius.size()*sizeof(float)));

  getTriangleContainer(0)->setCenterSize(static_cast<int>(d->centers.size()));
  getEdgeContainer(0)->setCenterSize(static_cast<int>(d->centers.size()));
  getTriangleContainer(0)->setFlatDataSize(static_cast<int>(d->vertices.size()));
  getEdgeContainer(0)->setFlatDataSize(static_cast<int>(d->edges.size()));
  setBuffersFilled(true);
}

void Scene_spheres_item::compute_bbox() const
{
}
