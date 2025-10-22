#include "Variational_medial_axis_skeleton_item.h"
#include <QApplication>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/constraint.h>
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLShaderProgram>
#include <QCursor>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Variational_medial_axis_sampling.h>
#include <sstream>

using namespace CGAL::Three;

typedef Scene_surface_mesh_item::Face_graph Mesh;

struct Variational_medial_axis_skeleton_item_priv
{
  CGAL::Variational_medial_axis_sampling<Mesh, CGAL::Parallel_if_available_tag> vmas;

  Variational_medial_axis_skeleton_item_priv(const Mesh& m)
    : vmas(m)
  {}
};


Variational_medial_axis_skeleton_item::Variational_medial_axis_skeleton_item(CGAL::Three::Scene_interface* scene_interface,
                                                                             const Scene_surface_mesh_item *sm_item,
                                                                             std::size_t nb_spheres)
{
  d = new Variational_medial_axis_skeleton_item_priv(*sm_item->face_graph());

  connect(Three::mainWindow(), SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));

  scene = scene_interface;

  d->vmas.sample(
      CGAL::parameters::number_of_spheres(nb_spheres));
}

void Variational_medial_axis_skeleton_item::fill_subitems()
{
  auto skeleton = d->vmas.export_skeleton();

  spheres = new Scene_spheres_item(this, 0, false, false);
  for (std::size_t i=0; i<skeleton.vertices().size(); ++i)
    spheres->add_sphere(skeleton.vertices()[i], i);

  spheres->setName("Medial Spheres");
  spheres->setRenderingMode(Gouraud);
  scene->addItem(spheres);
  spheres->computeElements();
  scene->changeGroup(spheres, this);
  lockChild(spheres);


  edges = new Scene_polylines_item();
  edges->setName("Skeleton Edges");
  for(const std::pair<std::size_t, std::size_t>& p : skeleton.edges())
  {
    edges->polylines.emplace_back();
    edges->polylines.back().push_back(skeleton.vertices()[p.first].center());
    edges->polylines.back().push_back(skeleton.vertices()[p.second].center());
  }
  scene->addItem(edges);
  edges->computeElements();
  scene->changeGroup(edges, this);
  lockChild(edges);

  faces = new Scene_polygon_soup_item();
  faces->setName("Skeleton Faces");
  faces->init_polygon_soup(skeleton.vertices().size(), skeleton.faces().size());
  for (std::size_t i=0; i<skeleton.vertices().size(); ++i)
    faces->new_vertex(skeleton.vertices()[i].center().x(),
                      skeleton.vertices()[i].center().y(),
                      skeleton.vertices()[i].center().z());
  for (const std::array<std::size_t, 3>& a : skeleton.faces())
    faces->new_triangle(a[0], a[1], a[2]);

  scene->addItem(faces);
  faces->computeElements();
  scene->changeGroup(faces, this);
  lockChild(faces);
}

Variational_medial_axis_skeleton_item::~Variational_medial_axis_skeleton_item()
{
  delete d;
}
