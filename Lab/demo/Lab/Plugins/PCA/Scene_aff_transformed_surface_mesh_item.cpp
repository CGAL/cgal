#include "Scene_aff_transformed_surface_mesh_item.h"

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Edge_container.h>

#include <QApplication>

using namespace CGAL::Three;

Scene_aff_transformed_surface_mesh_item::
Scene_aff_transformed_surface_mesh_item(Scene_surface_mesh_item* sm_item,
                                    const CGAL::qglviewer::Vec& pos)
  : Scene_aff_transformed_item(pos)
{
  d = new Scene_aff_transformed_surface_mesh_item_priv(sm_item, pos);
  setEdgeContainer(0, new Edge_container(Viewer_interface::PROGRAM_NO_SELECTION, false));
  compute_bbox();
  invalidateOpenGLBuffers();

  connect(sm_item, &Scene_surface_mesh_item::itemChanged ,
          this, &Scene_aff_transformed_surface_mesh_item::updateCache);
}

Scene_aff_transformed_surface_mesh_item::
~Scene_aff_transformed_surface_mesh_item()
{
  delete d;
}

void
Scene_aff_transformed_surface_mesh_item::
updateCache()
{
  compute_bbox();
  invalidateOpenGLBuffers();
}

QString
Scene_aff_transformed_surface_mesh_item::
toolTip() const
{
  return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                      "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                      "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
          .arg(d->sm_item->name());
}

void
Scene_aff_transformed_surface_mesh_item::
invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getEdgeContainer(0)->reset_vbos(ALL);
}

void
Scene_aff_transformed_surface_mesh_item::
initializeBuffers(Viewer_interface *v) const
{
  getEdgeContainer(0)->initializeBuffers(v);
  getEdgeContainer(0)->setFlatDataSize(d->nb_lines);
  d->positions_lines.clear();
  d->positions_lines.shrink_to_fit();
}

void
Scene_aff_transformed_surface_mesh_item::
compute_bbox() const
{
  setBbox(d->sm_item->bbox());
}

void
Scene_aff_transformed_surface_mesh_item::
computeElements() const
{
  d->compute_elements();

  getEdgeContainer(0)->allocate(Edge_container::Vertices,
                                d->positions_lines.data(),
                                static_cast<int>(d->positions_lines.size()*sizeof(float)));
  d->nb_lines = d->positions_lines.size();

  setBuffersFilled(true);
}

void
Scene_aff_transformed_surface_mesh_item::
drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);

  if(getBuffersFilled() && !getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }

  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }

  Edge_container* ec = getEdgeContainer(0);
  ec->setColor(this->color());
  ec->setFrameMatrix(this->getFMatrix());
  ec->draw(viewer, true);
}