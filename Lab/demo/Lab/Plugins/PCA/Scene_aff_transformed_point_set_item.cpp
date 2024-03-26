#include "Scene_aff_transformed_point_set_item.h"

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Point_container.h>

#include <QApplication>

using namespace CGAL::Three;

// const std::size_t limit_fast_drawing = 300000; //arbitraty large value

Scene_aff_transformed_point_set_item::
Scene_aff_transformed_point_set_item(Scene_points_with_normal_item* pts_item,
                                 const CGAL::qglviewer::Vec& pos)
  : Scene_aff_transformed_item(pos)
{
  d = new Scene_aff_transformed_point_set_item_priv(pts_item, pos);
  setPointContainer(0, new Point_container(Viewer_interface::PROGRAM_NO_SELECTION, false));
  compute_bbox();
  invalidateOpenGLBuffers();

  connect(pts_item, &Scene_points_with_normal_item::itemChanged ,
          this, &Scene_aff_transformed_point_set_item::updateCache);
}

Scene_aff_transformed_point_set_item::
~Scene_aff_transformed_point_set_item()
{
  delete d;
}

void
Scene_aff_transformed_point_set_item::
updateCache()
{
  compute_bbox();
  invalidateOpenGLBuffers();
}

QString
Scene_aff_transformed_point_set_item::
toolTip() const
{
  return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                      "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                      "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
      .arg(d->pts_item->name());
}

void
Scene_aff_transformed_point_set_item::
invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getPointContainer(0)->reset_vbos(ALL);
}

void
Scene_aff_transformed_point_set_item::
initializeBuffers(Viewer_interface* v) const
{
  Point_container* pc = getPointContainer(0);
  pc->initializeBuffers(v);
  pc->setFlatDataSize(d->nb_points);
  d->points.clear();
  d->points.shrink_to_fit();
}

void
Scene_aff_transformed_point_set_item::
compute_bbox() const
{
  setBbox(d->pts_item->bbox());
}

void
Scene_aff_transformed_point_set_item::
computeElements() const
{
  d->compute_elements();

  Point_container* pc = getPointContainer(0);
  pc->allocate(Point_container::Vertices,
               d->points.data(),
               static_cast<int>(d->points.size()*sizeof(float)));
  d->nb_points = d->points.size();

  setBuffersFilled(true);
}

void
Scene_aff_transformed_point_set_item::
drawPoints(CGAL::Three::Viewer_interface *viewer) const
{
  GLfloat point_size;
  viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
  viewer->setGlPointSize(6.f);
  double ratio_displayed = 1.0;
  // @fixme...
  // if((viewer->inFastDrawing () || d->pts_item->d->isPointSliderMoving()) &&
  //    ((d->nb_points)/3 > limit_fast_drawing)) // arbitrary large value
  //   ratio_displayed = 3 * limit_fast_drawing / static_cast<double>(d->nb_points);

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

  Point_container* pc = getPointContainer(0);
  pc->setColor(this->color());
  pc->setFrameMatrix(this->getFMatrix());

  std::size_t real_size = pc->getFlatDataSize();
  pc->setFlatDataSize(ratio_displayed * real_size);
  pc->draw(viewer, true);
  pc->setFlatDataSize(real_size);
}
