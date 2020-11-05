#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <QSlider>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Sqrt_extension.h>
#include <QWidgetAction>
#include <QMenu>
#include <QApplication>
#include <QThreadPool>

using namespace CGAL::Three;

struct PRIV{
  PRIV(Scene_item_rendering_helper* item)
    : item(item),
      is_bbox_computed(false),
      is_diag_bbox_computed(false),
      _diag_bbox(0),
      alphaSlider(0),
      are_buffers_filled(false)
  {}

  ~PRIV()
  {
    if(alphaSlider)
      delete alphaSlider;

    for(std::size_t i = 0; i < triangle_containers.size(); ++i)
    {
      delete triangle_containers[i];
    }
    for(std::size_t i = 0; i < edge_containers.size(); ++i)
    {
      delete edge_containers[i];
    }
    for(std::size_t i = 0; i < point_containers.size(); ++i)
    {
      delete point_containers[i];
    }
  }

  Scene_item_rendering_helper* item;
  Scene_item::Bbox _bbox;
  bool is_bbox_computed;
  bool is_diag_bbox_computed;
  double _diag_bbox;
  QSlider* alphaSlider;
  bool are_buffers_filled;
  std::map<CGAL::Three::Viewer_interface*, bool> isinit;
  QMap<CGAL::Three::Viewer_interface*, bool> buffers_init;
  std::vector<CGAL::Three::Triangle_container*> triangle_containers;
  std::vector<CGAL::Three::Edge_container*> edge_containers;
  std::vector<CGAL::Three::Point_container*> point_containers;

  void compute_diag_bbox();
};

Scene_item_rendering_helper::Scene_item_rendering_helper()
  :Scene_item(0,0),
    priv(new PRIV(this)){}

Scene_item_rendering_helper::~Scene_item_rendering_helper()
{
  if(priv)
    delete priv;
}

double Scene_item_rendering_helper::diagonalBbox() const {
  if(!priv->is_diag_bbox_computed)
    priv->compute_diag_bbox();
  priv->is_diag_bbox_computed = true;
  return priv->_diag_bbox;
}

void PRIV::compute_diag_bbox()
{
  const Scene_item::Bbox& b_box = item->bbox();
  _diag_bbox = CGAL::approximate_sqrt(
          CGAL::square(b_box.xmax() - b_box.xmin())
        + CGAL::square(b_box.ymax() - b_box.ymin())
        + CGAL::square(b_box.zmax() - b_box.zmin())
        );
}

float Scene_item_rendering_helper::alpha() const
{
  if(!priv->alphaSlider)
    return 1.0f;
  return (float)priv->alphaSlider->value() / 255.0f;
}

void Scene_item_rendering_helper::initGL(CGAL::Three::Viewer_interface* viewer) const
{
  if(!priv->alphaSlider)
  {
    priv->alphaSlider = new QSlider(::Qt::Horizontal);
    priv->alphaSlider->setMinimum(0);
    priv->alphaSlider->setMaximum(255);
    priv->alphaSlider->setValue(255);
  }

  Q_FOREACH(Triangle_container* tc, priv->triangle_containers)
  {
    if(!tc->isGLInit(viewer))
      tc->initGL(viewer);
  }
  Q_FOREACH(Edge_container* ec, priv->edge_containers)
  {
    if(!ec->isGLInit(viewer))
      ec->initGL(viewer);
  }
  Q_FOREACH(Point_container* pc, priv->point_containers)
  {
    if(!pc->isGLInit(viewer))
      pc->initGL(viewer);
  }
  if(!getBuffersFilled())
  {
    Gl_data_names flags;
    flags = (ALL);
    processData(flags);
  }
  priv->isinit[viewer] = true;
}

void Scene_item_rendering_helper::processData(Gl_data_names )const
{
  computeElements();
  priv->are_buffers_filled = true;
}


QMenu* Scene_item_rendering_helper::contextMenu()
{
  QMenu* resMenu = Scene_item::contextMenu();
  bool prop = property("menu_changed").toBool();
  if(!prop)
  {
    QMenu *container = new QMenu(tr("Alpha value"));
    QWidgetAction *sliderAction = new QWidgetAction(0);

    sliderAction->setDefaultWidget(priv->alphaSlider);
    container->addAction(sliderAction);
    resMenu->addMenu(container);
    setProperty("menu_changed", true);
  }
  return resMenu;
}

void Scene_item_rendering_helper::setAlpha(int alpha)
{
  if(!priv->alphaSlider)
  {
    priv->alphaSlider = new QSlider(::Qt::Horizontal);
    priv->alphaSlider->setMinimum(0);
    priv->alphaSlider->setMaximum(255);
    priv->alphaSlider->setValue(255);
  }
  priv->alphaSlider->setValue(alpha);
  redraw();
}

Scene_item::Bbox Scene_item_rendering_helper::bbox() const {
  if(!priv->is_bbox_computed)
    compute_bbox();
  priv->is_bbox_computed = true;
  return priv->_bbox;
}

bool Scene_item_rendering_helper::isInit(CGAL::Three::Viewer_interface* viewer)const
{
  if(priv->isinit.find(viewer) != priv->isinit.end())
    return priv->isinit[viewer];
  return false;
}

QSlider* Scene_item_rendering_helper::alphaSlider() { return priv->alphaSlider; }

void Scene_item_rendering_helper::setBbox(Bbox b)const { priv->_bbox = b; }

Triangle_container* Scene_item_rendering_helper::getTriangleContainer(std::size_t id)const
{
  return priv->triangle_containers[id];
}

Edge_container* Scene_item_rendering_helper::getEdgeContainer(std::size_t id)const
{
  return priv->edge_containers[id];
}

Point_container* Scene_item_rendering_helper::getPointContainer(std::size_t id)const
{
  return priv->point_containers[id];
}


void Scene_item_rendering_helper::setTriangleContainer(std::size_t id,
                                                       Triangle_container* tc)
{
  if(priv->triangle_containers.size() <= id)
  {
    priv->triangle_containers.resize(id+1, nullptr);
  }
  if(priv->triangle_containers[id])
    delete priv->triangle_containers[id];
  priv->triangle_containers[id] = tc;
}

void Scene_item_rendering_helper::setEdgeContainer(std::size_t id,
                                                   Edge_container* ec)
{
  if(priv->edge_containers.size() <= id)
  {
    priv->edge_containers.resize(id+1, nullptr);
  }
  if(priv->edge_containers[id])
    delete priv->edge_containers[id];
  priv->edge_containers[id] = ec;
}

void Scene_item_rendering_helper::setPointContainer(std::size_t id,
                                                   Point_container* pc)
{
  if(priv->point_containers.size() <= id)
  {
    priv->point_containers.resize(id+1, nullptr);
  }
  if(priv->point_containers[id])
    delete priv->point_containers[id];
  priv->point_containers[id] = pc;
}


void Scene_item_rendering_helper::setBuffersFilled(bool b) const
{
  priv->are_buffers_filled = b;
}

bool Scene_item_rendering_helper::getBuffersFilled()const
{
  return priv->are_buffers_filled;
}

bool Scene_item_rendering_helper::getBuffersInit(Viewer_interface* viewer) const
{
  return priv->buffers_init[viewer];
}

void Scene_item_rendering_helper::setBuffersInit(Viewer_interface* viewer, bool val)const
{
  priv->buffers_init[viewer] = val;
}

void Scene_item_rendering_helper::removeViewer(Viewer_interface *viewer)
{
  Q_FOREACH(Triangle_container* tc, priv->triangle_containers)
  {
    tc->removeViewer(viewer);
  }
  Q_FOREACH(Edge_container* ec, priv->edge_containers)
  {
    ec->removeViewer(viewer);
  }
  Q_FOREACH(Point_container* pc, priv->point_containers)
  {
    pc->removeViewer(viewer);
  }
}

void Scene_item_rendering_helper::newViewer(Viewer_interface *viewer)
{
  viewer->makeCurrent();
  Q_FOREACH(Triangle_container* tc, priv->triangle_containers)
  {
    if(!tc->isGLInit(viewer))
      tc->initGL(viewer);
  }
  Q_FOREACH(Edge_container* ec, priv->edge_containers)
  {
    if(!ec->isGLInit(viewer))
      ec->initGL(viewer);
  }
  Q_FOREACH(Point_container* pc, priv->point_containers)
  {
    if(!pc->isGLInit(viewer))
      pc->initGL(viewer);
  }
  //call function that recomputes vao binding
  initializeBuffers(viewer);
}
