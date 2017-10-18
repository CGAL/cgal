#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <QSlider>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <QWidgetAction>
#include <QMenu>
#include <QApplication>
#include <QThreadPool>

using namespace CGAL::Three;

struct D{
  D(Scene_item_rendering_helper* item)
    : item(item),
      is_bbox_computed(false),
      is_diag_bbox_computed(false),
      _diag_bbox(0),
      alphaSlider(0),
      are_buffers_filled(false),
      isinit(false)
  {}

  ~D()
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
  }

  Scene_item_rendering_helper* item;
  Scene_item::Bbox _bbox;
  bool is_bbox_computed;
  bool is_diag_bbox_computed;
  double _diag_bbox;
  QSlider* alphaSlider;
  bool are_buffers_filled;
  bool isinit;
  QMap<CGAL::Three::Viewer_interface*, bool> buffers_init;
  std::vector<CGAL::Three::Triangle_container*> triangle_containers;
  std::vector<CGAL::Three::Edge_container*> edge_containers;

  void compute_diag_bbox();
};

Scene_item_rendering_helper::Scene_item_rendering_helper()
  :d(new D(this)){}

Scene_item_rendering_helper::~Scene_item_rendering_helper()
{
  if(d)
    delete d;
}

double Scene_item_rendering_helper::diagonalBbox() const {
  if(!d->is_diag_bbox_computed)
    d->compute_diag_bbox();
  d->is_diag_bbox_computed = true;
  return d->_diag_bbox;
}

void D::compute_diag_bbox()
{
  const Scene_item::Bbox& b_box = item->bbox();
  _diag_bbox = CGAL::sqrt(
        CGAL::square(b_box.xmax() - b_box.xmin())
        + CGAL::square(b_box.ymax() - b_box.ymin())
        + CGAL::square(b_box.zmax() - b_box.zmin())
        );
}

float Scene_item_rendering_helper::alpha() const
{
  if(!d->alphaSlider)
    return 1.0f;
  return (float)d->alphaSlider->value() / 255.0f;
}

void Scene_item_rendering_helper::initGL()
{
  if(!d->alphaSlider)
  {
    d->alphaSlider = new QSlider(::Qt::Horizontal);
    d->alphaSlider->setMinimum(0);
    d->alphaSlider->setMaximum(255);
    d->alphaSlider->setValue(255);

    connect(d->alphaSlider, &QSlider::valueChanged,
            [this](){redraw();});
  }

  Q_FOREACH(QGLViewer* v, QGLViewer::QGLViewerPool())
  {
    if(!v)
      continue;
    Viewer_interface* viewer = static_cast<Viewer_interface*>(v);
    Q_FOREACH(Triangle_container* tc, d->triangle_containers)
    {
      if(!tc->isGLInit(viewer))
        tc->initGL(viewer);
    }
    Q_FOREACH(Edge_container* ec, d->edge_containers)
    {
      if(!ec->isGLInit(viewer))
        ec->initGL(viewer);
    }
  }
  if(!d->isinit)
  {
    Gl_data_names flags;
    flags = (ALL);
    processData(flags);
  }
  d->isinit = true;
}

void Scene_item_rendering_helper::newViewer(Viewer_interface *viewer)
{
  processData(ALL); //newViewer shouldn't be called when item is locked;
  Q_FOREACH(Triangle_container* tc, d->triangle_containers)
  {
    if(!tc->isGLInit(viewer))
      tc->initGL(viewer);
  }
  Q_FOREACH(Edge_container* ec, d->edge_containers)
  {
    if(!ec->isGLInit(viewer))
      ec->initGL(viewer);
  }
}


void Scene_item_rendering_helper::removeViewer(Viewer_interface *viewer)
{
  Q_FOREACH(Triangle_container* tc, d->triangle_containers)
  {
    tc->removeViewer(viewer);
  }
  Q_FOREACH(Edge_container* ec, d->edge_containers)
  {
    ec->removeViewer(viewer);
  }
}

QMenu* Scene_item_rendering_helper::contextMenu()
{
  QMenu* resMenu = Scene_item::contextMenu();
  bool prop = property("menu_changed").toBool();
  if(!prop)
  {
    QMenu *container = new QMenu(tr("Alpha value"));
    QWidgetAction *sliderAction = new QWidgetAction(0);

    sliderAction->setDefaultWidget(d->alphaSlider);
    container->addAction(sliderAction);
    resMenu->addMenu(container);
    setProperty("menu_changed", true);
  }
  return resMenu;
}

void Scene_item_rendering_helper::processData(Gl_data_names name)
{
  writing();

  QApplication::setOverrideCursor(::Qt::BusyCursor);
  WorkerThread *workerThread = new WorkerThread(this, name);
  connect(workerThread, &WorkerThread::finished, [this, workerThread]()
  {
    doneWriting();
    redraw();
    dataProcessed();
    QApplication::restoreOverrideCursor();
  });
  QThreadPool::globalInstance()->start(workerThread);
}

void Scene_item_rendering_helper::setAlpha(int alpha)
{
  if(!d->alphaSlider)
    initGL();
  d->alphaSlider->setValue(alpha);
  redraw();
}

Scene_item::Bbox Scene_item_rendering_helper::bbox() const {
  if(!d->is_bbox_computed)
    compute_bbox();
  d->is_bbox_computed = true;
  return d->_bbox;
}

bool Scene_item_rendering_helper::isInit()const { return d->isinit; }

QSlider* Scene_item_rendering_helper::alphaSlider() { return d->alphaSlider; }

void Scene_item_rendering_helper::setBbox(Bbox b) { d->_bbox = b; }

Triangle_container* Scene_item_rendering_helper::getTriangleContainer(std::size_t id)const
{
  return d->triangle_containers[id];
}

Edge_container* Scene_item_rendering_helper::getEdgeContainer(std::size_t id)const
{
  return d->edge_containers[id];
}

void Scene_item_rendering_helper::setTriangleContainer(std::size_t id,
                                                       Triangle_container* tc)
{
  if(d->triangle_containers.size() <= id)
  {
    d->triangle_containers.resize(id+1);
  }
  d->triangle_containers[id] = tc;
}

void Scene_item_rendering_helper::setEdgeContainer(std::size_t id,
                                                   Edge_container* ec)
{
  if(d->edge_containers.size() <= id)
  {
    d->edge_containers.resize(id+1);
  }
  d->edge_containers[id] = ec;
}

void Scene_item_rendering_helper::setBuffersFilled(bool b)
{
  d->are_buffers_filled = b;
}

bool Scene_item_rendering_helper::getBuffersFilled()const
{
  return d->are_buffers_filled;
}

bool Scene_item_rendering_helper::getBuffersInit(Viewer_interface* viewer) const
{
  return d->buffers_init[viewer];
}

void Scene_item_rendering_helper::setBuffersInit(Viewer_interface* viewer, bool val)
{
  d->buffers_init[viewer] = val;
}
