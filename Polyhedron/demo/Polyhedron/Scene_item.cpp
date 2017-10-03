#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <QMenu>
#include <QWidgetAction>
#include <iostream>
#include <QDebug>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <QSlider>

namespace CT = CGAL::Three;
const QColor CT::Scene_item::defaultColor = QColor(100, 100, 255);
CGAL::Three::Scene_item::Scene_item()
  : name_("unnamed"),
    color_(defaultColor),
    visible_(true),
    are_buffers_filled(false),
    isinit(false),
    rendering_mode(FlatPlusEdges),
    defaultContextMenu(0),
    is_locked(false),
    is_reading(0)
{
  is_bbox_computed = false;
  is_diag_bbox_computed = false;
  nb_isolated_vertices = 0;
  has_group = 0;
  parent_group = 0;
  is_selected = false;
  alphaSlider = NULL;
}

CGAL::Three::Scene_item::~Scene_item() {
  if(defaultContextMenu)
    defaultContextMenu->deleteLater();
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



void CGAL::Three::Scene_item::itemAboutToBeDestroyed(CGAL::Three::Scene_item* item) {
    if(this == item)
    {
      Q_EMIT aboutToBeDestroyed();
    }
}


QString modeName(RenderingMode mode) {
    switch(mode)
    {
    case Points:
        return QObject::tr("points");
    case ShadedPoints:
        return QObject::tr("shaded points");
    case Wireframe:
        return QObject::tr("wire");
    case Flat:
        return QObject::tr("flat");
    case FlatPlusEdges:
        return QObject::tr("flat+edges");
    case Gouraud:
        return QObject::tr("Gouraud");
    case PointsPlusNormals:
        return QObject::tr("pts+normals");
    case Splatting:
        return QObject::tr("splats");
    default:
        Q_ASSERT(false);
        return QObject::tr("unknown");
    }
}

const char* slotName(RenderingMode mode) {
    switch(mode)
    {
    case Points:
        return SLOT(setPointsMode());
    case ShadedPoints:
      return SLOT(setShadedPointsMode());
    case Wireframe:
        return SLOT(setWireframeMode());
    case Flat:
        return SLOT(setFlatMode());
    case FlatPlusEdges:
        return SLOT(setFlatPlusEdgesMode());
    case Gouraud:
        return SLOT(setGouraudMode());
    case PointsPlusNormals:
        return SLOT(setPointsPlusNormalsMode());
    case Splatting:
        return SLOT(setSplattingMode());
    default:
        Q_ASSERT(false);
        return "";
    }
}

// Rendering mode as a human readable string
QString CGAL::Three::Scene_item::renderingModeName() const
{
    return modeName(renderingMode());
}
QMenu* CGAL::Three::Scene_item::contextMenu()
{
    if(defaultContextMenu) {
        defaultContextMenu->setTitle(name());
        return defaultContextMenu;
    }

    defaultContextMenu = new QMenu(name());
    for(unsigned int mode = 0; mode < NumberOfRenderingMode;
        ++mode)
    {
        if(!supportsRenderingMode(RenderingMode(mode))) continue;
        QString mName = modeName(RenderingMode(mode));
        defaultContextMenu->addAction(tr("Set %1 Mode")
                                      .arg(mName),
                                      this,
                                      slotName(RenderingMode(mode)));
    }
    QMenu *container = new QMenu(tr("Alpha value"));
    QWidgetAction *sliderAction = new QWidgetAction(0);
    connect(alphaSlider, &QSlider::valueChanged, this, &Scene_item::redraw);

    sliderAction->setDefaultWidget(alphaSlider);
    container->addAction(sliderAction);
    defaultContextMenu->addMenu(container);
    return defaultContextMenu;
}

CGAL::Three::Scene_group_item* CGAL::Three::Scene_item::parentGroup() const {
  return has_group > 0 ? parent_group : NULL;
}

void CGAL::Three::Scene_item::
moveToGroup(CGAL::Three::Scene_group_item* group) {
  parent_group = group;
  if(group)
    has_group = group->has_group + 1;
  else
    has_group = 0;
}

void CGAL::Three::Scene_item::invalidateOpenGLBuffers() {}

void CGAL::Three::Scene_item::selection_changed(bool) {}
void CGAL::Three::Scene_item::setVisible(bool b)
{
  visible_ = b;
  if(b)
    itemVisibilityChanged();
}


void CGAL::Three::Scene_item::select(double /*orig_x*/,
                        double /*orig_y*/,
                        double /*orig_z*/,
                        double /*dir_x*/,
                        double /*dir_y*/,
                        double /*dir_z*/)
{
}

// set-up the uniform attributes of the shader programs.
void CGAL::Three::Scene_item::attribBuffers(CGAL::Three::Viewer_interface* viewer,
                                             int program_name) const
{
    viewer->attribBuffers(program_name);
    viewer->getShaderProgram(program_name)->bind();
    if(is_selected)
        viewer->getShaderProgram(program_name)->setUniformValue("is_selected", true);
    else
        viewer->getShaderProgram(program_name)->setUniformValue("is_selected", false);

    viewer->getShaderProgram(program_name)->setUniformValue("alpha", alpha());
    QColor c = this->color();
    if(program_name == Scene_item::PROGRAM_WITH_TEXTURE)
    {
       if(is_selected) c = c.lighter(120);
       viewer->getShaderProgram(program_name)->setAttributeValue
         ("color_facets",
          c.redF(),
          c.greenF(),
          c.blueF());
    }
    else if(program_name == PROGRAM_WITH_TEXTURED_EDGES)
    {
        if(is_selected) c = c.lighter(50);
        viewer->getShaderProgram(program_name)->setUniformValue
          ("color_lines",
           QVector3D(c.redF(), c.greenF(), c.blueF()));
    }
    viewer->getShaderProgram(program_name)->release();
}


QOpenGLShaderProgram* CGAL::Three::Scene_item::getShaderProgram(int name, CGAL::Three::Viewer_interface * viewer) const
{
    return viewer->getShaderProgram(name);
}

CGAL::Three::Scene_item::Header_data CGAL::Three::Scene_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  return data;
}

QString CGAL::Three::Scene_item::computeStats(int )
{
  return QString();
}

#include <CGAL/double.h>

void CGAL::Three::Scene_item::compute_diag_bbox()const
{
 const Bbox& b_box = bbox();
  _diag_bbox = CGAL::sqrt(
        CGAL::square(b_box.xmax() - b_box.xmin())
        + CGAL::square(b_box.ymax() - b_box.ymin())
        + CGAL::square(b_box.zmax() - b_box.zmin())
        );

}

void CT::Scene_item::processData()const
{
  Scene_item* cthis = const_cast<Scene_item*>(this);
  cthis->writing();

  WorkerThread *workerThread = new WorkerThread(cthis);
  connect(workerThread, &WorkerThread::finished, [cthis, workerThread]()
  {
    cthis->doneWriting();
    cthis->redraw();
    cthis->dataProcessed();
  });
  connect(workerThread, &WorkerThread::finished, workerThread, &WorkerThread::deleteLater);
  workerThread->start();
}

void CT::Scene_item::initGL() const
{
  if(!alphaSlider)
  {
    alphaSlider = new QSlider(::Qt::Horizontal);
    alphaSlider->setMinimum(0);
    alphaSlider->setMaximum(255);
    alphaSlider->setValue(255);
  }
  Q_FOREACH(QGLViewer* v, QGLViewer::QGLViewerPool())
  {
    if(!v)
      continue;
    CT::Viewer_interface* viewer = static_cast<CT::Viewer_interface*>(v);
    Q_FOREACH(CT::Triangle_container* tc, triangle_containers)
    {
      if(!tc->is_gl_init[viewer])
        tc->initGL(*this, viewer);
    }
    Q_FOREACH(CT::Edge_container* ec, edge_containers)
    {
      if(!ec->is_gl_init[viewer])
        ec->initGL(*this, viewer);
    }
  }
  if(!isinit)
  {
    processData();
  }
  isinit = true;
}

void Scene_item::newViewer(Viewer_interface *viewer)
{
  processData();
  Q_FOREACH(CT::Triangle_container* tc, triangle_containers)
  {
    if(!tc->is_gl_init[viewer])
      tc->initGL(*this, viewer);
  }
  Q_FOREACH(CT::Edge_container* ec, edge_containers)
  {
    if(!ec->is_gl_init[viewer])
      ec->initGL(*this, viewer);
  }
}

void Scene_item::removeViewer(Viewer_interface *viewer)
{
  Q_FOREACH(CT::Triangle_container* tc, triangle_containers)
  {
    tc->removeViewer(viewer);
  }
  Q_FOREACH(CT::Edge_container* ec, edge_containers)
  {
    ec->removeViewer(viewer);
  }
}

float Scene_item::alpha() const
{
  if(!alphaSlider)
    return 1.0f;
  return (float)alphaSlider->value() / 255.0f;
}
