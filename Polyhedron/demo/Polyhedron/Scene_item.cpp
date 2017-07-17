#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <QMenu>
#include <iostream>
#include <QDebug>
#include <CGAL/Three/Viewer_interface.h>
const QColor CGAL::Three::Scene_item::defaultColor = QColor(100, 100, 255);

CGAL::Three::Scene_item::Scene_item(int buffers_size, int vaos_size)
  : name_("unnamed"),
    color_(defaultColor),
    visible_(true),
    are_buffers_filled(false),
    rendering_mode(FlatPlusEdges),
    defaultContextMenu(0),
    buffersSize(buffers_size),
    vaosSize(vaos_size),
    vaos(vaos_size)
{
  is_bbox_computed = false;
  is_diag_bbox_computed = false;
  for(int i=0; i<vaosSize; i++)
  {
    addVaos(i);
    vaos[i]->create();
  }

  buffers.reserve(buffersSize);
  for(int i=0; i<buffersSize; i++)
  {
    QOpenGLBuffer n_buf;
    buffers.push_back(n_buf);
    buffers[i].create();
  }
  nb_isolated_vertices = 0;
  has_group = 0;
  parent_group = 0;
  is_selected = false;
}

CGAL::Three::Scene_item::~Scene_item() {
  if(defaultContextMenu)
    defaultContextMenu->deleteLater();
  for(int i=0; i<buffersSize; i++)
  {
    buffers[i].destroy();
  }
  for(int i=0; i<vaosSize; i++)
  {
    delete vaos[i];
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
    // defaultContextMenu->addAction(name());
    // defaultContextMenu->addSeparator();
    // QMenu* modeMenu = new QMenu(QObject::tr("Rendering mode"),
    //                             defaultContextMenu);
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
    // defaultContextMenu->addAction(modeMenu->menuAction());
    return defaultContextMenu;
}

CGAL::Three::Scene_group_item* CGAL::Three::Scene_item::parentGroup() const {
  return parent_group;
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
    if(viewer == 0)
        viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
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
