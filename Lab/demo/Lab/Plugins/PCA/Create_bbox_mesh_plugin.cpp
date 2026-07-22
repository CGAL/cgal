#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>

#include <CGAL/boost/graph/generators.h>

#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QApplication>

#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/CGAL_Lab_plugin_interface.h>

using namespace CGAL::Three;

class Create_bbox_mesh_plugin :
  public QObject,
  public CGAL::Three::CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*);

  QList<QAction*> actions() const;

  bool applicable(QAction*) const
  {
    bool at_least_one_non_empty = false;
    for(int index : scene->selectionIndices())
    {
      Scene_item* item = scene->item(index);
      if(!item->isFinite())
        return false;
      if(!item->isEmpty())
        at_least_one_non_empty = true;
    }

    return at_least_one_non_empty;
  }

protected:
  bool bbox(bool extended = false);

public Q_SLOTS:
  void createBbox()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const bool res = bbox();
    QApplication::restoreOverrideCursor();

    if(!res)
      QMessageBox::warning(mw, "Error", "Failed to compute the bounding box.");
  }

  void createExtendedBbox()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const bool res = bbox(true);
    QApplication::restoreOverrideCursor();

    if(!res)
      QMessageBox::warning(mw, "Error", "Failed to compute the extended bounding box");
  }

private:
  Scene_interface* scene;
  QMainWindow* mw;

  QAction* actionBbox;
  QAction* actionExtendedBbox;
};

void
Create_bbox_mesh_plugin::
init(QMainWindow* mainWindow,
     Scene_interface* scene_interface,
     Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;

  actionBbox = new QAction(tr("Create &Bounding Box"), mainWindow);
  actionBbox->setObjectName("createBboxMeshAction");
  connect(actionBbox, SIGNAL(triggered()),
          this, SLOT(createBbox()));

  actionExtendedBbox = new QAction(tr("Create &Extended Bounding Box"), mainWindow);
  actionExtendedBbox->setObjectName("createExtendedBboxMeshAction");
  connect(actionExtendedBbox, SIGNAL(triggered()),
          this, SLOT(createExtendedBbox()));
}

QList<QAction*>
Create_bbox_mesh_plugin::
actions() const
{
  return QList<QAction*>() << actionBbox << actionExtendedBbox;
}

bool
Create_bbox_mesh_plugin::
bbox(bool extended)
{
  Scene_interface::Bbox bbox;
  int item_count = 0;
  QString name;

  for(int index : scene->selectionIndices())
  {
    Scene_item* item = scene->item(index);
    if(item->isFinite() && !item->isEmpty())
    {
      if(item_count > 0)
      {
        bbox = bbox + item->bbox();
        if(item_count == 1)
          name = name + " and others";
      }
      else
      {
        bbox = item->bbox();
        name = item->name();
      }
      ++item_count;
    }
  }

  if(item_count == 0)
    return false;

  std::cout << "bbox: " << bbox.xmin() << " | "  << bbox.xmax()
            << "\n      " << bbox.ymin() << " | " << bbox.ymax()
            << "\n      " << bbox.zmin() << " | " << bbox.zmax()
            << std::endl;

  if(extended)
  {
    const double delta_x = ( bbox.xmax() - bbox.xmin() ) / 20.;
    const double delta_y = ( bbox.ymax() - bbox.ymin() ) / 20.;
    const double delta_z = ( bbox.zmax() - bbox.zmin() ) / 20.;
    bbox = Scene_interface::Bbox(bbox.xmin() - delta_x,
                                 bbox.ymin() - delta_y,
                                 bbox.zmin() - delta_z,
                                 bbox.xmax() + delta_x,
                                 bbox.ymax() + delta_y,
                                 bbox.zmax() + delta_z);

    std::cout << "extended bbox: " << bbox.xmin() << " | "  << bbox.xmax()
              << "\n               " << bbox.ymin() << " | " << bbox.ymax()
              << "\n               " << bbox.zmin() << " | " << bbox.zmax()
              << std::endl;
  }

  if((bbox.min)(0) > (bbox.max)(0) ||
     (bbox.min)(1) > (bbox.max)(1) ||
     (bbox.min)(2) > (bbox.max)(2))
  {
    return false;
  }

  Scene_item* item;
  EPICK::Iso_cuboid_3 ic(bbox);
  SMesh* p = new SMesh;
  CGAL::make_hexahedron(ic, *p);

  item = new Scene_surface_mesh_item(p);
  item->setName(name + (extended ? " (Extended Bbox)" : " (Bbox)"));
  item->setRenderingMode(Wireframe);
  item->setColor(Qt::black);
  scene->addItem(item);

  return true;
}

#include "Create_bbox_mesh_plugin.moc"
