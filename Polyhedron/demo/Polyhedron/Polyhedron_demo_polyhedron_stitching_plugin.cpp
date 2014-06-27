#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Polyhedron_stitching.h>


class Polyhedron_demo_polyhedron_stitching_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionDetectBorders;
  QAction* actionStitchBorders;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionDetectBorders << actionStitchBorders; }
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* /* m */)
  {
    actionDetectBorders= new QAction(tr("Detect polyhedron boundaries"), mainWindow);
    actionStitchBorders= new QAction(tr("Stitch polyhedron duplicated boundaries"), mainWindow);
    actionDetectBorders->setObjectName("actionDetectBorders");
    actionStitchBorders->setObjectName("actionStitchBorders");
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable() const {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(index)) )
        return true;
    }
    return false;
  }

public slots:
  void on_actionDetectBorders_triggered();
  void on_actionStitchBorders_triggered();

}; // end Polyhedron_demo_polyhedron_stitching_plugin

void Polyhedron_demo_polyhedron_stitching_plugin::on_actionDetectBorders_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Scene_polylines_item* new_item = new Scene_polylines_item();

      Polyhedron* pMesh = item->polyhedron();
      pMesh->normalize_border();

      for (Polyhedron::Halfedge_iterator
              it=pMesh->border_halfedges_begin(), it_end=pMesh->halfedges_end();
              it!=it_end; ++it)
      {
        if (!it->is_border()) continue;
        /// \todo build cycles and graph with nodes of valence 2.
        new_item->polylines.push_back( Scene_polylines_item::Polyline() );
        new_item->polylines.back().push_back( it->opposite()->vertex()->point() );
        new_item->polylines.back().push_back( it->vertex()->point() );
      }
      if (new_item->polylines.empty())
      {
        delete new_item;
      }
      else
      {
        new_item->setName(tr("Boundary of %1").arg(item->name()));
        new_item->setColor(Qt::red);
        scene->addItem(new_item);
      }
    }
  }
}

void Polyhedron_demo_polyhedron_stitching_plugin::on_actionStitchBorders_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Polyhedron* pMesh = item->polyhedron();
      CGAL::polyhedron_stitching(*pMesh);
      scene->itemChanged(item);
    }
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_polyhedron_stitching_plugin, Polyhedron_demo_polyhedron_stitching_plugin)

#include "Polyhedron_demo_polyhedron_stitching_plugin.moc"
