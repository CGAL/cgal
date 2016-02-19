#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>

#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Messages_interface.h"
using namespace CGAL::Three;
class Polyhedron_demo_orient_soup_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m);

  bool applicable(QAction* action) const {
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices()) {
      if(qobject_cast<Scene_polygon_soup_item*>(scene->item(index)))
        return true;
      else
        if (action==actionShuffle && 
            qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
            return true;
    }
    return false;
  }

  QList<QAction*> actions() const;

public Q_SLOTS:
  void orient();
  void shuffle();
  void displayNonManifoldEdges();

private:
  CGAL::Three::Scene_interface* scene;
  Messages_interface* messages;
  QMainWindow* mw;
  QAction* actionOrient;
  QAction* actionShuffle;
  QAction* actionDisplayNonManifoldEdges;

}; // end Polyhedron_demo_orient_soup_plugin

void Polyhedron_demo_orient_soup_plugin::init(QMainWindow* mainWindow,
                                              CGAL::Three::Scene_interface* scene_interface,
                                              Messages_interface* m)
{
  scene = scene_interface;
  mw = mainWindow;
  messages = m;
  actionOrient = new QAction(tr("&Orient Polygon Soup"), mainWindow);
  actionOrient->setObjectName("actionOrient");
  actionOrient->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionOrient, SIGNAL(triggered()),
          this, SLOT(orient()));

  actionShuffle = new QAction(tr("&Shuffle Polygon Soup"), mainWindow);
  actionShuffle->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionShuffle, SIGNAL(triggered()),
          this, SLOT(shuffle()));

  actionDisplayNonManifoldEdges = new QAction(tr("Display Non Manifold Edges"),
                                              mainWindow);
  actionDisplayNonManifoldEdges->setProperty("subMenuName", "View");
  connect(actionDisplayNonManifoldEdges, SIGNAL(triggered()),
          this, SLOT(displayNonManifoldEdges()));
}

QList<QAction*> Polyhedron_demo_orient_soup_plugin::actions() const {
  return QList<QAction*>() << actionOrient
                           << actionShuffle
                           << actionDisplayNonManifoldEdges;
}

void Polyhedron_demo_orient_soup_plugin::orient()
{
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polygon_soup_item* item = 
      qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(item)
    {
      //     qDebug()  << tr("I have the item %1\n").arg(item->name());
      QApplication::setOverrideCursor(Qt::WaitCursor);
      if(!item->orient()) {
         QMessageBox::information(mw, tr("Not orientable without self-intersections"),
                                      tr("The polygon soup \"%1\" is not directly orientable."
                                         " Some vertices have been duplicated and some self-intersections"
                                         " have been created.")
                                      .arg(item->name()));
      }

      Scene_polyhedron_item* poly_item = new Scene_polyhedron_item();
      if(item->exportAsPolyhedron(poly_item->polyhedron())) {
        poly_item->setName(item->name());
        poly_item->setColor(item->color());
        poly_item->setRenderingMode(item->renderingMode());
        poly_item->setVisible(item->visible());
        poly_item->invalidateOpenGLBuffers();
        poly_item->setProperty("source filename", item->property("source filename"));
        scene->replaceItem(index, poly_item);
        delete item;
      } else {
        item->invalidateOpenGLBuffers();
        scene->itemChanged(item);
      }

      QApplication::restoreOverrideCursor();
    }
    else{
      messages->warning(tr("This function is only applicable on polygon soups."));
    }
  }
}

void Polyhedron_demo_orient_soup_plugin::shuffle()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polygon_soup_item* item = 
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(item) {
    item->shuffle_orientations();
    //scene->itemChanged(item);
  }
  else {
    Scene_polyhedron_item* poly_item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(poly_item) {
      item = new Scene_polygon_soup_item();
      item->setName(poly_item->name());
      item->setRenderingMode(poly_item->renderingMode());
      item->setVisible(poly_item->visible());
      item->setProperty("source filename", poly_item->property("source filename"));
      item->load(poly_item);
      item->shuffle_orientations();
      item->setColor(poly_item->color());
      scene->replaceItem(index, item);
      delete poly_item;
    }
  }
}

void Polyhedron_demo_orient_soup_plugin::displayNonManifoldEdges()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polygon_soup_item* item = 
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(item)
  {
    item->setDisplayNonManifoldEdges(!item->displayNonManifoldEdges());
    scene->itemChanged(item);
  }
}

#include "Orient_soup_plugin.moc"
