#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>

#include "Scene_polygon_soup.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "Messages_interface.h"

class Polyhedron_demo_orient_soup_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface* m);
  QList<QAction*> actions() const;

public slots:
  void orient();
  void shuffle();
  void displayNonManifoldEdges();

private:
  Scene_interface* scene;
  Messages_interface* messages;
  QMainWindow* mw;
  QAction* actionOrient;
  QAction* actionShuffle;
  QAction* actionDisplayNonManifoldEdges;

}; // end Polyhedron_demo_orient_soup_plugin

void Polyhedron_demo_orient_soup_plugin::init(QMainWindow* mainWindow,
                                              Scene_interface* scene_interface,
                                              Messages_interface* m)
{
  scene = scene_interface;
  mw = mainWindow;
  messages = m;
  actionOrient = new QAction(tr("Orient polygon soup"), mainWindow);
  connect(actionOrient, SIGNAL(triggered()),
          this, SLOT(orient()));

  actionShuffle = new QAction(tr("Shuffle polygon soup"), mainWindow);
  connect(actionShuffle, SIGNAL(triggered()),
          this, SLOT(shuffle()));

  actionDisplayNonManifoldEdges = new QAction(tr("Display non manifold edges"),
                                              mainWindow);
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
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polygon_soup* item = 
    qobject_cast<Scene_polygon_soup*>(scene->item(index));

  if(item)
  {
//     qDebug()  << tr("I have the item %1\n").arg(item->name());
    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(!item->orient())
      messages->warning(tr("The polygon soup \"%1\" is not orientable.")
                       .arg(item->name()));
 //      QMessageBox::information(mw, tr("Not orientable"),
//                                tr("The polygon soup \"%1\" is not orientable.")
//                                .arg(item->name()));

    scene->itemChanged(item);
    
    QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_orient_soup_plugin::shuffle()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polygon_soup* item = 
    qobject_cast<Scene_polygon_soup*>(scene->item(index));

  if(item)
    item->shuffle_orientations();
  scene->itemChanged(item);
}

void Polyhedron_demo_orient_soup_plugin::displayNonManifoldEdges()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polygon_soup* item = 
    qobject_cast<Scene_polygon_soup*>(scene->item(index));

  if(item)
  {
    item->setDisplayNonManifoldEdges(!item->displayNonManifoldEdges());
    item->changed();
    scene->itemChanged(item);
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_orient_soup_plugin, Polyhedron_demo_orient_soup_plugin);

#include "Polyhedron_demo_orient_soup_plugin.moc"
