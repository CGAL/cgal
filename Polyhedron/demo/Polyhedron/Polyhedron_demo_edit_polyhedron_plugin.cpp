#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item.h"
#include <QAction>
#include <QMainWindow>

#include "Polyhedron_demo_plugin_interface.h"

class Polyhedron_demo_edit_polyhedron_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionEdit;
  }

public slots:
  void edit();

private:
  Scene_interface* scene;
  QAction* actionEdit;

}; // end Polyhedron_demo_edit_polyhedron_plugin

void Polyhedron_demo_edit_polyhedron_plugin::init(QMainWindow* mainWindow, 
                                                  Scene_interface* scene_interface)
{
  scene = scene_interface;
  actionEdit = new QAction(tr("Toggle &edition of item(s)"), mainWindow);
  connect(actionEdit, SIGNAL(triggered()),
          this, SLOT(edit()));
}

void Polyhedron_demo_edit_polyhedron_plugin::edit() {
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    if(Scene_polyhedron_item* poly_item = 
       qobject_cast<Scene_polyhedron_item*>(scene->item(i))) {
      Scene_edit_polyhedron_item* edit_poly = 
        new Scene_edit_polyhedron_item(poly_item);
      edit_poly->setColor(poly_item->color());
      edit_poly->setName(QString("%1 (edit)").arg(poly_item->name()));
      scene->replaceItem(i, edit_poly);
    } else if(Scene_edit_polyhedron_item* poly_item = 
              qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i))) {
      scene->replaceItem(i, poly_item->to_polyhedron_item());
      delete poly_item;
    }
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin, Polyhedron_demo_edit_polyhedron_plugin)

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
