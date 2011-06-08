#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item.h"
#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>

#include "Polyhedron_demo_plugin_helper.h"

class Polyhedron_demo_edit_polyhedron_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  Polyhedron_demo_edit_polyhedron_plugin() 
    : Polyhedron_demo_plugin_helper(), size(0)
  {}

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionToggleEdit;
  }

public slots:
  void on_actionToggleEdit_triggered();
  void edition();

private:
  QAction* actionToggleEdit;
  int size;
}; // end Polyhedron_demo_edit_polyhedron_plugin

void Polyhedron_demo_edit_polyhedron_plugin::init(QMainWindow* mainWindow, 
                                                  Scene_interface* scene_interface)
{
  actionToggleEdit = new QAction(tr("Toggle &edition of item(s)"), mainWindow);
  actionToggleEdit->setObjectName("actionToggleEdit");
  Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
}

void Polyhedron_demo_edit_polyhedron_plugin::on_actionToggleEdit_triggered() {
  bool found_polyhedron = false;
  bool edit_needed = false;
  QList<Scene_item*> changed_items;
  Q_FOREACH(Scene_interface::Item_id i, scene->selectionIndices())
  {
    if(Scene_polyhedron_item* poly_item = 
       qobject_cast<Scene_polyhedron_item*>(scene->item(i))) 
    {
      found_polyhedron = true;
      Scene_edit_polyhedron_item* edit_poly = 
        new Scene_edit_polyhedron_item(poly_item);
      edit_poly->setColor(poly_item->color());
      edit_poly->setName(QString("%1 (edit)").arg(poly_item->name()));
      scene->replaceItem(i, edit_poly);
      changed_items.push_back(scene->item(i));
      edit_needed = true;
      connect(edit_poly, SIGNAL(modified()),
              this, SLOT(edition()));
    } else if(Scene_edit_polyhedron_item* poly_item = 
              qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i))) 
    {
      found_polyhedron = true;
      scene->replaceItem(i, poly_item->to_polyhedron_item());
      delete poly_item;
    }
  }
  if(found_polyhedron == false) {
    QMessageBox::warning(mw, tr("Warning"),
                         tr("No polyhedron was selected"));
  }
  if(edit_needed) {
    size = QInputDialog::getInt(mw, 
                                tr("Polyhedron edition zone"),
                                tr("Size of edition zone:"),
                                size /* default value */ , 
                                0 /* min */ );
    std::cerr << "size = " << size << std::endl;
    Q_FOREACH(Scene_item* item, changed_items)
    {
      Scene_edit_polyhedron_item* poly_edit_item = 
        qobject_cast<Scene_edit_polyhedron_item*>(item);
      if(poly_edit_item) poly_edit_item->setZoneSize(size);
    }
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::edition() {
  QObject* obj = sender();
  Scene_edit_polyhedron_item* edit_item = 
    qobject_cast<Scene_edit_polyhedron_item*>(obj);
  if(!edit_item) {
    std::cerr << "ERROR" << __FILE__ << ":" << __LINE__ 
              << " : " << "unknown object type" << std::endl;
    return;
  }

  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  typedef Polyhedron::Vertex_handle Vertex_handle;

  const Point& orig = edit_item->original_position();
  const Vector move_vector = edit_item->current_position() - orig;

  Q_FOREACH(Vertex_handle vh, edit_item->selected_vertices())
  {
    vh->point() = vh->point() + move_vector;
  }
  edit_item->changed(); // that reset the original_position()
  scene->itemChanged(edit_item);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin, Polyhedron_demo_edit_polyhedron_plugin)

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
