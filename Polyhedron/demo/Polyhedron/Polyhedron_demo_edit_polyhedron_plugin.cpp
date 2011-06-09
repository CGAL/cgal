#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item.h"
#include <QAction>
#include <QKeySequence>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>
#include <QSettings>

#include <QDockWidget>
#include "ui_Deform_mesh.h"

#include "Polyhedron_type.h"  // defines the Polyhedron type

//#include <CGAL/Deform_mesh_BGL.h> // TODO: uncomment that! Do not compile
                                    // on Linux at the moment.

// Fake declaration needed because I cannot include
// <CGAL/Deform_mesh_BGL.h>.  TODO: replace that be the actual #include
namespace CGAL {
  template <typename Polyhedron>
  struct Deform_mesh_BGL {
    Deform_mesh_BGL(Polyhedron&) {}
  };
}

typedef CGAL::Deform_mesh_BGL<Polyhedron> Deform_mesh;

struct Polyhedron_deformation_data {
  Deform_mesh* deform_mesh;
  Polyhedron* polyhedron_copy; // For a possible undo operation. To be written.
};

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

  ~Polyhedron_demo_edit_polyhedron_plugin();

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionToggleEdit;
  }

public slots:
  void on_actionToggleEdit_triggered();
  void edition();

  void item_destroyed();

private:
  typedef std::map<QObject*, Polyhedron_deformation_data> Deform_map;
  Deform_map deform_map;

  Ui::DeformMesh deform_mesh_widget;
  QDockWidget* widget;

  QAction* actionToggleEdit;
  int size;
}; // end Polyhedron_demo_edit_polyhedron_plugin

Polyhedron_demo_edit_polyhedron_plugin::
~Polyhedron_demo_edit_polyhedron_plugin()
{
  QSettings settings;
  settings.beginGroup("Polyhedron edition");
  settings.setValue("Deform_mesh widget area", 
                    this->mw->dockWidgetArea(widget));
  settings.endGroup();
}

void Polyhedron_demo_edit_polyhedron_plugin::init(QMainWindow* mainWindow, 
                                                  Scene_interface* scene_interface)
{
  actionToggleEdit = new QAction(tr("Toggle &edition of item(s)"), mainWindow);
  actionToggleEdit->setObjectName("actionToggleEdit");
  actionToggleEdit->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_E));
  
  QSettings settings;
  settings.beginGroup("Polyhedron edition");
  int i = settings.value("Deform_mesh widget area", 
                         Qt::LeftDockWidgetArea).toInt();
  Qt::DockWidgetArea area = static_cast<Qt::DockWidgetArea>(i);
  settings.endGroup();
  widget = new QDockWidget();
  deform_mesh_widget.setupUi(widget);
  mainWindow->addDockWidget(area, widget);

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
      connect(edit_poly, SIGNAL(destroyed()),
              this, SLOT(item_destroyed()));
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

// Remove from 'deform_map' the metadata that corresponds to the deleted
// item.
void Polyhedron_demo_edit_polyhedron_plugin::item_destroyed() {
  QObject* obj = sender(); // the item that is destroyed
  Deform_map::iterator it = deform_map.find(obj);
  if(it != deform_map.end()) {
    delete it->second.polyhedron_copy;
    //    delete it->second.deform_mesh;  // TODO: uncomment that!
    deform_map.erase(it);
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

  Polyhedron* polyhedron = edit_item->polyhedron();

  Deform_mesh* deform = 0;  // Will be initialized below...

  Deform_map::iterator deform_it = deform_map.find(edit_item);
  if(deform_it == deform_map.end()) {
    // First time. Need to create the Deform_mesh object.
    deform = new Deform_mesh(*polyhedron);

    // create a new entry in the map
    Polyhedron_deformation_data& data = deform_map[edit_item]; 

    data.deform_mesh = deform;
    data.polyhedron_copy = new Polyhedron(*polyhedron); // copy

  } else {
    // In that case deform_it->second is the data associated to 'polyhedron'.
    deform = deform_it->second.deform_mesh;
  }

  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  typedef Polyhedron::Vertex_handle Vertex_handle;

  const Point& orig = edit_item->original_position();
  Q_UNUSED(orig)
  const Point& last_position = edit_item->last_position();
  const Vector translation = edit_item->current_position() - last_position;

  // ACTUAL DEFORMATION
  // This should be modified to use Deform_mesh instead.
  Q_FOREACH(Vertex_handle vh, edit_item->selected_vertices())
  {
    vh->point() = vh->point() + translation;
  }
  // END OF ACTUAL DEFORMATION

  // signal to the item that it needs to recompute its internal structures
  edit_item->changed(); // that reset the last_position()

  // signal to the scene that the item needs to be redrawn.
  scene->itemChanged(edit_item);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin, Polyhedron_demo_edit_polyhedron_plugin)

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
