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
    : Polyhedron_demo_plugin_helper(), size(0), edit_mode(false)
  {}

  ~Polyhedron_demo_edit_polyhedron_plugin();

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionToggleEdit;
  }

public slots:
  void on_actionToggleEdit_triggered(bool);
  void edition();

  void item_destroyed();
  void new_item_created(int item_id);

  void update_handlesRegionSize(int interestRegionSizeValue) {
    if(deform_mesh_widget.handlesRegionSize->value() > interestRegionSizeValue)
    {
      deform_mesh_widget.handlesRegionSize->setValue(interestRegionSizeValue);
    }
  }
  void update_interestRegionSize(int handlesRegionSizeValue) {
    if(deform_mesh_widget.interestRegionSize->value() < handlesRegionSizeValue)
    {
      deform_mesh_widget.interestRegionSize->setValue(handlesRegionSizeValue);
    }
  }

private:
  typedef Scene_interface::Item_id Item_id;

  Scene_edit_polyhedron_item* 
  convert_to_edit_polyhedron(Item_id, Scene_polyhedron_item*);

  Scene_polyhedron_item* 
  convert_to_plain_polyhedron(Item_id, Scene_edit_polyhedron_item*);

  typedef std::map<QObject*, Polyhedron_deformation_data> Deform_map;
  Deform_map deform_map;

  Ui::DeformMesh deform_mesh_widget;
  QDockWidget* widget;

  QAction* actionToggleEdit;
  int size;
  bool edit_mode;
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
  actionToggleEdit->setCheckable(true);
  
  QSettings settings;
  settings.beginGroup("Polyhedron edition");
  int i = settings.value("Deform_mesh widget area", 
                         Qt::LeftDockWidgetArea).toInt();
  Qt::DockWidgetArea area = static_cast<Qt::DockWidgetArea>(i);
  settings.endGroup();
  widget = new QDockWidget();
  deform_mesh_widget.setupUi(widget);
  mainWindow->addDockWidget(area, widget);

  // bind states of actionToggleEdit and editModeCb
  connect(actionToggleEdit, SIGNAL(triggered(bool)),
          deform_mesh_widget.editModeCb, SLOT(setChecked(bool)));
  connect(deform_mesh_widget.editModeCb, SIGNAL(clicked(bool)),
          actionToggleEdit, SLOT(setChecked(bool)));

  // make editModeCb actually trigger the slot
  connect(deform_mesh_widget.editModeCb, SIGNAL(clicked(bool)),
          this, SLOT(on_actionToggleEdit_triggered(bool)));

  // Connect Scene::newItem so that, if edit_mode==true, convert
  // automatically polyhedron items to "edit polyhedron" items.
  QObject* scene = dynamic_cast<QObject*>(scene_interface);
  if(scene) {
    connect(scene, SIGNAL(newItem(int)),
            this, SLOT(new_item_created(int)));
  } else {
    std::cerr << "ERROR " << __FILE__ << ":" << __LINE__ << " :"
              << " cannot convert scene_interface to scene!\n"; 
  }

  // Make sure handlesRegionSize->value() is always smaller than 
  // interestRegionSize->value()
  connect(deform_mesh_widget.handlesRegionSize, SIGNAL(valueChanged(int)),
          this, SLOT(update_interestRegionSize(int)));
  connect(deform_mesh_widget.interestRegionSize, SIGNAL(valueChanged(int)),
          this, SLOT(update_handlesRegionSize(int)));

  Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
}

void
Polyhedron_demo_edit_polyhedron_plugin::new_item_created(int item_id)
{
  if(edit_mode) {
    Scene_polyhedron_item* poly_item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
    if(poly_item) {
      convert_to_edit_polyhedron(item_id, poly_item);
    }
  }
}

Scene_edit_polyhedron_item*
Polyhedron_demo_edit_polyhedron_plugin::
convert_to_edit_polyhedron(Item_id i,
                           Scene_polyhedron_item* poly_item)
{
  QString poly_item_name = poly_item->name();
  Scene_edit_polyhedron_item* edit_poly = 
    new Scene_edit_polyhedron_item(poly_item);
  edit_poly->setColor(poly_item->color());
  edit_poly->setName(QString("%1 (edit)").arg(poly_item->name()));

  poly_item->setName(poly_item_name); // Because it is changed when the
                                      // name of edit_poly is changed.

  edit_poly->setVisible(poly_item->visible());
  edit_poly->setHandlesRegionSize(deform_mesh_widget.handlesRegionSize->value());
  edit_poly->setInterestRegionSize(deform_mesh_widget.interestRegionSize->value());
  connect(edit_poly, SIGNAL(modified()),
          this, SLOT(edition()));
  connect(edit_poly, SIGNAL(destroyed()),
          this, SLOT(item_destroyed()));
  connect(deform_mesh_widget.handlesRegionSize, SIGNAL(valueChanged(int)),
          edit_poly, SLOT(setHandlesRegionSize(int)));
  connect(deform_mesh_widget.interestRegionSize, SIGNAL(valueChanged(int)),
          edit_poly, SLOT(setInterestRegionSize(int)));
  scene->replaceItem(i, edit_poly);
  return edit_poly;
}

Scene_polyhedron_item*
Polyhedron_demo_edit_polyhedron_plugin::
convert_to_plain_polyhedron(Item_id i,
                            Scene_edit_polyhedron_item* edit_item) 
{
  Scene_polyhedron_item* poly_item = edit_item->to_polyhedron_item();
  scene->replaceItem(i, poly_item);
  delete edit_item;
  return poly_item;
}

void Polyhedron_demo_edit_polyhedron_plugin::on_actionToggleEdit_triggered(bool edit) {
  this->edit_mode = edit;
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    Scene_polyhedron_item* poly_item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(i));
    Scene_edit_polyhedron_item* edit_item = 
      qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(edit && poly_item) {
      convert_to_edit_polyhedron(i, poly_item);
    } else if(!edit && edit_item) {
      convert_to_plain_polyhedron(i, edit_item);
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

  // -- ACTUAL DEFORMATION --

  // This should be modified to use Deform_mesh instead.
  Q_FOREACH(Vertex_handle vh, edit_item->selected_vertices())
  {
    vh->point() = vh->point() + translation;
  }
  // could also use the list 'edit_item->vertices_in_region_of_interest()'

  // -- END OF ACTUAL DEFORMATION --

  // signal to the item that it needs to recompute its internal structures
  edit_item->changed(); // that reset the last_position()

  // signal to the scene that the item needs to be redrawn.
  scene->itemChanged(edit_item);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin, Polyhedron_demo_edit_polyhedron_plugin)

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
