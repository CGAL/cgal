#define  CGAL_DEFORM_ROTATION

#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item.h"
#include <QAction>
#include <QKeySequence>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>
#include <QSettings>

#include <QGLViewer/qglviewer.h>

#include <QDockWidget>
#include "ui_Deform_mesh.h"

#include "Polyhedron_type.h"  // defines the Polyhedron type

#include <CGAL/Deform_mesh.h> 
#include <CGAL/Taucs_solver_traits.h>       

typedef Polyhedron_vertex_deformation_index_map<Polyhedron> Vertex_index_map;
typedef Polyhedron_edge_deformation_index_map<Polyhedron> Edge_index_map;

typedef CGAL::Deform_mesh<Polyhedron, CGAL::Taucs_solver_traits<double>, Vertex_index_map, Edge_index_map> Deform_mesh;

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Polyhedron::Vertex_handle Vertex_handle;

struct Polyhedron_deformation_data {
  Deform_mesh* deform_mesh;
  bool preprocessed;                            // specify whether preprocessed or not
  std::map<Vertex_handle, Vector> handle_vectors;  // record transform vectors of all handles, 
                                                   // only for multiple handle region scenario 
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
  void start_deform();
  void clear_handles();
  void preprocess(Scene_edit_polyhedron_item* edit_item);
  void complete_deform(Scene_edit_polyhedron_item* edit_item);
  void edition();
  void usage_scenario_0(Scene_edit_polyhedron_item* edit_item);
  void usage_scenario_1(Scene_edit_polyhedron_item* edit_item);

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
  connect(deform_mesh_widget.startDeformPb, SIGNAL(clicked(bool)),
          this, SLOT(start_deform()));
  connect(deform_mesh_widget.clearHandlesPb, SIGNAL(clicked(bool)),
          this, SLOT(clear_handles()));
  connect(deform_mesh_widget.usageScenarioCb, SIGNAL(currentIndexChanged(int)),
          this, SLOT(clear_handles()));

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
  edit_poly->setGeodesicCircle(deform_mesh_widget.geodesicCircleCb->isChecked());
  edit_poly->setSharpFeature(deform_mesh_widget.sharpFeatureCb->isChecked());
  edit_poly->setUsageScenario(deform_mesh_widget.usageScenarioCb->currentIndex());
  edit_poly->setSelectedHandlesMoved(false);
  edit_poly->setSelectedVertexChanged(false);
  connect(edit_poly, SIGNAL(modified()),
          this, SLOT(edition()));
  connect(edit_poly, SIGNAL(destroyed()),
          this, SLOT(item_destroyed()));
  connect(deform_mesh_widget.handlesRegionSize, SIGNAL(valueChanged(int)),
          edit_poly, SLOT(setHandlesRegionSize(int)));
  connect(deform_mesh_widget.interestRegionSize, SIGNAL(valueChanged(int)),
          edit_poly, SLOT(setInterestRegionSize(int)));
  connect(deform_mesh_widget.geodesicCircleCb, SIGNAL(clicked(bool)),
          edit_poly, SLOT(setGeodesicCircle(bool)));
  connect(deform_mesh_widget.sharpFeatureCb, SIGNAL(clicked(bool)),
          edit_poly, SLOT(setSharpFeature(bool)));
  connect(deform_mesh_widget.usageScenarioCb, SIGNAL(currentIndexChanged(int)),
          edit_poly, SLOT(setUsageScenario(int)));
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
    delete it->second.deform_mesh;  // TODO: uncomment that!
    deform_map.erase(it);
  }
}

// Clear all the existing handles and ROI and recover mesh to original one
void Polyhedron_demo_edit_polyhedron_plugin::clear_handles() {

  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item =
    qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  Polyhedron* polyhedron = edit_item->polyhedron();
  Deform_map::iterator deform_it = deform_map.find(edit_item);
  if(deform_it != deform_map.end())
  {
    Polyhedron_deformation_data& data = deform_map[edit_item];
    Deform_mesh* deform = data.deform_mesh;
    deform->clear_handles();
    deform->clear_roi();
    data.preprocessed = false;
  }
  edit_item->clear_selected_handles();
  edit_item->clear_non_selected_handles();
  edit_item->clear_selected_roi();
  edit_item->clear_non_selected_roi();
  edit_item->clear_selected_vectors();
  edit_item->clear_non_selected_vectors();

  // signal to the item that it needs to recompute its internal structures
  edit_item->changed(); // that reset the last_position()

  // signal to the scene that the item needs to be redrawn.
  scene->itemChanged(edit_item);


  
}


// Start deformation: preprocess or complete deformation
void Polyhedron_demo_edit_polyhedron_plugin::start_deform() {

  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item =
    qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  // do deformation only when handles are selected
  if (edit_item->selected_vertices().isEmpty()) return;

  //if ( edit_item->usage_scenario() == 0 )
  {
    preprocess(edit_item);
  }
  /*else
  {
    complete_deform(edit_item);
  }*/
  

}

// Pre-processing of deformation, adaptive with usage_scenario_1
void Polyhedron_demo_edit_polyhedron_plugin::preprocess(Scene_edit_polyhedron_item* edit_item) {

  Polyhedron* polyhedron = edit_item->polyhedron();
  Deform_map::iterator deform_it = deform_map.find(edit_item);
  if(deform_it == deform_map.end())  // First time. Need to create the Deform_mesh object.
  {
    Polyhedron_deformation_data& new_data = deform_map[edit_item];
    Deform_mesh* new_deform = new Deform_mesh(*polyhedron, Vertex_index_map(), Edge_index_map());
    new_data.deform_mesh = new_deform;
    new_data.preprocessed = false;

    new_deform->clear_handles();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_handles())
      new_deform->insert_handle(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_handles())
      new_deform->insert_handle(vh);
    new_deform->clear_roi();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_roi())
      new_deform->insert_roi(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_roi())
      new_deform->insert_roi(vh);
  }

  Polyhedron_deformation_data& data = deform_map[edit_item];
  Deform_mesh* deform = data.deform_mesh;

  // precomputation of Laplacian matrix
  deform->preprocess();
  data.preprocessed = true;

}


// Complete deformation, adaptive with usage_scenario_2
void Polyhedron_demo_edit_polyhedron_plugin::complete_deform(Scene_edit_polyhedron_item* edit_item) {

  Polyhedron* polyhedron = edit_item->polyhedron();

  Deform_map::iterator deform_it = deform_map.find(edit_item);
  if(deform_it == deform_map.end()) return;

  Polyhedron_deformation_data& data = deform_map[edit_item];
  Deform_mesh* deform = data.deform_mesh;

  // precomputation of Laplacian matrix
  deform->preprocess();
  data.preprocessed = true;
  deform->deform();
  

  // signal to the item that it needs to recompute its internal structures
  edit_item->changed(); // that reset the last_position()

  // signal to the scene that the item needs to be redrawn.
  scene->itemChanged(edit_item);

}


// classic ROI + single handle paradigm
void Polyhedron_demo_edit_polyhedron_plugin::usage_scenario_0(Scene_edit_polyhedron_item* edit_item) {

  Polyhedron* polyhedron = edit_item->polyhedron();
 
  const Point& orig = edit_item->original_position();
  const Vector translation_origin = edit_item->current_position() - orig;
  const Point& last = edit_item->last_position();
  const Vector translation_last = edit_item->current_position() - last;

  Polyhedron_deformation_data& data = deform_map[edit_item];
  Deform_mesh* deform = data.deform_mesh;
  if ( translation_origin == Vector(0, 0, 0) && translation_last == Vector(0, 0, 0) )  // vertex selection: reset deform class
  { 
    deform->clear_handles();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_handles())
      deform->insert_handle(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_handles())
      deform->insert_handle(vh);
    deform->clear_roi();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_roi())
      deform->insert_roi(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_roi())
      deform->insert_roi(vh);
    data.preprocessed = false;
  }
  else                  // moving frame: actual deformation
  {
    if ( !data.preprocessed )    // need preprocessing
    {
      // precomputation of Laplacian matrix
      deform->preprocess();
      data.preprocessed = true;
    }

    // -- ACTUAL DEFORMATION --

    (*deform)(translation_origin);
    deform->deform();

    // -- END OF ACTUAL DEFORMATION --
  }
}

// point-based multiple handle manipulation
void Polyhedron_demo_edit_polyhedron_plugin::usage_scenario_1(Scene_edit_polyhedron_item* edit_item) {

  Polyhedron* polyhedron = edit_item->polyhedron();



  const Point& orig = edit_item->original_position();

  Point poi; // AF: This works only if we have a single ROI

  const Vector translation_origin = edit_item->current_position() - orig;
  const Point& last = edit_item->last_position();
  const Vector translation_last = edit_item->current_position() - last;
  Polyhedron_deformation_data& data = deform_map[edit_item];
  Deform_mesh* deform = data.deform_mesh;
  if ( translation_origin == Vector(0, 0, 0) && translation_last == Vector(0, 0, 0) )  // handle selection
  { 
    std::cerr << "reset something" << std::endl;
    deform->clear_handles();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_handles())
      deform->insert_handle(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_handles())
      deform->insert_handle(vh);
    deform->clear_roi();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_roi())
      deform->insert_roi(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_roi())
      deform->insert_roi(vh);
    data.preprocessed = false;
  }
  else                  // moving frame: move new handles
  {
    if ( !data.preprocessed )
    {
      edit_item->setSelectedVector(translation_last);
      edit_item->setSelectedHandlesMoved(true);

      Q_FOREACH(Vertex_handle vh, edit_item->selected_handles())
        data.handle_vectors[vh] = translation_origin;
    }
    else
    {
      // AF: for the rotation+ translation, we have to translate handles by  ORIGIN-poi
      //     make the rotation, translate back, and apply the additional translation.
      Vector vec = edit_item->selected_vector().second - edit_item->selected_vector().first;
      poi =  edit_item->selected_vector().first;
      double scalar = translation_origin*vec / vec.squared_length() /3.0;
      if (scalar > 1) scalar = 1;
      if (scalar < 0) scalar = 0;
   
      std::map<Vertex_handle, Vector>::iterator it = data.handle_vectors.begin();
      while ( it != data.handle_vectors.end() )   // apply scalar factor to each handle region
      {
        qglviewer::Quaternion quat(qglviewer::Vec(scalar * it->second.x(),
                                                  scalar * it->second.y(), 
                                                  scalar * it->second.z()), scalar * 3.14);

        qglviewer::Vec disp(scalar * it->second.x(), scalar * it->second.y(), scalar *  it->second.z());   

        (*deform)( it->first, poi, quat, disp);// 0 for the match
        //        (*deform)( it->first, scalar*(it->second) );
        it++;
      }
      deform->deform();
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

  // do deformation only when handles are selected
  if (edit_item->selected_vertices().isEmpty()) return;

  Polyhedron* polyhedron = edit_item->polyhedron();
  Deform_map::iterator deform_it = deform_map.find(edit_item);
  if(deform_it == deform_map.end())  // First time. Need to create the Deform_mesh object.
  {
    Polyhedron_deformation_data& new_data = deform_map[edit_item];
    Deform_mesh* new_deform = new Deform_mesh(*polyhedron, Vertex_index_map(), Edge_index_map());
    new_data.deform_mesh = new_deform;
    new_data.preprocessed = false;

    new_deform->clear_handles();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_handles())
      new_deform->insert_handle(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_handles())
      new_deform->insert_handle(vh);
    new_deform->clear_roi();
    Q_FOREACH(Vertex_handle vh, edit_item->selected_roi())
      new_deform->insert_roi(vh);
    Q_FOREACH(Vertex_handle vh, edit_item->non_selected_roi())
      new_deform->insert_roi(vh);
  }

  if ( edit_item->usage_scenario() == 0 )
  {
    usage_scenario_0(edit_item);
  }
  else
  {
    usage_scenario_1(edit_item);
  }

  // signal to the item that it needs to recompute its internal structures
  edit_item->changed(); // that reset the last_position()

  // signal to the scene that the item needs to be redrawn.
  scene->itemChanged(edit_item);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin, Polyhedron_demo_edit_polyhedron_plugin)

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
