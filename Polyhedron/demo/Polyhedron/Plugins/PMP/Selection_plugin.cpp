 #include <QtCore/qglobal.h>
#include <QMessageBox>
#include <QInputDialog>

#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include "ui_Selection_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include <QMultipleInputDialog.h>
#include <QDoubleSpinBox>
#include <QCheckBox>

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <Scene.h>
typedef Scene_surface_mesh_item Scene_face_graph_item;

typedef Scene_face_graph_item::Face_graph Face_graph;
typedef boost::property_map<Face_graph,CGAL::vertex_point_t>::type VPmap;

struct Is_terminal
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor , const Graph& )
  {
    return false; // degree(vd,g) != 2; is a bad test in case of parallel edges
  }
};


template <typename Graph>
struct Polyline_visitor
{
  Scene_polylines_item* item;
  const Graph& points_pmap;

  Polyline_visitor(Scene_polylines_item* item_,
                   const Graph& points_property_map)
    : item(item_),
      points_pmap(points_property_map)
  {}

  void start_new_polyline()
  {
    item->polylines.push_back( Scene_polylines_item::Polyline() );
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd)
  {
    item->polylines.back().push_back(points_pmap[vd]);
  }
  void end_polyline(){}
};
using namespace CGAL::Three;
class Polyhedron_demo_selection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "selection_plugin.json")
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90")
public:
  QString nameFilters() const override { return "Selection files(*.selection.txt)"; }
  QString name() const override { return "selection_sm_plugin"; }

  bool canLoad(QFileInfo) const override {
    Scene_item * item = CGAL::Three::Three::scene()->item(
          CGAL::Three::Three::scene()->mainSelectionIndex());
    Scene_facegraph_item* fg_item = qobject_cast<Scene_facegraph_item*>(item);
    if(fg_item)
      return true;
    Scene_polyhedron_selection_item* sel_item =
        qobject_cast<Scene_polyhedron_selection_item*>(item);
    if (sel_item)
      return true;
    return false;
  }

  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override {
      if(fileinfo.suffix().toLower() != "txt")
      {
        ok = false;
        return QList<Scene_item*>();
      }
      // There will be no actual loading at this step.
      Scene_polyhedron_selection_item* item = new Scene_polyhedron_selection_item();
      if(!item->load(fileinfo.filePath().toStdString())) {
          delete item;
        ok = false;
        return QList<Scene_item*>();
      }
      item->setName(fileinfo.baseName());
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>()<<item;
  }

  bool canSave(const CGAL::Three::Scene_item* scene_item) override {
      return qobject_cast<const Scene_polyhedron_selection_item*>(scene_item);
  }
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items) override {
    Scene_item* scene_item = items.front();
      const Scene_polyhedron_selection_item* item = qobject_cast<const Scene_polyhedron_selection_item*>(scene_item);
      if(item == NULL) { return false; }

      bool res = item->save(fileinfo.filePath().toStdString());
      if(res)
        items.pop_front();
      return res;
  }

  bool isDefaultLoader(const Scene_item* item) const override{
    if(qobject_cast<const Scene_polyhedron_selection_item*>(item))
      return true;
    return false;
  }

  bool applicable(QAction* action) const override {
    if(action == actionSelfIntersection)
      return qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()));
    else if(action == actionSelection)
      return qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()))
          || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    return false;
  }
  void print_message(QString message) { CGAL::Three::Three::information(message); }

  QList<QAction*> actions() const override {
    return QList<QAction*>() << actionSelection
                             << actionSelfIntersection;
  }

  using Polyhedron_demo_io_plugin_interface::init;
  virtual void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) override{
    mw = mainWindow;
    scene = scene_interface;
    messages = m;

    actionSelfIntersection = new QAction(tr("Self-&Intersection Test"), mw);
    actionSelfIntersection->setObjectName("actionSelfIntersection");
    actionSelfIntersection->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionSelfIntersection, SIGNAL(triggered()), this, SLOT(on_actionSelfIntersection_triggered()));

    actionSelection = new QAction(
          QString("Surface Mesh Selection")
          , mw);
    actionSelection->setObjectName("actionSelection");
    connect(actionSelection, SIGNAL(triggered()), this, SLOT(selection_action()));
    last_mode = 0;
    dock_widget = new QDockWidget(
                "Surface Mesh Selection"
          , mw);
    dock_widget->setVisible(false);
    ui_widget.setupUi(dock_widget);
    dock_widget->setWindowTitle(tr(
                                  "Surface Mesh Selection"
                                  ));
    connect(dock_widget, &QDockWidget::visibilityChanged,
            this, [this](bool b){
        this->set_highlighting(b);
      if(!b)
        this->set_operation_mode(-1);
    });

    addDockWidget(dock_widget);

    connect(ui_widget.hl_checkBox, SIGNAL(toggled(bool)), this, SIGNAL(set_highlighting(bool)));
    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Select_all_NTButton,  SIGNAL(clicked()), this, SLOT(on_Select_all_NTButton_clicked()));
    connect(ui_widget.Select_boundaryButton,  SIGNAL(clicked()), this, SLOT(on_Select_boundaryButton_clicked()));
    connect(ui_widget.Add_to_selection_button,  SIGNAL(clicked()), this, SLOT(on_Add_to_selection_button_clicked()));
    connect(ui_widget.Regularize_button,  SIGNAL(clicked()), this, SLOT(on_Regularize_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Clear_all_button,  SIGNAL(clicked()), this, SLOT(on_Clear_all_button_clicked()));
    connect(ui_widget.Inverse_selection_button,  SIGNAL(clicked()), this, SLOT(on_Inverse_selection_button_clicked()));
    connect(ui_widget.Select_isolated_components_button,  SIGNAL(clicked()), this, SLOT(on_Select_isolated_components_button_clicked()));
    connect(ui_widget.Get_minimum_button,  SIGNAL(clicked()), this, SLOT(on_Get_minimum_button_clicked()));
    connect(ui_widget.Create_selection_item_button,  SIGNAL(clicked()), this, SLOT(on_Create_selection_item_button_clicked()));
    connect(ui_widget.Selection_type_combo_box, SIGNAL(currentIndexChanged(int)),
            this, SLOT(on_Selection_type_combo_box_changed(int)));
    connect(ui_widget.lassoCheckBox, &QCheckBox::toggled,
            this, &Polyhedron_demo_selection_plugin::on_LassoCheckBox_changed);
    connect(ui_widget.Insertion_radio_button, SIGNAL(toggled(bool)), this, SLOT(on_Insertion_radio_button_toggled(bool)));
    connect(ui_widget.Brush_size_spin_box, SIGNAL(valueChanged(int)), this, SLOT(on_Brush_size_spin_box_changed(int)));
    connect(ui_widget.validateButton, SIGNAL(clicked()), this, SLOT(on_validateButton_clicked()));
    connect(ui_widget.Expand_reduce_button, SIGNAL(clicked()), this, SLOT(on_Expand_reduce_button_clicked()));
    connect(ui_widget.Select_sharp_edges_button, SIGNAL(clicked()), this, SLOT(on_Select_sharp_edges_button_clicked()));
    connect(ui_widget.selectionOrEuler, SIGNAL(currentChanged(int)), this, SLOT(on_SelectionOrEuler_changed(int)));
    connect(ui_widget.editionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_editionBox_changed(int)));

    ui_widget.Add_to_selection_button->hide();
    ui_widget.Select_all_NTButton->hide();
    QObject* scene = dynamic_cast<QObject*>(scene_interface);
    if(scene) {
      connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(item_about_to_be_destroyed(CGAL::Three::Scene_item*)));
      connect(scene, SIGNAL(newItem(int)), this, SLOT(new_item_created(int)));
      connect(scene, SIGNAL(selectionChanged(int)), this, SLOT(filter_operations()));
    }

    //Fill operations combo box.
    operations_strings = {
      "Create Point Set from Selected Vertices"           ,
      "Create Polyline from Selected Edges"               ,
      "Create Facegraph from Selected Facets"            ,
      "Erase Selected Facets"             ,
      "Keep Connected Components of Selected Facets"           ,
      "Expand Face Selection to Stay Manifold After Removal"   ,
      "Convert from Edge Selection to Facets Selection"        ,
      "Convert from Edge Selection to Point Selection"         ,
      "Convert from Facet Selection to Boundary Edge Selection",
      "Convert from Facet Selection to Point Selection"
    };

    operations_map[operations_strings[0]] = 0;
    operations_map[operations_strings[1]] = 1;
    operations_map[operations_strings[2]] = 2;
    operations_map[operations_strings[3]] = 3;
    operations_map[operations_strings[4]] = 4;
    operations_map[operations_strings[5]] = 5;
    operations_map[operations_strings[6]] = 6;
    operations_map[operations_strings[7]] = 7;
    operations_map[operations_strings[8]] = 8;
    operations_map[operations_strings[9]] = 9;
  }
  virtual void closure() override
  {
    dock_widget->hide();
  }
Q_SIGNALS:
  void save_handleType();
  void set_operation_mode(int);
  void set_highlighting(bool);

public Q_SLOTS:

  void on_actionSelfIntersection_triggered();

  void connectItem(Scene_polyhedron_selection_item* new_item)
  {
    connect(this, SIGNAL(save_handleType()),new_item, SLOT(save_handleType()));
    connect(new_item, SIGNAL(updateInstructions(QString)), this, SLOT(setInstructions(QString)));
    connect(this, SIGNAL(set_operation_mode(int)),new_item, SLOT(set_operation_mode(int)));
    connect(this, SIGNAL(set_highlighting(bool)),new_item, SLOT(set_highlighting(bool)));
    this->set_highlighting(ui_widget.hl_checkBox->isChecked());
    int item_id = scene->addItem(new_item);
   // QObject* scene_ptr = dynamic_cast<QObject*>(scene);
   // if (scene_ptr)
   //   connect(new_item,SIGNAL(simplicesSelected(CGAL::Three::Scene_item*)), scene_ptr, SLOT(setSelectedItem(CGAL::Three::Scene_item*)));
    connect(new_item,SIGNAL(isCurrentlySelected(Scene_facegraph_item_k_ring_selection*)), this, SLOT(isCurrentlySelected(Scene_facegraph_item_k_ring_selection*)));
    connect(new_item,SIGNAL(simplicesSelected(CGAL::Three::Scene_item*)), this, SLOT(filter_operations()));
    scene->setSelectedItem(item_id);
    on_SelectionOrEuler_changed(ui_widget.selectionOrEuler->currentIndex());
    if(last_mode == 0)
      on_Selection_type_combo_box_changed(ui_widget.Selection_type_combo_box->currentIndex());
    filter_operations();
  }
  // If the selection_item or the polyhedron_item associated to the k-ring_selector is currently selected,
  // set the k-ring_selector as currently selected. (A k-ring_selector tha tis not "currently selected" will
  // not process selection events)
  void isCurrentlySelected(Scene_facegraph_item_k_ring_selection* item)
  {
    if(scene->item_id(selection_item_map.find(item->poly_item)->second) == scene->mainSelectionIndex() ||
       scene->item_id(item->poly_item)== scene->mainSelectionIndex() )
      item->setCurrentlySelected(true);
    else
      item->setCurrentlySelected(false);
  }

  void setInstructions(QString s)
  {
    ui_widget.instructionsLabel->setText(s);
  }

  void printMessage(QString s)
  {
    print_message(s);
  }

  void selection_action() {
    dock_widget->show();
    dock_widget->raise();
    Scene_face_graph_item* poly_item = getSelectedItem<Scene_face_graph_item>();
    if(!poly_item || selection_item_map.find(poly_item) != selection_item_map.end()) { return; }
    Scene_polyhedron_selection_item* new_item = new Scene_polyhedron_selection_item(poly_item, mw);
    new_item->setName(QString("%1 (selection)").arg(poly_item->name()));
    connectItem(new_item);
  }

  Scene_polyhedron_selection_item* onTheFlyItem() {
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()));
    if(!poly_item)
      return NULL;
    Scene_polyhedron_selection_item* new_item = new Scene_polyhedron_selection_item(poly_item, mw);
    new_item->setName(QString("%1 (selection)").arg(poly_item->name()));
    connectItem(new_item);
    return new_item;

  }
  // Select all
  void on_Select_all_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      selection_item = onTheFlyItem();
    if(!selection_item)
    {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->select_all();
    filter_operations();
  }

  void on_Select_all_NTButton_clicked()
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      selection_item = onTheFlyItem();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->select_all_NT();
    filter_operations();
  }

  // Select Boundary
  void on_Select_boundaryButton_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      selection_item = onTheFlyItem();
    if(!selection_item)
    {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->select_boundary();
    filter_operations();
  }

  void on_Add_to_selection_button_clicked()
  {

    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->add_to_selection();
    filter_operations();
  }
  // Regularize selection using graph cut
  void on_Regularize_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    QMultipleInputDialog dialog ("Regularize Selection Border", mw);
    QDoubleSpinBox* weight = dialog.add<QDoubleSpinBox> ("Weight:");
    weight->setDecimals (6);
    weight->setRange (0., 0.999999); // [0;1[
    weight->setValue (0.5);

    QCheckBox* prevent_unselection = dialog.add<QCheckBox> ("Prevent unselection");

    if (dialog.exec() != QDialog::Accepted)
      return;

    boost::unordered_map<fg_face_descriptor, bool> is_selected_map;
    for(fg_face_descriptor fh : faces(*selection_item->polyhedron()))
    {
      if(selection_item->selected_facets.find(fh)
         == selection_item->selected_facets.end())
        is_selected_map[fh]=false;
      else
      {
        is_selected_map[fh]=true;
      }
    }

    auto border_length =
      [&]() -> double
      {
        double out = 0.;
        for(Scene_polyhedron_selection_item::fg_edge_descriptor ed : edges(*selection_item->polyhedron()))
        {
          fg_face_descriptor f0 = face (halfedge (ed, *selection_item->polyhedron()),
                                        *selection_item->polyhedron());
          fg_face_descriptor f1 = face (opposite(halfedge (ed, *selection_item->polyhedron()),
                                                 *selection_item->polyhedron()),
                                        *selection_item->polyhedron());
          if (is_selected_map[f0] == is_selected_map[f1])
            continue;

          fg_vertex_descriptor esource = source(ed, *selection_item->polyhedron());
          fg_vertex_descriptor etarget = target(ed, *selection_item->polyhedron());

          out += std::sqrt(CGAL::squared_distance (get (get(CGAL::vertex_point,*selection_item->polyhedron()), esource),
                                                   get (get(CGAL::vertex_point,*selection_item->polyhedron()), etarget)));
        }
        return out;
      };

    std::cerr << "[Selection Regularization] Weight = " << weight->value() << std::endl;

    std::cerr << "Length of border before regularization = " << border_length() << std::endl;

    CGAL::regularize_face_selection_borders (*selection_item->polyhedron(),
                                             boost::make_assoc_property_map(is_selected_map),
                                             weight->value(),
                                             CGAL::parameters::prevent_unselection (prevent_unselection->isChecked()));

    std::cerr << "Length of border after regularization = " << border_length() << std::endl;

    selection_item->selected_facets.clear();

    for(fg_face_descriptor fh : faces(*selection_item->polyhedron()))
    {
      if (is_selected_map[fh])
        selection_item->selected_facets.insert(fh);
    }
    selection_item->invalidateOpenGLBuffers();
    selection_item->itemChanged();
  }
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }
    if(ui_widget.tabWidget->currentIndex() == 0)
      selection_item->clear();
    else //case current tab is Components, then active simplices are edges.
    {
      selection_item->selected_edges.clear();
      selection_item->invalidateOpenGLBuffers();
      selection_item->itemChanged();
    }
    filter_operations();
  }
  void on_Clear_all_button_clicked(){
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->clear_all();
    filter_operations();
  }
  void on_Inverse_selection_button_clicked()
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }
    selection_item->inverse_selection();
    filter_operations();
  }
  // Isolated component related functions
  void on_Select_isolated_components_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      selection_item = onTheFlyItem();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    boost::optional<std::size_t> minimum =
      selection_item->select_isolated_components(ui_widget.Threshold_size_spin_box->value());
    if(minimum) {
      ui_widget.Threshold_size_spin_box->setValue((int) *minimum);
    }
    filter_operations();
  }
  void on_Get_minimum_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      selection_item = onTheFlyItem();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }
    boost::optional<std::size_t> minimum = selection_item->get_minimum_isolated_component();
    if(minimum) {
      ui_widget.Threshold_size_spin_box->setValue((int) *minimum);
    }
    filter_operations();
  }
  // Create selection item for selected polyhedron item
  void on_Create_selection_item_button_clicked() {

    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()));
    if(!poly_item) {
      print_message("Error: there is no selected "
                    "Surface_mesh "
                    "item!");
      return;
    }
    if(selection_item_map.find(poly_item) != selection_item_map.end())
      return;
    // all other arrangements (putting inside selection_item_map), setting names etc,
    // other params (e.g. k_ring) will be set inside new_item_created
    from_plugin = true;
    Scene_polyhedron_selection_item* new_item = new Scene_polyhedron_selection_item(poly_item, mw);
    new_item->setName(QString("%1 (selection)").arg(poly_item->name()));
    ui_widget.selectionOrEuler->setCurrentIndex(last_mode);
    connectItem(new_item);
    filter_operations();
  }
  void on_LassoCheckBox_changed(bool b)
  {
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it)
    {
      it->second->set_lasso_mode(b);
    }
  }
  void on_Selection_type_combo_box_changed(int index) {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_active_handle_type(static_cast<Active_handle::Type>(index));
      Q_EMIT save_handleType();
      switch(index)
      {
      case 0:
      case 1:
      case 2:
        ui_widget.lassoCheckBox->show();
        break;
      default:
        ui_widget.lassoCheckBox->hide();
        ui_widget.lassoCheckBox->setChecked(false);
        it->second->set_lasso_mode(false);
        break;
      }
      if(index == 1)
      {
        ui_widget.Select_all_NTButton->show();
        ui_widget.Add_to_selection_button->hide();
        ui_widget.Select_boundaryButton->hide();
        Q_EMIT set_operation_mode(-1);
      }
      else if(index == 2)
      {
        ui_widget.Select_all_NTButton->hide();
        ui_widget.Add_to_selection_button->hide();
        ui_widget.Select_boundaryButton->show();
        Q_EMIT set_operation_mode(-1);
      }
      else if(index == 4)
      {
        it->second->setPathSelection(true);
        ui_widget.Select_all_NTButton->hide();
        ui_widget.Add_to_selection_button->show();
        ui_widget.Select_boundaryButton->show();
        Q_EMIT set_operation_mode(-2);
      }
      else
      {
        ui_widget.Add_to_selection_button->hide();
        ui_widget.Select_all_NTButton->hide();
        ui_widget.Select_boundaryButton->hide();
        it->second->setPathSelection(false);
        Q_EMIT set_operation_mode(-1);
      }
    }
    filter_operations();
  }

  void on_Insertion_radio_button_toggled(bool toggle){
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_is_insert(toggle);
    }
  }
  void on_Brush_size_spin_box_changed(int value) {
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_k_ring(value);
    }
  }

  void on_validateButton_clicked() {
    switch(operations_map[ui_widget.operationsBox->currentText()])
    {
    //Create Point Set Item from Selected Vertices
    case 0:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      if(selection_item->selected_vertices.empty()) {
        print_message("Error: there is no selected vertex in polyhedron selection item!");
        return;
      }
      Scene_points_with_normal_item* point_item = new Scene_points_with_normal_item();
      point_item->setName(QString("%1-points").arg(selection_item->name()));
      Face_graph *poly = selection_item->polyhedron();
      VPmap vpm = get(CGAL::vertex_point,*poly);

      for(Scene_polyhedron_selection_item::Selection_set_vertex::iterator begin = selection_item->selected_vertices.begin();
         begin != selection_item->selected_vertices.end(); ++begin) {
        point_item->point_set()->insert(get(vpm,*begin));
      }
      scene->setSelectedItem( scene->addItem(point_item) );
      break;
    }
      //Create Polyline Item from Selected Edges
    case 1:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      if(selection_item->selected_edges.empty()) {
        print_message("Error: there is no selected edge in polyhedron selection item!");
        return;
      }
      Scene_polylines_item* polyline_item = new Scene_polylines_item();
      polyline_item->setName(QString("%1-edges").arg(selection_item->name()));

      typedef boost::adjacency_list < boost::listS,
          boost::vecS,
          boost::undirectedS,
          Kernel::Point_3 > Edge_graph;

      Edge_graph edge_graph;
      std::map<fg_vertex_descriptor, Edge_graph::vertex_descriptor> p2vd;
      std::map<fg_vertex_descriptor, Edge_graph::vertex_descriptor>::iterator it_find;
      bool insert_OK;

      Face_graph * poly = selection_item->polyhedron();
      VPmap vpm = get(CGAL::vertex_point,*poly);
      for(Scene_polyhedron_selection_item::Selection_set_edge::iterator begin = selection_item->selected_edges.begin();
          begin != selection_item->selected_edges.end(); ++begin)
      {
        fg_vertex_descriptor source = target(opposite(halfedge(*begin,*poly),*poly),*poly);
        boost::tie(it_find, insert_OK)
            = p2vd.insert(std::make_pair(source, Edge_graph::vertex_descriptor()));
        if (insert_OK)
        {
          it_find->second = add_vertex(edge_graph);
          edge_graph[it_find->second] = get(vpm,source);
        }
        Edge_graph::vertex_descriptor src=it_find->second;

        fg_vertex_descriptor targ = target(halfedge(*begin,*poly),*poly);
        boost::tie(it_find, insert_OK)
            = p2vd.insert(std::make_pair(targ, Edge_graph::vertex_descriptor()));
        if (insert_OK)
        {
          it_find->second = add_vertex(edge_graph);
          edge_graph[it_find->second] = get(vpm,targ);
        }
        Edge_graph::vertex_descriptor tgt=it_find->second;
        boost::add_edge(src, tgt, edge_graph);
      }


      Polyline_visitor<Edge_graph> polyline_visitor(polyline_item, edge_graph);
      CGAL::split_graph_into_polylines( edge_graph,
                                        polyline_visitor,
                                        Is_terminal() );
      scene->setSelectedItem( scene->addItem(polyline_item) );
      scene->itemChanged(polyline_item);
      break;
    }
      //Create Polyhedron Item from Selected Facets
    case 2:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }

      Scene_face_graph_item* poly_item = new Scene_face_graph_item();
      if(selection_item->export_selected_facets_as_polyhedron(poly_item->polyhedron())) {
        poly_item->setName(QString("%1-facets").arg(selection_item->name()));
        poly_item->invalidateOpenGLBuffers(); // for init()
        scene->setSelectedItem( scene->addItem(poly_item) );
        scene->itemChanged(poly_item);
      }
      else {
        delete poly_item;
        print_message("Error: polyhedron item is not created!");
      }
      break;
    }
      //Erase Selected Facets from Polyhedron Item
    case 3:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }

      selection_item->erase_selected_facets();
      break;
    }
      //Keep connected components of Selected Facets
    case 4:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if (!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      selection_item->keep_connected_components();
      break;
    }
      //Expand face selection
    case 5:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if (!selection_item ||
          selection_item->selected_facets.empty())
      {
        print_message("Error: Please select a selection item with a selection of faces.");
        return;
      }
      boost::unordered_map<fg_face_descriptor, bool> is_selected_map;
      int index = 0;
      for(fg_face_descriptor fh : faces(*selection_item->polyhedron()))
      {
        if(selection_item->selected_facets.find(fh)
           == selection_item->selected_facets.end())
          is_selected_map[fh]=false;
        else
        {
          is_selected_map[fh]=true;
        }
        ++index;
      }
      CGAL::expand_face_selection_for_removal(selection_item->selected_facets,
                                              *selection_item->polyhedron(),
                                        boost::make_assoc_property_map(is_selected_map));

      for(fg_face_descriptor fh : faces(*selection_item->polyhedron()))
      {
        if (is_selected_map[fh])
          selection_item->selected_facets.insert(fh);
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->itemChanged();
      break;
    }
      //Convert from Edge Selection to Facet Selection
    case 6:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      if(selection_item->selected_edges.empty()) {
        print_message("Error: there is no selected edge in this polyhedron selection item!");
        return;
      }
      const Face_graph& poly = *selection_item->polyhedron();
      for(Scene_polyhedron_selection_item::fg_edge_descriptor ed : selection_item->selected_edges)
      {
        if(!is_border(halfedge(ed,poly), poly)){
          selection_item->selected_facets.insert(face(halfedge(ed, poly), poly));
        }
        if(!is_border(opposite(halfedge(ed,poly), poly), poly)){
          selection_item->selected_facets.insert(face(opposite(halfedge(ed, poly), poly), poly));
        }
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->itemChanged();
      break;
    }
      //Convert from Edge Selection to Point Selection
    case 7:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      if(selection_item->selected_edges.empty()) {
        print_message("Error: there is no selected edge in this polyhedron selection item!");
        return;
      }
      const Face_graph& poly = *selection_item->polyhedron();

      for(Scene_polyhedron_selection_item::fg_edge_descriptor ed : selection_item->selected_edges)
      {
        selection_item->selected_vertices.insert(target(halfedge(ed, poly), poly));
        selection_item->selected_vertices.insert(source(halfedge(ed, poly), poly));
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->itemChanged();
      break;
    }
      //Convert from Facet Selection to Bounding Edge Selection
    case 8:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      if(selection_item->selected_facets.empty()) {
        print_message("Error: there is no selected facet in this polyhedron selection item!");
        return;
      }
      const Face_graph& poly = *selection_item->polyhedron();
      std::vector<Scene_polyhedron_selection_item::fg_halfedge_descriptor> boundary_edges;
      CGAL::Polygon_mesh_processing::border_halfedges(selection_item->selected_facets, poly, std::back_inserter(boundary_edges));
      for(Scene_polyhedron_selection_item::fg_halfedge_descriptor h : boundary_edges)
      {
        selection_item->selected_edges.insert(edge(h, poly));
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->itemChanged();
      break;
    }
      //Convert from Facet Selection to Points Selection
    case 9:
    {
      Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
      if(!selection_item) {
        print_message("Error: there is no selected polyhedron selection item!");
        return;
      }
      if(selection_item->selected_facets.empty()) {
        print_message("Error: there is no selected facet in this polyhedron selection item!");
        return;
      }
      const Face_graph& poly = *selection_item->polyhedron();
      for(Scene_polyhedron_selection_item::fg_face_descriptor fh : selection_item->selected_facets)
      {
        for(Scene_polyhedron_selection_item::fg_halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(fh,poly), poly) )
        {
          selection_item->selected_vertices.insert(target(h, poly));
        }
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->itemChanged();
      break;
    }
    default :
      break;
    }
    filter_operations();
    return;
  }

  void on_SelectionOrEuler_changed(int index)
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      return;
    else
      selection_item->on_Ctrlz_pressed();
    last_mode = index;
    switch(index)
    {
    //Selection mode
    case 0:
      Q_EMIT set_operation_mode(-1);
      on_Selection_type_combo_box_changed(ui_widget.Selection_type_combo_box->currentIndex());
      break;
      //Edition mode
    case 1:
    {
      bool is_valid = true;
      for(boost::graph_traits<Face_graph>::face_descriptor fd : faces(*selection_item->polyhedron()))
      {
        if (CGAL::is_triangle(halfedge(fd, *selection_item->polyhedron()), *selection_item->polyhedron())
            && CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(fd, *selection_item->polyhedron()))
        {
          is_valid = false;
          break;
        }
      }
      if(!is_valid)
      {
        QMessageBox::warning(mw,
                             tr("Degenerated Face_graph"),
                             tr("Degenerated faces have been detected. Problems may occur "
                                "for operations other than \"Move point\". "));
      }
      //remove lasso mode
      selection_item->set_lasso_mode(false);
      ui_widget.lassoCheckBox->setChecked(false);
      Q_EMIT save_handleType();
      on_editionBox_changed(ui_widget.editionBox->currentIndex());
      break;
    }
    }
  }

  void on_editionBox_changed(int mode )
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(selection_item)
    {
      selection_item->on_Ctrlz_pressed();
    }
    if(ui_widget.selectionOrEuler->currentIndex() == 0)
    {
      Q_EMIT set_operation_mode(-1);
    }
    else
    {
      Q_EMIT set_operation_mode(mode);
    }
    switch(mode)
    {
    //Join vertex
    case 0:
    //Split vertex
    case 1:
    {
      QPixmap pm(":/cgal/Polyhedron_3/resources/euler_vertex.png");
      ui_widget.docImage_Label->setPixmap(pm);
      break;
    }
    //Join face
    case 3:
    //Split face
    case 4:
    {
      QPixmap pm(":/cgal/Polyhedron_3/resources/euler_facet.png");
      ui_widget.docImage_Label->setPixmap(pm);
      break;
    }
    //Collapse edge
    case 5:
    {
      QPixmap pm(":/cgal/Polyhedron_3/resources/general_collapse.png");
      ui_widget.docImage_Label->setPixmap(pm);
      break;
    }
    //Add center vertex
    case 7:
    //Remove center vertex
    case 8:
    {
      QPixmap pm(":/cgal/Polyhedron_3/resources/euler_center.png");
      ui_widget.docImage_Label->setPixmap(pm);
      break;
    }
    //Add vertex and face to border
    case 9:
    {
      QPixmap pm(":/cgal/Polyhedron_3/resources/add_facet1.png");
      ui_widget.docImage_Label->setPixmap(pm);
      break;
    }
    //add facet to border
    case 10:
    {
      QPixmap pm(":/cgal/Polyhedron_3/resources/add_facet2.png");
      ui_widget.docImage_Label->setPixmap(pm);
      break;
    }
    default:
      ui_widget.docImage_Label->clear();
      break;
    }
    on_LassoCheckBox_changed(ui_widget.lassoCheckBox->isChecked());
  }
  void on_Select_sharp_edges_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item)
      selection_item = onTheFlyItem();
    if (!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    double angle = ui_widget.Sharp_angle_spinbox->value();
    selection_item->select_sharp_edges(angle);
    scene->itemChanged(selection_item);
    filter_operations();
  }

  void on_Expand_reduce_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    int steps = ui_widget.Expand_reduce_spin_box->value();
    selection_item->expand_or_reduce(steps);
    filter_operations();
  }
  // To handle empty selection items coming from loader
  void new_item_created(int item_id) {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(item_id));
    if(!selection_item) { return; }

    Scene_face_graph_item* poly_item = NULL;
    if(selection_item->polyhedron_item() == NULL) { //coming from selection_io loader
      bool found = false;
      for(int i = 0; i<scene->numberOfEntries(); ++i){
        poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(i));
        if(!poly_item)
          continue;
        if(!selection_item->actual_load(poly_item, mw)) {
          continue;
        }
        found = true;
        selection_item->invalidateOpenGLBuffers();
        scene->itemChanged(selection_item);
        break;
      }
      if(!found)
      {
        print_message("Error: loading selection item is not successful!");
        scene->erase(item_id);
        return;
      }
    } else {
      poly_item = getSelectedItem<Scene_face_graph_item>();
    }
    if(!poly_item) {
      print_message("Error: please select corresponding polyhedron item from Geometric Objects list.");
      scene->erase(item_id);
      return;
    }

    // now set default params both for selection items coming from selection_io, or on_Create_selection_item_button_clicked
    Active_handle::Type aht = static_cast<Active_handle::Type>(ui_widget.Selection_type_combo_box->currentIndex());
    bool is_insert = ui_widget.Insertion_radio_button->isChecked();
    int k_ring = ui_widget.Brush_size_spin_box->value();

    selection_item->set_active_handle_type(aht);
    selection_item->set_is_insert(is_insert);
    selection_item->set_k_ring(k_ring);
    selection_item->setRenderingMode(Flat);
    if(selection_item->name() == "unamed") {
      selection_item->setName(tr("%1 (selection)").arg(poly_item->name()));
    }

    selection_item_map.insert(std::make_pair(poly_item, selection_item));
    connect(this, SIGNAL(save_handleType()),selection_item, SLOT(save_handleType()));
    connect(selection_item, SIGNAL(updateInstructions(QString)), this, SLOT(setInstructions(QString)));
    connect(selection_item, SIGNAL(printMessage(QString)), this, SLOT(printMessage(QString)));
    connect(this, SIGNAL(set_operation_mode(int)),selection_item, SLOT(set_operation_mode(int)));
    connect(this, SIGNAL(set_highlighting(bool)),selection_item, SLOT(set_highlighting(bool)));
    this->set_highlighting(ui_widget.hl_checkBox->isChecked());
    //QObject* scene_ptr = dynamic_cast<QObject*>(scene);
    //if (scene_ptr)
    //  connect(selection_item,SIGNAL(simplicesSelected(CGAL::Three::Scene_item*)), scene_ptr, SLOT(setSelectedItem(CGAL::Three::Scene_item*)));
    connect(selection_item,SIGNAL(isCurrentlySelected(Scene_facegraph_item_k_ring_selection*)), this, SLOT(isCurrentlySelected(Scene_facegraph_item_k_ring_selection*)));
    on_LassoCheckBox_changed(ui_widget.lassoCheckBox->isChecked());

    if(!from_plugin){
      ui_widget.selectionOrEuler->setCurrentIndex(0);
    }
    else
      from_plugin = false;
    on_SelectionOrEuler_changed(ui_widget.selectionOrEuler->currentIndex());
    if(last_mode == 0)
      on_Selection_type_combo_box_changed(ui_widget.Selection_type_combo_box->currentIndex());
  }
  void item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item) {
    // if polyhedron item
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene_item);
    if(poly_item) {
      std::pair<Selection_item_map::iterator, Selection_item_map::iterator> res =
        selection_item_map.equal_range(poly_item);

      for(Selection_item_map::iterator begin = res.first; begin != res.second; ) {
        Scene_polyhedron_selection_item* selection_item = begin->second;
        selection_item_map.erase(begin++); // first erase from map, because scene->erase will cause a call to this function
        scene->erase( scene->item_id(selection_item) );
      }
    }
    // if polyhedron selection item
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene_item);
    if(selection_item) {
      Scene_face_graph_item* poly_item = selection_item->polyhedron_item();
      std::pair<Selection_item_map::iterator, Selection_item_map::iterator> res =
        selection_item_map.equal_range(poly_item);
      for(Selection_item_map::iterator begin = res.first; begin != res.second; ++begin) {
        if(begin->second == selection_item) {
          selection_item_map.erase(begin); break;
        }
      }
    }
  }
void filter_operations()
{
  Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(
                                                       scene->mainSelectionIndex()));
  if (!selection_item)
    return;
  QString current_op = ui_widget.operationsBox->currentText();
  ui_widget.operationsBox->clear();

  bool has_v(!selection_item->selected_vertices.empty() ||
             ui_widget.Selection_type_combo_box->currentIndex() == 0),
      has_e(!selection_item->selected_edges.empty()||
            ui_widget.Selection_type_combo_box->currentIndex() == 2||
            ui_widget.Selection_type_combo_box->currentIndex() == 4),
      has_f(!selection_item->selected_facets.empty()||
            ui_widget.Selection_type_combo_box->currentIndex() == 1||
            ui_widget.Selection_type_combo_box->currentIndex() == 3);

  if(has_v)
  {
    ui_widget.operationsBox->addItem(operations_strings[0]);
  }
  if(has_e)
  {
    ui_widget.operationsBox->addItem(operations_strings[1]);
    ui_widget.operationsBox->addItem(operations_strings[6]);
    ui_widget.operationsBox->addItem(operations_strings[7]);
  }
  if(has_f)
  {
    ui_widget.operationsBox->addItem(operations_strings[2]);
    ui_widget.operationsBox->addItem(operations_strings[3]);
    ui_widget.operationsBox->addItem(operations_strings[4]);
    ui_widget.operationsBox->addItem(operations_strings[5]);
    ui_widget.operationsBox->addItem(operations_strings[8]);
    ui_widget.operationsBox->addItem(operations_strings[9]);
  }
  if(!current_op.isEmpty())
    ui_widget.operationsBox->setCurrentText(current_op);
}
private:
  Messages_interface* messages;
  QAction* actionSelection;
  QAction *actionSelfIntersection;

  QDockWidget* dock_widget;
  Ui::Selection ui_widget;
  std::map<QString, int> operations_map;
  std::vector<QString> operations_strings;
typedef boost::unordered_map<Scene_face_graph_item*, Scene_polyhedron_selection_item*> Selection_item_map;
  Selection_item_map selection_item_map;
  int last_mode;
  bool from_plugin;
}; // end Polyhedron_demo_selection_plugin


template<class Mesh>
bool selfIntersect(Mesh* mesh, std::vector<std::pair<typename boost::graph_traits<Mesh>::face_descriptor,typename boost::graph_traits<Mesh>::face_descriptor> > &faces)
{
  if(!CGAL::is_triangle_mesh(*mesh))
  {
    CGAL::Three::Three::warning("%1 skipped because not triangulated.");
    return false;
  }
  // compute self-intersections
  CGAL::Polygon_mesh_processing::self_intersections
    (*mesh, std::back_inserter(faces),
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *mesh)));

  std::cout << "ok (" << faces.size() << " triangle pair(s))" << std::endl;
  return !faces.empty();
}

void Polyhedron_demo_selection_plugin::on_actionSelfIntersection_triggered()
{
  typedef boost::graph_traits<Face_graph>::face_descriptor Face_descriptor;
  typedef boost::graph_traits<Face_graph>::halfedge_descriptor halfedge_descriptor;
  QCursor tmp_cursor(Qt::WaitCursor);
  CGAL::Three::Three::CursorScopeGuard guard(tmp_cursor);
  bool found = false;
  std::vector<Scene_face_graph_item*> selected_polys;
  Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_face_graph_item* poly_item =
        qobject_cast<Scene_face_graph_item*>(scene->item(index));
    if(poly_item)
    {
      selected_polys.push_back(poly_item);
    }
  }
  Q_FOREACH(Scene_face_graph_item* poly_item, selected_polys)
  {
    Face_graph* mesh = poly_item->face_graph();
    std::vector<std::pair<Face_descriptor, Face_descriptor> > faces;
    // add intersecting triangles to a new Surface_mesh.
    if(selfIntersect(mesh, faces))
    {
      Scene_polyhedron_selection_item* selection_item = nullptr;
      bool should_add = true;
      if(selection_item_map.find(poly_item) != selection_item_map.end())
      {
        QApplication::restoreOverrideCursor();
        if(QMessageBox::question(mw, "Question", "Only one Selection Item can be associated to an item at once, "
                                 "and one already exists. Would you like to replace it ? (If not, this item will be skipped.)")
           == QMessageBox::Yes)
        {
          selection_item = selection_item_map.find(poly_item)->second;
          selection_item->clear();
          should_add = false;
        }
        else
          continue;
      }
      else
      {
        selection_item = new Scene_polyhedron_selection_item(poly_item, mw);
      }
      QApplication::setOverrideCursor(Qt::WaitCursor);
      //add the faces
      for(std::vector<std::pair<Face_descriptor, Face_descriptor> >::iterator fb = faces.begin();
          fb != faces.end(); ++fb) {
        selection_item->selected_facets.insert(fb->first);
        selection_item->selected_facets.insert(fb->second);

        //add the edges
        for(halfedge_descriptor he_circ : halfedges_around_face( halfedge(fb->first, *mesh), *mesh))
        {
          selection_item->selected_edges.insert(edge(he_circ, *mesh));
        }
        for(halfedge_descriptor he_circ : halfedges_around_face( halfedge(fb->second, *mesh), *mesh))
        {
          selection_item->selected_edges.insert(edge(he_circ, *mesh));
        }
      }

      selection_item->invalidateOpenGLBuffers();
      selection_item->setName(tr("%1 (selection) (intersecting triangles)").arg(poly_item->name()));
      if(should_add)
        connectItem(selection_item);
      poly_item->setRenderingMode(Wireframe);

      scene->itemChanged(poly_item);

      found = true;
    }
  }
  QApplication::restoreOverrideCursor();
  if(!found)
    QMessageBox::information(mw, tr("No self intersection"),
                             tr("None of the selected surfaces self-intersect."));
}
//Q_EXPORT_PLUGIN2(Polyhedron_demo_selection_plugin, Polyhedron_demo_selection_plugin)

#include "Selection_plugin.moc"
