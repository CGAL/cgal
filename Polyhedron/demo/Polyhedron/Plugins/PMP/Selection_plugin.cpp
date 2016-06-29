#include <QtCore/qglobal.h>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "ui_Selection_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <Scene.h>

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
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
        || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex())); 
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionSelection; }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionSelection = new QAction(tr("Selection"), mw);
    connect(actionSelection, SIGNAL(triggered()), this, SLOT(selection_action()));
    last_mode = 0;
    dock_widget = new QDockWidget("Selection", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);
    Scene* true_scene = static_cast<Scene*>(scene);
    if(true_scene)
      connect(true_scene, SIGNAL(selectionChanged(int)), this, SLOT(selectionChanged(int)));

    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Select_all_NTButton,  SIGNAL(clicked()), this, SLOT(on_Select_all_NTButton_clicked()));
    connect(ui_widget.Add_to_selection_button,  SIGNAL(clicked()), this, SLOT(on_Add_to_selection_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Clear_all_button,  SIGNAL(clicked()), this, SLOT(on_Clear_all_button_clicked()));
    connect(ui_widget.Inverse_selection_button,  SIGNAL(clicked()), this, SLOT(on_Inverse_selection_button_clicked()));
    connect(ui_widget.Select_isolated_components_button,  SIGNAL(clicked()), this, SLOT(on_Select_isolated_components_button_clicked()));
    connect(ui_widget.Get_minimum_button,  SIGNAL(clicked()), this, SLOT(on_Get_minimum_button_clicked()));
    connect(ui_widget.Create_selection_item_button,  SIGNAL(clicked()), this, SLOT(on_Create_selection_item_button_clicked()));    
    connect(ui_widget.Selection_type_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_type_combo_box_changed(int)));
    connect(ui_widget.Insertion_radio_button, SIGNAL(toggled(bool)), this, SLOT(on_Insertion_radio_button_toggled(bool)));
    connect(ui_widget.Brush_size_spin_box, SIGNAL(valueChanged(int)), this, SLOT(on_Brush_size_spin_box_changed(int)));
    connect(ui_widget.validateButton, SIGNAL(clicked()), this, SLOT(on_validateButton_clicked()));
    connect(ui_widget.Expand_reduce_button, SIGNAL(clicked()), this, SLOT(on_Expand_reduce_button_clicked()));
    connect(ui_widget.Select_sharp_edges_button, SIGNAL(clicked()), this, SLOT(on_Select_sharp_edges_button_clicked()));
    connect(ui_widget.modeBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_ModeBox_changed(int)));
    connect(ui_widget.editionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_editionBox_changed(int)));

    ui_widget.Add_to_selection_button->hide();
    ui_widget.Select_all_NTButton->hide();
    QObject* scene = dynamic_cast<QObject*>(scene_interface);
    if(scene) { 
      connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(item_about_to_be_destroyed(CGAL::Three::Scene_item*)));
      connect(scene, SIGNAL(newItem(int)), this, SLOT(new_item_created(int)));
    } 
  }
  virtual void closure()
  {
    dock_widget->hide();
  }
Q_SIGNALS:
  void save_handleType();
  void set_operation_mode(int);
public Q_SLOTS:
  //If the mainSelectedItem is a selection_item, disable the picking item selection. Else, enable it.
  void selectionChanged(int i)
  {
    QGLViewer* v = *QGLViewer::QGLViewerPool().begin();
    CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
    if(!viewer)
        return;
    Scene_polyhedron_selection_item* current_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(i));
    if(!current_item)
      viewer->setBindingSelect();
    else
      viewer->setNoBinding();
  }

  void connectItem(Scene_polyhedron_selection_item* new_item)
  {
    connect(this, SIGNAL(save_handleType()),new_item, SLOT(save_handleType()));
    connect(new_item, SIGNAL(updateInstructions(QString)), this, SLOT(setInstructions(QString)));
    connect(this, SIGNAL(set_operation_mode(int)),new_item, SLOT(set_operation_mode(int)));
    int item_id = scene->addItem(new_item);
    QObject* scene_ptr = dynamic_cast<QObject*>(scene);
    if (scene_ptr)
      connect(new_item,SIGNAL(simplicesSelected(CGAL::Three::Scene_item*)), scene_ptr, SLOT(setSelectedItem(CGAL::Three::Scene_item*)));
    connect(new_item,SIGNAL(isCurrentlySelected(Scene_polyhedron_item_k_ring_selection*)), this, SLOT(isCurrentlySelected(Scene_polyhedron_item_k_ring_selection*)));
    scene->setSelectedItem(item_id);
    on_ModeBox_changed(ui_widget.modeBox->currentIndex());
    if(last_mode == 0)
      on_Selection_type_combo_box_changed(ui_widget.Selection_type_combo_box->currentIndex());
  }
  // If the selection_item or the polyhedron_item associated to the k-ring_selector is currently selected,
  // set the k-ring_selector as currently selected. (A k-ring_selector tha tis not "currently selected" will
  // not process selection events)
  void isCurrentlySelected(Scene_polyhedron_item_k_ring_selection* item)
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

  void selection_action() {
    dock_widget->show();
    dock_widget->raise();
    Scene_polyhedron_item* poly_item = getSelectedItem<Scene_polyhedron_item>();
    if(!poly_item || selection_item_map.find(poly_item) != selection_item_map.end()) { return; }
    Scene_polyhedron_selection_item* new_item = new Scene_polyhedron_selection_item(poly_item, mw);
    connectItem(new_item);
  }

  Scene_polyhedron_selection_item* onTheFlyItem() {
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
    if(!poly_item)
      return NULL;
    Scene_polyhedron_selection_item* new_item = new Scene_polyhedron_selection_item(poly_item, mw);
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
  }

  void on_Add_to_selection_button_clicked()
  {

    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->add_to_selection();
  }
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }

    selection_item->clear();
  }
  void on_Clear_all_button_clicked(){
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }

    selection_item->clear_all();
  }
  void on_Inverse_selection_button_clicked()
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return;
    }
    selection_item->inverse_selection();
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
  }
  // Create selection item for selected polyhedron item
  void on_Create_selection_item_button_clicked() {
    Scene_polyhedron_item* poly_item = getSelectedItem<Scene_polyhedron_item>();
    if(!poly_item) {
      print_message("Error: there is no selected polyhedron item!");
      return; 
    }
    // all other arrangements (putting inside selection_item_map), setting names etc,
    // other params (e.g. k_ring) will be set inside new_item_created
    Scene_polyhedron_selection_item* new_item = new Scene_polyhedron_selection_item(poly_item, mw);
    ui_widget.modeBox->setCurrentIndex(last_mode);
    connectItem(new_item);
  }
  void on_Selection_type_combo_box_changed(int index) {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_active_handle_type(static_cast<Active_handle::Type>(index));
      Q_EMIT save_handleType();
      if(index == 1)
      {
        ui_widget.Select_all_NTButton->show();
        ui_widget.Add_to_selection_button->hide();
        Q_EMIT set_operation_mode(-1);
      }
      else if(index == 4)
      {
        it->second->setPathSelection(true);
        ui_widget.Add_to_selection_button->show();
        ui_widget.Select_all_NTButton->hide();
        Q_EMIT set_operation_mode(-2);
      }
      else
      {
        ui_widget.Add_to_selection_button->hide();
        ui_widget.Select_all_NTButton->hide();
        it->second->setPathSelection(false);
        Q_EMIT set_operation_mode(-1);
      }
    }
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
    switch(ui_widget.operationsBox->currentIndex())
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
      for(Scene_polyhedron_selection_item::Selection_set_vertex::iterator begin = selection_item->selected_vertices.begin();
         begin != selection_item->selected_vertices.end(); ++begin) {
         point_item->point_set()->push_back((*begin)->point());
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
      typedef Polyhedron::Vertex_handle Vertex_handle;
      Edge_graph edge_graph;
      std::map<Vertex_handle, Edge_graph::vertex_descriptor> p2vd;
      std::map<Vertex_handle, Edge_graph::vertex_descriptor>::iterator it_find;
      bool insert_OK;

      for(Scene_polyhedron_selection_item::Selection_set_edge::iterator begin = selection_item->selected_edges.begin();
          begin != selection_item->selected_edges.end(); ++begin)
      {
        Vertex_handle source = begin->halfedge()->opposite()->vertex();
        boost::tie(it_find, insert_OK)
            = p2vd.insert(std::make_pair(source, Edge_graph::vertex_descriptor()));
        if (insert_OK)
        {
          it_find->second = add_vertex(edge_graph);
          edge_graph[it_find->second] = source->point();
        }
        Edge_graph::vertex_descriptor src=it_find->second;

        Vertex_handle target = begin->halfedge()->vertex();
        boost::tie(it_find, insert_OK)
            = p2vd.insert(std::make_pair(target, Edge_graph::vertex_descriptor()));
        if (insert_OK)
        {
          it_find->second = add_vertex(edge_graph);
          edge_graph[it_find->second] = target->point();
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

      Scene_polyhedron_item* poly_item = new Scene_polyhedron_item();
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
    default :
      break;
    }
    return;
  }

  void on_ModeBox_changed(int index)
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(selection_item)
      selection_item->on_Ctrlz_pressed();
    last_mode = index;
    switch(index)
    {
    //Selection mode
    case 0:
      ui_widget.selection_groupBox->setVisible(true);
      ui_widget.edition_groupBox->setVisible(false);
      Q_EMIT set_operation_mode(-1);
      on_Selection_type_combo_box_changed(ui_widget.Selection_type_combo_box->currentIndex());
      break;
      //Edition mode
    case 1:
      ui_widget.selection_groupBox->setVisible(false);
      ui_widget.edition_groupBox->setVisible(true);
      Q_EMIT save_handleType();
      on_editionBox_changed(ui_widget.editionBox->currentIndex());
      break;
    }
  }

  void on_editionBox_changed(int mode )
  {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(selection_item)
      selection_item->on_Ctrlz_pressed();
    if(ui_widget.modeBox->currentIndex() == 0)
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
  }

  void on_Expand_reduce_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = getSelectedItem<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }

    int steps = ui_widget.Expand_reduce_spin_box->value();
    selection_item->expand_or_reduce(steps);
  }
  // To handle empty selection items coming from loader
  void new_item_created(int item_id) {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    Scene_polyhedron_selection_item* selection_item = 
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(item_id));
    if(!selection_item) { return; }

    Scene_polyhedron_item* poly_item = getSelectedItem<Scene_polyhedron_item>();
    if(!poly_item) {
      CGAL_assertion(selection_item->polyhedron_item() == NULL); // which means it is coming from selection_io loader
      print_message("Error: please select corresponding polyhedron item from Geometric Objects list.");
      scene->erase(item_id);
      return;
    }

    if(selection_item->polyhedron_item() == NULL) { //coming from selection_io loader
      if(!selection_item->actual_load(poly_item, mw)) {
        print_message("Error: loading selection item is not successful!");
        scene->erase(item_id);
        return;
      }
      selection_item->invalidateOpenGLBuffers();
      scene->itemChanged(selection_item);
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
    connect(this, SIGNAL(set_operation_mode(int)),selection_item, SLOT(set_operation_mode(int)));
    QObject* scene_ptr = dynamic_cast<QObject*>(scene);
    if (scene_ptr)
      connect(selection_item,SIGNAL(simplicesSelected(CGAL::Three::Scene_item*)), scene_ptr, SLOT(setSelectedItem(CGAL::Three::Scene_item*)));
    connect(selection_item,SIGNAL(isCurrentlySelected(Scene_polyhedron_item_k_ring_selection*)), this, SLOT(isCurrentlySelected(Scene_polyhedron_item_k_ring_selection*)));
  }
  void item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item) {
    // if polyhedron item
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene_item);
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
      Scene_polyhedron_item* poly_item = selection_item->polyhedron_item();
      std::pair<Selection_item_map::iterator, Selection_item_map::iterator> res =
        selection_item_map.equal_range(poly_item);
      for(Selection_item_map::iterator begin = res.first; begin != res.second; ++begin) {
        if(begin->second == selection_item) {
          selection_item_map.erase(begin); break;
        }
      }
    }
  }

private:
  Messages_interface* messages;
  QAction* actionSelection;

  QDockWidget* dock_widget;
  Ui::Selection ui_widget;
typedef std::multimap<Scene_polyhedron_item*, Scene_polyhedron_selection_item*> Selection_item_map;
  Selection_item_map selection_item_map;
  int last_mode;
}; // end Polyhedron_demo_selection_plugin

//Q_EXPORT_PLUGIN2(Polyhedron_demo_selection_plugin, Polyhedron_demo_selection_plugin)

#include "Selection_plugin.moc"
