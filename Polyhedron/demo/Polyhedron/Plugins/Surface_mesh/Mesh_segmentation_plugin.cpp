#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "ui_Mesh_segmentation_widget.h"
#include "Scene_surface_mesh_item.h"
#include "Scene.h"
#include "Color_map.h"

#include <CGAL/mesh_segmentation.h>
#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QElapsedTimer>
#include <QAction>
#include <QDebug>
#include <QObject>
#include <QDockWidget>
#include <QMessageBox>
//#include <QtConcurrentRun>
#include <map>
#include <algorithm>
#include <vector>
#include <CGAL/property_map.h>


template<class FaceGraphWithId, class ValueType>
struct FaceGraph_with_id_to_vector_property_map
    : public boost::put_get_helper<ValueType&,
             FaceGraph_with_id_to_vector_property_map<FaceGraphWithId, ValueType> >
{
public:
    typedef typename boost::graph_traits<FaceGraphWithId>::face_descriptor key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

  typedef typename boost::property_map<FaceGraphWithId, CGAL::face_index_t>::type Fidmap;

    FaceGraph_with_id_to_vector_property_map() : internal_vector(NULL) { }
    FaceGraph_with_id_to_vector_property_map(std::vector<ValueType>* internal_vector,
                                             Fidmap idmap)
         : internal_vector(internal_vector), idmap(idmap)
    {}

    reference operator[](key_type key) const
    {
      return (*internal_vector)[get(idmap, key)];
    }
private:
    std::vector<ValueType>* internal_vector;
    Fidmap idmap;

};

using namespace CGAL::Three;
class Polyhedron_demo_mesh_segmentation_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
    typedef std::map<Scene_surface_mesh_item*, std::vector<double> > SM_item_sdf_map;
public:

    QList<QAction*> actions() const {
        return QList<QAction*>() << actionSegmentation;
    }

    bool applicable(QAction*) const {
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    }

    void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
        this->scene = scene_interface;
        this->mw = mainWindow;
        actionSegmentation = new QAction("Mesh Segmentation", mw);
        actionSegmentation->setProperty("subMenuName", "Triangulated Surface Mesh Segmentation");
        actionSegmentation->setObjectName("actionSegmentation");

        // adding slot for itemAboutToBeDestroyed signal, aim is removing item from item-functor map.

        if( Scene* scene = dynamic_cast<Scene*>(scene_interface) ) {
            connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)));
        }
        init_color_map_sdf();
        init_color_map_segmentation();
        autoConnectActions();

        dock_widget = new QDockWidget("Mesh segmentation parameters", mw);
        dock_widget->setVisible(false); // do not show at the beginning
        ui_widget.setupUi(dock_widget);
        ui_widget.Smoothness_spin_box->setMaximum(10.0);
        mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

        connect(ui_widget.Partition_button,  SIGNAL(clicked()), this, SLOT(on_Partition_button_clicked()));
        connect(ui_widget.SDF_button,  SIGNAL(clicked()), this, SLOT(on_SDF_button_clicked()));
    }
    virtual void closure()
    {
      dock_widget->hide();
    }
    void check_and_set_ids(SMesh* face_graph)
    {
      face_graph->collect_garbage();
    }

    template<class SDFPropertyMap, class SceneFacegraphItem>
    void colorize_sdf(SceneFacegraphItem* item, SDFPropertyMap sdf_values, std::vector<QColor>& color_vector);
    template<class SegmentPropertyMap, class SceneFacegraphItem>
    void colorize_segmentation(SceneFacegraphItem* item, SegmentPropertyMap segment_ids, std::vector<QColor>& color_vector);

    void init_color_map_sdf();
    void init_color_map_segmentation();
    template<class FacegraphItem>
    void apply_SDF_button_clicked(FacegraphItem* item);
    template<class FacegraphItem>
    void apply_Partition_button_clicked(FacegraphItem* item);

    public Q_SLOTS:
        void on_actionSegmentation_triggered();
        void on_Partition_button_clicked();
        void on_SDF_button_clicked();
        void itemAboutToBeDestroyed(CGAL::Three::Scene_item*);
private:
    QAction*                      actionSegmentation;
    QDockWidget*                  dock_widget;
    Ui::Mesh_segmentation         ui_widget;

    std::vector<QColor>  color_map_sdf;
    std::vector<QColor>  color_map_segmentation;
    SM_item_sdf_map      sm_item_sdf_map;

    SM_item_sdf_map& get_sdf_map(Scene_surface_mesh_item*)
    {
      return sm_item_sdf_map;
    }
};

void Polyhedron_demo_mesh_segmentation_plugin::init_color_map_sdf()
{
    color_map_sdf = std::vector<QColor>(256);
    int r = 0, g = 0, b = 255;
    for(int i = 0; i <= 255; ++i)
    {
        if(i > 128 && i <= 192) { r = static_cast<int>( ((i - 128) / (192.0 - 128)) * 255 ); }
        if(i > 0 && i <= 98)    { g = static_cast<int>( ((i) / (98.0)) * 255 ); }
        if(i > 191 && i <=255)  { g = 255 - static_cast<int>( ((i - 191) / (255.0 - 191)) * 255 ); }
        if(i > 64 && i <= 127)  { b = 255 - static_cast<int>( ((i - 64) / (127.0 - 64)) * 255 ); }
        color_map_sdf[i] = QColor(r, g, b);
    }
}

void Polyhedron_demo_mesh_segmentation_plugin::init_color_map_segmentation()
{

    color_map_segmentation.push_back(QColor( 173, 35, 35));
    color_map_segmentation.push_back(QColor( 87, 87, 87));
    color_map_segmentation.push_back(QColor( 42, 75, 215));
    color_map_segmentation.push_back(QColor( 29, 105, 20));
    color_map_segmentation.push_back(QColor( 129, 74, 25));
    color_map_segmentation.push_back(QColor( 129, 38, 192));
    color_map_segmentation.push_back(QColor( 160, 160, 160));
    color_map_segmentation.push_back(QColor( 129, 197, 122));
    color_map_segmentation.push_back(QColor( 157, 175, 255));
    color_map_segmentation.push_back(QColor( 41, 208, 208));
    color_map_segmentation.push_back(QColor( 255, 146, 51));
    color_map_segmentation.push_back(QColor( 255, 238, 51));
    color_map_segmentation.push_back(QColor( 233, 222, 187));
    color_map_segmentation.push_back(QColor( 255, 205, 243));

}

void Polyhedron_demo_mesh_segmentation_plugin::itemAboutToBeDestroyed(CGAL::Three::Scene_item* scene_item)
{
    if(Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene_item)) {
      sm_item_sdf_map.erase(sm_item);
    }
}

void Polyhedron_demo_mesh_segmentation_plugin::on_actionSegmentation_triggered()
{ dock_widget->show(); dock_widget->raise();}

template<class FacegraphItem>
void Polyhedron_demo_mesh_segmentation_plugin::apply_SDF_button_clicked(FacegraphItem* item)
{
  typedef typename FacegraphItem::Face_graph Facegraph;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::size_t number_of_rays = ui_widget.Number_of_rays_spin_box->value();
  double cone_angle = (ui_widget.Cone_angle_spin_box->value()  / 180.0) * CGAL_PI;
  bool create_new_item = ui_widget.New_item_check_box->isChecked();

  typedef std::map<FacegraphItem*, std::vector<double> > Pair_map;
  typedef typename Pair_map::value_type Pair;
  typename Pair_map::iterator pair;

  FacegraphItem* active_item = item;

  if(create_new_item) {
      active_item = new FacegraphItem(*item->face_graph());
      active_item->setFlatPlusEdgesMode();
  }

  check_and_set_ids(active_item->face_graph());
  pair = get_sdf_map(active_item).insert(
          Pair(active_item, std::vector<double>()) ).first;

  pair->second.resize(num_faces(*item->face_graph()), 0.0);
  typename boost::property_map<Facegraph, CGAL::face_index_t>::type fidmap =
      get(CGAL::face_index, *pair->first->face_graph());
  FaceGraph_with_id_to_vector_property_map<Facegraph, double> sdf_pmap(&pair->second, fidmap);
  QElapsedTimer time;
  time.start();
  std::pair<double, double> min_max_sdf = sdf_values(*(pair->first->face_graph()), sdf_pmap, cone_angle, number_of_rays);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  std::cout << "SDF computation is completed. Min-SDF : " << min_max_sdf.first << " " "Max-SDF : " << min_max_sdf.second << std::endl;
  colorize_sdf(pair->first, sdf_pmap, pair->first->color_vector());

  pair->first->setName(tr("(SDF-%1-%2)").arg(number_of_rays).arg(ui_widget.Cone_angle_spin_box->value()));

  if(create_new_item) {
      Scene::Item_id index = scene->addItem(pair->first);
      item->setVisible(false);
      scene->itemChanged(item);
      pair->first->invalidateOpenGLBuffers();
      scene->itemChanged(pair->first);
      scene->setSelectedItem(index);
  }
  else {
    item->invalidateOpenGLBuffers();
    scene->itemChanged(scene->item_id(item));
  }
  QApplication::restoreOverrideCursor();
}
void Polyhedron_demo_mesh_segmentation_plugin::on_SDF_button_clicked()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(sm_item)
    apply_SDF_button_clicked(sm_item);
}

template<class FacegraphItem>
void Polyhedron_demo_mesh_segmentation_plugin::apply_Partition_button_clicked(FacegraphItem* item)
{

  typedef typename FacegraphItem::Face_graph Facegraph;
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::size_t number_of_clusters = ui_widget.Number_of_clusters_spin_box->value();
  double smoothness = ui_widget.Smoothness_spin_box->value();
  std::size_t number_of_rays = ui_widget.Number_of_rays_spin_box->value();
  double cone_angle = (ui_widget.Cone_angle_spin_box->value()  / 180.0) * CGAL_PI;
  bool create_new_item = ui_widget.New_item_check_box->isChecked();
  bool extract_segments = ui_widget.Extract_segments_check_box->isChecked();

  typename std::map<FacegraphItem*, std::vector<double> >::iterator pair;
  if(create_new_item)
  {
      // create new item
      FacegraphItem* new_item = new FacegraphItem(*item->face_graph());
      new_item->setFlatPlusEdgesMode();

      // copy SDF values of existing poly to new poly
      typename std::map<FacegraphItem*, std::vector<double> >::iterator it = get_sdf_map(item).find(item);
      const std::vector<double>& sdf_data = it == get_sdf_map(new_item).end() ?
                                            std::vector<double>() : it->second;
      pair = get_sdf_map(new_item).insert(std::make_pair(new_item, sdf_data) ).first;
  }
  else
  {
       std::pair<typename std::map<FacegraphItem*, std::vector<double> >::iterator, bool> res =
        get_sdf_map(item).insert(std::make_pair(item, std::vector<double>()) );
      pair = res.first;
  }
  bool isClosed = is_closed(*pair->first->face_graph());
  if(!isClosed)
  {
    QApplication::restoreOverrideCursor();
    QMessageBox::warning(mw, "Warning", "This mesh has boundaries, therefore the results may be unreliable or meaningless.");
    QApplication::setOverrideCursor(Qt::WaitCursor);
  }
  check_and_set_ids(pair->first->face_graph());
  QElapsedTimer time;
  time.start();
  typename boost::property_map<Facegraph, CGAL::face_index_t>::type fidmap =
      get(CGAL::face_index, *pair->first->face_graph());
  if(pair->second.empty()) { // SDF values are empty, calculate
    pair->second.resize(num_faces(*pair->first->face_graph()), 0.0);
    FaceGraph_with_id_to_vector_property_map<Facegraph, double> sdf_pmap(&pair->second, fidmap);
    sdf_values(*(pair->first->face_graph()), sdf_pmap, cone_angle, number_of_rays);
  }



  std::vector<std::size_t> internal_segment_map(num_faces(*pair->first->face_graph()));
  FaceGraph_with_id_to_vector_property_map<Facegraph, std::size_t> segment_pmap(&internal_segment_map, fidmap);
  FaceGraph_with_id_to_vector_property_map<Facegraph, double> sdf_pmap(&pair->second, fidmap);

  if(!isClosed)
  {
    bool has_sdf_values = false;
    for(typename boost::graph_traits<Facegraph>::face_descriptor f :
                  faces(*pair->first->face_graph()))
    {
      if(sdf_pmap[f] != -1
         && sdf_pmap[f] != (std::numeric_limits<double>::max)())
      {
        has_sdf_values = true;
        break;
      }
    }
    if(!has_sdf_values)
    {
      QApplication::restoreOverrideCursor();
      QMessageBox::warning(mw, "Error", "No SDF value could be computed, aborting...");
      return;
    }
  }
  std::size_t nb_segments = segmentation_from_sdf_values(*(pair->first->face_graph())
      ,sdf_pmap, segment_pmap, number_of_clusters, smoothness, extract_segments);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  std::cout << "Segmentation is completed. Number of segments : " << nb_segments << std::endl;

  colorize_segmentation(pair->first, segment_pmap, pair->first->color_vector());
  pair->first->setName(tr("(Segmentation-%1-%2)").arg(number_of_clusters).arg(smoothness));

  if(create_new_item) {
      Scene::Item_id index = scene->addItem(pair->first);
      item->setVisible(false);
      scene->itemChanged(item);
      pair->first->invalidateOpenGLBuffers();
      scene->itemChanged(pair->first);
      scene->setSelectedItem(index);
  }
  else {
    item->invalidateOpenGLBuffers();
    scene->itemChanged(scene->item_id(item));
  }

  QApplication::restoreOverrideCursor();
}
void Polyhedron_demo_mesh_segmentation_plugin::on_Partition_button_clicked()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(sm_item)
    apply_Partition_button_clicked(sm_item);
}

template<class SDFPropertyMap, class SceneFacegraphItem>
void Polyhedron_demo_mesh_segmentation_plugin::colorize_sdf(
     SceneFacegraphItem* item,
     SDFPropertyMap sdf_values,
     std::vector<QColor>& color_vector)
{
    typedef typename SceneFacegraphItem::Face_graph Facegraph;
    typedef typename boost::graph_traits<Facegraph>::face_iterator face_iterator;
    typedef typename boost::property_map<Facegraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
    typename SceneFacegraphItem::Face_graph* face_graph = item->face_graph();
    color_vector.clear();
    std::size_t patch_id = 0;
    Patch_id_pmap pidmap = get(CGAL::face_patch_id_t<int>(), *item->face_graph());

    for(face_iterator facet_it = faces(*face_graph).begin();
        facet_it != faces(*face_graph).end(); ++facet_it, ++patch_id)
    {
        double sdf_value = sdf_values[*facet_it];
        int gray_color = static_cast<int>(255 * sdf_value);
        if(gray_color < 0 || gray_color >= 256) {
          color_vector.push_back(QColor::fromRgb(0,0,0));
        }
        else {
          color_vector.push_back(color_map_sdf[gray_color]);
        }
        put(pidmap, *facet_it, static_cast<int>(patch_id));
    }
    item->setItemIsMulticolor(true);
    item->computeItemColorVectorAutomatically(false);
}


template<class SegmentPropertyMap, class SceneFacegraphItem>
void Polyhedron_demo_mesh_segmentation_plugin::colorize_segmentation(
     SceneFacegraphItem* item,
     SegmentPropertyMap segment_ids,
     std::vector<QColor>& color_vector)
{
    typedef typename SceneFacegraphItem::Face_graph Facegraph;
    typedef typename boost::graph_traits<Facegraph>::face_iterator face_iterator;
    typedef typename boost::property_map<Facegraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
    Facegraph* face_graph = item->face_graph();
    color_vector.clear();
    std::size_t max_segment = 0;
    Patch_id_pmap pidmap = get(CGAL::face_patch_id_t<int>(), *item->face_graph());
    for(face_iterator facet_it = faces(*face_graph).begin();
        facet_it != faces(*face_graph).end(); ++facet_it)
    {
        std::size_t segment_id = segment_ids[*facet_it];
        put(pidmap, *facet_it, static_cast<int>(segment_id));
        max_segment = (std::max)(max_segment, segment_id);
    }
    for(std::size_t i = 0; i <= max_segment; ++i)
    {
        QColor aColor = color_map_segmentation[(max_segment - i) % color_map_segmentation.size()];
        color_vector.push_back(aColor);
    }
    item->setItemIsMulticolor(true);
    item->computeItemColorVectorAutomatically(true);
    item->setProperty("NbPatchIds", static_cast<int>(max_segment + 1)); //for join_and_split plugin
}

#include "Mesh_segmentation_plugin.moc"
