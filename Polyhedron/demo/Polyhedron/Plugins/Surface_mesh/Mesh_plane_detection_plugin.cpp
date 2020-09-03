#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

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
//#include <QtConcurrentRun>
#include <map>
#include <algorithm>
#include <vector>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

#include "ui_Mesh_plane_detection_dialog.h"

namespace PMP = CGAL::Polygon_mesh_processing;

using namespace CGAL::Three;
class Polyhedron_demo_mesh_plane_detection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  public:

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionPlaneDetection;
  }

  bool applicable(QAction*) const {
    return
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionPlaneDetection = new QAction("Mesh Plane Detection", mw);
    actionPlaneDetection->setProperty("subMenuName", "Triangulated Surface Mesh Segmentation");
    actionPlaneDetection->setObjectName("actionPlaneDetection");

    // adding slot for itemAboutToBeDestroyed signal, aim is removing item from item-functor map.
    if( Scene* scene = dynamic_cast<Scene*>(scene_interface) ) {
      connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)));
    }

    autoConnectActions();
  }
  virtual void closure()
  {
  }

  template<class SegmentPropertyMap>
  void colorize_segmentation(Scene_surface_mesh_item* item, SegmentPropertyMap& segment_ids, std::vector<QColor>& color_vector);

  void check_and_set_ids(SMesh *sm);

public Q_SLOTS:
  void on_actionPlaneDetection_triggered();
  void itemAboutToBeDestroyed(CGAL::Three::Scene_item*);
private:
  QAction*                      actionPlaneDetection;


  //template <typename OutputIterator>
  void detect_planes_in_mesh (SMesh& polygon_mesh,
                              const double max_distance_to_plane,
                              const double max_accepted_angle,
                              const std::size_t min_region_size,
                              std::vector<std::size_t>& indices)
  {

    typedef EPICK::FT FT;
    using Polygon_mesh = SMesh;
    using Face_range   = typename Polygon_mesh::Face_range;
    using Neighbor_query = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Polygon_mesh>;
    using Region_type    = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<EPICK, Polygon_mesh>;
    using Sorting        = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<EPICK, Polygon_mesh, Neighbor_query>;

    std::cerr << "Detecting planes with:" << std::endl
              << " * Max distance to plane = " << max_distance_to_plane << std::endl
              << " * Max angle = " << max_accepted_angle << std::endl
              << " * Minimum region size = " << min_region_size << std::endl;


    using Region  = std::vector<std::size_t>;
    using Regions = std::vector<Region>;
    using Vertex_to_point_map = typename Region_type::Vertex_to_point_map;
    using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;
    const Face_range face_range = faces(polygon_mesh);

    // Create instances of the classes Neighbor_query and Region_type.
    Neighbor_query neighbor_query(polygon_mesh);
    const Vertex_to_point_map vertex_to_point_map(
          get(CGAL::vertex_point, polygon_mesh));
    Region_type region_type(
          polygon_mesh,
          max_distance_to_plane, max_accepted_angle, min_region_size,
          vertex_to_point_map);
    // Sort face indices.
    Sorting sorting(
          polygon_mesh, neighbor_query,
          vertex_to_point_map);
    sorting.sort();
    // Create an instance of the region growing class.
    Region_growing region_growing(
          face_range, neighbor_query, region_type,
          sorting.seed_map());
    // Run the algorithm.
    Regions regions;
    region_growing.detect(std::back_inserter(regions));
    std::cerr << regions.size()<< " planes detected" << std::endl;
    indices.resize(polygon_mesh.number_of_faces());
    for (std::size_t id = 0; id<regions.size(); ++id) {
     const auto& region = regions[id];
      for (const auto index : region){
        indices[index] = id;
      }
    }

    //std::copy (regions.begin(), regions.end(), output);
  }

};

void Polyhedron_demo_mesh_plane_detection_plugin::itemAboutToBeDestroyed(CGAL::Three::Scene_item*)
{
}

void Polyhedron_demo_mesh_plane_detection_plugin::on_actionPlaneDetection_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* poly_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if (poly_item)
    {
      SMesh& pmesh =*poly_item->polyhedron();

      QDialog dialog(mw);
      Ui::Mesh_plane_detection_dialog ui;
      ui.setupUi(&dialog);
      ui.maxDistanceToPlane->setMinimum(0.00001);
      // check user cancellation
      if(dialog.exec() == QDialog::Rejected)
        return;

      QElapsedTimer time;
      time.start();
      std::cerr << "Detecting planes... ";
      QApplication::setOverrideCursor(Qt::WaitCursor);
      QApplication::processEvents();

      //check_and_set_ids (&pmesh);
      std::vector<std::size_t> indices;
      detect_planes_in_mesh (pmesh,
                             ui.maxDistanceToPlane->value(),
                             std::fabs(ui.maximumDeviationFromNormalSpinBox->value()),
                             ui.minRegionSize->value(),
                             indices);

      //poly_item->set_color_vector_read_only(true);
      colorize_segmentation (poly_item, indices, poly_item->color_vector());

      std::cerr << "ok (" << time.elapsed() << " ms, "
                << pmesh.number_of_halfedges() / 2 << " edges)" << std::endl;

      poly_item->invalidateOpenGLBuffers();
      scene->itemChanged(index);
      QApplication::restoreOverrideCursor();
    }

}

template<class SegmentPropertyMap>
void Polyhedron_demo_mesh_plane_detection_plugin::colorize_segmentation(
                                                                        Scene_surface_mesh_item* item,
                                                                        SegmentPropertyMap& segment_ids,
                                                                        std::vector<QColor>& color_vector)
{
  item->setItemIsMulticolor(true);
  item->computeItemColorVectorAutomatically(true);
  SMesh* sm = item->face_graph();
  color_vector.clear();
  std::size_t max_segment = 0;
  for(SMesh::Face_iterator facet_it = sm->faces_begin();
      facet_it != sm->faces_end(); ++facet_it)
    {
      std::size_t segment_id = segment_ids[static_cast<std::size_t>(*facet_it)];
      sm->property_map<face_descriptor, int>("f:patch_id").first[*facet_it] = static_cast<int>(segment_id);
      max_segment = (std::max)(max_segment, segment_id);
    }
  color_vector.push_back(QColor(0, 0, 0));
  for(std::size_t i = 1; i <= max_segment; ++i)
    color_vector.push_back (QColor (CGAL::get_default_random().get_int (41, 255),
                                    CGAL::get_default_random().get_int (35, 238),
                                    CGAL::get_default_random().get_int (35, 255)));
}

#include "Mesh_plane_detection_plugin.moc"
