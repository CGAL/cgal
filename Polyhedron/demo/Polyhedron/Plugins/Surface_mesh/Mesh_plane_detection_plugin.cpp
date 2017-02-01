#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include "Scene.h"
#include "Color_map.h"

#include <CGAL/mesh_segmentation.h>
#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>
#include <QDebug>
#include <QObject>
#include <QDockWidget>
//#include <QtConcurrentRun>
#include <map>
#include <algorithm>
#include <vector>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

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
      qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
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
  void colorize_segmentation(Scene_polyhedron_item* item, SegmentPropertyMap& segment_ids, std::vector<QColor>& color_vector);
  
  void check_and_set_ids(Polyhedron* polyhedron);

public Q_SLOTS:
  void on_actionPlaneDetection_triggered();
  void itemAboutToBeDestroyed(CGAL::Three::Scene_item*);
private:
  QAction*                      actionPlaneDetection;


  template <typename OutputIterator>
  void detect_planes_in_mesh (Polyhedron& mesh, const double area_min,
                              const double angle_max,
                              OutputIterator output)
  {
    std::cerr << "Detecting planes with:" << std::endl
              << " * Area min = " << area_min << std::endl
              << " * Min cos angle = " << angle_max << std::endl;
    std::vector<std::size_t> label_region(mesh.size_of_facets(), 0);
    int class_index = 0;

    std::vector<typename Polyhedron::Facet_handle> facets (mesh.size_of_facets());
    for (typename Polyhedron::Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
      facets[f->id()] = f;
    
    for (typename Polyhedron::Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
      {
        if (label_region[f->id()] != 0)
          continue;
        class_index++;
        label_region[f->id()] = class_index;
        double area = PMP::face_area (f, mesh, PMP::parameters::geom_traits(Kernel()));
            
        //characteristics of the seed
        Kernel::Vector_3 normal_seed = PMP::compute_face_normal (f, mesh, PMP::parameters::geom_traits(Kernel()));
        Kernel::Point_3 pt_seed = f->halfedge()->vertex()->point();
        Kernel::Plane_3 optimal_plane(pt_seed, normal_seed);
                   //        Kernel::Plane_3 optimal_plane = f->plane();
              
        //initialization containers
        std::vector<std::size_t> index_container (1, f->id());
        std::vector<std::size_t> index_container_former_ring (1, f->id());
        std::list<std::size_t> index_container_current_ring;

        //propagation
        bool propagation = true;
        do{

          propagation = false;

          for (std::size_t k = 0; k < index_container_former_ring.size(); k++)
            {
              typename Polyhedron::Halfedge_around_facet_circulator
                circ = facets[index_container_former_ring[k]]->facet_begin(), start = circ;

              do
                {
                  if (circ->is_border_edge())
                    continue;
                  
                  typename Polyhedron::Facet_handle
                    neighbor = circ->opposite()->facet();
                  std::size_t neighbor_index = neighbor->id();
                  if (label_region[neighbor_index] == 0)
                    {
                      Kernel::Vector_3 normal
                        = PMP::compute_face_normal (neighbor, mesh, PMP::parameters::geom_traits(Kernel()));

                      if (std::fabs(normal * optimal_plane.orthogonal_vector()) > angle_max)
                        {
                          label_region[neighbor_index] = class_index;
                          propagation = true;
                          index_container_current_ring.push_back(neighbor_index);
                          area += PMP::face_area (neighbor, mesh, PMP::parameters::geom_traits(Kernel()));
                        }
                    }
                }
              while (++ circ != start);
            }
			
          //update containers
          index_container_former_ring.clear();
          for (std::list<std::size_t>::iterator it = index_container_current_ring.begin();
               it != index_container_current_ring.end(); ++it)
            {
              index_container_former_ring.push_back(*it);
              index_container.push_back(*it);
            }
          index_container_current_ring.clear();

        } while (propagation);

        //Test the number of inliers -> reject if inferior to Nmin
        if (area < area_min)
          {
            class_index--;
            label_region[f->id()] = 0;
            for (std::size_t k = 0; k < index_container.size(); k++)
              label_region[index_container[k]] = 0;
          }
      }
    std::cerr << class_index << " planes detected" << std::endl;

    std::copy (label_region.begin(), label_region.end(), output);
  }

};

void Polyhedron_demo_mesh_plane_detection_plugin::itemAboutToBeDestroyed(CGAL::Three::Scene_item*)
{
}

void Polyhedron_demo_mesh_plane_detection_plugin::on_actionPlaneDetection_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (poly_item)
    {
      Polyhedron& pmesh =*poly_item->polyhedron();

      QDialog dialog(mw);
      Ui::Mesh_plane_detection_dialog ui;
      ui.setupUi(&dialog);
  
      // check user cancellation
      if(dialog.exec() == QDialog::Rejected)
        return;

      QTime time;
      time.start();
      std::cerr << "Detecting planes... ";
      QApplication::setOverrideCursor(Qt::WaitCursor);
      QApplication::processEvents();

      check_and_set_ids (&pmesh);
      std::vector<int> indices;
      detect_planes_in_mesh (pmesh,
                             ui.minimumAreaDoubleSpinBox->value(),
                             std::fabs(std::cos (CGAL_PI * ui.maximumDeviationFromNormalSpinBox->value() / 180.)),
                             std::back_inserter (indices));

      poly_item->set_color_vector_read_only(true);
      colorize_segmentation (poly_item, indices, poly_item->color_vector());
      
      std::cerr << "ok (" << time.elapsed() << " ms, " 
                << pmesh.size_of_halfedges() / 2 << " edges)" << std::endl;

      poly_item->invalidateOpenGLBuffers();
      scene->itemChanged(index);
      QApplication::restoreOverrideCursor();
    }
  
}


void Polyhedron_demo_mesh_plane_detection_plugin::check_and_set_ids(Polyhedron* polyhedron)
{
  std::size_t facet_id = 0;
  for(Polyhedron::Facet_iterator facet_it = polyhedron->facets_begin();
      facet_it != polyhedron->facets_end(); ++facet_it, ++facet_id)
    facet_it->id() = facet_id;
}

template<class SegmentPropertyMap>
void Polyhedron_demo_mesh_plane_detection_plugin::colorize_segmentation(
                                                                        Scene_polyhedron_item* item,
                                                                        SegmentPropertyMap& segment_ids,
                                                                        std::vector<QColor>& color_vector)
{
  item->setItemIsMulticolor(true);
  Polyhedron* polyhedron = item->polyhedron();
  color_vector.clear();
  std::size_t max_segment = 0;
  for(Polyhedron::Facet_iterator facet_it = polyhedron->facets_begin(); 
      facet_it != polyhedron->facets_end(); ++facet_it)   
    {
      std::size_t segment_id = segment_ids[facet_it->id()];
      facet_it->set_patch_id(static_cast<int>(segment_id));
      max_segment = (std::max)(max_segment, segment_id);      
    }
  color_vector.push_back(QColor(0, 0, 0));
  for(std::size_t i = 1; i <= max_segment; ++i)   
    color_vector.push_back (QColor (CGAL::get_default_random().get_int (41, 255),
                                    CGAL::get_default_random().get_int (35, 238),
                                    CGAL::get_default_random().get_int (35, 255)));
}

#include "Mesh_plane_detection_plugin.moc"
