#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"
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


  template <typename OutputIterator>
  void detect_planes_in_mesh (SMesh& mesh, const double area_min,
                              const double angle_max,
                              OutputIterator output)
  {
    typedef SMesh::size_type size_type;
    
    std::cerr << "Detecting planes with:" << std::endl
              << " * Area min = " << area_min << std::endl
              << " * Min cos angle = " << angle_max << std::endl;
    std::vector<int> label_region(mesh.number_of_faces(), 0);
    int class_index = 0;

    
    for (typename SMesh::Face_iterator f = mesh.faces_begin(); f != mesh.faces_end(); ++f)
      {
        if (label_region[*f] != 0)
          continue;
        class_index++;
        label_region[*f] = class_index;
        double area = PMP::face_area (*f, mesh, PMP::parameters::geom_traits(EPICK()));
            
        //characteristics of the seed
        EPICK::Vector_3 normal_seed = PMP::compute_face_normal (*f, mesh, PMP::parameters::geom_traits(EPICK()));
        EPICK::Point_3 pt_seed = mesh.point(target(halfedge(*f, mesh), mesh));
        EPICK::Plane_3 optimal_plane(pt_seed, normal_seed);
                   //        Kernel::Plane_3 optimal_plane = f->plane();
              
        //initialization containers
        std::vector<size_type> index_container (1,*f);
        std::vector<size_type> index_container_former_ring (1, *f);
        std::list<size_type> index_container_current_ring;

        //propagation
        bool propagation = true;
        do{

          propagation = false;

          for (size_type k = 0; k < index_container_former_ring.size(); k++)
            {
              typename SMesh::Halfedge_around_face_circulator
                circ( mesh.halfedge(SMesh::Face_index(index_container_former_ring[k])), mesh)
                , start = circ;

              do
                {
                  if (is_border(*circ, mesh))
                    continue;
                  
                  typename SMesh::Face_index
                    neighbor = mesh.face(opposite(*circ, mesh));
                  size_type neighbor_index = neighbor;
                  if (label_region[neighbor_index] == 0)
                    {
                      EPICK::Vector_3 normal
                        = PMP::compute_face_normal (neighbor, mesh, PMP::parameters::geom_traits(EPICK()));

                      if (std::fabs(normal * optimal_plane.orthogonal_vector()) > angle_max)
                        {
                          label_region[neighbor_index] = class_index;
                          propagation = true;
                          index_container_current_ring.push_back(neighbor_index);
                          area += PMP::face_area (neighbor, mesh, PMP::parameters::geom_traits(EPICK()));
                        }
                    }
                }
              while (++ circ != start);
            }
			
          //update containers
          index_container_former_ring.clear();
          for (std::list<size_type>::iterator it = index_container_current_ring.begin();
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
            label_region[*f] = 0;
            for (size_type k = 0; k < index_container.size(); k++)
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
  
  Scene_surface_mesh_item* poly_item = 
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if (poly_item)
    {
      SMesh& pmesh =*poly_item->polyhedron();

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

      //check_and_set_ids (&pmesh);
      std::vector<int> indices;
      detect_planes_in_mesh (pmesh,
                             ui.minimumAreaDoubleSpinBox->value(),
                             std::fabs(std::cos (CGAL_PI * ui.maximumDeviationFromNormalSpinBox->value() / 180.)),
                             std::back_inserter (indices));

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
