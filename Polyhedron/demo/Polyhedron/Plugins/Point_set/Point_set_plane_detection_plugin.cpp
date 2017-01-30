#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

#include <CGAL/Random.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>

#include "ui_Point_set_plane_detection_plugin.h"


using namespace CGAL::Three;
class Polyhedron_demo_point_set_plane_detection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  QAction* actionDetect;
  QAction* actionEstimateParameters;

  typedef Point_set_3<Kernel>::Point_map PointPMap;
  typedef Point_set_3<Kernel>::Vector_map NormalPMap;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    actionDetect = new QAction(tr("Point Set Plane Detection"), mainWindow);
    actionDetect->setObjectName("actionDetect");
    actionEstimateParameters = new QAction(tr("Point Set Plane Detection (parameter estimation)"), mainWindow);
    actionEstimateParameters->setObjectName("actionEstimateParameters");
    autoConnectActions();
  }

  bool applicable(QAction*) const {
    Scene_points_with_normal_item* item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if (item && item->has_normals())
      return true;
    return false;
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionDetect << actionEstimateParameters;
  }

public Q_SLOTS:
  void on_actionDetect_triggered();
  void on_actionEstimateParameters_triggered();

private:

  typedef Kernel::Plane_3 Plane_3;
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;

  template <typename PointSet, typename OutputIterator>
  void region_growing(PointSet& point_set,
                      double cluster_epsilon,
                      double epsilon,
                      std::size_t Nmin,
                      double maximal_deviation_of_normals,
                      OutputIterator output)
  {
    typedef typename PointSet::Index Index;
    
    typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
    typedef CGAL::Search_traits_adapter <Index,
                                         typename PointSet::Point_map, SearchTraits_3> Search_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Neighbor_search::Distance Distance;

    // build kdtree
    Tree tree(point_set.begin(),
              point_set.end(),
              typename Tree::Splitter(),
              Search_traits (point_set.point_map())
              );
    Distance tr_dist(point_set.point_map());
    
    //Initialization structures
    std::vector <int> label_subregion (point_set.size(), -1);
    int class_index = -1;

    for (typename PointSet::iterator it = point_set.begin(); it != point_set.end(); ++ it)
      {
        if (label_subregion[*it] != -1)
          continue;

        label_subregion[*it] = class_index++;

        int conti = 0; 	//for accelerate least_square fitting 

        Plane_3 optimal_plane(point_set.point(*it),
                              point_set.normal(*it));

        //initialization containers
        std::vector<Index> index_container (1, *it);
        std::vector<Index> index_container_former_ring (1, *it);
        std::list<Index> index_container_current_ring;

        //propagation
        bool propagation = true;
        do
          {
            propagation = false;

            for (std::size_t k = 0; k < index_container_former_ring.size(); k++)
              {
                Index point_index = index_container_former_ring[k];
                
                Neighbor_search search(tree, point_set.point(point_index), 10, 0, true, tr_dist);
                
                for (typename Neighbor_search::iterator nit = search.begin();
                     nit != search.end(); ++ nit)
                  {
                    if (nit->second > cluster_epsilon * cluster_epsilon)
                      break;
                    
                    Index neighbor_index = nit->first;

                    if (label_subregion[neighbor_index] != -1)
                      continue;

                    const Point_3& neighbor = point_set.point(neighbor_index);
                    Point_3 neighbor_projection = optimal_plane.projection(neighbor);
                    double distance = CGAL::squared_distance(neighbor, neighbor_projection);

                    if (distance > epsilon * epsilon
                        || std::fabs(point_set.normal(neighbor_index)
                                     * optimal_plane.orthogonal_vector()) < maximal_deviation_of_normals)
                      continue;

                    label_subregion[neighbor_index] = class_index;
                    propagation = true;
                    index_container_current_ring.push_back(neighbor_index);
                    conti++;

                    if ((conti<50 && conti % 10 == 0) || (conti>50 && conti % 500 == 0))
                      {
                        std::list<Point_3> listp;
                        for (std::size_t pm = 0; pm < index_container.size(); pm++)
                          listp.push_back(point_set.point(index_container[pm]));

                        Plane_3 reajusted_plane;
                        CGAL::linear_least_squares_fitting_3(listp.begin(),
                                                             listp.end(),
                                                             reajusted_plane,
                                                             CGAL::Dimension_tag<0>());
                        optimal_plane = reajusted_plane;
                      }
                  }
              }

            //update containers
            index_container_former_ring.clear();
            for (typename std::list<Index>::iterator lit = index_container_current_ring.begin();
                 lit != index_container_current_ring.end(); ++lit)
              {
                index_container_former_ring.push_back(*lit);
                index_container.push_back(*lit);
              }
            index_container_current_ring.clear();
          }
        while (propagation);

        if (index_container.size() >= Nmin)
          {
            std::vector<Point_3> out;
            for (std::size_t k = 0; k < index_container.size(); k++)
              out.push_back (point_set.point(index_container[k]));
            *(output ++) = out;
          }
        else
          {
            class_index--;
            label_subregion[*it] = -1;
            for (std::size_t k = 0; k < index_container.size(); k++)
              label_subregion[index_container[k]] = -1;
          }
      }
  }


  Kernel::Point_2 to_2d (const Point_3& centroid,
                         const Vector_3& base1,
                         const Vector_3& base2,
                         const Point_3& query)
  {
    Vector_3 v (centroid, query);
    return Kernel::Point_2 (v * base1, v * base2);
  }

  Point_3 to_3d (const Point_3& centroid,
                 const Vector_3& base1,
                 const Vector_3& base2,
                 const Kernel::Point_2& query)
  {
    return centroid + query.x() * base1 + query.y() * base2;
  }

  void build_alpha_shape (Point_set& points, const Point_3& centroid, const Plane_3& plane,
                          Scene_polyhedron_item* item, double epsilon);

}; // end Polyhedron_demo_point_set_plane_detection_plugin

class Point_set_demo_point_set_plane_detection_dialog : public QDialog, private Ui::PointSetPlaneDetectionDialog
{
  Q_OBJECT
public:
  Point_set_demo_point_set_plane_detection_dialog(QWidget * /*parent*/ = 0)
  {
    setupUi(this);
  }

  double cluster_epsilon() const { return m_cluster_epsilon_field->value(); }
  double epsilon() const { return m_epsilon_field->value(); }
  unsigned int min_points() const { return m_min_pts_field->value(); }
  double normal_tolerance() const { return m_normal_tolerance_field->value(); }
  bool generate_alpha() const { return m_generate_alpha->isChecked(); }
  bool generate_subset() const { return !(m_do_not_generate_subset->isChecked()); }
};

void Polyhedron_demo_point_set_plane_detection_plugin::on_actionDetect_triggered() {

  CGAL::Random rand(time(0));
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
    {
      // Gets point set
      Point_set* points = item->point_set();

      if(points == NULL)
        return;

      // Gets options
      Point_set_demo_point_set_plane_detection_dialog dialog;
      if(!dialog.exec())
        return;
    
      Scene_group_item *planes_item = new Scene_group_item(QString("%1 (point subsets)").arg(item->name()));
      if (dialog.generate_subset())
        scene->addItem(planes_item);
      planes_item->setExpanded(false);
      Scene_group_item *alpha_item = new Scene_group_item(QString("%1 (alpha shapes)").arg(item->name()));
      if (dialog.generate_alpha())
        scene->addItem(alpha_item);
      alpha_item->setExpanded(false);
    
      QApplication::setOverrideCursor(Qt::WaitCursor);

      std::vector<std::vector<Point_3> > planes;
      
      region_growing(*points,
                     dialog.cluster_epsilon(),
                     dialog.epsilon(),
                     dialog.min_points(),
                     dialog.normal_tolerance(),
                     std::back_inserter (planes));

      std::cerr << planes.size() << " shape(s) found" << std::endl;
      
      for (std::size_t i = 0; i < planes.size(); ++ i)
        {
          Scene_points_with_normal_item *point_item = new Scene_points_with_normal_item;
          point_item->point_set()->add_normal_map();

          for (std::size_t n = 0; n < planes[i].size(); ++ n)
            point_item->point_set()->insert(planes[i][n]);
      
          Point_3 centroid;
          Plane_3 plane;
          CGAL::linear_least_squares_fitting_3(planes[i].begin(),
                                               planes[i].end(),
                                               plane, centroid,
                                               CGAL::Dimension_tag<0>());
          unsigned char r, g, b;
          r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          b = static_cast<unsigned char>(64 + rand.get_int(0, 192));

          point_item->setRbgColor(r, g, b);

          point_item->setName (QString("Plane #%1").arg(item->name()));
          if (dialog.generate_alpha ())
            {
              Scene_polyhedron_item* poly_item = new Scene_polyhedron_item;

              build_alpha_shape (*(point_item->point_set()), centroid, plane,
                                 poly_item, dialog.cluster_epsilon());
          
              poly_item->setColor(point_item->color ());
              poly_item->setName (QString("Alpha shape #%1").arg(item->name()));
              poly_item->setRenderingMode (Flat);

              scene->addItem(poly_item);
              scene->changeGroup(poly_item, alpha_item);
            }

          for(Point_set::iterator it = point_item->point_set()->begin(); it != point_item->point_set()->end(); ++it)
            point_item->point_set()->normal(*it) = plane.orthogonal_vector();

          if (dialog.generate_subset())
            {
              scene->addItem (point_item);
              scene->changeGroup(point_item, planes_item);
            }
          else
            delete point_item;
        }
      if (!(dialog.generate_subset()))
        delete planes_item;
    
      if (!(dialog.generate_alpha()))
        delete alpha_item;

      scene->itemChanged(index);

      QApplication::restoreOverrideCursor();

      item->setVisible(false);
    }
}

void Polyhedron_demo_point_set_plane_detection_plugin::on_actionEstimateParameters_triggered() {

  CGAL::Random rand(time(0));
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
    {
      // Gets point set
      Point_set* points = item->point_set();

      if(points == NULL)
        return;

      if (points->nb_selected_points() == 0)
        {
          QMessageBox::information(NULL,
                                   tr("Warning"),
                                   tr("Selection is empty.\nTo estimate parameters, please select a planar section."));
          return;
        }
      
      QApplication::setOverrideCursor(Qt::WaitCursor);

      typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
      typedef CGAL::Search_traits_adapter <Point_set::Index,
                                           Point_set::Point_map, SearchTraits_3> Search_traits;
      typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
      typedef typename Neighbor_search::Tree Tree;
      typedef typename Neighbor_search::Distance Distance;

      // build kdtree
      Tree tree(points->first_selected(),
                points->end(),
                typename Tree::Splitter(),
                Search_traits (points->point_map())
                );
      Distance tr_dist(points->point_map());

      Plane_3 plane;
      CGAL::linear_least_squares_fitting_3(boost::make_transform_iterator
                                           (points->first_selected(),
                                            CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                            (points->point_map())),
                                           boost::make_transform_iterator
                                           (points->end(),
                                            CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                            (points->point_map())),
                                           plane,
                                           CGAL::Dimension_tag<0>());

      std::vector<double> epsilon, dispersion, cluster_epsilon;

      Vector_3 norm = plane.orthogonal_vector();
      norm = norm / std::sqrt (norm * norm);
      for (Point_set::iterator it = points->first_selected(); it != points->end(); ++ it)
        {
          double dist = CGAL::squared_distance (plane, points->point(*it));
          epsilon.push_back(dist);

          double disp = std::fabs (norm * points->normal(*it));
          dispersion.push_back (disp);

          Neighbor_search search(tree, points->point(*it), 2, 0, true, tr_dist);
          typename Neighbor_search::iterator nit = search.begin();
          ++ nit;
          double eps = nit->second;
          cluster_epsilon.push_back(eps);
        }

      std::sort (epsilon.begin(), epsilon.end());
      std::sort (dispersion.begin(), dispersion.end());
      std::sort (cluster_epsilon.begin(), cluster_epsilon.end());
      
      QApplication::restoreOverrideCursor();

      
      QMessageBox::information(NULL,
                               tr("Estimated Parameters"),
                               tr("Epsilon = [%1 ; %2 ; %3 ; %4 ; %5]\nNormal Tolerance = [%6 ; %7 ; %8 ; %9 ; %10]\nMinimum Number of Points = %11\nConnectivity Epsilon = [%12 ; %13 ; %14 ; %15 ; %16]")
                               .arg(epsilon.front())
                               .arg(epsilon[epsilon.size() / 10])
                               .arg(epsilon[epsilon.size() / 2])
                               .arg(epsilon[9 * epsilon.size() / 10])
                               .arg(epsilon.back())
                               .arg(dispersion.back())
                               .arg(dispersion[9 * dispersion.size() / 10])
                               .arg(dispersion[dispersion.size() / 2])
                               .arg(dispersion[dispersion.size() / 10])
                               .arg(dispersion.front())
                               .arg(points->nb_selected_points())
                               .arg(cluster_epsilon.front())
                               .arg(cluster_epsilon[cluster_epsilon.size() / 10])
                               .arg(cluster_epsilon[cluster_epsilon.size() / 2])
                               .arg(cluster_epsilon[9 * cluster_epsilon.size() / 10])
                               .arg(cluster_epsilon.back()));
    }
}

void Polyhedron_demo_point_set_plane_detection_plugin::build_alpha_shape
(Point_set& points,  const Point_3& centroid, const Plane_3& plane,
 Scene_polyhedron_item* item, double epsilon)
{
  typedef Kernel::Point_2  Point_2;
  typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
  typedef CGAL::Alpha_shape_face_base_2<Kernel>  Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> Triangulation_2;
  typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;

  Vector_3 base1 = plane.base1();
  Vector_3 base2 = plane.base2();
  base1 = base1 / std::sqrt(base1 * base1);
  base2 = base2 / std::sqrt(base2 * base2);
  
  std::vector<Point_2> projections;
  projections.reserve (points.size ());

  for (Point_set::const_iterator it = points.begin(); it != points.end(); ++ it)
    projections.push_back (to_2d(centroid, base1, base2, points.point(*it)));

  Alpha_shape_2 ashape (projections.begin (), projections.end (), epsilon);
  
  std::map<Alpha_shape_2::Vertex_handle, std::size_t> map_v2i;

  Scene_polygon_soup_item *soup_item = new Scene_polygon_soup_item;
  
  soup_item->init_polygon_soup(points.size(), ashape.number_of_faces ());
  std::size_t current_index = 0;

  for (Alpha_shape_2::Finite_faces_iterator it = ashape.finite_faces_begin ();
       it != ashape.finite_faces_end (); ++ it)
    {
      if (ashape.classify (it) != Alpha_shape_2::INTERIOR)
        continue;

      for (int i = 0; i < 3; ++ i)
        {
          if (map_v2i.find (it->vertex (i)) == map_v2i.end ())
            {
              map_v2i.insert (std::make_pair (it->vertex (i), current_index ++));
              Point_3 p = to_3d (centroid, base1, base2, it->vertex (i)->point ());
              soup_item->new_vertex (p.x (), p.y (), p.z ());
            }
        }
      soup_item->new_triangle (map_v2i[it->vertex (0)],
                               map_v2i[it->vertex (1)],
                               map_v2i[it->vertex (2)]);
    }

  soup_item->orient();
  soup_item->exportAsPolyhedron (item->polyhedron());

  if (soup_item->isEmpty ())
    {
      std::cerr << "POLYGON SOUP EMPTY" << std::endl;
      for (std::size_t i = 0; i < projections.size (); ++ i)
        std::cerr << projections[i] << std::endl;
      
    }
  
  delete soup_item;
}




#include <QtPlugin>

#include "Point_set_plane_detection_plugin.moc"
