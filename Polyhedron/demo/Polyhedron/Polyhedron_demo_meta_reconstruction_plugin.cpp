#include <QTime>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QObject>

#include <fstream>

#include <CGAL/array.h>
#include <CGAL/centroid.h>
#include <CGAL/PCA_util.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>

#include "Scene_polygon_soup_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "ui_Polyhedron_demo_meta_reconstruction_plugin.h"

namespace MetaReconstruction
{
  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  // types for K nearest neighbors search
  typedef CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  template <typename OutputIterator>
  void generate_scales (const unsigned int scale_min, double factor,
			OutputIterator out)
  {
    for (unsigned int scale = scale_min; scale < 400;
	 scale = static_cast<unsigned int>(scale * factor))
      *(out ++) = scale;
  }
  
  unsigned int scale_of_anisotropy (const Point_set& points)
  {
    // Tree tree(points.begin(), points.end());
    
    // double ratio_kept = (points.size() < 1000)
    //   ? 1. : 1000. / (points.size());
    
    // std::vector<Point> subset;
    // for (std::size_t i = 0; i < points.size (); ++ i)
    //   if (rand() / (double)RAND_MAX < ratio_kept)
    // 	subset.push_back (points[i]);
    
    // std::vector<unsigned int> scales;
    // generate_scales (6, 1.5, std::back_inserter (scales));

    // std::vector<double> scores (scales.size(), 0.);

    // Point ideal (0.0, 0.5, 0.5);
    
    // for (std::size_t i = 0; i < subset.size (); ++ i)
    //   {
    // 	Neighbor_search search(tree, subset[i],scales.back());
    // 	std::vector<Point> neighbors;

    // 	unsigned int nb = 0;
    // 	std::size_t index = 0;
    // 	for (Search_iterator search_iterator = search.begin();
    // 	     search_iterator != search.end (); ++ search_iterator, ++ nb)
    // 	  {
    // 	    neighbors.push_back (search_iterator->first);

    // 	    if (nb + 1 == scales[index])
    // 	      {
    // 		Point centroid = CGAL::centroid (neighbors.begin (), neighbors.end ());
    // 		CGAL::cpp11::array<double, 6> covariance  = {{ 0., 0., 0., 0., 0., 0. }};
    // 		CGAL::internal::assemble_covariance_matrix_3 (neighbors.begin (),
    // 							      neighbors.end (),
    // 							      covariance,
    // 							      centroid,
    // 							      Kernel(),
    // 							      NULL,
    // 							      CGAL::Dimension_tag<0>());
    // 		CGAL::cpp11::array<double, 3> eigenvalues = {{ 0., 0., 0. }};
    // 		CGAL::Default_diagonalize_traits<double, 3>::diagonalize_selfadjoint_covariance_matrix
    // 		  (covariance, eigenvalues);

    // 		Vector v (eigenvalues[0], eigenvalues[1], eigenvalues[2]);
    // 		v = v / (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
    // 		scores[index] += CGAL::squared_distance (ideal, CGAL::ORIGIN + v);
								
		
    // 		++ index;
    // 		if (index == scales.size ())
    // 		  break;
    // 	      }
    // 	  }
    //   }

    // std::ofstream f ("test.plot");
    // for (std::size_t i = 0; i < scores.size(); ++ i)
    //   f << scales[i] << " " << scores[i] / subset.size() << std::endl;
    // f.close();
    
    return 6;
  }

  
  unsigned int scale_of_noise (const Point_set& points, unsigned int scale_min = 6)
  {
    Tree tree(points.begin(), points.end());
    
    double ratio_kept = (points.size() < 1000)
      ? 1. : 1000. / (points.size());
    
    std::vector<Point> subset;
    for (std::size_t i = 0; i < points.size (); ++ i)
      if (rand() / (double)RAND_MAX < ratio_kept)
    	subset.push_back (points[i]);
    
    std::vector<unsigned int> scales;
    generate_scales (6, 1.5, std::back_inserter (scales));

    std::vector<unsigned int> chosen;

    for (std::size_t i = 0; i < subset.size (); ++ i)
      {
    	Neighbor_search search(tree, subset[i],scales.back());
	double current = 0.;
    	unsigned int nb = 0;
    	std::size_t index = 0;
	double minimum = (std::numeric_limits<double>::max)();
	unsigned int c = 0;
	
    	for (Search_iterator search_iterator = search.begin();
    	     search_iterator != search.end (); ++ search_iterator, ++ nb)
    	  {
	    current += search_iterator->second;

    	    if (nb + 1 == scales[index])
    	      {
		double score = std::sqrt (current / scales[index])
		  / std::pow (scales[index], 0.375); // NB ^ (5/12)

		if (score < minimum)
		  {
		    minimum = score;
		    c = scales[index];
		  }

    		++ index;
    		if (index == scales.size ())
    		  break;
    	      }
    	  }
	chosen.push_back (c);
      }

    std::sort (chosen.begin (), chosen.end());
    std::cerr << "Scale = " << chosen[chosen.size () / 2] << std::endl;
    
    return chosen[chosen.size() / 2];
  }

  void simplify_point_set (const Point_set& points, unsigned int scale)
  {

  }

  void smooth_point_set (const Point_set& points, unsigned int scale)
  {

  }
}

class Polyhedron_demo_meta_reconstruction_plugin_dialog : public QDialog, private Ui::MetaReconstructionOptionsDialog
{
  Q_OBJECT
public:
  Polyhedron_demo_meta_reconstruction_plugin_dialog(QWidget* /*parent*/ = 0)
  {
    setupUi(this);
  }

  bool boundaries() const { return m_boundaries->isChecked(); }
  bool interpolate() const { return m_interpolate->isChecked(); }
};

#include <CGAL/Scale_space_surface_reconstruction_3.h>

class Polyhedron_demo_meta_reconstruction_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionScaleSpaceReconstruction;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionScaleSpaceReconstruction = new QAction(tr("Meta surface reconstruction"), mainWindow);
    actionScaleSpaceReconstruction->setObjectName("actionMetaReconstruction");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionScaleSpaceReconstruction;
  }

public Q_SLOTS:
  void on_actionMetaReconstruction_triggered();
}; // end class Polyhedron_meta_reconstruction_plugin


void Polyhedron_demo_meta_reconstruction_plugin::on_actionMetaReconstruction_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

      //generate the dialog box to set the options
      Polyhedron_demo_meta_reconstruction_plugin_dialog dialog;
      if(!dialog.exec())
	return;

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      std::cout << "Meta surface reconstruction with the following requirements:" << std::endl
		<< (dialog.boundaries() ? " * Output shape has boundaries" : " * Output shape is closed") << std::endl
		<< (dialog.interpolate() ? " * Output shape passes through input points"
		    : " * Output shape approximates input points") << std::endl;


      std::cout << "Analysis of input point set:" << std::endl;
      time.start();

      // TODO: analyse if point set is isotropic
      unsigned int aniso_scale = MetaReconstruction::scale_of_anisotropy (*points);
      bool isotropic = (aniso_scale > 6);

      // TODO: analyse if point set is noisy
      unsigned int noise_scale = MetaReconstruction::scale_of_noise (*points, aniso_scale);
      bool noisy = (noise_scale > 6);
      
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      
      if (!(dialog.interpolate()) && (noisy || isotropic))
	{
	  points = new Point_set();
	  std::copy (pts_item->point_set()->begin(), pts_item->point_set()->end(),
		     std::back_inserter (*points));
	  
	  std::cout << "Preprocessing:" << std::endl;
	  time.restart();

	  if (isotropic)
	    {
	      // TODO: simplify point set
	      MetaReconstruction::simplify_point_set (*points, aniso_scale);
	    }
	  if (noisy)
	    {
	      // TODO: smooth point set
	      MetaReconstruction::smooth_point_set (*points, noise_scale);
	    }

	  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	}
      

      if (dialog.interpolate())
	{
	  if (noisy)
	    { 
	      std::cout << "Scale space reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      std::cout << "Advancing front reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }

	}
      else
	{
	  if (dialog.boundaries())
	    {
	      std::cout << "Scale space reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      std::cout << "Poisson reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	}

      if (!(dialog.interpolate()) && (noisy || isotropic))
	delete points;

      // default cursor
      QApplication::restoreOverrideCursor();
    }
}

#include "Polyhedron_demo_meta_reconstruction_plugin.moc"
