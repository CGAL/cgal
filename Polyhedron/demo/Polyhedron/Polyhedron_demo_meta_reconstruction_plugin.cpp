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
#include <CGAL/compute_average_spacing.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>

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
  void generate_scales (OutputIterator out,
			const unsigned int scale_min = 6,
			const double factor = 1.15,
			unsigned int nb_scales = 30)
  {
    unsigned int prev = -1;
    
    for (unsigned int i = 0; i < nb_scales; ++ i)
      {
	unsigned int current = static_cast<unsigned int>(scale_min * std::pow (factor, i));
	if (current != prev)
	  {
	    *(out ++) = current;
	    prev = current;
	  }
	else
	  ++ nb_scales;
      }
  }
  
  unsigned int scale_of_anisotropy (const Point_set& points, double& size)
  {
    Tree tree(points.begin(), points.end());
    
    double ratio_kept = (points.size() < 1000)
      ? 1. : 1000. / (points.size());
    
    std::vector<Point> subset;
    for (std::size_t i = 0; i < points.size (); ++ i)
      if (rand() / (double)RAND_MAX < ratio_kept)
    	subset.push_back (points[i]);
    
    std::vector<unsigned int> scales;
    generate_scales (std::back_inserter (scales));

    std::vector<unsigned int> chosen;

    for (std::size_t i = 0; i < subset.size (); ++ i)
      {
    	Neighbor_search search(tree, subset[i],scales.back());
	double current = 0.;
    	unsigned int nb = 0;
    	std::size_t index = 0;
	double maximum = 0.;
	unsigned int c = 0;
	
    	for (Search_iterator search_iterator = search.begin();
    	     search_iterator != search.end (); ++ search_iterator, ++ nb)
    	  {
	    current += search_iterator->second;

    	    if (nb + 1 == scales[index])
    	      {
		double score = std::sqrt (current / scales[index])
		  / std::pow (scales[index], 0.75); // NB ^ (3/4)

		if (score > maximum)
		  {
		    maximum = score;
		    c = scales[index];
		  }

    		++ index;
    		if (index == scales.size ())
    		  break;
    	      }
    	  }
	chosen.push_back (c);
      }

    double mean = 0.;
    for (std::size_t i = 0; i < chosen.size(); ++ i)
      mean += chosen[i];
    mean /= chosen.size();
    unsigned int aniso_scale = static_cast<unsigned int>(mean);

    size = 0.;
    for (std::size_t i = 0; i < subset.size (); ++ i)
      {
    	Neighbor_search search(tree, subset[i], aniso_scale);
	size += std::sqrt ((-- search.end())->second);
      }
    size /= subset.size();
    
    return mean;
  }

  
  unsigned int scale_of_noise (const Point_set& points)
  {
    Tree tree(points.begin(), points.end());
    
    double ratio_kept = (points.size() < 1000)
      ? 1. : 1000. / (points.size());
    
    std::vector<Point> subset;
    for (std::size_t i = 0; i < points.size (); ++ i)
      if (rand() / (double)RAND_MAX < ratio_kept)
    	subset.push_back (points[i]);
    
    std::vector<unsigned int> scales;
    generate_scales (std::back_inserter (scales));

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
    
    return chosen[chosen.size() / 2];
  }

  void simplify_point_set (Point_set& points, double size)
  {
    points.erase (CGAL::grid_simplify_point_set (points.begin (), points.end (), size),
		  points.end ());
  }

  void smooth_point_set (Point_set& points, unsigned int scale)
  {
    CGAL::jet_smooth_point_set(points.begin(), points.end(), scale);
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

      Scene_points_with_normal_item* new_item = NULL;
      if (!(dialog.interpolate()))
	{
	  new_item = new Scene_points_with_normal_item();
	  new_item->setName(QString("%1 (preprocessed)").arg(pts_item->name()));
	  new_item->set_has_normals (pts_item->has_normals());
	  new_item->setColor(pts_item->color());
	  new_item->setRenderingMode(pts_item->renderingMode());
	  new_item->setVisible(pts_item->visible());
	  new_item->resetSelection();
	  new_item->invalidate_buffers();

	  points = new_item->point_set();
	  std::copy (pts_item->point_set()->begin(), pts_item->point_set()->end(),
		     std::back_inserter (*points));
	}

      std::cerr << "Analysing isotropy of point set... ";
      time.start();
      double aniso_size;
      unsigned int aniso_scale = MetaReconstruction::scale_of_anisotropy (*points, aniso_size);
      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;

      bool isotropic = (aniso_scale == 6);
      std::cerr << (isotropic ? " -> Point set is isotropic" : " -> Point set is anisotropic") << std::endl;
      if (!(dialog.interpolate()) && !isotropic)
	{
	  std::cerr << "Correcting anisotropy of point set... ";
	  time.restart();
	  std::size_t prev_size = points->size ();
	  MetaReconstruction::simplify_point_set (*points, aniso_size);
	  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	  std::cerr << " -> " << prev_size - points->size() << " point(s) removed ("
		    << 100. * (prev_size - points->size()) / (double)(prev_size)
		    << "%)" << std::endl;
	}

      std::cerr << "Analysing noise of point set... ";
      time.restart();
      unsigned int noise_scale = MetaReconstruction::scale_of_noise (*points);
      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
      bool noisy = (noise_scale > 6);
      std::cerr << (noisy ? " -> Point set is noisy" : " -> Point set is noise-free") << std::endl;
      
      if (!(dialog.interpolate()) && noisy)
	{
	  std::cerr << "Denoising point set... ";
	  time.restart();
	  MetaReconstruction::smooth_point_set (*points, noise_scale);
	  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	}

      if (dialog.interpolate())
	{
	  if (noisy)
	    { 
	      std::cerr << "Scale space reconstruction... ";
	      time.restart();
	      // TODO

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      std::cerr << "Advancing front reconstruction... ";
	      time.restart();
	      // TODO

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }

	}
      else
	{
	  if (dialog.boundaries())
	    {
	      std::cerr << "Scale space reconstruction... ";
	      time.restart();
	      // TODO

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      if (!(new_item->has_normals()))
		{
		  std::cerr << "Estimation of normal vectors... ";
		  time.restart();
		  // TODO

		  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
		}
	      
	      std::cerr << "Poisson reconstruction... ";
	      time.restart();
	      // TODO

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	}

      if (!(dialog.interpolate()))
	{
	  if (noisy || !isotropic)
	    scene->addItem(new_item);
	  else
	    delete new_item;
	}

      // default cursor
      QApplication::restoreOverrideCursor();
    }
}

#include "Polyhedron_demo_meta_reconstruction_plugin.moc"
