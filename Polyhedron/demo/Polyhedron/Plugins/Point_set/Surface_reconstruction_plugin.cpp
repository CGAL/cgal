#include <QTime>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QObject>

#include <fstream>

#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

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
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>

#include "ui_Surface_reconstruction_plugin.h"

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
Polyhedron* poisson_reconstruct(const Point_set& points,
                                Kernel::FT sm_angle, // Min triangle angle (degrees). 
                                Kernel::FT sm_radius, // Max triangle size w.r.t. point set average spacing. 
                                Kernel::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                                const QString& solver, // solver name
                                bool use_two_passes,
				bool do_not_fill_holes);


struct Radius {

  double bound;

  Radius(double bound)
    : bound(bound)
  {}

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (const AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return adv.smallest_radius_delaunay_sphere (c, index);
    }

    // If radius > bound, return infinity so that facet is not used
    double d  = 0;
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                              c->vertex((index+2)%4)->point()));
    if(d>bound) return adv.infinity();
    d = sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();

    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere (c, index);
  }

};


namespace SurfaceReconstruction
{
  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  // types for K nearest neighbors search
  typedef CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  typedef Neighbor_search::iterator Search_iterator;

  typedef CGAL::Scale_space_surface_reconstruction_3<Kernel> ScaleSpace;
  
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
    
    return aniso_scale;
  }

  
  unsigned int scale_of_noise (const Point_set& points, double& size)
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

    unsigned int noise_scale = chosen[chosen.size() / 2];
    size = 0.;
    for (std::size_t i = 0; i < subset.size (); ++ i)
      {
    	Neighbor_search search(tree, subset[i], noise_scale);
	size += std::sqrt ((-- search.end())->second);
      }
    size /= subset.size();

    
    return noise_scale;
  }

  void simplify_point_set (Point_set& points, double size)
  {
    points.erase (CGAL::grid_simplify_point_set (points.begin (), points.end (), size),
		  points.end ());
  }

  void smooth_point_set (Point_set& points, unsigned int scale)
  {
    CGAL::jet_smooth_point_set<Concurrency_tag>(points.begin(), points.end(),
                                                scale);
  }

  template <typename ItemsInserter>
  void scale_space (const Point_set& points, ItemsInserter items,
		    unsigned int scale, bool generate_smooth = false,
                    bool separate_shells = false, bool force_manifold = true,
                    unsigned int samples = 300, unsigned int iterations = 4)
  {
    ScaleSpace reconstruct (scale, samples);
    reconstruct.reconstruct_surface(points.begin (), points.end (), iterations,
                                    separate_shells, force_manifold);

    for( unsigned int sh = 0; sh < reconstruct.number_of_shells(); ++sh )
      {
        Scene_polygon_soup_item* new_item
          = new Scene_polygon_soup_item ();
        new_item->setColor(Qt::magenta);
        new_item->setRenderingMode(FlatPlusEdges);
        new_item->init_polygon_soup(points.size(), reconstruct.number_of_triangles ());

        Scene_polygon_soup_item* smooth_item = NULL;
        if (generate_smooth)
          {
            smooth_item = new Scene_polygon_soup_item ();
            smooth_item->setColor(Qt::magenta);
            smooth_item->setRenderingMode(FlatPlusEdges);
            smooth_item->init_polygon_soup(points.size(), reconstruct.number_of_triangles ());
          }

        std::map<unsigned int, unsigned int> map_i2i;
        unsigned int current_index = 0;
    
        for (ScaleSpace::Triple_iterator it = reconstruct.shell_begin (sh);
             it != reconstruct.shell_end (sh); ++ it)
          {
            for (unsigned int ind = 0; ind < 3; ++ ind)
              {
                if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
                  {
                    map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
                    Point p = points[(*it)[ind]].position();
                    new_item->new_vertex (p.x (), p.y (), p.z ());
                    
                    if (generate_smooth)
                      {
                        p = *(reconstruct.points_begin() + (*it)[ind]);
                        smooth_item->new_vertex (p.x (), p.y (), p.z ());
                      }
                  }
              }
            new_item->new_triangle( map_i2i[(*it)[0]],
                                    map_i2i[(*it)[1]],
                                    map_i2i[(*it)[2]] );
            if (generate_smooth)
              smooth_item->new_triangle( map_i2i[(*it)[0]],
                                         map_i2i[(*it)[1]],
                                         map_i2i[(*it)[2]] );
              
          }

        *(items ++) = new_item;
        if (generate_smooth)
          *(items ++) = smooth_item;
      }

    if (force_manifold)
      {
        std::ptrdiff_t num = std::distance( reconstruct.garbage_begin(  ),
                                            reconstruct.garbage_end(  ) );

        Scene_polygon_soup_item* new_item
          = new Scene_polygon_soup_item ();
        new_item->setColor(Qt::blue);
        new_item->setRenderingMode(FlatPlusEdges);
        new_item->init_polygon_soup(points.size(), num);

        Scene_polygon_soup_item* smooth_item = NULL;
        if (generate_smooth)
          {
            smooth_item = new Scene_polygon_soup_item ();
            smooth_item->setColor(Qt::blue);
            smooth_item->setRenderingMode(FlatPlusEdges);
            smooth_item->init_polygon_soup(points.size(), num);
          }

        std::map<unsigned int, unsigned int> map_i2i;

        unsigned int current_index = 0;
        for (ScaleSpace::Triple_iterator it=reconstruct.garbage_begin(),
               end=reconstruct.garbage_end();it!=end;++it)
          {
            for (unsigned int ind = 0; ind < 3; ++ ind)
              {
                if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
                  {
                    map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
                    Point p = points[(*it)[ind]].position();
                    new_item->new_vertex (p.x (), p.y (), p.z ());
                    
                    if (generate_smooth)
                      {
                        p = *(reconstruct.points_begin() + (*it)[ind]);
                        smooth_item->new_vertex (p.x (), p.y (), p.z ());
                      }
                  }

              }
            new_item->new_triangle( map_i2i[(*it)[0]],
                                    map_i2i[(*it)[1]],
                                    map_i2i[(*it)[2]] );
            if (generate_smooth)
              smooth_item->new_triangle( map_i2i[(*it)[0]],
                                         map_i2i[(*it)[1]],
                                         map_i2i[(*it)[2]] );
          }

        *(items ++) = new_item;
        if (generate_smooth)
          *(items ++) = smooth_item;

      }
  }
  
  void advancing_front (const Point_set& points, Scene_polyhedron_item* new_item, double size,
                        double radius_ratio_bound = 5., double beta = 0.52)
  {
    Polyhedron& P = * const_cast<Polyhedron*>(new_item->polyhedron());
    Radius filter (size);

    CGAL::advancing_front_surface_reconstruction (points.begin (), points.end (), P, filter,
                                                  radius_ratio_bound, beta);
						  
  }

  void compute_normals (Point_set& points, unsigned int neighbors)
  {
    CGAL::jet_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
                                                CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                                                2 * neighbors);
    
    points.erase (CGAL::mst_orient_normals (points.begin(), points.end(),
					    CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
					    2 * neighbors),
		  points.end ());
  }
  
}



class Polyhedron_demo_surface_reconstruction_plugin_dialog : public QDialog, private Ui::SurfaceReconstructionDialog
{
  Q_OBJECT
public:
  Polyhedron_demo_surface_reconstruction_plugin_dialog(QWidget* /*parent*/ = 0)
  {
    setupUi(this);
    
#ifdef CGAL_EIGEN3_ENABLED
    m_inputSolver->addItem("Eigen - built-in CG");
    m_inputSolver->addItem("Eigen - built-in simplicial LDLt");
#endif
  }

  unsigned int method () const
  {
    if (buttonAuto->isChecked ())       return 0;
    if (buttonAdvancing->isChecked ())  return 1;
    if (buttonScaleSpace->isChecked ()) return 2;
    if (buttonPoisson->isChecked ())    return 3;
    return -1;
  }
  bool boundaries () const { return m_boundaries->isChecked (); }
  bool interpolate () const { return m_interpolate->isChecked (); }
  double longest_edge () const { return m_longestEdge->value (); }
  double radius_ratio_bound () const { return m_radiusRatioBound->value (); }
  double beta_angle () const { return m_betaAngle->value (); }
  unsigned int neighbors () const { return m_neighbors->value (); }
  unsigned int samples () const { return m_samples->value (); }
  unsigned int iterations () const { return m_iterations->value (); }
  bool separate_shells () const { return m_genShells->isChecked (); }
  bool force_manifold () const { return m_forceManifold->isChecked (); }
  bool generate_smoothed () const { return m_genSmooth->isChecked (); }
  double angle () const { return m_inputAngle->value (); }
  double radius () const { return m_inputRadius->value (); }
  double distance () const { return m_inputDistance->value (); }
  bool two_passes () const { return m_inputTwoPasses->isChecked (); }
  bool do_not_fill_holes () const { return m_doNotFillHoles->isChecked (); }
  QString solver () const { return m_inputSolver->currentText (); }
};

#include <CGAL/Scale_space_surface_reconstruction_3.h>

class Polyhedron_demo_surface_reconstruction_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  QAction* actionSurfaceReconstruction;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface) {

    actionSurfaceReconstruction = new QAction(tr("Surface Reconstruction"), mainWindow);
    actionSurfaceReconstruction->setObjectName("actionSurfaceReconstruction");

    CGAL::Three::Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  void automatic_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void advancing_front_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void scale_space_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void poisson_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  
  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSurfaceReconstruction;
  }

public Q_SLOTS:
  void on_actionSurfaceReconstruction_triggered();
}; // end class Polyhedron_surface_reconstruction_plugin


void Polyhedron_demo_surface_reconstruction_plugin::on_actionSurfaceReconstruction_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      //generate the dialog box to set the options
      Polyhedron_demo_surface_reconstruction_plugin_dialog dialog;
      if(!dialog.exec())
	return;

      unsigned int method = dialog.method ();
      switch (method)
        {
        case 0:
          automatic_reconstruction (dialog);
          break;
        case 1:
          advancing_front_reconstruction (dialog);
          break;
        case 2:
          scale_space_reconstruction (dialog);
          break;
        case 3:
          poisson_reconstruction (dialog);
          break;
        default:
          std::cerr << "Error: unkown method." << std::endl;
          return;
        }
      

    }
}

void Polyhedron_demo_surface_reconstruction_plugin::automatic_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

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
	  new_item->invalidateOpenGLBuffers();

	  points = new_item->point_set();
	  std::copy (pts_item->point_set()->begin(), pts_item->point_set()->end(),
		     std::back_inserter (*points));
	}

      std::cerr << "Analysing isotropy of point set... ";
      time.start();
      double aniso_size;
      unsigned int aniso_scale = SurfaceReconstruction::scale_of_anisotropy (*points, aniso_size);
      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;

      std::cerr << " -> Scale / size = " << aniso_scale << " / " << aniso_size << std::endl;
      
      bool isotropic = (aniso_scale == 6);
      std::cerr << (isotropic ? " -> Point set is isotropic" : " -> Point set is anisotropic") << std::endl;
      if (!(dialog.interpolate()) && !isotropic)
	{
	  std::cerr << "Correcting anisotropy of point set... ";
	  time.restart();
	  std::size_t prev_size = points->size ();
	  SurfaceReconstruction::simplify_point_set (*points, aniso_size);
	  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	  std::cerr << " -> " << prev_size - points->size() << " point(s) removed ("
		    << 100. * (prev_size - points->size()) / (double)(prev_size)
		    << "%)" << std::endl;
	}

      std::cerr << "Analysing noise of point set... ";
      time.restart();
      double noise_size;
      unsigned int noise_scale = SurfaceReconstruction::scale_of_noise (*points, noise_size);
      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
      
      std::cerr << " -> Scale / size = " << noise_scale << " / " << noise_size << std::endl;
      
      bool noisy = (noise_scale > 6);
      std::cerr << (noisy ? " -> Point set is noisy" : " -> Point set is noise-free") << std::endl;
      
      if (!(dialog.interpolate()) && noisy)
	{
	  std::cerr << "Denoising point set... ";
	  time.restart();
	  SurfaceReconstruction::smooth_point_set (*points, noise_scale);
          new_item->set_has_normals (false);
	  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	}

      if (dialog.interpolate())
	{
	  if (noisy)
	    { 
	      std::cerr << "Scale space reconstruction... ";
	      time.restart();

              std::vector<Scene_polygon_soup_item*> reco_items;

	      SurfaceReconstruction::scale_space (*points, std::back_inserter (reco_items),
                                                  (std::max)(noise_scale, aniso_scale));
	      
              for (std::size_t i = 0; i < reco_items.size (); ++ i)
                {
                  reco_items[i]->setName(tr("%1 (scale space)").arg(scene->item(index)->name()));
                  scene->addItem (reco_items[i]);
                }

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      std::cerr << "Advancing front reconstruction... ";
	      time.restart();

	      Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
	      SurfaceReconstruction::advancing_front (*points, reco_item, 10. * (std::max)(noise_size, aniso_size));
	      
	      reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
	      reco_item->setColor(Qt::magenta);
	      reco_item->setRenderingMode(FlatPlusEdges);
	      scene->addItem(reco_item);

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }

	}
      else
	{
	  if (dialog.boundaries())
	    {
	      std::cerr << "Advancing front reconstruction... ";
	      time.restart();

	      Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
	      SurfaceReconstruction::advancing_front (*points, reco_item, 10. * (std::max)(noise_size, aniso_size));
	      
	      reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
	      reco_item->setColor(Qt::magenta);
	      reco_item->setRenderingMode(FlatPlusEdges);
	      scene->addItem(reco_item);

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      if (!(new_item->has_normals()))
		{
		  std::cerr << "Estimation of normal vectors... ";
		  time.restart();

		  SurfaceReconstruction::compute_normals (*points, noise_scale);
		  
		  new_item->set_has_normals (true);
		  new_item->setRenderingMode(PointsPlusNormals);
		  
		  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
		}
	      
	      std::cerr << "Poisson reconstruction... ";
	      time.restart();

	      Polyhedron* pRemesh = poisson_reconstruct(*points, 20,
							100 * (std::max)(noise_size, aniso_size),
							(std::max)(noise_size, aniso_size),
							QString ("Eigen - built-in CG"), false, false);
	      if(pRemesh)

		{
		  // Add polyhedron to scene
		  Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(pRemesh);
		  reco_item->setName(tr("%1 (poisson)").arg(pts_item->name()));
		  reco_item->setColor(Qt::lightGray);
		  scene->addItem(reco_item);
		}

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


void Polyhedron_demo_surface_reconstruction_plugin::advancing_front_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      std::cerr << "Advancing front reconstruction... ";

      Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
      SurfaceReconstruction::advancing_front (*points, reco_item,
                                              dialog.longest_edge (),
                                              dialog.radius_ratio_bound (),
                                              CGAL_PI * dialog.beta_angle () / 180.);
	      
      reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
      reco_item->setColor(Qt::magenta);
      reco_item->setRenderingMode(FlatPlusEdges);
      scene->addItem(reco_item);

      QApplication::restoreOverrideCursor();
    }
}


void Polyhedron_demo_surface_reconstruction_plugin::scale_space_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      std::cout << "Scale scape surface reconstruction...";
      
      std::vector<Scene_polygon_soup_item*> reco_items;

      SurfaceReconstruction::scale_space (*points, std::back_inserter (reco_items),
                                          dialog.neighbors (),
                                          dialog.generate_smoothed (),
                                          dialog.separate_shells (),
                                          dialog.force_manifold (),
                                          dialog.samples (),
                                          dialog.iterations ());

      for (std::size_t i = 0; i < reco_items.size (); ++ i)
        {
          if (dialog.force_manifold () && i > reco_items.size () - 3)
            {
              if (dialog.generate_smoothed () && i % 2)
                reco_items[i]->setName(tr("%1 (scale space smooth garbage)").arg(scene->item(index)->name()));
              else
                reco_items[i]->setName(tr("%1 (scale space garbage)").arg(scene->item(index)->name()));
            }
          else
            {
              if (dialog.generate_smoothed ())
                {
                  if (i % 2)
                    reco_items[i]->setName(tr("%1 (scale space smooth shell %2)").arg(scene->item(index)->name()).arg((i+1)/2));
                  else
                    reco_items[i]->setName(tr("%1 (scale space shell %2)").arg(scene->item(index)->name()).arg(((i+1)/2)+1));
                }
              else
                reco_items[i]->setName(tr("%1 (scale space shell %2)").arg(scene->item(index)->name()).arg(i+1));
            }
          scene->addItem (reco_items[i]);
        }

      QApplication::restoreOverrideCursor();
    }
}


void Polyhedron_demo_surface_reconstruction_plugin::poisson_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(point_set_item)
    {
      // Gets point set
      Point_set* points = point_set_item->point_set();
      if(!points) return;

      const double sm_angle     = dialog.angle ();
      const double sm_radius    = dialog.radius ();
      const double sm_distance  = dialog.distance ();
      const QString sm_solver   = dialog.solver();
      bool use_two_passes = dialog.two_passes();
      bool do_not_fill_holes = dialog.do_not_fill_holes();

      QApplication::setOverrideCursor(Qt::WaitCursor);

      if (!(point_set_item->has_normals()))
        {
          std::cerr << "Estimation of normal vectors... ";

          SurfaceReconstruction::compute_normals (*points, 12);
		  
          point_set_item->set_has_normals (true);
          point_set_item->setRenderingMode(PointsPlusNormals);

        }


      // Reconstruct point set as a polyhedron
      Polyhedron* pRemesh = poisson_reconstruct(*points, sm_angle, sm_radius, sm_distance, sm_solver, use_two_passes,
                                                do_not_fill_holes);
      if(pRemesh)
        {
          // Add polyhedron to scene
          Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pRemesh);
          new_item->setName(tr("%1 Poisson (%2 %3 %4)")
                            .arg(point_set_item->name())
                            .arg(sm_angle)
                            .arg(sm_radius)
                            .arg(sm_distance));
          new_item->setColor(Qt::lightGray);
          scene->addItem(new_item);


          // Hide point set
          point_set_item->setVisible(false);
          scene->itemChanged(index);
        }

      QApplication::restoreOverrideCursor();
    }
}


#include "Surface_reconstruction_plugin.moc"
