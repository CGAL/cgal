#include <QTime>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QObject>

#include <fstream>

#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include "SMesh_type.h"
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
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Alpha_shape_mesher.h>
#include <CGAL/Scale_space_reconstruction_3/Weighted_PCA_smoother.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/structure_point_set.h>

#include "ui_Surface_reconstruction_plugin.h"
#include "CGAL/Kernel_traits.h"

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
// types for K nearest neighbors search
typedef CGAL::Search_traits_3<Kernel> Tree_traits;
typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Neighbor_search::iterator Search_iterator;

typedef CGAL::Scale_space_surface_reconstruction_3<Kernel> ScaleSpace;
typedef CGAL::Scale_space_reconstruction_3::Advancing_front_mesher<Kernel> ScaleSpaceAFM;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher<Kernel> ScaleSpaceASM;
typedef CGAL::Scale_space_reconstruction_3::Jet_smoother<Kernel> ScaleSpaceJS;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother<Kernel> ScaleSpaceWPS;

typedef CGAL::cpp11::array<std::size_t,3> Facet;
template<class Mesh, typename Traits>
struct Construct{
  typedef CGAL::cpp11::array<std::size_t,3> Facet;
  typedef typename Traits::Point_3  Point_3;
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VPmap;
  Mesh& mesh;
  VPmap vpmap;
  std::vector<typename boost::graph_traits<Mesh>::vertex_descriptor> vertices;
  template < typename PointIterator>
  Construct(Mesh& mesh,PointIterator b, PointIterator e)
    : mesh(mesh)
  {
    vpmap = get(boost::vertex_point, mesh);
    for(; b!=e; ++b){
      typename boost::graph_traits<Mesh>::vertex_descriptor v;
      v = add_vertex(mesh);
      vertices.push_back(v);
      put(vpmap, v,  *b);
    }
  }

  Construct& operator=(const Facet f)
  {
    typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
    std::vector<vertex_descriptor> facet;
    facet.resize(3);
    facet[0]=vertices[f[0]];
    facet[1]=vertices[f[1]];
    facet[2]=vertices[f[2]];
    CGAL::Euler::add_face(facet, mesh);
    return *this;
  }
  Construct&
  operator*() { return *this; }
  Construct&
  operator++() { return *this; }
  Construct
  operator++(int) { return *this; }
};


// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
Polyhedron* poisson_reconstruct_polyhedron(Point_set& points,
                                Kernel::FT sm_angle, // Min triangle angle (degrees).
                                Kernel::FT sm_radius, // Max triangle size w.r.t. point set average spacing.
                                Kernel::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                                const QString& solver_name, // solver name
                                bool use_two_passes,
                                bool do_not_fill_holes);
// Reconstructs a surface mesh from a point set and returns it.
SMesh* poisson_reconstruct_sm(Point_set& points,
                                Kernel::FT sm_angle, // Min triangle angle (degrees).
                                Kernel::FT sm_radius, // Max triangle size w.r.t. point set average spacing.
                                Kernel::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                                const QString& solver_name, // solver name
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


template <typename Structuring>
struct Priority_with_structure_coherence {

  Structuring& structuring;
  double bound;
  
  Priority_with_structure_coherence(Structuring& structuring,
                                    double bound)
    : structuring (structuring), bound (bound)
  {}

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // If perimeter > bound, return infinity so that facet is not used
    if (bound != 0)
      {
        double d  = 0;
        d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                  c->vertex((index+2)%4)->point()));
        if(d>bound) return adv.infinity();
        d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                                   c->vertex((index+3)%4)->point()));
        if(d>bound) return adv.infinity();
        d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                   c->vertex((index+3)%4)->point()));
        if(d>bound) return adv.infinity();
      }

    Facet f = {{ c->vertex ((index + 1) % 4)->info (),
                 c->vertex ((index + 2) % 4)->info (),
                 c->vertex ((index + 3) % 4)->info () }};

    double weight = 100. * (5 - structuring.facet_coherence (f));

    return weight * adv.smallest_radius_delaunay_sphere (c, index);
  }

};


struct On_the_fly_pair{
  const Point_set& points;
  typedef std::pair<Point, std::size_t> result_type;

  On_the_fly_pair(const Point_set& points) : points(points) {}
  
  result_type
  operator()(const Point_set::Index& i) const
  {
    return result_type(points.point(i), (std::size_t)i);
  }

};


namespace SurfaceReconstruction
{
  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  // types for K nearest neighbors search
  typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
  typedef CGAL::Search_traits_adapter <Point_set::Index, Point_set::Point_map, SearchTraits_3> Search_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  typedef Neighbor_search::Distance Distance;
  typedef Neighbor_search::iterator Search_iterator;

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
    Tree tree(points.begin_or_selection_begin(), points.end(),
              Tree::Splitter(), Search_traits (points.point_map()));
    
    double ratio_kept = (points.size() < 1000)
      ? 1. : 1000. / (points.size());
    
    std::vector<Point> subset;
    for (Point_set::const_iterator it = points.begin(); it != points.end(); ++ it)
      if (CGAL::get_default_random().get_double() < ratio_kept)
    	subset.push_back (points.point(*it));
    
    std::vector<unsigned int> scales;
    generate_scales (std::back_inserter (scales));

    std::vector<unsigned int> chosen;
    Distance tr_dist (points.point_map());
    
    for (std::size_t i = 0; i < subset.size (); ++ i)
      {
    	Neighbor_search search(tree, subset[i],scales.back(), 0, true, tr_dist);
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
    	Neighbor_search search(tree, subset[i], aniso_scale, 0, true, tr_dist);
	size += std::sqrt ((-- search.end())->second);
      }
    size /= subset.size();
    
    return aniso_scale;
  }

  
  unsigned int scale_of_noise (const Point_set& points, double& size)
  {
    Tree tree(points.begin_or_selection_begin(), points.end(),
              Tree::Splitter(), Search_traits (points.point_map()));
    Distance tr_dist (points.point_map());
    
    double ratio_kept = (points.size() < 1000)
      ? 1. : 1000. / (points.size());
    
    std::vector<Point> subset;
    for (Point_set::const_iterator it = points.begin(); it != points.end(); ++ it)
      if (CGAL::get_default_random().get_double() < ratio_kept)
    	subset.push_back (points.point(*it));
    
    std::vector<unsigned int> scales;
    generate_scales (std::back_inserter (scales));

    std::vector<unsigned int> chosen;
    
    for (std::size_t i = 0; i < subset.size (); ++ i)
      {
    	Neighbor_search search(tree, subset[i],scales.back(), 0, true, tr_dist);
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
    	Neighbor_search search(tree, subset[i], noise_scale, 0, true, tr_dist);
	size += std::sqrt ((-- search.end())->second);
      }
    size /= subset.size();

    
    return noise_scale;
  }

  void simplify_point_set (Point_set& points, double size)
  {
    points.set_first_selected (CGAL::grid_simplify_point_set (points.begin (), points.end (), points.point_map(), size));
    points.delete_selection();
  }

  void smooth_point_set (Point_set& points, unsigned int scale)
  {
    CGAL::jet_smooth_point_set<Concurrency_tag>(points.begin(), points.end(), points.point_map(),
                                                scale);
  }

  template <typename ItemsInserter>
  void scale_space (const Point_set& points, ItemsInserter items,
                    bool jet_smoother, 
                    unsigned int iterations,
                    unsigned int neighbors, unsigned int fitting, unsigned int monge,
                    unsigned int neighborhood_size, unsigned int samples,
                    bool advancing_front_mesher,
                    bool generate_smooth,
                    double longest_edge, double radius_ratio_bound, double beta_angle,
                    bool separate_shells, bool force_manifold)
  {
    ScaleSpace reconstruct (points.points().begin(), points.points().end());

    double squared_radius = 0.;
    if (jet_smoother)
    {
      ScaleSpaceJS smoother(neighbors, fitting, monge);
      reconstruct.increase_scale(iterations, smoother);
      if (!advancing_front_mesher)
        squared_radius = CGAL::compute_average_spacing<Concurrency_tag> (points.points().begin(),
                                                                         points.points().end(), neighbors);
    }
    else
    {
      ScaleSpaceWPS smoother(neighborhood_size, samples);
      reconstruct.increase_scale(iterations, smoother);
      squared_radius = smoother.squared_radius();
      
    }

    if (advancing_front_mesher)
    {
      ScaleSpaceAFM mesher (longest_edge, radius_ratio_bound, beta_angle);
      reconstruct.reconstruct_surface (mesher);

      Scene_polygon_soup_item* new_item
        = new Scene_polygon_soup_item ();
      new_item->setColor(Qt::lightGray);
      new_item->setRenderingMode(FlatPlusEdges);
      new_item->init_polygon_soup(points.size(), reconstruct.number_of_facets ());

      Scene_polygon_soup_item* smooth_item = NULL;
      if (generate_smooth)
      {
        smooth_item = new Scene_polygon_soup_item ();
        smooth_item->setColor(Qt::lightGray);
        smooth_item->setRenderingMode(FlatPlusEdges);
        smooth_item->init_polygon_soup(points.size(), reconstruct.number_of_facets ());
      }

      std::map<std::size_t, std::size_t> map_i2i;
      std::size_t current_index = 0;
    
      for (ScaleSpace::Facet_iterator it = reconstruct.facets_begin();
           it != reconstruct.facets_end(); ++ it)
      {
        for (unsigned int ind = 0; ind < 3; ++ ind)
        {
          if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
          {
            map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
            Point p = points.point(*(points.begin_or_selection_begin() + (*it)[ind]));
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
    else
    {
      ScaleSpaceASM mesher (squared_radius, separate_shells, force_manifold);
      reconstruct.reconstruct_surface (mesher);

      for( unsigned int sh = 0; sh < mesher.number_of_shells(); ++sh )
      {
        Scene_polygon_soup_item* new_item
          = new Scene_polygon_soup_item ();
        new_item->setColor(Qt::lightGray);
        new_item->setRenderingMode(FlatPlusEdges);
        new_item->init_polygon_soup(points.size(), mesher.number_of_triangles ());

        Scene_polygon_soup_item* smooth_item = NULL;
        if (generate_smooth)
        {
          smooth_item = new Scene_polygon_soup_item ();
          smooth_item->setColor(Qt::lightGray);
          smooth_item->setRenderingMode(FlatPlusEdges);
          smooth_item->init_polygon_soup(points.size(), mesher.number_of_triangles ());
        }

        std::map<unsigned int, unsigned int> map_i2i;
        unsigned int current_index = 0;
    
        for (ScaleSpaceASM::Facet_iterator it = mesher.shell_begin (sh);
             it != mesher.shell_end (sh); ++ it)
        {
          for (unsigned int ind = 0; ind < 3; ++ ind)
          {
            if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
            {
              map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
              Point p = points.point(*(points.begin_or_selection_begin() + (*it)[ind]));
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
        std::ptrdiff_t num = std::distance( mesher.garbage_begin(  ),
                                            mesher.garbage_end(  ) );

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

        std::map<std::size_t, std::size_t> map_i2i;

        std::size_t current_index = 0;
        for (ScaleSpaceASM::Facet_iterator it=mesher.garbage_begin(),
               end=mesher.garbage_end();it!=end;++it)
        {
          for (unsigned int ind = 0; ind < 3; ++ ind)
          {
            if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
            {
              map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
              Point p = points.point(*(points.begin_or_selection_begin() + (*it)[ind]));
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
  }

  struct Point_set_make_pair_point_index
    : public std::unary_function<const Point_set::Index&, std::pair<Kernel::Point_3, std::size_t> >
  {
    const Point_set& point_set;
    Point_set_make_pair_point_index (const Point_set& point_set) : point_set (point_set) { }
    std::pair<Kernel::Point_3, std::size_t> operator() (const Point_set::Index& i) const
    {
      return std::make_pair (point_set.point (i), i);
    }
  };
  

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

  template<class FaceGraphItem>
  void advancing_front (const Point_set& points, FaceGraphItem* new_item, double size,
                        double radius_ratio_bound = 5., double beta = 0.52)
  {

    // TODO: build DT with indices
    typedef typename FaceGraphItem::Face_graph FaceGraph;
    typedef typename CGAL::Kernel_traits<typename boost::property_traits<typename boost::property_map<FaceGraph,
        typename boost::vertex_point_t>::type>::value_type>::Kernel Traits;

    FaceGraph& P = * const_cast<FaceGraph*>(new_item->face_graph());
    Radius filter(size);
    Construct<FaceGraph, Traits> construct(P,points.points().begin(),points.points().end());
    CGAL::advancing_front_surface_reconstruction(points.points().begin(),
                                                 points.points().end(),
                                                 construct,
                                                 filter,
                                                 radius_ratio_bound,
                                                 beta);
						  
  }

  void compute_normals (Point_set& points, unsigned int neighbors)
  {
    CGAL::jet_estimate_normals<Concurrency_tag>(points.begin_or_selection_begin(), points.end(),
                                                points.point_map(),
                                                points.normal_map(),
                                                2 * neighbors);

    points.set_first_selected (CGAL::mst_orient_normals (points.begin(), points.end(),
                                                         points.point_map(),
                                                         points.normal_map(),
                                                         2 * neighbors));
    points.delete_selection();
  }
  
  struct build_from_pair
  {
    Point_set& m_pts;

    build_from_pair (Point_set& pts) : m_pts (pts) { }

    void operator() (const std::pair<Point_set::Point, Point_set::Vector>& pair)
    {
      m_pts.insert (pair.first, pair.second);
    }


  };
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
    return tabWidget->currentIndex();
  }
  
  // Auto
  bool boundaries () const { return m_boundaries->isChecked (); }
  bool interpolate () const { return m_interpolate->isChecked (); }

  // Advancing front
  double longest_edge () const { return m_longestEdge->value (); }
  double radius_ratio_bound () const { return m_radiusRatioBound->value (); }
  double beta_angle () const { return m_betaAngle->value (); }

  // Scale Space
  bool scalespace_js() const { return m_scalespace_jet->isChecked(); }
  unsigned int iterations () const { return m_iterations->value (); }
  unsigned int neighbors () const { return m_neighbors->value(); }
  unsigned int fitting () const { return m_fitting->value(); }
  unsigned int monge () const { return m_monge->value(); }
  unsigned int neighborhood_size () const { return m_neighborhood_size->value (); }
  unsigned int samples () const { return m_samples->value (); }
  bool scalespace_af() const { return m_scalespace_af->isChecked(); }
  bool generate_smoothed () const { return m_genSmooth->isChecked (); }
  double longest_edge_2 () const { return m_longestEdge_2->value (); }
  double radius_ratio_bound_2 () const { return m_radiusRatioBound_2->value (); }
  double beta_angle_2 () const { return m_betaAngle_2->value (); }
  bool separate_shells () const { return m_genShells->isChecked (); }
  bool force_manifold () const { return m_forceManifold->isChecked (); }


  // Poisson
  double angle () const { return m_inputAngle->value (); }
  double radius () const { return m_inputRadius->value (); }
  double distance () const { return m_inputDistance->value (); }
  bool two_passes () const { return m_inputTwoPasses->isChecked (); }
  bool do_not_fill_holes () const { return m_doNotFillHoles->isChecked (); }

  // RANSAC
  double connectivity_tolerance () const { return m_connectivityTolerance->value (); }
  double noise_tolerance () const { return m_noiseTolerance->value (); }
  unsigned int min_size_subset () const { return m_minSizeSubset->value (); }
  bool generate_structured () const { return m_generateStructured->isChecked (); }
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
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionSurfaceReconstruction = new QAction(tr("Surface Reconstruction"), mainWindow);
    actionSurfaceReconstruction->setObjectName("actionSurfaceReconstruction");
    autoConnectActions();

  }

  void automatic_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void advancing_front_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void scale_space_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void poisson_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void ransac_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  
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
      dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);
      if(!dialog.exec())
	return;

      unsigned int method = dialog.method ();
      switch (method)
        {
        case 0:
          advancing_front_reconstruction (dialog);
          break;
        case 1:
          poisson_reconstruction (dialog);
          break;
        case 2:
          scale_space_reconstruction (dialog);
          break;
        case 3:
          ransac_reconstruction (dialog);
          break;
        case 4:
          automatic_reconstruction (dialog);
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
	  new_item = new Scene_points_with_normal_item(*pts_item);
	  new_item->setName(QString("%1 (preprocessed)").arg(pts_item->name()));
	  new_item->resetSelection();
	  new_item->invalidateOpenGLBuffers();

	  points = new_item->point_set();
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
          new_item->point_set()->remove_normal_map();
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
                                                  true,
                                                  4,
                                                  (std::max)(noise_scale, aniso_scale), 2, 2,
                                                  0, 0,
                                                  true,
                                                  false,
                                                  10. * (std::max)(noise_size, aniso_size), 5., 0.52,
                                                  false, false);
	      
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

              if(mw->property("is_polyhedron_mode").toBool())
              {
                Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
                SurfaceReconstruction::advancing_front (*points, reco_item, 10. * (std::max)(noise_size, aniso_size));

                reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
                reco_item->setColor(Qt::lightGray);
                reco_item->setRenderingMode(FlatPlusEdges);
                scene->addItem(reco_item);
              }
              else
              {
                Scene_surface_mesh_item* reco_item = new Scene_surface_mesh_item(SMesh());
                SurfaceReconstruction::advancing_front (*points, reco_item, 10. * (std::max)(noise_size, aniso_size));

                reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
                reco_item->setColor(Qt::lightGray);
                reco_item->setRenderingMode(FlatPlusEdges);
                scene->addItem(reco_item);
              }
	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }

	}
      else
	{
	  if (dialog.boundaries())
	    {
	      std::cerr << "Advancing front reconstruction... ";
	      time.restart();

              if(mw->property("is_polyhedron_mode").toBool())
              {
                Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
                SurfaceReconstruction::advancing_front (*points, reco_item, 10. * (std::max)(noise_size, aniso_size));

                reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
                reco_item->setColor(Qt::lightGray);
                reco_item->setRenderingMode(FlatPlusEdges);
                scene->addItem(reco_item);
              }
              else
              {
                Scene_surface_mesh_item* reco_item = new Scene_surface_mesh_item(SMesh());
                SurfaceReconstruction::advancing_front (*points, reco_item, 10. * (std::max)(noise_size, aniso_size));

                reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
                reco_item->setColor(Qt::lightGray);
                reco_item->setRenderingMode(FlatPlusEdges);
                scene->addItem(reco_item);
              }

	      std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      if (!(new_item->has_normals()))
		{
		  std::cerr << "Estimation of normal vectors... ";
		  time.restart();

		  SurfaceReconstruction::compute_normals (*points, noise_scale);
		  
		  new_item->point_set()->add_normal_map();
		  new_item->setRenderingMode(PointsPlusNormals);
		  
		  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
		}
	      
	      std::cerr << "Poisson reconstruction... ";
              time.restart();
              Polyhedron* pRemesh = NULL;
              SMesh* smRemesh = NULL;
              if(mw->property("is_polyhedron_mode").toBool())
                pRemesh = poisson_reconstruct_polyhedron(*points,
                                                         20,
                                                         100 * (std::max)(noise_size, aniso_size),
                                                         (std::max)(noise_size, aniso_size),
                                                         QString ("Eigen - built-in CG"), false, false);
              else
                smRemesh = poisson_reconstruct_sm(*points,
                                                  20,
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
              else if(smRemesh)
              {
                // Add polyhedron to scene
                Scene_surface_mesh_item* reco_item = new Scene_surface_mesh_item(smRemesh);
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

      if(mw->property("is_polyhedron_mode").toBool())
      {
        Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
        SurfaceReconstruction::advancing_front (*points, reco_item,
                                                dialog.longest_edge (),
                                                dialog.radius_ratio_bound (),
                                                CGAL_PI * dialog.beta_angle () / 180.);

        reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
        reco_item->setColor(Qt::lightGray);
        reco_item->setRenderingMode(FlatPlusEdges);
        scene->addItem(reco_item);
      }
      else
      {
        Scene_surface_mesh_item* reco_item = new Scene_surface_mesh_item(SMesh());
        SurfaceReconstruction::advancing_front (*points, reco_item,
                                                dialog.longest_edge (),
                                                dialog.radius_ratio_bound (),
                                                CGAL_PI * dialog.beta_angle () / 180.);

        reco_item->setName(tr("%1 (advancing front)").arg(scene->item(index)->name()));
        reco_item->setColor(Qt::lightGray);
        reco_item->setRenderingMode(FlatPlusEdges);
        scene->addItem(reco_item);
      }

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
                                          dialog.scalespace_js(),
                                          dialog.iterations(),
                                          dialog.neighbors(), dialog.fitting(), dialog.monge(),
                                          dialog.neighborhood_size (), dialog.samples(),
                                          dialog.scalespace_af(),
                                          dialog.generate_smoothed (),
                                          dialog.longest_edge_2(), dialog.radius_ratio_bound_2(),
                                          CGAL_PI * dialog.beta_angle_2 () / 180.,
                                          dialog.separate_shells (), dialog.force_manifold ());

      for (std::size_t i = 0; i < reco_items.size (); ++ i)
        {
          if (!(dialog.scalespace_af()))
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
          }
          else
          {
            if (dialog.generate_smoothed ())
            {
              if (i % 2)
                reco_items[i]->setName(tr("%1 (scale space smooth)").arg(scene->item(index)->name()));
              else
                reco_items[i]->setName(tr("%1 (scale space)").arg(scene->item(index)->name()));
            }
            else
              reco_items[i]->setName(tr("%1 (scale space)").arg(scene->item(index)->name()));
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
          points->add_normal_map();
          SurfaceReconstruction::compute_normals (*points, 12);
		  

          point_set_item->setRenderingMode(PointsPlusNormals);

        }


      // Reconstruct point set as a polyhedron
      Polyhedron* pRemesh = NULL;
      SMesh* smRemesh= NULL;
      if(mw->property("is_polyhedron_mode").toBool())
        pRemesh = poisson_reconstruct_polyhedron(*points, sm_angle, sm_radius, sm_distance, sm_solver, use_two_passes,
                                                 do_not_fill_holes);
      else
        smRemesh = poisson_reconstruct_sm(*points, sm_angle, sm_radius, sm_distance, sm_solver, use_two_passes,
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
      else if(smRemesh)
      {
        // Add polyhedron to scene
        Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(smRemesh);
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

void Polyhedron_demo_surface_reconstruction_plugin::ransac_reconstruction
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

      QApplication::setOverrideCursor(Qt::WaitCursor);

      CGAL::Timer global_timer;
      global_timer.start();

      CGAL::Timer local_timer;

      if (!(point_set_item->has_normals()))
        {
          local_timer.start();
                
          std::cerr << "Estimation of normal vectors... ";
          points->add_normal_map();
          CGAL::jet_estimate_normals<Concurrency_tag>(points->begin(), points->end(),
                                                      points->point_map(),
                                                      points->normal_map(),
                                                      12);
          local_timer.stop();
          point_set_item->setRenderingMode(PointsPlusNormals);

          std::cerr << "done in " << local_timer.time() << " second(s)" << std::endl;
          local_timer.reset();
        }

      typedef Point_set::Point_map PointMap;
      typedef Point_set::Vector_map NormalMap;

      typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Kernel, Point_set, PointMap, NormalMap> Traits;
      typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Shape_detection;

      local_timer.start();
      Shape_detection shape_detection;
      shape_detection.set_input(*points, points->point_map(), points->normal_map());

      shape_detection.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();

      Shape_detection::Parameters op;
      op.probability = 0.05;
      op.min_points = dialog.min_size_subset();
      op.epsilon = dialog.noise_tolerance();
      op.cluster_epsilon = dialog.connectivity_tolerance();
      op.normal_threshold = 0.7;

      shape_detection.detect(op);
      local_timer.stop();
      std::cout << shape_detection.shapes().size() << " plane(s) found in "
                << local_timer.time() << " second(s)" << std::endl;
      local_timer.reset();
      
      std::cout << "Structuring point set... " << std::endl;
      typedef CGAL::Point_set_with_structure<Traits> Structuring;

      local_timer.start();
      Structuring structuring (points->begin (), points->end (),
                               points->point_map(), points->normal_map(),
                               shape_detection,
                               op.cluster_epsilon);

      Scene_points_with_normal_item *structured = new Scene_points_with_normal_item;
      structured->point_set()->add_normal_map();
      for (std::size_t i = 0; i < structuring.size(); ++ i)
        structured->point_set()->insert (structuring.point(i), structuring.normal(i));

      local_timer.stop ();
      std::cerr << structured->point_set()->size() << " point(s) generated in "
                << local_timer.time() << std::endl;
      local_timer.reset();

      std::cerr << "Reconstructing... ";
      local_timer.start();

      Priority_with_structure_coherence<Structuring> priority (structuring, 10. * op.cluster_epsilon);

      if(mw->property("is_polyhedron_mode").toBool())
      {
        Scene_polyhedron_item* reco_item = new Scene_polyhedron_item(Polyhedron());
        Polyhedron& P = * const_cast<Polyhedron*>(reco_item->polyhedron());
        Construct<Polyhedron, Traits> construct(P,structured->point_set()->points().begin(),structured->point_set()->points().end());
        CGAL::advancing_front_surface_reconstruction(structured->point_set()->points().begin(),
                                                     structured->point_set()->points().end(),
                                                     construct,
                                                     priority,
                                                     5.,
                                                     0.52);
        local_timer.stop();
        std::cerr << "done in " << local_timer.time() << " second(s)" << std::endl;

        reco_item->setName(tr("%1 (RANSAC-based reconstruction)").arg(scene->item(index)->name()));
        reco_item->setColor(Qt::magenta);
        reco_item->setRenderingMode(FlatPlusEdges);
        scene->addItem(reco_item);
      }
      else
      {
        Scene_surface_mesh_item* reco_item = new Scene_surface_mesh_item(SMesh());
        SMesh& P = * const_cast<SMesh*>(reco_item->polyhedron());
        Construct<SMesh, Traits> construct(P,structured->point_set()->points().begin(),structured->point_set()->points().end());
        CGAL::advancing_front_surface_reconstruction(structured->point_set()->points().begin(),
                                                     structured->point_set()->points().end(),
                                                     construct,
                                                     priority,
                                                     5.,
                                                     0.52);
        local_timer.stop();
        std::cerr << "done in " << local_timer.time() << " second(s)" << std::endl;
        reco_item->setName(tr("%1 (RANSAC-based reconstruction)").arg(scene->item(index)->name()));
        reco_item->setColor(Qt::magenta);
        reco_item->setRenderingMode(FlatPlusEdges);
        scene->addItem(reco_item);
      }
      if (dialog.generate_structured ())
      {
        structured->setName(tr("%1 (structured)").arg(point_set_item->name()));
        structured->setRenderingMode(PointsPlusNormals);
        structured->setColor(Qt::blue);
        scene->addItem (structured);
      }
      else
        delete structured;

      std::cerr << "All done in " << global_timer.time() << " seconds." << std::endl;
      
      QApplication::restoreOverrideCursor();
    }
}


#include "Surface_reconstruction_plugin.moc"
