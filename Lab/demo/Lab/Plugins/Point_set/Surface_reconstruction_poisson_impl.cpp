//----------------------------------------------------------
// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
//----------------------------------------------------------

// CGAL
#include <CGAL/AABB_tree.h> // must be included before kernel
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Eigen_solver_traits.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <math.h>

#include "Kernel_type.h"
#include "SMesh_type.h"
#include "Scene_points_with_normal_item.h"

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

template <typename Triangulation>
class Marching_tets
{
private:
  typedef typename Triangulation::FT FT;
  typedef typename Triangulation::Point Point_3;
  typedef typename Triangulation::Vector Vector_3;
  typedef std::array<std::size_t, 3> Facet;

  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Triangulation::Finite_cells_iterator Finite_cells_iterator;

  const Triangulation& m_tr;

public:

  Marching_tets (const Triangulation& tr) : m_tr (tr) { }

  template <typename VertexOutputIterator,
            typename FacetOutputIterator>
  void output_to_polygon_soup (VertexOutputIterator vertices, FacetOutputIterator facets, FT value = 0.) const
  {
    std::vector<std::array<Point_3, 3> > cell_facets;

    std::size_t nb_points = 0;
    std::map<Point_3, std::size_t> map_p2i;

    for (Finite_cells_iterator cit = m_tr.finite_cells_begin();
         cit != m_tr.finite_cells_end(); ++ cit)
    {
      contour (cit, value, cell_facets);

      for (const std::array<Point_3, 3>& f : cell_facets)
      {
        std::array<std::size_t, 3> facet;
        std::size_t idx = 0;
        for (const Point_3& p : f)
        {
          typename std::map<Point_3, std::size_t>::iterator iter;
          bool inserted = false;
          std::tie (iter, inserted) = map_p2i.insert (std::make_pair (p, nb_points));
          if (inserted)
          {
            *(vertices ++) = p;
            ++ nb_points;
          }
          facet[idx ++] = iter->second;
        }
        *(facets ++) = facet;
      }

      cell_facets.clear();
    }
  }

private:

  void contour (Cell_handle cell, const FT value, std::vector<std::array<Point_3, 3> >& facets) const
  {
    std::vector<Point_3> cell_points;
    Vector_3 direction;

    if(!extract_level_set_points (cell, value, cell_points, direction))
      return;

    if(cell_points.size() == 3)
    {
      const Point_3& a = cell_points[0];
      const Point_3& b = cell_points[1];
      const Point_3& c = cell_points[2];

      Vector_3 n = CGAL::cross_product((b - a), (c - a));

      if(n * direction >= 0)
        facets.push_back ({a, b, c});
      else
        facets.push_back ({a, c, b});
    }
    else if(cell_points.size() == 4)
    {
      // compute normal
      Vector_3 u = cell_points[1] - cell_points[0];
      Vector_3 v = cell_points[2] - cell_points[0];
      Vector_3 n = CGAL::cross_product(u, v);

      if(n * direction <= 0)
      {
        facets.push_back ({cell_points[0], cell_points[2], cell_points[3]});
        facets.push_back ({cell_points[0], cell_points[3], cell_points[1]});
      }
      else
      {
        facets.push_back ({cell_points[0], cell_points[1], cell_points[3]});
        facets.push_back ({cell_points[0], cell_points[3], cell_points[2]});
      }
    }
  }

  bool extract_level_set_points(Cell_handle cell, const FT value, std::vector<Point_3>& points, Vector_3& direction) const
  {
    Point_3 point;
    if(level_set(cell, value, 0, 1, point, direction)) points.push_back(point);
                if(level_set(cell, value, 0, 2, point, direction)) points.push_back(point);
                if(level_set(cell, value, 0, 3, point, direction)) points.push_back(point);
                if(level_set(cell, value, 1, 2, point, direction)) points.push_back(point);
                if(level_set(cell, value, 1, 3, point, direction)) points.push_back(point);
                if(level_set(cell, value, 2, 3, point, direction)) points.push_back(point);
                return points.size() != 0;
  }

  bool level_set(Cell_handle cell, const FT value, const int i1, const int i2, Point_3& p, Vector_3& direction) const
  {
    const Point_3& p1 = cell->vertex(i1)->point();
    const Point_3& p2 = cell->vertex(i2)->point();
    double v1 = cell->vertex(i1)->f();
    double v2 = cell->vertex(i2)->f();

    if(v1 <= value && v2 >= value)
    {
      double ratio = (value - v1) / (v2 - v1);
      p = p1 + ratio * (p2 - p1);
      direction = p2 - p1;
      return true;
    }
    else if(v2 <= value && v1 >= value)
    {
      double ratio = (value - v2) / (v1 - v2);
      p = p2 + ratio * (p1 - p2);
      direction = p1 - p2;
      return true;
    }
    return false;
  }

};


// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
SMesh* poisson_reconstruct(Point_set& points,
                           bool marching_tets,
                           Kernel::FT sm_angle, // Min triangle angle (degrees).
                           Kernel::FT sm_radius, // Max triangle size w.r.t. point set average spacing.
                           Kernel::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                           bool conjugate_gradient,
                           bool use_two_passes,
                           bool do_not_fill_holes)
{
  // Poisson implicit function
  typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

  // Surface mesher
  typedef CGAL::Surface_mesh_default_triangulation_3 STr;
  typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
  typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

  // AABB tree
  typedef CGAL::AABB_face_graph_triangle_primitive<SMesh> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits> AABB_tree;
  CGAL::Timer task_timer; task_timer.start();

  //***************************************
  // Checks requirements
  //***************************************

  if (points.size() == 0)
  {
    std::cerr << "Error: empty point set" << std::endl;
    return nullptr;
  }

  bool points_have_normals = points.has_normal_map();
  if ( ! points_have_normals )
  {
    std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
    return nullptr;
  }

  CGAL::Timer reconstruction_timer; reconstruction_timer.start();

  //***************************************
  // Computes implicit function
  //***************************************


  std::cerr << "Computes Poisson implicit function "
            << "using " << (conjugate_gradient ? "Conjugate Gradient" : "Simplicial LDLT") << "..." << std::endl;


  // Creates implicit function from the point set.
  // Note: this method requires an iterator over points
  // + property maps to access each point's position and normal.
  Poisson_reconstruction_function function(points.begin_or_selection_begin(), points.end(),
                                           points.point_map(), points.normal_map());

  bool ok = false;
  if(conjugate_gradient)
  {
    CGAL::Eigen_solver_traits<Eigen::SimplicialCholesky<CGAL::Eigen_sparse_matrix<double>::EigenType> > solver;
    ok = function.compute_implicit_function(solver, use_two_passes);
  }
  else
  {
    CGAL::Eigen_solver_traits<Eigen::ConjugateGradient<CGAL::Eigen_sparse_matrix<double>::EigenType> > solver;
    solver.solver().setTolerance(1e-6);
    solver.solver().setMaxIterations(1000);
    ok = function.compute_implicit_function(solver, use_two_passes);
  }

  // Computes the Poisson indicator function f()
  // at each vertex of the triangulation.
  if ( ! ok )
  {
    std::cerr << "Error: cannot compute implicit function" << std::endl;
    return nullptr;
  }

  // Prints status
  std::cerr << "Total implicit function (triangulation+refinement+solver): " << task_timer.time() << " seconds\n";
  task_timer.reset();

  SMesh* mesh = new SMesh;

  if (marching_tets)
  {

    std::cerr << "Marching tetrahedra..." << std::endl;

    std::vector<Point_3> vertices;
    std::vector<std::array<std::size_t, 3> > facets;
    Marching_tets<typename Poisson_reconstruction_function::Triangulation> marching_tets (function.tr());
    marching_tets.output_to_polygon_soup (std::back_inserter(vertices), std::back_inserter(facets));

    if (!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(facets))
    {
      std::cerr << "Error: result is not oriented" << std::endl;
      delete mesh;
      return nullptr;
    }
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(vertices, facets, *mesh);
  }
  else
  {
    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Computes average spacing
    Kernel::FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(points.all_or_selection_if_not_empty(),
                                                                                6 /* knn = 1 ring */,
                                                                                points.parameters());

    // Gets one point inside the implicit surface
    Kernel::Point_3 inner_point = function.get_inner_point();
    Kernel::FT inner_point_value = function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      delete mesh;
      return nullptr;
    }

    // Gets implicit function's radius
    Kernel::Sphere_3 bsphere = function.bounding_sphere();
    Kernel::FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    Kernel::FT sm_sphere_radius = 5.0 * radius;
    Kernel::FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Kernel::Sphere_3(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*average_spacing,  // Max triangle size
                                                        sm_distance*average_spacing); // Approximation error

    CGAL_TRACE_STREAM << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius<<",\n"
                      << "                    angle="<<sm_angle << " degrees,\n"
                      << "                    triangle size="<<sm_radius<<" * average spacing="<<sm_radius*average_spacing<<",\n"
                      << "                    distance="<<sm_distance<<" * average spacing="<<sm_distance*average_spacing<<",\n"
                      << "                    dichotomy error=distance/"<<sm_distance*average_spacing/sm_dichotomy_error<<",\n"
                      << "                    Manifold_with_boundary_tag)\n";

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    // Prints status
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
              << tr.number_of_vertices() << " output vertices"
              << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0)
    {
      delete mesh;
      return nullptr;
    }

    // Converts to polyhedron
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, *mesh);

    // Prints total reconstruction duration
    std::cerr << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";
  }

  //***************************************
  // Computes reconstruction error
  //***************************************

  // Constructs AABB tree and computes internal KD-tree
  // data structure to accelerate distance queries
  AABB_tree tree(faces(*mesh).first, faces(*mesh).second, *mesh);

  // Computes distance from each input point to reconstructed mesh
  double max_distance = DBL_MIN;
  double avg_distance = 0;

  std::set<typename boost::graph_traits<SMesh>::face_descriptor> faces_to_keep;

  for (Point_set::const_iterator p=points.begin_or_selection_begin(); p!=points.end(); p++)
  {
    typename AABB_traits::Point_and_primitive_id pap = tree.closest_point_and_primitive (points.point (*p));
    double distance = std::sqrt(CGAL::squared_distance (pap.first, points.point(*p)));

    max_distance = (std::max)(max_distance, distance);
    avg_distance += distance;

    typename boost::graph_traits<SMesh>::face_descriptor f = pap.second;
    faces_to_keep.insert (f);
  }
  avg_distance /= double(points.size());

  std::cerr << "Reconstruction error:" << std::endl
            << "  max = " << max_distance << std::endl
            << "  avg = " << avg_distance << std::endl;

  if (do_not_fill_holes)
  {
    typename boost::graph_traits<SMesh>::face_iterator it = faces(*mesh).begin ();
    while (it != faces(*mesh).end ())
          {
      typename boost::graph_traits<SMesh>::face_iterator current = it ++;

      if (faces_to_keep.find (*current) == faces_to_keep.end ())
      {
        CGAL::Euler::remove_face(halfedge (*current, *mesh), *mesh);
      }

          }

  }
  return mesh;
}

