#include <CGAL/AABB_polyhedral_oracle.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"
#include "Scene_item.h"
#include <qgl.h>
#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Simple_cartesian.h>

#include "C2t3_type.h"
#include "Scene_c2t3_item.h"

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>

#include <CGAL/Timer.h>

#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>

#include <algorithm>
#include <sstream>

#include <CGAL/array.h>

#include "Scene_item.h"

#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QApplication>
#include <QThread>
#include <QMessageBox>

template <class Tr>
class Surface_mesh_modified_criteria_3
{
public:
  typedef Tr Triangulation;
  typedef typename Tr::Geom_traits::FT FT;

  typedef typename CGAL::array<FT, 3> Quality;
  typedef typename Tr::Facet Facet;
  const FT sq_distance_bound;

  Surface_mesh_modified_criteria_3(const FT angle_bound,
                                   const FT radius_bound,
                                   const FT distance_bound)
    : sq_distance_bound(distance_bound*distance_bound/100),
      curvature_size_criterion(distance_bound),
      uniform_size_criterion(radius_bound),
      aspect_ratio_criterion(angle_bound)
      
  {
  }

  bool is_bad (const Facet& f, Quality& q) const
  {
    const typename Tr::Point& pa = f.first->vertex((f.second+1)%4)->point();
    const typename Tr::Point& pb = f.first->vertex((f.second+2)%4)->point();
    const typename Tr::Point& pc = f.first->vertex((f.second+3)%4)->point();
    if( squared_distance(pa, pb) < sq_distance_bound )
      return false;
    if( squared_distance(pc, pb) < sq_distance_bound )
      return false;
    if( squared_distance(pa, pc) < sq_distance_bound )
      return false;
    if(aspect_ratio_criterion.is_bad(f, q[0]))
      return true;
    else {
      q[0] = 1;
      if(uniform_size_criterion.is_bad(f, q[1]))
        return true;
      else {
        q[1] = 1;
        if(curvature_size_criterion.is_bad(f, q[2]))
          return true;
      }
    }
    return false;
  }
private:
  CGAL::Surface_mesher::Curvature_size_criterion<Tr> curvature_size_criterion;
  // bound on Hausdorff distance does not play any role if bigger than
  // the square of the Uniform_size_criterion

  CGAL::Surface_mesher::Uniform_size_criterion<Tr> uniform_size_criterion;
  // bound on radii of surface Delaunay balls
  
  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr> aspect_ratio_criterion;
  // lower bound on minimum angle in degrees

}; // end class Surface_mesh_default_criteria_3


namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

//
// Types for meshing
//
typedef CGAL::Simple_cartesian<double> Simple_cartesian_kernel;
// input surface
// typedef CGAL::Mesh_3::Robust_intersection_traits_3<Kernel> IGT;
typedef CGAL::AABB_polyhedral_oracle<Polyhedron,Kernel,Simple_cartesian_kernel> Input_surface;

#include "Mesher_base.h"

// Class template, refines Mesher_base
// That allows to create meshers with different criteria or manifold tag,
// and thread them with the API of Mesher_base (mesh/one_step/stop).
template <typename Criteria, typename Manifold_tag>
class Mesher : public Mesher_base
{
  typedef typename CGAL::Surface_mesher_generator<
    C2t3,
    Input_surface,
    Criteria,
    Manifold_tag,
    CGAL_SURFACE_MESHER_VERBOSITY >::type MyMesher;

  MyMesher mesher;
  const C2t3& c2t3;
  const Input_surface& surface;
  CGAL::Null_mesh_visitor visitor;
public:
  Mesher(QObject* parent, 
         C2t3& c2t3, 
         const Input_surface& surface,
         const Criteria& criteria)
    : Mesher_base(parent), 
      mesher(c2t3, surface, surface, criteria),
      c2t3(c2t3),
      surface(surface)
  {
    typename Input_surface::Construct_initial_points get_initial_points =
      surface.construct_initial_points_object();

    get_initial_points(surface,
                       CGAL::inserter(c2t3.triangulation()),
                       20);
    mesher.init();
  }
        
  void mesh()
  {
    int global_nbsteps = 0;
    int nbsteps = 0;
    CGAL::Timer timer;
    timer.start();
    is_stopped = false;

    std::cerr << "Legende of the following line: "
              << "(#vertices,#steps," << mesher.debug_info_header()
              << ")\n";

    while(!is_stopped && !mesher.is_algorithm_done())
    {
      one_step();
      ++nbsteps;
      ++global_nbsteps;
      if(timer.time() > 1)
      {
        std::cerr 
	  << boost::format("\r             \r"
			   "(%1%,%2%,%3%) (%|4$.1f| vertices/s)")
	  % c2t3.triangulation().number_of_vertices()
	  % global_nbsteps % mesher.debug_info()
	  % (nbsteps / timer.time());
        qApp->processEvents();
        nbsteps = 0;
        timer.reset();
      }
    }
  }

  void one_step()
  {
    mesher.one_step(visitor);
  }
};

// That thread takes a Mesher_base* as parent. It just launches the meshing
// process.
struct Meshing_thread : QThread 
{
  Mesher_base* mesher;

  Meshing_thread(Mesher_base* parent) 
    : QThread(parent), mesher(parent)
  {
  }

  void run() {
    mesher->mesh();
    mesher->moveToThread(QApplication::instance()->thread());
  }
};

typedef Tr::Geom_traits GT;
typedef Tr::Geom_traits::FT FT;

Scene_item* cgal_code_remesh(QWidget* parent, 
                             Polyhedron* pMesh,
                             const double angle,
                             const double sizing,
                             const double approx,
                             int tag) {
// };

// class Mesh_process : public QObject {
//   Q_OBJECT

//   QWidget* parent;
//   Polyhedron* pMesh;
//   const double angle;
//   const double sizing;
//   const double approx;
//   int tag;

// public:
//   Mesh_process(QWidget* parent, 
//                Polyhedron* pMesh,
//                const double angle,
//                const double sizing,
//                const double approx,
//                int tag)
//     : parent(parent),
//       pMesh(pMesh),
//       angle(angle),
//       sizing(sizing),
//       approx(approx)
//       tag(tag)
//   {
//   }

//   Scene_item* launch() 

  if(!pMesh) return 0;

  CGAL::set_error_behaviour(CGAL::ABORT);
  // remesh


  Tr& triangulation = * new Tr;; // 3D-Delaunay triangulation
  C2t3& c2t3 = *(new C2t3(triangulation));
  // C2t3 c2t3(triangulation); // 2D-complex in 3D-Delaunay triangulation

  // meshing parameters
  typedef Surface_mesh_modified_criteria_3<Tr> Criteria;
  const Criteria facets_criteria(angle,sizing,approx);

  // const Criteria new_facets_criteria(facets_criteria);

  // AABB tree
  CGAL::Timer timer;
  timer.start();
  std::cerr << "Build AABB tree...";
  Input_surface input(*pMesh);
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // initial point set
  timer.reset();
  std::cerr << "Insert initial point set... ";

  { // new scope for the initialization, so that the vector
    // polyhedron_points is destroyed as soon as the initialization is
    // finished
    typedef Kernel::Point_3 Point;

    std::vector<Point> polyhedron_points;
    polyhedron_points.reserve(pMesh->size_of_vertices());
    std::copy(pMesh->points_begin(), pMesh->points_end(), 
              std::back_inserter(polyhedron_points));

    typedef std::vector<Point>::size_type size_type;
    size_type nb_initial_points = 10;
    nb_initial_points = (std::min)(nb_initial_points, polyhedron_points.size());
    for(size_type n = 0;
        n < nb_initial_points || (n < 10 * nb_initial_points && 
                                  triangulation.dimension() < 3 );
        n = triangulation.number_of_vertices())
    {
      const int pos = CGAL::default_random.get_int(0, (int)polyhedron_points.size());
      triangulation.insert(polyhedron_points[pos]);
    }
  }
  if(triangulation.dimension() < 3)
    return 0;

  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // return new Scene_c2t3_item(c2t3);
  // remesh
  timer.reset();
  Mesher_base* mesher;
  std::cerr << "Remesh...";
  QMessageBox* message_box = new QMessageBox(QMessageBox::NoIcon,
                                            "Remeshing...",
                                            "Meshing process is running...",
                                            QMessageBox::Cancel,
                                            parent);
  switch(tag) {
  case 0: 
    mesher = new Mesher<Criteria, 
      CGAL::Non_manifold_tag>(0, c2t3, input, facets_criteria);
        ;
    break;
  case 1:
    mesher = new Mesher<Criteria, 
      CGAL::Manifold_tag>(0, c2t3, input, facets_criteria);
    break;
  default:
    mesher = new Mesher<Criteria, 
      CGAL::Manifold_with_boundary_tag>(0, c2t3, input, facets_criteria);
  }
  QObject::connect(message_box, SIGNAL(buttonClicked( QAbstractButton *)),
                   mesher, SLOT(stop()));
  message_box->show();
  qApp->processEvents();

  
  Meshing_thread* thread = new Meshing_thread(mesher);
  mesher->moveToThread(thread);
  thread->start();
  while(!thread->isFinished())
  {
    qApp->processEvents();
    thread->wait(200);
  }
  delete message_box;
  delete mesher;
  std::cerr << "done (" << timer.time() << " ms, " << triangulation.number_of_vertices() << " vertices)" << std::endl;

  if(triangulation.number_of_vertices() > 0)
  {
    // add remesh as new polyhedron
    Polyhedron *pRemesh = new Polyhedron;
    CGAL::Complex_2_in_triangulation_3_polyhedron_builder<C2t3, Polyhedron> builder(c2t3);
    pRemesh->delegate(builder);
    if(c2t3.number_of_facets() != pRemesh->size_of_facets())
    {
      delete pRemesh;
      std::stringstream temp_file;
      if(!CGAL::output_surface_facets_to_off(temp_file, c2t3))
      {
        std::cerr << "Cannot write the mesh to an off file!\n";
        return 0;
      }
      Scene_polygon_soup_item* soup = new Scene_polygon_soup_item();
      if(!soup->load(temp_file))
      {
        std::cerr << "Cannot reload the mesh from an off file!\n";
        return 0;
      }
      else
        return soup;
    } else {
      return new Scene_polyhedron_item(pRemesh);
    }
  }
  else
    return 0;
}
