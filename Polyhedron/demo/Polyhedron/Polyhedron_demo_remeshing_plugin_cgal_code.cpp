#include <CGAL/AABB_polyhedral_oracle.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"
#include "Scene_item.h"
#include <qgl.h>
#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup.h"

#include <CGAL/Simple_cartesian.h>

#include "C2t3_type.h"

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

#include <Qt/qglobal.h>
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

  Surface_mesh_modified_criteria_3(const FT angle_bound,
                                   const FT radius_bound,
                                   const FT distance_bound)
    : curvature_size_criterion(distance_bound),
      uniform_size_criterion(radius_bound),
      aspect_ratio_criterion(angle_bound)
      
  {
  }

  bool is_bad (const Facet& f, Quality& q) const
  {
    const typename Tr::Point& pa = f.first->vertex((f.second+1)%4)->point();
    const typename Tr::Point& pb = f.first->vertex((f.second+2)%4)->point();
    const typename Tr::Point& pc = f.first->vertex((f.second+3)%4)->point();
    int nb_protecting_balls = 0;
    if(pa.weight() != FT(0)) ++nb_protecting_balls;
    if(pb.weight() != FT(0)) ++nb_protecting_balls;
    if(pc.weight() != FT(0)) ++nb_protecting_balls;
    if(nb_protecting_balls == 0)
    { 
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
    }
    else {
      q[0] = 1;
      if(nb_protecting_balls == 3)
        return false;
      else if(uniform_size_criterion.is_bad(f, q[1]))
        return true;
      // else {
      //   q[1] = 1;
      //   if(nb_protecting_balls == 1 && curvature_size_criterion.is_bad(f, q[2]))
      //     return true;
      // }
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
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

//
// Types for meshing
//
typedef CGAL::Simple_cartesian<double> Simple_cartesian_kernel;
// input surface
typedef CGAL::Mesh_3::Robust_intersection_traits_3<Kernel> IGT;
typedef CGAL::AABB_polyhedral_oracle<Polyhedron,IGT,Simple_cartesian_kernel> Input_surface;


// A base non-templated class, to allow 
class Mesher_base : public QObject {
  Q_OBJECT
protected:
  bool is_stopped;
public:
  Mesher_base(QObject* parent) : QObject(parent) {
    is_stopped = true;
  };
  virtual ~Mesher_base() {}
public slots:
  virtual void mesh() = 0;
  virtual void one_step() = 0;

  void stop() {
    std::cerr << "STOP!\n";
    is_stopped = true;
  }
};

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

class Q_DECL_EXPORT Scene_c2t3_item : public Scene_item
{
  Q_OBJECT
public:
  Scene_c2t3_item(const C2t3& c2t3)
    : sphere_display_list(0), quadric(0), c2t3_(c2t3)
  {
  }

  ~Scene_c2t3_item()
  {
    if(quadric != 0)
      gluDeleteQuadric(quadric);
    if(sphere_display_list  != 0)
      glDeleteLists(sphere_display_list, 1);
  }

  const C2t3& c2t3() const {
    return c2t3_;
  }
  bool isFinite() const { return true; }
  bool isEmpty() const {
    return c2t3().triangulation().number_of_vertices() == 0;
  }

  Bbox bbox() const {
    if(isEmpty())
      return Bbox();
    else {
      bool first = true;
      CGAL::Bbox_3 result;
      for(Tr::Finite_vertices_iterator
            vit = ++c2t3().triangulation().finite_vertices_begin(),
            end = c2t3().triangulation().finite_vertices_end();
          vit != end; ++vit)
      {
        if(vit->point().weight() > 0) {
          if(first) {
            result = vit->point().bbox();
            first = false;
          } else { 
            result = result + vit->point().bbox();
          }
        }
      }
      return Bbox(result.xmin(), result.ymin(), result.zmin(),
                  result.xmax(), result.ymax(), result.zmax());
    }
  }


  Scene_c2t3_item* clone() const {
    return 0;
  }

  QString toolTip() const {
    return tr("<p><b>2D complex in a 3D triangulation</b></p>"
              "<p>Number of vertices: %1<br />"
              "Number of surface facets: %2<br />")
      .arg(c2t3().triangulation().number_of_vertices())
      .arg(c2t3().number_of_facets());
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud); // CHECK THIS!
  }

  void draw() const {
    ::glBegin(GL_TRIANGLES);
    for(C2t3::Facet_iterator
          fit = c2t3().facets_begin(),
          end = c2t3().facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Tr::Geom_traits::Point_3& pa = cell->vertex((index+1)&3)->point();
      const Tr::Geom_traits::Point_3& pb = cell->vertex((index+2)&3)->point();
      const Tr::Geom_traits::Point_3& pc = cell->vertex((index+3)&3)->point();
      draw_triangle(pa, pb, pc);
    }
    ::glEnd();
    
    GLenum gl_error = ::glGetError();
    if(gl_error != GL_NO_ERROR)
      std::cerr << "GL error: " << gluErrorString(gl_error) << std::endl;

    if(!draw_spheres)
      return;

    // force wireframe for protecting spheres
    GLint polygon_mode[2];
    ::glGetIntegerv(GL_POLYGON_MODE, &polygon_mode[0]);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    for(Tr::Finite_vertices_iterator 
          vit = c2t3().triangulation().finite_vertices_begin(),
          end =  c2t3().triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      draw_sphere(vit->point());
    }
    ::glPolygonMode(GL_FRONT, polygon_mode[0]);
    ::glPolygonMode(GL_BACK, polygon_mode[1]);
  }

  void draw_sphere(const Tr::Point p) const 
  {
    if(p.weight() > 0) {
      if(sphere_display_list == 0) {
        sphere_display_list = glGenLists(1);
        if(sphere_display_list == 0)
          std::cerr << "ERROR: Cannot create display list!\n";
        if(quadric == 0)
          quadric = gluNewQuadric();
        if(quadric == 0)
          std::cerr << "ERROR: Cannot create GLU quadric!\n";
        glNewList(sphere_display_list, GL_COMPILE);
        gluSphere(quadric, 1., 10, 10);
        glEndList();
        if(glGetError() != GL_NO_ERROR)
          std::cerr << gluErrorString(glGetError());
      }
      glPushMatrix();
      glTranslated(CGAL::to_double(p.point().x()),
                   CGAL::to_double(p.point().y()),
                   CGAL::to_double(p.point().z()));
      const GLdouble r = CGAL::to_double(CGAL_NTS sqrt(p.weight()));
      glScaled(r, r, r);
      glCallList(sphere_display_list);
      glPopMatrix();
    }
  }

public slots:
  void show_spheres(bool b) {
    draw_spheres = b;
  }

private:
  static void draw_triangle(const Tr::Point& pa,
                            const Tr::Point& pb,
                            const Tr::Point& pc) {
    Tr::Geom_traits::Vector_3 n = cross_product(pb - pa, pc -pa);
    n = n / CGAL::sqrt(n*n);

    ::glNormal3d(n.x(),n.y(),n.z());

    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glVertex3d(pc.x(),pc.y(),pc.z());
  }

private:
  mutable GLuint sphere_display_list;
  mutable GLUquadric* quadric;
  C2t3 c2t3_;
  bool draw_spheres;
};

typedef Tr::Geom_traits GT;

bool insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const GT::FT size);

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

  std::cerr << "Insert protecting balls... ";
  timer.reset();
  insert_spheres(c2t3, pMesh, sizing/1.5);
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // initial point set
  timer.reset();
  std::cerr << "Insert initial point set... ";
  typedef CGAL::Cartesian_converter<Kernel,GT> Converter;
  Converter convert;

  { // new scope for the initialization, so that the vector
    // polyhedron_points is destroyed as soon as the initialization is
    // finished
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
      const int pos = CGAL::default_random.get_int(0, polyhedron_points.size());
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
    // // add remesh as new polyhedron
    // Polyhedron *pRemesh = new Polyhedron;
    // CGAL::Complex_2_in_triangulation_3_polyhedron_builder<C2t3,
    // Polyhedron> builder(c2t3);
    return new Scene_c2t3_item(c2t3);

    // try {
    //   CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
    //   pRemesh->delegate(builder);
    // } catch(CGAL::Failure_exception)
    // {
    // }
    // CGAL::set_error_behaviour(CGAL::ABORT);

    // if(c2t3.number_of_facets() != pRemesh->size_of_facets())
    // {
    //   delete pRemesh;
    //   std::stringstream temp_file;
    //   namespace sm = CGAL::Surface_mesher;
    //   if(!CGAL::output_surface_facets_to_off(temp_file, c2t3, 
    //                                          sm::NO_OPTION | sm::IO_VERBOSE))
    //   {
    //     std::cerr << "Cannot write the mesh to an off file!\n";
    //     std::cerr << temp_file.str() << std::endl;
    //     return 0;
    //   }
    //   Scene_polygon_soup* soup = new Scene_polygon_soup();
    //   if(!soup->load(temp_file))
    //   {
    //     std::cerr << "Cannot reload the mesh from an off file!\n";
    //     return 0;
    //   }
    //   else
    //     return soup;
    // } else {
    //   return new Scene_polyhedron_item(pRemesh);
    // }
  }
  else
    return 0;
}

#include "Polyhedron_demo_remeshing_plugin_cgal_code.moc"
