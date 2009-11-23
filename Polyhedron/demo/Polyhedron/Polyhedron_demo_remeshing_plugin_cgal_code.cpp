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

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>

#include <CGAL/Timer.h>

#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>

#include <algorithm>
#include <sstream>


#include "Scene_item.h"
#include <Qt/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

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
    // draw_sphere(c2t3().triangulation().finite_vertices_begin()->point());
    for(Tr::Finite_vertices_iterator 
          vit = c2t3().triangulation().finite_vertices_begin(),
          end =  c2t3().triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      draw_sphere(vit->point());
    }
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

private:
  mutable GLuint sphere_display_list;
  mutable GLUquadric* quadric;
  C2t3 c2t3_;
};

typedef Tr::Geom_traits GT;

bool insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const GT::FT size);

Scene_item* cgal_code_remesh(Polyhedron* pMesh,
                             const double angle,
                             const double sizing,
                             const double approx,
                             int tag) {
  CGAL::set_error_behaviour(CGAL::ABORT);
  if(!pMesh) return 0;

  // remesh


  Tr& triangulation = * new Tr;; // 3D-Delaunay triangulation
  C2t3& c2t3 = *(new C2t3(triangulation));
  // C2t3 c2t3(triangulation); // 2D-complex in 3D-Delaunay triangulation

  // meshing parameters
  CGAL::Surface_mesh_default_criteria_3<Tr> facets_criteria(angle,sizing,approx);

  // AABB tree
  CGAL::Timer timer;
  timer.start();
  std::cerr << "Build AABB tree...";
  typedef CGAL::Simple_cartesian<double> Simple_cartesian_kernel;
  // input surface
  typedef CGAL::AABB_polyhedral_oracle<Polyhedron,Kernel,Simple_cartesian_kernel> Input_surface;
  Input_surface input(*pMesh);
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // initial point set
  timer.reset();
  std::cerr << "Insert initial point set...";
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

  insert_spheres(c2t3, pMesh, sizing);
  std::cerr << c2t3.number_of_facets() << std::endl;
  return new Scene_c2t3_item(c2t3);
  // remesh
  timer.reset();
  std::cerr << "Remesh...";
  switch(tag) {
  case 0: 
    CGAL::make_surface_mesh(c2t3, input, input, facets_criteria, CGAL::Non_manifold_tag());
    break;
  case 1:
    CGAL::make_surface_mesh(c2t3, input, input, facets_criteria, CGAL::Manifold_tag());
    break;
  default:
    CGAL::make_surface_mesh(c2t3, input, input, facets_criteria, CGAL::Manifold_with_boundary_tag());
  }

  std::cerr << "done (" << timer.time() << " ms, " << triangulation.number_of_vertices() << " vertices)" << std::endl;

  if(triangulation.number_of_vertices() > 0)
  {
    // add remesh as new polyhedron
    Polyhedron *pRemesh = new Polyhedron;
    CGAL::Complex_2_in_triangulation_3_polyhedron_builder<C2t3, Polyhedron> builder(c2t3);
    try {
      CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
      pRemesh->delegate(builder);
    } catch(CGAL::Failure_exception)
    {
    }
    CGAL::set_error_behaviour(CGAL::ABORT);

    if(c2t3.number_of_facets() != pRemesh->size_of_facets())
    {
      delete pRemesh;
      std::stringstream temp_file;
      if(!CGAL::output_surface_facets_to_off(temp_file, c2t3))
      {
        std::cerr << "Cannot write the mesh to an off file!\n";
        return 0;
      }
      Scene_polygon_soup* soup = new Scene_polygon_soup();
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

#include "Polyhedron_demo_remeshing_plugin_cgal_code.moc"
