#include <CGAL/AABB_intersections.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <fstream>

#include <CGAL/Timer.h>


// @TODO: Is that the right kernel?!
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Point Point_3;

#include "Scene_item.h"
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

class Q_DECL_EXPORT Scene_c3t3_item : public Scene_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item(const C3t3& c3t3)
    : c3t3_(c3t3), frame(new ManipulatedFrame())

  {}

  ~Scene_c3t3_item()
  {
    delete frame;
  }

  const C3t3& c3t3() const {
    return c3t3_;
  }

  bool manipulatable() const {
    return true;
  }
  ManipulatedFrame* manipulatedFrame() {
    return frame;
  }

  void setPosition(float x, float y, float z) {
    frame->setPosition(x, y, z);
  }

  void setNormal(float x, float y, float z) {
    frame->setOrientation(x, y, z, 0.f);
  }

  Kernel::Plane_3 plane() const {
    const qglviewer::Vec& pos = frame->position();
    const qglviewer::Vec& n =
      frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
    return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
  }

  bool isFinite() const { return true; }
  bool isEmpty() const {
    return c3t3().triangulation().number_of_vertices() == 0;
  }

  Bbox bbox() const {
    if(isEmpty())
      return Bbox();
    else {
      CGAL::Bbox_3 result = c3t3().triangulation().finite_vertices_begin()->point().bbox();
      for(Tr::Finite_vertices_iterator
            vit = ++c3t3().triangulation().finite_vertices_begin(),
            end = c3t3().triangulation().finite_vertices_end();
          vit != end; ++vit)
      {
        result = result + vit->point().bbox();
      }
      return Bbox(result.xmin(), result.ymin(), result.zmin(),
                  result.xmax(), result.ymax(), result.zmax());
    }
  }

  Scene_c3t3_item* clone() const {
    return 0;
  }

  QString toolTip() const {
    int number_of_tets = 0;
    for(Tr::Finite_cells_iterator
          cit = c3t3().triangulation().finite_cells_begin(),
          end = c3t3().triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      if( c3t3().is_in_complex(cit) )
        ++number_of_tets;
    }
    return tr("<p><b>3D complex in a 3D triangulation</b></p>"
              "<p>Number of vertices: %1<br />"
              "Number of surface facets: %2<br />"
              "Number of volume tetrahedra: %3</p>")
      .arg(c3t3().triangulation().number_of_vertices())
      .arg(c3t3().number_of_facets())
      .arg(number_of_tets);
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud && m!=PointsPlusNormals); // CHECK THIS!
  }

  void draw() const {
    ::glPushMatrix();
    ::glMultMatrixd(frame->matrix());
    QGLViewer::drawGrid((float)complex_diag());
    ::glPopMatrix();

    if(isEmpty())
      return;

    GLboolean two_side;
    ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &two_side);
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    const Kernel::Plane_3& plane = this->plane();
    GLdouble clip_plane[4];
    clip_plane[0] = -plane.a();
    clip_plane[1] = -plane.b();
    clip_plane[2] = -plane.c();
    clip_plane[3] = -plane.d();

    ::glClipPlane(GL_CLIP_PLANE0, clip_plane);
    ::glEnable(GL_CLIP_PLANE0);
    ::glBegin(GL_TRIANGLES);
    for(C3t3::Facet_iterator
          fit = c3t3().facets_begin(),
          end = c3t3().facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Kernel::Point_3& pa = cell->vertex((index+1)&3)->point();
      const Kernel::Point_3& pb = cell->vertex((index+2)&3)->point();
      const Kernel::Point_3& pc = cell->vertex((index+3)&3)->point();
      typedef Kernel::Oriented_side Side;
      using CGAL::ON_ORIENTED_BOUNDARY;
      const Side sa = plane.oriented_side(pa);
      const Side sb = plane.oriented_side(pb);
      const Side sc = plane.oriented_side(pc);
      if( sa != ON_ORIENTED_BOUNDARY &&
          sb != ON_ORIENTED_BOUNDARY &&
          sc != ON_ORIENTED_BOUNDARY &&
          sb == sa && sc == sa )
      {
        draw_triangle(pa, pb, pc);
      }
    }
    ::glEnd();
    ::glDisable(GL_CLIP_PLANE0);

    ::glBegin(GL_TRIANGLES);
// workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
    CGALglcolor(this->color().darker(150));
#undef darker
    for(Tr::Finite_cells_iterator
          cit = c3t3().triangulation().finite_cells_begin(),
          end = c3t3().triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      if(! c3t3().is_in_complex(cit) )
        continue;

        const Kernel::Point_3& pa = cit->vertex(0)->point();
        const Kernel::Point_3& pb = cit->vertex(1)->point();
        const Kernel::Point_3& pc = cit->vertex(2)->point();
        const Kernel::Point_3& pd = cit->vertex(3)->point();
        typedef Kernel::Oriented_side Side;
        using CGAL::ON_ORIENTED_BOUNDARY;
        const Side sa = plane.oriented_side(pa);
        const Side sb = plane.oriented_side(pb);
        const Side sc = plane.oriented_side(pc);
        const Side sd = plane.oriented_side(pd);

        if( sa == ON_ORIENTED_BOUNDARY ||
            sb == ON_ORIENTED_BOUNDARY ||
            sc == ON_ORIENTED_BOUNDARY ||
            sd == ON_ORIENTED_BOUNDARY ||
            sb != sa || sc != sa || sd != sa)
        {
          draw_triangle(pa, pb, pc);
          draw_triangle(pa, pb, pd);
          draw_triangle(pa, pc, pd);
          draw_triangle(pb, pc, pd);
        }

//       for(int i = 0; i < 4; ++i) {
//         if(c3t3().is_in_complex(cit, i)) continue;
//         const Point_3& pa = cit->vertex((i+1)&3)->point();
//         const Point_3& pb = cit->vertex((i+2)&3)->point();
//         const Point_3& pc= cit->vertex((i+3)&3)->point();
//         typedef Kernel::Oriented_side Side;
//         using CGAL::ON_ORIENTED_BOUNDARY;
//         const Side sa = plane.oriented_side(pa);
//         const Side sb = plane.oriented_side(pb);
//         const Side sc = plane.oriented_side(pc);

//         if( sa == ON_ORIENTED_BOUNDARY ||
//             sb == ON_ORIENTED_BOUNDARY ||
//             sc == ON_ORIENTED_BOUNDARY ||
//             sb != sa || sc != sa )
//         {
//           draw_triangle(pa, pb, pc);
//         }
//       }
    }
    ::glEnd();
    if(!two_side)
      ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  };

private:
  static void draw_triangle(const Kernel::Point_3& pa,
                            const Kernel::Point_3& pb,
                            const Kernel::Point_3& pc) {
    Kernel::Vector_3 n = cross_product(pb - pa, pc -pa);
    n = n / CGAL::sqrt(n*n);

    ::glNormal3d(n.x(),n.y(),n.z());

    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glVertex3d(pc.x(),pc.y(),pc.z());
  }

  double complex_diag() const {
    const Bbox& bbox = this->bbox();
    const double& xdelta = bbox.xmax-bbox.xmin;
    const double& ydelta = bbox.ymax-bbox.ymin;
    const double& zdelta = bbox.zmax-bbox.zmin;
    const double diag = std::sqrt(xdelta*xdelta +
                                  ydelta*ydelta +
                                  zdelta*zdelta);
    return diag * 0.7;
  }

  C3t3 c3t3_;

  qglviewer::ManipulatedFrame* frame;
};

Scene_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                             const QString filename,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing)
{
  if(!pMesh) return 0;

  // remesh

  // Set mesh criteria
  Facet_criteria facet_criteria(angle, sizing, approx); // angle, size, approximation
  Cell_criteria cell_criteria(4, tets_sizing); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  CGAL::Timer timer;
  timer.start();
  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << angle << std::endl
            << "  facets size bound: " << sizing << std::endl
            << "  approximation bound: " << approx << std::endl
            << "  tetrahedra size bound: " << tets_sizing << std::endl;
  std::cerr << "Build AABB tree...";
  // Create domain
  Mesh_domain domain(*pMesh);
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // Meshing
  std::cerr << "Mesh...";
  Scene_c3t3_item* new_item = 
    new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(domain, criteria));

  std::cerr << "done (" << timer.time() << " ms, " << new_item->c3t3().triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(new_item->c3t3().triangulation().number_of_vertices() > 0)
  {
    std::ofstream medit_out("out.mesh");
    new_item->c3t3().output_to_medit(medit_out);
    const Scene_item::Bbox& bbox = new_item->bbox();
    new_item->setPosition((float)(bbox.xmin + bbox.xmax)/2.f,
                          (float)(bbox.ymin + bbox.ymax)/2.f,
                          (float)(bbox.zmin + bbox.zmax)/2.f);
    return new_item;
  }
  else {
    delete new_item;
    return 0;
  }
}

#include "Scene_c3t3_item.moc"
