#include "Scene_c3t3_item.h"

#include <QVector>
#include <QColor>
#include <map>
#include <CGAL/gl.h>

#include "Scene_item_with_display_list.h"
#include "Scene_interface.h"
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

void draw_triangle(const Kernel::Point_3& pa,
                          const Kernel::Point_3& pb,
                          const Kernel::Point_3& pc) {
  Kernel::Vector_3 n = cross_product(pb - pa, pc -pa);
  n = n / CGAL::sqrt(n*n);

  ::glNormal3d(n.x(),n.y(),n.z());

  ::glVertex3d(pa.x(),pa.y(),pa.z());
  ::glVertex3d(pb.x(),pb.y(),pb.z());
  ::glVertex3d(pc.x(),pc.y(),pc.z());
}

double complex_diag(const Scene_item* item) {
  const Scene_item::Bbox& bbox = item->bbox();
  const double& xdelta = bbox.xmax-bbox.xmin;
  const double& ydelta = bbox.ymax-bbox.ymin;
  const double& zdelta = bbox.zmax-bbox.zmin;
  const double diag = std::sqrt(xdelta*xdelta +
                                ydelta*ydelta +
                                zdelta*zdelta);
  return diag * 0.7;
}


struct Scene_c3t3_item_priv {
  Scene_c3t3_item_priv(const C3t3& c3t3_) : c3t3(c3t3_) {}

  C3t3 c3t3;
  QVector<QColor> colors;
};

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3)
  : d(new Scene_c3t3_item_priv(c3t3)), frame(new ManipulatedFrame())
{
  connect(frame, SIGNAL(modified()),
          this, SLOT(changed()));
  typedef std::set<int> Indices;
  typedef Indices::size_type size_type;
  Indices indices;
  int max = 0;
  for(Tr::Finite_cells_iterator
        cit = this->c3t3().triangulation().finite_cells_begin(),
        end = this->c3t3().triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    indices.insert(cit->subdomain_index());
  }
  d->colors.resize(max+1);
  size_type nb_domains = indices.size();
  size_type i = 0;
  for(Indices::iterator
        it = indices.begin(),
        end = indices.end();
      it != end; ++it, ++i) 
  {
    d->colors[*it] = QColor::fromHsvF( 1. / nb_domains * i, 1., 0.8);
  }
}

Scene_c3t3_item::~Scene_c3t3_item()
{
  delete frame;
}

const C3t3& 
Scene_c3t3_item::c3t3() const {
  return d->c3t3;
}

Kernel::Plane_3 
Scene_c3t3_item::plane() const {
  const qglviewer::Vec& pos = frame->position();
  const qglviewer::Vec& n =
    frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
}

Scene_item::Bbox 
Scene_c3t3_item::bbox() const {
  if(isEmpty())
    return Bbox();
  else {
    CGAL::Bbox_3 result = c3t3().triangulation().vertices_begin()->point().bbox();
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

QString 
Scene_c3t3_item::toolTip() const {
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

void
Scene_c3t3_item::direct_draw() const {
  ::glPushMatrix();
  ::glMultMatrixd(frame->matrix());
  QGLViewer::drawGrid((float)complex_diag(this));
  ::glPopMatrix();

  if(isEmpty())
    return;

  std::cerr << "Direct_draw\n";
  GLboolean lighting = ::glIsEnabled(GL_LIGHTING);
  GLboolean two_side;
  ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &two_side);
  // ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  if(!lighting)
    ::glDisable(GL_LIGHTING);

  const Kernel::Plane_3& plane = this->plane();
  GLdouble clip_plane[4];
  clip_plane[0] = -plane.a();
  clip_plane[1] = -plane.b();
  clip_plane[2] = -plane.c();
  clip_plane[3] = -plane.d();

  GLint i;
  ::glGetIntegerv(GL_POLYGON_MODE, &i);
  std::cerr << "Polygon mode: " << i << std::endl;
  std::cerr << "Lighting: " 
            << std::boolalpha << (bool)::glIsEnabled(GL_LIGHTING) << std::endl;
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
    if(cell->subdomain_index() != 0 &&
       cell->neighbor(index)->subdomain_index() != 0)
    {
      continue;
    }

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
      if(cell->subdomain_index() == 0) {
        CGALglcolor(d->colors[cell->neighbor(index)->subdomain_index()]);
      }
      else {
        CGALglcolor(d->colors[cell->subdomain_index()]);
      }
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
      false && cit != end; ++cit)
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
      CGALglcolor(d->colors[cit->subdomain_index()].darker(150));
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
  if(lighting)
    ::glEnable(GL_LIGHTING);
  else
    ::glDisable(GL_LIGHTING);
};

#include "Scene_c3t3_item.moc"
