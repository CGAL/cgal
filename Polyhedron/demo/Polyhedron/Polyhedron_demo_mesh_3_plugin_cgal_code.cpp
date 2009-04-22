#include <CGAL/AABB_tree/Triangle_3_segment_3_intersection.h>
#include <CGAL/AABB_tree/Triangle_3_ray_3_intersection.h>

#include "Polyhedron_type.h"

#include <CGAL/AABB_tree/AABB_polyhedral_oracle.h>
#include <CGAL/AABB_tree/AABB_tree.h>

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Weighted_point_with_surface_index_geom_traits.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Volume_mesher_cell_base_3.h>


#include <CGAL/Volume_mesher_3.h>
#include <CGAL/Surface_mesher/Vertices_on_the_same_surface_criterion.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/Slivers_exuder.h>


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/Timer.h>

// Geometric traits
typedef CGAL::Regular_triangulation_filtered_traits_3<Kernel> Regular_traits;
struct My_traits : public Regular_traits
{
  typedef CGAL::Weighted_point_with_surface_index<Regular_traits::Point_3> Weighted_point_3;
  //     typedef CGAL::Point_with_psc_localisation<My_traits1::Weighted_point_3> Weighted_point_3;
  typedef Weighted_point_3 Point_3;
};

// Vertex base
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits> Vb1;
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, My_traits, Vb1> Vb;

// Cell base
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Surface_mesh_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Volume_mesher_cell_base_3<My_traits, Cb2> Cb;

typedef CGAL::Triangulation_cell_base_with_circumcenter_3<My_traits,
                                                          Cb> Cb_with_circumcenter;

// TDS
typedef CGAL::Triangulation_data_structure_3<Vb, Cb_with_circumcenter> Tds;
// Triangulation
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;

// 3D complex
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef My_traits::Weighted_point_3 Point_3;

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

class Q_DECL_EXPORT Scene_c3t3_item : public Scene_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item()
    : tr_(), c2t3_(tr_), frame(new ManipulatedFrame())

  {}

  ~Scene_c3t3_item()
  {
    delete frame;
  }

  const Tr& triangulation() const {
    return tr_;
  }

  Tr& triangulation() {
    return tr_;
  }

  const C2t3& c2t3() const {
    return c2t3_;
  }

  C2t3& c2t3() {
    return c2t3_;
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

  Plane_3 plane() const {
    const qglviewer::Vec& pos = frame->position();
    const qglviewer::Vec& n = 
      frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
    return Plane_3(n[0], n[1],  n[2], - n * pos);
  }

  bool isFinite() const { return true; }
  bool isEmpty() const {
    return triangulation().number_of_vertices() == 0;
  }

  Bbox bbox() const { 
    if(isEmpty())
      return Bbox(); 
    else {
      CGAL::Bbox_3 result = triangulation().vertices_begin()->point().bbox();
      for(Tr::Finite_vertices_iterator
            vit = ++triangulation().finite_vertices_begin(),
            end = triangulation().finite_vertices_end();
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
          cit = triangulation().finite_cells_begin(),
          end = triangulation().finite_cells_end();
        cit != end; ++cit) 
    {
      if( cit->is_in_domain() )
        ++number_of_tets;
    }
    return tr("<p><b>3D complex in a 3D triangulation</b></p>"
              "<p>Number of vertices: %1<br />"
              "Number of surface facets: %2<br />"
              "Number of volume tetrahedra: %3</p>")
      .arg(triangulation().number_of_vertices())
      .arg(c2t3().number_of_facets())
      .arg(number_of_tets);
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const { 
    return (m != Gouraud); // CHECK THIS!
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

    const My_traits::Plane_3& plane = this->plane();
    GLdouble clip_plane[4];
    clip_plane[0] = -plane.a();
    clip_plane[1] = -plane.b();
    clip_plane[2] = -plane.c();
    clip_plane[3] = -plane.d();

    ::glClipPlane(GL_CLIP_PLANE0, clip_plane);
    ::glEnable(GL_CLIP_PLANE0);
    ::glBegin(GL_TRIANGLES);
    for(C2t3::Facet_iterator 
          fit = c2t3().facets_begin(),
          end = c2t3().facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Point_3& pa = cell->vertex((index+1)&3)->point();
      const Point_3& pb = cell->vertex((index+2)&3)->point();
      const Point_3& pc = cell->vertex((index+3)&3)->point();
      typedef My_traits::Oriented_side Side;
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
    CGALglcolor(this->color().darker(150));
    for(Tr::Finite_cells_iterator
          cit = triangulation().finite_cells_begin(),
          end = triangulation().finite_cells_end();
        cit != end; ++cit) 
    {
      if(! cit->is_in_domain() )
        continue;

        const Point_3& pa = cit->vertex(0)->point();
        const Point_3& pb = cit->vertex(1)->point();
        const Point_3& pc = cit->vertex(2)->point();
        const Point_3& pd = cit->vertex(3)->point();
        typedef My_traits::Oriented_side Side;
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
//         if(c2t3().is_in_complex(cit, i)) continue;
//         const Point_3& pa = cit->vertex((i+1)&3)->point();
//         const Point_3& pb = cit->vertex((i+2)&3)->point();
//         const Point_3& pc= cit->vertex((i+3)&3)->point();
//         typedef My_traits::Oriented_side Side;
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
  static void draw_triangle(const Point_3& pa, 
                            const Point_3& pb, 
                            const Point_3& pc) {
    My_traits::Vector_3 n = cross_product(pb - pa, pc -pa);
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

  Tr tr_; // 3D-Delaunay triangulation
  C2t3 c2t3_;

  qglviewer::ManipulatedFrame* frame;
};

Scene_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing)
{
  if(!pMesh) return 0;

  // remesh

  typedef Tr::Geom_traits GT;

  Scene_c3t3_item* new_item = new Scene_c3t3_item();

  Tr& triangulation = new_item->triangulation();
  C2t3& c2t3 = new_item->c2t3();

  // meshing parameters
  typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
  typedef CGAL::Surface_mesher::Standard_criteria <Criterion> Facets_criteria;
  CGAL::Surface_mesher::Curvature_size_criterion<Tr>
    curvature_size_criterion (approx);
  CGAL::Surface_mesher::Uniform_size_criterion<Tr>
    uniform_size_criterion (sizing);
  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
    aspect_ratio_criterion (angle);
  CGAL::Surface_mesher::Vertices_on_the_same_surface_criterion<Tr>
    vertices_on_the_same_surface_criterion;

  std::vector<Criterion*> criterion_vector;
  criterion_vector.push_back(&aspect_ratio_criterion);
  criterion_vector.push_back(&uniform_size_criterion);
  criterion_vector.push_back(&curvature_size_criterion);
//   criterion_vector.push_back(&vertices_on_the_same_surface_criterion);
  Facets_criteria facets_criteria (criterion_vector);

  typedef CGAL::Mesh_criteria_3<Tr> Tets_criteria;
  Tets_criteria tets_criteria(4., tets_sizing);

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
  unsigned int nb_initial_points = 10;
  Polyhedron::Point_const_iterator it;
  typedef CGAL::Cartesian_converter<Kernel,GT> Converter;
  Converter convert;
  unsigned int i = 0;
  for(it = pMesh->points_begin();
      it != pMesh->points_end(), i < nb_initial_points;
      it++, i++)
    triangulation.insert(convert(*it));
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  std::cerr << "Number of initial points: "
            << new_item->triangulation().number_of_vertices() << "\n";

  typedef CGAL::Volume_mesher_3<
    C2t3, 
    Input_surface, 
    Facets_criteria,
    Tets_criteria,
    Input_surface
      > Mesher;

  Mesher mesher(c2t3, input, facets_criteria, tets_criteria, input);
    

  // remesh
  timer.reset();
  std::cerr << "Remesh...";
  mesher.refine_surface();
  mesher.refine_mesh();
  std::cerr << "done (" << timer.time() << " ms, " << triangulation.number_of_vertices() << " vertices)" << std::endl;

  using CGAL::output_to_medit;

  if(triangulation.number_of_vertices() > 0)
  {
    std::ofstream medit_out("out.mesh");
    CGAL::output_to_medit(medit_out, c2t3);
    const Scene_item::Bbox& bbox = new_item->bbox();
    new_item->setPosition((bbox.xmin + bbox.xmax)/2.f,
                          (bbox.ymin + bbox.ymax)/2.f,
                          (bbox.zmin + bbox.zmax)/2.f);
    return new_item;
  }
  else
    return 0;
}

#include "Scene_c3t3_item.moc"
