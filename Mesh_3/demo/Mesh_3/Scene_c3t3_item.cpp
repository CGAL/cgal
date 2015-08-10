#include "config.h"

#include "Scene_c3t3_item.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QPainter>

#include <map>
#include <vector>
#include <CGAL/gl.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

#include <CGAL_demo/Scene_interface.h>
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
struct Scene_c3t3_item_priv {
  Scene_c3t3_item_priv() : c3t3() {}
  Scene_c3t3_item_priv(const C3t3& c3t3_) : c3t3(c3t3_) {}

  C3t3 c3t3;
  QVector<QColor> colors;
};
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



/**************************************************
****************SHADER FUNCTIONS******************/

void Scene_c3t3_item::compile_shaders()
{
    //The mesh
    for(int i=0; i< vboSize; i++)
        buffers[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();

    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec3 normal;\n"
        "attribute highp vec3 inColor; \n"

        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 mv_matrix; \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "varying highp vec4 color; \n"
        "void main(void)\n"
        "{\n"
        "   color = vec4(inColor, 1.0); \n"
        "   fP = mv_matrix * vertex; \n"
        "   fN = mat3(mv_matrix)* normal; \n"
        "   gl_Position = mvp_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 330 \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "varying vec4 color; \n"
        "uniform highp vec4 light_pos;  \n"
        "uniform highp vec4 light_diff; \n"
        "uniform highp vec4 light_spec; \n"
        "uniform highp vec4 light_amb;  \n"
        "uniform float spec_power ; \n"

        "void main(void) { \n"

        "   vec3 L = light_pos.xyz - fP.xyz; \n"
        "   vec3 V = -fP.xyz; \n"
        "   vec3 N; \n"
        "   if(fN == vec3(0.0,0.0,0.0)) \n"
        "       N = vec3(0.0,0.0,0.0); \n"
        "   else \n"
        "       N = normalize(fN); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"
        "   vec3 R = reflect(-L, N); \n"
        "   vec4 diffuse = max(abs(dot(N,L)),0) * light_diff*color; \n"
        "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

         "gl_FragColor = color*light_amb + diffuse + specular; \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }

    //The grid
    //Vertex source code
    const char vertex_source_grid[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"

        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 f_matrix; \n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix * f_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source_grid[] =
    {
        "#version 120 \n"
        "uniform vec4 color; \n"
        "void main(void) { \n"
        "gl_FragColor = color; \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader_grid = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_grid ->compileSourceCode(vertex_source_grid))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader_grid = new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_grid ->compileSourceCode(fragment_source_grid ))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program_grid .addShader(vertex_shader_grid ))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_grid .addShader(fragment_shader_grid ))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_grid .link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }

}

void Scene_c3t3_item::compute_elements()
{
    v_poly.resize(0);
    normal.resize(0);
    color_triangles.resize(0);

    draw_grid((float)complex_diag(this),&v_grid);
    if(isEmpty())
      return;

    const Kernel::Plane_3& plane = this->plane();
//TRIANGLES
    for(C3t3::Facets_in_complex_iterator
          fit = c3t3().facets_in_complex_begin(),
          end = c3t3().facets_in_complex_end();
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
      using CGAL::ON_NEGATIVE_SIDE;
      const Side sa = plane.oriented_side(pa);
      const Side sb = plane.oriented_side(pb);
      const Side sc = plane.oriented_side(pc);
      if(sa == ON_NEGATIVE_SIDE &&
         sb == ON_NEGATIVE_SIDE &&
         sc == ON_NEGATIVE_SIDE)
      {
  #ifdef SHOW_REMAINING_BAD_ELEMENT_IN_RED

          Tr::Facet mirror_facet = c3t3().triangulation().mirror_facet(*fit);
          bool blueOrRed = false;
          if(cell->mark == index || mirror_facet.first->mark == mirror_facet.second)
          {
            std::cerr << "================== BAD TRIANGLE =================" << std::endl;
            blueOrRed = true;

            if(cell->mark2 != -1)
            {
              const Kernel::Point_3& pa2 = cell->vertex((cell->mark2+1)&3)->point();
              const Kernel::Point_3& pb2 = cell->vertex((cell->mark2+2)&3)->point();
              const Kernel::Point_3& pc2 = cell->vertex((cell->mark2+3)&3)->point();

              std::cerr << "================== BLUE =================" << std::endl;

              Kernel::Vector_3 n = cross_product(pb2 - pa2, pc2 -pa2);
              n = n / CGAL::sqrt(n*n);

              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              v_poly.push_back(pa2.x()); v_poly.push_back(pa2.y()); v_poly.push_back(pa2.z());
              v_poly.push_back(pb2.x()); v_poly.push_back(pb2.y()); v_poly.push_back(pb2.z());
              v_poly.push_back(pc2.x()); v_poly.push_back(pc2.y()); v_poly.push_back(pc2.z());
              color_triangles.push_back(0.0); color_triangles.push_back(0.0); color_triangles.push_back(1.0);
              color_triangles.push_back(0.0); color_triangles.push_back(0.0); color_triangles.push_back(1.0);
              color_triangles.push_back(0.0); color_triangles.push_back(0.0); color_triangles.push_back(1.0);

              const Tr::Facet f_blue(cell, cell->mark2);
              Tr::Facet mirror_f_blue = c3t3().triangulation().mirror_facet(f_blue);
              const Kernel::Point_3& dual_edge_pa = c3t3().triangulation().dual(f_blue.first);
              const Kernel::Point_3& dual_edge_pb = c3t3().triangulation().dual(mirror_f_blue.first);
              const Kernel::Point_3& dual_edge_pc = dual_edge_pa + Kernel::Vector_3(0.001, 0., 0.);
              Kernel::Vector_3 n = cross_product(dual_edge_pb - dual_edge_pa, dual_edge_pc -dual_edge_pa);
              n = n / CGAL::sqrt(n*n);

              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

              v_poly.push_back(dual_edge_pa.x()); v_poly.push_back(dual_edge_pa.y()); v_poly.push_back(dual_edge_pa.z());
              v_poly.push_back(dual_edge_pb.x()); v_poly.push_back(dual_edge_pb.y()); v_poly.push_back(dual_edge_pb.z());
              v_poly.push_back(dual_edge_pc.x()); v_poly.push_back(dual_edge_pc.y()); v_poly.push_back(dual_edge_pc.z());

              color_triangles.push_back(1.0); color_triangles.push_back(1.0); color_triangles.push_back(0.0);
              color_triangles.push_back(1.0); color_triangles.push_back(1.0); color_triangles.push_back(0.0);
              color_triangles.push_back(1.0); color_triangles.push_back(1.0); color_triangles.push_back(0.0);

            else if(mirror_facet.first->mark2 != -1)
            {
              const Kernel::Point_3& pa2 = mirror_facet.first->vertex((mirror_facet.first->mark2+1)&3)->point();
              const Kernel::Point_3& pb2 = mirror_facet.first->vertex((mirror_facet.first->mark2+2)&3)->point();
              const Kernel::Point_3& pc2 = mirror_facet.first->vertex((mirror_facet.first->mark2+3)&3)->point();


              std::cerr << "================== BLUE =================" << std::endl;
              Kernel::Vector_3 n = cross_product(pb2 - pa2, pc2 -pa2);
              n = n / CGAL::sqrt(n*n);

              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              v_poly.push_back(pa2.x()); v_poly.push_back(pa2.y()); v_poly.push_back(pa2.z());
              v_poly.push_back(pb2.x()); v_poly.push_back(pb2.y()); v_poly.push_back(pb2.z());
              v_poly.push_back(pc2.x()); v_poly.push_back(pc2.y()); v_poly.push_back(pc2.z());
              color_triangles.push_back(0.0); color_triangles.push_back(0.0); color_triangles.push_back(1.0);
              color_triangles.push_back(0.0); color_triangles.push_back(0.0); color_triangles.push_back(1.0);
              color_triangles.push_back(0.0); color_triangles.push_back(0.0); color_triangles.push_back(1.0);


              const Tr::Facet f_blue(mirror_facet.first, mirror_facet.first->mark2);
              Tr::Facet mirror_f_blue = c3t3().triangulation().mirror_facet(f_blue);
              const Kernel::Point_3& dual_edge_pa = c3t3().triangulation().dual(f_blue.first);
              const Kernel::Point_3& dual_edge_pb = c3t3().triangulation().dual(mirror_f_blue.first);
              const Kernel::Point_3& dual_edge_pc = dual_edge_pa + Kernel::Vector_3(0.001, 0., 0.);

              Kernel::Vector_3 n = cross_product(dual_edge_pb - dual_edge_pa, dual_edge_pc -dual_edge_pa);
              n = n / CGAL::sqrt(n*n);

              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
              normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

              v_poly.push_back(dual_edge_pa.x()); v_poly.push_back(dual_edge_pa.y()); v_poly.push_back(dual_edge_pa.z());
              v_poly.push_back(dual_edge_pb.x()); v_poly.push_back(dual_edge_pb.y()); v_poly.push_back(dual_edge_pb.z());
              v_poly.push_back(dual_edge_pc.x()); v_poly.push_back(dual_edge_pc.y()); v_poly.push_back(dual_edge_pc.z());

              color_triangles.push_back(1.0); color_triangles.push_back(1.0); color_triangles.push_back(0.0);
              color_triangles.push_back(1.0); color_triangles.push_back(1.0); color_triangles.push_back(0.0);
              color_triangles.push_back(1.0); color_triangles.push_back(1.0); color_triangles.push_back(0.0);
            }
          }
          else
          {
                Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
                n = n / CGAL::sqrt(n*n);

                normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
                normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
                normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

                v_poly.push_back(pa.x()); v_poly.push_back(pa.y()); v_poly.push_back(pa.z());
                v_poly.push_back(pb.x()); v_poly.push_back(pb.y()); v_poly.push_back(pb.z());
                v_poly.push_back(pc.x()); v_poly.push_back(pc.y()); v_poly.push_back(pc.z());

            if(cell->subdomain_index() == 0) {
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].blueF());

            }
            else {
              CGALglcolor(d->colors[cell->subdomain_index()]);
              color_triangles.push_back(d->colors[cell->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].blueF());
            }
          }
        }


  #else
          Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
          n = n / CGAL::sqrt(n*n);

          normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
          normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
          normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

          v_poly.push_back(pa.x()); v_poly.push_back(pa.y()); v_poly.push_back(pa.z());
          v_poly.push_back(pb.x()); v_poly.push_back(pb.y()); v_poly.push_back(pb.z());
          v_poly.push_back(pc.x()); v_poly.push_back(pc.y()); v_poly.push_back(pc.z());

          if(cell->subdomain_index() == 0) {
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->neighbor(index)->subdomain_index()].blueF());

          }
          else {
              color_triangles.push_back(d->colors[cell->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].blueF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].redF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].greenF());
              color_triangles.push_back(d->colors[cell->subdomain_index()].blueF());
          }
  #endif
      }
    }

//end TRIANGLES
    //re TRIANGLES
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
          QColor temp = d->colors[cit->subdomain_index()].darker(150);
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());

          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());

          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());

          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());
          color_triangles.push_back(temp.redF());
          color_triangles.push_back(temp.greenF());
          color_triangles.push_back(temp.blueF());


      Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
      n = n / CGAL::sqrt(n*n);

      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

      v_poly.push_back(pa.x()); v_poly.push_back(pa.y()); v_poly.push_back(pa.z());
      v_poly.push_back(pb.x()); v_poly.push_back(pb.y()); v_poly.push_back(pb.z());
      v_poly.push_back(pc.x()); v_poly.push_back(pc.y()); v_poly.push_back(pc.z());

      n = cross_product(pb - pa, pd - pa);
      n = n / CGAL::sqrt(n*n);

      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

      v_poly.push_back(pa.x()); v_poly.push_back(pa.y()); v_poly.push_back(pa.z());
      v_poly.push_back(pb.x()); v_poly.push_back(pb.y()); v_poly.push_back(pb.z());
      v_poly.push_back(pd.x()); v_poly.push_back(pd.y()); v_poly.push_back(pd.z());

      n = cross_product(pc - pa, pd - pa);
      n = n / CGAL::sqrt(n*n);

      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

      v_poly.push_back(pa.x()); v_poly.push_back(pa.y()); v_poly.push_back(pa.z());
      v_poly.push_back(pc.x()); v_poly.push_back(pc.y()); v_poly.push_back(pc.z());
      v_poly.push_back(pd.x()); v_poly.push_back(pd.y()); v_poly.push_back(pd.z());

      n = cross_product(pc - pb, pd - pb);
      n = n / CGAL::sqrt(n*n);

      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());
      normal.push_back(n.x()); normal.push_back(n.y()); normal.push_back(n.z());

      v_poly.push_back(pb.x()); v_poly.push_back(pb.y()); v_poly.push_back(pb.z());
      v_poly.push_back(pc.x()); v_poly.push_back(pc.y()); v_poly.push_back(pc.z());
      v_poly.push_back(pd.x()); v_poly.push_back(pd.y()); v_poly.push_back(pd.z());
      }
    }
//re end TRIANGLES

}

void Scene_c3t3_item::initialize_buffers() const
{
    rendering_program.bind();

        vao[0].bind();
        buffers[0].bind();
        buffers[0].allocate(v_poly.data(), static_cast<int>(v_poly.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(normal.data(), static_cast<int>(normal.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[1].release();

        buffers[2].bind();
        buffers[2].allocate(color_triangles.data(), static_cast<int>(color_triangles.size()*sizeof(float)));
        colorLocation[0] = rendering_program.attributeLocation("inColor");
        rendering_program.enableAttributeArray(colorLocation[0]);
        rendering_program.setAttributeBuffer(colorLocation[0],GL_FLOAT,0,3);
        buffers[2].release();

        vao[0].release();

    rendering_program.release();

    rendering_program_grid.bind();

        vao[1].bind();
        buffers[3].bind();
        buffers[3].allocate(v_grid.data(), static_cast<int>(v_grid.size()*sizeof(float)));
        poly_vertexLocation[1] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[3].release();
        vao[1].release();

    rendering_program_grid.release();
    are_buffers_initialized = true;

}

void Scene_c3t3_item::attrib_buffers(Viewer* viewer) const
{
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    QMatrix4x4 fMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvMatrix.data()[i] = (float)mat[i];
    }
    frame->getMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        fMatrix.data()[i] = (float)mat[i];
    }


    QVector4D	position(0.0f,0.0f,1.0f,1.0f );
    GLboolean isTwoSide;
    viewer->glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE,&isTwoSide);
    // define material
     QVector4D	ambient;
     QVector4D	diffuse;
     QVector4D	specular;
     GLfloat      shininess ;
    // Ambient
    ambient[0] = 0.29225f;
    ambient[1] = 0.29225f;
    ambient[2] = 0.29225f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.50754f;
    diffuse[1] = 0.50754f;
    diffuse[2] = 0.50754f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.0f;
    specular[1] = 0.0f;
    specular[2] = 0.0f;
    specular[3] = 0.0f;
    // Shininess
    shininess = 51.2f;


    rendering_program.bind();
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    mvLocation[0] = rendering_program.uniformLocation("mv_matrix");
    lightLocation[0] = rendering_program.uniformLocation("light_pos");
    lightLocation[1] = rendering_program.uniformLocation("light_diff");
    lightLocation[2] = rendering_program.uniformLocation("light_spec");
    lightLocation[3] = rendering_program.uniformLocation("light_amb");
    lightLocation[4] = rendering_program.uniformLocation("spec_power");

    rendering_program.release();
    rendering_program_grid.bind();

    mvpLocation[1] = rendering_program_grid.uniformLocation("mvp_matrix");
    fmatLocation = rendering_program_grid.uniformLocation("f_matrix");
    colorLocation[1] = rendering_program_grid.uniformLocation("color");
    rendering_program_grid.release();

    rendering_program.bind();
    rendering_program.setUniformValue(lightLocation[0], position);
    rendering_program.setUniformValue(twosideLocation, isTwoSide);
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);
    rendering_program.setUniformValue(mvLocation[0], mvMatrix);
    rendering_program.setUniformValue(lightLocation[1], diffuse);
    rendering_program.setUniformValue(lightLocation[2], specular);
    rendering_program.setUniformValue(lightLocation[3], ambient);
    rendering_program.setUniformValue(lightLocation[4], shininess);
    rendering_program.release();

    rendering_program_grid.bind();
    rendering_program_grid.setUniformValue(mvpLocation[1], mvpMatrix);
    rendering_program_grid.setUniformValue(fmatLocation, fMatrix);
    rendering_program_grid.release();

}

enum { DRAW = 0, DRAW_EDGES = 1 };

void
Scene_c3t3_item::draw(Viewer* viewer)const {
    if(!are_buffers_initialized)
        initialize_buffers();
     vao[0].bind();
    attrib_buffers(viewer);
    rendering_program.bind();
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(v_poly.size()/3));
    rendering_program.release();
    vao[0].release();
}

void
Scene_c3t3_item::draw_edges(Viewer* viewer) const {
    if(!are_buffers_initialized)
        initialize_buffers();
    vao[1].bind();
    attrib_buffers(viewer);
    rendering_program_grid.bind();
    QColor color;
    color.setRgbF(this->color().redF(), this->color().greenF(), this->color().blueF());
    rendering_program_grid.setUniformValue(colorLocation[1], color);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_grid.size()/3));
    rendering_program_grid.release();
    vao[1].release();
}



template<typename C3t3>
std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value);


Scene_c3t3_item::
Scene_c3t3_item()
  : d(new Scene_c3t3_item_priv())
  , frame(new ManipulatedFrame())
  , histogram_()
  , data_item_(NULL)
  , indices_()
{
  compile_shaders();
  connect(frame, SIGNAL(modified()), this, SLOT(changed()));
  c3t3_changed();
}

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3)
  : d(new Scene_c3t3_item_priv(c3t3)), frame(new ManipulatedFrame())
  , histogram_(), data_item_(NULL), indices_()
{
  compile_shaders();
  connect(frame, SIGNAL(modified()), this, SLOT(changed()));
  c3t3_changed();
}

Scene_c3t3_item::~Scene_c3t3_item()
{
  for(int i=0; i<vboSize; i++)
      buffers[i].destroy();
  for(int i=0; i<vaoSize; i++)
      vao[i].destroy();
  delete frame;
  delete d;
}

const C3t3& 
Scene_c3t3_item::c3t3() const {
  return d->c3t3;
}

C3t3& 
Scene_c3t3_item::c3t3()
{
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
  return tr("<p>3D complex in a 3D triangulation:<br />"
            "<b>%4</b></p>"
            "<p>Number of vertices: %1<br />"
            "Number of surface facets: %2<br />"
            "Number of volume tetrahedra: %3</p>")
    .arg(c3t3().triangulation().number_of_vertices())
    .arg(c3t3().number_of_facets_in_complex())
    .arg(number_of_tets)
    .arg(this->name());
}

QPixmap
Scene_c3t3_item::graphicalToolTip() const
{
  if ( ! histogram_.isNull() )
  {
    return histogram_;
  }
  else
  {
    const_cast<Scene_c3t3_item&>(*this).build_histogram();
    return histogram_;
  }
}
 
void
Scene_c3t3_item::build_histogram()
{
#ifdef CGAL_MESH_3_DEMO_BIGGER_HISTOGRAM_WITH_WHITE_BACKGROUNG
  // Create an histogram_ and display it
  const int height = 280;
  const int top_margin = 5;
  const int left_margin = 20;
  const int drawing_height = height-top_margin*2;
  const int width = 804;
  const int cell_width = 4;
  const int text_margin = 3;
  const int text_height = 34;
  
  histogram_ = QPixmap(width,height+text_height);
  histogram_.fill(QColor(255,255,255));
#else
  // Create an histogram_ and display it
  const int height = 140;
  const int top_margin = 5;
  const int left_margin = 20;
  const int drawing_height = height-top_margin*2;
  const int width = 402;
  const int cell_width = 2;
  const int text_margin = 3;
  const int text_height = 20;
  
  histogram_ = QPixmap(width,height+text_height);
  histogram_.fill(QColor(192,192,192));
#endif  

  QPainter painter(&histogram_);
  painter.setPen(Qt::black);
  painter.setBrush(QColor(128,128,128));
  //painter.setFont(QFont("Arial", 30));
  
  // Build histogram_ data
  double min_value, max_value;
  std::vector<int> histo_data = create_histogram(c3t3(),min_value,max_value);
  
  // Get maximum value (to normalize)
  int max_size = 0;
  for ( std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end() ;
       it != end ; ++it )
  {
    max_size = (std::max)(max_size,*it);
  }
  
  // colored histogram
  int j = 0;
  
  // draw
  int i=left_margin;
  for ( std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end() ;
       it != end ; ++it, i+=cell_width )
  {
    int line_height = static_cast<int>( std::ceil(static_cast<double>(drawing_height) *
      static_cast<double>(*it)/static_cast<double>(max_size)) + .5);
    
    painter.fillRect(i,
                     drawing_height+top_margin-line_height,
                     cell_width,
                     line_height,
                     get_histogram_color(j++));
  }
  
  // draw bottom horizontal line
  painter.setPen(Qt::blue);
  
  painter.drawLine(QPoint(left_margin, drawing_height + top_margin),
                   QPoint(left_margin + static_cast<int>(histo_data.size())*cell_width, 
                          drawing_height + top_margin));

  
  // draw min value and max value
  const int min_tr_width = static_cast<int>( 2*(std::floor(min_value)*cell_width + left_margin) );
  const int max_tr_width = static_cast<int>( 
    2*((histo_data.size()-std::floor(max_value))*cell_width + left_margin) );
  const int tr_y = drawing_height + top_margin + text_margin;
  
  painter.setPen(get_histogram_color(min_value));
  QRect min_text_rect (0, tr_y, min_tr_width, text_height);
  painter.drawText(min_text_rect, Qt::AlignCenter, tr("%1").arg(min_value,0,'f',1));
  
  painter.setPen(get_histogram_color(max_value));           
  QRect max_text_rect (width - max_tr_width, tr_y, max_tr_width, text_height);
  painter.drawText(max_text_rect, Qt::AlignCenter, tr("%1").arg(max_value,0,'f',1));
}

template<typename C3t3>
std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value)
{
  typedef typename C3t3::Triangulation::Point Point_3;
  
	std::vector<int> histo(181,0);
  
  min_value = 180.;
  max_value = 0.;
  
	for (typename C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit)
	{
		if( !c3t3.is_in_complex(cit))
			continue;
		
#ifdef CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM
    if (c3t3.in_dimension(cit->vertex(0)) <= 1
     || c3t3.in_dimension(cit->vertex(1)) <= 1
     || c3t3.in_dimension(cit->vertex(2)) <= 1
     || c3t3.in_dimension(cit->vertex(3)) <= 1)
      continue;
#endif //CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM
		
		const Point_3& p0 = cit->vertex(0)->point();
		const Point_3& p1 = cit->vertex(1)->point();
		const Point_3& p2 = cit->vertex(2)->point();
		const Point_3& p3 = cit->vertex(3)->point();
		
		double a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0,p1,p2,p3)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p2, p1, p3)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p3, p1, p2)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p2, p0, p3)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p3, p0, p2)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p2, p3, p0, p1)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
	}
  
	return histo;	
}


QColor
Scene_c3t3_item::get_histogram_color(const double v) const
{
  if ( v < 5 )            { return Qt::red; }
  else if ( v < 10 )      { return QColor(215,108,0); }
  else if ( v < 15 )      { return QColor(138,139,0); }
  else if ( v < 165 )     { return QColor(60,136,64); }
  else if ( v < 170 )     { return QColor(138,139,1); }
  else if ( v < 175 )     { return QColor(215,108,0); }
  else /* 175<v<=180 */   { return Qt::red; }
}


void
Scene_c3t3_item::setColor(QColor c)
{
  color_ = c;
  compute_color_map(c);
}

void
Scene_c3t3_item::c3t3_changed()
{
  // Update colors
  // Fill indices map and get max subdomain value
  indices_.clear();
  
  int max = 0;
  for(C3t3::Cells_in_complex_iterator cit = this->c3t3().cells_in_complex_begin(),
      end = this->c3t3().cells_in_complex_end() ; cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    indices_.insert(cit->subdomain_index());
  }

  d->colors.resize(max+1);
  compute_color_map(color_);
  
  // Rebuild histogram
  build_histogram();
  compute_elements();
  are_buffers_initialized = false;
}

void
Scene_c3t3_item::compute_color_map(const QColor& c)
{
  typedef Indices::size_type size_type;

  size_type nb_domains = indices_.size();
  size_type i = 0;
  for(Indices::iterator it = indices_.begin(), end = indices_.end();
      it != end; ++it, ++i)
  {
    double hue = c.hueF() + 1./nb_domains * i;
    if ( hue > 1 ) { hue -= 1.; }
    d->colors[*it] = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
  }
}

void Scene_c3t3_item::draw_grid(float diag, std::vector<float> *positions_grid)
{
    positions_grid->resize(0);
           float x = (2*(float)diag)/10.0;
           float y = (2*(float)diag)/10.0;
           for(int u = 0; u < 11; u++)
           {

               positions_grid->push_back(-(float)diag + x* u);
               positions_grid->push_back(-(float)diag);
               positions_grid->push_back(0.0);

               positions_grid->push_back(-(float)diag + x* u);
               positions_grid->push_back((float)diag);
               positions_grid->push_back(0.0);
           }
           for(int v=0; v<11; v++)
           {

               positions_grid->push_back(-(float)diag);
               positions_grid->push_back(-(float)diag + v * y);
               positions_grid->push_back(0.0);

               positions_grid->push_back((float)diag);
               positions_grid->push_back(-(float)diag + v * y);
               positions_grid->push_back(0.0);
           }
       }

void Scene_c3t3_item::contextual_changed()
{
    if(frame->isInMouseGrabberPool())
        c3t3_changed();
}

