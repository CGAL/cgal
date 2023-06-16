
#include "Geodesic_arcs.h"

#include <iostream>
#include <iterator>
#include <vector>

#include <qvector3d.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Vector_3.h>

#include "arr_print.h"


using Kernel        = CGAL::Exact_predicates_exact_constructions_kernel;
using Geom_traits   = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Point         = Geom_traits::Point_2;
using Curve         = Geom_traits::Curve_2;
using Topol_traits  = CGAL::Arr_spherical_topology_traits_2<Geom_traits>;
using Arrangement   = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;


using Dir3 =  Kernel::Direction_3	;
std::ostream& operator << (std::ostream& os, const Dir3& d)
{
  os << d.dx() << ", " << d.dy() << ", " << d.dz();
  return os;
}

using Approximate_point_2 = Geom_traits::Approximate_point_2;
std::ostream& operator << (std::ostream& os, const Approximate_point_2& d)
{
  os << d.dx() << ", " << d.dy() << ", " << d.dz();
  return os;
}

using Approximate_number_type = Geom_traits::Approximate_number_type;
using Approximate_kernel      = Geom_traits::Approximate_kernel;
using Approximate_Vector_3    = CGAL::Vector_3<Approximate_kernel>;
using Approximate_Direction_3 = Approximate_kernel::Direction_3;
using Direction_3             = Kernel::Direction_3;


std::ostream& operator << (std::ostream& os, const Approximate_Vector_3& v)
{
  os << v.x() << ", " << v.y() << ", " << v.z();
  //os << v.hx() << ", " << v.hy() << ", " << v.hz() << ", " << v.hw();
  return os;
}


Geodesic_arcs::Geodesic_arcs()
{
  initializeOpenGLFunctions();

  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Arrangement arr(&traits);


  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  
  std::vector<Curve>  xcvs;
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 0, 1)));
  xcvs.push_back(ctr_cv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  //xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0), Dir3(0, 0, -1)));
  //xcvs.push_back(ctr_cv(Dir3(0, 0, -1)));

  auto approx = traits.approximate_2_object();

  const double error = 0.001;
  
  std::vector<QVector3D> vertex_data;

  m_arc_offsets.clear();
  m_arc_offsets.push_back(0);
  for (const auto& xcv : xcvs)
  {
    std::vector<Approximate_point_2> v;
    auto oi2 = approx(xcv, error, std::back_insert_iterator(v));

    for (const auto& p : v)
    {
      const QVector3D arc_point(p.dx(), p.dy(), p.dz());
      vertex_data.push_back(arc_point);
    }
    const auto current_vertex_data_size = vertex_data.size();
    m_arc_offsets.push_back(current_vertex_data_size);
  }
  //std::cout << "offset count = " << m_arc_offsets.size() << std::endl;


  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);


  // Vertex Buffer
  glGenBuffers(1, &m_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  auto vertex_buffer_size = sizeof(QVector3D) * vertex_data.size();
  auto vertex_buffer_data = reinterpret_cast<const void*>(vertex_data.data());
  glBufferData(GL_ARRAY_BUFFER,
    vertex_buffer_size,
    vertex_buffer_data,
    GL_STATIC_DRAW);

  // Position Vertex-Attribute
  GLint position_attrib_index = 0;
  const void* position_offset = 0;
  GLsizei stride = 0;
  glVertexAttribPointer(position_attrib_index,
                        3,
                        GL_FLOAT, GL_FALSE,
                        stride,
                        position_offset);
  glEnableVertexAttribArray(position_attrib_index);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void Geodesic_arcs::draw()
{
  glBindVertexArray(m_vao);
  {
    for (int i = 1; i < m_arc_offsets.size(); i++)
    {
      const auto first = m_arc_offsets[i - 1];
      const auto count = m_arc_offsets[i] - first;
      glDrawArrays(GL_LINE_STRIP, first, count);
    }
  }
  glBindVertexArray(0);
}
