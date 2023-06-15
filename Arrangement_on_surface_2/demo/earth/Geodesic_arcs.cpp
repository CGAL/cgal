
#include "Geodesic_arcs.h"

#include <qvector3d.h>

#include <iostream>
#include <iterator>
#include <vector>
using namespace std;


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include "arr_print.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel         Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>         Geom_traits;
typedef Geom_traits::Point_2                                      Point;
typedef Geom_traits::Curve_2                                      Curve;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits>        Topol_traits;
typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits> Arrangement;


typedef Kernel::Direction_3	Dir3;
ostream& operator << (ostream& os, const Dir3& d)
{
  os << d.dx() << ", " << d.dy() << ", " << d.dz();
  return os;
}

typedef Geom_traits::Approximate_point_2	Approximate_point_2;
ostream& operator << (ostream& os, const Approximate_point_2& d)
{
  os << d.dx() << ", " << d.dy() << ", " << d.dz();
  return os;
}

#include <CGAL/Vector_3.h>
typedef Geom_traits::Approximate_number_type	Approximate_number_type;
typedef Geom_traits::Approximate_kernel			Approximate_kernel;
typedef CGAL::Vector_3<Approximate_kernel>		Approximate_Vector_3;
typedef Approximate_kernel::Direction_3         Approximate_Direction_3;

typedef Kernel::Direction_3                Direction_3;


ostream& operator << (ostream& os, const Approximate_Vector_3& v)
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

  
  vector<Curve>  xcvs;
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0), Dir3(0, 0, -1)));
  //xcvs.push_back(ctr_cv(Dir3(0, 0, -1)));

  auto approx = traits.approximate_2_object();

  const double error = 0.001;
  std::vector<Approximate_point_2> v;
    
  const auto& xcv = xcvs[1];
  //for (const auto& xcv : xcvs)
  {
    auto oi2 = approx(xcv, error, std::back_insert_iterator(v));

    for (auto it = v.begin(); it != v.end(); ++it)
      cout << *it << endl;
    cout << "num points output = " << v.size() << endl;
  }
 

  std::vector<QVector3D> vertex_data;
  for (const auto& p : v)
  {
    const QVector3D arc_point(p.dx(), p.dy(), p.dz());
    vertex_data.push_back(arc_point);
  }
  m_num_arc_points = v.size(); // CAREFUL: not size of vertex_data!!!


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
    glDrawArrays(GL_LINE_STRIP, 0, m_num_arc_points);
  }
  glBindVertexArray(0);
}