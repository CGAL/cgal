#ifndef _VOLUME_H
#define _VOLUME_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>

#include "surface.h"
#include "binary_image.h"

#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/format.hpp>

#include <queue>
#include <vector>
#include <iterator> // std::back_inserter

#include <QString>

// kernel
// #include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Segment_3 Segment_3;

typedef CBinary_image_3<FT,Point> Binary_image;

#include <mc/MarchingCubes.h>

// surface mesher
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>
typedef CGAL::Surface_mesh_vertex_base_3<Kernel> Vb;
typedef CGAL::Surface_mesh_cell_base_3<Kernel> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Binary_image> Surface_3;
typedef CGAL::Surface_mesh_traits_generator_3<Surface_3>::Type Oracle;
typedef std::pair<Triangle_3,Vector> Facet;

class QMainWindow;
class QDoubleSpinBox;

class Volume : public Surface
{
  Q_OBJECT
public:
  Volume(QObject* parent);
  ~Volume() {}

private:
  Binary_image m_image;

  // options
  FT m_sm_angle;
  FT m_sm_radius;
  FT m_sm_distance;
  double m_isovalue;
  double m_relative_precision;

  // visualization
  bool m_view_surface;
  bool m_view_mc;

  std::vector<Facet> m_surface;
  std::vector<Facet> m_surface_mc;
  MarchingCubes mc ;

  QMainWindow* parent;
  bool m_inverse_normals;
  bool two_sides;
  bool draw_triangles_edges;
  bool use_gouraud;

  QDoubleSpinBox* spinBox_isovalue;
  QDoubleSpinBox* spinBox_radius_bound;
  QDoubleSpinBox* spinBox_distance_bound;
private:
  template <typename Iterator>
  void gl_draw_surface(Iterator begin, Iterator end);

  template <typename PointsOutputIterator>
  void search_for_connected_components(PointsOutputIterator);

public:
  void gl_draw_surface();
  void gl_draw_surface_mc();
  void gl_draw_marchingcube();
private:
  void gl_draw_one_marching_cube_vertex(int);

signals:
void new_bounding_box(double, double, double, double, double, double);

public slots:
  void set_inverse_normals(const bool);
  void set_two_sides(const bool);
  void set_draw_triangles_edges(const bool);
  void set_use_gouraud(const bool);
  void open(const QString& filename);
  void draw();
  void get_bbox(float& /*xmin*/, float& /*ymin*/, float& /*zmin*/,
		float& /*xmax*/, float& /*ymax*/, float& /*zmax*/) {}
  void close() {}
  void display_marchin_cube();
  void display_surface_mesher_result();
  void set_isovalue(double);
  void set_radius_bound(double);
  void set_distance_bound(double);
private:
  void status_message(QString);
  void busy() const;
  void not_busy() const;
  void changed_parameters();
};

template <typename PointsOutputIterator>
void Volume::search_for_connected_components(PointsOutputIterator it)
{
  const unsigned int nx = m_image.xdim();
  const unsigned int ny = m_image.ydim();
  const unsigned int nz = m_image.zdim();

  typedef unsigned char Marker;

  static const Marker outside_mark 
    = ( std::numeric_limits<Marker>::is_bounded ? 
        std::numeric_limits<Marker>::max() :
        std::numeric_limits<Marker>::infinity() );

  boost::multi_array<Marker, 3> visited(boost::extents[nx][ny][nz]);
  typedef boost::tuple<int, int, int> Indices;
  typedef std::set<Indices> Zone;
  typedef std::queue<Indices> Indices_queue;

  int number_of_connected_components = 0;
  for(unsigned int i=0;i<nx;i++)
    for(unsigned int j=0;j<ny;j++)
      for(unsigned int k=0;k<nz;k++)
      {
//         Zone zone;
        if(visited[i][j][k]>0)
          continue;
        if(m_image.value(i, j, k) <= m_isovalue) {
          visited[i][j][k] = outside_mark;
          continue;
        }

        // if we reach here, (i, j, k) is a new connected component
        ++number_of_connected_components;
        std::cerr << boost::format("Found new connected component (#%5%) "
                                   "at voxel (%1%, %2%, %3%), value=%4%\n")
          % i % j % k % m_image.value(i, j, k) % number_of_connected_components;

        int nb_voxels = 0;

        Indices_queue queue;
        Indices indices(i, j ,k);
        queue.push(indices);
//         zone.insert(indices);

        bool seed_found = false;

        Indices bbox_min = indices;
        Indices bbox_max = indices;

        while(!queue.empty()) // walk through the connected component
        {
          Indices indices = queue.front();
          queue.pop();

          // warning: local indices i, j and k.
          const int i = boost::get<0>(indices);
          const int j = boost::get<1>(indices);
          const int k = boost::get<2>(indices);

          if(visited[i][j][k]>0)
            continue;
          visited[i][j][k] = number_of_connected_components;
          ++nb_voxels;

          boost::get<0>(bbox_min) = std::min(i, boost::get<0>(bbox_min));
          boost::get<0>(bbox_max) = std::max(i, boost::get<0>(bbox_max));
          boost::get<1>(bbox_min) = std::min(j, boost::get<1>(bbox_min));
          boost::get<1>(bbox_max) = std::max(j, boost::get<1>(bbox_max));
          boost::get<2>(bbox_min) = std::min(k, boost::get<2>(bbox_min));
          boost::get<2>(bbox_max) = std::max(k, boost::get<2>(bbox_max));

          int nb_neighbors = 0;

          static const int neighbors_offset[6][3] = { { +1,  0,  0 },
                                                      { -1,  0,  0 },
                                                      {  0, +1,  0 },
                                                      {  0, -1,  0 },
                                                      {  0,  0, +1 },
                                                      {  0,  0, -1 } };
          // Visit neighbors.
          // (i_n, j_n, k_n) are indices of neighbors.
          for(int n = 0; n < 6; ++n)
          {
            const int i_n = i + neighbors_offset[n][0];
            const int j_n = j + neighbors_offset[n][1];
            const int k_n = k + neighbors_offset[n][2];
            if(i_n < 0 || i_n >= static_cast<int>(nx)) {
              ++nb_neighbors; // fake neighbor
              continue;
            }
            if(j_n < 0 || j_n >= static_cast<int>(ny)) {
              ++nb_neighbors; // fake neighbor
              continue;
            }
            if(k_n < 0 || k_n >= static_cast<int>(nz)) {
              ++nb_neighbors; // fake neighbor
              continue;
            }
            if(m_image.value(i_n, j_n, k_n) <= m_isovalue)
              visited[i_n][j_n][k_n] = outside_mark;
            else 
            {
              ++nb_neighbors;
              if(!visited[i_n][j_n][k_n]) {
                Indices indices(i_n, j_n, k_n);
                queue.push(indices);
              }
            }
          } // end for neighbors

          if(!seed_found && nb_neighbors == 6)
          {
            *it++ = m_image.point(i, j, k);
            std::cerr << boost::format("Found seed %5%, which is voxel (%1%, %2%, %3%), value=%4%\n")
              % i % j % k %  m_image.value(i, j, k) % m_image.point(i, j, k);
            seed_found = true;
          }
        } // end while !queue.empty() (with local indices i, j, k)

        // if no seed has been found, take the first voxel of the connected
        // component
        if(!seed_found) {
          *it++ = m_image.point(i, j, k);
          std::cerr << boost::format("No seed found. Choose %1%\n") %  m_image.point(i, j, k);
        }
        std::cerr << boost::format("There was %1% voxels in that component.\n"
                                   "The bounding box is (%2% %3% %4%, %5% %6% %7%).\n")
          % nb_voxels
          % boost::get<0>(bbox_min) % boost::get<1>(bbox_min) % boost::get<2>(bbox_min)
          % boost::get<0>(bbox_max) % boost::get<1>(bbox_max) % boost::get<2>(bbox_max);
      } // end for i,j,k
} // end function Volume::search_for_connected_components()

template <typename Iterator>
void Volume::gl_draw_surface(Iterator begin, Iterator end)
{
  ::glBegin(GL_TRIANGLES);
  unsigned int counter = 0;
  for(Iterator it = begin; it != end; ++it)
  {
    const Facet& f = *it;

    const Vector& n = f.second;

    if(m_inverse_normals)
      ::glNormal3d(-n.x(),-n.y(),-n.z());
    else
      ::glNormal3d(n.x(),n.y(),n.z());

    const Triangle_3& t = f.first;
    const Point& a = t[0];
    const Point& b = t[1];
    const Point& c = t[2];
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());
    ++counter;
  }
  ::glEnd();
}

#endif // _VOLUME_H
