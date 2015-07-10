#ifndef _VOLUME_H
#define _VOLUME_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/Timer.h>


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
#include <QFileInfo>

#ifdef CGAL_USE_VTK
class vtkImageReader;
class vtkImageData;
class vtkDICOMImageReader;
class vtkDemandDrivenPipeline;
class vtkImageGaussianSmooth;
#endif // CGAL_USE_VTK

class QTreeWidgetItem;

// kernel
// #include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel1;

#include <CGAL/Point_with_psc_localisation.h>
struct Kernel : public Kernel1 {
  typedef CGAL::Point_with_psc_localisation<Kernel::Point_3,
                                            const QTreeWidgetItem*> Point_3;
};

typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Segment_3 Segment_3;

// typedef CGAL::Triple<Triangle_3,Vector,const QTreeWidgetItem*> Facet;

typedef boost::tuple<Triangle_3,Vector,const QTreeWidgetItem*> Facet;

typedef CBinary_image_3<FT,Point> Binary_image;

class QTreeWidgetItem;

// surface mesher
// #define CGAL_MESHES_NO_OUTPUT
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>
typedef CGAL::Surface_mesh_vertex_base_3<Kernel> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<unsigned char, Kernel> Cb1;
typedef CGAL::Surface_mesh_cell_base_3<Kernel, Cb1> Cb2;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<Kernel, Cb2> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Binary_image> Surface_3;

#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
#include <mc/MarchingCubes.h>
#endif

class MainWindow;
class QDoubleSpinBox;
class Viewer;
class Values_list;

class Volume : public Surface
{
  Q_OBJECT
public:
  Volume(MainWindow* mw);
  ~Volume();

private:
  Binary_image m_image;

  // options
  FT m_sm_angle;
  FT m_sm_radius;
  FT m_sm_distance;
  double m_relative_precision;

  // visualization
  bool m_view_surface;
  bool m_draw_triangulation;
  QColor m_triangulation_color;
  bool m_inverse_normals;
  bool two_sides;
  bool draw_triangles_edges;
  bool use_gouraud;
  bool show_bbox;

  std::vector<Facet> m_surface;
  Tr del;            // 3D-Delaunay triangulation
  C2t3 c2t3;         // 2D complex in 3D triangulation

  MainWindow* mw;
  QFileInfo fileinfo;
  Values_list* values_list;
  QDoubleSpinBox* spinBox_radius_bound;
  QDoubleSpinBox* spinBox_distance_bound;

  bool direct_draw; // do not use display lists
  std::vector<GLuint> lists_draw_surface;
  bool lists_draw_surface_is_valid;
  GLuint list_draw_marching_cube;
  bool list_draw_marching_cube_is_valid;

  CGAL::Timer sm_timer;
  int sm_total_time;

#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  std::vector<Facet> m_surface_mc;
  MarchingCubes<unsigned char> mc ;
  std::vector<int> nbs_of_mc_triangles;
  std::vector<GLuint> lists_draw_surface_mc;
  bool lists_draw_surface_mc_is_valid;
  CGAL::Timer mc_timer;
  int mc_total_time;
public:
  void gl_draw_surface_mc();
  void gl_draw_marchingcube();
private:
  void gl_draw_one_marching_cube_vertex(int);

#endif // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE

  bool m_view_mc; // that boolean is here even with if
		  // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
                  // is not defined.

#ifdef CGAL_USE_VTK
  vtkImageReader* vtk_reader;
  vtkImageData* vtk_image;
  vtkDICOMImageReader* dicom_reader;
  vtkDemandDrivenPipeline* executive;
  vtkImageGaussianSmooth* smoother;
#endif // CGAL_USE_VTK

public Q_SLOTS:
void display_marchin_cube();

private:
  template <typename Iterator>
  void gl_draw_surface(Iterator begin, Iterator end, const QTreeWidgetItem* = 0);

  template <typename PointsOutputIterator,
	    typename DomainsOutputIterator,
	    typename TransformOperator>
  void search_for_connected_components(PointsOutputIterator,
				       DomainsOutputIterator, 
				       TransformOperator);

public:
  void gl_draw_surface();

Q_SIGNALS:

  void new_bounding_box(double, double, double, double, double, double);

public Q_SLOTS:
  void only_in();
  void set_inverse_normals(const bool);
  void set_two_sides(const bool);
  void set_draw_triangles_edges(const bool);
  void set_triangulation_edges_color();
  void set_draw_triangulation(const bool);
  void set_use_gouraud(const bool);
  void set_show_bbox(const bool);
  bool open(const QString& filename);
#ifdef CGAL_USE_VTK
  bool open_vtk(const QString& filename);
#endif
  bool open_xt(const QString& filename);
  bool opendir(const QString& dirname);
  void finish_open();
  void export_off();
  void save_image_to_inr();
  void check_can_export_off();
  void draw();
  void get_bbox(float& /*xmin*/, float& /*ymin*/, float& /*zmin*/,
		float& /*xmax*/, float& /*ymax*/, float& /*zmax*/) {}
  void close() {}
  void display_surface_mesher_result();
  void set_radius_bound(double);
  void set_distance_bound(double);
  void changed_parameters();

  void labellizedToogled(bool);

  void save_image_settings(QString);
  void load_image_settings(QString);
private:
  void status_message(QString);
  void busy() const;
  void not_busy() const;
};

template <typename PointsOutputIterator,
	  typename DomainsOutputIterator,
	  typename TransformOperator>
void Volume::search_for_connected_components(PointsOutputIterator it,
					     DomainsOutputIterator dom_it,
					     TransformOperator transform)
{
  const std::size_t nx = m_image.xdim();
  const std::size_t ny = m_image.ydim();
  const std::size_t nz = m_image.zdim();

  const double max_v = (std::max)((std::max)(m_image.vx(),
                                             m_image.vy()),
                                  m_image.vz());

  typedef unsigned char Marker;
  typedef typename TransformOperator::result_type Label;

  boost::multi_array<Marker, 3> visited(boost::extents[nx][ny][nz]);
  typedef boost::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
                                Indices;
  typedef std::queue<Indices>   Indices_queue;
  typedef std::vector<Indices>  Border_vector;

  int number_of_connected_components = 0;
  for(std::size_t i=0; i<nx; i++)
    for(std::size_t j=0; j<ny; j++)
      for(std::size_t k=0; k<nz; k++)
      {
        if(visited[i][j][k]>0)
          continue;
        const Label current_label = transform(m_image.value(i, j, k));
	*dom_it++ = current_label;
        if(current_label == Label()) {
          visited[i][j][k] = 3;
          continue;
        }

        // if we reach here, (i, j, k) is a new connected component
        ++number_of_connected_components;
        std::cerr << boost::format("Found new connected component (#%5%) "
                                   "at voxel (%1%, %2%, %3%), value=%4%, volume id=%6%\n")
          % i % j % k
          % m_image.value(i, j, k) 
          % number_of_connected_components
          % (int)current_label;

        int nb_voxels = 0;

        Indices_queue queue;
        Indices indices(i, j ,k, 0);
        queue.push(indices);

        Border_vector border;

        /*
         * First pass is a BFS to retrieve all the connected component, and
         * its border.
         * Second pass is a BFS initialized with all voxel of the border.
         * The last voxel of that BFS is used as the seed.
         */
        int pass = 1; // pass will be equal to 2 in second pass

        Indices bbox_min = indices;
        Indices bbox_max = indices;

        while(!queue.empty()) // walk through the connected component
        {
          Indices indices = queue.front();
          queue.pop();

          // warning: those indices i, j and k are local to the while loop
          const std::size_t i = boost::get<0>(indices);
          const std::size_t j = boost::get<1>(indices);
          const std::size_t k = boost::get<2>(indices);
          const std::size_t depth = boost::get<3>(indices);

          if(visited[i][j][k] < pass)
          {
            visited[i][j][k] = pass;
            if(pass == 1 )
            {
              ++nb_voxels;
              boost::get<0>(bbox_min) = (std::min)(i, boost::get<0>(bbox_min));
              boost::get<0>(bbox_max) = (std::max)(i, boost::get<0>(bbox_max));
              boost::get<1>(bbox_min) = (std::min)(j, boost::get<1>(bbox_min));
              boost::get<1>(bbox_max) = (std::max)(j, boost::get<1>(bbox_max));
              boost::get<2>(bbox_min) = (std::min)(k, boost::get<2>(bbox_min));
              boost::get<2>(bbox_max) = (std::max)(k, boost::get<2>(bbox_max));
            }

            static const int neighbors_offset[6][3] = { { +1,  0,  0 },
                                                        { -1,  0,  0 },
                                                        {  0, +1,  0 },
                                                        {  0, -1,  0 },
                                                        {  0,  0, +1 },
                                                        {  0,  0, -1 } };
            bool voxel_is_on_border = false;

            // Visit neighbors.
            // (i_n, j_n, k_n) are indices of neighbors.
            for(int n = 0; n < 6; ++n)
            {
              const ptrdiff_t i_n = i + neighbors_offset[n][0];
              const ptrdiff_t j_n = j + neighbors_offset[n][1];
              const ptrdiff_t k_n = k + neighbors_offset[n][2];
              if(i_n < 0 || i_n >= static_cast<ptrdiff_t>(nx) ||
                 j_n < 0 || j_n >= static_cast<ptrdiff_t>(ny) ||
                 k_n < 0 || k_n >= static_cast<ptrdiff_t>(nz))
              {
                voxel_is_on_border = true;
                continue;
              }
              else
              {
                if(transform(m_image.value(i_n, j_n, k_n)) == current_label)
                {
                  if(visited[i_n][j_n][k_n] < pass) {
                    Indices indices(i_n, j_n, k_n, depth+1);
                    queue.push(indices);
                  }
                }
                else
                  voxel_is_on_border = true;
              }
            } // end for neighbors

            if(pass == 1 && voxel_is_on_border)
              border.push_back(indices);
          } // end if voxel not already visited

          if(queue.empty()) {
            if(pass == 1)
            { // End of first pass. Begin second pass with the voxels of
              // the border.
              for(typename Border_vector::const_iterator
                    border_it = border.begin(), border_end = border.end();
                  border_it != border_end; ++border_it)
                queue.push(*border_it);
              pass = 2;
            }
            else // end of second pass, return the last visited voxel
            {
// 	      if(nb_voxels >= 100)
	      {
		*it++ = std::make_pair(m_image.point(i, j, k), (depth+1)*max_v);
		std::cerr << boost::format("Found seed %5%, which is voxel (%1%, %2%, %3%), value=%4%\n")
		  % i % j % k %  m_image.value(i, j, k) % m_image.point(i, j, k);
	      }
            }
          } // end if queue.empty()
        } // end while !queue.empty() (with local indices i, j, k)

        std::cerr << boost::format("There was %1% voxel(s) in that component.\n"
                                   "The bounding box is (%2% %3% %4%, %5% %6% %7%).\n"
                                   "%8% voxel(s) on border\n")
          % nb_voxels
          % boost::get<0>(bbox_min) % boost::get<1>(bbox_min) % boost::get<2>(bbox_min)
          % boost::get<0>(bbox_max) % boost::get<1>(bbox_max) % boost::get<2>(bbox_max)
          % border.size();
      } // end for i,j,k
} // end function Volume::search_for_connected_components()

#endif // _VOLUME_H
