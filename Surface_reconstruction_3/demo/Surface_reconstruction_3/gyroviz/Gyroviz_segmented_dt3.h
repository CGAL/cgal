// Author     : Nader Salman

#ifndef _Gyroviz_segmented_dt3_
#define _Gyroviz_segmented_dt3_


// Delaunay triangulation 3 && Intersections
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/triangulation_ds_cell_base_3.h>
#include <CGAL/intersections.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/surface_reconstruction_assertions.h>


#include <list>
#include <vector>
#include <sstream>

#include <CImg.h>
using namespace cimg_library;


template <class Gt, class Tds>
class Gyroviz_segmented_dt3 : public CGAL::Delaunay_triangulation_3<Gt, Tds>
{

private:

  typedef CGAL::Delaunay_triangulation_3 <Gt, Tds> Base;

public:

  typedef Tds Triangulation_data_structure;
  typedef Gt Geom_traits;

  typedef typename Geom_traits::FT                FT;
  typedef typename Geom_traits::Ray_3             Ray_3;
  typedef typename Geom_traits::Plane_3           Plane_3;
  typedef typename Geom_traits::Point_3           Point_3;
  typedef Gyroviz_point_dt3<Geom_traits>          Gyroviz_point_dt3;
  typedef typename Geom_traits::Line_3            Line_3;
  typedef typename Geom_traits::Vector_3          Vector_3;
  typedef typename Geom_traits::Segment_3         Segment_3;
  typedef typename Geom_traits::Triangle_3        Triangle_3;
  typedef typename Geom_traits::Tetrahedron_3     Tetrahedron;
  typedef typename Geom_traits::Sphere_3          Sphere;
  typedef typename Geom_traits::Iso_cuboid_3      Iso_cuboid_3; 

  typedef typename Base::Vertex                   Vertex;
  typedef typename Base::Facet                    Facet;
  typedef typename Base::Edge                     Edge;
  typedef typename Base::Cell                     Cell;

  typedef typename Base::Vertex_handle            Vertex_handle;
  typedef typename Base::Cell_handle              Cell_handle;
  typedef typename Base::Locate_type              Locate_type;

  typedef typename Base::Cell_circulator          Cell_circulator;
  typedef typename Base::Facet_circulator         Facet_circulator;
  typedef typename Base::Cell_iterator            Cell_iterator;
  typedef typename Base::Facet_iterator           Facet_iterator;
  typedef typename Base::Edge_iterator            Edge_iterator;
  typedef typename Base::Vertex_iterator          Vertex_iterator;
  typedef typename Base::Point_iterator           Point_iterator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Base::All_cells_iterator       All_cells_iterator;
  
  // tetrahedron type (in case we will scult the tetrahedrization)
  static const unsigned char UNKNOWN  = -1;
  static const unsigned char VIDE  =  0;
  static const unsigned char PLEIN   =  1;

  // Data members
private:

  Iso_cuboid_3 m_bounding_box;       // Triangulation's bounding box
  Point_3      m_barycenter;         // Triangulation's barycenter
  FT           m_standard_deviation; // Triangulation's standard deviation
  
  
  // points near detected borders 
  bool on_border_2D_vertices(Point_2 p2, CImg <unsigned char> image)
  {

    // if pixel any of 9x9 surrounding pixels is on border 
    // keep vertex(white color is used in frei-chen gradient operator)

    if (image((unsigned int)p2.x(),(unsigned int)p2.y(),0,0)==255)
    {
      return true;
    }

    else if((unsigned int)p2.x() == 0 && (unsigned int)p2.y() == 0) //upper left pixel
    {
      if (image((unsigned int)p2.x()+1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x(),(unsigned int)p2.y()+1,0,0)  == 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()+1,0,0)== 255)
      {  
        return true;
      }

    }
    
    else if((unsigned int)p2.x() == image.dimx() && (unsigned int)p2.y() == 0) //upper right pixel
    {
      if (image((unsigned int)p2.x()-1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x(),(unsigned int)p2.y()+1,0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()+1,0,0)== 255) 
      {
        return true;
      }
    }

    else if((unsigned int)p2.x() == 0 && (unsigned int)p2.y() == image.dimy()) //lower left pixel
    {
      if (image((unsigned int)p2.x(),(unsigned int)p2.y()-1,0,0)  == 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()-1,0,0)== 255) 
      {
        return true;
      }

    }

    else if((unsigned int)p2.x() == image.dimx() && (unsigned int)p2.y() == image.dimy()) //lower right pixel
    {
      if (image((unsigned int)p2.x(),(unsigned int)p2.y()-1,0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()-1,0,0)== 255) 
      {
        return true;
      }
    }

    else if((unsigned int)p2.x() > 0 && (unsigned int)p2.x() < image.dimx() && (unsigned int)p2.y() == 0) // upper band
    {
      if (image((unsigned int)p2.x()-1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()+1,0,0)== 255 ||
        image((unsigned int)p2.x(),  (unsigned int)p2.y()+1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()+1,0,0)== 255) 
      {
        return true;
      }      
    }

    else if((unsigned int)p2.x() > 0 && (unsigned int)p2.x() < image.dimx() && (unsigned int)p2.y() == image.dimy()) // lower band
    {
      if (image((unsigned int)p2.x()-1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x(),  (unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y(),0,0)  == 255) 
      {
        return true;
      }
    }

    else if((unsigned int)p2.x() == 0 && (unsigned int)p2.y() > 0  && (unsigned int)p2.y() < image.dimy()) // left band
    {
      if (image((unsigned int)p2.x(),  (unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()+1,0,0)== 255 ||
        image((unsigned int)p2.x(),  (unsigned int)p2.y()+1,0,0)== 255) 
      {
        return true;
      }
    }

    else if((unsigned int)p2.x() == image.dimx() && (unsigned int)p2.y() > 0  && (unsigned int)p2.y() < image.dimy()) // right band
    {
      if (image((unsigned int)p2.x(),  (unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()+1,0,0)== 255 ||
        image((unsigned int)p2.x(),  (unsigned int)p2.y()+1,0,0)== 255) 
      {
        return true;
      }
    }

    else // middle of the image corner and bands excluded
    {
      if (image((unsigned int)p2.x()-1,(unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()-1,(unsigned int)p2.y()+1,0,0)== 255 ||
        image((unsigned int)p2.x(),  (unsigned int)p2.y()+1,0,0)== 255 ||
        image((unsigned int)p2.x(),  (unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()-1,0,0)== 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y(),0,0)  == 255 ||
        image((unsigned int)p2.x()+1,(unsigned int)p2.y()+1,0,0)== 255)
      {
        return true;
      }
    }
    return false;
  }



 bool read_pnt_and_image(char *pFilename, CImg <unsigned char> image, int image_number)
  {
    //extract from pFilename the image number  
    std::string temp ( pFilename );
    std::string filename_without_path = temp.substr(temp.size()-15);
    std::string extract_number = filename_without_path.substr(7,filename_without_path.size()-4);
    
    int image_number = atoi(extract_number.c_str());
    
    //assert(image_number == );
    
    FILE *pFile = fopen(pFilename,"r");
    if(pFile == NULL)
      return false;

    //scan vertices and add them to triangulation with corresponding info
    int lineNumber = 0;
    char pLine[512];

    while( fgets(pLine, 512, pFile))
    {
      lineNumber++;

      // read 2D/3D coordinates 
      //(on suppose avoir fait du traitement 
      // des fichiers pnt a l'etape precedente)
      int  unused1, is_reconstructed;
      double p,q,x,y,z;

      if (sscanf(pLine,"%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf", &p,&q,&unused1,&is_reconstructed,&x,&y,&z) == 7)
      {
        if (is_reconstructed == 2)
        {
          Point_2 point_2(p,q);
          Point_3 point_3(x,y,z);

         // test if point_2 is on one of the borders of the segmented image
         if(on_border_2D_vertices(point_2, image))
         {
           Vertex_handle vh = this->insert(point_3);
           vh->info() = image_number;
         }

        }
      }
    }
    fclose(pFile);
    update_bounding_box();
    return (this->number_of_vertices() > 0);
  }


bool read_set_of_pnt(char *pFilename, CImg <unsigned char> image, int image_number)
  {





  /// Get the region of interest, ignoring the outliers.
  /// This method is used to define the OpenGL arcball sphere.
  Sphere region_of_interest() const
  {
    // A good candidate is a sphere containing the dense region of the point cloud:
    // - center point is barycenter
    // - Radius is 2 * standard deviation
    float radius = 2.f * (float)m_standard_deviation;
    return Sphere(m_barycenter, radius*radius);
  }

  /// Update region of interest.
  /// Owner is responsible to call this function after modifying the triangulation.
  void update_bounding_box()
  {
    // Update bounding box and barycenter.
    // TODO: we should use the functions in PCA component instead.
    FT xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = ymin = zmin =  1e38;
    xmax = ymax = zmax = -1e38;
    Vector_3 v = CGAL::NULL_VECTOR;
    FT norm = 0;
    assert(points_begin() != points_end());
    Finite_vertices_iterator fv = this->finite_vertices_begin();
    for(; fv != this->finite_vertices_end(); ++fv)
    {
      const Point_3& p = fv->point();

      // update bbox
      xmin = (std::min)(p.x(),xmin);
      ymin = (std::min)(p.y(),ymin);
      zmin = (std::min)(p.z(),zmin);
      xmax = (std::max)(p.x(),xmax);
      ymax = (std::max)(p.y(),ymax);
      zmax = (std::max)(p.z(),zmax);

      // update barycenter
      v = v + (p - CGAL::ORIGIN);
      norm += 1;
    }
    //
    Point_3 p(xmin,ymin,zmin);
    Point_3 q(xmax,ymax,zmax);
    m_bounding_box = Iso_cuboid_3(p,q);
    //
    m_barycenter = CGAL::ORIGIN + v / norm;

    /// Compute standard deviation
    Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    /*Finite_vertices_iterator*/ fv = this->finite_vertices_begin();
    for(; fv != this->finite_vertices_end(); ++fv)
    {
      const Point_3& p = fv->point();
      sq_radius += sqd(p, m_barycenter);
    }
    sq_radius /= number_of_vertices();
    m_standard_deviation = CGAL::sqrt(sq_radius);
  }















#endif // _Gyroviz_segmented_dt3_