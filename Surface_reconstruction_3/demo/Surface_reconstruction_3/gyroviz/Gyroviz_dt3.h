// Author     : Nader Salman

#ifndef _Gyroviz_dt3_
#define _Gyroviz_dt3_


// the idea in here is to read a specific file format structured as follows
// 3 columns corresponding to the 3D coordinates of a track followed by the
// list of the images in which he has been "seen".

// Kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Delaunay triangulation 3 && Intersections
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/triangulation_ds_cell_base_3.h>
#include <CGAL/intersections.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include "Gyroviz_info_for_dt3.h"
#include "Gyroviz_point_dt3.h"

#include <list>
#include <vector>
#include <sstream>

template <class Gt, class Tds>
class Gyroviz_dt3 : public CGAL::Delaunay_triangulation_3<Gt, Tds>
{

  // private types 
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



  // tetrahedron type
  static const unsigned char UNKNOWN  = -1;
  static const unsigned char OUTSIDE  =  0;
  static const unsigned char INSIDE   =  1;

  // Data members
private:

  Iso_cuboid_3 m_bounding_box;       // Triangulation's bounding box
  Point_3      m_barycenter;         // Triangulation's barycenter
  FT           m_standard_deviation; // Triangulation's standard deviation


  // Public methods
public:

  // Default constructor, copy constructor and operator =() are fine.

  // If some vertex belongs to the final surface then it should be visible
  // in the views it comes from.
  // All the tetrahedra intersected by a ray emanating from the vertex to
  // the camera center of one of these views should be labbelled as outside
  // (0) the tetrahedron behind the vertex should be labelled as inside (1).


  /// Insert points to the triangulation using a spatial sort.
  /// Return the number of inserted points.
  Base::insert;
  int insert(typename std::vector<Gyroviz_point_dt3>::iterator first,
    typename std::vector<Gyroviz_point_dt3>::iterator last)
  {
    // spatial sorting
    std::random_shuffle(first,last);
    //typedef typename CGAL::Gyroviz_point_dt3_spatial_sort_traits<Geom_traits> Spatial_sort_traits;
    //CGAL::spatial_sort(first,last,Spatial_sort_traits());
    CGAL::spatial_sort(first,last,this->geom_traits());

    Cell_handle cell;
    unsigned int index = 0;
    typename std::vector<Gyroviz_point_dt3>::const_iterator it;
    for(it = first; it != last; it++, index++)
    {
      const Gyroviz_point_dt3& pwc = *it;
      Vertex_handle v = insert(pwc,cell);
      v->info() = Gyroviz_info_for_dt3(pwc.cameras());
      cell = v->cell();
    }

    update_bounding_box();

    return index;
  }


  bool save_pnt(char *pFilename)
  {
  }


  // Read 3D points of an OFF file.
  bool read_pwc(const char* pFilename)
  {
    FILE *pFile = fopen(pFilename,"rt");
    if(pFile == NULL)
    {
      std::cerr << "Error: cannot open " << pFilename;
      return false;
    }

    long pointsCount = 0, camCount = 0; // number of vertices and number of cameras in the file
    int lineNumber = 0; // current line number

    char pLine[4096]; // current line buffer
    //char rest_of_line [512]; // current camera list buffer

    std::vector<Point_3> list_of_camera_coordinates; // container of cameras
    std::vector<Gyroviz_point_dt3> gpts; // container of points to be read
    double x,y,z;
    double Cx,Cy,Cz;   


    while(fgets(pLine,sizeof(pLine),pFile))
    {
      lineNumber++;

      // Read file signature on first line
      if (lineNumber == 1)
      {
        char signature[512];

        if ( (sscanf(pLine,"%s",signature) != 1) || (strcmp(signature, "CamOFF") != 0) )
        {
          // unsupported file format
          std::cerr << "Incorrect file format line " << lineNumber << " of " << pFilename;
          return false;
        }
      }    

      // Read number of vertices and cameras on 2nd line    
      else if (lineNumber == 2)
      {
        if (sscanf(pLine,"%ld %ld",&camCount,&pointsCount) != 2)
        {
          std::cerr << "Error line " << lineNumber << " of " << pFilename;
          return false;
        }
      }       

      // Read 3D points on next lines
      else if (list_of_camera_coordinates.size() < camCount)
      {
        if(sscanf(pLine,"%lg\t%lg\t%lg",&Cx,&Cy,&Cz) == 3)
        {
          Point_3 cam_coord(Cx,Cy,Cz);
          list_of_camera_coordinates.push_back(cam_coord);
        }
        // ...or skip comment line
      }

      // Read 3D points + camera indices on next lines
      else {

        std::istringstream iss(pLine);
        Point_3 position;

        if (iss >> position)
        {
          //TRACE("Line number : %d, List_cameras : ", lineNumber);
          int value;
          std::vector<Point_3> list_of_cameras;
          while (iss >> value)
          {
            Point_3 cam_coord = list_of_camera_coordinates[value-1];
            list_of_cameras.push_back(cam_coord);
            // TRACE("%d ", value);
          } 

          // TRACE("\n", lineNumber);

          // vertex contains a 3D point followed by a vector of 3D
          // points corresponding to the camera positions.


          // TEMPORARY? Skip 3D points with no cameras.
          if (list_of_cameras.begin() != list_of_cameras.end())
          {
            //Point_3 point_3(x,y,z);
            //Vertex_handle vh = this->insert(position);
            //vh->info()       = Gyroviz_info_for_dt3(list_of_cameras);  
            //// TEST Each cell corresponding to the current vertex is flagged as UNKNOWN
            //Cell_handle cell;
            //cell = vh->cell();

            //if(!is_infinite(cell)){cell->info() = INSIDE;}
            //else cell->info() = OUTSIDE;
            gpts.push_back(Gyroviz_point_dt3(position,list_of_cameras.begin(),
              list_of_cameras.end()));
          }
        }
      }
    }
    fclose(pFile);

    // Insert points
    insert(gpts.begin(), gpts.end());
    update_bounding_box();

    this->inside_outside(); // label correctly all the cells.
    return true;
  }



  // Get indices different from i and j
  void other_two_indices(int i, int j, int* k, int* l)
  {
    CGAL_surface_reconstruction_assertion(i != j);
    bool k_done = false;
    bool l_done = false;
    for(int index=0;index<4;index++)
    {
      if(index != i && index != j)
      {
        if(!k_done)
        {
          *k = index;
          k_done = true;
        }
        else
        {
          *l = index;
          l_done = true;
        }
      }
    }
    CGAL_surface_reconstruction_assertion(k_done);
    CGAL_surface_reconstruction_assertion(l_done);
  }



  // Get last indice different from the 3 indices in the vector
  void other_indice(std::vector<int> list_vertices, int* l)
  {
    CGAL_surface_reconstruction_assertion(list_vertices.size() == 3);
    bool l_done = false;
    int sum_of_vertices;


    for(int i = 0; i<list_vertices.size(); ++i)
    {
      sum_of_vertices += list_vertices[i];
    }

    if(sum_of_vertices == 3) { *l=3; l_done = true; } 
    else if(sum_of_vertices == 4) { *l=2; l_done = true; }
    else if(sum_of_vertices == 5) { *l=1; l_done = true; }
    else { *l=0; l_done = true; }

    CGAL_surface_reconstruction_assertion(l_done);
  }


  enum Type{NOTHING, VERTEX, EDGE, FACET};

  void do_intersect_from_facet(Facet f, Segment_3 S, Type& out_type, int& out_i, int& out_j)
  {
    // initialization
    out_type = NOTHING;
    out_i = out_j = -1;

    std::vector<int> intersected_facets;

    if(is_infinite(f.first)) return;

    // count number of inttersection
    for(int k = 0; k<4; ++k) 
    {
      if(k != f.second)
      {
        Triangle_3 cell_facet = this->triangle(f.first,k);   

        if(do_intersect(cell_facet,S))
        {
          intersected_facets.push_back(k);
        }
      }
    }

    int number_intersected_facets = intersected_facets.size();

    if(number_intersected_facets == 0) 
      out_type = NOTHING;

    else if(number_intersected_facets == 1)
    {
      out_type = FACET;
      out_i = intersected_facets[0];
    }
    else if(number_intersected_facets == 2)
    {
      other_two_indices(intersected_facets[0], intersected_facets[1], &out_i, &out_j);
      out_type = EDGE;
    }
    else
    {
      out_type = VERTEX;
      out_i = f.second; 
    }
  }



  void do_intersect_from_vertex(Cell_handle c, int index, Segment_3 S, Type& out_type, int& out_i, int& out_j)
  {
    // initialization
    out_type = NOTHING;
    out_i = out_j = -1;
    Line_3 L = S.supporting_line();

    if(is_infinite(c)) return;

    // opposite facet to the input vertex(index)
    Triangle_3 cell_facet = this->triangle(c,index);

    if(!do_intersect(cell_facet,S))
    {
      out_type = NOTHING;
      return;
    }

    for(int k = 0; k<4; ++k) 
    {
      if(k != index)
      {
        // tests if S intersects one of the remaining vertices
        if(S.has_on(c->vertex(k)->point()))
        {
          out_type = VERTEX;
          out_i = k;
          return;
        }
      }
    }

    for(int k = 0; k<4; ++k) 
    {
      if(k != index)
      {
        // tests if S is coplanar to one of the 3 incident facets to the vertex
        Plane_3 plane(c->vertex( (k+1)&3 )->point(), 
          c->vertex( (k+2)&3 )->point(),
          c->vertex( (k+3)&3 )->point());

        if(plane.has_on(L))
        {
          out_type = EDGE;
          other_two_indices(index, k, &out_i, &out_j);
          return;
        }
      }
    }

    // else we do intersect the facet index
    out_type = FACET;
    out_i = index;    
  }



  void do_intersect_from_edge(Edge e, Segment_3 S, Type& out_type, int& out_i, int& out_j)
  {

    // initialization
    out_type = NOTHING;
    out_i = out_j = -1;

    Cell_handle c = e.first;
    int in_i = e.second, in_j = e.third;
    int k,l;

    Triangle_3 cell_facet_i = this->triangle(c,in_i);
    Triangle_3 cell_facet_j = this->triangle(c,in_j);

    other_two_indices(in_i,in_j,&k,&l);

    if(is_infinite(c)) return;

    // this code can be optimized
    if(!do_intersect(cell_facet_i,S) && !do_intersect(cell_facet_j,S))
    {
      out_type = NOTHING;
      return;
    }

    // tests if S intersects one of the remaining vertices
    else if(S.has_on(c->vertex(k)->point()))
    {
      out_type = VERTEX;
      out_i = k;
    }

    // tests if S intersects one of the remaining vertices
    else if(S.has_on(c->vertex(l)->point()))
    {
      out_type = VERTEX;
      out_i = l;
    }

    else if(do_intersect(cell_facet_i,S) && do_intersect(cell_facet_j,S))
    {
      out_type = EDGE;
      out_i    = k;
      out_j    = l; 
    }

    else if(do_intersect(cell_facet_i,S))
    {
      out_type = FACET;
      out_i = in_i; 
    }

    else
    {
      out_type = FACET;
      out_i = in_j; 
    }
  }



  void inside_outside()
  {
    
    // lets initialize all the cells to INSIDE
    for (Cell_iterator cell = cells_begin(); cell != cells_end(); cell++)
    {
      // TEST Each cell corresponding to the current vertex is flagged as "UNKNOWN"
      if(!is_infinite(cell)){cell->info() = INSIDE;}
      else cell->info() = OUTSIDE;
      
    }

    // TEST 
    int counter_vertex = 0;

    Finite_vertices_iterator fv = this->finite_vertices_begin();

    for(;fv != this->finite_vertices_end(); ++fv)
    {
      counter_vertex++;

      std::vector<Point_3> list_of_cams = fv->info().get_list_of_cameras();

      TRACE("Vertex : %d, List_of_cameras_size : %d\n",counter_vertex,fv->info().get_list_of_cameras().size());

      // construct segments emanating from the vertex (fv) to cameras.
      for(int j = 0; j<fv->info().get_list_of_cameras().size(); ++j)
      {

        Segment_3  segment (fv->point(), list_of_cams[j]);

        //////////////////////////////////
        // input variable initialization
        Cell_handle in_Tetra = NULL;
        Type        in_type  = NOTHING;
        int in_i = -1, in_j = -1;
        /////////////////////////////////

        // TODO : Opposite tetra must be labbelled as INSIDE!

        std::list<Cell_handle> incident_c;
        incident_cells(static_cast<Vertex_handle>(fv),back_inserter(incident_c));

        // create an iterator for these cells
        typename std::list<Cell_handle>::iterator cell_it = incident_c.begin();


        //// TEST
        //Cell_handle INSIDE_cell, OUTSIDE_cell;
        //
        //for( ; cell_it != incident_c.end(); ++cell_it) 
        //{
        //  if(is_infinite(*cell_it)) continue;

        //  int index_in_cell_it;
        //  if((*cell_it)->has_vertex(fv,index_in_cell_it))
        //  {
        //    // opposite facet to the input vertex(index)
        //    Triangle_3 opposite_cell_facet = this->triangle(*cell_it,index_in_cell_it);

        //    if(do_intersect(opposite_cell_facet,segment.opposite()))
        //    {
        //     
        //      (*cell_it)->info() = INSIDE;
        //      
        //      // TEST
        //      TRACE("INSIDE Cell : %x\n",&*cell_it);
        //      INSIDE_cell = *cell_it;
        //      
        //      break;
        //    }
        //  }
        //}

        cell_it = incident_c.begin();

        // lets find the first Cell

        for( ; cell_it != incident_c.end(); ++cell_it) 
        {
          int index_in_cell_it;
          if((*cell_it)->has_vertex(fv,index_in_cell_it))
          {
            do_intersect_from_vertex(*cell_it, index_in_cell_it, segment, in_type, in_i, in_j);

            if(in_type != NOTHING)
            {
              in_Tetra = (*cell_it);

              ////TEST
              //TRACE("OUTSIDE Cell : %x\n",&in_Tetra);
              //OUTSIDE_cell = in_Tetra;

              break;
            }
          }
        }

        ////TEST
        //assert(&OUTSIDE_cell != &INSIDE_cell);


        // (in_Tetra, in_type, in_i, in_j) variables are initialized for the upcoming tests
        while(in_type != NOTHING)
        {

          // compute the output (out_Tetra, out_type, out_i, out_j) variables of in_Tetra

          // output variable initialization
          Cell_handle out_Tetra = in_Tetra;
          Type        out_type  = NOTHING;
          int out_i = -1, out_j = -1;

          switch(in_type)
          {
          case FACET:
            {
              //TRACE("Case in_type == FACET\n");
              Facet F(in_Tetra,in_i);
              do_intersect_from_facet(F, segment, out_type, out_i, out_j);
              break;
            }

          case EDGE:
            {
              //TRACE("Case in_type == EDGE\n");
              Edge E(in_Tetra,in_i,in_j);
              do_intersect_from_edge(E, segment, out_type, out_i, out_j);
              break;
            }

          case VERTEX:
            {
              //TRACE("Case in_type == VERTEX\n");
              do_intersect_from_vertex(in_Tetra, in_i, segment, out_type, out_i, out_j);
              break;
            }
          }

          // TODO : All the tetrahedra intersected by a segment emanating from the vertex to
          //        the camera center of one of these views should be labbelled as outside.
          //        Exception : segment is tangent to the tetrahedra

          //TRACE("OUTSIDE\n");
          in_Tetra->info() = OUTSIDE;

          // compute the next (in_Tetra, in_type, in_i, in_j) variables
          in_type = out_type;

          switch(out_type)
          {

          case NOTHING : 
            break;

          case FACET:
            {
              //TRACE("Case out_type == FACET\n");
              in_Tetra = out_Tetra->neighbor(out_i);
              in_type  = FACET;
              in_i     = this->mirror_index(out_Tetra,out_i);

              break;
            }

          case EDGE:
            {
              //TRACE("Case out_type == EDGE\n");
              Edge E(out_Tetra, out_i, out_j);

              Cell_circulator ccir = this->incident_cells(E, out_Tetra);
              Cell_circulator cdone = ccir;

              int index_in_ccir_i, index_in_ccir_j; 

              do{

                if((Cell_handle)ccir != out_Tetra)
                {
                  if((*ccir).has_vertex(out_Tetra->vertex(out_i),index_in_ccir_i)){}

                  if((*ccir).has_vertex(out_Tetra->vertex(out_j),index_in_ccir_j)){}

                  Edge edge_ccir(&ccir, index_in_ccir_i, index_in_ccir_j);

                  do_intersect_from_edge(edge_ccir, segment, in_type, in_i, in_j);
                  if(in_type != NOTHING)
                  {
                    break;
                  }
                  ++ccir;
                }
                ++ccir;
              }while ( ccir != cdone );

              break;
            }

          case VERTEX:
            {
              //TRACE("Case out_type == VERTEX\n");
              Vertex_handle v = out_Tetra->vertex(out_i);

              std::list<Cell_handle> incident_t;
              incident_cells(v,back_inserter(incident_t));

              typename std::list<Cell_handle>::iterator tetra_it = incident_t.begin();

              for( ; tetra_it != incident_t.end(); ++tetra_it) 
              {
                if(*tetra_it != out_Tetra)
                {
                  int index_in_tetra_it;
                  if((*tetra_it)->has_vertex(v,index_in_tetra_it))
                  {
                    do_intersect_from_vertex((*tetra_it), index_in_tetra_it, segment, in_type, in_i, in_j);
                    if(in_type != NOTHING)
                    {
                      break;
                    }
                  }
                }
              }         

              break;
            }   
          }//end::switch(out_type)
        }//end::while(in_type != NOTHING)             
      }//end::list_of_cameras
    }//end::for_all_vertices
  }//end::inside_outside





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



  void gl_draw_3D_vertices(const unsigned char r, const unsigned char g,
    const unsigned char b, const float size)
  {
    // Draw input points
    ::glPointSize(size);
    ::glColor3ub(r,g,b);
    ::glBegin(GL_POINTS);

    Finite_vertices_iterator fv = this->finite_vertices_begin();
    for(; fv != this->finite_vertices_end(); ++fv)
    {
      Point_3 p = fv->point();
      ::glVertex3d(p.x(),p.y(),p.z());
    }

    ::glEnd();
  }


  void gl_draw_3D_delaunay_triangulation(const unsigned char r, const unsigned char g,
    const unsigned char b, const float width)
  {
    // Draw 3D delaunay triangulation
    ::glLineWidth(width);    
    ::glColor3ub(r,g,b);
    ::glBegin(GL_LINES);

    Finite_edges_iterator fe = this->finite_edges_begin();
    for(; fe != this->finite_edges_end(); ++fe)
    {
      Segment s = segment(*fe);
      Point_3 p1 = s.source();
      Point_3 p2 = s.target();
      ::glVertex3d(p1.x(),p1.y(),p1.z());
      ::glVertex3d(p2.x(),p2.y(),p2.z());
    } 

    ::glEnd();
  } 


  void gl_draw_3D_rays(const unsigned char r, const unsigned char g,
    const unsigned char b)
  {

    // Draw 3D delaunay triangulation
    ::glColor3ub(r,g,b);
    ::glBegin(GL_LINES);

    Finite_vertices_iterator fv = this->finite_vertices_begin();
    for(; fv != this->finite_vertices_end(); ++fv)
    {

      Point_3  vertex_source       = fv->point();
      std::vector<Point_3> list_of_cams = fv->info().get_list_of_cameras();

      for(int i=0; i<list_of_cams.size(); ++i)
      {
        Point_3 camera  = list_of_cams[i];
        ::glVertex3d(vertex_source.x(), vertex_source.y(), vertex_source.z());
        ::glVertex3d(camera.x(), camera.y(), camera.z());
      }
    }

    ::glEnd();
  } 


  void gl_draw_3D_inside_tetrahedrons(const unsigned char r, const unsigned char g,
    const unsigned char b)
  {

    // Draw 3D delaunay triangulation
    ::glColor3ub(r,g,b);
    ::glBegin(GL_TRIANGLES);

    Facet_iterator ff = this->facets_begin();
    for(; ff != this->facets_end(); ++ff)
    {
      Cell_handle c  = ff->first;
      int index_in_c = ff->second;

      Cell_handle neighbor_c = c->neighbor(index_in_c);

      if(!is_infinite(c)&& !is_infinite(neighbor_c) && c->info() == INSIDE && neighbor_c->info() == OUTSIDE
        || !is_infinite(c) && !is_infinite(neighbor_c) && c->info() == OUTSIDE && neighbor_c->info() == INSIDE)
      {
        int k = ff->second;
        Point_3 p1 = c->vertex( (k+1)&3 )->point();
        Point_3 p2 = c->vertex( (k+2)&3 )->point();
        Point_3 p3 = c->vertex( (k+3)&3 )->point();
        ::glVertex3d(p1.x(), p1.y(), p1.z());
        ::glVertex3d(p2.x(), p2.y(), p2.z());
        ::glVertex3d(p3.x(), p3.y(), p3.z());
      }
    } 

    ::glEnd();
  }    

};


#endif // _Gyroviz_dt3_