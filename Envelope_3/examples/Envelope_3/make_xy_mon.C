
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>



#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h> 

#include <queue>
#include <list>
#include <vector>
#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL/Timer.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/IO/Arr_iostream.h>

//typedef double                                    NT;
typedef CGAL::Quotient<CGAL::MP_Float>            NT;
typedef CGAL::Simple_cartesian<NT>                Kernel;
typedef Kernel::RT                                RT;
typedef Kernel::Point_3                           Point_3; 
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Plane_3                           Plane_3;
typedef Kernel::Line_2                            Line_2;
typedef CGAL::Polyhedron_3<Kernel>                Polyhedron;
typedef Polyhedron::Halfedge_const_handle         Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator          Facet_const_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator 
  Halfedge_around_facet_const_circulator;

typedef Kernel::Oriented_side_2                   Oriented_side_2;
typedef Kernel::Construct_line_2                  Construct_line_2;
typedef Kernel::Construct_projected_xy_point_2    
  Construct_projected_xy_point_2;
typedef Kernel::Construct_plane_3                 Construct_plane_3;
typedef CGAL::Arr_segment_traits_2<Kernel>        Traits_2;
typedef Traits_2::X_monotone_curve_2              X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>             Arrangement_2;
typedef Arrangement_2::Halfedge_handle            Arr_halfedge;


class Make_xy_monotone
{
protected:
  typedef std::list<Arrangement_2*>      XY_surfaces_list;

  std::queue<Facet_const_iterator>                             m_queue;
  CGAL::Unique_hash_map<Facet_const_iterator, bool>            m_status_map;
  CGAL::Unique_hash_map<Halfedge_const_handle, Arr_halfedge>   m_edges_map;

public:
 
  void operator()(const Polyhedron& p, XY_surfaces_list& xy_surfaces)
  {
    m_status_map = 
      CGAL::Unique_hash_map<Facet_const_iterator, bool>(false, p.size_of_facets());

    m_edges_map = CGAL::Unique_hash_map<Halfedge_const_handle,
                                        Arr_halfedge>(Arr_halfedge());

    for (Facet_const_iterator f = p.facets_begin(); f != p.facets_end(); ++f)
    {
      if (! m_status_map[f])
        bfs(f, xy_surfaces);
    }
  }

  void bfs(Facet_const_iterator f, XY_surfaces_list& xy_surfaces)
  {
    m_edges_map.clear(); // clear the previous data of the edges map
    
    m_status_map[f] = true;
    if(is_vertical(f))
      return;

    Arrangement_2* arr = new Arrangement_2;
    xy_surfaces.push_back(arr);
    m_queue.push(f);
    while(!m_queue.empty())
    {
      f = m_queue.front();
      m_queue.pop();

      // add f to current xy_surface
      add_facet_to_arr(*arr, f);
      Halfedge_around_facet_const_circulator haf = f->facet_begin();
      Halfedge_around_facet_const_circulator curr = haf;
      do
      {
        if(curr->opposite()->is_border())
        {
          ++curr;
          continue;
        }
        Facet_const_iterator neighbour_f = curr->opposite()->face();

        bool& is_done = m_status_map[neighbour_f];
        if (!is_done)
        {
          if(is_vertical(neighbour_f))
          {
            is_done = true;
          }
          else
            if(are_xy_mon(curr, f, neighbour_f))
            {
              is_done = true;
              m_queue.push(neighbour_f);
            }
        }
            
        ++curr;
      } while(curr != haf);
    }
  }

  bool is_vertical(Facet_const_iterator f)
  {
    CGAL_assertion(f->facet_degree() > 2);
    Halfedge_const_handle he = f->halfedge();
    const Point_3& p1 = he->vertex()->point();
    const Point_3& p2 = he->next()->vertex()->point();
    const Point_3& p3 = he->next()->next()->vertex()->point();
    Kernel k;
    return (k.collinear_2_object()(project_xy(p1),
                                   project_xy(p2),
                                   project_xy(p3)));
  }
  bool are_xy_mon(Halfedge_const_handle he, Facet_const_iterator f1, Facet_const_iterator f2)
  {
    /*
     *                    t
     *                   *|*
     *                  * | *
     *                 * ^|  *
     *                *  ||   *
     *            p1 * he||    * p2
     *               * f1|| f2 *
     *                *  ||   *
     *                 *  |  *
     *                  * | *
     *                   *|*
     *                    s
     */

    CGAL_assertion(!is_vertical(f1) && !is_vertical(f2));
    
    const Point_2& s = project_xy(he->opposite()->vertex()->point());
    const Point_2& t = project_xy(he->vertex()->point());

    const Point_2& p1 = project_xy(he->next()->vertex()->point());
    const Point_2& p2 = project_xy(he->opposite()->next()->vertex()->point());

    Kernel k;
    Line_2 l = k.construct_line_2_object()(s, t); 
    CGAL::Oriented_side os1, os2;

    os1 = k.oriented_side_2_object()(l, p1);
    os2 = k.oriented_side_2_object()(l, p2);

    CGAL_assertion(os1 != CGAL::ON_ORIENTED_BOUNDARY &&
                   os2 != CGAL::ON_ORIENTED_BOUNDARY);

    return (os1 != os2);
  }


  /*Point_2  project_xy (const Point_3& p)
  {
    Kernel k;
    Plane_3 xy_plane = k.construct_plane_3_object()(RT(0), RT(0), RT(1), RT(0));
    return (k.construct_projected_xy_point_2_object()(xy_plane, p));
  }*/

  Point_2  project_xy (const Point_3& p)
  {
    return (Point_2(p.x(), p.y()));
  }

  bool is_back_facet(Facet_const_iterator f)
  {
    Kernel k;
    Plane_3 supp_plane =
      k.construct_plane_3_object()(f->halfedge()->vertex()->point(),
                                   f->halfedge()->next()->vertex()->point(),
                                   f->halfedge()->next()->next()->vertex()->point());

    return (supp_plane.c() > 0);
  }

  void add_facet_to_arr(Arrangement_2& arr, Facet_const_iterator f)
  {
    //static int i=1;
    bool is_back = is_back_facet(f);
    Halfedge_around_facet_const_circulator haf = f->facet_begin();
    Halfedge_around_facet_const_circulator curr = haf;

    Kernel k;
    do
    {
      Halfedge_const_handle he = curr;
      
      // look for that polyhedron's halfedge at the edges map
      Arr_halfedge arr_he = m_edges_map[he->opposite()];
      Arr_halfedge prev_he; 

      if(arr_he == Arr_halfedge())
      {
        // this edge was not inserted to the arrangement yet

        const Point_2& s = project_xy(curr->opposite()->vertex()->point());
        const Point_2& t = project_xy(curr->vertex()->point());

        X_monotone_curve_2 cv(s, t);

        Arr_halfedge next_he = find_halfedge_before_cw(he, is_back);
        if(prev_he == Arr_halfedge())
          prev_he = find_halfedge_before_cw(he->opposite(), is_back);

        bool next_exist = (next_he!=Arr_halfedge());
        bool prev_exist = (prev_he!=Arr_halfedge());

        CGAL::Comparison_result res = k.compare_xy_2_object()(s, t);

        if(!next_exist && !prev_exist)
        {
          //std::cout<<"("<<i++<<")"<<"insert_in_face_interior: " << cv<<"\n";
          prev_he = arr.insert_in_face_interior(cv, arr.unbounded_face());
          if(prev_he->direction() != res)
            prev_he = prev_he->twin();
        }
        else
          if(!next_exist && prev_exist)
          {
            if(res == CGAL::SMALLER)
            {
              //std::cout<<"("<<i++<<")"<<"insert_from_left_vertex: " << cv<<"\n";
              prev_he = arr.insert_from_left_vertex(cv, prev_he);
              if(prev_he->direction() != res)
                prev_he = prev_he->twin();
            }
            else
            {
              //std::cout<<"("<<i++<<")"<<"insert_from_right_vertex: " << cv<<"\n";
              prev_he = arr.insert_from_right_vertex(cv, prev_he);
              if(prev_he->direction() != res)
                prev_he = prev_he->twin();
            }
          }
          else
            if(next_exist && !prev_exist)
            {
              if(res == CGAL::SMALLER)
              {
               // std::cout<<"("<<i++<<")"<<"insert_from_right_vertex: " << cv<<"\n";
                prev_he = arr.insert_from_right_vertex(cv, next_he);
                if(prev_he->direction() != res)
                  prev_he = prev_he->twin();
              }
              else
              {
                //std::cout<<"("<<i++<<")"<<"insert_from_left_vertex: " << cv<<"\n";
                prev_he = arr.insert_from_left_vertex(cv, next_he);
                if(prev_he->direction() != res)
                  prev_he = prev_he->twin();
              }
            }
            else
            {
              //std::cout<<"("<<i++<<")"<<"insert_at_vertices: " << cv<<"\n";
              //std::cout<<"prev_he: " <<prev_he->source()->point()<<" --> "<< prev_he->target()->point()<<"\n";
              //std::cout<<"next_he: " <<next_he->source()->point()<<" --> "<< next_he->target()->point()<<"\n";
              //TODO: can be done much more efficient using the Arr_accessor !!
              prev_he = arr.insert_at_vertices(cv, prev_he, next_he);
              if(prev_he->direction() != res)
                prev_he = prev_he->twin();
            }
        //map 'he' of the polyhderon to 'prev_he'
        m_edges_map[he] = prev_he;
      }
      else
      {
        prev_he = arr_he->twin();
      }

      CGAL_assertion(prev_he->source()->point() ==
                     project_xy(he->opposite()->vertex()->point()));
      CGAL_assertion(prev_he->target()->point() ==
                     project_xy(he->vertex()->point()));

     
      ++curr;
    } while(curr != haf);
    

  }

  Arr_halfedge find_halfedge_before_cw(Halfedge_const_handle he, bool is_back_face)
  {
    // iterate the halfedges around the target vertex of 'he' counter-clockwise
    // and find the first halfedge that was inserterd to the arrangement.
    Halfedge_const_handle curr;
    if(is_back_face)
      curr = he->prev_on_vertex();
    else
      curr = he->next_on_vertex();
    do
    {
      Arr_halfedge res = m_edges_map[curr];
      if(res != Arr_halfedge())
        return res;

      res = m_edges_map[curr->opposite()];
      if(res != Arr_halfedge())
        return res->twin();

      if(is_back_face)
        curr = curr->prev_on_vertex();
      else
        curr = curr->next_on_vertex();
    } while(curr != he);

    return Arr_halfedge();
  }

};

Point_2  project_xy (const Point_3& p)
{
  return (Point_2(p.x(), p.y()));
}

bool test_xy_mon(std::vector<Facet_const_iterator>& surf)
{
  CGAL::Unique_hash_map<Halfedge_const_handle, bool> edges_hash(false);
  typedef CGAL::Arr_segment_traits_2<Kernel>   Traits;
  typedef  Traits::Curve_2  Curve_2;
  std::list<Curve_2>  curves;
  for(unsigned int i=0; i< surf.size(); ++i)
  {
    Facet_const_iterator curr_f = surf[i];
    Halfedge_around_facet_const_circulator haf = curr_f->facet_begin();
    Halfedge_around_facet_const_circulator curr = haf;
    do
    {
      Halfedge_const_handle he = curr;
      if(!edges_hash[he] && !edges_hash[he->opposite()])
      {
        Curve_2 cv(project_xy(he->vertex()->point()),
                   project_xy(he->opposite()->vertex()->point()));
        curves.push_back(cv);
        edges_hash[he] = edges_hash[he->opposite()] = true;
      }
      
      ++curr;
    }while(curr != haf);
  }

  bool do_x = CGAL::do_curves_intersect(curves.begin(), curves.end());
  return (!do_x);
}
int main(int argc, char** argv)
{
  if(argc<2)
  {
    std::cout<<"Missing off file\n";
    return (-1);
  }

  std::ifstream off_file(argv[1]);
  if(!off_file.is_open())
  {
    std::cout<<"Failed to open file\n";
    return (-1);
  }

  CGAL::Timer t;
  Polyhedron p;
  //off_file >> p;
  // reads a polyhedron from `in' and appends it to P.
  t.start();
  CGAL::scan_OFF(off_file, p, true);
  t.stop();

  if ( (off_file.rdstate() & std::ios::badbit ) != 0 )
    std::cerr << "Error opening "<<argv[1] <<"\n";

  std::cout<<"Reading OFF file took: " <<t.time()<<" seconds\n";
  t.reset();
  if(p.empty())
    return 0;

  Facet_const_iterator f = p.facets_begin();

  std::cout << "|V| = " <<p.size_of_vertices()<<std::endl;
  std::cout << "|E| = " <<p.size_of_halfedges()<<std::endl;
  std::cout << "|F| = " <<p.size_of_facets()<<std::endl;

  std::list<Arrangement_2*> xy_surfs;
  typedef std::list<Arrangement_2*>::iterator Xy_surf_itr;
  Make_xy_monotone make_xy;
  t.start();
  make_xy(p, xy_surfs);
  t.stop();
  std::cout<<"Making the surface xy-monotone took: " <<t.time()<<" seconds\n";


  std::cout<<"There are: " << xy_surfs.size() << " xy-surfaces\n";
  unsigned int i=1;
  for(Xy_surf_itr itr = xy_surfs.begin(); itr != xy_surfs.end(); ++itr, ++i)
  {
    std::cout<<"#" <<i<<" surface : "<<std::endl;
    std::cout<<"  |V| = "<<(*itr)->number_of_vertices()<<std::endl;
    std::cout<<"  |E| = "<<(*itr)->number_of_edges()<<std::endl;
    std::cout<<"  |F| = "<<(*itr)->number_of_faces()<<std::endl;

    char buffer[33];
    buffer[0] = '1' + i-1;
    buffer[1] = '\0';
    std::string prog_name (argv[1]);
    std::string filename("surface");
    std::string prefix (buffer);
    std::string temp(".arr");
    prefix += temp;
    filename += prefix; 
    std::cout << filename<<"\n";
    std::ofstream arr_file(filename.c_str());
    if(!arr_file.is_open())
      std::cout<<"Failed to open: " << filename<<"\n";
    arr_file << **itr;

  }

  return (0);
}
