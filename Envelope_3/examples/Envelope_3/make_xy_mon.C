
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <limits>



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
#include <CGAL/Arr_accessor.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


////typedef double                                    NT;
typedef CGAL::Quotient<CGAL::MP_Float>            NT;
typedef CGAL::Simple_cartesian<NT>                Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::RT                                RT;
typedef Kernel::Point_3                           Point_3; 
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Plane_3                           Plane_3;
typedef Kernel::Line_2                            Line_2;
typedef Kernel::Segment_2                         Segment_2;
typedef Kernel::Iso_rectangle_2                   Iso_rectangle_2;

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
typedef CGAL::Arr_accessor<Arrangement_2>         Arr_accessor;
typedef Arrangement_2::Vertex_handle              Vertex_handle;
typedef Arrangement_2::Halfedge_handle            Arr_halfedge;

void  my_failure_function (const char *type,
                           const char *expression,
                           const char *file,
                           int line,
                           const char *explanation) 
{
  std::cout<<explanation<<"\n";
}
CGAL::Timer t_boundary;
CGAL::Timer t_do_x;

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

    bool is_back_face = is_back_facet(f);
    //std::cout<<"bfs!!!\n";
    Arrangement_2* arr = new Arrangement_2;
    xy_surfaces.push_back(arr);
    m_queue.push(f);
    while(!m_queue.empty())
    {
      f = m_queue.front();
      m_queue.pop();

      // add f to current xy_surface
      bool _was_added = add_facet_to_arr(*arr, f, is_back_face);
      if(!_was_added)
        continue;
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
            if(are_xy_mon(curr))
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
  bool are_xy_mon(Halfedge_const_handle he)
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

    CGAL_assertion(!is_vertical(he->face()) &&
                   !is_vertical(he->opposite()->face()));
    
    const Point_2& p1 = project_xy(he->next()->vertex()->point());
    const Point_2& p2 = project_xy(he->opposite()->next()->vertex()->point());

    Kernel k;

    CGAL_assertion(m_edges_map[he] != Arr_halfedge());
    const Line_2 & l = m_edges_map[he]->curve().line();

    CGAL::Oriented_side os1 = k.oriented_side_2_object()(l, p1);
    CGAL::Oriented_side os2 = k.oriented_side_2_object()(l, p2);

    CGAL_assertion(os1 != CGAL::ON_ORIENTED_BOUNDARY &&
                   os2 != CGAL::ON_ORIENTED_BOUNDARY);

    return (os1 != os2);
  }

  // when we add a new facet to the xy-surface, we have to check that the 
  // projected edges of the faces doesn't intersect the projected boundary
  // of the xy-surface (just cheking adjacent facet is not good enough)
  // Halfedge_const_handle he : shared halfegde between the current 
  // xy-monotone surface and  facet f, such that, he->opposite()->face() == f
  bool does_projected_facet_intersect_boundary(Facet_const_iterator f,
                                               bool is_back_face)
  {
    //std::cout<<"does_projected_facet_intersect_boundary ...\n";
    t_boundary.start();
    // find that shared halfedge
    Halfedge_const_handle he;
    Arr_halfedge arr_he;
    Halfedge_around_facet_const_circulator haf = f->facet_begin();
    Halfedge_around_facet_const_circulator curr = haf;
    do
    {
      if((arr_he = m_edges_map[curr->opposite()]) != Arr_halfedge())
      {
        he = curr->opposite();
        break;
      }
      curr++;
    }while(curr != haf);

    if(arr_he == Arr_halfedge())
    {
      // when we add he first facet of the xy-monotone surface, thare aren't
      // any arrangement edges we dont need to check for intersection 
      // with the boundary.
      t_boundary.stop();
      return false;
    }

    if(is_back_face)
      arr_he = arr_he->twin();

    CGAL_assertion(he->opposite()->face() == f);
    CGAL_assertion(! arr_he->twin()->face()->is_unbounded());

    Arr_halfedge curr_he;
    
    Halfedge_const_handle shared_he = he->opposite();
    Halfedge_const_handle curr_facet_he;
    
    Kernel k;
  
    /*int c =0;
    for(curr_he = arr_he->next(); curr_he!=arr_he; curr_he=curr_he->next(), c++)
      ;
    std::cout<<"There are : " << c << " edges at the boundary\n";*/
    for(curr_he = arr_he->next(); curr_he!=arr_he; curr_he=curr_he->next())
    {
      for(curr_facet_he = shared_he->next(); curr_facet_he!=shared_he; curr_facet_he=curr_facet_he->next())
      {
        if(m_edges_map[curr_facet_he->opposite()] != Arr_halfedge())
          continue;
        
        bool _cont = false;
        Halfedge_const_handle hh = curr_facet_he;
        do
        {
          Arr_halfedge temp_he = m_edges_map[hh];
          if(temp_he == Arr_halfedge())
          {
            temp_he = m_edges_map[hh->opposite()];
            if(temp_he != Arr_halfedge())
              temp_he = temp_he->twin();
            else 
            {
              hh = hh->next()->opposite();
              continue;
            }
          }
          if(temp_he->target() == curr_he->source() ||
             temp_he->target() == curr_he->target())
          {
            //std::cout<<"_cont(1) !!!\n";
            _cont = true;
            break;
          }
          hh = hh->next()->opposite();
        }while(hh != curr_facet_he);
         
        if(_cont)
          continue;

        //_cont = false;
        Halfedge_const_handle opp_curr_facet_he = curr_facet_he->opposite();
        hh = opp_curr_facet_he;
        do
        {
          Arr_halfedge temp_he = m_edges_map[hh];
          if(temp_he == Arr_halfedge())
          {
            temp_he = m_edges_map[hh->opposite()];
            if(temp_he != Arr_halfedge())
              temp_he = temp_he->twin();
            else
            {
              hh = hh->next()->opposite();
              continue;
            }
          }
          if(temp_he->target() == curr_he->source() ||
            temp_he->target() == curr_he->target())
          {
            //std::cout<<"_cont(1) !!!\n";
            _cont = true;
            break;
          }
          hh = hh->next()->opposite();
        }while(hh != opp_curr_facet_he);

        if(_cont)
          continue;
            
        Segment_2 seg(project_xy(curr_facet_he->opposite()->vertex()->point()),
                      project_xy(curr_facet_he->vertex()->point()));
        t_do_x.start();
        bool do_x = k.do_intersect_2_object()(static_cast<Segment_2>(curr_he->curve()), seg);
        t_do_x.stop();
        if(do_x)
        {
          /*Halfedge_around_facet_const_circulator haf = f->facet_begin();
          Halfedge_around_facet_const_circulator curr = haf;
          do
          {
            std::cout<< curr->opposite()->vertex()->point() <<" --> " << curr->vertex()->point()<<"\n";
            curr++;
          }while(curr != haf);

          std::ofstream file("boundary.arr");
          if(!file.is_open())
          {
            std::cout<<"failed to open file!!!\n";
            abort();
          }
          Arr_halfedge hh = arr_he;
          do
          {
            file << hh->curve() << std::endl;
            hh = hh->next();
          }while(hh!=arr_he);
          std::cout<<"bug: " << curr_facet_he->opposite()->vertex()->point() <<" --> " << curr_facet_he->vertex()->point()<<"\n";
          std::cout<<"bug2: " << curr_he->curve() << "\n";
          std::cout<<"is_back_facet? " << is_back_facet(f)<<"\n";
          t_boundary.stop();
          std::cout<<"return (true)\n";*/
          t_boundary.stop();
          return true;
        }
      }
    }

     
    t_boundary.stop();
    return false;
  }
       
  /*Iso_rectangle_2 construct_iso_rectangle_from_facet(Halfedge_const_handle he)
  {
    Kernel k;

    Halfedge_const_handle curr_facet_he = he;
    Point_2 left, right, bottom, top;

    left = right = bottom = top = curr_facet_he->target()->point();
    curr_facet_he = curr_facet_he->next();
    do
    {
      const Point_2& pt = curr_facet_he->target()->point();
      if(k.compare_x_2_object(pt, left) == CGAL::SMALLER)
        left = pt;
      else
        if(k.compare_x_2_object(pt, right) == CGAL::LARGER)
          right = pt;

      if(k.compare_y_2_object(pt, bottom) == CGAL::SMALLER)
        bottom = pt;
      else
        if(k.compare_y_2_object(pt, top) == CGAL::LARGER)
          top = pt;

      curr_facet_he = curr_facet_he->next();
    }while(curr_facet_he != he);
   

    return (Iso_rectangle_2 rect(left, right, bottom, top));
  }*/
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

  bool add_facet_to_arr(Arrangement_2& arr, Facet_const_iterator f, bool is_back)
  {
    //std::cout<<"add_facet_to_arr!!!\n";
    if(does_projected_facet_intersect_boundary(f, is_back))
    {
      //std::cout<<"found intersection of current facet with boundary!!!\n";
      m_status_map[f] = false;
      return false;
    }
    //static int i=1;
    Halfedge_around_facet_const_circulator haf = f->facet_begin();
    Halfedge_around_facet_const_circulator curr = haf;
    do
    {
      std::cout<< curr->opposite()->vertex()->point() <<" --> " << curr->vertex()->point()<<"\n";
      curr++;
    }while(curr != haf);
    std::cout<<"*******************************\n";

    haf = f->facet_begin();
    curr = haf;

    do
    {
      Halfedge_const_handle he = curr;
      
      /*std::cout<<"current polyhedron edge: " << he->opposite()->vertex()->point() <<" --> " << curr->vertex()->point()<<"\n";
      std::cout<<"is_back? : " << is_back << "\n";*/
      // look for that polyhedron's halfedge at the edges map
      Arr_halfedge arr_he = m_edges_map[he->opposite()];
      CGAL_assertion(m_edges_map[he] == Arr_halfedge());
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

        //CGAL::Comparison_result res = k.compare_xy_2_object()(s, t);
        CGAL::Comparison_result res = cv.is_directed_right() ? CGAL::SMALLER : CGAL::LARGER;

        if(!next_exist && !prev_exist)
        {          
          Arr_accessor arr_accessor(arr);
          Vertex_handle v1 = arr_accessor.create_vertex(s);
          Vertex_handle v2 = arr_accessor.create_vertex(t);
          prev_he = 
            arr_accessor.insert_in_face_interior_ex(cv,
                                                    arr.unbounded_face(),
                                                    v1, v2, res);

          //std::cout<<"("<<i++<<")"<<"insert_in_face_interior: " << cv<<"\n";
          //prev_he = arr.insert_in_face_interior(cv, arr.unbounded_face());
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
                //std::cout<<"("<<i++<<")"<<"insert_from_right_vertex: " << cv<<"\n";
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
              /*std::cout<<"("<<i++<<")"<<"insert_at_vertices: " << cv<<"\n";
              std::cout<<"prev_he: " <<prev_he->source()->point()<<" --> "<< prev_he->target()->point()<<"\n";
              std::cout<<"next_he: " <<next_he->source()->point()<<" --> "<< next_he->target()->point()<<"\n";*/
              //TODO: can be done much more efficient using the Arr_accessor !!
              bool new_face = false;
              Arr_accessor arr_accessor(arr);
              if(is_back && !is_opp_face_get_closed(he)  ||  !is_back && is_opp_face_get_closed(he))
                prev_he = arr_accessor.insert_at_vertices_ex(cv, prev_he, next_he, res, new_face);
              /*if(is_back && is_curr_facet_get_closed(he)  ||  !is_back && !is_curr_facet_get_closed(he))
                prev_he = arr_accessor.insert_at_vertices_ex(cv, prev_he, next_he, res, new_face);*/
              else
              {
                if(res == CGAL::SMALLER)
                  prev_he = arr_accessor.insert_at_vertices_ex(cv, next_he, prev_he, CGAL::LARGER, new_face);
                else
                  prev_he = arr_accessor.insert_at_vertices_ex(cv, next_he, prev_he, CGAL::SMALLER, new_face);
              }
              if(prev_he->direction() != res)
                prev_he = prev_he->twin();
              /*prev_he = arr.insert_at_vertices(cv, prev_he, next_he);
              if(prev_he->direction() != res)
                prev_he = prev_he->twin();*/
            }
        //map 'he' of the polyhderon to 'prev_he'

            std::ofstream arr_file("temp.arr");
            if(!arr_file.is_open())
            {
              std::cout<<"cannot open file !!!\n";
              exit(-1);
            }
            arr_file << arr;
            if(!CGAL::is_valid(arr))
              abort();
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

    return true;
  }

  bool is_opp_face_get_closed(Halfedge_const_handle he)
  {
    Halfedge_const_handle first = he->opposite();
    for(first = first->next(); first != he->opposite(); first = first->next())
    {
      if(m_edges_map[first->opposite()] == Arr_halfedge())
        return false;
    }
    return true;
  }

  bool is_curr_facet_get_closed(Halfedge_const_handle he)
  {
    Halfedge_const_handle first = he->next();
    for(; first != he; first = first->next())
    {
      if(m_edges_map[first] == Arr_halfedge() || 
         m_edges_map[first->opposite()] == Arr_halfedge())
        return false;
    }
    return true;

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
  CGAL::set_warning_handler(my_failure_function);
  //CGAL::set_error_behaviour (CGAL::CONTINUE);

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
  std::cout<<"Intersecting with boundary took : " <<t_boundary.time()<<" seconds\n";
  std::cout<<"do intersect took : " <<t_do_x.time()<<" seconds\n";
  


  std::cout<<"There are: " << xy_surfs.size() << " xy-surfaces\n";
  unsigned int i=1;
  for(Xy_surf_itr itr = xy_surfs.begin(); itr != xy_surfs.end(); ++itr, ++i)
  {
    std::cout<<"#" <<i<<" surface : "<<std::endl;
    std::cout<<"  |V| = "<<(*itr)->number_of_vertices()<<std::endl;
    std::cout<<"  |E| = "<<(*itr)->number_of_edges()<<std::endl;
    std::cout<<"  |F| = "<<(*itr)->number_of_faces()<<std::endl;

    std::cout<<"is_valid? " << CGAL::is_valid(**itr)<<std::endl;
    if(i<10)
    {
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

  }

  return (0);
}
