#ifndef PROTECTED_SETS_H
#define PROTECTED_SETS_H

#include <CGAL/Cartesian_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include <algorithm>
#include <list>
#include <chrono>
#include <random>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::FT FT;
typedef K::Point_d Point_d;
typedef K::Vector_d Vector_d;
typedef K::Oriented_side_d Oriented_side_d;
typedef K::Has_on_positive_side_d Has_on_positive_side_d;
typedef K::Sphere_d Sphere_d;
typedef K::Hyperplane_d Hyperplane_d;

typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;
typedef Delaunay_triangulation::Facet Facet;
typedef Delaunay_triangulation::Vertex_handle Delaunay_vertex;
typedef Delaunay_triangulation::Full_cell_handle Full_cell_handle;

typedef std::vector<Point_d> Point_Vector;
typedef CGAL::Euclidean_distance<K> Euclidean_distance;

FT  _sfty = pow(10,-14);

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// AUXILLARY FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Insert a point in Delaunay triangulation. If you are working in a flat torus, the procedure adds all the 3^d copies in adjacent cubes as well 
 *  
 *  W is the initial point vector
 *  chosen_landmark is the index of the chosen point in W
 *  landmarks_ind is the vector of indices of already chosen points in W
 *  delaunay is the Delaunay triangulation
 *  landmark_count is the current number of chosen vertices
 *  torus is true iff you are working on a flat torus [-1,1]^d
 *  OUT: Vertex handle to the newly inserted point
 */
Delaunay_vertex insert_delaunay_landmark_with_copies(Point_Vector& W, int chosen_landmark, std::vector<int>& landmarks_ind, Delaunay_triangulation& delaunay, int& landmark_count, bool torus)
{
  if (!torus)
    {
      Delaunay_vertex v =delaunay.insert(W[chosen_landmark]);
      landmarks_ind.push_back(chosen_landmark);
      landmark_count++;
      return v;
    }
  else
    {
      int D = W[0].size();
      int nb_cells = pow(3, D);
      Delaunay_vertex v;
      for (int i = 0; i < nb_cells; ++i)
        {
          std::vector<FT> point;
          int cell_i = i;
          for (int l = 0; l < D; ++l)
            {
              point.push_back(W[chosen_landmark][l] + 2.0*(cell_i%3-1));
              cell_i /= 3;
            }
          v = delaunay.insert(point);
        }
      landmarks_ind.push_back(chosen_landmark);
      landmark_count++;
      return v;
    }
}

/** Small check if the vertex v is in the full cell fc
 */

bool vertex_is_in_full_cell(Delaunay_triangulation::Vertex_handle v, Full_cell_handle fc)
{
  for (auto v_it = fc->vertices_begin(); v_it != fc->vertices_end(); ++v_it)
    if (*v_it == v)
      return true;
  return false;
}

/** Fill chosen point vector from indices with copies if you are working on a flat torus
 *  
 *  IN:  W is the point vector
 *  OUT: landmarks is the output vector
 *  IN:  landmarks_ind is the vector of indices
 *  IN:  torus is true iff you are working on a flat torus [-1,1]^d
 */

void fill_landmarks(Point_Vector& W, Point_Vector& landmarks, std::vector<int>& landmarks_ind, bool torus)
{
  if (!torus)
    for (unsigned j = 0; j < landmarks_ind.size(); ++j)
      landmarks.push_back(W[landmarks_ind[j]]);
  else
    {
      int D = W[0].size();
      int nb_cells = pow(3, D);
      int nbL = landmarks_ind.size();
      // Fill landmarks
      for (int i = 0; i < nb_cells-1; ++i)
        for (int j = 0; j < nbL; ++j)
          {
            int cell_i = i;
            Point_d point;
            for (int l = 0; l < D; ++l)
              {
                point.push_back(W[landmarks_ind[j]][l] + 2.0*(cell_i-1));
                cell_i /= 3;
              }
            landmarks.push_back(point);
          }
    }
}

/** Fill a vector of all simplices in the Delaunay triangulation giving integer indices to vertices
 *
 *  IN: t is the Delaunay triangulation
 *  OUT: full_cells is the output vector
 */

void fill_full_cell_vector(Delaunay_triangulation& t, std::vector<std::vector<int>>& full_cells)
{
  // Store vertex indices in a map
  int ind = 0; //index of a vertex
  std::map<Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
  for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
    if (t.is_infinite(v_it))
      continue;
    else
      index_of_vertex[v_it] = ind++;
  // Write full cells as vectors in full_cells
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    {
      if (t.is_infinite(fc_it))
        continue;
      Point_Vector vertices;
      for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
        vertices.push_back((*fc_v_it)->point());
      Sphere_d cs( vertices.begin(), vertices.end());
      Point_d csc = cs.center();
      bool in_cube = true; 
      for (auto xi = csc.cartesian_begin(); xi != csc.cartesian_end(); ++xi)
        if (*xi > 1.0 || *xi < -1.0)
          {
            in_cube = false; break;
          }
      if (!in_cube)
        continue;
      std::vector<int> cell;
      for (auto v_it = fc_it->vertices_begin(); v_it != fc_it->vertices_end(); ++v_it)
        cell.push_back(index_of_vertex[*v_it]);
      full_cells.push_back(cell);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IS VIOLATED TEST
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Check if a newly created cell is protected from old vertices
 *  
 *  t is the Delaunay triangulation
 *  vertices is the vector containing the point to insert and a facet f in t
 *  v1 is the vertex of t, such that f and v1 form a simplex
 *  v2 is the vertex of t, such that f and v2 form another simplex
 *  delta is the protection constant
 *  power_protection is true iff the delta-power protection is used
 */

bool new_cell_is_violated(Delaunay_triangulation& t, std::vector<Point_d>& vertices, const Delaunay_vertex& v1, const Delaunay_vertex v2, FT delta, bool power_protection)
{
  assert(vertices.size() == vertices[0].size() ||
         vertices.size() == vertices[0].size() + 1); //simplex size = d | d+1
  assert(v1 != v2);
  if (vertices.size() == vertices[0].size() + 1)
  // FINITE CASE
    {
      Sphere_d cs(vertices.begin(), vertices.end());
      Point_d center_cs = cs.center();
      FT r = sqrt(Euclidean_distance().transformed_distance(center_cs, vertices[0]));
      /*
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        if (!t.is_infinite(v_it)) 
          {
            //CGAL::Oriented_side side = Oriented_side_d()(cs, (v_it)->point()); 
            if (std::find(vertices.begin(), vertices.end(), v_it->point()) == vertices.end())
              {
                FT dist2 = Euclidean_distance().transformed_distance(center_cs, (v_it)->point());
                if (!power_protection)
                  if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
                    return true;
                if (power_protection)
                  if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
                    return true;                 
              }
          }
      */
      if (!t.is_infinite(v1))
        {
          FT dist2 = Euclidean_distance().transformed_distance(center_cs, v1->point());
          if (!power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
              return true;
          if (power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
              return true;
        }
      if (!t.is_infinite(v2))
        {
          FT dist2 = Euclidean_distance().transformed_distance(center_cs, v2->point());
          if (!power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
              return true;
          if (power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
              return true;
        }
    }
  else
  // INFINITE CASE
    {
      Delaunay_triangulation::Vertex_iterator v = t.vertices_begin();
      while (t.is_infinite(v) || std::find(vertices.begin(), vertices.end(), v->point()) == vertices.end())
        v++;
      Hyperplane_d facet_plane(vertices.begin(), vertices.end(), v->point(), CGAL::ON_POSITIVE_SIDE);
      Vector_d orth_v = facet_plane.orthogonal_vector();
      /*
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        if (!t.is_infinite(v_it))
          if (std::find(vertices.begin(), vertices.end(), v_it->point()) == vertices.end())
            {
              std::vector<FT> coords;
              Point_d p = v_it->point();
              auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
              for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
                coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
              Point_d p_delta = Point_d(coords);
              bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
              bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);
              if (!p_is_inside && p_delta_is_inside)
                return true;
            }
      */
      if (!t.is_infinite(v1))
        {
          std::vector<FT> coords;
          Point_d p = v1->point();
          auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
          for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
            coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
          Point_d p_delta = Point_d(coords);
          bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
          bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);
          if (!p_is_inside && p_delta_is_inside)
            return true;
        }
      if (!t.is_infinite(v2))
        {
          std::vector<FT> coords;
          Point_d p = v2->point();
          auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
          for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
            coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
          Point_d p_delta = Point_d(coords);
          bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
          bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);
          if (!p_is_inside && p_delta_is_inside)
            return true;
        }
    }
  return false;
}
  
/** Auxillary recursive function to check if the point p violates the protection of the cell c and
 *  if there is a violation of an eventual new cell
 *
 *  p is the point to insert
 *  t is the current triangulation
 *  c is the current cell (simplex)
 *  parent_cell is the parent cell (simplex)
 *  index is the index of the facet between c and parent_cell from parent_cell's point of view
 *  D is the dimension of the triangulation
 *  delta is the protection constant
 *  marked_cells is the vector of all visited cells containing p in their circumscribed ball
 *  power_protection is true iff you are working with delta-power protection
 *
 *  OUT: true iff inserting p hasn't produced any violation so far
 */

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, Full_cell_handle c, Full_cell_handle parent_cell, int index, int D, FT delta, std::vector<Full_cell_handle>& marked_cells, bool power_protection)
{
  Euclidean_distance ed;
  std::vector<Point_d> vertices;
  if (!t.is_infinite(c))
    {
      // if the cell is finite, we look if the protection is violated
      for (auto v_it = c->vertices_begin(); v_it != c->vertices_end(); ++v_it)
        vertices.push_back((*v_it)->point());
      Sphere_d cs( vertices.begin(), vertices.end());
      Point_d center_cs = cs.center();
      FT r = sqrt(ed.transformed_distance(center_cs, vertices[0]));
      FT dist2 = ed.transformed_distance(center_cs, p);
      // if the new point is inside the protection ball of a non conflicting simplex
      if (!power_protection)
        if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
          return true;
      if (power_protection)
        if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
          return true;    
      // if the new point is inside the circumscribing ball : continue violation searching on neighbours
      //if (dist2 < r*r)
      //if (dist2 < (5*r+delta)*(5*r+delta))
      if (dist2 < r*r)
        {
          c->tds_data().mark_visited();
          marked_cells.push_back(c);
          for (int i = 0; i < D+1; ++i)
            {
              Full_cell_handle next_c = c->neighbor(i);
              if (next_c->tds_data().is_clear() &&
                  is_violating_protection(p, t, next_c, c, i, D, delta, marked_cells, power_protection))
                return true;
            }
        }
      // if the new point is outside the protection sphere
      else
        {
          // facet f is on the border of the conflict zone : check protection of simplex {p,f}
          // the new simplex is guaranteed to be finite
          vertices.clear(); vertices.push_back(p);
          for (int i = 0; i < D+1; ++i)
            if (i != index)
              vertices.push_back(parent_cell->vertex(i)->point());
          Delaunay_vertex vertex_to_check = t.infinite_vertex();
          for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
            if (!vertex_is_in_full_cell(*vh_it, parent_cell))
              {
                vertex_to_check = *vh_it; break;
              }
          if (new_cell_is_violated(t, vertices, vertex_to_check, parent_cell->vertex(index), delta, power_protection)) 
          //if (new_cell_is_violated(t, vertices, vertex_to_check->point(), delta)) 
            return true;            
        }
    } 
  else
    {
      // Inside of the convex hull is + side. Outside is - side.
      for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
        if (!t.is_infinite(*vh_it))
          vertices.push_back((*vh_it)->point());
      Delaunay_triangulation::Vertex_iterator v_it = t.vertices_begin();
      while (t.is_infinite(v_it) || vertex_is_in_full_cell(v_it, c))
        v_it++;
      Hyperplane_d facet_plane(vertices.begin(), vertices.end(), v_it->point(), CGAL::ON_POSITIVE_SIDE);
      //CGAL::Oriented_side outside = Oriented_side_d()(facet_plane, v_it->point());
      Vector_d orth_v = facet_plane.orthogonal_vector();
      std::vector<FT> coords;
      auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
      for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
        coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
      Point_d p_delta = Point_d(coords);
      bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p) && (Oriented_side_d()(facet_plane, p) != CGAL::ZERO);
      bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);

      if (!p_is_inside && p_delta_is_inside)
        return true;
      //if the cell is infinite we look at the neighbours regardless
      if (p_is_inside)
        {
          c->tds_data().mark_visited();
          marked_cells.push_back(c);    
          for (int i = 0; i < D+1; ++i)
            {
              Full_cell_handle next_c = c->neighbor(i);
              if (next_c->tds_data().is_clear() &&
                  is_violating_protection(p, t, next_c, c, i, D, delta, marked_cells, power_protection))
                return true;
            }
        }
      else
        {
          // facet f is on the border of the conflict zone : check protection of simplex {p,f}
          // the new simplex is finite if the parent cell is finite
          vertices.clear(); vertices.push_back(p);
          for (int i = 0; i < D+1; ++i)
            if (i != index)
              if (!t.is_infinite(parent_cell->vertex(i)))
                vertices.push_back(parent_cell->vertex(i)->point());
          Delaunay_vertex vertex_to_check = t.infinite_vertex();
          for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
            if (!vertex_is_in_full_cell(*vh_it, parent_cell))
              {
                vertex_to_check = *vh_it; break;
              }
          if (new_cell_is_violated(t, vertices, vertex_to_check, parent_cell->vertex(index), delta, power_protection)) 
          //if (new_cell_is_violated(t, vertices, vertex_to_check->point(), delta)) 
            return true;
        }
    }
  //c->tds_data().clear_visited();
  //marked_cells.pop_back();
  return false;
}

/** Checks if inserting the point p in t will make conflicts
 *
 *  p is the point to insert
 *  t is the current triangulation
 *  D is the dimension of triangulation
 *  delta is the protection constant
 *  power_protection is true iff you are working with delta-power protection
 *  OUT: true iff inserting p produces a violation of delta-protection.
 */

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, int D, FT delta, bool power_protection)
{
  Euclidean_distance ed;
  Delaunay_triangulation::Vertex_handle v;
  Delaunay_triangulation::Face f(t.current_dimension()); 
  Delaunay_triangulation::Facet ft; 
  Delaunay_triangulation::Full_cell_handle c; 
  Delaunay_triangulation::Locate_type lt;
  std::vector<Full_cell_handle> marked_cells;
  c = t.locate(p, lt, f, ft, v);
  bool violation_existing_cells = is_violating_protection(p, t, c, c, 0, D, delta, marked_cells, power_protection);
  for (Full_cell_handle fc : marked_cells)
    fc->tds_data().clear();
  return violation_existing_cells;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!! THE INTERFACE FOR LANDMARK CHOICE IS BELOW !!!!!!!!!!//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// LANDMARK CHOICE PROCEDURE
///////////////////////////////////////////////////////////////////////

/** Procedure to compute a maximal protected subset from a point cloud. All OUTs should be empty at call.
 *
 *  IN:  W is the initial point cloud having type Epick_d<Dynamic_dimension_tag>::Point_d
 *  IN:  nbP is the size of W
 *  OUT: landmarks is the output vector for the points
 *  OUT: landmarks_ind is the output vector for the indices of the selected points in W
 *  IN:  delta is the constant of protection
 *  OUT: full_cells is the output vector of the simplices in the final Delaunay triangulation
 *  IN:  torus is true iff you are working on a flat torus [-1,1]^d
 */ 

void landmark_choice_protected_delaunay(Point_Vector& W, int nbP, Point_Vector& landmarks, std::vector<int>& landmarks_ind, FT delta, std::vector<std::vector<int>>& full_cells, bool torus, bool power_protection)
{
  unsigned D = W[0].size();
  Euclidean_distance ed;
  Delaunay_triangulation t(D);
  CGAL::Random rand;
  int landmark_count = 0;
  std::list<int> index_list;
  //      shuffle the list of indexes (via a vector)
  {
    std::vector<int> temp_vector;
    for (int i = 0; i < nbP; ++i)
      temp_vector.push_back(i);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(temp_vector.begin(), temp_vector.end(), std::default_random_engine(seed));
    //CGAL::spatial_sort(temp_vector.begin(), temp_vector.end());
    for (std::vector<int>::iterator it = temp_vector.begin(); it != temp_vector.end(); ++it)
      index_list.push_front(*it);
  }
  if (!torus)
    for (unsigned pos1 = 0; pos1 < D+1; ++pos1)
      {
        std::vector<FT> point;
        for (unsigned i = 0; i < pos1; ++i)
          point.push_back(-1);
        if (pos1 != D)
          point.push_back(1);
        for (unsigned i = pos1+1; i < D; ++i)
          point.push_back(0);
        assert(point.size() == D);
        W[index_list.front()] = Point_d(point);
        insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count, torus);
        index_list.pop_front();
      }
  else if (D == 2)
    {
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2; ++j)
          {
            W[index_list.front()] = Point_d(std::vector<FT>{i*0.5, j*1.0});
            insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count, torus);
            index_list.pop_front();
            W[index_list.front()] = Point_d(std::vector<FT>{0.25+i*0.5, 0.5+j});
            insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count, torus);
            index_list.pop_front();
          }
    }
  else
    std::cout << "No torus starter available for dim>2\n";
  std::list<int>::iterator list_it = index_list.begin();
  while (list_it != index_list.end())
    {
      if (!is_violating_protection(W[*list_it], t, D, delta, power_protection))
        {
          // If no conflicts then insert in every copy of T^3

          insert_delaunay_landmark_with_copies(W, *list_it, landmarks_ind, t, landmark_count, torus);
          index_list.erase(list_it);
          list_it = index_list.begin();
          /*
          // PIECE OF CODE FOR DEBUGGING PURPOSES
          
          Delaunay_vertex inserted_v = insert_delaunay_landmark_with_copies(W, *list_it, landmarks_ind, t, landmark_count);
          if (triangulation_is_protected(t, delta))
            {
              index_list.erase(list_it);
              list_it = index_list.begin();
            }
          else
            { //THAT'S WHERE SOMETHING'S WRONG
              t.remove(inserted_v);
              landmarks_ind.pop_back();
              landmark_count--;
              write_delaunay_mesh(t, W[*list_it], is2d);
              is_violating_protection(W[*list_it], t_old, D, delta); //Called for encore
            }
          */
          //std::cout << "index_list_size() = " << index_list.size() << "\n";
        }
      else
        {
          list_it++;
          //std::cout << "!!!!!WARNING!!!!! A POINT HAS BEEN OMITTED!!!\n";
        }
      //if (list_it != index_list.end())
      //  write_delaunay_mesh(t, W[*list_it], is2d);
    }
  fill_landmarks(W, landmarks, landmarks_ind, torus);
  fill_full_cell_vector(t, full_cells);
  /*
  if (triangulation_is_protected(t, delta))
    std::cout << "Triangulation is ok\n";
  else
    {
      std::cout << "Triangulation is BAD!! T_T しくしく!\n";
    }
  */
  //write_delaunay_mesh(t, W[0], is2d);
  //std::cout << t << std::endl;
}

#endif
