#ifndef CGAL_MESH_3_TET_SOUP_TO_C3T3_H
#define CGAL_MESH_3_TET_SOUP_TO_C3T3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/IO/File_medit.h>

#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

// The purpose of this file is to rebuild manually a triangulation from a c3t3
// output file (because generating a new triangulation every time is expensive)

namespace CGAL
{
template<class Tr>
void build_vertices(Tr& tr,
                    const std::vector<typename Tr::Point>& points,
                    std::vector<typename Tr::Vertex_handle>& vertex_handle_vector)
{
  typedef typename Tr::Vertex_handle            Vertex_handle;

  vertex_handle_vector[0] = tr.tds().create_vertex(); // creates the infinite vertex
  tr.infinite_vertex() = vertex_handle_vector[0];

  // build vertices
  for(std::size_t i=0; i<points.size(); ++i)
  {
    Vertex_handle vh = tr.tds().create_vertex();
    vertex_handle_vector[i+1] = vh;
    vh->set_point(points[i]);
  }
}

template<class Tr>
void add_facet_to_incident_cells_map(const typename Tr::Cell_handle c, int i,
                                     std::map<std::set<typename Tr::Vertex_handle>,
                                     std::vector<std::pair<typename Tr::Cell_handle,
                                     int> > >& incident_cells_map)
{
  typedef typename Tr::Vertex_handle                            Vertex_handle;
  typedef typename Tr::Cell_handle                              Cell_handle;
  typedef std::set<Vertex_handle>                               Facet;
  typedef std::pair<Cell_handle, int>                           Incident_cell;
  typedef std::map<Facet, std::vector<Incident_cell> >          Incident_cells_map;

  // the opposite vertex of f in c is i
  Facet f;
  f.insert(c->vertex((i + 1) % 4));
  f.insert(c->vertex((i + 2) % 4));
  f.insert(c->vertex((i + 3) % 4));
  CGAL_precondition(f.size() == 3);

  Incident_cell e = std::make_pair(c, i);
  std::vector<Incident_cell> vec;
  vec.push_back(e);
  std::pair<typename Incident_cells_map::iterator, bool> is_insert_successful =
      incident_cells_map.insert(std::make_pair(f, vec));
  if(!is_insert_successful.second) // the entry already exists in the map
  {
    // a facet must have exactly two incident cells
    CGAL_assertion(is_insert_successful.first->second.size() == 1);
    is_insert_successful.first->second.push_back(e);
  }
}

template<class Tr>
void build_finite_cells(Tr& tr,
                        const std::vector<boost::array<int,5> >& finite_cells,
                        const std::vector<typename Tr::Vertex_handle>& vertex_handle_vector,
                        std::map<std::set<typename Tr::Vertex_handle>,
                        std::vector<std::pair<typename Tr::Cell_handle,
                        int> > >& incident_cells_map,
                        const std::map<boost::array<int,3>, int>& border_facets)
{
  typedef boost::array<int, 5>              Tet_with_ref; // 4 ids + 1 reference

  typedef typename Tr::Vertex_handle                            Vertex_handle;
  typedef typename Tr::Cell_handle                              Cell_handle;

  // build the finite cells
  for(std::size_t i=0; i<finite_cells.size(); ++i)
  {
    const Tet_with_ref& tet = finite_cells[i];
    boost::array<Vertex_handle, 4> vs;

    for(int j=0; j<4; ++j)
    {
      CGAL_precondition(static_cast<std::size_t>(tet[j]) < tr.number_of_vertices() &&
                        tet[j] >= 0);
      vs[j] = vertex_handle_vector[tet[j] + 1];
      CGAL_postcondition(vs[j] != Vertex_handle());
    }

    // this assertion also tests for degeneracy
    CGAL_assertion(CGAL::orientation(vs[0]->point().point(), vs[1]->point().point(),
        vs[2]->point().point(), vs[3]->point().point()) == POSITIVE);

    Cell_handle c = tr.tds().create_cell(vs[0], vs[1], vs[2], vs[3]);
    c->info() = tet[4]; // the reference encodes the interior/exterior info

    // the reference must either be -1 for finite interior or > 0 for subdomains
    CGAL_precondition(tet[4] != 0);
    // assign cells to vertices
    for(int j=0; j<4; ++j)
    {
      if(vs[j]->cell() == Cell_handle())
        vs[j]->set_cell(c);
    }

    // build the map used for adjacency later
    for(int j=0; j<4; ++j)
    {
      add_facet_to_incident_cells_map<Tr>(c, j, incident_cells_map);
      if(border_facets.size() !=0)
      {
        boost::array<int,3> facet;
        facet[0]=tet[(j+ 1) % 4];
        facet[1]=tet[(j+ 2) % 4];
        facet[2]=tet[(j+ 3) % 4];
        //find the circular permutation that puts the smallest index in the first place.
        int n0 = (std::min)((std::min)(facet[0],facet[1]), facet[2]);
        int k=0;
        boost::array<int,3> f;
        do
        {
          f[0]=facet[(0+k)%3];
          f[1]=facet[(1+k)%3];
          f[2]=facet[(2+k)%3];
          ++k;
        }while(f[0] != n0);
        if(border_facets.find(f) != border_facets.end())
          c->set_surface_patch_index(j, border_facets.at(f));
        else
        {
          int temp = f[2];
          f[2]=f[1];
          f[1]=temp;
          if(border_facets.find(f) != border_facets.end())
            c->set_surface_patch_index(j, border_facets.at(f));
          else
            c->set_surface_patch_index(j, 0);
        }
      }
    }
  }
}

template<class Tr>
void add_infinite_facets_to_incident_cells_map(typename Tr::Cell_handle c,
                                               int inf_vert_pos,
                                               std::map<std::set<typename Tr::Vertex_handle>,
                                               std::vector<std::pair<typename Tr::Cell_handle,
                                               int> > >& incident_cells_map)
{
  int l = (inf_vert_pos + 1) % 4;
  add_facet_to_incident_cells_map<Tr>(c, l, incident_cells_map);
  l = (inf_vert_pos + 2) % 4;
  add_facet_to_incident_cells_map<Tr>(c, l, incident_cells_map);
  l = (inf_vert_pos + 3) % 4;
  add_facet_to_incident_cells_map<Tr>(c, l, incident_cells_map);
}

template<class Tr>
void build_infinite_cells(Tr& tr,
                          const std::vector<boost::array<int,4> >& infinite_cells,
                          const std::vector<typename Tr::Vertex_handle>& vertex_handle_vector,
                          std::map<std::set<typename Tr::Vertex_handle>,
                          std::vector<std::pair<typename Tr::Cell_handle,
                          int> > >& incident_cells_map)
{
  typedef boost::array<int, 4>                            Tet; // 4 ids

  typedef typename Tr::Vertex_handle                      Vertex_handle;
  typedef typename Tr::Cell_handle                        Cell_handle;
  typedef std::set<Vertex_handle>                         Facet;
  typedef std::pair<Cell_handle, int>                     Incident_cell;
  typedef std::map<Facet, std::vector<Incident_cell> >    Incident_cells_map;

  // build the infinite cells if provided
  for(std::size_t i=0; i<infinite_cells.size(); ++i)
  {
    const Tet& tet = infinite_cells[i];

    int inf_pos = -1; // position of the infinite vertex
    boost::array<Vertex_handle, 4> vs;
    for(int j=0; j<4; ++j)
    {
      if(tet[j] == -1)
      {
        vs[j] = tr.infinite_vertex();
        inf_pos = j;
      }
      else
        vs[j] = vertex_handle_vector[tet[j] + 1];
    }
    CGAL_precondition(inf_pos != -1);

    Cell_handle c = tr.tds().create_cell(vs[0], vs[1], vs[2], vs[3]);

    // could simply be 'if(!i)', but this is clearer
    if(tr.infinite_vertex()->cell() == Cell_handle())
      tr.infinite_vertex()->set_cell(c);

    // add the cell to the incident cells map
    Facet f;
    for(int j=1; j<=3; ++j)
      f.insert(vertex_handle_vector[tet[(inf_pos + j) % 4] + 1]);

    // we're adding infinite cells, and have already inserted the finite cells
    // in the map, so there has to be an entry in the map for this face already
    typename Incident_cells_map::iterator it = incident_cells_map.find(f);
    CGAL_assertion(it != incident_cells_map.end());
    (it->second).push_back(std::make_pair(c, inf_pos));

    // the three infinite facets
    add_infinite_facets_to_incident_cells_map<Tr>(c, inf_pos, incident_cells_map);
  }

  // manually build the infinite cells if they were not provided
  if(infinite_cells.empty())
  {
    // check the incident cells map for facets who only have one incident cell
    // and build the infinite cell on the opposite side
    typename Incident_cells_map::iterator it = incident_cells_map.begin();
    typename Incident_cells_map::iterator end = incident_cells_map.end();
    for(; it != end; ++it)
    {
      if(it->second.size() == 2) // facet already has both its incident cells
        continue;
      CGAL_assertion(it->second.size() == 1);

      Cell_handle c = it->second[0].first;
      int i = it->second[0].second;

      Cell_handle opp_c;
      // the infinite cell that we are creating needs to be well oriented...
      int inf_vert_position_in_opp_c = 0;
      if(i == 0 || i == 2)
        opp_c = tr.tds().create_cell(tr.infinite_vertex(),
                                     c->vertex((i+1)%4),
                                     c->vertex((i+2)%4),
                                     c->vertex((i+3)%4));
      else
        opp_c = tr.tds().create_cell(tr.infinite_vertex(),
                                     c->vertex((i+1)%4),
                                     c->vertex((i+3)%4),
                                     c->vertex((i+2)%4));

      // set the infinite_vertex's incident cell
      if(tr.infinite_vertex()->cell() == Cell_handle())
        tr.infinite_vertex()->set_cell(opp_c);

      // add the facets to the incident cells map

      // the only finite facet
      it->second.push_back(std::make_pair(opp_c, inf_vert_position_in_opp_c));
      CGAL_assertion(it->second.size() == 2);
    }
  }
}

template<class Tr>
bool assign_neighbors(Tr& tr,
                      const std::map<std::set<typename Tr::Vertex_handle>,
                      std::vector<std::pair<typename Tr::Cell_handle,
                      int> > >& incident_cells_map)
{
  typedef typename Tr::Vertex_handle                              Vertex_handle;
  typedef typename Tr::Cell_handle                                Cell_handle;
  typedef std::set<Vertex_handle>                                 Facet;
  typedef std::pair<Cell_handle, int>                             Incident_cell;
  typedef std::map<Facet, std::vector<Incident_cell> >            Incident_cells_map;

  typename Incident_cells_map::const_iterator icit = incident_cells_map.begin();
  for(; icit!=incident_cells_map.end(); ++icit)
  {
    const std::vector<Incident_cell>& adjacent_cells = icit->second;
    if(adjacent_cells.size() != 2)
      return false;

    Cell_handle c0 = adjacent_cells[0].first;
    int i0 = adjacent_cells[0].second;
    Cell_handle c1 = adjacent_cells[1].first;
    int i1 = adjacent_cells[1].second;

    tr.tds().set_adjacency(c0, i0, c1, i1);
  }
  return true;
}

template<class Tr, bool c3t3_loader_failed>
bool build_triangulation(Tr& tr,
                         const std::vector<typename Tr::Point>& points,
                         const std::vector<boost::array<int,5> >& finite_cells,
                         const std::vector<boost::array<int,4> >& infinite_cells,
                         const std::map<boost::array<int,3>, int>& border_facets)
{
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;
  typedef std::set<Vertex_handle>               Facet;

  // associate to a face the two (at most) incident tets and the id of the face in the cell
  typedef std::pair<Cell_handle, int>                   Incident_cell;
  typedef std::map<Facet, std::vector<Incident_cell> >  Incident_cells_map;

  Incident_cells_map incident_cells_map;
  std::vector<Vertex_handle> vertex_handle_vector(points.size() + 1); // id to vertex_handle

  CGAL_precondition(!points.empty() && !finite_cells.empty());

  tr.tds().clear(); // not tr.clear() since it calls tr.init() which we don't want

  build_vertices<Tr>(tr, points, vertex_handle_vector);
  build_finite_cells<Tr>(tr, finite_cells, vertex_handle_vector, incident_cells_map, border_facets);
  build_infinite_cells<Tr>(tr, infinite_cells, vertex_handle_vector, incident_cells_map);
  tr.tds().set_dimension(3);
  if(!assign_neighbors<Tr>(tr, incident_cells_map))
    return false;

  std::cout << "built triangulation : " << std::endl;
  std::cout << tr.number_of_vertices() << " vertices" << std::endl;
  std::cout << tr.number_of_cells() << " cells" << std::endl;

  BOOST_FOREACH(Vertex_handle vh, vertex_handle_vector)
  {
    vh->set_dimension(3);
  }
  if(c3t3_loader_failed)
  {
    return true;
  }
  else
    return tr.is_valid(true);
}

template<class Tr, bool c3t3_loader_failed>
bool build_triangulation_from_file(std::istream& is,
                                   Tr& tr)
{
  typedef typename Tr::Point                                  Point_3;

  typedef boost::array<int, 3> Facet; // 3 = id
  typedef boost::array<int, 4> Tet; // 4 = id
  typedef boost::array<int, 5> Tet_with_ref; // first 4 = id, fifth = reference

  std::size_t finite_exterior_cells_counter = 0;
  std::vector<Tet_with_ref> finite_cells;
  std::vector<Tet> infinite_cells;
  std::vector<Point_3> points;
  std::map<Facet, int> border_facets;

  // grab the vertices
  int dim;
  int nv, nf, ntet, ref;
  std::string word;

  is >> word >> dim; // MeshVersionFormatted 1
  is >> word >> dim; // Dimension 3

  CGAL_assertion(dim == 3);

  while(is >> word && word != "End")
  {
    if(word == "Vertices")
    {
      is >> nv;
      for(int i=0; i<nv; ++i)
      {
        double x,y,z;
        is >> x >> y >> z>>ref;
        points.push_back(Point_3(x,y,z));
      }
    }

    if(word == "Triangles")
    {
      is >> nf;
      for(int i=0; i<nf; ++i)
      {
        int n1, n2, n3, surface_patch_id;
        is >> n1 >> n2 >> n3 >> surface_patch_id;
        // no use for boundary facets for now...
        Facet facet;
        facet[0] =n1 -1;
        facet[1] =n2 -1;
        facet[2] =n3 -1;
        //find the circular permutation that puts the smallest index in the first place.
        int n0 = (std::min)((std::min)(facet[0],facet[1]), facet[2]);
        int k=0;
        Facet f;
        do
        {
          f[0]=facet[(0+k)%3];
          f[1]=facet[(1+k)%3];
          f[2]=facet[(2+k)%3];
          ++k;
        }while(f[0] != n0);
        border_facets.insert(std::make_pair(f, surface_patch_id));
      }
    }

    if(word == "Tetrahedra")
    {
      is >> ntet;
      for(int i=0; i<ntet; ++i)
      {
        int n0, n1, n2, n3, reference;
        is >> n0 >> n1 >> n2 >> n3 >> reference;
        Tet_with_ref t;
        t[0] = n0 - 1;
        t[1] = n1 - 1;
        t[2] = n2 - 1;
        t[3] = n3 - 1;
        t[4] = reference;
        finite_cells.push_back(t);
      }
    }
  }



  if(!finite_exterior_cells_counter)
  {
    // The finite interior cells MUST be the convex hulls of the point
    std::cerr << "WARNING: no finite exterior cell provided..." << std::endl;
  }

  CGAL_precondition(!finite_cells.empty());

  bool is_well_built = build_triangulation<Tr, c3t3_loader_failed>(tr, points, finite_cells, infinite_cells, border_facets);
  return is_well_built;
}

}  // namespace CGAL

#endif // CGAL_MESH_3_TET_SOUP_TO_C3T3_H
