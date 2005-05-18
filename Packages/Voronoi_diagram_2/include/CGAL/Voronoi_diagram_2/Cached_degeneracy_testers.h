#ifndef CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <cstdlib>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================


template<class Edge_tester_t, class Project_t>
class Cached_edge_degeneracy_tester
{
 public:
  typedef Edge_tester_t    Edge_degeneracy_tester;
  typedef Project_t        Data_structure_project;
  typedef typename Edge_degeneracy_tester::Dual_graph Dual_graph;
  typedef typename Edge_degeneracy_tester::Edge Edge;
  typedef typename Edge_degeneracy_tester::Face_handle Face_handle;
  typedef typename Edge_degeneracy_tester::Edge_circulator Edge_circulator;

  typedef typename Edge_degeneracy_tester::All_edges_iterator
  All_edges_iterator;

  typedef typename Edge_degeneracy_tester::Finite_edges_iterator
  Finite_edges_iterator;

 private:
  class Edge_hash_function
    : public Handle_hash_function
  {
  private:
    typedef Handle_hash_function     Base;

  public:
    typedef Base::result_type        result_type;

    template<class Edge>
    result_type operator()(const Edge& e) const
    {
      return (Base::operator()(e.first)) << e.second;
    }
  };

  // true if degenerate, false otherwise
  typedef Unique_hash_map<Edge,bool,Edge_hash_function>  Edge_map;

  Edge opposite(const Edge& e) const {
    int i_mirror =
      Data_structure_project()(dg_).mirror_index(e.first, e.second);
    return Edge( e.first->neighbor(e.second), i_mirror );
  }

 public:
  Cached_edge_degeneracy_tester(const Dual_graph* dual = NULL)
    : dg_(dual), e_tester(dual) {}

  bool operator()(const Edge& e) const {
    if ( emap.is_defined(e) ) { return emap[e]; }
    bool b = e_tester(e);
    emap[e] = b;
    emap[opposite(e)] = b;
    return b;
  }

  bool operator()(const Face_handle& f, int i) const {
    return operator()(Edge(f,i));
  }

  bool operator()(const Edge_circulator& ec) const {
    return operator()(*ec);
  }

  bool operator()(const All_edges_iterator& eit) const {
    return operator()(*eit);
  }

  bool operator()(const Finite_edges_iterator& eit) const {
    return operator()(*eit);
  }

 private:
  const Dual_graph* dg_;
  Edge_degeneracy_tester e_tester;
  mutable Edge_map emap;
};


//=========================================================================
//=========================================================================


template<class Face_degeneracy_t, class Edge_degeneracy_t>
class Cached_face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef Face_degeneracy_t                     Face_degeneracy_tester;
  typedef Edge_degeneracy_t                     Edge_degeneracy_tester;

  typedef typename Face_degeneracy_tester::Dual_graph     Dual_graph;
  typedef typename Face_degeneracy_tester::Vertex_handle  Vertex_handle;

  typedef typename Face_degeneracy_tester::Finite_vertices_iterator
  Finite_vertices_iterator;

 private:
  typedef Unique_hash_map<Vertex_handle,bool>   Vertex_map;

 public:
  Cached_face_degeneracy_tester(const Dual_graph* dual = NULL)
    : f_tester(dual) {}

  Cached_face_degeneracy_tester(const Dual_graph* dual,
				const Edge_degeneracy_tester* e_tester)
    : dg_(dual), f_tester(dual, e_tester) {}

  bool operator()(const Vertex_handle& v) const {
    if ( vmap.is_defined(v) ) { return vmap[v]; }
    bool b = f_tester(v);
    vmap[v] = b;
    return b;
  }
 
  bool operator()(const Finite_vertices_iterator& vit) const {
    return operator()(Vertex_handle(vit));
  }

 private:
  const Dual_graph* dg_;
  Face_degeneracy_tester f_tester;
  mutable Vertex_map vmap;
};


//=========================================================================
//=========================================================================


template<class Tester_handle_t>
struct Handle_to_tester_adaptor
{
  typedef Tester_handle_t                       Tester_handle;
  typedef typename Tester_handle::element_type  Tester;
  typedef typename Tester::result_type          result_type;

  typedef typename Tester::Dual_graph           Dual_graph;

  Handle_to_tester_adaptor() {}
  Handle_to_tester_adaptor(const Dual_graph* dg = NULL)
  {
    Tester rc_tester(dg);
    h_.initialize_with(rc_tester);
  }

  template<class Edge_tester>
  Handle_to_tester_adaptor(const Dual_graph* dg,
			   const Edge_tester& e_tester)
  {
    Tester rc_tester(dg, e_tester);
    h_.initialize_with(rc_tester);
  }

  template<typename argument_type>
  result_type operator()(const argument_type& arg) const {
    return h_.Ptr()->operator()(arg);
  }

  template<typename argument_type_1, typename argument_type_2>
  result_type operator()(const argument_type_1& arg1,
			 const argument_type_2& arg2) const {
    return h_.Ptr()->operator()(arg1, arg2);
  }

 private:
  Tester_handle h_;
};


//=========================================================================
//=========================================================================

template<class Edge_tester_t, class Project_t>
class Ref_counted_edge_degeneracy_tester_base
  : public Cached_edge_degeneracy_tester<Edge_tester_t,Project_t>,
    public Ref_counted_virtual
{
 private:
  typedef Cached_edge_degeneracy_tester<Edge_tester_t,Project_t>  Base;

 public:
  typedef bool                           result_type;
  typedef typename Base::Dual_graph      Dual_graph;

  Ref_counted_edge_degeneracy_tester_base(const Dual_graph* dg = NULL)
    : Base(dg) {}

  ~Ref_counted_edge_degeneracy_tester_base() {}
};

//=========================================================================

template<class Face_tester_t, class Edge_tester_t>
class Ref_counted_face_degeneracy_tester_base
  : public Cached_face_degeneracy_tester<Face_tester_t,Edge_tester_t>,
    public Ref_counted_virtual
{
 private:
  typedef Cached_face_degeneracy_tester<Face_tester_t,Edge_tester_t> Base;

 public:
  typedef bool                         result_type;
  typedef typename Base::Dual_graph    Dual_graph;

  Ref_counted_face_degeneracy_tester_base(const Dual_graph* dg = NULL)
    : Base(dg) {}

  Ref_counted_face_degeneracy_tester_base(const Dual_graph* dg,
					  const Edge_tester_t* e_tester)
    : Base(dg,e_tester) {}

  ~Ref_counted_face_degeneracy_tester_base() {}
};


//=========================================================================
//=========================================================================


template<class Edge_tester_t, class Project_t>
class Ref_counted_edge_degeneracy_tester
  : public Handle_to_tester_adaptor
  <Handle_for_virtual
   <Ref_counted_edge_degeneracy_tester_base<Edge_tester_t,
					    Project_t>
   > >
{
 private:
  typedef
  Ref_counted_edge_degeneracy_tester_base<Edge_tester_t,Project_t>
  Ref_counted_tester_base;

  typedef Handle_for_virtual<Ref_counted_tester_base>
  Ref_counted_tester_base_handle;

  typedef Handle_to_tester_adaptor<Ref_counted_tester_base_handle>
  Base;

 public:
  typedef typename Base::Dual_graph Dual_graph;

  Ref_counted_edge_degeneracy_tester(const Dual_graph* dg = NULL)
    : Base(dg) {}
};

//=========================================================================

template<class Face_tester_t, class Edge_tester_t>
class Ref_counted_face_degeneracy_tester
  : public Handle_to_tester_adaptor
  <Handle_for_virtual
   <Ref_counted_face_degeneracy_tester_base<Face_tester_t,
					    Edge_tester_t>
   > >
{
 private:
  typedef
  Ref_counted_face_degeneracy_tester_base<Face_tester_t,Edge_tester_t>
  Ref_counted_tester_base;

  typedef Handle_for_virtual<Ref_counted_tester_base>
  Ref_counted_tester_base_handle;

  typedef Handle_to_tester_adaptor<Ref_counted_tester_base_handle>
  Base;

 public:
  typedef typename Base::Dual_graph Dual_graph;

  Ref_counted_face_degeneracy_tester(const Dual_graph* dg = NULL)
    : Base(dg) {}

  Ref_counted_face_degeneracy_tester(const Dual_graph* dg,
  				     const Edge_tester_t* e_tester)
    : Base(dg, e_tester) {}
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
