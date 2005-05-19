#ifndef CGAL_DELAUNAY_GRAPH_CONCEPT_H
#define CGAL_DELAUNAY_GRAPH_CONCEPT_H 1

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <iostream>
#include <fstream>


CGAL_BEGIN_NAMESPACE


template<class K>
class Delaunay_graph_concept
{
 public:
  typedef unsigned int   size_type;
  typedef K              Geom_traits;

  typedef typename Geom_traits::Point_2  Site_2;
  typedef typename Geom_traits::Point_2  Point_2;

 private:
  struct Dummy_iterator : public Emptyset_iterator
  {
    typedef unsigned int size_type;

    Dummy_iterator() {}
    Dummy_iterator(const Dummy_iterator&) {}

    template< class T >
    Dummy_iterator& operator=(const T&) { return *this; }

    Dummy_iterator& operator++()        { return *this; }
    Dummy_iterator& operator++(int)     { return *this; }

    Dummy_iterator& operator--()        { return *this; }
    Dummy_iterator& operator--(int)     { return *this; }

    Dummy_iterator& operator*()         { return *this; }
    Dummy_iterator* operator->()        { return this; }

    const Dummy_iterator& operator*() const   { return *this; }
    const Dummy_iterator* operator->() const  { return this; }

    bool operator==(const Dummy_iterator&) const {
      return true;
    }

    bool operator!=(const Dummy_iterator&) const {
      return false;
    }

    bool operator<(const Dummy_iterator& other) const {
      return this < &other;
    }

  };

  //  static Dummy_iterator dummy_iterator;
  static const Dummy_iterator& dummy_it() {
    static Dummy_iterator dummy_it_static;
    return dummy_it_static;
  }

 public:

  //  struct Vertex_circulator {};

  struct Vertex_handle : public Dummy_iterator {
    typedef Vertex_handle&  reference;
    typedef Vertex_handle*  pointer;

    Vertex_handle() {}
    Vertex_handle(const Dummy_iterator&) {}

    Vertex_handle& operator*()         { return *this; }
    Vertex_handle* operator->()        { return this; }

    const Vertex_handle& operator*() const   { return *this; }
    const Vertex_handle* operator->() const  { return this; }

    Site_2 site() const { return Site_2(); }
  };

  struct Face_handle : public Dummy_iterator {
    typedef Face_handle&  reference;
    typedef Face_handle*  pointer;

    Face_handle() {}
    Face_handle(const Dummy_iterator&) {}

    Face_handle& operator*()         { return *this; }
    Face_handle* operator->()        { return this; }

    const Face_handle& operator*() const   { return *this; }
    const Face_handle* operator->() const  { return this; }

    Vertex_handle vertex(int) const { return Vertex_handle(); }

    Face_handle neighbor(int) const { return Face_handle(); }
    //      int mirror_index(int) const { return 0; }
  };

  typedef Vertex_handle Vertex;
  typedef Face_handle   Face;

  template<class V, class F>
  struct Data_structure_t {
    typedef F Face_handle;

    int mirror_index(Face_handle, int) const { return 0; }
  };

  typedef Data_structure_t<Vertex_handle,Face_handle> Data_structure;
  typedef Data_structure        Triangulation_data_structure;

  //  typedef typename Data_structure::Face_handle     Face_handle;
  //  typedef typename Data_structure::Vertex_handle   Vertex_handle;

  static const Vertex_handle infinite_vertex() {
    static Vertex_handle inf_v;
    return inf_v;
  }

  static const Data_structure& data_structure() {
    static Data_structure ds;
    return ds;
  }

  static const Data_structure& tds() {
    return data_structure();
  }

  typedef std::pair<Face_handle,int>  Edge;

 private:
  class Dummy_edge_iterator : public Emptyset_iterator
  {
  public:
    typedef unsigned int size_type;

    typedef Edge   value_type;
    typedef const Edge*  pointer;
    typedef Edge&  reference;

    Dummy_edge_iterator() {}
    Dummy_edge_iterator(const Dummy_iterator&) {}

    template< class T >
    Dummy_edge_iterator& operator=(const T&) { return *this; }

    Dummy_edge_iterator& operator++()        { return *this; }
    Dummy_edge_iterator& operator++(int)     { return *this; }

    Edge& operator*()         { return e; }
    Edge* operator->()        { return &e; }

    const Edge& operator*() const   { return e; }
    const Edge* operator->() const  { return (&e); }

    bool operator==(const Dummy_edge_iterator&) const {
      return true;
    }

    bool operator!=(const Dummy_edge_iterator&) const {
      return false;
    }

    bool operator<(const Dummy_edge_iterator& other) const {
      return this < &other;
    }

  protected:
    Edge e;
  };

  static const Dummy_edge_iterator& dummy_edge_it() {
    static Dummy_edge_iterator dummy_edge_it_static;
    return dummy_edge_it_static;
  }

 public:
  typedef Dummy_edge_iterator  All_edges_iterator;
  struct Finite_edges_iterator : public Dummy_edge_iterator {
    Finite_edges_iterator() {}
    Finite_edges_iterator(const Dummy_edge_iterator&) {}
  };

  typedef Dummy_iterator  All_faces_iterator;
  typedef Dummy_iterator  Finite_faces_iterator;

  typedef Dummy_iterator  All_vertices_iterator;
  typedef Dummy_iterator  Finite_vertices_iterator;

  typedef Dummy_iterator  Vertex_circulator;
  typedef Dummy_iterator  Face_circulator;

  struct Edge_circulator : public Dummy_edge_iterator
  {
    Edge_circulator() {}

    template< class T >
    Edge_circulator& operator=(const T&) { return *this; }

    Edge_circulator& operator++()        { return *this; }
    Edge_circulator& operator++(int)     { return *this; }

    Edge& operator*()         { return this->e; }
    Edge* operator->()        { return &this->e; }

    const Edge& operator*() const   { return this->e; }
    const Edge* operator->() const  { return &this->e; }

    bool operator==(const Edge_circulator&) const {
      return true;
    }

    bool operator!=(const Edge_circulator&) const {
      return false;
    }
  };

  

 private:
  static const Edge_circulator& dummy_edge_circulator() {
    static Edge_circulator dummy_edge_circulator_static;
    return dummy_edge_circulator_static;
  }

 public:
  void insert(const Site_2&) {}

  size_type number_of_vertices() const { return 0; }
  size_type number_of_faces() const { return 0; }
  size_type number_of_edges() const { return 0; }

  size_type degree(Vertex_handle v) const { return 0; }

  int dimension() const { return -1; }

  bool is_infinite(const Dummy_iterator&) const { return false; }
  bool is_infinite(const Dummy_edge_iterator&) const { return false; }
  bool is_infinite(const Vertex_handle&) const { return false; }
  bool is_infinite(const Face_handle&) const { return false; }
  bool is_infinite(const Edge&) const { return false; }
  bool is_infinite(const Face_handle&,int) const { return false; }

  All_edges_iterator all_edges_begin() const {
    return dummy_edge_it();
  }

  All_edges_iterator all_edges_end() const {
    return dummy_edge_it();
  }

  Finite_edges_iterator finite_edges_begin() const {
    return dummy_edge_it();
  }

  Finite_edges_iterator finite_edges_end() const {
    return dummy_edge_it();
  }

  All_vertices_iterator all_vertices_begin() const {
    return dummy_it();
  }

  All_vertices_iterator all_vertices_end() const {
    return dummy_it();
  }

  Finite_vertices_iterator finite_vertices_begin() const {
    return dummy_it();
  }

  Finite_vertices_iterator finite_vertices_end() const {
    return dummy_it();
  }

  All_faces_iterator all_faces_begin() const {
    return dummy_it();
  }

  All_faces_iterator all_faces_end() const {
    return dummy_it();
  }

  Finite_faces_iterator finite_faces_begin() const {
    return dummy_it();
  }

  Finite_faces_iterator finite_faces_end() const {
    return dummy_it();
  }

  Edge_circulator incident_edges(const Vertex_handle&) const {
    return dummy_edge_circulator();
  }

  Edge_circulator incident_edges(const Vertex_handle&,
				 const Face_handle&) const {
    return dummy_edge_circulator();
  }

  Vertex_circulator incident_vertices(const Vertex_handle&) const {
    return dummy_it();
  }

  Face_circulator incident_faces(const Vertex_handle&) const {
    return dummy_it();
  }


  Point_2 dual(const Face_handle&) const { return Point_2(); }

  bool is_valid(bool = true,int = 0) const { return true; }
};


template<class K>
std::istream&
operator>>(std::istream& is,
	   const typename Delaunay_graph_concept<K>::Site_2&)
{
  double d;
  is >> d;
  return is;
}

template<class K>
std::ostream&
operator<<(std::ostream& os,
	   const typename Delaunay_graph_concept<K>::Site_2&)
{
  return os;
}

CGAL_END_NAMESPACE


#endif // CGAL_DELAUNAY_GRAPH_CONCEPT_H
