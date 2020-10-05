#ifndef _REFINER_H_
#define _REFINER_H_

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <queue>

template <class  Kernel,  class  Polyhedron>
class CEdge
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

private:
  FT m_sqlen;
  Halfedge_handle m_he;

public:
  // life cycle
  CEdge(const Halfedge_handle& he)
  {
    m_sqlen = squared_len(he);
    m_he = he;
  }
  CEdge(const CEdge& edge)
    : m_sqlen(edge.sqlen()),
    m_he(edge.he())
  {
  }
  ~CEdge() {}

public:
  // squared length of an edge
  static FT squared_len(Halfedge_handle he)
  {
    return CGAL::squared_distance(he->vertex()->point(),
      he->opposite()->vertex()->point());
  }

public:
  // accessors
  FT& sqlen() { return m_sqlen; }
  const FT sqlen() const { return m_sqlen; }

  Halfedge_handle he() { return m_he; }
  const Halfedge_handle he() const { return m_he; }
};

// functor for priority queue
template<class Edge>
struct less // read more priority
{
  bool operator()(const Edge& e1,
    const Edge& e2) const
  {
    return e1.sqlen() < e2.sqlen();
  }
};

template <class  Kernel,  class  Polyhedron>
class Refiner
{
  // types
  typedef typename Kernel::FT FT;
  typedef CEdge<Kernel, Polyhedron> Edge;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Edge_iterator   Edge_iterator;
  typedef std::priority_queue<Edge,
    std::vector<Edge>,
    ::less<Edge> > PQueue;
  // data
  PQueue m_queue;
  Polyhedron* m_pMesh;
public :
  // life cycle
  Refiner(Polyhedron* pMesh)
  {
    m_pMesh = pMesh;
  }
  ~Refiner() {}

public :

  void fill_queue(const FT& max_sqlen)
  {
    for(Edge_iterator he = m_pMesh->edges_begin();
      he != m_pMesh->edges_end();
      he++)
      if(Edge::squared_len(he) > max_sqlen)
        m_queue.push(Edge(he));
  }

  void fill_queue()
  {
    for(Edge_iterator he = m_pMesh->edges_begin();
      he != m_pMesh->edges_end();
      he++)
      m_queue.push(Edge(he));
  }

  // run
  void run_nb_splits(const unsigned int nb_splits)
  {
    // fill queue
    fill_queue();

    unsigned int nb = 0;
    while(nb < nb_splits)
    {
      // get a copy of the candidate edge
      Edge edge = m_queue.top();
      m_queue.pop();

      Halfedge_handle he = edge.he();
      // update point
      Halfedge_handle hnew = m_pMesh->split_edge(he);
      hnew->vertex()->point() = CGAL::midpoint(he->vertex()->point(),
        he->opposite()->vertex()->point());

      // hit has been split into two edges
      m_queue.push(Edge(hnew));
      m_queue.push(Edge(he));

      // split facet if possible
      if(!hnew->is_border())
      {
        m_pMesh->split_facet(hnew,hnew->next()->next());
        m_queue.push(Edge(hnew->next()));
      }

      // split facet if possible
      if(!hnew->opposite()->is_border())
      {
        m_pMesh->split_facet(hnew->opposite()->next(),
          hnew->opposite()->next()->next()->next());
        m_queue.push(Edge(hnew->opposite()->prev()));
      }

      nb++;
    } // end while
  } // end run



  // run
  unsigned int operator()(const FT& max_sqlen)
  {
    // fill queue
    fill_queue(max_sqlen);

    unsigned int nb_split = 0;
    while(!m_queue.empty())
    {
      // get a copy of the candidate edge
      Edge edge = m_queue.top();
      m_queue.pop();

      Halfedge_handle he = edge.he();
      FT sqlen = Edge::squared_len(he);
      if(sqlen > max_sqlen)
      {
        // update point
        Halfedge_handle hnew = m_pMesh->split_edge(he);
        hnew->vertex()->point() = CGAL::midpoint(he->vertex()->point(),
          he->opposite()->vertex()->point());

        // hit has been split into two edges
        m_queue.push(Edge(hnew));
        m_queue.push(Edge(he));

        // split facet if possible
        if(!hnew->is_border())
        {
          m_pMesh->split_facet(hnew,hnew->next()->next());
          m_queue.push(Edge(hnew->next()));
        }

        // split facet if possible
        if(!hnew->opposite()->is_border())
        {
          m_pMesh->split_facet(hnew->opposite()->next(),
            hnew->opposite()->next()->next()->next());
          m_queue.push(Edge(hnew->opposite()->prev()));
        }

        nb_split++;
      } // end if(sqlen > max_sqlen)
    } // end while(!m_queue.empty())
    return nb_split;
  } // end run
};

#endif // _REFINER_H_


