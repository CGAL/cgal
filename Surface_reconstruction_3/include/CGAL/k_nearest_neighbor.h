#ifndef CGAL_K_NEIGHBOR_NEIGHBOR_H
#define CGAL_K_NEIGHBOR_NEIGHBOR_H

#include <CGAL/basic.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <list>

CGAL_BEGIN_NAMESPACE


template <class Vertex_handle>
struct KVertex
{
  double m_coord[3];
	Vertex_handle m_vertex_handle;

  KVertex()
  {
    m_coord[0] =
    m_coord[1] =
    m_coord[2] = 0;
    m_vertex_handle = NULL;
  }

  KVertex(const KVertex& v)
  {
    m_coord[0] = v.x();
    m_coord[1] = v.y();
    m_coord[2] = v.z();
    m_vertex_handle = v.vertex_handle();
  }

  KVertex(double x,
          double y,
          double z,
          Vertex_handle vertex_handle)
  {
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
    m_vertex_handle = vertex_handle;
  }

	Vertex_handle& vertex_handle() { return m_vertex_handle; }
	const Vertex_handle& vertex_handle() const { return m_vertex_handle; }

  const double x() const { return m_coord[0]; }
  const double y() const { return m_coord[1]; }
  const double z() const { return m_coord[2]; }

  double& x() { return m_coord[0]; }
  double& y() { return m_coord[1]; }
  double& z() { return m_coord[2]; }

  bool operator==(const KVertex& p) const
  {
    return (x() == p.x()) &&
           (y() == p.y()) &&
           (z() == p.z());
  }

  bool  operator!=(const KVertex& p) const { return ! (*this == p); }
}; // end of class KVertex



template <class Vertex_handle>
struct Kernel_traits< KVertex<Vertex_handle> > {
  struct Kernel {
    typedef double FT;
    typedef double RT;
  };
};


template <class Vertex_handle>
struct Construct_coord_iterator {

	typedef typename CGAL::KVertex<Vertex_handle> KVertex;

  const double* operator()(const KVertex& p) const
  { return static_cast<const double*>(p.m_coord); }

  const double* operator()(const KVertex& p, int)  const
  { return static_cast<const double*>(p.m_coord+3); }
};

// We have put the glue layer in this file as well, that is a class that
// allows to iterate over the Cartesian coordinates of the KVertex, and a
// class to construct such an iterator for a KVertex.
// We next need a distance class

template <class Vertex_handle>
struct Distance {
	typedef typename CGAL::KVertex<Vertex_handle> KVertex;
  typedef KVertex Query_item;

  double transformed_distance(const KVertex& p1, const KVertex& p2) const {
    double distx = p1.x()-p2.x();
    double disty = p1.y()-p2.y();
    double distz = p1.z()-p2.z();
    return distx * distx +
			     disty * disty +
					 distz * distz;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const KVertex& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    double distance(0.0), h = p.x();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.y();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.z();
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }

  template <class TreeTraits>
  double max_distance_to_rectangle(const KVertex& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    double h = p.x();

    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);

    h=p.y();
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.z();
    double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }

  double new_distance(double& dist, double old_off, double new_off,
		      int cutting_dimension)  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return sqrt(d); }

}; // end of struct Distance


template <class Gt, class Vertex_handle>
class K_nearest_neighbor
{
public:
  typedef Gt  Geom_traits;

	typedef typename Geom_traits::FT FT;
	typedef typename Geom_traits::Point_3 Point;
	typedef typename CGAL::KVertex<Vertex_handle> KVertex;
	typedef typename CGAL::Construct_coord_iterator<Vertex_handle> KConstruct_coord_iterator;
  typedef typename CGAL::Search_traits<double,
                                       KVertex,
                                       const double*,
                                       KConstruct_coord_iterator> Traits;
	typedef typename CGAL::Distance<Vertex_handle> KDistance;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Traits,KDistance> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Iterator;

private:
  Tree m_tree;

public:
  K_nearest_neighbor() {}
  ~K_nearest_neighbor() {}

  void init(const std::list<KVertex>& vertices)
  {
    m_tree = Tree(vertices.begin(),vertices.end());
  }

	void clear()
	{
		m_tree = Tree();
	}

  bool k_nearest_neighbors(const KVertex& query,
		                       const unsigned int nb,
													 std::list<KVertex>& kvertices)
  {
    Neighbor_search search(m_tree,query,nb); // only n nearest neighbors
		Iterator it = search.begin();
		for(unsigned int i=0;i<nb;i++,it++)
		{
			if(it == search.end())
				return false;
			kvertices.push_back(it->first);
		}
		return true;
  }

};


CGAL_END_NAMESPACE

#endif // CGAL_K_NEIGHBOR_NEIGHBOR_H

