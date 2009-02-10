#ifndef CGAL_AOS3_COMBINATORIAL_ARRANGEMENT_H
#define CGAL_AOS3_COMBINAOTIRAL_ARRANGEMENT_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

class Combinatorial_arrangement {
public:
  typedef Sphere_key Sphere_3_key;
  /*struct Surface {
    enum Part {X,Y,S};
    Surface(){}
    static Surface sphere(Sphere_3_key k) {
      return Surface(k, S);
    }
    static Surface rule(Sphere_3_key k, Coordinate_index r){
      if (r== plane_coordiante(0)) return Surface(k, X);
      else return Surface(k,Y);
    }
  private:
    Surface(Sphere_3_key k, int p): k_(k), p_(p){}
    Sphere_3_key k_;
    int p_;
    };*/
  typedef Combinatorial_curve Surface;
  struct Point {
    /*
      Point types:
      S1, S2,S3 2 of them per tuple
      E[SR][SR] 2 of them: front/back of sphere hitting a vertex
      S[SR][SR]
      hmmm, just maintain an unordered tuple

      ((K,[SR])^3, [01])-- three curves and whether it is the first or second
     */
  private:
    Surface surfaces_[3];
    int index_;
  };
  typedef Combinatorial_vertex Curve;
  
  struct Cellular_DS_traits {
    struct Point_2{};// do I really need this?
  };

  struct Cellular_DS_items {
    typedef CellularDS_default_items P;
    template <class Refs, class Traits>
    struct Vertex_wrapper {
      struct Vertex: public CellularDS_vertex_base<Refs> {
	Vertex(){}
	Vertex(const Point &pt): pt_(pt){}
	const Point &point() const {
	  return pt_;
	}
	typedef Point Description;
      private:
	Point pt_;
      };
    };


    template <class Refs, class CDS>
    struct Dart_wrapper {
      struct Dart: public CellularDS_dart_base<Refs, CDS> {
	Dart(){}
	Dart(const Curve &c): c_(c){}
	const Curve &curve() const {
	  return c_;
	}
	typedef Curve Description;
      private:
	Curve c_;
      };
    };

    template < class Refs, class Traits, class CDS> 
    struct Halffacet_wrapper { 
      struct Halffacet: public CellularDS_halffacet_base<Refs, CDS> {
	Halffacet(){}
	Halffacet(const Surface &c): c_(c){}
	const Slice_curve surface() const {
	  return c_;
	}
	typedef Surface Description;
      private:
	Surface c_;
      };
    }; 
    template < class Refs, class Traits, class HDS> 
    struct Cell_wrapper { 
      struct Cell: public CellularDS_cell_base<Refs, HDS> {
	
      private:
	std::vector<Sphere_3_key> containing_;
      }; 
    }; 
    
  };

  typedef CellularDS_list<, Items> CDS;

  

  /*
    basic operations:
    -- insert vertex in edge
        - add a vertex to each cell, 
	- insert in each edges
	- connect new edge across cells
    -- connect two vertices with edge
        - split face in each cell
	- connect the new edges
	- connect the new facets
    -- merge two cells through a face
        - copy faces from one cell into the other, build map
	- connect opposites
    -- create a new cell which cuts then end of a given cell along some edges
        - move end facet forward
	- create new cell with two copies of end facet
	- create side facets
	- connect up
    -- end a facet
        - add new vertex, connect two edges to vertex
	- 

    Key operations:
    -- cut off end of cell parallel to sweep plane
          -- create new cell with two copies of end face 
    -- star a face of a cell
    -- end a facet in a vertex
    -- replace a vertex in the sweep plane with a hds and connect back to the vertex
    -- replace a hds in the sweep plane with a vertex
    -- unstar some faces

    
   */

  
private:
  CDS cds_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#endif
