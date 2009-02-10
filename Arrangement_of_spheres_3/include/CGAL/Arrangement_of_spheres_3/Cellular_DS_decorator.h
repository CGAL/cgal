#ifndef CGAL_AOS3_CELLULAR_DS_DECORATOR_H
#define CGAL_AOS3_CELLULAR_DS_DECORATOR_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

template <class CDS>
class Cellular_DS_decorator {
  CDS cds_;
public:
  
  typedef typename CDS::Vertex_handle Vertex_handle;
  typedef typename CDS::Dart_handle Dart_handle;
  typedef typename CDS::Halfface_handle Halffacet_handle;
  typedef typename CDS::Cell_handle Cell_handle;

  typedef typename CDS::Vertex Vertex;
  typedef typename CDS::Dart Dart;
  typedef typename CDS::Halfface Halffacet;
  typedef typename CDS::Cell Cell;
  typedef typename CDS::Vertex::Description Point;
  typedef typename CDS::Dart::Description Curve;
  typedef typename CDS::Halfface::Description Surface;


  typedef Sphere_key Key;

  Cellular_DS_decorator(CDS &cds): cds_(cds){}

  // create cube
  void create_cube() {
    Cell_handle outside= cds_.cells.push_back(Cell());
    Cell_handle inside= cds_.cells.push_back(Cell());

    Surface sides[2][3];
    sides[0][0]=Surface(Key::lb_key(), Surface::FRONT_FACE);
    sides[0][1]=Surface(Key::lb_key(), Surface::R_RULE);
    sides[0][2]=Surface(Key::lb_key(), Surface::T_RULE);
    sides[1][0]=Surface(Key::ub_key(), Surface::BACK_FACE);
    sides[1][1]=Surface(Key::ub_key(), Surface::L_RULE);
    sides[1][2]=Surface(Key::ub_key(), Surface::B_RULE);

    Halfface_handle facets[2][3];
    //for (int i=0; i< 2; ++i) {
    for (int j=0; j<3; ++j) {
      facets[0][j]= cds_.halffacets.push_back(Halffacet(sides[0][j]),
					      Halffacet(sides[0][j].opposite()));
      facets[1][j]= cds_.halffacets.push_back(Halffacet(sides[1][j].opposite()),
					      Halffacet(sides[1][j]));
      facets[0][j]->set_cell(inside);
      facets[1][j]->set_cell(outside);
    }

    inside->set_facet(facets[0][0]);
    outside->set_facet(facets[1][0]);
      //}
    
    Vertex_handle vertices[2][2][2];
    Point points[2][2][2];
    for (int i=0; i<2; ++i) {
      for (int j=0; j< 2; ++j) {
	for (int k=0; k<2; ++k) {
	  points[i][j][k]= Point(sides[i][0],
				 sides[j][1],
				 sides[k][2]);
	  vertices[i][j][k]=cds_.vertices_push_back(Vertex(points[i][j][k]));
	}
      }
    }
    
    new_edge(facets[0][0], facets[0][1], 
	     vertices[0][0][0], vertices[0][1][0]);
    new_edge(facets[0][0], facets[1][2], 
	     vertices[0][1][0], vertices[1][1][0]);
    new_edge(facets[0][0], facets[1][1], 
	     vertices[1][1][0], vertices[1][0][0]);
    new_edge(facets[0][0], facets[0][2], 
	     vertices[1][0][0], vertices[0][0][0]);

    new_edge(facets[1][1], facets[0][2], 
	     vertices[1][0][1], vertices[1][0][0]);
    new_edge(facets[0][2], facets[0][1], 
	     vertices[0][0][1], vertices[0][0][0]);
    new_edge(facets[0][1], facets[1][2], 
	     vertices[1][1][0], vertices[0][1][0]);
    new_edge(facets[1][2], facets[1][1], 
	     vertices[1][1][1], vertices[1][1][0]);

    new_edge(facets[1][0], facets[1][1], 
	     vertices[1][0][1], vertices[1][1][1]);
    new_edge(facets[1][0], facets[0][2], 
	     vertices[0][0][1], vertices[1][0][1]);
    new_edge(facets[1][0], facets[0][1], 
	     vertices[1][1][0], vertices[0][0][1]);
    new_edge(facets[1][2], facets[1][0], 
	     vertices[1][1][0], vertices[1][1][1]);


    
  }

  void audit() {
    for (CDS::Cell_iterator cit= cds_.cells_begin(); cit != cds_.cells_end();
	 ++cit) {
      // use 2D audit
      CGAL::HalfedgeDS_const_decorator<CDS::Cell> chds(*cit);
      if (!chds.is_valid(true, 1)) {
	CGAL_assertion(0);
      }
      
      // all cell pointers point to right thing
      for (CDS::Cell::Face_iterator it= cit->faces_begin();
	   it != cit->faces_end(); ++it) {
	CGAL_assertion(it->cell() == cit);
      }
      
    }

    // both sides of a halffacet have same edges
    for (CDS::Halffacet_iterator cit= cds_.halffacets_begin();
	 cit != cds_.halffacets_end();
	 ++cit) {
      Halffacet_handle fh=cit;
      Halffacet_handle ofh=fh->opposite();
      Dart_handle dh= fh->dart();
      Dart_handle dhm= dh->mirror();
      do {
	CGAL_assertion(dh->mirror() == dhm);
	dh=dh->next();
	dhm=dhm->prev();
      } while (dh != fh->dart());
    }

    // all pointers non-null
      
    // labels correct!


  }

  // close cell to a point
  // end cells with a new cell
  // insert HDS in cell and extrude
  // cross edges
  // edge across vertex
  // open face in surface
  // verify
};

CGAL_AOS3_END_INTERNAL_NAMESPACE
#endif
