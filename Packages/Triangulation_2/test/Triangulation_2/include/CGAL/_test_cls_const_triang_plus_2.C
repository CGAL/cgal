#include <CGAL/_test_cls_const_Del_triangulation_2.C>

template <class TriangPlus>
void
_test_cls_const_triang_plus_2( const TriangPlus & )
{
  typedef TriangPlus                           TrP;
  typedef typename TrP::Geom_traits            Gt;
  typedef typename Gt::Point_2                 Point;

  typedef typename TrP::Vertex_handle          Vertex_handle;
  typedef typename TrP::Constraint             Constraint;
  typedef typename TrP::Constraint_hierarchy   Hierarchy;
  typedef typename TrP::Context                Context;
  typedef typename TrP::Context_iterator       Context_iterator;
  typedef typename TrP::Vertices_in_constraint Vertices_in_constraint;

  // _test_cls_const_Del_triangulation( TrP() );

  Point pt[12] = {
     Point(0,0), Point(0,4), 
     Point(1,0), Point(1,4), 
     Point(2,0), Point(2,4),
     Point(-1,1), Point(3,1),
     Point(-1,2), Point(3,2),
     Point(0.5,1), Point(2.5,1)
   };

  TrP trp;
  Vertex_handle vh[12];
  for(int i=0; i<12; i++){
    vh[i] = trp.insert(pt[i]);
  }
  for(int j=0; j<11; j+=2){
    trp.insert(vh[j],vh[j+1]);
  }

  trp.insert(Point(4,4), Point(4,5));
  trp.push_back(Point(4,6));
  trp.push_back(Constraint(Point(4,3), Point(3,4)));


  Vertices_in_constraint vit = trp.vertices_in_constraint_begin(vh[10],vh[11]);
  assert (*vit == vh[10]);
  Vertex_handle va = *++vit;
  Vertex_handle vb = *++vit;
  assert (*++vit == vh[11]);
  assert (++vit == trp.vertices_in_constraint_end(vh[10],vh[11]));
  assert(trp.number_of_enclosing_constraints(va,vb) == 2);
  Context_iterator cit1 = trp.contexts_begin(va,vb);
  Context_iterator cit2 = cit1++;
  assert( cit1->number_of_vertices() == 4  || cit1->number_of_vertices() == 7);
  Vertices_in_constraint vit1 = cit1->first(); 
  Vertices_in_constraint vit2 = cit2->first();
  if ( cit1->number_of_vertices() == 4 ) {
    assert(*vit1 == vh[10]);
    assert(*vit2 == vh[6] );
    assert(*--(cit1->past()) == vh[11]);
    assert( *--(cit2->past()) == vh[7]);
  }
  else {
    assert(*vit1 == vh[6]);
    assert(*vit2 == vh[10]);
    assert(*--(cit1->past()) == vh[7]);
    assert(*--(cit2->past()) == vh[11]);
  }
  assert(*(cit1->current()) == va);
  assert( *(cit2->current()) == va);
  return;
  
  
}


