#include <CGAL/_test_cls_const_Del_triangulation_2.h>
#include <CGAL/use.h>
template <class TrP>
void
_test_cls_const_triang_plus_2( const TrP & )
{
  //typedef TriangPlus                           TrP;
  typedef typename TrP::Geom_traits            Gt;
  typedef typename Gt::Point_2                 Point;

  typedef typename TrP::Vertex_handle          Vertex_handle;
  typedef typename TrP::Constraint             Constraint;
  typedef typename TrP::Constraint_hierarchy   Hierarchy;
  typedef typename TrP::Context                Context;
  typedef typename TrP::Context_iterator       Context_iterator;
  typedef typename TrP::Vertices_in_constraint Vertices_in_constraint;

  CGAL_USE_TYPE(Hierarchy);
  CGAL_USE_TYPE(Context);
  std::cout << " call test of constrained triangulations" <<std::endl;
  _test_cls_const_Del_triangulation( TrP() );

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
     trp.insert_constraint(vh[j],vh[j+1]);
  }

  trp.insert(Point(4,4), Point(4,5));
  trp.push_back(Point(4,6));
  trp.push_back(Constraint(Point(4,3), Point(3,4)));

  // test access to the hierarchy
  std::cout << " test acces to the constraint hierarchy" << std::endl;
  Vertices_in_constraint vit = trp.vertices_in_constraint_begin(vh[10],vh[11]);
  assert (*vit == vh[10] || *vit == vh[11] );
  Vertex_handle va = *++vit;
  Vertex_handle vb = *++vit;
  assert (*++vit == vh[11] || *vit == vh[10]);
  assert (++vit == trp.vertices_in_constraint_end(vh[10],vh[11]));
  assert(trp.number_of_enclosing_constraints(va,vb) == 2);
  Context_iterator cit1 = trp.contexts_begin(va,vb);
  Context_iterator cit2 = cit1++;
  //trp.print_hierarchy();
  assert( cit1->number_of_vertices() == 4  || cit1->number_of_vertices() == 7);
  Vertices_in_constraint firstin1 = cit1->vertices_begin();
  Vertices_in_constraint lastin1 = --(cit1->vertices_end());
  Vertices_in_constraint currentin1 = cit1->current();
  Vertices_in_constraint firstin2 = cit2->vertices_begin();
  Vertices_in_constraint lastin2 = --(cit2->vertices_end());
  Vertices_in_constraint currentin2 = cit2->current();
  if ( cit1->number_of_vertices() == 4) {
    assert( (*firstin1 == vh[10] &&  *lastin1 == vh[11]) ||
	    (*firstin1 == vh[11] &&  *lastin1 == vh[10]));
    assert( (*firstin2 == vh[6] &&  *lastin2 == vh[7]) ||
	    (*firstin2 == vh[7] &&  *lastin2 == vh[6]));
  }
  else {
    assert( (*firstin1 == vh[6] &&  *lastin1 == vh[7]) ||
	    (*firstin1 == vh[7] &&  *lastin1 == vh[6]));
    assert( (*firstin2 == vh[10] &&  *lastin2 == vh[11]) ||
	    (*firstin2 == vh[11] &&  *lastin2 == vh[10]));
  }
  assert( (*currentin1 == va &&  *++currentin1 == vb) ||
	  (*currentin1 == vb &&  *++currentin1 == va));
  assert( (*currentin2 == va &&  *++currentin2 == vb) ||
	  (*currentin2 == vb &&  *++currentin2 == va));

  //test copy and swap
  std::cout << "test copy and swap" << std::endl;
  TrP  trp2(trp);
  TrP  trp3 = trp2;
  // trp.print_hierarchy();
  // trp2.print_hierarchy();
  // trp3.print_hierarchy();
  assert(trp3.number_of_constraints() == trp.number_of_constraints());
  assert(trp3.number_of_subconstraints() == trp.number_of_subconstraints());
  trp2.clear();
  //trp2.print_hierarchy();
  trp2.swap(trp3);
  assert(trp2.number_of_constraints() == trp.number_of_constraints());
  assert(trp2.number_of_subconstraints() == trp.number_of_subconstraints());
  assert(trp3.number_of_constraints() == 0);
  // trp2.print_hierarchy();
  // trp3.print_hierarchy();

  //test remove_constraint
  std::cout << " test removal of constraint" << std::endl;
  trp.remove_constraint(vh[10],vh[11]);
  trp.remove_constraint(vh[6],vh[7]);

  std::cerr << " test a special configuration" << std::endl;
  trp.clear();
  Vertex_handle v1 = trp.insert(Point(0, 0));
  Vertex_handle v2 = trp.insert(Point(-1, 0));
  Vertex_handle v3 = trp.insert(Point(1, 0));
  trp.insert_constraint(v1, v2);
  trp.insert_constraint(v2, v3);
  trp.insert_constraint(v3, v1);
  trp.remove_constraint(v1, v2);
  trp.remove_constraint(v2, v3);
  trp.remove_constraint(v3, v1);

  std::cout << std::endl;
  return;
}


