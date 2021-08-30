#include <CGAL/_test_cls_const_Del_triangulation_2.h>
#include <CGAL/use.h>
template <class TrP>
void
_test_cls_const_triang_plus_2( const TrP & )
{
  //typedef TriangPlus                           TrP;
  typedef typename TrP::Geom_traits                     Gt;
  typedef typename Gt::Point_2                          Point;

  typedef typename TrP::Vertex_handle                   Vertex_handle;
  typedef typename TrP::Constraint                      Constraint;
  typedef typename TrP::Constraint_iterator             Constraint_iterator;
  typedef typename TrP::Constraint_hierarchy            Hierarchy;
  typedef typename TrP::Context                         Context;
  typedef typename TrP::Context_iterator                Context_iterator;
  typedef typename TrP::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
  typedef typename TrP::Constraint_id                   Constraint_id;

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
  Constraint_id cid;
  for(int j=0; j<11; j+=2){
     cid = trp.insert_constraint(vh[j],vh[j+1]);
  }

  trp.insert(Point(4,4), Point(4,5));
  trp.push_back(Point(4,6));
  trp.push_back(Constraint(Point(4,3), Point(3,4)));

  // test access to the hierarchy
  std::cout << " test acces to the constraint hierarchy" << std::endl;
  Vertices_in_constraint_iterator vit = trp.vertices_in_constraint_begin(cid);
  assert (*vit == vh[10] || *vit == vh[11] );
  Vertex_handle va = *++vit;
  Vertex_handle vb = *++vit;
  assert (*++vit == vh[11] || *vit == vh[10]);
  assert (++vit == trp.vertices_in_constraint_end(cid));
  assert(trp.number_of_enclosing_constraints(va,vb) == 2);
  Context_iterator cit1 = trp.contexts_begin(va,vb);
  Context_iterator cit2 = cit1++;
  //trp.print_hierarchy();
  assert( cit1->number_of_vertices() == 4  || cit1->number_of_vertices() == 7);
  Vertices_in_constraint_iterator firstin1 = cit1->vertices_begin();
  Vertices_in_constraint_iterator lastin1 = --(cit1->vertices_end());
  Vertices_in_constraint_iterator currentin1 = cit1->current();
  Vertices_in_constraint_iterator firstin2 = cit2->vertices_begin();
  Vertices_in_constraint_iterator lastin2 = --(cit2->vertices_end());
  Vertices_in_constraint_iterator currentin2 = cit2->current();
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
  trp.remove_constraint(cid);

  std::cerr << " test a special configuration" << std::endl;
  trp.clear();
  Vertex_handle v1 = trp.insert(Point(0, 0));
  Vertex_handle v2 = trp.insert(Point(-1, 0));
  Vertex_handle v3 = trp.insert(Point(1, 0));
  Constraint_id cid1 = trp.insert_constraint(v1, v2);
  Constraint_id cid2 = trp.insert_constraint(v2, v3);
  Constraint_id cid3 = trp.insert_constraint(v3, v1);
  trp.remove_constraint(cid1);
  trp.remove_constraint(cid2);
  trp.remove_constraint(cid3);

  std::cerr << " test the configuration of bug #2999" << std::endl;
  trp.clear();
  {
    Point const p1(942455,   2306674.6);
    Point const p2(942452.3, 2306671.9);
    Point const p3(942441.1, 2306667.8);
    Point const p4(942453.8, 2306673.4);
    Point const p6(942447.7, 2306686  );
    trp.insert_constraint(p1,p2);
    trp.insert_constraint(p3,p4); // that insert an intersection point p5,
                                  // that should be snapped to p4
    trp.insert_constraint(p4,p6);
  }

  std::cout << std::endl;
  std::cout << "test IO" << std::endl;

  trp.clear();
  {
    Point c1[2] = { Point(0,0), Point(1,0) };
    Point c2[3] = { Point(0,1), Point(1,1), Point(2,1) };
    Point c3[4] = { Point(0,0), Point(1,0), Point(1,1), Point(0,1) };

    trp.insert_constraint(c1, c1 + 2);
    trp.insert_constraint(c2, c2 + 3);
    trp.insert_constraint(c3, c3 + 4);
    std::ofstream out("cdtplus.txt");
    out << trp;
    out.close();
    trp.clear();
    std::ifstream in("cdtplus.txt");
    in >> trp;
    assert(trp.number_of_constraints() == 3);
    std::size_t n = 0;
    for(Constraint_iterator cit = trp.constraints_begin(); cit != trp.constraints_end(); ++cit){
      Constraint_id  cid = *cit;
      n += cid.second->all_size();
    }
    assert( n == 9);
  }
  return;
}


