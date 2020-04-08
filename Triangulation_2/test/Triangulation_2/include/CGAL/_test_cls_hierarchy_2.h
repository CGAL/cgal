#include <CGAL/Triangulation_2/internal/Constraint_hierarchy_2.h>

struct Less
{
  bool operator()(int i, int j) const
  {
    return i < j;
  }
};

typedef int Vh;
typedef bool Data;
typedef CGAL::Constraint_hierarchy_2<Vh,
                                     Less,
                                     Data>      Hierarchy;
typedef Hierarchy::H_constraint                 H_constraint;
typedef Hierarchy::H_vertex                     H_vertex;
typedef Hierarchy::H_vertex_it                  H_vertex_it;
typedef Hierarchy::H_constraint_list            H_constraint_list;
typedef Hierarchy::H_constraint_list::iterator  H_constraint_it;
typedef Hierarchy::H_c_iterator                 H_c_iterator;
typedef Hierarchy::H_sc_iterator                H_sc_iterator;
typedef Hierarchy::H_vertex_map                 H_vertex_map;
typedef Hierarchy::H_context                    H_context;
typedef Hierarchy::H_context_iterator           H_context_iterator;

void
_test_cls_hierarchy_2()
{
  Vh v[10];
  for(int i=0; i <10; i++) { v[i] = i;}

  Hierarchy h;
  h.insert_constraint(v[1],v[2]);
  h.insert_constraint(v[4],v[3]);
  h.split_constraint(v[1],v[2],v[5]);
  h.split_constraint(v[4],v[3],v[6]);
  h.split_constraint(v[3],v[6],v[7]);
  H_constraint c;
  h.enclosing_constraint(std::make_pair(v[1], v[5]), c);
  assert(c ==  std::make_pair(v[1],v[2]));
  h.enclosing_constraint(std::make_pair(v[3], v[7]), c);
  assert(c ==  std::make_pair(v[3],v[4]));

  // result should be 152 3764
  // h.print();
  H_context co = h.context(v[1],v[5]);
    H_vertex_it vit = co.vertices_begin();
  assert( *vit++ == v[1] &&
          *vit++ == v[5] &&
          *vit++ == v[2] &&
          vit == co.vertices_end() &&
          *(co.current()) == v[1]) ;
  co = h.context(v[7],v[6]);
  vit = co.vertices_begin();
  assert( *vit++ == v[3]);
  assert( *vit++ == v[7]);
  assert( *vit++ == v[6]);
  assert( *vit++ == v[4]);
  assert( vit == co.vertices_end());
  assert(*(co.current()) == v[7]) ;
  co = h.context(v[6],v[7]);
  assert(*(co.current()) == v[7]) ;
  H_vertex_it v_in_c = h.vertices_in_constraint_begin(v[4],v[3]);
  assert(*v_in_c == v[3]);
  assert(*++v_in_c == v[7]);
  assert(*++v_in_c == v[6]);
  assert(*++v_in_c == v[4]);
  assert(++v_in_c == h.vertices_in_constraint_end(v[4],v[3]) );

  h.constrain_vertex(v[6]);
  h.set_data(v[6], true);
  assert(h.get_data(v[6]));
  assert(h.is_constrained_vertex(v[6]));
  h.unconstrain_vertex(v[6]);
  h.remove_Steiner(v[6],v[4],v[7]);

  //h. print();
  // result should be 152  374

  Vh w;
  h.oriented_end(v[1],v[5],w);  assert(w == v[2]);
  h.oriented_end(v[4],v[7],w);  assert(w == v[3]);
  h.next_along_sc(v[4],v[7],w); assert(w == v[3]);
  assert( !h.next_along_sc(v[7],v[4],w));

  h.insert_constraint(v[1],v[8]);
  h.split_constraint(v[1],v[8],v[5]);
  assert(h.is_subconstrained_edge(v[1],v[5]));
  assert(h.is_constrained_edge(v[1],v[2]));
  H_constraint_list hcl;
  h.enclosing_constraints(v[1], v[5], hcl);
  assert(hcl.size() == 2);

  //h.print();
  // result should be 152 374 158
  assert(h.number_of_enclosing_constraints(v[1],v[5]) == 2);
  H_context_iterator co_it = h.contexts_begin(v[1],v[5]);
  assert(co_it->number_of_vertices() == 3);
  vit = co_it->vertices_begin();
  assert ( *vit++ == v[1] &&
           *vit++ == v[5] &&
           *vit++ == v[2]);
  co_it++;
  vit = co_it->vertices_begin();
  assert ( *vit++ == v[1] &&
           *vit++ == v[5] &&
           *vit++ == v[8]);
  co_it++;
  assert (co_it  == h.contexts_end(v[1],v[5]));

  //test clear() and copy() and swap()
  Hierarchy ch(h);
  assert( ch.number_of_constraints() == h.number_of_constraints());
  assert( ch.number_of_subconstraints() == h.number_of_subconstraints());
  ch.clear();
  ch = h;
  assert( ch.number_of_constraints() == h.number_of_constraints());
  assert( ch.number_of_subconstraints() == h.number_of_subconstraints());
  ch.clear();
  std::size_t nc = h.number_of_constraints();
  std::size_t nsc = h.number_of_subconstraints();
  ch.swap(h);
  assert (ch.number_of_constraints() == nc);
  assert( ch.number_of_subconstraints() == nsc);
  assert (h.number_of_constraints() == 0);
  assert( h.number_of_subconstraints() == 0);
  ch.swap(h);

  //test remove constraint
  //h.print();
  h.remove_constraint(v[1],v[2]);
  // h.print();
  h.remove_constraint(v[3],v[4]);
  // h.print();

  return;
}




