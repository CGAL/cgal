#ifndef CGAL_SM_DECORATOR_H 
#define CGAL_SM_DECORATOR_H

#include <CGAL/basic.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>

CGAL_BEGIN_NAMESPACE

/*{\Moptions print_title=yes }*/ 
/*{\Moptions outfile=SM_decorator.man }*/
/*{\Manpage {SM_decorator}{Sphere_map,Kernel}
{Topological sphere map decorator}{D}}*/

template <typename Sphere_map_, typename Kernel_>
class SM_decorator : public SM_const_decorator<Sphere_map_,Kernel_>
{ 
public:
typedef SM_const_decorator<Sphere_map_,Kernel_> Base;
typedef SM_decorator<Sphere_map_,Kernel_> Self;

/*{\Mdefinition ...}*/

/*{\Mtypes 5}*/

typedef Sphere_map_ Sphere_map;
typedef Kernel_ Kernel;

typedef typename Kernel_::Sphere_point Sphere_point;
/*{\Mtypemember embedding vertices.}*/

typedef typename Kernel_::Sphere_segment Sphere_segment;
/*{\Mtypemember embedding edges.}*/

typedef typename Kernel_::Sphere_circle Sphere_circle;
/*{\Mtypemember embedding loops.}*/

typedef typename Sphere_map_::Mark   Mark;
/*{\Mtypemember attributes of objects (vertices, edges, faces).}*/

typedef size_t Size_type;
/*{\Mtypemember size type.}*/

enum { BEFORE = -1, AFTER = 1 };
/*{\Menum insertion order labels.}*/

typedef void*  GenPtr;

#define USING(t) typedef typename Sphere_map_::t t
USING(Vertex); 
USING(Vertex_const_handle);
USING(Vertex_const_iterator);
USING(Halfedge); 
USING(Halfedge_const_handle); 
USING(Halfedge_const_iterator);
USING(Face);
USING(Face_const_handle);
USING(Face_const_iterator);
USING(Halfloop);
USING(Halfloop_const_handle);
USING(Halfloop_const_iterator);
USING(Vertex_handle);
USING(Vertex_iterator);
USING(Halfedge_handle); 
USING(Halfedge_iterator);
USING(Face_handle);
USING(Face_iterator);
USING(Halfloop_handle);
USING(Halfloop_iterator);
USING(Object_handle);
#undef USING
#define USING(t) typedef typename Base::t t
USING(Halfedge_around_vertex_const_circulator);
USING(Halfedge_around_face_const_circulator);
USING(Face_cycle_const_iterator);
#undef USING


/*{\Mtext Local types are handles, iterators and circulators of the
following kind: |Vertex_handle|, |Vertex_iterator|, |Halfedge_handle|,
|Halfedge_iterator|, |Halfloop_handle|, |Halfloop_iterator|,
|Face_handle|, |Face_iterator|.  Additionally the following
circulators are defined.}*/

typedef CircFromIt<
        Halfedge_iterator, 
        move_edge_around_vertex<Halfedge_iterator> > 
        Halfedge_around_vertex_circulator;
/*{\Mtypemember circulating the adjacency list of an vertex |v|.}*/

typedef CircFromIt<
        Halfedge_iterator, 
        move_edge_around_face<Halfedge_iterator> > 
        Halfedge_around_face_circulator;
/*{\Mtypemember circulating the face cycle of an face |f|.}*/

typedef typename Sphere_map_::Face_cycle_iterator 
  Face_cycle_iterator;
/*{\Mtypemember iterating all face cycles of an face |f|.
The iterator has method |bool is_vertex()|, |bool is_halfedge()|,
|bool is_halfloop()|, and can be converted to the corresponding
handles |Vertex_handle|, |Halfedge_handle|, or 
|Halfloop_handle|.}*/

public:

/*{\Mcreation 3}*/
SM_decorator() : Base() {}
SM_decorator(const Self& D) : Base(D) {}

SM_decorator(Sphere_map_& M) : Base(M) {}
/*{\Mcreate constructs a plane map decorator exploring |M|.}*/
 

/*{\Moperations 4 4}*/

Sphere_map& sm() { return *psm_; }

Vertex_const_handle source(Halfedge_const_handle e) const
{ return e->source_; }
Vertex_const_handle target(Halfedge_const_handle e) const
{ return e->twin_->source_; }
Halfedge_const_handle twin(Halfedge_const_handle e) const 
{ return e->twin_; }
Halfloop_const_handle twin(Halfloop_const_handle l) const 
{ return l->twin_; }

Vertex_handle source(Halfedge_handle e) const
/*{\Mop returns the source of |e|.}*/
{ return e->source_; }

Vertex_handle target(Halfedge_handle e) const
/*{\Mop returns the target of |e|.}*/
{ return e->twin_->source_; }

Halfedge_handle twin(Halfedge_handle e) const 
/*{\Mop returns the twin of |e|.}*/
{ return e->twin_; }

Halfloop_handle twin(Halfloop_handle l) const 
/*{\Mop returns the twin of |l|.}*/
{ return l->twin_; }

bool is_isolated(Vertex_const_handle v) const
{ return (v->edge_ == Halfedge_handle()); }

bool is_isolated(Vertex_handle v) const
/*{\Mop returns |true| when |v| is linked to the interior of a face.}*/
{ return (v->edge_ == Halfedge_handle()); }

Halfedge_const_handle first_out_edge(Vertex_const_handle v) const
{ return v->edge_; }
Halfedge_const_handle last_out_edge(Vertex_const_handle v) const
{ return cap(v->edge_); }

Halfedge_handle first_out_edge(Vertex_handle v) const
/*{\Mop returns one edge with source |v|. It's the starting point for
  the circular iteration over the edges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->edge_; }

Halfedge_handle last_out_edge(Vertex_handle v) const
/*{\Mop returns one edge with source |v|. \precond |!is_isolated(v)|.}*/
{ return cap(v->edge_); }

Halfedge_const_handle cyclic_adj_succ(Halfedge_const_handle e) const
{ return e->prev_->twin_; }
Halfedge_const_handle cyclic_adj_pred(Halfedge_const_handle e) const
{ return e->twin_->next_; }

Halfedge_handle cyclic_adj_succ(Halfedge_handle e) const
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->prev_->twin_; }

Halfedge_handle cyclic_adj_pred(Halfedge_handle e) const
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->twin_->next_; }

Halfedge_const_handle next(Halfedge_const_handle e) const
{ return e->next_; }
Halfedge_const_handle previous(Halfedge_const_handle e) const
{ return e->prev_; }

Halfedge_handle next(Halfedge_handle e) const
/*{\Mop returns the next edge in the face cycle containing |e|.}*/
{ return e->next_; }

Halfedge_handle previous(Halfedge_handle e) const
/*{\Mop returns the previous edge in the face cycle containing |e|.}*/
{ return e->prev_; }

Face_const_handle face(Halfedge_const_handle e) const
{ return e->face_; }
Face_const_handle face(Halfloop_const_handle l) const
{ return l->face_; }
Face_const_handle face(Vertex_const_handle v) const
{ return v->face_; }

Face_handle face(Halfedge_handle e) const
/*{\Mop returns the face incident to |e|.}*/
{ return e->face_; }

Face_handle face(Halfloop_handle l) const
/*{\Mop returns the face incident to |l|.}*/
{ return l->face_; }

Face_handle face(Vertex_handle v) const
/*{\Mop returns the face incident to |v|. \precond |is_isolated(v)|.}*/
{ return v->face_; }

/*{\Mtext \headerline{Iteration} \setopdims{3.3cm}{0cm}}*/
  
Vertex_iterator vertices_begin() const
{ return psm_->vertices_begin(); }
Vertex_iterator vertices_end() const
{ return psm_->vertices_end(); }
Halfedge_iterator halfedges_begin() const
{ return psm_->halfedges_begin(); }
Halfedge_iterator halfedges_end() const
{ return psm_->halfedges_end(); }
Face_iterator faces_begin() const
{ return psm_->faces_begin(); }
Face_iterator faces_end() const
{ return psm_->faces_end(); }
Halfloop_iterator halfloops_begin() const
{ return psm_->halfloops_begin(); }
Halfloop_iterator halfloops_end() const
{ return psm_->halfloops_end(); }

Halfloop_handle halfloop() const
{ return psm_->loops_; }

Halfedge_around_vertex_const_circulator 
  out_edges(Vertex_const_handle v) const
{ return Halfedge_around_vertex_circulator(first_out_edge(v)); }
Face_cycle_const_iterator face_cycles_begin(Face_const_handle f) const
{ return f->bounday_.begin(); }
Face_cycle_const_iterator face_cycles_end(Face_const_handle f) const
{ return f->boundary_.end(); }

Face_cycle_iterator face_cycles_begin(Face_handle f) const
/*{\Mop returns an iterator for all bounding face cycles of |f|.
The iterator is is convertable to |Vertex_handle|, 
|Halfloop_handle|, or |Halfedge_handle|.}*/
{ return f->boundary_.begin(); }

Face_cycle_iterator face_cycles_end(Face_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->boundary_.end(); }

Halfedge_around_vertex_circulator 
  out_edges(Vertex_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.
\precond the adjacency list is not empty.}*/
{ return Halfedge_around_vertex_circulator(first_out_edge(v)); }

/*{\Mtext \headerline{Update Operations}}*/

void clear() const
/*{\Mop reintializes |P| to the empty map.}*/
{ psm_->clear(); }

bool is_closed_at_source(Halfedge_handle e) const
/*{\Mop returns |true| when |prev(e) == twin(e)|.}*/
{ return previous(e) == twin(e); }

bool is_closed_at_target(Halfedge_handle e) const
/*{\Mop returns |true| when |next(e) == twin(e)|.}*/
{ return next(e) == twin(e); }

Halfedge_handle cas(Halfedge_handle e) const 
{ return cyclic_adj_succ(e); } 
Halfedge_handle cap(Halfedge_handle e) const
{ return cyclic_adj_pred(e); }

Vertex_handle new_vertex(const Sphere_point& p) const
/*{\Mop creates a new vertex.}*/
{ return psm_->new_vertex(p); }

Halfedge_handle new_edge_pair() const
/*{\Xop creates a new edge pair. No connectivity is provided.}*/
{ return psm_->new_halfedge_pair(); }

Halfloop_handle new_loop_pair() const
/*{\Mop creates a new loop pair.
\precond No sloop pair exists in the local graph.}*/ 
{ CGAL_nef_assertion(!has_loop());
  return psm_->new_halfloop_pair(); }

Face_handle new_face() const
/*{\Mop creates a new face.}*/
{ return psm_->new_face(); }

void delete_vertex_only(Vertex_handle v) const
/*{\Mop deletes |v| without any connectivity update.}*/
{ psm_->delete_vertex(v); }

void delete_edge_pair_only(Halfedge_handle e) const
/*{\Mop deletes |e| and its twin without any connectivity update.}*/
{ psm_->delete_halfedge_pair(e); }

void delete_halfedge_only(Halfedge_handle e) const
/*{\Mop deletes |e| without its twin and without any connectivity update.}*/
{ psm_->delete_halfedge(e); }

void delete_face_only(Face_handle f) const
/*{\Mop deletes |f| without any connectivity update.}*/
{ psm_->delete_face(f); }

void delete_loop_only() const
/*{\Mop deletes the loop and its twin without any connectivity update.}*/ 
{ CGAL_nef_assertion(psm_->loops_);
  psm_->delete_halfloop_pair(psm_->loops_); }

template <typename H>
bool is_boundary_object(H h) const
{ return psm_->is_boundary_object(h); }

template <typename H>
void store_boundary_object(H h, Face_handle f) const
{ f->boundary_.push_back(make_object(h));
  psm_->store_boundary_item(h, --(psm_->face_cycles_end(f)));
}

template <typename H>
void undo_boundary_object(H h, Face_handle f) const
{ CGAL_nef_assertion(psm_->is_boundary_object(h));
  Face_cycle_iterator it = psm_->boundary_item(h);
  psm_->undef_boundary_item(h);
  f->boundary_.erase(it);
}

void link_as_face_cycle(Halfedge_handle e, Face_handle f) const
/*{\Mop creates a new face cycle of |f| and 
   makes |e| the entry point of it.}*/
{
  Halfedge_around_face_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) hfc->face_ = f;
  store_boundary_object(e,f);
} 

void link_as_loop(Halfloop_handle l, Face_handle f) const
/*{\Mop creates a new trivial face cycle of |f| and 
   makes |l| the singular object of it.}*/
{ store_boundary_object(l,f); l->face_ = f; } 

void link_as_isolated_vertex(Vertex_handle v, Face_handle f) const
/*{\Mop creates a new trivial face cycle of |f|.
   (makes |v| an isolated vertex within |f|).}*/
{ store_boundary_object(v,f); v->face_ = f; } 

void unlink_as_face_cycle(Halfedge_handle e) const
/*{\Mop removes the face cycle defined by |e| from |face(e)|.
    Does not update the face links of the corresponding face cycle
    edges. \precond |e| is the entry object of the face cycle.}*/
{ undo_boundary_object(e,face(e)); }
  
void unlink_as_loop(Halfloop_handle l) const
/*{\Mop removes the trivial face cycle defined by |l| from
   |face(l)|. Does not update |l|'s face link.}*/
{ undo_boundary_object(l,face(l)); }

void unlink_as_isolated_vertex(Vertex_handle v) const
/*{\Mop removes the trivial face cycle defined by |v| from
   |face(v)|. Does not update |v|'s face link.
   \precond |v| is a trivial face cycle of |face(v)|.}*/
{ undo_boundary_object(v,face(v)); }

void clear_face_cycle_entries(Face_handle f) const
{ psm_->reset_object_list(f->boundary_);
  // removes entries of list and the hashed membership
}

Halfedge_handle new_edge_pair(Vertex_handle v1, 
			      Vertex_handle v2) const
/*{\Mop creates a new pair of edges |(e1,e2)| representing |(v1,v2)| 
  by appending the |ei| to |A(vi)| $(i=1,2)$.}*/
{ Halfedge_handle e1 = new_edge_pair();
  Halfedge_handle e2 = twin(e1);
  if ( ! is_isolated(v1) ) 
    set_adjacency_at_source_between(cap(first_out_edge(v1)),e1,
                                    first_out_edge(v1));
  else
    close_tip_at_source(e1,v1);
  if ( ! is_isolated(v2) )
    set_adjacency_at_source_between(cap(first_out_edge(v2)),e2,
                                    first_out_edge(v2));
  else 
    close_tip_at_source(e2,v2);
  return e1;
}


Halfedge_handle new_edge_pair(Halfedge_handle e1, 
			      Halfedge_handle e2,
			      int pos1 = AFTER, int pos2 = AFTER) const
/*{\Mop creates a new pair of edges |(es1,es2)| representing the uedge
  |\{source(e1),source(e2)\}| by inserting the |esi| before or after |ei| 
  into the cyclic adjacency list of |source(ei)| depending on |posi| 
  $(i=1,2)$ from |\Mname::BEFORE|, |\Mname::AFTER|.}*/
{ 
  Halfedge_handle er = new_edge_pair();
  Halfedge_handle ero = twin(er);
  if (pos1 < 0) { // before e1
    set_adjacency_at_source_between(cap(e1),er,e1);
    if ( e1 == first_out_edge(source(e1)) )
      set_first_out_edge(source(e1),er);
  } else { // after e1
    set_adjacency_at_source_between(e1,er,cas(e1));
  }
  if (pos2 < 0) { // before e2
    set_adjacency_at_source_between(cap(e2),ero,e2);
    if ( e2 == first_out_edge(source(e2)) )
      set_first_out_edge(source(e2),ero);
  } else { // after e2
    set_adjacency_at_source_between(e2,ero,cas(e2));
  }
  return er;
}

Halfedge_handle new_edge_pair(Halfedge_handle e, Vertex_handle v,
                               int pos = AFTER) const
/*{\Mop creates a new pair of edges  |(e1,e2)| representing the uedge
  |\{source(e),v\}| by inserting |e1| before or after |e| into cyclic 
  adjacency list of |source(e)| depending on |pos| from |\Mname::BEFORE|,
  |\Mname::AFTER| and appending |e2| at |A(v)|.}*/
{
  Halfedge_handle e_new = new_edge_pair();
  Halfedge_handle e_opp = twin(e_new);
  if (pos < 0) { // before e
    set_adjacency_at_source_between(cap(e),e_new,e);
    if ( e == first_out_edge(source(e)) )
      set_first_out_edge(source(e),e_new);
  } else  // after e
    set_adjacency_at_source_between(e,e_new,cas(e));
  
  if ( ! is_isolated(v) ) {
    Halfedge_handle e_first = first_out_edge(v);
    set_adjacency_at_source_between(cap(e_first),e_opp,e_first);
  } else
    close_tip_at_source(e_opp,v);
  return e_new;
}


Halfedge_handle new_edge_pair(Vertex_handle v, Halfedge_handle e,
                           int pos = AFTER) const
/*{\Mop symmetric to the previous one.}*/
{ return twin(new_edge_pair(e,v,pos)); }

void delete_edge_pair(Halfedge_handle e) const
/*{\Mop deletes |e| and its twin and maintains the adjacency at its source 
        and its target.}*/
{ remove_from_adj_list_at_source(e);
  remove_from_adj_list_at_source(twin(e));
  delete_edge_pair_only(e);
}

void delete_vertex(Vertex_handle v) const
/*{\Mop deletes |v| and all outgoing edges |A(v)| as well as their twins. 
   Updates the adjacency at the targets of the edges in |A(v)|.}*/
{ 
  if ( ! is_isolated(v) ) {
    Halfedge_handle e = first_out_edge(v);
    while ( e != cap(e) ) 
      delete_edge_pair(cap(e));  
    delete_edge_pair(e); 
  }
  delete_vertex_only(v);
}

void delete_face(Face_handle f) const
/*{\Mop deletes the face |f| without consideration of topological 
   linkage.}*/
{ clear_face_cycle_entries(f); delete_face_only(f); }

bool has_outdeg_two(Vertex_handle v) const
/*{\Mop return true when |v| has outdegree two.}*/
// does this work for looping edges?
{ if (is_isolated(v)) return false;
  Halfedge_handle e1 = first_out_edge(v);
  Halfedge_handle e2 = last_out_edge(v);
  return (e1!=e2 && e2==cas(e1));
}

void link_as_prev_next_pair(Halfedge_handle e1, Halfedge_handle e2) const 
/*{\Xop makes |e1| and |e2| adjacent in the face cycle 
   $\ldots-|e1-e2|-\ldots$.
   Afterwards |e1 = previous(e2)| and |e2 = next(e1)|.}*/
{ e1->next_ = e2; e2->prev_ = e1; }

void merge_edge_pairs_at_target(Halfedge_handle e) const
/*{\Mop merges the edge pairs at |v = target(e)|. |e| and |twin(e)| 
  are preserved, |next(e)|, |twin(next(e))| and |v| are deleted
  in the merger. \precond |v| has outdegree two. The adjacency at 
  |source(e)| and |target(next(e))| is kept consistent.
  If |next(e)| was entry point of |face(e)| then |e| takes this role.
  The same holds for |twin(next(e))| and |face(twin(e))|.}*/
{
  TRACEN("merge_edge_pairs_at_target "<<PH(e));
  Halfedge_handle en = next(e), eno = twin(en), enn, enno,
               eo = twin(e) ;
  if ( is_closed_at_target(en) ) { enn = eo; enno=e; }
  else { enn = next(en), enno = previous(eno); }
  Vertex_handle v = target(e), vn = target(en);
  CGAL_nef_assertion(has_outdeg_two(v));
  Face_handle f1 = face(en), f2 = face(eno);
  // transfer the opposite face cycles e-en-enn to e-enn
  if ( enn != eno ) {
    link_as_prev_next_pair(e,enn);
    link_as_prev_next_pair(enno,eo);
  } else {
    link_as_prev_next_pair(e,eo);
  }
  // set vertex of e and deal with vertex-halfedge incidence
  eo->source_ = vn;

  if ( first_out_edge(vn) == eno ) set_first_out_edge(vn,eo);
  if ( is_boundary_object(en) ) 
  { undo_boundary_object(en,f1); store_boundary_object(e,f1); }
  if ( is_boundary_object(eno) )
  { undo_boundary_object(eno,f2); store_boundary_object(eo,f2); }
  delete_vertex_only(v);
  delete_edge_pair_only(en);
  TRACEN("END "<<PH(previous(e))<<PH(e)<<PH(next(e)));
}

void convert_edge_to_loop(Halfedge_handle e) const
/*{\Mop converts the edge at |v = target(e)| to the unique
  loop |l| of |\Mvar|. |e|, |twin(e)| and |v| are deleted
  in the conversion. \precond there was no loop in |\Mvar|. 
  As |e| was entry point of |face(e)| then |l| takes this role.}*/
{ TRACEN("convert_edge_to_loop "<<PH(e));
  CGAL_nef_assertion( !has_loop() );
  Halfloop_handle l = new_loop_pair();
  Vertex_handle v = target(e);
  Face_handle f1 = face(e), f2 = face(twin(e));
  CGAL_nef_assertion( is_boundary_object(e) && 
                  is_boundary_object(twin(e)) );
  undo_boundary_object(e,f1); undo_boundary_object(twin(e),f2);
  link_as_loop(l,f1), link_as_loop(twin(l),f2);
  circle(l) = circle(e); circle(twin(l)) = circle(twin(e));
  mark(l) = mark(e);
  delete_vertex_only(v);
  delete_edge_pair_only(e);
}

void flip_diagonal(Halfedge_handle e) const
{ Halfedge_handle r = twin(e);
  Halfedge_handle en = next(e), enn= next(en);
  Halfedge_handle rn = next(r), rnn= next(rn);
  TRACEN(PH(e)<<PH(en)<<PH(enn));
  TRACEN(PH(r)<<PH(rn)<<PH(rnn));
  CGAL_nef_assertion( next(enn)==e && next(rnn)==r );
  remove_from_adj_list_at_source(e);
  remove_from_adj_list_at_source(r);
  set_adjacency_at_source_between(enn,e,twin(en));
  set_adjacency_at_source_between(rnn,r,twin(rn));
}
     

/*{\Mtext \headerline{Incomplete topological update primitives}}*/

Halfedge_handle new_edge_pair_at_source
  (Halfedge_handle e, int pos = AFTER) const
/*{\Xop creates a new pair of edges  |(e1,e2)| representing |(source(e),())| 
  by inserting |e1| before or after |e| into cyclic adjacency list of
  |source(e)| depending on |pos| from |\Mname::BEFORE, \Mname::AFTER|.}*/
{ Halfedge_handle e_new = new_edge_pair();
  if (pos < 0) set_adjacency_at_source_between(cap(e),e_new,e);
  else         set_adjacency_at_source_between(e,e_new,cas(e));
  return e_new;
}

Halfedge_handle new_edge_pair_at_source
  (Vertex_handle v, int pos = AFTER) const
/*{\Mop creates a new pair of edges  |(e1,e2)| representing |(v,())| 
  by inserting |e1| at the beginning (BEFORE) or end (AFTER)
  of adjacency list of |v|.}*/
{ Halfedge_handle e1 = new_edge_pair();
  if ( ! is_isolated(v) ) {
    Halfedge_handle ef = first_out_edge(v);
    if (pos < 0) { // before e1
      set_adjacency_at_source_between(cap(ef),e1,ef);
      set_first_out_edge(v,e1);
    } else {
      set_adjacency_at_source_between(ef,e1,cas(ef));
    }
  } else
    close_tip_at_source(e1,v);
  return e1;
}

void delete_edge_pair_at_source(Halfedge_handle e) const
/*{\Mop deletes |e| and its twin and maintains the adjacency at its 
  source.}*/
{ remove_from_adj_list_at_source(e);
  delete_edge_pair_only(e);
}

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) const
/*{\Mop makes |v| the target of |e| appends |twin(e)| to the adjacency list
   of |v|.}*/
{ if ( !is_isolated(v) ) 
    set_adjacency_at_source_between(cap(first_out_edge(v)),twin(e),
      first_out_edge(v));
  else
    close_tip_at_target(e,v);
}

void link_as_source_of(Halfedge_handle e, Vertex_handle v) const
/*{\Mop makes |source(e) = v| and sets |e| as the first
        out edge if |v| was isolated before.}*/
{ e->source_ = v;
  if (v->edge_ == Halfedge_handle()) v->edge_ = e; }

void link_as_target_of(Halfedge_handle e, Vertex_handle v) const
/*{\Mop makes |target(e) = v| and sets |e| as the first
  in edge if |v| was isolated before.}*/
{ link_as_source_of(twin(e),v); }

void set_adjacency_at_source_between(Halfedge_handle e, Halfedge_handle en) 
  const 
/*{\Mop makes |e| and |en| neigbors in the cyclic ordered adjacency list 
    around |v=source(e)|. \precond |source(e)==source(en)|.}*/
{ CGAL_nef_assertion(source(e)==source(en));
  link_as_prev_next_pair(twin(en),e);
}

void set_adjacency_at_source_between(Halfedge_handle e1, 
                                     Halfedge_handle e_between, 
                                     Halfedge_handle e2) const 
/*{\Mop inserts |e_between| into the adjacency list around |source(e1)| 
  between |e1| and |e2| and makes |source(e1)| the source of |e_between|. 
  \precond |source(e1)==source(e2)|.}*/
{ 
  e_between->source_ = source(e1);
  set_adjacency_at_source_between(e1,e_between);
  set_adjacency_at_source_between(e_between,e2);
}

void close_tip_at_source(Halfedge_handle e, Vertex_handle v) const 
/*{\Mop sets |v| as source of |e| and closes the tip by setting the 
  corresponding pointers such that |prev(e) == twin(e)| and
  |next(twin(e)) == e|.}*/
{ link_as_source_of(e,v); 
  link_as_prev_next_pair(twin(e),e); }

void close_tip_at_target(Halfedge_handle e, Vertex_handle v) const 
/*{\Mop sets |v| as target of |e| and closes the tip by setting the 
  corresponding pointers such that |prev(twin(e)) == e| and 
  |next(e) == twin(e)|.}*/
{ link_as_target_of(e,v);
  link_as_prev_next_pair(e,twin(e)); }


void remove_from_adj_list_at_source(Halfedge_handle e) const
/*{\Mop removes a halfedge pair |(e,twin(e)| from the adjacency list
of |source(e)|. Afterwards |next(prev(e))==next(twin(e))| and
|first_out_edge(source(e))| is valid if |degree(source(v))>1| before
the operation.}*/
{
  Vertex_handle v = source(e);
  if ( is_closed_at_source(e) ) { // last outgoing
    v->edge_ = Halfedge_handle();
  } else {
    if (e == first_out_edge(v)) v->edge_ = cap(e);
    set_adjacency_at_source_between(cap(e),cas(e));
  }
}

void set_face(Halfedge_handle e, Face_handle f) const
{ e->face_ = f; }
void set_face(Vertex_handle v, Face_handle f) const
{ v->face_ = f; }
void set_first_out_edge(Vertex_handle v, Halfedge_handle e) const
{ v->edge_ = e; }
void set_prev(Halfedge_handle e, Halfedge_handle ep) const
{ e->prev_ = ep; }
void set_next(Halfedge_handle e, Halfedge_handle en) const
{ e->next_ = en; }
void set_source(Halfedge_handle e, Vertex_handle v) const
{ e->source_ = v; }



/*{\Mtext \headerline{Associated Information}\restoreopdims}*/

Sphere_point& point(Vertex_handle v) const
/*{\Mop returns the embedding of |v|.}*/
{ return v->point_; }

Sphere_circle& circle(Halfedge_handle e) const
/*{\Mop returns the plane supporting |e|.}*/
{ return e->circle_; }

Sphere_circle& circle(Halfloop_handle l) const
/*{\Mop returns the plane supporting |e|.}*/
{ return l->circle_; }

Mark& mark(Vertex_handle v) const
/*{\Mop returns the mark of |v|.}*/
{ return v->mark_; }

Mark& mark(Halfedge_handle e) const
/*{\Mop returns the mark of |e|.}*/
{ return ( &*e < &*twin(e) ) ? e->mark_ : twin(e)->mark_; }

Mark& mark(Halfloop_handle l) const
/*{\Mop returns the mark of |l|.}*/
{ return ( &*l < &*twin(l) ) ? l->mark_ : twin(l)->mark_; }

Mark& mark(Face_handle f) const
/*{\Mop returns the mark of |f|.}*/
{ return f->mark_; }

void unify_marks(Halfedge_handle e) const
{ if ( &*e < &*twin(e) ) twin(e)->mark_ = e->mark_; 
  else                   e->mark_ = twin(e)->mark_;
}

const Sphere_point& point(Vertex_const_handle v) const
{ return v->point_; }
const Sphere_circle& circle(Halfedge_const_handle e) const
{ return e->circle_; }
const Sphere_circle& circle(Halfloop_const_handle l) const
{ return l->circle_; }

const Mark& mark(Vertex_const_handle v) const
{ return v->mark_; }
const Mark& mark(Halfedge_const_handle e) const
{ return ( &*e < &*twin(e) ) ? e->mark_ : twin(e)->mark_; }
const Mark& mark(Halfloop_const_handle l) const
{ return ( &*l < &*twin(l) ) ? l->mark_ : twin(l)->mark_; }
const Mark& mark(Face_const_handle f) const
{ return f->mark_; }

void unify_tmp_marks(Halfedge_handle e) const
{ if (&*e < &*twin(e)) twin(e)->mark_ = e->mark_; 
  else e->mark_ = twin(e)->mark_; }

void set_marks_in_face_cycle(Halfedge_handle e, Mark m) const
{ Halfedge_around_face_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) mark(hfc) = mark(target(hfc)) = m;
}

Mark& mark_of_halfsphere(int i) const
{ CGAL_nef_assertion(i);
  if (i<0) return psm_->m_neg_; 
  return psm_->m_pos_; }

GenPtr& info(Vertex_handle v) const
{ return v->info_; }
GenPtr& info(Halfedge_handle e) const
{ return e->info_; }
GenPtr& info(Halfloop_handle l) const
{ return l->info_; }
GenPtr& info(Face_handle f) const
{ return f->info_; }

const GenPtr& info(Vertex_const_handle v) const
{ return v->info_; }
const GenPtr& info(Halfedge_const_handle e) const
{ return e->info_; }
const GenPtr& info(Halfloop_const_handle l) const
{ return l->info_; }
const GenPtr& info(Face_const_handle f) const
{ return f->info_; }

 
/*{\Mtext \headerline{Iteration}}*/
/*{\Mtext The list of all objects can be accessed via iterator ranges.
For comfortable iteration we also provide iterations macros. 
The iterator range access operations are of the following kind:\\
|Vertex_iterator vertices_begin()/vertices_end()|\\
|Halfedge_iterator halfedges_begin()/halfedges_end()|\\
|Face_iterator faces_begin()/faces_end()|

The macros are then |CGAL_forall_vertices_of(v,V)|,
|CGAL_forall_halfedges_of(e,V)|, |CGAL_forall_edges_of(e,V)|, 
|CGAL_forall_faces_of(f,V)|, |CGAL_forall_face_cycles_of(fc,F)|.}*/


}; // SM_decorator


CGAL_END_NAMESPACE
#endif // CGAL_SM_DECORATOR_H


