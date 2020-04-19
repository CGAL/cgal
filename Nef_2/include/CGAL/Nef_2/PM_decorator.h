// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_PM_DECORATOR_H
#define CGAL_PM_DECORATOR_H

#include <CGAL/license/Nef_2.h>

#include <CGAL/Nef_2/PM_const_decorator.h>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif
#include <CGAL/Unique_hash_map.h>
#include <vector>

namespace CGAL {

/*{\Moptions outfile=PMDecorator.man }*/
/*{\Msubst
PM_decorator#PMDecorator
PM_const_decorator#PMConstDecorator
}*/
/*{\Manpage {PMDecorator}{}{Plane map manipulation}{D}}*/

template <typename HDS>
class PM_decorator : public PM_const_decorator<HDS>
{
/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to examine and modify a plane map. |\Mvar| inherits from
|PM_const_decorator| but provides additional manipulation operations.}*/

/*{\Mgeneralization PM_const_decorator}*/

public:
/*{\Mtypes 3}*/
  typedef PM_decorator<HDS>        Self;
  typedef PM_const_decorator<HDS>  Base;
  typedef HDS                       Plane_map;
  typedef typename HDS::Vertex_base Vertex_base;
  typedef typename HDS::Halfedge_base Halfedge_base;
  typedef typename HDS::Face_base Face_base;
  typedef typename HDS::Vertex Vertex;
  typedef typename HDS::Halfedge Halfedge;
  typedef typename HDS::Face Face;
  typedef typename HDS::Vertex_handle Vertex_handle;
  typedef typename HDS::Vertex_iterator Vertex_iterator;
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename HDS::Halfedge_iterator Halfedge_iterator;
  typedef typename HDS::Face_handle Face_handle;
  typedef typename HDS::Face_iterator Face_iterator;
  typedef typename HDS::Vertex_const_handle Vertex_const_handle;
  typedef typename HDS::Halfedge_const_handle Halfedge_const_handle;
  typedef typename HDS::Face_const_handle Face_const_handle;
  typedef typename HDS::Vertex_const_iterator Vertex_const_iterator;
  typedef typename HDS::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename HDS::Face_const_iterator Face_const_iterator;
  typedef typename Base::Hole_const_iterator Hole_const_iterator ;
  typedef typename Base::Isolated_vertex_const_iterator Isolated_vertex_const_iterator;
  typedef typename Base::Point_const_iterator Point_const_iterator;
  typedef typename Base::Mark Mark;
typedef typename Base::Point Point;
typedef typename Base::GenPtr GenPtr;

/*{\Mtext Local types are handles, iterators and circulators of the following
kind: |Vertex_handle|, |Vertex_iterator|, |Halfedge_handle|,
|Halfedge_iterator|, |Face_handle|, |Face_iterator|.  Additionally the
following circulators are defined. The |circulators| can be constructed from
the corresponding halfedge handles or iterators.}*/

typedef CircFromIt<
        Halfedge_iterator,
        move_halfedge_around_vertex<Halfedge_iterator> >
        Halfedge_around_vertex_circulator;
/*{\Mtypemember circulating the outgoing halfedges in $A(v)$.}*/

typedef CircFromIt<
        Halfedge_iterator,
        move_halfedge_around_face<Halfedge_iterator> >
        Halfedge_around_face_circulator;
/*{\Mtypemember circulating the halfedges in the face cycle of a
face |f|.}*/

typedef typename Face::Hole_iterator Hole_iterator;
/*{\Mtypemember iterating all holes of a face |f|. The type is
convertible to |Halfedge_handle|.}*/

typedef typename Face::Isolated_vertex_iterator Isolated_vertex_iterator;
/*{\Mtypemember iterating all isolated vertices of a face |f|.
The type generalizes |Vertex_handle|.}*/


/* note: originally I had the mhavs, mhafs hardwired to Halfedge
   in this class scope. egcs 290.60 reacted with an internal compiler
   error; this recursive instatiation scheme works however!
   what a shitty world */

enum { BEFORE = -1, AFTER = 1 };
/*{\Menum insertion order labels.}*/




  using Base::phds;

/*{\Mcreation 3}*/

PM_decorator() : Base() {}
PM_decorator(const PM_decorator& D) :
  Base((PM_const_decorator<HDS>)D) {}
PM_decorator& operator=(const PM_decorator& D)
{ Base::phds = ((PM_const_decorator<HDS>)D).phds; return *this; }

PM_decorator(Plane_map& p) : Base(p) {}
/*{\Mcreate constructs a decorator working on |P|.}*/

#define BASE(t) { return Base::t; }
Hole_const_iterator  holes_begin(Face_const_handle f) const
BASE(holes_begin(f))
Hole_const_iterator  holes_end(Face_const_handle f) const
BASE(holes_end(f))
Isolated_vertex_const_iterator
isolated_vertices_begin(Face_const_handle f) const
BASE(isolated_vertices_begin(f))
Isolated_vertex_const_iterator
isolated_vertices_end(Face_const_handle f) const
BASE(isolated_vertices_end(f))
Vertex_const_handle source(Halfedge_const_handle e) const
BASE(source(e))
Vertex_const_handle target(Halfedge_const_handle e) const
BASE(target(e))
Halfedge_const_handle twin(Halfedge_const_handle e) const
BASE(twin(e))
bool is_isolated(Vertex_const_handle v) const
BASE(is_isolated(v))
Halfedge_const_handle first_out_edge(Vertex_const_handle v) const
BASE(first_out_edge(v))
Halfedge_const_handle last_out_edge(Vertex_const_handle v) const
BASE(last_out_edge(v))
Halfedge_const_handle cyclic_adj_succ(
  Halfedge_const_handle e) const
BASE(cyclic_adj_succ(e))
Halfedge_const_handle cyclic_adj_pred(
  Halfedge_const_handle e) const
BASE(cyclic_adj_pred(e))
Halfedge_const_handle next(Halfedge_const_handle e) const
BASE(next(e))
Halfedge_const_handle previous(Halfedge_const_handle e) const
BASE(previous(e))
Face_const_handle face(Halfedge_const_handle e) const
BASE(face(e))
Face_const_handle face(Vertex_const_handle v) const
BASE(face(v))
const Mark& mark(Vertex_const_handle v) const
BASE(mark(v))
const Mark& mark(Halfedge_const_handle e) const
BASE(mark(e))
const Mark& mark(Face_const_handle f) const
BASE(mark(f))
const Point& point(Vertex_const_handle v) const
BASE(point(v))
const GenPtr& info(Vertex_const_handle v) const
BASE(info(v))
const GenPtr& info(Halfedge_const_handle e) const
BASE(info(e))
const GenPtr& info(Face_const_handle f) const
BASE(info(f))
#undef BASE

/*{\Moperations 3 3}*/

Plane_map& plane_map() const
/*{\Mop returns the plane map decorated.}*/
{ return *Base::phds; }

void clear() const
/*{\Mop reinitializes |P| to the empty map.}*/
{ this->phds->clear(); }

Vertex_handle source(Halfedge_handle e) const
/*{\Mop returns the source of |e|.}*/
{ return e->opposite()->vertex(); }

Vertex_handle target(Halfedge_handle e) const
/*{\Mop returns the target of |e|.}*/
{ return e->vertex(); }

Halfedge_handle twin(Halfedge_handle e) const
/*{\Mop returns the twin of |e|.}*/
{ return e->opposite(); }

bool is_isolated(Vertex_handle v) const
/*{\Mop returns |true| iff |v| is linked to the interior of a face.
This is equivalent to the condition that $A(v) = \emptyset$.}*/
{ return v->is_isolated(); }

bool is_closed_at_source(Halfedge_handle e) const
/*{\Mop returns |true| when |prev(e) == twin(e)|.}*/
{ return e->prev() == e->opposite(); }

Halfedge_handle first_out_edge(Vertex_handle v) const
/*{\Mop returns a halfedge with source |v|. It's the starting point for
  the circular iteration over the halfedges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->halfedge()->opposite(); }

Halfedge_handle last_out_edge(Vertex_handle v) const
/*{\Mop returns a the halfedge with source |v| that is the last
  in the circular iteration before encountering |first_out_edge(v)|
  again. \precond |!is_isolated(v)|.}*/
{ return v->halfedge()->next(); }

Halfedge_handle cas(Halfedge_handle e) const
{ return e->prev()->opposite(); }

Halfedge_handle cap(Halfedge_handle e) const
{ return e->opposite()->next(); }

Halfedge_handle cyclic_adj_succ(Halfedge_handle e) const
{ return cas(e); }
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/

Halfedge_handle cyclic_adj_pred(Halfedge_handle e) const
{ return cap(e); }
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/

Halfedge_handle adj_succ_at_source(Halfedge_handle e) const
{ if (e==last_out_edge(source(e))) return Halfedge_handle();
  return cas(e); }
// the edge of source(e) is the first in the adj list

Halfedge_handle adj_pred_at_source(Halfedge_handle e) const
{ if (e==first_out_edge(source(e))) return Halfedge_handle();
  return cap(e); }
// the edge of source(e) is the first in the adj list

Halfedge_handle next(Halfedge_handle e) const
/*{\Mop returns the next edge in the face cycle containing |e|.}*/
{ return e->next(); }

Halfedge_handle previous(Halfedge_handle e) const
/*{\Mop returns the previous edge in the face cycle containing |e|.}*/
{ return e->prev(); }

Face_handle face(Halfedge_handle e) const
/*{\Mop returns the face incident to |e|.}*/
{ return e->face(); }

Face_handle face(Vertex_handle v) const
/*{\Mop returns the face incident to |v|.
   \precond |is_isolated(v)|.}*/
{ return v->face(); }

Halfedge_handle halfedge(Face_handle f) const
/*{\Mop returns a halfedge in the bounding face cycle of |f|
(|Halfedge_handle()| if there is no bounding face cycle).}*/
{ return f->halfedge(); }

Vertex_iterator   vertices_begin() const
{ return this->phds->vertices_begin(); }
Halfedge_iterator halfedges_begin() const
{ return this->phds->halfedges_begin(); }
Face_iterator     faces_begin() const
{ return this->phds->faces_begin(); }
Vertex_iterator   vertices_end() const
{ return this->phds->vertices_end(); }
Halfedge_iterator halfedges_end() const
{ return this->phds->halfedges_end(); }
Face_iterator     faces_end() const
{ return this->phds->faces_end(); }

/*{\Mtext \headerline{Iteration} \setopdims{6.5cm}{0cm}}*/

Halfedge_around_vertex_circulator
  out_edges(Vertex_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.}*/
{ return Halfedge_around_vertex_circulator(first_out_edge(v)); }

Halfedge_around_face_circulator
  face_cycle(Face_handle f) const
/*{\Mop returns a circulator for the outer face cycle of |f|.}*/
{ return Halfedge_around_face_circulator(f->halfedge()); }

Hole_iterator  holes_begin(Face_handle f) const
/*{\Mop returns an iterator for all holes in the interior of |f|.
   A |Hole_iterator| can be assigned to a
   |Halfedge_around_face_circulator|.}*/
{ return f->fc_begin(); }

Hole_iterator  holes_end(Face_handle f) const
/*{\Mop returns the past-the-end iterator of |f|.}*/
{ return f->fc_end(); }

Isolated_vertex_iterator isolated_vertices_begin(Face_handle f) const
/*{\Mop returns an iterator for all isolated vertices in the interior
of |f|.}*/
{ return f->iv_begin(); }

Isolated_vertex_iterator isolated_vertices_end(Face_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->iv_end(); }

/*{\Mtext \headerline{Update Operations} \restoreopdims}*/

Vertex_handle new_vertex(const Vertex_base& vb = Vertex_base()) const
/*{\Mop creates a new vertex.}*/
{ Vertex_handle v = this->phds->vertices_push_back(vb);
  v->set_halfedge(Halfedge_handle());
  return v;
}

Face_handle new_face(const Face_base& fb = Face_base()) const
/*{\Mop creates a new face.}*/
{ Face_handle f = this->phds->faces_push_back(fb);
  return f;
}


void link_as_outer_face_cycle(Face_handle f, Halfedge_handle e) const
/*{\Mop makes |e| the entry point of the outer face cycle of |f| and
makes |f| the face of all halfedges in the face cycle of |e|.}*/
{
  Halfedge_around_face_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) hfc->set_face(f);
  f->set_halfedge(e);
}

void link_as_hole(Face_handle f, Halfedge_handle e) const
/*{\Mop makes |e| the entry point of a hole face cycle of |f| and
    makes |f| the face of all halfedges in the face cycle of |e|.}*/
{
  Halfedge_around_face_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) hfc->set_face(f);
  f->store_fc(e);
}

void link_as_isolated_vertex(Face_handle f, Vertex_handle v) const
/*{\Mop makes |v| an isolated vertex within |f|.}*/
{  f->store_iv(v); v->set_face(f); }

void clear_face_cycle_entries(Face_handle f) const
/*{\Mop removes all isolated vertices and halfedges that
are entrie points into face cycles from the lists of |f|.}*/
{ f->clear_all_entries(); }


Halfedge_handle new_halfedge_pair(Vertex_handle v1, Vertex_handle v2,
                                  Halfedge_base hb = Halfedge_base()) const
/*{\Mop creates a new pair of edges |(e1,e2)| representing |(v1,v2)|
  by appending the |ei| to |A(vi)| $(i=1,2)$.}*/
{ Halfedge_handle e1 = this->phds->edges_push_back(hb,hb);
  Halfedge_handle e2 = e1->opposite();
  e1->set_face(Face_handle()); e2->set_face(Face_handle());
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


Halfedge_handle new_halfedge_pair(Halfedge_handle e1,
                                  Halfedge_handle e2,
                                  Halfedge_base hb = Halfedge_base(),
                                  int pos1 = AFTER, int pos2 = AFTER) const
/*{\Mop creates a new pair of edges |(h1,h2)| representing
  |(source(e1),source(e2))| by inserting the |hi| before or after |ei|
  into the cyclic adjacency list of |source(ei)| depending on
  |posi| $(i=1,2)$ from |\Mname::BEFORE, \Mname::AFTER|.}*/
{
  Halfedge_handle er = this->phds->edges_push_back(hb,hb);
  Halfedge_handle ero = er->opposite();
  er->set_face(Face_handle()); ero->set_face(Face_handle());
  if (pos1 < 0) { // before e1
    set_adjacency_at_source_between(cap(e1),er,e1);
    if ( e1 == first_out_edge(source(e1)) )
      make_first_out_edge(er); // added 22/8/00
  } else { // after e1
    set_adjacency_at_source_between(e1,er,cas(e1));
  }
  if (pos2 < 0) { // before e2
    set_adjacency_at_source_between(cap(e2),ero,e2);
    if ( e2 == first_out_edge(source(e2)) )
      make_first_out_edge(ero);
  } else { // after e2
    set_adjacency_at_source_between(e2,ero,cas(e2));
  }
  return er;
}

Halfedge_handle new_halfedge_pair(Halfedge_handle e, Vertex_handle v,
                                  Halfedge_base hb = Halfedge_base(),
                                  int pos = AFTER) const
/*{\Mop creates a new pair of edges  |(e1,e2)| representing |(source(e),v)|
  by inserting |e1| before or after |e| into the cyclic adjacency list of
  |source(e)| depending on |pos| from |\Mname::BEFORE, \Mname::AFTER|
  and appending |e2| to |A(v)|.}*/
{
  Halfedge_handle e_new = this->phds->edges_push_back(hb,hb);
  Halfedge_handle e_opp = e_new->opposite();
  e_new->set_face(Face_handle()); e_opp->set_face(Face_handle());

  if (pos < 0) { // before e
    set_adjacency_at_source_between(cap(e),e_new,e);
    if ( e == first_out_edge(source(e)) )
      make_first_out_edge(e_new);
  } else  // after e
    set_adjacency_at_source_between(e,e_new,cas(e));

  if ( ! is_isolated(v) ) {
    Halfedge_handle e_first = first_out_edge(v);
    set_adjacency_at_source_between(cap(e_first),e_opp,e_first);
  } else
    close_tip_at_source(e_opp,v);

  return e_new;
}


Halfedge_handle new_halfedge_pair(Vertex_handle v, Halfedge_handle e,
                                  Halfedge_base hb = Halfedge_base(),
                                  int pos = AFTER) const
/*{\Mop symmetric to the previous one.}*/
{ return new_halfedge_pair(e,v,hb,pos)->opposite(); }




void delete_halfedge_pair(Halfedge_handle e) const
/*{\Mop deletes |e| and its twin and updates the adjacency at its source
        and its target.}*/
{ remove_from_adj_list_at_source(e);
  remove_from_adj_list_at_source(e->opposite());
  this->phds->edges_erase(e);
}

void delete_vertex(Vertex_handle v) const
/*{\Mop deletes |v| and all outgoing edges |A(v)| as well as their twins.
        Updates the adjacency at the targets of the edges in |A(v)|.}*/
{
  if ( ! is_isolated(v) ) {
    Halfedge_handle e = first_out_edge(v);
    while ( e != cap(e) )
      delete_halfedge_pair(cap(e));
    delete_halfedge_pair(e);
  }
  this->phds->vertices_erase(v);
}

void delete_face(Face_handle f) const
/*{\Mop deletes the face |f| without consideration of topological linkage.}*/
{ this->phds->faces_erase(f); }


bool has_outdeg_two(Vertex_handle v) const
/*{\Mop return true when |v| has outdegree two.}*/
{ if (v->is_isolated()) return false;
  Halfedge_handle e1 = v->halfedge();
  Halfedge_handle e2 = e1->next()->opposite();
  return (e1!=e2 && e2->next()->opposite()==e1);
}

void merge_halfedge_pairs_at_target(Halfedge_handle e) const
/*{\Mop merges the halfedge pairs at |v = target(e)|. |e| and
  |twin(e)| are preserved, |next(e)|, |twin(next(e))| and |v| are deleted
  in the merger. \precond |v| has outdegree two. The adjacency at |source(e)|
  and |target(next(e))| is kept consistent.}*/
{
  CGAL_NEF_TRACEN("merge_halfedge_pairs_at_target "<<PE(e));
  Halfedge_handle eo = e->opposite(),
                  en = e->next(), eno = en->opposite(),
                  enn = en->next(), enno = eno->prev();
  Vertex_handle v = e->vertex(), vn = en->vertex();
  CGAL_assertion(has_outdeg_two(v));
  Face_handle f1 = en->face(), f2 = eno->face();
  // transfer the opposite face cycles e-en-enn to e-enn
  if ( enn != eno ) {
    e->set_next(enn); enn->set_prev(e);
    eo->set_prev(enno); enno->set_next(eo);
  } else {
    e->set_next(eo); eo->set_prev(e);
  }
  // set vertex of e and deal with vertex-halfedge incidence
  e->set_vertex(vn);
  if (vn->halfedge()==en) vn->set_halfedge(e);
  if (en->is_hole_entry())
  { f1->remove_fc(en); f1->store_fc(e); }
  if (eno->is_hole_entry())
  { f2->remove_fc(eno); f2->store_fc(eo); }
  if (f1->halfedge() == en) f1->set_halfedge(e);
  if (f2->halfedge() == eno) f2->set_halfedge(eo);
  this->phds->vertices_erase(v);
  this->phds->edges_erase(en);
}

void flip_diagonal(Halfedge_handle e) const
{ Halfedge_handle r = twin(e);
  Halfedge_handle en = e->next(), enn= en->next();
  Halfedge_handle rn = r->next(), rnn= rn->next();
  CGAL_assertion( enn->next()==e && rnn->next()==r );
  remove_from_adj_list_at_source(e);
  remove_from_adj_list_at_source(r);
  set_adjacency_at_source_between(enn,e,twin(en));
  set_adjacency_at_source_between(rnn,r,twin(rn));
}


/*{\Mtext \headerline{Incomplete topological update primitives}}*/

Halfedge_handle new_halfedge_pair_at_source
  (Halfedge_handle e, int pos = AFTER, Halfedge_base hb =
   Halfedge_base()) const
/*{\Xop creates a new pair of edges  |(e1,e2)| representing |(source(e),())|
  by inserting |e1| before or after |e| into cyclic adjacency list of
  |source(e)| depending on |pos| from |\Mname::BEFORE, \Mname::AFTER|.}*/
{
  Halfedge_handle e_new = this->phds->edges_push_back(hb,hb);
  if (pos < 0) // before e
    set_adjacency_at_source_between(cap(e),e_new,e);
  else  // after e
    set_adjacency_at_source_between(e,e_new,cas(e));
  return e_new;
}

Halfedge_handle new_halfedge_pair_at_source
  (Vertex_handle v, int pos = AFTER, Halfedge_base hb = Halfedge_base()) const
/*{\Mop creates a new pair of edges  |(e1,e2)| representing |(v,())|
  by inserting |e1| at the beginning (BEFORE) or end (AFTER)
  of adjacency list of |v|.}*/
{ Halfedge_handle e1 = this->phds->edges_push_back(hb,hb);
  Halfedge_handle e2 = e1->opposite();
  e1->set_face(Face_handle()); e2->set_face(Face_handle());
  if ( ! is_isolated(v) ) {
    Halfedge_handle ef = first_out_edge(v);
    set_adjacency_at_source_between(cap(ef),e1,ef);
    if ( pos == BEFORE ) v->set_halfedge(e2);
  } else
    close_tip_at_source(e1,v);
  return e1;
}

void delete_halfedge_pair_at_source(Halfedge_handle e) const
/*{\Mop deletes |e| and its twin and updates the adjacency at its
  source.}*/
{ remove_from_adj_list_at_source(e);
  this->phds->edges_erase(e);
}

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) const
/*{\Mop makes |v| the target of |e| and appends |twin(e)| to $A(v)$.}*/
{ if ( ! is_isolated(v) )
    set_adjacency_at_source_between(cap(first_out_edge(v)),twin(e),
      first_out_edge(v));
  else
    close_tip_at_target(e,v);
}

Halfedge_handle new_halfedge_pair_without_vertices() const
/*{\Mop inserts an open edge pair, and inits all link slots to their default
    handles.}*/
{
  Halfedge_handle e_new = this->phds->edges_push_back(Halfedge(),Halfedge());
  return e_new;
}

void delete_vertex_only(Vertex_handle v) const
/*{\Mop deletes |v| without consideration of adjacency.}*/
{ this->phds->vertices_erase(v); }

void delete_halfedge_pair_only(Halfedge_handle e) const
/*{\Mop deletes |e| and its twin without consideration of adjacency.}*/
{ this->phds->edges_erase(e); }

void link_as_target_of(Halfedge_handle e, Vertex_handle v) const
/*{\Mop makes |target(e) = v| and sets |e| as the first
        in-edge if |v| was isolated before.}*/
{ e->set_vertex(v);
  if (v->halfedge() == Halfedge_handle()) v->set_halfedge(e); }

void link_as_source_of(Halfedge_handle e, Vertex_handle v) const
/*{\Mop makes |source(e) = v| and sets |e| as the first
        out-edge if |v| was isolated before.}*/
{ link_as_target_of(e->opposite(),v); }

void make_first_out_edge(Halfedge_handle e) const
/*{\Mop makes |e| the first outgoing halfedge in the cyclic adjacency
    list of |source(e)|.}*/
{ source(e)->set_halfedge(e->opposite()); }


void set_adjacency_at_source_between(Halfedge_handle e, Halfedge_handle en)
  const
/*{\Mop makes |e| and |en| neigbors in the cyclic ordered adjacency list
  around |v=source(e)|. \precond |source(e)==source(en)|.}*/
{ CGAL_assertion(source(e)==source(en));
  link_as_prev_next_pair(en->opposite(),e);
}

void set_adjacency_at_source_between(Halfedge_handle e1,
                                     Halfedge_handle e_between,
                                     Halfedge_handle e2) const
/*{\Mop inserts |e_between| into the adjacency list around |source(e1)|
  between |e1| and |e2| and makes |source(e1)| the source of |e_between|.
  \precond |source(e1)==source(e2)|.}*/
{ e_between->opposite()->set_vertex(source(e1));
  set_adjacency_at_source_between(e1,e_between);
  set_adjacency_at_source_between(e_between,e2);
}

void close_tip_at_target(Halfedge_handle e, Vertex_handle v) const
/*{\Mop sets |v| as target of |e| and closes the tip by setting the
  corresponding pointers such that |prev(twin(e)) == e| and
  |next(e) == twin(e)|.}*/
{ link_as_target_of(e,v);
  link_as_prev_next_pair(e,e->opposite()); }

void close_tip_at_source(Halfedge_handle e, Vertex_handle v) const
/*{\Mop sets |v| as source of |e| and closes the tip by setting the
  corresponding pointers such that |prev(e) == twin(e)| and
  |next(twin(e)) == e|.}*/
{ close_tip_at_target(e->opposite(),v); }


void remove_from_adj_list_at_source(Halfedge_handle e) const
/*{\Mop removes a halfedge pair |(e,twin(e)| from the adjacency list
        of |source(e)|. Afterwards |next(prev(e))==next(twin(e))|
        and |first_out_edge(source(e))| is valid if
        |degree(source(v))>1| before the operation.}*/
{
  Vertex_handle vs = source(e);
  if ( is_closed_at_source(e) ) { // last outgoing
    vs->set_halfedge(Halfedge_handle());
  } else {
    if (e == first_out_edge(vs))
      vs->set_halfedge(e->prev());
    set_adjacency_at_source_between(cap(e),cas(e));
  }
}


void unlink_as_hole(Halfedge_handle e) const
/*{\Mop removes |e|'s existence as an face cycle entry point of |face(e)|.
    Does not update the face links of the corresponding face cycle
    halfedges.}*/
{ e->face()->remove_fc(e); }

void unlink_as_isolated_vertex(Vertex_handle v) const
/*{\Mop removes |v|'s existence as an isolated vertex in |face(v)|.
    Does not update |v|'s face link.}*/
{ v->face()->remove_iv(v); }

void link_as_prev_next_pair(Halfedge_handle e1, Halfedge_handle e2) const
/*{\Mop makes |e1| and |e2| adjacent in the face cycle $\ldots-|e1-e2|-\ldots$.
    Afterwards |e1 = previous(e2)| and |e2 = next(e1)|.}*/
{ e1->set_next(e2); e2->set_prev(e1); }

void set_face(Halfedge_handle e, Face_handle f) const
/*{\Mop makes |f| the face of |e|.}*/
{ e->set_face(f); }

void set_face(Vertex_handle v, Face_handle f) const
/*{\Mop makes |f| the face of |v|.}*/
{ v->set_face(f); }

void set_halfedge(Face_handle f, Halfedge_handle e) const
/*{\Mop makes |e| entry edge in the outer face cycle of |f|.}*/
{ f->set_halfedge(e); }

void set_hole(Face_handle f, Halfedge_handle e) const
/*{\Mop makes |e| entry edge in a hole face cycle of |f|.}*/
{ f->store_fc(e); }

void set_isolated_vertex(Face_handle f, Vertex_handle v) const
/*{\Mop makes |v| an isolated vertex of |f|.}*/
{ f->store_iv(v); }


/*{\Mtext \headerline{Cloning}\setopdims{2cm}{1cm}}*/

void clone(const Plane_map& H) const;
/*{\Mop clones |H| into |P|. Afterwards |P| is a copy of |H|.\\
  \precond |H.check_integrity_and_topological_planarity()| and
  |P| is empty.}*/

template <typename LINKDA>
void clone_skeleton(const Plane_map& H, const LINKDA& L) const;
/*{\Mop clones the skeleton of |H| into |P|. Afterwards |P| is a copy
of |H|. The link data accessor allows to transfer information from
the old to the new objects. It needs the function call operators:\\
|void operator()(Vertex_handle vn, Ver\-tex_\-const_\-handle vo) const|\\
|void operator()(Halfedge_handle hn, Half\-edge_\-const_\-handle ho) const|\\
where |vn,hn| are the cloned objects and |vo,ho| are the original
objects.\\
\precond |H.check_integrity_and_topological_planarity()| and
|P| is empty.}*/


void reflecting_inversion()
/*{\Xop inverts the topological links corresponding to a reflecting
inversion. Assume that the plane map is embedded into the x-y plane
and one looks at it from the tip of the positive z-axis in space. Now
change your view point to a point on the negative z-axis. As a
consequence faces are right of edges and adjacency list are clockwise
order-preserving. This operation recreates our embedding invariant
(faces are left of edges and adjacency lists are counterclockwise
order-preserving).}*/
{
  // swap faces:
  Halfedge_iterator e;
  for (e = halfedges_begin(); e != halfedges_end(); ++(++e)) {
    Face_handle f1 = face(e), f2 = face(twin(e));
    e->set_face(f2); twin(e)->set_face(f1);
  }
  // reverse adjacency lists:
  std::vector<Halfedge_handle> A;
  Vertex_iterator v;
  for (v = vertices_begin(); v != vertices_end(); ++v) {
    if ( is_isolted(v) ) continue;
    Halfedge_around_vertex_circulator h = out_edges(v), hend(h);
    CGAL_For_all(h,hend) A.push_back(h);
    int n = A.size();
    for (int i=0; i<n; ++i)
      set_adjacency_at_source_between(A[(i+1)%n],A[i],A[(i-1)%n]);
  }
  CGAL_error_msg("test this");
}

/*{\Mtext \headerline{Associated Information}\restoreopdims}*/

Point& point(Vertex_handle v) const
/*{\Mop returns the embedding of |v|.}*/
{ return v->point(); }


Mark& mark(Vertex_handle v) const
/*{\Mop returns the mark of |v|.}*/
{ return v->mark(); }

Mark& mark(Halfedge_handle e) const
/*{\Mop returns the mark of |e|.}*/
{ if (&*e < &*(e->opposite())) return e->mark();
  else return e->opposite()->mark();
} // we store the mark in the container with smaller memory address !

Mark& mark(Face_handle f) const
/*{\Mop returns the mark of |f|.}*/
{ return f->mark(); }

void set_marks_in_face_cycle(Halfedge_handle e, Mark m) const
{
  Halfedge_around_face_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) {
    mark(hfc) = mark(target(hfc)) = m;
  }
}

GenPtr& info(Vertex_handle v) const
/*{\Mop returns a generic information slot.}*/
{ return v->info(); }

GenPtr& info(Halfedge_handle e) const
/*{\Mop returns a generic information slot.}*/
{ return e->info(); }

GenPtr& info(Face_handle f) const
/*{\Mop returns a generic information slot.}*/
{ return f->info(); }


}; // PM_decorator<HDS>

template <typename  HDS>
void PM_decorator<HDS>::clone(const HDS& H) const
{
  CGAL_assertion(this->number_of_vertices()==0&&
                 this->number_of_halfedges()==0&&
                 this->number_of_faces()==0);

  PM_const_decorator<HDS> DC(H);
  CGAL_assertion((DC.check_integrity_and_topological_planarity(),1));
  CGAL::Unique_hash_map<Vertex_const_iterator,Vertex_handle>     Vnew;
  CGAL::Unique_hash_map<Halfedge_const_iterator,Halfedge_handle> Hnew;
  CGAL::Unique_hash_map<Face_const_iterator,Face_handle>         Fnew;

  /* First clone all objects and store correspondance in three maps.*/
  Vertex_const_iterator vit, vend = H.vertices_end();
  for (vit = H.vertices_begin(); vit!=vend; ++vit)
    Vnew[vit] = this->phds->vertices_push_back(Vertex_base());
  Halfedge_const_iterator eit, eend = H.halfedges_end();
  for (eit = H.halfedges_begin(); eit!=eend; ++(++eit)) {
    Hnew[eit] = this->phds->edges_push_back(Halfedge_base(),Halfedge_base());
    Hnew[eit->opposite()] = Hnew[eit]->opposite();
  }
  Face_const_iterator fit, fend = H.faces_end();
  for (fit = H.faces_begin(); fit!=fend; ++fit) {
    Fnew[fit] = this->phds->faces_push_back(Face_base());
  }

  /* Now copy topology.*/
  Vertex_iterator vit2, vend2 = this->phds->vertices_end();
  for (vit = H.vertices_begin(), vit2 = vertices_begin();
       vit2!=vend2; ++vit, ++vit2) {
    mark(vit2) = DC.mark(vit);
    point(vit2) = DC.point(vit);
    if ( DC.is_isolated(vit) ) vit2->set_face(Fnew[vit->face()]);
    else vit2->set_halfedge(Hnew[vit->halfedge()]);
  }
  Halfedge_iterator eit2, eend2 = this->phds->halfedges_end();
  for (eit = H.halfedges_begin(), eit2 = halfedges_begin();
       eit2!=eend2; ++eit, ++eit2) {
    eit2->set_prev(Hnew[eit->prev()]);
    eit2->set_next(Hnew[eit->next()]);
    eit2->set_vertex(Vnew[eit->vertex()]);
    eit2->set_face(Fnew[eit->face()]);
    mark(eit2) = DC.mark(eit);
  }

  Face_iterator fit2, fend2 = faces_end();
  for (fit = H.faces_begin(), fit2 = faces_begin();
       fit2!=fend2; ++fit, ++fit2) {
    fit2->set_halfedge(Hnew[fit->halfedge()]);
      // outer face cycle
    Hole_const_iterator fcit, fcend = holes_end(fit);
    for (fcit = holes_begin(fit); fcit!=fcend; ++fcit) {
      fit2->store_fc(Hnew[fcit]);
    } // hole face cycles
    Isolated_vertex_const_iterator ivit;
    for (ivit = isolated_vertices_begin(fit);
         ivit != isolated_vertices_end(fit); ++ivit)
      fit2->store_iv(Vnew[ivit]);
      // isolated vertices in the interior
    mark(fit2) = DC.mark(fit);
  }
  CGAL_assertion((this->check_integrity_and_topological_planarity(),1));
}


template <typename HDS>
template <typename LINKDA>
void PM_decorator<HDS>::
clone_skeleton(const HDS& H, const LINKDA& L) const
{
  CGAL_assertion(this->number_of_vertices()==0&&
                 this->number_of_halfedges()==0&&
                 this->number_of_faces()==0);

  PM_const_decorator<HDS> DC(H);
  CGAL_assertion((DC.check_integrity_and_topological_planarity(),1));
  CGAL::Unique_hash_map<Vertex_const_iterator,Vertex_handle>     Vnew;
  CGAL::Unique_hash_map<Halfedge_const_iterator,Halfedge_handle> Hnew;

  /* First clone all objects and store correspondance in the two maps.*/
  Vertex_const_iterator vit, vend = H.vertices_end();
  for (vit = H.vertices_begin(); vit!=vend; ++vit) {
    Vertex_handle v = this->phds->vertices_push_back(Vertex_base());
    Vnew[vit] = v;
  }
  Halfedge_const_iterator eit, eend = H.halfedges_end();
  for (eit = H.halfedges_begin(); eit!=eend; ++(++eit)) {
    Halfedge_handle e = this->phds->edges_push_back(Halfedge_base(),Halfedge_base());
    Hnew[eit] = e; Hnew[eit->opposite()] = e->opposite();
  }

  /* Now copy topology.*/
  Vertex_iterator vit2, vend2 = vertices_end();
  for (vit = H.vertices_begin(), vit2 = vertices_begin();
       vit2!=vend2; ++vit, ++vit2) {
    mark(vit2) = DC.mark(vit);
    point(vit2) = DC.point(vit);
    if ( !DC.is_isolated(vit) )
      vit2->set_halfedge(Hnew[vit->halfedge()]);
    L(vit2,vit);
  }
  Halfedge_iterator eit2, eend2 = this->phds->halfedges_end();
  for (eit = H.halfedges_begin(), eit2 = halfedges_begin();
       eit2!=eend2; ++eit, ++eit2) {
    eit2->set_prev(Hnew[eit->prev()]);
    eit2->set_next(Hnew[eit->next()]);
    eit2->set_vertex(Vnew[eit->vertex()]);
    mark(eit2) = DC.mark(eit);
    // eit2->set_face(Face_handle((Face*)&*(eit->face())));
    L(eit2,eit);
    // link to face of original
  }
}

} //namespace CGAL

#endif //CGAL_PM_DECORATOR_H
