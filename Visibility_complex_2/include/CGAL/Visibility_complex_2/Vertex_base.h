// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_VERTEX_BASE_H
#define CGAL_VISIBILITY_COMPLEX_2_VERTEX_BASE_H

#include <CGAL/basic.h>
#include <CGAL/Visibility_complex_2/ccw_cw_traits.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {


template<class Vc_>
class Vertex_base;

// ------------- Vertex_base class --------------------------

template<class Vc_>
class Vertex_base 
    : public Vc_::Gt::Bitangent_2
{
public:

    typedef typename Vc_::Gt                   Gt;

    typedef typename Gt::Point_2               Point_2;
    typedef typename Gt::Bitangent_2           Bitangent_2;
    typedef Bitangent_2                        BT;
    typedef typename Gt::Arc_2                 Arc_2;
    typedef typename Gt::Disk                  Disk;

public:

    typedef typename BT::Disk_handle           Disk_handle;
    typedef typename BT::Type                  Type;
    typedef typename BT::Type_util             Type_util;

    typedef typename Vc_::Vertex               Vertex;
    typedef typename Vc_::Vertex_handle        Vertex_handle;
    typedef typename Vc_::Edge                 Edge;
    typedef typename Vc_::Edge_handle          Edge_handle;
    typedef typename Vc_::Face                 Face;
    typedef typename Vc_::Face_handle          Face_handle;
    typedef typename Vc_::Vertex_const_handle  Vertex_const_handle;
    typedef typename Vc_::Edge_const_handle    Edge_const_handle;
    typedef typename Vc_::Face_const_handle    Face_const_handle;

private:

    Edge_handle              edge_[6];
    Vertex_handle            pi_;
    Face_handle              sup_;
    Face_handle              inf_;

public:    
  using BT::target_object;
  using BT::source_object;
  using BT::is_xx_left;
  using BT::is_left_xx;
  using BT::type;
  using BT::LL;
  using BT::LR;
  using BT::RL;
  using BT::RR;
  using BT::reverse;

    // CONSTRUCTORS
    Vertex_base() 
	: BT() , 
	  pi_(0) , sup_(0) , inf_(0)
    { 
	for (int i=0;i<6;i++) edge_[i] = 0; 
    }

    Vertex_base(Type t , Disk_handle start , 
					    Disk_handle finish)
	: BT(t,start,finish) ,
	  pi_(0) , sup_(0) , inf_(0)
    {
	for (int i=0;i<6;i++) edge_[i] = 0; 
    }
    Vertex_base(Edge_handle start , Edge_handle finish) 
    {
	*this = Bitangent_2(Type_util()(start->sign(),finish->sign()),
			    start->object(),finish->object(),
			    *start,*finish);
	for (int i=0;i<6;i++) edge_[i] = 0; 
	pi_ = 0; sup_ = 0; inf_ = 0; 
    }

    Vertex_base(const Bitangent_2& b) 
	: Bitangent_2(b)
    { 
	for (int i=0;i<6;i++) edge_[i] = 0; 
	pi_ = 0; sup_ = 0; inf_ = 0;
    }
    Vertex_base(const Vertex_base&sibling,bool reverse,Type t)
      : Bitangent_2(sibling,reverse,t),
        pi_(0) , sup_(0) , inf_(0) {
      for (int i=0;i<6;i++) edge_[i] = 0; 
  }

  // -------------------------------------------------------------------------
    // DESTRUCTOR
    ~Vertex_base() { 
	if (pi_ != 0 && pi_->pi_ == (static_cast<Vertex_handle>(this))) {
	  pi_->pi_ = 0;
        }
	if (sup_ != 0 && sup_->inf() == 
            (static_cast<Vertex_handle>(this))) sup_->set_inf(0);
	if (inf_ != 0 && inf_->sup() == 
            (static_cast<Vertex_handle>(this))) inf_->set_sup(0);
        Edge_handle edges[8];
        for (int i=0;i<6;++i) edges[i]=edge_[i];
        if (pi_) {
          edges[6]=pi_->edge_[4];
          edges[7]=pi_->edge_[5];
        } else {
          edges[6]=0;
          edges[7]=0;
        }
	for (int i = 0 ; i < 8 ; i++) {
	    if (edges[i] != 0 && edges[i]->inf() == 
                (static_cast<Vertex_handle>(this))) edges[i]->set_inf(0);
	    if (edges[i] != 0 && edges[i]->sup() == 
                (static_cast<Vertex_handle>(this))) edges[i]->set_sup(0);
	}
    } 
    // ---------- obstacle bitangents
    bool is_constraint() const { return (edge_[4] != 0 && edge_[5] != 0); }
    // -------- return the opposite oriented bitangent
    Vertex_handle pi();
    Vertex_const_handle pi() const;
private:
    void set_pi_aux(Vertex_handle v);
public:
    void set_pi(Vertex_handle v);

    Edge_handle ccw_target_edge() { return edge_[0]; }
    Edge_handle ccw_source_edge() { return edge_[1]; }
    Edge_handle ccw_edge(Disk_handle p) {
      CGAL_precondition(p == source_object() || p == target_object());
      return (p == source_object()) ? ccw_source_edge() : ccw_target_edge();
    }

    Edge_handle cw_target_edge()  { return edge_[2]; }
    Edge_handle cw_source_edge()  { return edge_[3]; }
    Edge_handle cw_edge(Disk_handle p) {
      CGAL_precondition(p == source_object() || p == target_object());
      return (p == source_object()) ? cw_source_edge() : cw_target_edge();
    }

    Edge_const_handle ccw_target_edge() const { return edge_[0]; }
    Edge_const_handle ccw_source_edge() const { return edge_[1]; }
    Edge_const_handle ccw_edge(Disk_handle p) const {
      CGAL_precondition(p == source_object() || p == target_object());
      return (p == source_object()) ? ccw_source_edge() : ccw_target_edge();
    }
    Edge_const_handle cw_target_edge()  const { return edge_[2]; }
    Edge_const_handle cw_source_edge()  const { return edge_[3]; }
    Edge_const_handle cw_edge(Disk_handle p)  const {
      CGAL_precondition(p == source_object() || p == target_object());
      return (p == source_object()) ? cw_source_edge() : cw_target_edge();
    }



    void set_ccw_edge(const Edge_handle& e);
    void set_cw_edge (const Edge_handle& e);
    void set_ccw_edge(const Edge_handle& e , const Disk_handle& p);
    void set_cw_edge (const Edge_handle& e , const Disk_handle& p);

    Vertex_handle ccR() 
    { return (ccw_target_edge() == 0) ? 0 : ccw_target_edge()->sup(); }
    Vertex_handle ccL() 
    { return (ccw_source_edge() == 0) ? 0 : ccw_source_edge()->sup(); }
    Vertex_handle cwR() 
    { return (cw_source_edge()  == 0) ? 0 :  cw_source_edge()->inf();  }
    Vertex_handle cwL() 
    { return (cw_target_edge()  == 0) ? 0 : cw_target_edge()->inf();  }

    Vertex_const_handle ccR() const 
    { return (ccw_target_edge() == 0) ? 0 : ccw_target_edge()->sup(); }
    Vertex_const_handle ccL() const 
    { return (ccw_source_edge() == 0) ? 0 : ccw_source_edge()->sup(); }
    Vertex_const_handle cwR() const 
    { return (cw_source_edge()  == 0) ? 0 :  cw_source_edge()->inf();  }
    Vertex_const_handle cwL() const 
    { return (cw_target_edge()  == 0) ? 0 : cw_target_edge()->inf();  }


    Edge_handle target_cusp_edge()  { return edge_[4]; }
    Edge_handle source_cusp_edge()  { return edge_[5]; }
    Edge_handle cusp_edge(Disk_handle p) {
	return (p == target_object()) ? edge_[4] : edge_[5];
    }

    Edge_const_handle target_cusp_edge()  const { return edge_[4]; }
    Edge_const_handle source_cusp_edge()  const { return edge_[5]; }
    Edge_const_handle cusp_edge(Disk_handle p) const {
	return (p == target_object()) ? edge_[4] : edge_[5];
    }

    void set_target_cusp_edge(Edge_handle e) { edge_[4] = e; }
    void set_source_cusp_edge(Edge_handle e) { edge_[5] = e; }
    void set_cusp_edge(Edge_handle e, Disk_handle p) {
	if (p == target_object()) edge_[4] = e;
	else edge_[5] = e;
    }

    Face_handle inf()             { return inf_; }
    Face_const_handle inf() const { return inf_; }
    void        set_inf(Face_handle f) { inf_ = f; }
    Face_handle sup()             { return sup_; }
    Face_const_handle sup() const { return sup_; }
    void        set_sup(Face_handle f) { sup_ = f; }

    // The two degenerate faces with sink this
    Face_handle target_cusp_face() 
    {
	if (target_cusp_edge() == 0) return 0;
	return (is_xx_left()) ? target_cusp_edge()->dl(): 
				target_cusp_edge()->dr(); 
    }
    Face_handle source_cusp_face()
    { 
	if (source_cusp_edge() == 0) return 0;
	return (is_left_xx()) ? source_cusp_edge()->ul() : 
				source_cusp_edge()->ur(); 
    }
    Face_handle cusp_face(Disk_handle p)
    {
	CGAL_precondition(p == source_object() || p == target_object());
	return 
          (p == source_object()) ? source_cusp_face() : target_cusp_face();
    }
    Face_const_handle target_cusp_face() const 
    {
	if (target_cusp_edge() == 0) return 0;
	return (is_xx_left()) ? target_cusp_edge()->dl(): 
				target_cusp_edge()->dr(); 
    }
    Face_const_handle source_cusp_face() const
    { 
	if (source_cusp_edge() == 0) return 0;
	return (is_left_xx()) ? source_cusp_edge()->ul() : 
				source_cusp_edge()->ur(); 
    }
    Face_const_handle cusp_face(Disk_handle p) const
    {
	CGAL_precondition(p == source_object() || p == target_object());
	return 
          (p == source_object()) ? source_cusp_face() : target_cusp_face();
    }

    Vertex_handle raw_pi() {
      return pi_;
    }
    Vertex_const_handle raw_pi() const {
      return pi_;
    }
};

template< class Vc_ >
inline void
Vertex_base<Vc_>::set_ccw_edge(const Edge_handle& e , 
						  const Disk_handle& p)
{
    if (target_object() == p) edge_[0] = e;
    else                      edge_[1] = e;
}

template< class Vc_ >
inline void
Vertex_base<Vc_>::set_cw_edge(const Edge_handle& e , 
						 const Disk_handle& p)
{
    if (target_object() == p) edge_[2] = e;
    else                      edge_[3] = e;
}

template< class Vc_ >
inline void
Vertex_base<Vc_>::set_ccw_edge(const Edge_handle& e)
{
    CGAL_precondition(e == 0 || e->object() == source_object() || 
				e->object() == target_object());
    if (e == 0) { edge_[0] = 0; edge_[1] = 0; }
    else if (e->object() == target_object()) edge_[0] = e;
    else edge_[1] = e;
}

template< class Vc_ >
inline void
Vertex_base<Vc_>::set_cw_edge(const Edge_handle& e)
{
    CGAL_precondition(e == 0 || e->object() == source_object() || 
			        e->object() == target_object());
    if (e == 0) { edge_[2] = 0; edge_[3] = 0; }
    else if (e->object() == target_object()) edge_[2] = e;
    else edge_[3] = e;
}



template< class Vc_ >
typename Vertex_base<Vc_>::Vertex_handle
Vertex_base<Vc_>::pi()
{
    if (pi_ == 0) {
        pi_ = new Vertex(*this, true, type()); 
	pi_->pi_ = static_cast<Vertex_handle>(this);
    }
    return pi_;
}

template< class Vc_ >
typename Vertex_base<Vc_>::Vertex_const_handle
Vertex_base<Vc_>::pi() const {
  return const_cast<Vertex_base*>(this)->pi();
}

template< class Vc_ >
void
Vertex_base<Vc_>::set_pi_aux(Vertex_handle v)
{
    CGAL_precondition(v != 0);
    if (pi_ != 0 && pi_ != v) {
	if (pi_->cw_target_edge() != 0)
          v->set_cw_edge(pi_->cw_target_edge());
	if (pi_->cw_source_edge() != 0)
          v->set_cw_edge(pi_->cw_source_edge());
	if (pi_->ccw_target_edge() != 0)
          v->set_ccw_edge(pi_->ccw_target_edge());
	if (pi_->ccw_source_edge() != 0)
          v->set_ccw_edge(pi_->ccw_source_edge());
	if (pi_->cwL() != 0) pi_->cw_target_edge()->set_sup(v);
	if (pi_->cwR() != 0) pi_->cw_source_edge()->set_sup(v);
	if (pi_->ccR() != 0) pi_->ccw_target_edge()->set_inf(v);
	if (pi_->ccL() != 0) pi_->ccw_source_edge()->set_inf(v);
	if (pi_->inf() != 0) {
	    v->set_inf(pi_->inf());
	    pi_->inf()->set_sup(v);
	}
	if (pi_->sup() != 0) {
	    v->set_sup(pi_->sup());
	    pi_->sup()->set_inf(v);
	}
	delete pi_; 
    }
    pi_=v;
}

template< class Vc_ >
void
Vertex_base<Vc_>::set_pi(Vertex_handle v) {
    CGAL_precondition(v != 0);
    CGAL_precondition(v->target_object()==source_object());
    CGAL_precondition(v->source_object()==target_object());
    CGAL_precondition(v->type()==reverse(type()));
    if (pi_!=v) set_pi_aux(v);
    if (v->pi_!=static_cast<Vertex_handle>(this)) 
      v->set_pi_aux(static_cast<Vertex_handle>(this));
}

}
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_VERTEX_BASE_H
