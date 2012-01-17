// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_ITERATOR_ADAPTORS_H
#define CGAL_VORONOI_DIAGRAM_2_ITERATOR_ADAPTORS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <iterator>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class VDA, class Iterator, class Base_iterator, class Value_type>
class Iterator_adaptor_base
{
 private:
  typedef Iterator_adaptor_base<VDA,Iterator,Base_iterator,Value_type> Self;

 public:
  typedef Value_type                                 value_type;
  typedef Value_type&                                reference;
  typedef Value_type*                                pointer;
  typedef typename Base_iterator::size_type          size_type;
  typedef typename Base_iterator::difference_type    difference_type;
  typedef typename Base_iterator::iterator_category  iterator_category;

  Iterator_adaptor_base(const VDA* vda = NULL) : vda_(vda), value_() {}
  Iterator_adaptor_base(const VDA* vda, const Base_iterator& cur)
    : vda_(vda), value_(), cur_(cur) {}
  Iterator_adaptor_base(const Self& other) { copy_from(other); }

  Self& operator=(const Self& other) {
    copy_from(other);
    return *this;
  }

  pointer operator->() const {
    const Iterator* ptr = static_cast<const Iterator*>(this);
    ptr->eval_pointer();
    return &value_;
  }

  reference operator*() const {
    const Iterator* ptr = static_cast<const Iterator*>(this);
    ptr->eval_reference();
    return value_;
  }

  Iterator& operator++() {
    Iterator* ptr = static_cast<Iterator*>(this);
    ptr->increment();
    return *ptr;
  }

  Iterator& operator--() {
    Iterator* ptr = static_cast<Iterator*>(this);
    ptr->decrement();
    return *ptr;
  }

  Iterator operator++(int) {
    Iterator tmp(*this);
    (static_cast<Iterator*>(this))->operator++();
    return tmp;
  }

  Iterator operator--(int) {
    Iterator tmp(*this);
    (static_cast<Iterator*>(this))->operator--();
    return tmp;
  }

  bool operator==(const Iterator& other) const {
    if ( vda_ == NULL ) { return other.vda_ == NULL; }
    if ( other.vda_ == NULL ) { return vda_ == NULL; }
    return ( vda_ == other.vda_ && cur_ == other.cur_ );
  }

  bool operator!=(const Iterator& other) const {
    return !(*this == other);
  }

 protected:
  void copy_from(const Self& other) {
    vda_ = other.vda_;
    value_ = other.value_;
    if ( vda_ != NULL ) {
      cur_ = other.cur_;
    }
  }

 protected:
  const VDA* vda_;
  mutable value_type value_;
  Base_iterator cur_;
};

//=========================================================================
//=========================================================================


template<class VDA, class Base_it, class Tag = CGAL::Tag_true>
class Edge_iterator_adaptor
  : public Iterator_adaptor_base<VDA,
				 Edge_iterator_adaptor<VDA,Base_it,Tag>,
				 Base_it,
				 typename VDA::Halfedge>
{
 protected:
  typedef Triangulation_cw_ccw_2   CW_CCW_2;

 private:
  typedef Edge_iterator_adaptor<VDA,Base_it,Tag>        Self;
  // Base_it is essentially VDA::Non_degenerate_edges_iterator
  typedef Base_it                                       Base_iterator;
  typedef typename VDA::Halfedge                        Halfedge;
  typedef typename VDA::Halfedge_handle                 Halfedge_handle;
  typedef Halfedge                                      Value;

  typedef Iterator_adaptor_base<VDA,Self,Base_iterator,Value> Base;

  friend class Iterator_adaptor_base<VDA,Self,Base_iterator,Value>;

 public:
  Edge_iterator_adaptor(const VDA* vda = NULL) : Base(vda) {}
  Edge_iterator_adaptor(const VDA* vda, const Base_iterator& cur)
    : Base(vda, cur) {}

  operator Halfedge_handle() const {
    eval_reference();
    return Halfedge_handle(this->value_);
  }

 private:
  Edge_iterator_adaptor(const Base& base) : Base(base) {}

  void eval_pointer() const {
    eval_pointer(Tag());
  }

  void eval_pointer(const Tag_false&) const {
    this->value_ = *this->cur_;
  }

  void eval_pointer(const Tag_true&) const {
    if ( this->vda_->dual().dimension() == 1 ) {
      int cw_i  = CW_CCW_2::cw( this->cur_->second );
      int ccw_i = CW_CCW_2::ccw( this->cur_->second );
      this->value_ =
	Halfedge(this->vda_, this->cur_->first->vertex(ccw_i),
		 this->cur_->first->vertex(cw_i));
      return;
    }

    this->value_ =
      Halfedge(this->vda_, this->cur_->first, this->cur_->second);

    typename VDA::Delaunay_graph::Edge e = this->value_.dual();

    int j = CW_CCW_2::ccw( e.second );
    typename VDA::Delaunay_graph::Vertex_handle v = e.first->vertex(j);

    if ( this->vda_->face_rejector()(this->vda_->dual(), v) ) {
      this->value_ = *this->value_.opposite();
    }
  }

  void eval_reference() const {
    eval_pointer();
  }

  void increment() {
    ++this->cur_;
  }

  void decrement() {
    --this->cur_;
  }
};


//=========================================================================
//=========================================================================

template<class VDA, class Base_it>
class Halfedge_iterator_adaptor
  : public Iterator_adaptor_base<VDA,
				 Halfedge_iterator_adaptor<VDA,Base_it>,
				 Base_it,
				 typename VDA::Halfedge>
{
private:
  typedef Halfedge_iterator_adaptor<VDA,Base_it>   Self;

  typedef typename VDA::Halfedge_handle    Halfedge_handle;
  typedef typename VDA::Halfedge           Halfedge;
  typedef Halfedge                         Value;
  // this is essentially VDA::Edge_iterator
  typedef Base_it                          Base_iterator;

  typedef Iterator_adaptor_base<VDA,Self,Base_iterator,Value> Base;

  friend class Iterator_adaptor_base<VDA,Self,Base_iterator,Value>;

 public:
  Halfedge_iterator_adaptor(const VDA* vda = NULL) : Base(vda) {}
  Halfedge_iterator_adaptor(const VDA* vda, Base_iterator cur)
    : Base(vda, cur), is_first_(true) {}

  Halfedge_iterator_adaptor(const Self& other)
    : Base(other) { copy_from(other); }

  operator Halfedge_handle() const {
    eval_reference();
    return Halfedge_handle(this->value_);
  }

  Self& operator=(const Self& other) {
    copy_from(other);
    return *this;
  }

  bool operator==(const Self& other) const {
    return Base::operator==(other) && is_first_ == other.is_first_;
  }

  bool operator!=(const Self& other) const {
    return !(*this == other);
  }

 private:
  Halfedge_iterator_adaptor(const Base& base)
    : Base(base), is_first_(true) {}

 private:
  void increment() {
    if ( is_first_ ) {
      is_first_ = false;
      return;
    }

    is_first_ = true;
    ++this->cur_;
  }

  void decrement() {
    if ( is_first_ ) {
      is_first_ = false;
      --this->cur_;
      return;
    }

    is_first_ = true;
  }


  void eval_pointer() const {
    this->value_ = *this->cur_;
    if ( !is_first_ ) {
      this->value_ = *this->value_.opposite();
    }
  }

  void eval_reference() const {
    eval_pointer();
  }

  void copy_from(const Self& other) {
    this->vda_ = other.vda_;
    this->value_ = other.value_;
    is_first_ = other.is_first_;
    if ( this->vda_ != NULL ) {
      this->cur_ = other.cur_;
    }
  }

private:
  bool is_first_;
};


//=========================================================================
//=========================================================================

template<class VDA, class Base_it>
class Face_iterator_adaptor
  : public Iterator_adaptor_base<VDA,
				 Face_iterator_adaptor<VDA,Base_it>,
				 Base_it,
				 typename VDA::Face>
{
 private:
  typedef Face_iterator_adaptor<VDA,Base_it>    Self;

  // Base_it is essentially VDA::Delaunay_vertices_iterator
  typedef Base_it                              Base_iterator;
  typedef typename VDA::Face                   Face;
  typedef typename VDA::Face_handle            Face_handle;
  typedef Face                                 Value;

  typedef Iterator_adaptor_base<VDA,Self,Base_iterator,Value> Base;

  friend class Iterator_adaptor_base<VDA,Self,Base_iterator,Value>;

 public:
  Face_iterator_adaptor(const VDA* vda = NULL) : Base(vda) {}

  Face_iterator_adaptor(const VDA* vda,
			const Base_iterator& cur)
    : Base(vda, cur) {}

  operator Face_handle() const {
    eval_reference();
    return Face_handle(this->value_);
  }

 private:
  Face_iterator_adaptor(const Base& base) : Base(base) {}

  void increment() { ++this->cur_; }
  void decrement() { --this->cur_; }

  void eval_pointer() const {
    this->value_ = Face(this->vda_, this->cur_.base());
  }

  void eval_reference() const {
    eval_pointer();
  }
};

//=========================================================================
//=========================================================================

template<class VDA, class Base_it>
class Vertex_iterator_adaptor :
  public Iterator_adaptor_base<VDA,
			       Vertex_iterator_adaptor<VDA,Base_it>,
			       Base_it,
			       typename VDA::Vertex>
{
 private:
  typedef Vertex_iterator_adaptor<VDA,Base_it>    Self;
  // Base_it is essentially VDA::Non_degenerate_vertices_iterator
  typedef Base_it                                 Base_iterator;
  typedef typename VDA::Vertex                    Vertex;
  typedef typename VDA::Vertex_handle             Vertex_handle;
  typedef Vertex                                  Value;

  typedef Iterator_adaptor_base<VDA,Self,Base_iterator,Value> Base;

  friend class Iterator_adaptor_base<VDA,Self,Base_iterator,Value>;

 public:
  Vertex_iterator_adaptor(const VDA* vda = NULL) : Base(vda) {}

  Vertex_iterator_adaptor(const VDA* vda, const Base_iterator& cur)
    : Base(vda, cur) {}

  operator Vertex_handle() const {
    eval_reference();
    return Vertex_handle(this->value_);
  }

 private:
  Vertex_iterator_adaptor(const Base& base) : Base(base) {}

  void increment() { ++this->cur_; }
  void decrement() { --this->cur_; }

  void eval_pointer() const {
    this->value_ = Vertex(this->vda_, this->cur_.base());
  }

  void eval_reference() const {
    eval_pointer();
  }
};

//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_ITERATOR_ADAPTORS_H
