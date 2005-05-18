#ifndef CGAL_VORONOI_DIAGRAM_2_CIRCULATOR_ADAPTORS_H
#define CGAL_VORONOI_DIAGRAM_2_CIRCULATOR_ADAPTORS_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/circulator_bases.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class Halfedge_handle, class Circulator>
class Circulator_from_halfedge_adaptor
  : public Bidirectional_circulator_base<Halfedge_handle>
{
 private:
  typedef Bidirectional_circulator_base<Halfedge_handle>     Base;
  typedef Circulator_from_halfedge_adaptor<Halfedge_handle,Circulator>  Self;

 public:
  Circulator_from_halfedge_adaptor() : cur_() {}
  Circulator_from_halfedge_adaptor(const Halfedge_handle& he) : cur_(he) {}

  Circulator& operator++() {
    Circulator* ptr = static_cast<Circulator*>(this);
    ptr->increment();
    return *ptr;
  }

  Circulator& operator--() {
    Circulator* ptr = static_cast<Circulator*>(this);
    ptr->decrement();
    return *ptr;
  }

  Circulator operator++(int) {
    Circulator tmp(*this);
    (static_cast<Circulator*>(this))->operator++();
    return tmp;
  }

  Circulator operator--(int) {
    Circulator tmp(*this);
    (static_cast<Circulator*>(this))->operator--();
    return tmp;
  }

  typename Base::pointer   operator->() const { return &cur_; }
  typename Base::reference operator*() { return cur_; }

  bool operator==(const Circulator& other) const {
    return cur_ == other.cur_;

  }

  bool operator!=(const Circulator& other) const {
    return cur_ != other.cur_;
  }

 protected:
  Halfedge_handle cur_;
};

//=========================================================================
//=========================================================================

template<class Halfedge_handle>
class Halfedge_around_vertex_circulator_adaptor
  : public Circulator_from_halfedge_adaptor<Halfedge_handle,
      Halfedge_around_vertex_circulator_adaptor<Halfedge_handle> >
{
 private:
  typedef Halfedge_around_vertex_circulator_adaptor<Halfedge_handle> Self;
  typedef Circulator_from_halfedge_adaptor<Halfedge_handle,Self>     Base;

  friend class Circulator_from_halfedge_adaptor<Halfedge_handle,Self>;

 public:
  Halfedge_around_vertex_circulator_adaptor() : Base() {}
  Halfedge_around_vertex_circulator_adaptor(const Halfedge_handle& he)
    : Base(he) {}

 private:
  Halfedge_around_vertex_circulator_adaptor(const Base& b) : Base(b) {}

  void increment() {
    this->cur_ = this->cur_->next()->opposite();
  }

  void decrement() {
    //    this->cur_ = this->cur_->previous()->opposite();
    this->cur_ = this->cur_->opposite()->previous();
  }
};

//=========================================================================
//=========================================================================

template<class Halfedge_handle>
class Ccb_halfedge_circulator_adaptor
  : public Circulator_from_halfedge_adaptor<Halfedge_handle,
      Ccb_halfedge_circulator_adaptor<Halfedge_handle> >
{
 private:
  typedef Ccb_halfedge_circulator_adaptor<Halfedge_handle>        Self;
  typedef Circulator_from_halfedge_adaptor<Halfedge_handle,Self>  Base;

  friend class Circulator_from_halfedge_adaptor<Halfedge_handle,Self>;

 public:
  Ccb_halfedge_circulator_adaptor() : Base() {}
  Ccb_halfedge_circulator_adaptor(const Halfedge_handle& he) : Base(he) {}

 private:
  Ccb_halfedge_circulator_adaptor(const Base& b) : Base(b) {}

  void increment() {
    this->cur_ = this->cur_->next();
  }

  void decrement() {
    this->cur_ = this->cur_->previous();
  }
};

//=========================================================================
//=========================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_CIRCULATOR_ADAPTORS_H
