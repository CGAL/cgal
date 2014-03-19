
template < class I, class P > struct Substitute_iterator;

template < class I, class P >
bool operator==(const Substitute_iterator<I,P>&, const Substitute_iterator<I,P>&);

template < class I, class P >
struct Substitute_iterator {
  typedef I                                Iterator;
  typedef P                                Predicate;
  typedef Substitute_iterator<I,P>             Self;
  typedef std::iterator_traits<I>          ITI;
  typedef typename ITI::reference          reference;
  typedef typename ITI::pointer            pointer;
  typedef typename ITI::value_type         value_type;
  typedef typename ITI::difference_type    difference_type;
  typedef typename ITI::iterator_category  iterator_category;
  //  // Special for circulators.
  //typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
  //typedef typename  C_S_Traits::size_type               size_type;

protected:
  Iterator  c_;       // current position.
  reference o_;       // object to substitute
  Predicate p_;       // Substitute o_ if  p_(x).
public:

  Substitute_iterator() {}

  Substitute_iterator(Iterator c, const Predicate& p, reference o)
  : c_(c), p_(p), o_(o)
  {
  }

  Self& operator++() {
    ++c_;
    return *this;
  }

  Self& operator--() {
    --c_;
    return *this;
  }

  Self operator++(int) {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self operator--(int) {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  reference operator*() const { if (p_(*c_)) return o_ ; else  return *c_;  }
  pointer operator->() const  { if (p_(*c_)) return &o_; else  return &*c_; }
  const Predicate& predicate() const { return p_; }
  Iterator base() const { return c_; }

  //bool is_end() const { return (c_ == e_); }

  friend bool operator== <>(const Self&, const Self&);



  // OPERATIONS Random Access Category
  // ---------------------------------

  Self& operator+=( difference_type n) {
    c_ += n;
    return *this;
  }
  Self  operator+( difference_type n) const {
    Self tmp = *this;
    return tmp += n;
  }
  Self& operator-=( difference_type n) {
    return operator+=( -n);
  }
  Self  operator-( difference_type n) const {
    Self tmp = *this;
    return tmp += -n;
  }
  difference_type  operator-( const Self& i) const { return c_ - i.c_; }
  reference  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
  }
  bool operator< ( const Self& i) const { return ( c_ < i.c_); }
  bool operator> ( const Self& i) const { return i < *this; }
  bool operator<=( const Self& i) const { return !(i < *this); }
  bool operator>=( const Self& i) const { return !(*this < i); }

};



template < class I, class P >
inline
bool operator==(const Substitute_iterator<I,P>& it1,
                const Substitute_iterator<I,P>& it2)
{
  return it1.base() == it2.base();
}

template < class I, class P >
inline
bool operator!=(const Substitute_iterator<I,P>& it1,
                const Substitute_iterator<I,P>& it2)
{ return !(it1 == it2); }
