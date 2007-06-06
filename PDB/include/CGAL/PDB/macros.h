#ifndef CGAL_PDB_MACROS_H
#define CGAL_PDB_MACROS_H

#define CGAL_PDB_SUBSCRIPT(type, expr) type& operator[](unsigned int i){ expr;}\
  const type& operator[](unsigned int i) const { expr;}

#define CGAL_PDB_COPY_CONSTRUCTOR(TC) TC(const TC &o){copy_from(o);}\
  TC& operator=(const TC &o) {copy_from(o); return *this;}

#define CGAL_PDB_ACCESSOR(type, name, expr) const type &name() const{expr;}

#define CGAL_PDB_SET(type, name, expr) void set_##name(const type &k) {expr;}

#define CGAL_PDB_RWFIELD(type, name, var) \
  const type &name() const {return var;}\
  void set_##name(const type &k) {var=k;}


#define CGAL_PDB_OUTPUT(type)\
  inline std::ostream& operator<<(std::ostream&out, const type &t){	\
    return t.write(out);						\
  }

#define CGAL_PDB_ITERATOR(uc_name, lc_name, it_type, bexpr, eexpr)	\
  typedef it_type uc_name##_iterator;					\
  uc_name##_iterator lc_name##s_begin() {bexpr;}			\
  uc_name##_iterator lc_name##s_end() {eexpr;}

#define CGAL_PDB_CONST_ITERATOR(uc_name, lc_name, it_type, bexpr, eexpr) \
  typedef it_type uc_name##_const_iterator;				\
  uc_name##_const_iterator lc_name##s_begin() const {bexpr;}		\
  uc_name##_const_iterator lc_name##s_end() const {eexpr;}

#define CGAL_PDB_FIND(ucname, expr)				\
  ucname##_const_iterator find(ucname##_key k) const {expr;}	\
  ucname##_iterator find(ucname##_key k) {expr;}

#define CGAL_PDB_INSERT(ucname, expr)			\
  ucname##_iterator insert(ucname##_key k, const ucname &m) {expr;}

#define CGAL_PDB_SIZE(lcname, expr)		\
  size_t number_of_##lcname() const {expr;}

#define CGAL_PDB_SWAP(type)			\
  inline void swap(type &a, type &b) {		\
    a.swap_with(b);				\
  }

#define CGAL_PDB_IFNONEQUAL(a,b,cmp) if (a cmp b) return true;	\
  else if (b cmp a) return false;

#endif
