#ifndef CGAL_DEBUG_INTEGER_H
#define CGAL_DEBUG_INTEGER_H

#include <CGAL/leda_integer.h>

//CGAL_BEGIN_NAMESPACE

class leda_debug_integer{
  leda_integer li;
  double d;
  inline void synchronize() { d=li.todouble();}
public:
  leda_debug_integer() {li=leda_integer();synchronize();}
  leda_debug_integer(const int i) {li=leda_integer(i);synchronize();}  
  leda_debug_integer(const unsigned int i) {li=leda_integer(i);synchronize();}  
  leda_debug_integer(const long i) {li=leda_integer(i);synchronize();}  
  leda_debug_integer(const unsigned long l) {li=leda_integer(l);synchronize();}  
  leda_debug_integer(const double d) {li=leda_integer(d);synchronize();}  
  /*
  leda_debug_integer(const int i) {li=leda_integer(i);synchronize();}  
  leda_debug_integer(const int i) {li=leda_integer(i);synchronize();}  
  leda_debug_integer(const int i) {li=leda_integer(i);synchronize();}  
  */
  leda_debug_integer(const leda_integer& i) {li=i;synchronize();}
  /*
     integer(int sz, const digit* vec);
     integer(char* s);
     integer(const number_string& s) { PTR=NULL; (*this)=integer(s.cstring()); }
*/

  leda_debug_integer& operator=(const leda_integer& x) { li=x; synchronize(); return *this; }
  leda_debug_integer& operator=(const leda_debug_integer& x) { li=x.li; synchronize(); return *this; }
 
  leda_debug_integer operator-() const { return -li; }
  leda_debug_integer operator~() const { return ~li; }
  leda_debug_integer operator<<(long n) const { return li << n; }
  leda_debug_integer operator>>(long n) const { return li >> n; }
  leda_debug_integer operator+= (const leda_debug_integer& b) { li = li + b.li; synchronize(); return *this; }
  leda_debug_integer operator-= (const leda_debug_integer& b) {  li = li - b.li; synchronize(); return *this; }
  leda_debug_integer operator*= (const leda_debug_integer& b) {  li = li * b.li; synchronize(); return *this; }
  leda_debug_integer operator/= (const leda_debug_integer& b) { li = li / b.li; synchronize(); return *this; }
  leda_debug_integer operator%= (const leda_debug_integer& b) {  li = li % b.li; synchronize(); return *this; }
  leda_debug_integer operator&= (const leda_debug_integer& b) {  li = li & b.li; synchronize(); return *this; }
  leda_debug_integer operator|= (const leda_debug_integer& b) { li = li | b.li; synchronize(); return *this; }
  leda_debug_integer operator<<=(int n) {  li = li << n; synchronize(); return *this; }
  leda_debug_integer operator>>=(int n) {  li = li >> n; synchronize(); return *this; }
  leda_debug_integer operator++ () { li++; synchronize(); return *this; }
  leda_debug_integer operator++ (int) {li.operator++(); synchronize(); return *this; }
  leda_debug_integer operator-- () {li--; synchronize(); return *this; }
  leda_debug_integer operator-- (int) {li.operator--(); synchronize(); return *this; }

  friend leda_debug_integer operator + (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li+b.li; }
  friend leda_debug_integer operator - (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li-b.li; }
  friend leda_debug_integer operator * (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li*b.li; }
  friend leda_debug_integer operator / (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li/b.li; }
  friend leda_debug_integer operator % (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li%b.li; }
  friend leda_debug_integer operator & (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li&b.li; }
  friend leda_debug_integer operator | (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li|b.li; }

  /*  friend leda_debug_integer operator + (const leda_debug_integer& a, const double b) { return a.li+b; }
  friend leda_debug_integer operator - (const leda_debug_integer& a, const double b) { return a.li-b; }
  friend leda_debug_integer operator * (const leda_debug_integer& a, const double b) { return a.li*b; }
  friend leda_debug_integer operator / (const leda_debug_integer& a, const double b) { return a.li/b; }
  friend leda_debug_integer operator % (const leda_debug_integer& a, const double b) { return a.li%b; }
  friend leda_debug_integer operator & (const leda_debug_integer& a, const double b) { return a.li&b; }
  friend leda_debug_integer operator | (const leda_debug_integer& a, const double b) { return a.li|b; }
  */
  
  friend bool operator == (const leda_debug_integer& a, const leda_debug_integer& b)  { return a.li==b.li; }
  friend bool operator <  (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li<b.li; }
  friend bool operator >  (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li>b.li; }
  friend inline bool operator != (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li!=b.li; }
  friend inline bool operator >= (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li>=b.li; }
  friend inline bool operator <= (const leda_debug_integer& a, const leda_debug_integer& b) { return a.li<=b.li; }

  bool operator==(int n) const { return li==n; }
  bool operator< (int n) const { return li<n; }
  bool operator> (int n) const { return li>n; }
  bool operator!=(int n) const { return li!=n; }
  bool operator>=(int n) const { return li>=n; }
  bool operator<=(int n) const { return li<=n; }

  friend inline std::ostream& operator << (std::ostream & out, 
                                           const leda_debug_integer& a){
    return out << a.li;
  }
  friend inline std::istream& operator >> (std::istream & in, 
                                           leda_debug_integer& a){
    in >> a.li;
    a.synchronize();
    return in;
  }

  bool   is_long()   const { return li.is_long();}
  long to_long() const { return li.to_long();}
  double to_double() const { return li.to_double();}
  /*
  number_string to_string() const { return li.to_string();}
  leda_debug_integer& from_string(number_string s)    const { return li.from_string();}
  number_string tostring() const { return li.tostring();}
  leda_debug_integer& fromstring(number_string s)   const { return li.fromstring();}
  */
  double todouble() const { return li.todouble();}
  long   tolong()   const { return li.tolong();}
  bool   islong()   const { return li.islong();}
  /*
#ifndef CGAL_CFG_NO_NAMESPACE
  friend inline
  double
  to_double(const leda_debug_integer & i)
  { return to_double(i.li); }
#endif // CGAL_CFG_NO_NAMESPACE
  */
  friend inline
  Number_tag
  number_type_tag(const leda_debug_integer& )
  { return Number_tag(); }
  
  friend inline
  bool
  is_finite(const leda_debug_integer &)
  { return true; }
  
  friend inline
  bool
  is_valid(const leda_debug_integer &)
  { return true; }
  
  friend inline
  io_Operator
  io_tag(const leda_debug_integer &)
  { return io_Operator(); }
  
#ifndef CGAL_CFG_NO_NAMESPACE
  friend inline
  Sign
  sign(const leda_debug_integer& n)
  { return (Sign)::sign(n.li); }
#endif // CGAL_CFG_NO_NAMESPACE
  
};
  
//CGAL_END_NAMESPACE


#endif
