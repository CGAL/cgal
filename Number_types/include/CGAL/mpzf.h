#include <cstdlib>
#include <assert.h>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <gmp.h>

// FIXME:
// this code is experimental. It assumes there is an int64_t type, it
// may assume little endianness, it uses a gcc builtin, things get
// instantiated that should be in src/, others work only thanks to
// inlining, etc. And the aors function is complicated enough that it
// likely has bugs.

// On a dataset provided by Andreas, replacing Gmpq with this type in
// Epick reduced the running time of the construction of a Delaunay
// triangulation by a factor larger than 6 (more is possible, see
// the comment before mpzf::init()).

#if !defined(CGAL_HAS_THREADS)
#define CGAL_THREAD_LOCAL
#elif __cplusplus >= 201103L
#define CGAL_THREAD_LOCAL thread_local
#elif defined(_MSC_VER)
#define CGAL_THREAD_LOCAL __declspec(thread)
#else
#define CGAL_THREAD_LOCAL __thread
// Too bad for the others
#endif

namespace mpzf_impl {
template <class T, class = void> struct pool1 {
  static T pop() { T ret = data.back(); data.pop_back(); return ret; }
  static void push(T t) { data.push_back(t); }
  static bool empty() { return data.empty(); }
  static const int extra = 0;
  private:
  // thread_local would be fine, but not __declspec(thread) and for most
  // compilers not __thread, since data is not POD.
  static std::vector<T> data;
};
template <class T, class D> std::vector<T> pool1<T,D>::data;

// Use an intrusive single-linked list instead (allocate one more limb and use
// it to store the pointer to next), the difference isn't that noticable (still
// the list wins).  Neither is thread-safe (both can be with threadlocal, and
// the list can be with an atomic compare-exchange). With gcc, TLS affects the
// vector version (almost +20% when constructing a Delaunay triangulation), but
// for the list it seems basically free?!

// Warning: this isn't truly generic!
template <class T, class = void> struct pool2 {
  static T pop() { T ret = data; data = ptr(data); return ret; }
  static void push(T t) { ptr(t) = data; data = t; }
  static bool empty() { return data == 0; }
  static const int extra = 1; // TODO: handle the case where a pointer is larger than a mp_limb_t
  private:
  BOOST_STATIC_ASSERT(sizeof(T) >= sizeof(T*));
  static CGAL_THREAD_LOCAL T data;
  static T& ptr(T t) { t -= extra+1; return *reinterpret_cast<T*>(t); }
};
template <class T, class D> CGAL_THREAD_LOCAL T pool2<T,D>::data = 0;
}


struct mpzf {
  typedef mpzf_impl::pool2<mp_limb_t*,mpzf> pool;

  mp_limb_t* data;
  int size;
  int exp;

  struct allocate{};
  struct noalloc{};

  // We could also use a fixed-size allocation, since evaluating a
  // polynomial in double will never require more than roughly
  // (2100*degree) bits, or (33*degree) mp_limb_t, which is very
  // small. I checked by including an array of 150 limbs in every mpzf
  // and it shaved another 17% from the running-time.
  // BONUS: doing that would be thread-safe!
  void init(unsigned mini=2){
    if(!pool::empty()){
      data = pool::pop();
      if(data[-1] >= mini) return; // TODO: when mini==2, no need to check
      data -= pool::extra+1;
      delete[] data; // too small, useless
    }
    if(mini<2) mini=2;
    data = (new mp_limb_t[mini+(pool::extra+1)]) + (pool::extra+1);
    data[-1] = mini;
  }
  void clear(){
    while(*--data==0); ++data; // in case we skipped final zeroes
    pool::push(data);
  }
  ~mpzf(){ clear(); }
  mpzf(): size(0), exp(0) { init(); }
  mpzf(noalloc){}
  mpzf(allocate,int i) { init(i); }
  mpzf& operator=(mpzf const& x){
    unsigned asize=std::abs(x.size);
    if(asize==0) { exp=0; size=0; return *this; }
    if(this==&x) return *this;
    while(*--data==0); // factor that code somewhere?
    if(*data<asize){
      data -= pool::extra;
      delete[] data;
      init(asize);
    } else ++data;
    size=x.size;
    exp=x.exp;
    mpn_copyi(data,x.data,asize);
    return *this;
  }
  mpzf(mpzf const& x){
    int asize=std::abs(x.size);
    init(asize);
    size=x.size;
    exp=x.exp;
    if(size!=0) mpn_copyi(data,x.data,asize);
  }
#if __cplusplus >= 201103L
  mpzf(mpzf&& x):data(x.data),size(x.size),exp(x.exp){
    x.init(); // yes, that's a shame...
  }
  mpzf& operator=(mpzf&& x){
    // Should have 2 levels, one with just the members and one with
    // the algos, so that I could write a single std::swap.
    std::swap(size,x.size);
    std::swap(exp ,x.exp );
    std::swap(data,x.data);
    return *this;
  }
#endif
  mpzf(double d){
    init();
    union {
      struct { uint64_t man:52; uint64_t exp:11; uint64_t sig:1; } s;
      double d;
    } u;
    u.d = d;
    if(u.s.exp==0){
      if(d==0){ size=0; exp=0; return; }
      throw "denormal\n";
    }
    uint64_t m = (1L<<52)|u.s.man;
    int e1 = (int)u.s.exp+13;
    int e2 = e1 % 64;
    exp = e1 / 64 - 17;
    if(__builtin_ctzll(m)+e2>=64){
      data[0] = m >> (64-e2);
      size = 1;
      ++exp;
    }else{
      data[0] = m << e2;
      if(e2>11){
	data[1] = m >> (64-e2);
	size = 2;
      } else {
	size = 1;
      }
    }
    if(u.s.sig) size=-size;
    //print();
  }
  // For debug purposes only
  void print()const{
    //std::cout << "size: " << size << std::endl;
    if(size==0) { std::cout << "zero\n"; return; }
    if(size<0) std::cout << "- ";
    int asize = std::abs(size);
    std::cout << std::hex;
    while(--asize>=0) { std::cout << data[asize] << ' '; }
    std::cout << std::dec << "exp " << exp << ' ';
    asize = std::abs(size);
    std::cout << "double: " << ldexp(data[asize-1],64*(exp+asize-1))*((size<0)?-1:1) << '\n';
  }
  friend int abscmp(mpzf const&a, mpzf const&b){
    int asize=std::abs(a.size);
    int bsize=std::abs(b.size);
    int ah=asize+a.exp;
    int bh=bsize+b.exp;
    if(ah!=bh) return ah-bh;
    int minsize=std::min(asize,bsize);
    mp_limb_t* adata=a.data+(asize-1);
    mp_limb_t* bdata=b.data+(bsize-1);
    for(int i=0;i<minsize;++i,--adata,--bdata){
      mp_limb_t aa=*adata;
      mp_limb_t bb=*bdata;
      if(aa!=bb) return (aa<bb)?-1:1;
    }
    return asize-bsize; // this assumes that we get rid of trailing zeros...
  }
  friend bool operator<(mpzf const&a, mpzf const&b){
    if((a.size^b.size)<0) return a.size<0;
    return ((a.size<0)?abscmp(b,a):abscmp(a,b))<0;
  }
  friend bool operator>(mpzf const&a, mpzf const&b){
    return b<a;
  }
  friend bool operator>=(mpzf const&a, mpzf const&b){
    return !(a<b);
  }
  friend bool operator<=(mpzf const&a, mpzf const&b){
    return !(a>b);
  }
  static mpzf aors(mpzf const&a, mpzf const&b, int bsize){
    mpzf res=noalloc();
    if(bsize==0){
      int size=std::abs(a.size);
      res.init(size);
      res.exp=a.exp;
      res.size=a.size;
      if(size!=0) mpn_copyi(res.data,a.data,size);
      return res;
    }
    int asize=a.size;
    if(asize==0){
      int size=std::abs(bsize);
      res.init(size);
      res.exp=b.exp;
      res.size=bsize;
      mpn_copyi(res.data,b.data,size);
      return res;
    }
    if((asize^bsize)>=0){
      // Addition
      int absasize=std::abs(asize);
      int absbsize=std::abs(bsize);
      mp_limb_t* adata=a.data;
      mp_limb_t* bdata=b.data;
      int aexp=a.exp;
      int bexp=b.exp;
      if(aexp<bexp){ res.exp=a.exp; aexp=0; bexp=b.exp-a.exp; }
      else { res.exp=b.exp; aexp=a.exp-b.exp; bexp=0; }
      res.init(std::max(absasize+aexp,absbsize+bexp)+1);
      mp_limb_t* rdata=res.data;
      res.size=0;
      // TODO: if aexp>0, swap a and b so we don't repeat the code.
      if(0<bexp){
	if(absasize<=bexp){ // no overlap
	  mpn_copyi(rdata, adata, absasize);
	  rdata+=absasize;
	  rdata=std::fill_n(rdata,bexp-absasize,0);
	  mpn_copyi(rdata, bdata, absbsize);
	  res.size=absbsize+bexp;
	  if(bsize<0) res.size=-res.size;
	  return res;
	} else {
	  mpn_copyi(rdata, adata, bexp);
	  adata+=bexp;
	  absasize-=bexp;
	  rdata+=bexp;
	  res.size=bexp;
	}
      }
      else if(0<aexp){
	if(absbsize<=aexp){ // no overlap
	  mpn_copyi(rdata, bdata, absbsize);
	  rdata+=absbsize;
	  rdata=std::fill_n(rdata,aexp-absbsize,0);
	  mpn_copyi(rdata, adata, absasize);
	  res.size=absasize+aexp;
	  if(asize<0) res.size=-res.size;
	  return res;
	} else {
	  mpn_copyi(rdata, bdata, aexp);
	  bdata+=aexp;
	  absbsize-=aexp;
	  rdata+=aexp;
	  res.size=aexp;
	}
      }
      if(absasize>=absbsize){
	mp_limb_t carry=mpn_add(rdata,adata,absasize,bdata,absbsize);
	res.size+=absasize;
	if(carry!=0){
	  res.size++;
	  rdata[absasize]=carry;
	}
      } else {
	mp_limb_t carry=mpn_add(rdata,bdata,absbsize,adata,absasize);
	res.size+=absbsize;
	if(carry!=0){
	  res.size++;
	  rdata[absbsize]=carry;
	}
      }
      // unnecessary if a.exp != b.exp
      while(/*res.size>0&&*/res.data[0]==0){--res.size;++res.data;++res.exp;}
      if(bsize<0) res.size=-res.size;
    } else {
      // Subtraction
      const mpzf *x, *y;
      int xsize=a.size;
      int ysize=bsize;
      int cmp=abscmp(a,b);
      if(cmp==0){ res.init(); res.size=0; res.exp=0; return res; }
      if(cmp<0) { x=&b; y=&a; std::swap(xsize, ysize); }
      else { x=&a; y=&b; }
      int absxsize=std::abs(xsize);
      int absysize=std::abs(ysize);
      mp_limb_t* xdata=x->data;
      mp_limb_t* ydata=y->data;
      int xexp=x->exp;
      int yexp=y->exp;
      if(xexp<yexp){ res.exp=xexp; yexp-=xexp; xexp=0; }
      else { res.exp=yexp; xexp-=yexp; yexp=0; }
      res.init(std::max(absxsize+xexp,absysize+yexp)+1);
      mp_limb_t* rdata=res.data;
      res.size=0;
      bool carry1=false;
      if(0<yexp){ // must have overlap since x is larger
	mpn_copyi(rdata, xdata, yexp);
	xdata+=yexp;
	absxsize-=yexp;
	rdata+=yexp;
	res.size=yexp;
      }
      else if(0<xexp){
	if(absysize<=xexp){ // no overlap
	  mpn_neg(rdata, ydata, absysize); // assert that it returns 1
	  rdata+=absysize;
	  rdata=std::fill_n(rdata,xexp-absysize,-1);
	  mpn_sub_1(rdata, xdata, absxsize, 1);
	  res.size=absxsize+xexp;
	  if(res.data[res.size-1]==0) --res.size;
	  if(xsize<0) res.size=-res.size;
	  return res;
	} else {
	  mpn_neg(rdata, ydata, xexp); // assert that it returns 1
	  ydata+=xexp;
	  absysize-=xexp;
	  rdata+=xexp;
	  res.size=xexp;
	  carry1=true; // assumes no trailing zeros
	}
      }
      int carry=mpn_sub(rdata, xdata, absxsize, ydata, absysize);
      if(carry1) carry+=mpn_sub_1(rdata, rdata, absxsize, 1);
      assert(carry==0);
      res.size+=absxsize;
      while(/*res.size>0&&*/res.data[res.size-1]==0) --res.size;
      if(xsize<0) res.size=-res.size;
    }
    return res;
  }
  friend mpzf operator+(mpzf const&a, mpzf const&b){
    return aors(a,b,b.size);
  }
  friend mpzf operator-(mpzf const&a, mpzf const&b){
    return aors(a,b,-b.size);
  }
  friend mpzf operator*(mpzf const&a, mpzf const&b){
    int asize=std::abs(a.size);
    int bsize=std::abs(b.size);
    int siz=asize+bsize;
    mpzf res(allocate(),siz);
    if(asize==0||bsize==0){res.exp=0;res.size=0;return res;}
    res.exp=a.exp+b.exp;
    // TODO: call mpn_mul_123456 if possible?
    mp_limb_t high;
    if(asize>=bsize)
      high = mpn_mul(res.data,a.data,asize,b.data,bsize);
    else
      high = mpn_mul(res.data,b.data,bsize,a.data,asize);
    if(high==0) --siz;
    if(res.data[0]==0) { ++res.data; ++res.exp; --siz; }
    res.size=((a.size^b.size)>=0)?siz:-siz;
    //res.print();
    return res;
  }
  friend mpzf square(mpzf const&a){
    int asize=std::abs(a.size);
    int siz=2*asize;
    mpzf res(allocate(),siz);
    res.exp=2*a.exp;
    if(asize==0){res.size=0;return res;}
    mpn_sqr(res.data,a.data,asize);
    mp_limb_t high = res.data[siz-1];
    if(high==0) --siz;
    if(res.data[0]==0) { ++res.data; ++res.exp; --siz; }
    res.size=siz;
    //res.print();
    return res;
  }
};

#if 0
void f(double d){
  union { struct { uint64_t man:52; uint64_t exp:11; uint64_t sig:1; } s; double d; } u;
  u.d = d;
  uint64_t m = (1L<<52)|u.s.man;
  int e = (int)u.s.exp-1075;
  std::cout << std::hex << m << ' ';
  std::cout << std::dec << e << ' ';
  std::cout << (u.s.sig?"-":"+") << ldexp(m,e) << std::endl;
}

void g(double d){
  union { struct { uint64_t man:52; uint64_t exp:11; uint64_t sig:1; } s; double d; } u;
  u.d = d;
  if(u.s.exp==0){std::cout << "zero or denorm\n"; return;}
  uint64_t m = (1L<<52)|u.s.man;
  int e1 = (int)u.s.exp+13;
  int e2 = e1 % 64;
  int e = e1 / 64 - 17;
  uint64_t m2;
  uint64_t m1;
  __int128 mm;
  if(__builtin_ctzll(m)+e2>=64){
    m2 = 0;
    m1 = m >> (64-e2);
    mm = m1;
    ++e;
  }else{
    m1 = m << e2;
    m2 = e2 ? (m >> (64-e2)) : 0;
    mm = ((__int128)m2 << 64) | m1;
  }
  std::cout << std::hex << m2 << ' ' << m1 << "   ";
  std::cout << std::dec << e << ' ';
  std::cout << (u.s.sig?"-":"+") << ldexp(mm,64*e) << std::endl;
}

int main(){
//  g(ldexp(1,-999));
//  g(ldexp(1,999));
//  g(0);
//  mpzf x(-.000001);
  mpzf x=+ldexp(1,128);
  mpzf y=x+ldexp(1,65); y.print();
  mpzf z=+ldexp(1,132);
  //y*y;
  //y*y*y*y*y*y*y*y*y*y*y;
  (z+y).print();
  (y+z).print();
  (z-y).print();
  (y-z).print();
}
#endif

#if 1
namespace CGAL {
  template <> class Algebraic_structure_traits< mpzf >
    : public Algebraic_structure_traits_base< mpzf, Integral_domain_without_division_tag >  {
      public:
	typedef Tag_true            Is_exact;
	typedef Tag_false            Is_numerical_sensitive;
    };
  template <> class Real_embeddable_traits< mpzf >
    : public INTERN_RET::Real_embeddable_traits_base< mpzf , CGAL::Tag_true > {
      public:
	class Sgn
	  : public std::unary_function< Type, ::CGAL::Sign > {
	    public:
	      ::CGAL::Sign operator()( const Type& x ) const {
		return CGAL::Sign((x.size<0)?-1:(x.size>0));
	      }
	  };
    };
}
#endif
