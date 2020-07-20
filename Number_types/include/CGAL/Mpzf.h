// Copyright (c) 2013
// INRIA Saclay - Ile de France (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)        :  Marc Glisse

#ifndef CGAL_MPZF_H
#define CGAL_MPZF_H
#include <cstdlib>
#include <algorithm>
#include <climits>
#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#ifdef CGAL_USE_GMPXX
# include <CGAL/gmpxx.h>
#else
# include <CGAL/gmp.h>
#endif
#include <CGAL/enum.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
//#include <CGAL/Gmpzf.h>

#include <CGAL/Coercion_traits.h>

// The following is currently assumed in several places. I hope I am not
// making too many other assumptions.
// * limbs are 64 bits
// * if using gcc, sizeof(long long)==8
// * mpn_neg(_n) exists
// * IEEE double
// * not too fancy endianness
#if !defined(CGAL_DO_NOT_USE_MPZF) \
    && __GNU_MP_VERSION * 10 + __GNU_MP_VERSION_MINOR >= 43 \
    && GMP_NUMB_BITS == 64
#define CGAL_HAS_MPZF 1

// GMP-4.3.* has a different name for mpn_neg.
#ifndef mpn_neg
#define mpn_neg mpn_neg_n
#endif
// GMP-4.3.0 is missing mpn_sqr.
#ifndef mpn_sqr
#define mpn_sqr(dest,a,n) mpn_mul_n(dest,a,a,n)
#endif
// GMP before 5.0 doesn't provide mpn_copyi.
#ifndef mpn_copyi
#define mpn_copyi(dst, src, siz) std::copy((src), (src)+(siz), (dst))
#endif

#include <boost/cstdint.hpp>

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(_BitScanForward64)
#pragma intrinsic(_BitScanReverse64)
#endif

#ifdef __xlC__
#include <builtins.h>
#endif

#include <CGAL/assertions.h>
#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/version.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4146 4244 4267 4702 4800)
     // warning on - applied on unsigned number
     // conversion with loss of data
     // conversion with loss of data
     // unreachable code
     // int to bool performance
#endif


/*
#ifdef CGAL_MPZF_NO_USE_CACHE
# ifdef CGAL_MPZF_USE_CACHE
#  undef CGAL_MPZF_USE_CACHE
# endif
#else
# if !defined(CGAL_MPZF_USE_CACHE) \
     && defined(CGAL_HAS_THREADS) \
     && !defined(CGAL_I_PROMISE_I_WONT_USE_MANY_THREADS)
#  define CGAL_MPZF_USE_CACHE 1
# endif
#endif
*/
#define CGAL_MPZF_USE_CACHE 1

// On a dataset provided by Andreas, replacing Gmpq with this type in
// Epick reduced the running time of the construction of a Delaunay
// triangulation by a factor larger than 6

#if !defined(CGAL_HAS_THREADS)
#define CGAL_MPZF_THREAD_LOCAL
#define CGAL_MPZF_TLS
#elif defined(CGAL_CAN_USE_CXX11_THREAD_LOCAL)
#define CGAL_MPZF_THREAD_LOCAL thread_local
#define CGAL_MPZF_TLS thread_local
#elif defined(_MSC_VER)
#define CGAL_MPZF_THREAD_LOCAL __declspec(thread)
#define CGAL_MPZF_TLS
#else
#define CGAL_MPZF_THREAD_LOCAL __thread
#define CGAL_MPZF_TLS
// Too bad for the others
#endif
namespace CGAL {
namespace Mpzf_impl {
// Warning: these pools aren't generic at all!

// Not thread-safe
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
// it to store the pointer to next), the difference isn't that noticeable (still
// the list wins).  Neither is thread-safe (both can be with threadlocal, and
// the list can be with an atomic compare-exchange (never tried)).  With gcc,
// TLS has a large effect on classes with constructor/destructor, but is free
// for a simple pointer.  The slowdown is because of PR55812.

// Leaks at thread destruction
template <class T, class = void> struct pool2 {
  static T pop() { T ret = data(); data() = ptr(data()); return ret; }
  static void push(T t) { ptr(t) = data(); data() = t; }
  static bool empty() { return data() == 0; }
  static const int extra = 1; // TODO: handle the case where a pointer is larger than a mp_limb_t
  private:
  CGAL_static_assertion(sizeof(T) >= sizeof(T*));
  static T& data () {
    static CGAL_MPZF_TLS T data_ = 0;
    return data_;
  }
  static T& ptr(T t) { t -= extra+1; return *reinterpret_cast<T*>(t); }
};

#if defined(CGAL_CAN_USE_CXX11_THREAD_LOCAL)
template <class T, class = void> struct pool3 {
  static T pop() { T ret = data(); data() = ptr(data()); return ret; }
  static void push(T t) { ptr(t) = data(); data() = t; }
  static bool empty() { return data() == 0; }
  static const int extra = 1; // TODO: handle the case where a pointer is larger than a mp_limb_t
  private:
  CGAL_static_assertion(sizeof(T) >= sizeof(T*));
  struct cleaner {
    T data_ = 0;
    ~cleaner(){
      // Deallocate everything. As an alternative, we could store it in a
      // global location, for re-use by a later thread.
      while (!empty())
        delete[] (pop() - (extra + 1));
    }
  };
  static T& data () {
    static thread_local cleaner obj;
    return obj.data_;
  }
  static T& ptr(T t) { t -= extra+1; return *reinterpret_cast<T*>(t); }
};
#endif

// No caching
template <class T, class = void> struct no_pool {
  static T pop() { throw "Shouldn't be here!"; }
  static void push(T t) { delete [] (t - (extra+1)); }
  static bool empty() { return true; }
  static const int extra = 0;
};

// Only used with an argument known not to be 0.
inline int ctz (boost::uint64_t x) {
#if defined(_MSC_VER)
  unsigned long ret;
  _BitScanForward64(&ret, x);
  return (int)ret;
#elif defined(__xlC__)
  return __cnttz8 (x);
#else
  // Assume long long is 64 bits
  return __builtin_ctzll (x);
#endif
}
inline int clz (boost::uint64_t x) {
#if defined(_MSC_VER)
  unsigned long ret;
  _BitScanReverse64(&ret, x);
  return 63 - (int)ret;
#elif defined(__xlC__)
  // Macro supposedly not defined on z/OS.
  return __cntlz8 (x);
#else
  return __builtin_clzll (x);
#endif
}

// In C++11, std::fill_n returns a pointer to the end, but in C++03,
// it returns void.
inline mp_limb_t* fill_n_ptr(mp_limb_t* p, int n, int c) {
  return std::fill_n (p, n, c);
}
} // namespace Mpzf_impl

#undef CGAL_MPZF_THREAD_LOCAL
#undef CGAL_MPZF_TLS

// TODO:
// * make data==0 a valid state for number 0. Incompatible with the cache. I
//   tried, and it doesn't seem to help (may even hurt a bit).
struct Mpzf {
  private:
#ifdef CGAL_MPZF_USE_CACHE
  // More experiments to determine the best value would be good. It likely
  // depends on the usage. Note that pool2 is fast enough that a conditional
  // cache slows it down. A purely static cache (crash if it isn't large
  // enough) still wins by about 11% on the Delaunay_3 construction, but is
  // more complicated to handle.
  // Evaluating a polynomial in double will never require more than roughly
  // (2100*degree) bits, or (33*degree) mp_limb_t, which is very small. I
  // checked by including an array of 150 limbs in every Mpzf (that's where
  // the 11% number comes from).
  // BONUS: doing that is thread-safe!
  static const unsigned int cache_size = 8;
#endif
//#if !defined(CGAL_HAS_THREADS) || defined(CGAL_I_PROMISE_I_WONT_USE_MANY_THREADS)
//  typedef Mpzf_impl::pool2<mp_limb_t*,Mpzf> pool;
//#elif defined(CGAL_CAN_USE_CXX11_THREAD_LOCAL)
//  typedef Mpzf_impl::pool3<mp_limb_t*,Mpzf> pool;
//#else
  typedef Mpzf_impl::no_pool<mp_limb_t*,Mpzf> pool;
//#endif

  mp_limb_t* data_; /* data_[0] is never 0 (except possibly for 0). */
  inline mp_limb_t*& data() { return data_; }
  inline mp_limb_t const* data() const { return data_; }

#ifdef CGAL_MPZF_USE_CACHE
  mp_limb_t cache[cache_size + 1];
#endif
  int size; /* Number of relevant limbs in data_. */
  int exp; /* The number is data_ (an integer) * 2 ^ (64 * exp). */
  typedef int Exponent_type;
  typedef int Size_type;

  struct allocate{};
  struct noalloc{};

  void init(unsigned mini=2){
#ifdef CGAL_MPZF_USE_CACHE
    if (mini <= cache_size) {
      cache[0] = cache_size;
      data() = cache + 1;
      return;
    }
#endif
    if(!pool::empty()){
      data() = pool::pop();
      if(data()[-1] >= mini) return; // TODO: when mini==2, no need to check
      delete[] (data() - (pool::extra+1)); // too small, useless
    }
    if(mini<2) mini=2;
    data() = (new mp_limb_t[mini+(pool::extra+1)]) + (pool::extra+1);
    data()[-1] = mini;
  }
  void clear(){
    // while(*--data()==0);
    // This line gave a misscompilation by Intel Compiler 2019
    // (19.0.0.117). I replaced it by the following two lines:
    // -- Laurent Rineau, sept. 2018
    --data();
    while(*data()==0) { --data(); } // in case we skipped final zeroes

#ifdef CGAL_MPZF_USE_CACHE
    if (data() == cache) return;
#endif
    ++data();
    pool::push(data());
  }

  Mpzf(noalloc){}
  Mpzf(allocate,int i) { init(i); }

  public:

  static void clear_pool () {
    while (!pool::empty())
      delete[] (pool::pop() - (pool::extra + 1));
  }

  ~Mpzf(){
    clear();
  }
  Mpzf(): size(0), exp(0) {
    init();
  }
  Mpzf& operator=(Mpzf const& x){
    unsigned asize=std::abs(x.size);
    if(asize==0) { exp=0; size=0; return *this; }
    if(this==&x) return *this;
    while(*--data()==0); // factor that code somewhere?
    if(*data()<asize){
#ifdef CGAL_MPZF_USE_CACHE
      if (data() != cache)
#endif
        delete[] (data() - pool::extra);
      init(asize);
    } else ++data();
    size=x.size;
    exp=x.exp;
    mpn_copyi(data(),x.data(),asize);
    return *this;
  }
  Mpzf(Mpzf const& x){
    int asize=std::abs(x.size);
    init(asize);
    size=x.size;
    exp=x.exp;
    if(size!=0) mpn_copyi(data(),x.data(),asize);
  }
#if defined(CGAL_MPZF_USE_CACHE)
  Mpzf(Mpzf&& x)noexcept:size(x.size),exp(x.exp){
    auto xd = x.data();
    while(*--xd==0);
    if (xd != x.cache) {
      data() = x.data();
      x.init();
    } else {
      init();
      if(size!=0) mpn_copyi(data(),x.data(),std::abs(size));
    }
    x.size = 0;
  }
  Mpzf& operator=(Mpzf&& x)noexcept{
    if (this == &x) return *this; // is this needed?
    size = x.size;
    exp = x.exp;
    auto xd = x.data();
    auto td = data();
    while(*--xd==0);
    while(*--td==0);
    if (xd != x.cache) {
      data() = x.data();
      if (td != cache) {
        pool::push(td+1);
        // should we instead give it to x in case x is reused?
        // x.data() = td + 1;
      }
      x.init();
    } else {
      // In some cases data points in the middle of the buffer, reset it
      data() = td + 1;
      if(size!=0) mpn_copyi(data(),x.data(),std::abs(size));
    }
    x.size = 0;
    return *this;
  }
#else
  Mpzf(Mpzf&& x):data_(x.data()),size(x.size),exp(x.exp){
    x.init(); // yes, that's a shame...
    x.size = 0;
    x.exp = 0;
  }
  Mpzf& operator=(Mpzf&& x)noexcept{
    size = x.size;
    // In case something tries to read it, size needs to be smaller than data
    x.size = 0;
    exp = x.exp;
    std::swap(data(),x.data());
    return *this;
  }
  friend void swap(Mpzf&a, Mpzf&b)noexcept{
    std::swap(a.size, b.size);
    std::swap(a.exp, b.exp);
    std::swap(a.data(), b.data());
  }
  friend Mpzf operator-(Mpzf&& x){
    Mpzf ret = std::move(x);
    ret.size = -ret.size;
    return ret;
  }
#endif
  Mpzf(int i) : exp(0) {
    // assume that int is smaller than mp_limb_t
    init();
    if      (i == 0)    { size = 0; }
    else if (i >  0)    { size = 1; data()[0] = i; }
    else /* (i <  0) */ { size =-1; data()[0] = -(mp_limb_t)i; }
    // cast to mp_limb_t because -INT_MIN is undefined
  }
  Mpzf(unsigned int i) : exp(0) {
    // assume that int is smaller than mp_limb_t
    init();
    if      (i == 0)    { size = 0; }
    else /* (i >  0) */ { size = 1; data()[0] = i; }
  }
  Mpzf(long i) : exp(0) {
    // assume that long is smaller than mp_limb_t
    init();
    if      (i == 0)    { size = 0; }
    else if (i >  0)    { size = 1; data()[0] = i; }
    else /* (i <  0) */ { size =-1; data()[0] = -(mp_limb_t)i; }
    // cast to mp_limb_t because -LONG_MIN is undefined
  }
  Mpzf(unsigned long i) : exp(0) {
    // assume that long is smaller than mp_limb_t
    init();
    if      (i == 0)    { size = 0; }
    else /* (i >  0) */ { size = 1; data()[0] = i; }
  }
  Mpzf(double d){
    init();
    using boost::uint64_t;
    union {
#ifdef CGAL_LITTLE_ENDIAN
      struct { uint64_t man:52; uint64_t exp:11; uint64_t sig:1; } s;
#else /* CGAL_BIG_ENDIAN */
      //WARNING: untested!
      struct { uint64_t sig:1; uint64_t exp:11; uint64_t man:52; } s;
#endif
      double d;
    } u;
    u.d = d;
    uint64_t m;
    uint64_t dexp = u.s.exp;
    CGAL_assertion_msg(dexp != 2047, "Creating an Mpzf from infinity or NaN.");
    if (dexp == 0) {
      if (d == 0) { size=0; exp=0; return; }
      else { // denormal number
        m = u.s.man;
        ++dexp;
      }
    } else {
      m = (1LL<<52) | u.s.man;
    }
    int e1 = (int)dexp+13;
    // FIXME: make it more general! But not slower...
    CGAL_static_assertion(GMP_NUMB_BITS == 64);
    int e2 = e1 % 64;
    exp = e1 / 64 - 17;
    // 52+1023+13==17*64 ?
#if 0
    // This seems very slightly faster
    if(Mpzf_impl::ctz(m)+e2>=64){
      data()[0] = m >> (64-e2);
      size = 1;
      ++exp;
    }else{
      data()[0] = m << e2;
      if(e2>11){ // Wrong test for denormals
        data()[1] = m >> (64-e2);
        size = 2;
      } else {
        size = 1;
      }
    }
#else
    mp_limb_t d0 = (m << e2) & GMP_NUMB_MASK;
    mp_limb_t d1 = 0;
    if (e2 != 0) // shifting by 64 is UB
      d1 = m >> (GMP_NUMB_BITS - e2);
    if (d0 == 0) {
      data()[0] = d1;
      size = 1;
      ++exp;
    }
    else {
      data()[0] = d0;
      if (d1 == 0) {
        size = 1;
      }
      else {
        data()[1] = d1;
        size = 2;
      }
    }
#endif
    if(u.s.sig) size=-size;
    //CGAL_assertion(to_double()==IA_force_to_double(d));
  }

#ifdef CGAL_USE_GMPXX
  Mpzf(mpz_class const&z){
    init_from_mpz_t(z.get_mpz_t());
  }
#endif
  Mpzf(Gmpz const&z){
    init_from_mpz_t(z.mpz());
  }
  void init_from_mpz_t(mpz_t const z){
    exp=Exponent_type(mpz_scan1(z,0)/GMP_NUMB_BITS);
    size=Size_type(mpz_size(z)-exp);
    init(size);
    mpn_copyi(data(),z->_mp_d+exp,size);
  }

#if 0
  // For debug purposes only
  void print()const{
    //std::cout << "size: " << size << std::endl;
    if(size==0) { std::cout << "zero\n"; return; }
    if(size<0) std::cout << "- ";
    int asize = std::abs(size);
    std::cout << std::hex;
    while(--asize>=0) { std::cout << data()[asize] << ' '; }
    std::cout << std::dec << "exp " << exp << ' ';
    std::cout << std::dec << "size " << size << ' ';
    asize = std::abs(size);
    std::cout << "double: " << std::ldexp((double)data()[asize-1],64*(exp+asize-1))*((size<0)?-1:1) << '\n';
  }
#endif

  friend int Mpzf_abscmp(Mpzf const&a, Mpzf const&b){
    int asize=std::abs(a.size);
    int bsize=std::abs(b.size);
    // size==0 should mean exp==-infinity, like with double.
    // Since it doesn't, test for it explicitly.
    if (bsize == 0) return asize;
    if (asize == 0) return -1;
    int ah=asize+a.exp;
    int bh=bsize+b.exp;
    if(ah!=bh) return ah-bh;
    int minsize=(std::min)(asize,bsize);
    const mp_limb_t* adata=a.data()+(asize-1);
    const mp_limb_t* bdata=b.data()+(bsize-1);
    for(int i=0;i<minsize;++i,--adata,--bdata){
      const mp_limb_t aa=*adata;
      const mp_limb_t bb=*bdata;
      if(aa!=bb) return (aa<bb)?-1:1;
    }
    return asize-bsize; // this assumes that we get rid of trailing zeros...
  }
  friend int Mpzf_cmp (Mpzf const&a, Mpzf const&b){
    if ((a.size ^ b.size) < 0) return (a.size < 0) ? -1 : 1;
    int res = Mpzf_abscmp(a, b);
    return (a.size < 0) ? -res : res;
  }
  friend bool operator<(Mpzf const&a, Mpzf const&b){
    if((a.size ^ b.size) < 0) return a.size < 0;
    return ((a.size < 0) ? Mpzf_abscmp(b, a) : Mpzf_abscmp(a, b)) < 0;
  }
  friend bool operator>(Mpzf const&a, Mpzf const&b){
    return b<a;
  }
  friend bool operator>=(Mpzf const&a, Mpzf const&b){
    return !(a<b);
  }
  friend bool operator<=(Mpzf const&a, Mpzf const&b){
    return !(a>b);
  }
  friend bool operator==(Mpzf const&a, Mpzf const&b){
    if (a.exp != b.exp || a.size != b.size) return false;
    if (a.size == 0) return true;
    return mpn_cmp(a.data(), b.data(), std::abs(a.size)) == 0;
  }
  friend bool operator!=(Mpzf const&a, Mpzf const&b){
    return !(a==b);
  }
  friend Mpzf const& min BOOST_PREVENT_MACRO_SUBSTITUTION (Mpzf const&a, Mpzf const&b){
    return (b<a)?b:a;
  }
  friend Mpzf const& max BOOST_PREVENT_MACRO_SUBSTITUTION (Mpzf const&a, Mpzf const&b){
    return (a<b)?b:a;
  }
  private:
  static Mpzf aors(Mpzf const&a, Mpzf const&b, int bsize){
    Mpzf res=noalloc();
    if(bsize==0){
      int size=std::abs(a.size);
      res.init(size);
      res.exp=a.exp;
      res.size=a.size;
      if(size!=0) mpn_copyi(res.data(),a.data(),size);
      return res;
    }
    int asize=a.size;
    if(asize==0){
      int size=std::abs(bsize);
      res.init(size);
      res.exp=b.exp;
      res.size=bsize;
      mpn_copyi(res.data(),b.data(),size);
      return res;
    }
    if((asize^bsize)>=0){
      // Addition
      int absasize=std::abs(asize);
      int absbsize=std::abs(bsize);
      const mp_limb_t* adata=a.data();
      const mp_limb_t* bdata=b.data();
      int aexp=a.exp;
      int bexp=b.exp;
      if(aexp<bexp){ res.exp=a.exp; aexp=0; bexp=b.exp-a.exp; }
      else { res.exp=b.exp; aexp=a.exp-b.exp; bexp=0; }
      res.init((std::max)(absasize+aexp,absbsize+bexp)+1);
      mp_limb_t* rdata=res.data();
      res.size=0;
      // TODO: if aexp>0, swap a and b so we don't repeat the code.
      if(0<bexp){
        if(absasize<=bexp){ // no overlap
          mpn_copyi(rdata, adata, absasize);
          rdata+=absasize;
          rdata=Mpzf_impl::fill_n_ptr(rdata,bexp-absasize,0);
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
          rdata=Mpzf_impl::fill_n_ptr(rdata,aexp-absbsize,0);
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
      while(/*res.size>0&&*/res.data()[0]==0){--res.size;++res.data();++res.exp;}
      if(bsize<0) res.size=-res.size;
    } else {
      // Subtraction
      const Mpzf *x, *y;
      int xsize=a.size;
      int ysize=bsize;
      int cmp=Mpzf_abscmp(a,b);
      if(cmp==0){ res.init(); res.size=0; res.exp=0; return res; }
      if(cmp<0) { x=&b; y=&a; std::swap(xsize, ysize); }
      else { x=&a; y=&b; }
      int absxsize=std::abs(xsize);
      int absysize=std::abs(ysize);
      const mp_limb_t* xdata=x->data();
      const mp_limb_t* ydata=y->data();
      int xexp=x->exp;
      int yexp=y->exp;
      if(xexp<yexp){ res.exp=xexp; yexp-=xexp; xexp=0; }
      else { res.exp=yexp; xexp-=yexp; yexp=0; }
      res.init((std::max)(absxsize+xexp,absysize+yexp)+1);
      mp_limb_t* rdata=res.data();
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
          rdata=Mpzf_impl::fill_n_ptr(rdata,xexp-absysize,-1);
          mpn_sub_1(rdata, xdata, absxsize, 1);
          res.size=absxsize+xexp;
          while(/*res.size>0&&*/res.data()[res.size-1]==0) --res.size;
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
      CGAL_assertion_code( mp_limb_t carry= )
        mpn_sub(rdata, xdata, absxsize, ydata, absysize);
      if(carry1)
        CGAL_assertion_code( carry+= )
          mpn_sub_1(rdata, rdata, absxsize, 1);
      CGAL_assertion(carry==0);
      res.size+=absxsize;
      while(/*res.size>0&&*/res.data()[res.size-1]==0) --res.size;
      while(/*res.size>0&&*/res.data()[0]==0){--res.size;++res.data();++res.exp;}
      if(xsize<0) res.size=-res.size;
    }
    return res;
  }

  public:
  friend Mpzf operator+(Mpzf const&a, Mpzf const&b){
    return aors(a,b,b.size);
  }

  friend Mpzf operator-(Mpzf const&a, Mpzf const&b){
    return aors(a,b,-b.size);
  }

  friend Mpzf operator*(Mpzf const&a, Mpzf const&b){
    int asize=std::abs(a.size);
    int bsize=std::abs(b.size);
    int siz=asize+bsize;
    Mpzf res(allocate(),siz);
    if(asize==0||bsize==0){res.exp=0;res.size=0;return res;}
    res.exp=a.exp+b.exp;
    mp_limb_t high;
    if(asize>=bsize)
      high = mpn_mul(res.data(),a.data(),asize,b.data(),bsize);
    else
      high = mpn_mul(res.data(),b.data(),bsize,a.data(),asize);
    if(high==0) --siz;
    if(res.data()[0]==0) { ++res.data(); ++res.exp; --siz; }
    res.size=((a.size^b.size)>=0)?siz:-siz;
    return res;
  }

  friend Mpzf Mpzf_square(Mpzf const&a){
    int asize=std::abs(a.size);
    int siz=2*asize;
    Mpzf res(allocate(),siz);
    res.exp=2*a.exp;
    if(asize==0){res.size=0;return res;}
    mpn_sqr(res.data(),a.data(),asize);
    mp_limb_t high = res.data()[siz-1];
    if(high==0) --siz;
    if(res.data()[0]==0) { ++res.data(); ++res.exp; --siz; }
    res.size=siz;
    return res;
  }

  friend Mpzf operator/(Mpzf const&a, Mpzf const&b){
    // FIXME: Untested
    int asize=std::abs(a.size);
    int bsize=std::abs(b.size);
    int siz=asize+2-bsize;
    Mpzf res(allocate(),asize+2);
    if(bsize==0){throw std::range_error("Division by zero");}
    if(asize==0){res.exp=0;res.size=0;return res;}
    res.size=siz;
    res.exp=a.exp-b.exp;
    const mp_limb_t *adata = a.data();
    const mp_limb_t *bdata = b.data();
    mp_limb_t *qp = res.data();
    mp_limb_t *rp = qp + siz;
    if(Mpzf_impl::ctz(adata[0]) >= Mpzf_impl::ctz(bdata[0])){ // Easy case
      --res.size;
      mpn_tdiv_qr(qp, rp, 0, adata, asize, bdata, bsize);
      CGAL_assertion_code(
          for (int i=0; i<bsize; ++i)
            if (rp[i] != 0) throw std::logic_error("non exact Mpzf division");
      )
    }
    else if(adata[-1]==0){ // We are lucky
      --adata; ++asize; --res.exp;
      mpn_tdiv_qr(qp, rp, 0, adata, asize, bdata, bsize);
      CGAL_assertion_code(
          for (int i=0; i<bsize; ++i)
            if (rp[i] != 0) throw std::logic_error("non exact Mpzf division");
      )
    }
    else{
      --res.exp;
      Mpzf a2(allocate(),asize+1);
      a2.data()[0]=0;
      mpn_copyi(a2.data()+1,a.data(),asize);
      // No need to complete a2, we just want the buffer.
      //a2.size=(a.size<0)?(a.size-1):(a.size+1);
      //a2.exp = a.exp-1;
      mpn_tdiv_qr(qp, rp, 0, a2.data(), asize+1, bdata, bsize);
      CGAL_assertion_code(
          for (int i=0; i<bsize; ++i)
            if (rp[i] != 0) throw std::logic_error("non exact Mpzf division");
      )
    }
    while(/*res.size>0&&*/res.data()[res.size-1]==0) --res.size;
    //while(/*res.size>0&&*/res.data()[0]==0){--res.size;++res.data();++res.exp;}
    if((a.size^b.size)<0) res.size=-res.size;
    return res;
  }

  friend Mpzf Mpzf_gcd(Mpzf const&a, Mpzf const&b){
    // FIXME: Untested
    if (a.size == 0) return b;
    if (b.size == 0) return a;
    int asize=std::abs(a.size);
    int bsize=std::abs(b.size);
    int atz=Mpzf_impl::ctz(a.data()[0]);
    int btz=Mpzf_impl::ctz(b.data()[0]);
    int rtz=(std::min)(atz,btz);
    Mpzf tmp(allocate(), asize);
    Mpzf res(allocate(), bsize);
    if (atz != 0) {
      mpn_rshift(tmp.data(), a.data(), asize, atz);
      if(tmp.data()[asize-1]==0) --asize;
    }
    else { mpn_copyi(tmp.data(), a.data(), asize); }
    if (btz != 0) {
      mpn_rshift(res.data(), b.data(), bsize, btz);
      if(res.data()[bsize-1]==0) --bsize;
    }
    else { mpn_copyi(res.data(), b.data(), bsize); }
    res.exp = 0; // Pick b.exp? or the average? 0 helps return 1 more often.
    if (asize < bsize)
      res.size = Size_type(mpn_gcd(res.data(), res.data(), bsize, tmp.data(), asize));
    else
      res.size = Size_type(mpn_gcd(res.data(), tmp.data(), asize, res.data(), bsize));
    if(rtz!=0) {
      mp_limb_t c = mpn_lshift(res.data(), res.data(), res.size, rtz);
      if(c) { res.data()[res.size]=c; ++res.size; }
    }
    return res;
  }

  friend bool Mpzf_is_square(Mpzf const&x){
    if (x.size < 0) return false;
    if (x.size == 0) return true;
    // Assume that GMP_NUMB_BITS is even.
    return mpn_perfect_square_p (x.data(), x.size);
  }

  friend Mpzf Mpzf_sqrt(Mpzf const&x){
    // FIXME: Untested
    if (x.size < 0) throw std::range_error("Sqrt of negative number");
    if (x.size == 0) return 0;
    if (x.exp % 2 == 0) {
      Mpzf res(allocate(), (x.size + 1) / 2);
      res.exp = x.exp / 2;
      res.size = (x.size + 1) / 2;
      CGAL_assertion_code(mp_size_t rem=)
      mpn_sqrtrem(res.data(), 0, x.data(), x.size);
      CGAL_assertion(rem==0);
      return res;
    }
    else if (x.data()[-1] == 0) {
      Mpzf res(allocate(), (x.size + 2) / 2);
      res.exp = (x.exp - 1) / 2;
      res.size = (x.size + 2) / 2;
      CGAL_assertion_code(mp_size_t rem=)
      mpn_sqrtrem(res.data(), 0, x.data()-1, x.size+1);
      CGAL_assertion(rem==0);
      return res;
    }
    else {
      Mpzf res(allocate(), (x.size + 2) / 2);
      res.exp = (x.exp - 1) / 2;
      res.size = (x.size + 2) / 2;
      CGAL_assertion_code(mp_size_t rem=)
      mpn_sqrtrem(res.data(), 0, x.data(), x.size);
      CGAL_assertion(rem==0);
      mpn_lshift(res.data(), res.data(), res.size, GMP_NUMB_BITS / 2);
      return res;
    }
  }

  friend Mpzf operator+(Mpzf const&x){
    return x;
  }
  friend Mpzf operator-(Mpzf const&x){
    Mpzf ret = x;
    ret.size = -ret.size;
    return ret;
  }
  Mpzf& operator+=(Mpzf const&x){ *this=*this+x; return *this; }
  Mpzf& operator-=(Mpzf const&x){ *this=*this-x; return *this; }
  Mpzf& operator*=(Mpzf const&x){ *this=*this*x; return *this; }
  Mpzf& operator/=(Mpzf const&x){ *this=*this/x; return *this; }

  bool is_canonical () const {
    if (size == 0) return true;
    if (data()[0] == 0) return false;
    if (data()[std::abs(size)-1] == 0) return false;
    return true;
  }

  bool is_zero () const {
    return size==0;
  }

  bool is_one () const {
    return exp==0 && size==1 && data()[0]==1;
  }

  CGAL::Sign sign () const { return CGAL::sign(size); }

  double to_double () const {
    // Assumes GMP_NUMB_BITS == 64
    using std::ldexp;
    if(size==0) return 0;
    int asize = std::abs(size);
    mp_limb_t top = data()[asize-1];
    double dtop = (double)top;
    if(top >= (1LL<<53) || asize == 1) /* ok */ ;
    else { dtop += (double)data()[asize-2] * ldexp(1.,-GMP_NUMB_BITS); }
    return ldexp( (size<0) ? -dtop : dtop, (asize-1+exp) * GMP_NUMB_BITS);
  }

  std::pair<double, double> to_interval () const {
    // Assumes GMP_NUMB_BITS == 64
    if (size == 0) return std::make_pair(0., 0.);
    double dl, dh;
    int asize = std::abs(size);
    int e = 64 * (asize - 1 + exp);
    mp_limb_t x = data()[asize-1];
    int lz = Mpzf_impl::clz(x);
    if (lz <= 11) {
      if (lz != 11) {
        e += (11 - lz);
        x >>= (11 - lz);
      }
      dl = double(x);
      dh = double(x + 1);
      // Check for the few cases where dh=x works (asize==1 and the evicted
      // bits from x were 0s)
    }
    else if (asize == 1) {
      dl = dh = double(x); // conversion is exact
    }
    else {
      mp_limb_t y = data()[asize-2];
      e -= (lz - 11);
      x <<= (lz - 11);
      y >>= (75 - lz);
      x |= y;
      dl = double(x);
      dh = double(x + 1);
      // Check for the few cases where dh=x works (asize==2 and the evicted
      // bits from y were 0s)
    }
    typedef Interval_nt<> IA;
    IA res (dl, dh);
    res = ldexp (res, e);
    if (size < 0) res = -res;
    return CGAL::to_interval(res);
    // Use ldexp(Interval_nt,int) to delegate the hard thinking
    // about over/underflow.
  }

#ifdef CGAL_USE_GMPXX
  explicit
  operator mpq_class () const {
    mpq_class q;
    export_to_mpq_t(q.get_mpq_t());
    return q;
  }
#endif

  explicit
  operator Gmpq () const {
    Gmpq q;
    export_to_mpq_t(q.mpq());
    return q;
  }
  void export_to_mpq_t(mpq_t q) const {
    /* q must be 0/1 before this call */
    CGAL_precondition(mpq_cmp_ui(q,0,1)==0);
    if (size != 0) {
      mpz_import (mpq_numref (q),
                  std::abs(size),
                  -1, // data()[0] is the least significant part
                  sizeof(mp_limb_t),
                  0, // native endianness inside mp_limb_t
                  GMP_NAIL_BITS, // should be 0
                  data());
      if (exp > 0)
        mpq_mul_2exp(q, q, (sizeof(mp_limb_t) * CHAR_BIT *  exp));
      else if (exp < 0)
        mpq_div_2exp(q, q, (sizeof(mp_limb_t) * CHAR_BIT * -exp));

      if (size < 0)
        mpq_neg(q,q);
    }
  }
#if 0
  explicit
// This makes Mpzf==int ambiguous
  operator Gmpzf () const {
    mpz_t z;
    z->_mp_d=const_cast<mp_limb_t*>(data());
    z->_mp_size=size;
    Gmpzf m(z);
// Only works for a very limited range of exponents
    Gmpzf e(std::ldexp(1.,GMP_NUMB_BITS*exp));
    return m*e;
  }
#endif

  friend void simplify_quotient(Mpzf& a, Mpzf& b){
    // Avoid quotient(2^huge_a/2^huge_b)
    a.exp -= b.exp;
    b.exp = 0;
    // Simplify with gcd?
  }
};

// Copied from Gmpzf, not sure that's the best thing to do.
inline
std::ostream& operator<< (std::ostream& os, const Mpzf& a)
{
    return os << a.to_double();
}

inline
std::istream& operator>> (std::istream& is, Mpzf& a)
{
  double d;
  is >> d;
  if (is)
    a = d;
  return is;
}


  template <> struct Algebraic_structure_traits< Mpzf >
    : public Algebraic_structure_traits_base< Mpzf, Integral_domain_without_division_tag >  {
      typedef Tag_true            Is_exact;
      typedef Tag_false            Is_numerical_sensitive;

      struct Is_zero
        : public CGAL::cpp98::unary_function< Type, bool > {
          bool operator()( const Type& x ) const {
            return x.is_zero();
          }
        };

      struct Is_one
        : public CGAL::cpp98::unary_function< Type, bool > {
          bool operator()( const Type& x ) const {
            return x.is_one();
          }
        };

      struct Gcd
        : public CGAL::cpp98::binary_function< Type, Type, Type > {
          Type operator()(
              const Type& x,
              const Type& y ) const {
            return Mpzf_gcd(x, y);
          }
        };

      struct Square
        : public CGAL::cpp98::unary_function< Type, Type > {
          Type operator()( const Type& x ) const {
            return Mpzf_square(x);
          }
        };

      struct Integral_division
        : public CGAL::cpp98::binary_function< Type, Type, Type > {
          Type operator()(
              const Type& x,
              const Type& y ) const {
            return x / y;
          }
        };

      struct Sqrt
        : public CGAL::cpp98::unary_function< Type, Type > {
          Type operator()( const Type& x) const {
            return Mpzf_sqrt(x);
          }
        };

      struct Is_square
        : public CGAL::cpp98::binary_function< Type, Type&, bool > {
          bool operator()( const Type& x, Type& y ) const {
            // TODO: avoid doing 2 calls.
            if (!Mpzf_is_square(x)) return false;
            y = Mpzf_sqrt(x);
            return true;
          }
          bool operator()( const Type& x) const {
            return Mpzf_is_square(x);
          }
        };

    };
  template <> struct Real_embeddable_traits< Mpzf >
    : public INTERN_RET::Real_embeddable_traits_base< Mpzf , CGAL::Tag_true > {
      struct Sgn
        : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign > {
          ::CGAL::Sign operator()( const Type& x ) const {
            return x.sign();
          }
        };

      struct To_double
        : public CGAL::cpp98::unary_function< Type, double > {
            double operator()( const Type& x ) const {
              return x.to_double();
            }
        };

      struct Compare
        : public CGAL::cpp98::binary_function< Type, Type, Comparison_result > {
            Comparison_result operator()(
                const Type& x,
                const Type& y ) const {
              return CGAL::sign(Mpzf_cmp(x,y));
            }
        };

      struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
            std::pair<double, double> operator()( const Type& x ) const {
              return x.to_interval();
            }
        };

    };

CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Mpzf)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,Mpzf)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,Mpzf)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,Mpzf)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float    ,Mpzf)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double   ,Mpzf)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpz     ,Mpzf)
#ifdef CGAL_USE_GMPXX
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(mpz_class,Mpzf)
#endif

}

/* There isn't much Eigen can do with such a type,
 * mostly this is here for IsInteger to protect people.
 */
namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CGAL::Mpzf>
  {
    typedef CGAL::Mpzf Real;
    /* Should this be Quotient<Mpzf>? Gmpq?  */
    typedef CGAL::Mpzf NonInteger;
    typedef CGAL::Mpzf Nested;
    typedef CGAL::Mpzf Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      /* Only exact divisions are supported, close enough to an integer.
       * This way we get compilation failures instead of runtime.  */
      IsInteger = 1,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 30,
      MulCost = 50
    };
  };
}

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // GMP_NUMB_BITS == 64
#endif // CGAL_MPZF_H
