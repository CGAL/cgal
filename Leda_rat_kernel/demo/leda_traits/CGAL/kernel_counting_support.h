// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/kernel_counting_support.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================


#if !defined(CGAL_DO_COUNTING)
#define CGAL_DO_COUNTING

#include <CGAL/basic.h>
#include <map>
#include <cstring>
#include <typeinfo>

CGAL_BEGIN_NAMESPACE

// a map is used for counting ...
// why do we use a template ?
// we want to have different instantiations for different kernels ...

template<class T>
struct Do_counting {

static std::map<std::string,int> COUNTER;

template<class O1,class O2,class A1>
void operator()(const O1& func1, const O2& func2, const A1& arg1) const
{
 std::string str = typeid(func1).name();
 
#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else  COUNTER[str] = COUNTER[str] + 1;
}

template<class O1,class O2,class A1,class A2>
void operator()(const O1& func1, const O2& func2, const A1& arg1, const A2& arg2) const
{
 std::string str = typeid(func1).name();

#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1;

}

template<class O1,class O2,class A1,class A2,class A3>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3) const
{
 std::string str = typeid(func1).name();
 
#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1; 
}

template<class O1,class O2,class A1,class A2,class A3,class A4>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4) const
{
 std::string str = typeid(func1).name();

#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1;
}

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5) const
{
 std::string str = typeid(func1).name();
 
#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1; 
}

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6) const
{
 std::string str = typeid(func1).name();
 
#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1; 
}

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6,class A7>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6,
		const A7& arg7) const
{
 std::string str = typeid(func1).name();
 
#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1; 
}

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6,
		const A7& arg7, const A8& arg8) const
{
 std::string str = typeid(func1).name();
 
#if defined(DEBUG) 
  std::cout << str << "\n";
#endif 
 
 if (COUNTER.find(str) == COUNTER.end()) COUNTER[str] = 1;
 else COUNTER[str] = COUNTER[str] + 1; 
}

static void print_counters(std::ostream& os)
{
 std::map<std::string,int>::iterator iter = COUNTER.begin();
 for(;iter != COUNTER.end(); iter++){
   os << iter->first << ":    " << iter->second << "\n";
 } 
 
//#ifdef __GNUG__
// os <<  __PRETTY_FUNCTION__ << "\n";
//#endif 
}

static void clear_counters() { COUNTER.clear(); }

};

template<class T>
std::map<std::string,int> Do_counting<T>::COUNTER;

CGAL_END_NAMESPACE

#endif
