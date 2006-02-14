// a few useful functions for usage of GeoWin with the rat-kernel-traits ...

#ifndef CEP_GEOWIN_LEDA_RAT_KERNEL_H
#define CEP_GEOWIN_LEDA_RAT_KERNEL_H

#include <CGAL/basic.h>
#include <LEDA/geowin.h>
#include <LEDA/geowin_init.h>
#include <list>

#if !defined(LEDA_NAMESPACE)
#define LEDA_BEGIN_NAMESPACE
#define LEDA_END_NAMESPACE
#endif

LEDA_BEGIN_NAMESPACE

template<class T>
void geowin_generate_objects(GeoWin& gw, std::list<T>& L)
{
  leda_list<T> H;
  geowin_generate_objects(gw,H);
  T obj;
  forall(obj,H) L.push_back(obj);
}


template<class T>
leda_string geowin_info_fcn(const std::list<T>& L)
{ int n = L.size();
  if (n == 1)
    return leda_string("\\black\\tt %d %s", n, leda_tname((T*)0));
  else
    return leda_string("\\black\\tt %d %ss",n, leda_tname((T*)0));
}

LEDA_END_NAMESPACE



#endif
