// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 11
//
// file          : include/CGAL/R_Tree/examples/R_tree_key.h
// package       : ExternalMemoryStructures (0.631)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// chapter       : $CGAL_Chapter: Basic / External Data Structures $
// source        : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Peter Widmayer <widmayer@inf.ethz.ch>)
//
// Example classes that have the interfaces of Range_tree_interface.h
// ======================================================================
#ifndef __R_tree_key_H__
#define __R_tree_key_H__
#include <CGAL/basic.h>
#include <iostream>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

// a 1 dimensional key which is an interval
class R_tree_key_1
{
  int min(int a, int b) const { return (a<b) ? a : b; }
  int max(int a, int b) const { return (a<b) ? b : a; }
public:
  int xmin, xmax;
  size_t size(void) const {
    return sizeof(*this);
  }

  R_tree_key_1(){xmin=xmax=0;}
  R_tree_key_1(int x1, int x2){
    xmin=x1;
    xmax=x2;
  }
  
  R_tree_key_1(const R_tree_key_1 &t){
    xmin=t.xmin;
    xmax=t.xmax;
  }
      
  bool operator==(const R_tree_key_1 &p) const {
    return ((xmin==p.xmin) && (xmax == p.xmax));
  }

  R_tree_key_1 & operator=(const R_tree_key_1 &t)  {
    xmin=t.xmin;
    xmax=t.xmax;
    return *this;
  }

  void unify( const R_tree_key_1 & p, const R_tree_key_1 & q){
    xmin = min(p.xmin, q.xmin);
    xmax = max(p.xmax, q.xmax);
  }

  //returns true if *this includes y
  bool include(const R_tree_key_1& y) const {
    if ((xmin > y.xmin))
      return false;
    if ((xmax < y.xmax))
      return false;
    return true;
  }

  bool compare(const R_tree_key_1 & y) const{
    return intersect(y);
  }


  bool intersect( const R_tree_key_1& y) const {
    if ((xmax <= y.xmin) || (y.xmax <= xmin)) {
      return false;
    }
    return true;
  } 



  double cost() const {
    return (xmax - xmin);
  }


  void intersection(const R_tree_key_1 & p, const R_tree_key_1 & q){
    if(p.intersect(q))
      {
	xmin = max(p.xmin, q.xmin);
	xmax = min(p.xmax, q.xmax);
      }
  }


  // compute the distance between the centers of the keys
  double center_dist( const R_tree_key_1& q) const {
    double x=xmin + 0.5*(xmax-xmin);
    double qx=q.xmin + 0.5*(q.xmax-q.xmin);
    double dist = (x-qx)*(x-qx);
    return dist;
  }

  void read(char **s) {
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    xmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    xmax=*((int *)from_int);

    delete[] from_int;
  }
  void write(char **s)
    {
      int sint=(int)sizeof(int);
      char *from_int=new char[sint];
      int i,r=0;
      memcpy(from_int,(char *)(&xmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i] = from_int[i];
      r += sint;

      memcpy(from_int,(char *)(&xmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      delete[] from_int;
    }

  class C_Compare_1{
  public:
    bool operator()(const R_tree_key_1 &p, const R_tree_key_1 &q)
    {
      if(p.xmin < q.xmin)
        return true;
      else 
	if(p.xmin==q.xmin)
	  if(p.xmax < q.xmax)
	    return true;
      return false;
    }
  };

  typedef C_Compare_1 compare_1;

  void dump(int depth=0) const	{
    std::cerr << "(" << xmin << "),(" << xmax << ")";
  }

};

std::ostream& operator << (std::ostream& s, R_tree_key_1 m) {
  s  << "(" << m.xmin << ")(" << m.xmax << ")" << std::endl;
  return s;
}


// a 2 dimensional key which is a 2 dimensional rectangle.
class R_tree_key_2
{
  int min(int a, int b) const { return (a<b) ? a : b; }
  int max(int a, int b) const { return (a<b) ? b : a; }
public:
  int xmin, xmax, ymin, ymax;
  size_t size(void) const {
    return 4*sizeof(int)+4;
//    return sizeof(*this);
  }

  R_tree_key_2(){xmin=xmax=ymin=ymax=0;}
  R_tree_key_2(int x1, int x2, int y1, int y2){
    xmin=x1;
    xmax=x2;
    ymin=y1;
    ymax=y2;
  }
  
  R_tree_key_2(const R_tree_key_2 &t){
    xmin=t.xmin;
    xmax=t.xmax;
    ymin=t.ymin;
    ymax=t.ymax;
  }
      
  bool operator==(const R_tree_key_2 &p) const {
    return ((xmin==p.xmin) && (xmax == p.xmax) 
            && (ymin == p.ymin) &&(ymax == p.ymax));
  }

  R_tree_key_2 & operator=(const R_tree_key_2 &t)  {
    xmin=t.xmin;
    xmax=t.xmax;
    ymin=t.ymin;
    ymax=t.ymax;
    return *this;
  }

  class C_Compare_1{
  public:
    bool operator()(const R_tree_key_2 &p, const R_tree_key_2 &q)
    {
      if(p.xmin < q.xmin)
        return true;
      else 
	if(p.xmin==q.xmin)
	  if(p.xmax < q.xmax)
	    return true;
      return false;
    }
  };

  class C_Compare_2{
  public:
    bool operator()(const R_tree_key_2 &p, const R_tree_key_2 &q)
    {
      if(p.ymin < q.ymin)
        return true;
      else 
	if(p.ymin==q.ymin)
	  if(p.ymax < q.ymax)
	    return true;
      return false;
    }
  };
  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;

   void unify( const R_tree_key_2 & p, const R_tree_key_2 & q){
    xmin = min(p.xmin, q.xmin);
    xmax = max(p.xmax, q.xmax);
    ymin = min(p.ymin, q.ymin);
    ymax = max(p.ymax, q.ymax);
  }

  //returns true if *this includes y
  bool include(const R_tree_key_2& y) const {
    if ((xmin > y.xmin) || (ymin >  y.ymin))
      return false;
    if ((xmax < y.xmax) || (ymax < y.ymax))
      return false;
    return true;
  }

  bool compare(const R_tree_key_2 & y) const{
    return intersect(y);
  }


  bool intersect( const R_tree_key_2& y) const {
    if ((xmax <= y.xmin) || (y.xmax <= xmin)) {
      return false;
    }
    if ((ymax <= y.ymin) || (y.ymax <= ymin)) {
      return false;
    }
    return true;
  } 


  double cost() const {
    return (xmax - xmin) * (ymax - ymin);
  }

  void intersection( const R_tree_key_2 & p, const R_tree_key_2 & q){
    if(p.intersect(q))
      {
	xmin = max(p.xmin, q.xmin);
	xmax = min(p.xmax, q.xmax);
	ymin = max(p.ymin, q.ymin);
	ymax = min(p.ymax, q.ymax);
      }
  }


  // compute the distance between the centers of the keys
  double center_dist( const R_tree_key_2& q) const {
    double x=xmin + 0.5*(xmax-xmin);
    double y=ymin + 0.5*(ymax-ymin);
    double qx=q.xmin + 0.5*(q.xmax-q.xmin);
    double qy=q.ymin + 0.5*(q.ymax-q.ymin);
    double dist = (x-qx)*(x-qx) + (y-qy)*(y-qy); 
    return dist;
  }

  void read(char **s) {
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    xmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    xmax=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    ymin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    ymax=*((int *)from_int);

    delete[] from_int;
  }
  void write(char **s)
    {
      int sint=(int)sizeof(int);
      char *from_int=new char[sint];
      int i,r=0;
      memcpy(from_int,(char *)(&xmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i] = from_int[i];
      r += sint;

      memcpy(from_int,(char *)(&xmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&ymin),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&ymax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];

      delete[] from_int;
    }

  void dump(int depth=0) const	{
    std::cerr << "(" << xmin << "/" << ymin 
         << "),(" << xmax << "/" << ymax << ")";
  }

};

std::ostream& operator << (std::ostream& s, R_tree_key_2 m) {
  s  << "(" << m.xmin << "/" << m.ymin << ")("
     << m.xmax << "/" << m.ymax << ")" << std::endl;
  return s;
}





// a 3 dimensional key which is a 3 dimensional rectangle.
class R_tree_key_3
{
  int min(int a, int b) const { return (a<b) ? a : b; }
  int max(int a, int b) const { return (a<b) ? b : a; }
public:
  int xmin, xmax, ymin, ymax, zmin, zmax;
  size_t size(void) const {
    return sizeof(*this);
  }

  R_tree_key_3(){xmin=xmax=ymin=ymax=zmin=zmax=0;}
  R_tree_key_3(int x1, int x2, int y1, int y2, int z1, int z2){
    xmin=x1;
    xmax=x2;
    ymin=y1;
    ymax=y2;
    zmin=z1;
    zmax=z2;
  }
  
  R_tree_key_3(const R_tree_key_3 &t){
    xmin=t.xmin;
    xmax=t.xmax;
    ymin=t.ymin;
    ymax=t.ymax;
    zmin=t.zmin;
    zmax=t.zmax;
  }
      
  bool operator==(const R_tree_key_3 &p) const {
    return ((xmin==p.xmin) && (xmax == p.xmax) 
            && (ymin == p.ymin) &&(ymax == p.ymax)
	    && (zmin == p.zmin) &&(zmax == p.zmax));
  }

  R_tree_key_3 & operator=(const R_tree_key_3 &t)  {
    xmin=t.xmin;
    xmax=t.xmax;
    ymin=t.ymin;
    ymax=t.ymax;    
    zmin=t.zmin;
    zmax=t.zmax;
    return *this;
  }

  void unify( const R_tree_key_3 & p, const R_tree_key_3 & q){
    xmin = min(p.xmin, q.xmin);
    xmax = max(p.xmax, q.xmax);
    ymin = min(p.ymin, q.ymin);
    ymax = max(p.ymax, q.ymax);
    zmin = min(p.zmin, q.zmin);
    zmax = max(p.zmax, q.zmax);
  }

  //returns true if *this includes y
  bool include(const R_tree_key_3& y) const {
    if ((xmin > y.xmin) || (ymin >  y.ymin))
      return false;
    if ((xmax < y.xmax) || (ymax < y.ymax))
      return false;
    if ((zmax < y.zmax) || (zmax < y.zmax))
      return false;
    return true;
  }

  bool compare(const R_tree_key_3 & y) const{
    return intersect(y);
  }


  bool intersect( const R_tree_key_3& y) const {
    if ((xmax <= y.xmin) || (y.xmax <= xmin)) {
      return false;
    }
    if ((ymax <= y.ymin) || (y.ymax <= ymin)) {
      return false;
    }
    if ((zmax <= y.zmin) || (y.zmax <= zmin)) {
      return false;
    }
    return true;
  } 


  double cost() const {
    return (xmax - xmin) * (ymax - ymin)* (zmax - zmin);
  }

  void intersection( const R_tree_key_3 & p, const R_tree_key_3 & q){
    if(p.intersect(q))
      {
	xmin = max(p.xmin, q.xmin);
	xmax = min(p.xmax, q.xmax);
	ymin = max(p.ymin, q.ymin);
	ymax = min(p.ymax, q.ymax);
	zmin = max(p.zmin, q.zmin);
	zmax = min(p.zmax, q.zmax);
      }
  }


  // compute the distance between the centers of the keys
  double center_dist( const R_tree_key_3& q) const {
    double x=xmin + 0.5*(xmax-xmin);
    double y=ymin + 0.5*(ymax-ymin);
    double z=zmin + 0.5*(zmax-zmin);
    double qx=q.xmin + 0.5*(q.xmax-q.xmin);
    double qy=q.ymin + 0.5*(q.ymax-q.ymin);
    double qz=q.zmin + 0.5*(q.zmax-q.zmin);
    double dist = (x-qx)*(x-qx) + (y-qy)*(y-qy) + (z-qz)*(z-qz); 
    return dist;
  }

  void read(char **s) {
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    xmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    xmax=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    ymin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    ymax=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    zmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    zmax=*((int *)from_int);

    delete[] from_int;
  }

  void write(char **s)
    {
      int sint=(int)sizeof(int);
      char *from_int=new char[sint];
      int i,r=0;
      memcpy(from_int,(char *)(&xmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i] = from_int[i];
      r += sint;

      memcpy(from_int,(char *)(&xmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&ymin),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&ymax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&zmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&zmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];

      delete[] from_int;
    }

  class C_Compare_1{
  public:
    bool operator()(const R_tree_key_3 &p, const R_tree_key_3 &q)
    {
      if(p.xmin < q.xmin)
        return true;
      else 
	if(p.xmin==q.xmin)
	  if(p.xmax < q.xmax)
	    return true;
      return false;
    }
  };

  class C_Compare_2{
  public:
    bool operator()(const R_tree_key_3 &p, const R_tree_key_3 &q)
    {
      if(p.ymin < q.ymin)
        return true;
      else 
	if(p.ymin==q.ymin)
	  if(p.ymax < q.ymax)
	    return true;
      return false;
    }
  };

  class C_Compare_3{
  public:
    bool operator()(const R_tree_key_3 &p, const R_tree_key_3 &q)
    {
      if(p.ymin < q.ymin)
        return true;
      else 
	if(p.ymin==q.ymin)
	  if(p.ymax < q.ymax)
	    return true;
      return false;
    }
  };
  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;
  typedef C_Compare_3 compare_3;


  void dump(int depth=0) const	{
    std::cerr << "(" << xmin << "/" << ymin << "/" << zmin 
         << "),(" << xmax << "/" << ymax << "/" << zmax << ")";
  }

};

std::ostream& operator << (std::ostream& s, R_tree_key_3 m) {
  s  << "(" << m.xmin << "/" << m.ymin <<  "/" << m.zmin << ")("
     << m.xmax << "/" << m.ymax <<  "/" << m.zmax << ")" << std::endl;
  return s;
}


// a 4 dimensional key which is a 4 dimensional rectangle.
class R_tree_key_4
{
  int min(int a, int b) const { return (a<b) ? a : b; }
  int max(int a, int b) const { return (a<b) ? b : a; }
public:
  int xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
  size_t size(void) const {
    return sizeof(*this);
  }

  R_tree_key_4(){xmin=xmax=ymin=ymax=zmin=zmax=0;}
  R_tree_key_4(int x1, int x2, int y1, int y2, int z1, int z2, int w1, int w2){
    xmin=x1;
    xmax=x2;
    ymin=y1;
    ymax=y2;
    zmin=z1;
    zmax=z2;
    wmin=w1;
    wmax=w2;
  }
  
  R_tree_key_4(const R_tree_key_4 &t){
    xmin=t.xmin;
    xmax=t.xmax;
    ymin=t.ymin;
    ymax=t.ymax;
    zmin=t.zmin;
    zmax=t.zmax;
    wmin=t.wmin;
    wmax=t.wmax;
  }
      
  bool operator==(const R_tree_key_4 &p) const {
    return ((xmin==p.xmin) && (xmax == p.xmax) 
            && (ymin == p.ymin) &&(ymax == p.ymax)
	    && (zmin == p.zmin) &&(zmax == p.zmax)
	    && (wmin == p.wmin) &&(wmax == p.wmax));
  }

  R_tree_key_4 & operator=(const R_tree_key_4 &t)  {
    xmin=t.xmin;
    xmax=t.xmax;
    ymin=t.ymin;
    ymax=t.ymax;    
    zmin=t.zmin;
    zmax=t.zmax;
    wmin=t.wmin;
    wmax=t.wmax;
    return *this;
  }

  void unify( const R_tree_key_4 & p, const R_tree_key_4 & q){
    xmin = min(p.xmin, q.xmin);
    xmax = max(p.xmax, q.xmax);
    ymin = min(p.ymin, q.ymin);
    ymax = max(p.ymax, q.ymax);
    zmin = min(p.zmin, q.zmin);
    zmax = max(p.zmax, q.zmax);
    wmin = min(p.wmin, q.wmin);
    wmax = max(p.wmax, q.wmax);
  }

  //returns true if *this includes y
  bool include(const R_tree_key_4& y) const {
    if ((xmin > y.xmin) || (ymin >  y.ymin))
      return false;
    if ((xmax < y.xmax) || (ymax < y.ymax))
      return false;
    if ((zmax < y.zmax) || (zmax < y.zmax))
      return false;
    if ((wmax < y.wmax) || (wmax < y.wmax))
      return false;
    return true;
  }

  bool compare(const R_tree_key_4 & y) const{
    return intersect(y);
  }


  bool intersect( const R_tree_key_4& y) const {
    if ((xmax <= y.xmin) || (y.xmax <= xmin)) {
      return false;
    }
    if ((ymax <= y.ymin) || (y.ymax <= ymin)) {
      return false;
    }
    if ((zmax <= y.zmin) || (y.zmax <= zmin)) {
      return false;
    }
    if ((wmax <= y.wmin) || (y.wmax <= wmin)) {
      return false;
    }
    return true;
  } 


  double cost() const {
    return (xmax - xmin) * (ymax - ymin)* (zmax - zmin)* (wmax - wmin);
  }


  void intersection( const R_tree_key_4 & p, const R_tree_key_4 & q){
    if(p.intersect(q))
      {
	xmin = max(p.xmin, q.xmin);
	xmax = min(p.xmax, q.xmax);
	ymin = max(p.ymin, q.ymin);
	ymax = min(p.ymax, q.ymax);
	zmin = max(p.zmin, q.zmin);
	zmax = min(p.zmax, q.zmax);
	zmin = max(p.wmin, q.wmin);
	zmax = min(p.wmax, q.wmax);
      }
  }


  // compute the distance between the centers of the keys
  double center_dist( const R_tree_key_4& q) const {
    double x=xmin + 0.5*(xmax-xmin);
    double y=ymin + 0.5*(ymax-ymin);
    double z=zmin + 0.5*(zmax-zmin);
    double w=wmin + 0.5*(wmax-wmin);
    double qx=q.xmin + 0.5*(q.xmax-q.xmin);
    double qy=q.ymin + 0.5*(q.ymax-q.ymin);
    double qz=q.zmin + 0.5*(q.zmax-q.zmin);
    double qw=q.wmin + 0.5*(q.wmax-q.wmin);
    double dist = (x-qx)*(x-qx) + (y-qy)*(y-qy) + 
      (z-qz)*(z-qz) + (w-qw)*(w-qw); 
    return dist;
  }

  void read(char **s) {
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    xmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    xmax=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    ymin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    ymax=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    zmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    zmax=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    wmin=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    wmax=*((int *)from_int);

    delete[] from_int;
  }

  void write(char **s)
    {
      int sint=(int)sizeof(int);
      char *from_int=new char[sint];
      int i,r=0;
      memcpy(from_int,(char *)(&xmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i] = from_int[i];
      r += sint;

      memcpy(from_int,(char *)(&xmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&ymin),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&ymax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&zmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&zmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&wmin),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];
      r += sint;

      memcpy(from_int, (char *)(&wmax),sint);
      for (i=0; i<sint; i++)
	(*s)[i+r] = from_int[i];

      delete[] from_int;
    }


  class C_Compare_1{
  public:
    bool operator()(const R_tree_key_4 &p, const R_tree_key_4 &q)
    {
      if(p.xmin < q.xmin)
        return true;
      else 
	if(p.xmin==q.xmin)
	  if(p.xmax < q.xmax)
	    return true;
      return false;
    }
  };

  class C_Compare_2{
  public:
    bool operator()(const R_tree_key_4 &p, const R_tree_key_4 &q)
    {
      if(p.ymin < q.ymin)
        return true;
      else 
	if(p.ymin==q.ymin)
	  if(p.ymax < q.ymax)
	    return true;
      return false;
    }
  };

  class C_Compare_3{
  public:
    bool operator()(const R_tree_key_4 &p, const R_tree_key_4 &q)
    {
      if(p.ymin < q.ymin)
        return true;
      else 
	if(p.ymin==q.ymin)
	  if(p.ymax < q.ymax)
	    return true;
      return false;
    }
  };


  class C_Compare_4{
  public:
    bool operator()(const R_tree_key_4 &p, const R_tree_key_4 &q)
    {
      if(p.ymin < q.ymin)
        return true;
      else 
	if(p.ymin==q.ymin)
	  if(p.ymax < q.ymax)
	    return true;
      return false;
    }
  };
  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;
  typedef C_Compare_3 compare_3;
  typedef C_Compare_4 compare_4;



  void dump(int depth=0) const	{
    std::cerr << "(" << xmin << "/" << ymin << "/" << zmin << "/" << wmin 
         << "),(" << xmax << "/" << ymax << "/" << zmax << "/" << wmax << ")";
  }

};

std::ostream& operator << (std::ostream& s, R_tree_key_4 m) {
  s  << "(" << m.xmin << "/" << m.ymin <<  "/" << m.zmin << "/" << m.wmin
     << ")("
     << m.xmax << "/" << m.ymax <<  "/" << m.zmax << "/" << m.wmax << ")"
     << std::endl;
  return s;
}


//************************************************************
// Necessary classes for the split method of BKS    


template<class Container>
class X_Compare_1{
  typedef typename Container::Key Key;
  typedef typename Key::compare_1 comp_1;
  comp_1 ccc1;
public:
  bool operator()(const Container &p, const Container &q)
    {
      return ccc1(p.key,q.key);
    }
};

template<class Container>
class X_Compare_2{
  typedef typename Container::Key Key;
  typedef typename Key::compare_2 comp_2;
  comp_2 ccc2;
public:
  bool operator()(const Container &p, const Container &q)
    { 
      return ccc2(p.key,q.key);
    }
};

template<class Container>
class X_Compare_3{
  typedef typename Container::Key Key;
  typedef typename Key::compare_3 comp_3;
  comp_3 ccc3;
public:
  bool operator()(const Container &p, const Container &q)
    { 
      return ccc3(p.key,q.key);
    }
};

template<class Container>
class X_Compare_4{
  typedef typename Container::Key Key;
  typedef typename Key::compare_4 comp_4;
  comp_4 ccc4;
public:
  bool operator()(const Container &p, const Container &q)
    { 
      return ccc4(p.key,q.key);
    }
};


// return true if sort axis exists
template<class Container>
class sort_axis_key_1_dim
{
public:
  typedef typename Container::Key Key;
  typedef typename Key::compare_1 comp_1;

  typedef X_Compare_1<Container> compare_1;
  compare_1 cc1;
  bool operator()(int split_axis, Container *first, Container *last)
    {
      if(split_axis==0)
	{  
	  std::stable_sort(first,last,cc1);
	  return true;
	}
      else
	return false;

      return false;
    }
};




// return true if sort axis exists
template<class Container>
class sort_axis_key_2_dim
{
public:
  typedef typename Container::Key Key;
  typedef typename Key::compare_1 comp_1;
  typedef typename Key::compare_2 comp_2;

  typedef X_Compare_1<Container> compare_1;
  typedef X_Compare_2<Container> compare_2;
  compare_1 cc1;
  compare_2 cc2;
  bool operator()(int split_axis, Container *first, Container *last)
    {
      if(split_axis==0)
	{  
	  std::stable_sort(first,last,cc1);
	  return true;
	}
      else
	if(split_axis==1)
	  {
	    std::stable_sort(first,last,cc2);
	    return true;
	  }
	else
	  return false;

      return false;
    }
};


// return true if sort axis exists
template<class Container>
class sort_axis_key_3_dim
{
public:
  typedef typename Container::Key Key;
  typedef typename Key::compare_1 comp_1;
  typedef typename Key::compare_2 comp_2;
  typedef typename Key::compare_3 comp_3;

  typedef X_Compare_1<Container> compare_1;
  typedef X_Compare_2<Container> compare_2;
  typedef X_Compare_3<Container> compare_3;
  compare_1 cc1;
  compare_2 cc2;
  compare_3 cc3;
  bool operator()(int split_axis, Container *first, Container *last)
    {
      if(split_axis==0)
	{  
	  std::stable_sort(first,last,cc1);
	  return true;
	}
      else
	if(split_axis==1)
	  {
	    std::stable_sort(first,last,cc2);
	    return true;
	  }
	else
	  if(split_axis==2)
	    {
	      std::stable_sort(first,last,cc3);
	      return true;
	    }
	  else
	    return false;

      return false;
    }
};





// return true if sort axis exists
template<class Container>
class sort_axis_key_4_dim
{
public:
  typedef typename Container::Key Key;
  typedef typename Key::compare_1 comp_1;
  typedef typename Key::compare_2 comp_2;
  typedef typename Key::compare_3 comp_3;
  typedef typename Key::compare_4 comp_4;

  typedef X_Compare_1<Container> compare_1;
  typedef X_Compare_2<Container> compare_2;
  typedef X_Compare_3<Container> compare_3;
  typedef X_Compare_4<Container> compare_4;
  compare_1 cc1;
  compare_2 cc2;
  compare_3 cc3;
  compare_4 cc4;
  bool operator()(int split_axis, Container *first, Container *last)
    {
      if(split_axis==0)
	{  
	  std::stable_sort(first,last,cc1);
	  return true;
	}
      else
	if(split_axis==1)
	  {
	    std::stable_sort(first,last,cc2);
	    return true;
	  }
	else
	  if(split_axis==2)
	    {
	      std::stable_sort(first,last,cc3);
	      return true;
	    }
	  else
	    if(split_axis==3)
	      {
		std::stable_sort(first,last,cc4);
		return true;
	      }
	    else
	      return false;

      return false;
    }
};


CGAL_END_NAMESPACE
#endif




