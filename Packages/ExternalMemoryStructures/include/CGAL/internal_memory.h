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
// file          : include/CGAL/internal_memory.h
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
// Provide access structure for data in a file.
// ======================================================================
#ifndef __INTERNAL_MEMORY_H__
#define __INTERNAL_MEMORY_H__

#include <CGAL/basic.h>
#include <ctime>
#include <fstream>
#include <cstring>
#include <cassert>

CGAL_BEGIN_NAMESPACE

//first written element in the data file. It stores implicitly 
//a list of deleted elements by link pfirst. The element pfirst 
//points to then stores a link to the next deleted element.
//plast always keeps the end of file position.
class Theader{
public:
  typedef ptrdiff_t difference_type;
  int  pfirst; //position of first deleted element
  int  plast;  //end of file position
  int number_of_elements;
  size_t size(void) { return 3*sizeof(int); }
  Theader() : pfirst(-1), plast(0), number_of_elements(0) {}
  
  void read(char *s){
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (s)[i];
    r += sint;
    pfirst=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (s)[i+r];
    r += sint;
    plast=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (s)[i+r];
    number_of_elements=*((int *)from_int);

    delete[] from_int;
  }
   void write (char *s) {
     int sint=(int)sizeof(int);
     char *from_int=new char[sint];
     int i,r=0;
     memcpy(from_int,(char *)(&pfirst),sint);
     for (i=0; i<sint; i++){
       (s)[i] = from_int[i];
     }
     r += sint;
     memcpy(from_int, (char *)(&plast),sint);
     for (i=0; i<sint; i++)
       (s)[i+r] = from_int[i];
     r += sint;
     memcpy(from_int, (char *)(&number_of_elements),sint);
     for (i=0; i<sint; i++)
       (s)[i+r] = from_int[i];
     r += sint;
     delete[] from_int;
  }
};

//Every stored item is associated with such a struct. 
//Here we store if the element is still valid. If it is
//not valid a pointer to the next non valid element is 
//stored (pnext).
class Telem {
public:
  typedef ptrdiff_t difference_type;
  bool deleted;
  time_t t;
  int pnext;
  size_t size(void) { return 1*sizeof(bool) + sizeof(int) + sizeof(time_t);}
  void read(char *s){
    int sint=(int)sizeof(int);
    int sbool=(int)sizeof(bool);
    int stime_t=(int)sizeof(time_t);
    
    char *from_bool=new char[sbool];
    char *from_int=new char[sint];
    char *from_time_t=new char[stime_t];
    int i,r=0;

    for (i=0; i<sbool; i++)
      from_bool[i] = (s)[i+r];
    r += sbool;
    deleted=*(reinterpret_cast<bool *>(from_bool));

    for (i=0; i<sint; i++)
      from_int[i] = (s)[i+r];
    r += sint;
    pnext=*((int *)from_int);

    for (i=0; i<stime_t; i++)
      from_time_t[i] = (s)[i+r];
    r += stime_t;
    t=(int)from_time_t;

    delete[] from_int;
    delete[] from_bool;
    delete[] from_time_t;
  }

  void write (char *s) {   
    int sint=(int)sizeof(int);
    int sbool=(int)sizeof(bool);
    int stime_t=(int)sizeof(time_t);
    
    char *from_bool=new char[sbool];
    char *from_int=new char[sint];
    char *from_time_t=new char[stime_t];
    int i,r=0;
    
    memcpy(from_bool, (char *)(&deleted), sbool);
    for (i=0; i<sbool; i++)
      (s)[i+r] = from_bool[i];
    r += sbool;
    
    memcpy(from_int,(char *)(&pnext),sint);
    for (i=0; i<sint; i++)
      (s)[i+r] = from_int[i];
    r += sint;
    
    memcpy(from_time_t,(char *)(&t),stime_t);
    for (i=0; i<stime_t; i++)
      (s)[i+r] = from_time_t[i];
    
    delete[] from_int;
    delete[] from_bool;
    delete[] from_time_t;
  }
};


//Elements of arbitrary types T that provide a read and a 
//write function are stored in a file.
class internal_memory{
public:
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  Theader header;
  Telem      elem;
   int IO_pagesize;
   int number_of_pages;
   int used_pages;
  char *raw_data;
  //cache(){ }

  //the file is opend and initialized if it is empty. Otherwise
  //the header of the file which contains the necessary information
  //is retrieved.
  internal_memory( int pagesize)
  {
    IO_pagesize=pagesize;
    number_of_pages=100;
    used_pages=0;
    int to_alloc = number_of_pages*(IO_pagesize+elem.size())+header.size();
    raw_data = new char[to_alloc];
    header.plast = 0;
    header.pfirst = -1;
    header.write(raw_data);
    header.read(raw_data);
  }

  //the file is cleaned up if it is empty.
  ~internal_memory() {
    close();
  } 
  void close(){
    delete [] raw_data;    
  }

  //same as the constructor.
  void open()
  {
    header.plast = 0;
    header.pfirst = -1;
    header.number_of_elements = 0;
    header.write(raw_data);
    header.read(raw_data);
  }
  
    
  //an new element is inserted. Its position is returned.
  difference_type insert(char** x) {
    difference_type pos = header.pfirst;
    Telem oldelem;

    elem.deleted = false;
    elem.pnext=-1;
    elem.t=time(&elem.t);
    header.number_of_elements++;
    if (pos==-1){
      //make raw_data larger!
      if(header.number_of_elements-1>=number_of_pages)
      {
	char *tmp_data = new char[2*number_of_pages*(IO_pagesize+elem.size())
	                          +header.size()];
	memcpy(tmp_data,raw_data,
	       number_of_pages*(IO_pagesize+elem.size())+header.size());
	delete [] raw_data;
	raw_data=tmp_data;
	number_of_pages=2*number_of_pages;
      }
      pos = header.plast;
      header.plast++;
      header.write(raw_data);
    }
    else {
      oldelem.read(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      header.pfirst=oldelem.pnext;
      header.write(raw_data);
    }
    elem.write(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
    memcpy(raw_data+(pos)*(elem.size()+IO_pagesize)+header.size()+elem.size(),
           *x,IO_pagesize);
    return pos;
  }

  //an new element is inserted. Its position is returned.
  difference_type get_pos() {
    difference_type pos = header.pfirst;
    Telem oldelem;
    header.number_of_elements++;
    if (pos==-1){
      if(header.number_of_elements-1>=number_of_pages) //make raw_data larger!
      {
	char *tmp_data = new char[2*number_of_pages*(IO_pagesize+elem.size())
	                          +header.size()];
	memcpy(tmp_data,raw_data,
	       number_of_pages*(IO_pagesize+elem.size())+header.size());
	delete [] raw_data;
	raw_data=tmp_data;
	number_of_pages=2*number_of_pages;
      }
      pos = header.plast;
      header.plast++;
      header.write(raw_data);
    }
    else {
      oldelem.read(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      header.pfirst=oldelem.pnext;
      header.write(raw_data);
    }
    return pos;
  }



  //an new element is inserted. Its position is returned.
  bool insert(difference_type pos, char** x) {
    elem.deleted = false;
    elem.pnext=-1;
    elem.t=time(&elem.t);
    elem.write(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
    memcpy(raw_data+(pos)*(elem.size()+IO_pagesize)+header.size()+elem.size(),
           *x,IO_pagesize);
    return true;
  }



  //returns the number of valid elements in the file.
  difference_type number_of_elements(){
    return header.number_of_elements;
  }


  // the element associated with position pos is returned.
  // x has to provide size space for writing     
  bool get(difference_type pos, char **x)
  {
    if(pos <0 || pos > header.plast)
      return false;
    else{
      elem.read(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      if(elem.deleted)
	return false;
      memcpy(*x,raw_data+(pos)*(elem.size()+IO_pagesize)+header.size()
                +elem.size(),IO_pagesize);
    }
    return true;
  }

  // the element at position pos is updated with x.
  bool update(difference_type pos, char **x){
    bool ret = true;
    if(pos <0 || pos > header.plast)
      return false;
    else{
      elem.read(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      if(elem.deleted)
	{
	  elem.deleted=true;
	  ret = false;
	}
      elem.write(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      memcpy(raw_data+(pos)*(elem.size()+IO_pagesize)+header.size()+elem.size()
             ,*x,IO_pagesize);
    }
    return ret;
  }

  //the element at position pos is declared deleted. Its free place
  //is stored in the header.
  bool erase(difference_type pos)
  {
    if(pos <0 || pos > header.plast)
      return false;
    else{
      elem.read(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      if(elem.deleted)
	return false;
      if(!elem.deleted)
	header.number_of_elements--;
      elem.deleted = true;
      elem.pnext=header.pfirst;
      header.pfirst=pos;
      elem.write(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      header.write(raw_data);
    }
    return true;
  }
  
  //returns true if the element at position pos is deleted.
  bool deleted(difference_type pos)
  {
    if(pos <0 || pos > header.plast)
      return true;
    else{
      elem.read(raw_data+(pos)*(elem.size() + IO_pagesize) + header.size());
      if(elem.deleted)
	return true;
    }
    return false;
  }
};

CGAL_END_NAMESPACE
#endif

