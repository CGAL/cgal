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
// file          : include/CGAL/cache.h
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
#ifndef __CACHE_H__
#define __CACHE_H__
#ifdef CGAL_D_DOGN_Cache 
  #define CGAL_DOGN_Cache(cmd) cmd 
#else
  #define CGAL_DOGN_Cache(cmd) 
#endif 

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
  int pfirst; //position of first deleted element
  int plast;  //end of file position
  int number_of_elements;
  //  size_t size(void) { return 3*sizeof(int)+3; }
  size_t size(void) { return 3*sizeof(int); }
  Theader() : pfirst(-1), plast(0), number_of_elements(0) {}

  void read(char **s){
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    pfirst=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    plast=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    number_of_elements=*((int *)from_int);
    CGAL_DOGN_Cache(
               std::cerr << "Theader number of elements:"
	                 << number_of_elements << std::endl;
	       )
    delete[] from_int;
  }
   void write (char **s) {
     int sint=(int)sizeof(int);
     char *from_int=new char[sint];
     int i,r=0;
     memcpy(from_int,(char *)(&pfirst),sint);
     for (i=0; i<sint; i++){
       (*s)[i] = from_int[i];
     }
     r += sint;
     memcpy(from_int, (char *)(&plast),sint);
     for (i=0; i<sint; i++)
       (*s)[i+r] = from_int[i];
     r += sint;
     memcpy(from_int, (char *)(&number_of_elements),sint);
     for (i=0; i<sint; i++)
       (*s)[i+r] = from_int[i];
     r += sint;
     delete[] from_int;
     CGAL_DOGN_Cache(
		std::cerr << "Theader number of elements:"
		<< number_of_elements << std::endl;
		)
  }
};

//Every stored item is associated with such a struct. 
//Here we store if the element is still valid. If it is
//not valid a pointer to the next non valid element is 
//stored (pnext).
class Telem {
public:
  typedef ptrdiff_t difference_type;
  bool dirty; //true if the element has to be written on disk
  bool deleted;
  int pnext;
  time_t t;

  //  size_t size(void)
  //  { return 2*sizeof(bool) + sizeof(int) + sizeof(time_t)+4; }
  size_t size(void) { return 2*sizeof(bool) + sizeof(int) + sizeof(time_t);}
  void read(char **s){
    int sint=(int)sizeof(int);
    int sbool=(int)sizeof(bool);
    int stime_t=(int)sizeof(time_t);
    
    char *from_bool=new char[sbool];
    char *from_int=new char[sint];
    char *from_time_t=new char[stime_t];
    int i,r=0;
    for (i=0; i<sbool; i++)
      from_bool[i] = (*s)[i];
    r += sbool;
    dirty=*(reinterpret_cast<bool *>(from_bool));

    for (i=0; i<sbool; i++)
      from_bool[i] = (*s)[i+r];
    r += sbool;
    deleted=*(reinterpret_cast<bool *>(from_bool));

    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i+r];
    r += sint;
    pnext=*((int *)from_int);

    for (i=0; i<stime_t; i++)
      from_time_t[i] = (*s)[i+r];
    r += stime_t;
    t=(int)from_time_t;

    delete[] from_int;
    delete[] from_bool;
    delete[] from_time_t;
  }

  void write (char **s) {   
    int sint=(int)sizeof(int);
    int sbool=(int)sizeof(bool);
    int stime_t=(int)sizeof(time_t);
    
    char *from_bool=new char[sbool];
    char *from_int=new char[sint];
    char *from_time_t=new char[stime_t];
    int i,r=0;
    memcpy(from_bool,(char *)(&dirty),sbool);
    for (i=0; i<sbool; i++)
      (*s)[i] = from_bool[i];
    r += sbool;
    
    memcpy(from_bool, (char *)(&deleted), sbool);
    for (i=0; i<sbool; i++)
      (*s)[i+r] = from_bool[i];
    r += sbool;
    
    memcpy(from_int,(char *)(&pnext),sint);
    for (i=0; i<sint; i++)
      (*s)[i+r] = from_int[i];
    r += sint;
    
    memcpy(from_time_t,(char *)(&t),stime_t);
    for (i=0; i<stime_t; i++)
      (*s)[i+r] = from_time_t[i];
    
    delete[] from_int;
    delete[] from_bool;
    delete[] from_time_t;
  }

};


struct Buffer_item{
  typedef ptrdiff_t difference_type;
  typedef Buffer_item B_i;
  Telem elem;
  char *raw_data;
  difference_type diskpos;
  B_i *next;
  B_i *prev;
};
  


//Elements of arbitrary types T that provide a read and a 
//write function are stored in a file.
template <const  int PagesInMemory>
class cache{
public:
  std::fstream XX;
  int IO_pagesize;
  const static  int IO_pages=PagesInMemory;
  typedef Buffer_item  BufferItem;
  char **header_string;
  char **elem_string;
  int header_size;
  int elem_size;
  BufferItem *buffer_first; //pointer to the first item in the buffer
  BufferItem *buffer_last;  //pointer to the last item in the buffer
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  Theader header;
  Telem      elem;
  cache( int Pagesize){
    header_string = new (char *);
    *header_string= new char[header.size()];
    header_size=header.size();
    elem_string = new (char *);
    *elem_string= new char[elem.size()];
    elem_size=elem.size();
    IO_pagesize=Pagesize;
    //initializiation of the buffer as a doubly connected list.
    //the complete space for all list items is allocated.

    if(IO_pages >0)
      {
	BufferItem * bi1,*bi2;
	 int i;
	buffer_first= bi1 = new BufferItem;
	bi1->raw_data=new char[IO_pagesize];
	bi1->prev=0;
	bi1->elem.dirty=false;
	bi1->elem.deleted=true;
	bi1->diskpos=-1;
	for(i=1;i<IO_pages;i++)
	  {
	    bi2=new BufferItem;
	    bi1->next=bi2;
	    bi2->raw_data=new char[IO_pagesize];
	    bi2->prev=bi1;
	    bi1=bi2;
	    bi1->elem.dirty=false;
	    bi1->elem.deleted=true;
	    bi1->diskpos=-1;
	  }
	buffer_last=bi1;
	bi1->next=0;
      }
  }

  //the file is opend and initialized if it is empty. Otherwise
  //the header of the file which contains the necessary information
  //is retrieved.
  cache( int Pagesize, char* name, 
	 int m = (std::fstream::binary|std::fstream::in | std::fstream::out))
  { 
    header_string = new (char *);
    *header_string= new char[header.size()];
    header_size=header.size();
    elem_string = new (char *);
    *elem_string= new char[elem.size()];
    elem_size=elem.size();

    //fake the SGI CC - compiler: it does not want to open empty 
    //files for reading 
    std::fstream fake_open;
    fake_open.open(name,
                   std::fstream::app|std::fstream::binary|std::fstream::out);
    if(!fake_open.is_open())
      std::cerr << "fake open is not open";
    else
      std::cerr << "header fake open is open";
    fake_open.close();

    std::cerr << "cache open called with " << name << std::endl;
    IO_pagesize=Pagesize;
    XX.open(name, m);  
    if(!XX.is_open())
      {
	std::cerr << name << "did not open -- this should no more happen";
      }
    XX.seekp(0, std::fstream::end);
    int teller;
    teller = XX.tellp();
    if (teller <= 0) {
      header.plast = 0;
      header.pfirst = -1;
      XX.seekp(0);
      header.write(header_string);
      XX.write(*header_string,header_size);
      CGAL_DOGN_Cache(
		 if(!XX)
		 std::cerr << "\n\n 1 Failed for fstream \n\n";
		 else
		 std::cerr << "\n\n 1 Did not fail\n\n";
		 )
    }
    XX.seekg(0);
    XX.read(*header_string,header_size);
    header.read(header_string);
    //initializiation of the buffer as a doubly connected list.
    //the complete space for all list items is allocated.
    if(IO_pages >0)
      {
	BufferItem * bi1, *bi2;
        int i;
	buffer_first= bi1 = new BufferItem; 
	bi1->raw_data=new char[IO_pagesize];
	bi1->prev=0;
	bi1->elem.dirty=false;
	bi1->elem.deleted=true;
	bi1->diskpos=-1;
	for(i=1;i<IO_pages;i++)
	  {
	    bi2=new BufferItem;
	    bi1->next=bi2;
	    bi2->raw_data=new char[IO_pagesize];	
	    bi2->prev=bi1;
	    bi1=bi2;
	    bi1->elem.dirty=false;
	    bi1->elem.deleted=true;
	    bi1->diskpos=-1;
	  }
	buffer_last=bi1;
	bi1->next=0;
	CGAL_DOGN_Cache(
		   test_buffer_items();
		   )
      }
  }

  //the file is cleaned up if it is empty.
  ~cache() {
    close();
  } 
  
  void close()
  {
    if(header.number_of_elements==0)
      {
	header.pfirst=-1;
	header.plast=0;
	XX.seekp(0);
	XX.clear();
	XX.seekp(0);
	header.write(header_string);
	XX.write(*header_string,header_size);
	CGAL_DOGN_Cache(
		   if(!XX)
		   std::cerr << "\n\n 2 Failed for fstream \n\n";
		   else
		   std::cerr << "\n\n 2 Did not fail\n\n";
		   )
      }
    else
      {
	XX.seekp(0);
	header.write(header_string);
	XX.write(*header_string,header_size);
	CGAL_DOGN_Cache(
		   if(!XX)
		   std::cerr << "\n\n 3 Failed for fstream \n\n";
		   else
		   std::cerr << "\n\n 3 Did not fail\n\n";
		   )
	CGAL_DOGN_Cache(
		   XX.seekg(0);
		   XX.read(*header_string,header_size);
		   header.read(header_string);
		   )
      }
    if(IO_pages>0)
      {
	BufferItem *bi=buffer_first;
	BufferItem *bi2;
	while(bi!=buffer_last)
	  {
	    if(bi->elem.dirty) //write to disk
	      {
		bi->elem.dirty=false;
		XX.seekp((bi->diskpos)*(elem.size() + IO_pagesize)
		                                    + header.size());
		bi->elem.write(elem_string);
		XX.write(*elem_string,elem_size);
		XX.seekp((bi->diskpos)*(elem.size() + IO_pagesize)
		                               + header.size() + elem.size());
		XX.write(bi->raw_data,IO_pagesize);
		CGAL_DOGN_Cache(
			   if(!XX)
			   std::cerr << "\n\n 4 Failed for fstream \n\n";
			   else
			   std::cerr << "\n\n 4 Did not fail\n\n";
			   )
		  
		  }
	    bi2=bi;
	    bi=bi->next;
	    delete [] bi2->raw_data;
	    delete bi2;	
	  }
	if(bi->elem.dirty) //write to disk
	  {
	    bi->elem.dirty=false;
	    XX.seekp((bi->diskpos)*(elem.size() + IO_pagesize)
	                                        + header.size());
	    bi->elem.write(elem_string);
	    XX.write(*elem_string,elem_size);
	    XX.seekp((bi->diskpos)*(elem.size() + IO_pagesize)
	                                        + elem.size() + header.size());
	    XX.write(bi->raw_data,IO_pagesize);
	    CGAL_DOGN_Cache(
		       if(!XX)
		       std::cerr << "\n\n 5 Failed for fstream \n\n";
		       else
		       std::cerr << "\n\n 5 Did not fail\n\n";
		       )

	  }
	delete [] bi->raw_data;
	delete bi;
      }
    XX.close();
  } 




  //same as the constructor.
   void open(char* name,
             int m =(std::fstream::binary|std::fstream::in|std::fstream::out)){
    //close the opened file first.
    std::cerr << "cache open called with " << name << std::endl;
	
    if(XX.is_open())
      {
	cerr << "Cache file was already open. then it does not work???";
	close();
      }
    //fake the SGI CC - compiler: it does not want to open empty 
    //files for reading 
    std::fstream fake_open;
    fake_open.open(name,
                   std::fstream::app|std::fstream::binary|std::fstream::out);
    CGAL_DOGN_Cache(
	       if(!fake_open.is_open())
	       std::cerr << "fake open is not open";
	       else
	       std::cerr << "header fake open is open";
	       )
    fake_open.close();

    XX.open(name, m);  
    if(!XX.is_open()){
      std::cerr << name << "can not be opened";
    }
    XX.seekp(0, std::fstream::end);
    int teller;
    teller = XX.tellp();
    if (teller <= 0) {
      header.plast = 0;
      header.pfirst = -1;
      header.number_of_elements = 0;
      XX.seekp(0);
      header.write(header_string);
      XX.write(*header_string,header_size);
    }
    XX.seekg(0);
    XX.read(*header_string,header_size);
    header.read(header_string);
  }

  void test_buffer_items()
  {
    if(IO_pages >0)
      {
	BufferItem *bi=buffer_first;
	std::cerr << "BUFFERITEMS::::::::::::::::\n";
	while(bi!=0)
	  {
	    std::cerr << "bi:" << bi->diskpos
	              << " prev :" << bi->prev
		      << "next :" << bi->next << std::endl;
	    bi=bi->next;
	  }
	std::cerr << "blast:" << buffer_last->diskpos
	          << " prev :" << buffer_last->prev
		  << "next :" << buffer_last->next << std::endl;
	std::cerr << "BUFFERITEMS::::::END::::::::::\n";
      }
  }

  BufferItem *find_item(difference_type pos)
  {  
    BufferItem *bi=buffer_first;
    while(bi->diskpos!=pos && bi!=buffer_last)
      bi=bi->next;
    if(bi==buffer_last && bi->diskpos!=pos)
      return 0;
    else
      return bi;
  }
    
  //an new element is inserted. Its position is returned. 
  difference_type get_pos() {
    difference_type pos = header.pfirst;
    Telem oldelem;
    header.number_of_elements++;
    if (pos==-1){
      pos = header.plast;
      header.plast++;
      CGAL_DOGN_Cache(
		 if(!XX)
		 std::cerr << "\n\n 6a Failed for fstream \n\n";
		 else
		 std::cerr << "\n\n 6a Did not fail\n\n";
		 )
      if(IO_pages<=0)
	{
	  XX.seekp(0);
	  header.write(header_string);
	  XX.write(*header_string,header_size);
	  CGAL_DOGN_Cache(
		     if(!XX)
		     std::cerr << "\n\n 6b Failed for fstream \n\n";
		     else
		     std::cerr << "\n\n 6b Did not fail\n\n";
		     )
	}
    }
    else //write element on an old place -> 
      //find out what the pnext of that place ist
      { 
	if(IO_pages >0)
	  {
	    BufferItem *bi;
	    if((bi=find_item(pos))!=0){
	      header.pfirst=bi->elem.pnext;
	    }
	    else
	      {
		XX.seekg( (pos)*(elem.size() + IO_pagesize) + header.size());
		XX.read(*elem_string,elem_size);
		oldelem.read(elem_string);
		header.pfirst=oldelem.pnext;
	      } 
	  }
	else
	  {
	    XX.seekg( (pos)*(elem.size() + IO_pagesize) + header.size());
	    XX.read(*elem_string,elem_size);
	    oldelem.read(elem_string);
	    header.pfirst=oldelem.pnext;
	  }
	XX.seekp(0);
	header.write(header_string);
	XX.write(*header_string,header_size);
	CGAL_DOGN_Cache(
		   if(!XX)
		   std::cerr << "\n\n 7 Failed for fstream \n\n";
		   else
		   std::cerr << "\n\n 7 Did not fail\n\n";
		   )
      }
    return pos;
  }
  
  bool insert(difference_type pos, char** x) {
    if(IO_pages>0)
      {
	// write buffer_last on disk
	if(buffer_last->elem.dirty)
	  {
	    CGAL_DOGN_Cache(
		    std::cerr << "Write the last element of buffer to disk\n";
		    std::cerr << "Write elem at pos: "
		              << (buffer_last->diskpos)*(elem.size()
			         + IO_pagesize) + header.size()
			      << std::endl;

		       )
	    buffer_last->elem.dirty=false;
	    XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
	                                                 + header.size());
	    buffer_last->elem.write(elem_string);
	    XX.write(*elem_string,elem_size);
	    XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
	                                   + elem.size() + header.size());
	    XX.write(buffer_last->raw_data,IO_pagesize);
	  }
	//write element into buffer_last
	buffer_last->elem.dirty=true;
	buffer_last->elem.deleted=false;
	buffer_last->elem.t=elem.t;
	buffer_last->elem.pnext=elem.pnext;
	buffer_last->elem.dirty=true;
	memcpy(buffer_last->raw_data,*x,IO_pagesize);
    
	//change pointers such that the element is the first in buffer
	buffer_last->prev->next=0;
	buffer_last->next=buffer_first;
	buffer_first->prev=buffer_last;
	buffer_first=buffer_last;
	buffer_last=buffer_last->prev;
	buffer_first->prev=0;
	buffer_first->diskpos=pos;
	CGAL_DOGN_Cache(
		   test_buffer_items();
		   )
      }
    else
      {
	elem.deleted = false;
	elem.dirty=false;
	elem.pnext=-1;
	elem.t=time(&elem.t);
	XX.seekp((pos)*(elem.size() + IO_pagesize) + header.size());
	elem.write(elem_string);
	XX.write(*elem_string,elem_size);
	XX.seekp((pos)*(elem.size() + IO_pagesize) + elem.size()
	                                           + header.size());
	XX.write(*x,IO_pagesize);
	CGAL_DOGN_Cache(
		   if(!XX)
		   std::cerr << "\n\n 9 Failed for fstream \n\n";
		   else
		   std::cerr << "\n\n 9 Did not fail\n\n";
		   )
      }
    return true;
  }

  //returns the number of valid elements in the file.
  difference_type number_of_elements(){
    return header.number_of_elements;
  }


  // the element associated with position pos is returned.
  // x has to provide IO_pagesize space for writing     
  bool get(difference_type pos, char **x)
  {
    if(pos <0 || pos > header.plast)
      return false;
    else{
      if(IO_pages>0)
	{
	  BufferItem *bi;
	  if((bi=find_item(pos))!=0){ //element in internal memory
	    memcpy(*x,bi->raw_data,IO_pagesize);
	    if(bi!=buffer_first)
	      {
		if(bi!=buffer_last)
		  bi->next->prev=bi->prev;
		else
		  buffer_last=bi->prev;
		bi->prev->next=bi->next;
		bi->next=buffer_first;
		bi->prev=0;
		bi->next->prev=bi;
		buffer_first=bi;
	      }
	  }
	  else
	    {
	      CGAL_DOGN_Cache(
			 std::cerr << "Element in external memory\n";
			 )
	      //write buffer_last to disk and read the new element in
	      if(buffer_last->elem.dirty)
		{
		  buffer_last->elem.dirty=false;
		  XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
		                                 + header.size());
		  CGAL_DOGN_Cache(
			     std::cerr << "Write elem at pos: "
			               << (buffer_last->diskpos)*(elem.size()
				          + IO_pagesize) + header.size()
				       << std::endl;
			     std::cerr << "Write data at pos: "
			               << (buffer_last->diskpos)*(elem.size()
				          + IO_pagesize) + header.size()
					  + elem.size()
			               << std::endl;
			     )
		  buffer_last->elem.write(elem_string);
		  XX.write(*elem_string,elem_size);
		  XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
		           + header.size() + elem.size());
		  XX.write(buffer_last->raw_data,IO_pagesize);
		}
	      XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	      XX.read(*elem_string,elem_size);
	      buffer_last->elem.read(elem_string);

	      if(buffer_last->elem.deleted)
		return false;
	      XX.seekg((pos)*(elem.size() + IO_pagesize)
	               + header.size()+ elem.size());
	      XX.read(*x,IO_pagesize);
	      memcpy(buffer_last->raw_data,*x,IO_pagesize);
	      buffer_last->diskpos=pos;
	      buffer_last->prev->next=0;
	      buffer_last->next=buffer_first;
	      buffer_first->prev=buffer_last;
	      buffer_first=buffer_last;
	      buffer_last=buffer_last->prev;
	      buffer_first->prev=0;
	    }
	}
      else
	{
	  XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	  CGAL_DOGN_Data(
		    std::cerr << "*********Get elem at pos: "
		              << (pos)*(elem.size() + IO_pagesize)
			         + header.size()
			      << std::endl;
		    std::cerr << "******Get data at pos: "
		              << (pos)*(elem.size() + IO_pagesize)
			         + header.size()+ elem.size()
			      << std::endl;
		    )
	  XX.read(*elem_string,elem_size);
	  elem.read(elem_string);
	  if(elem.deleted)
	    return false;
	  XX.seekg((pos)*(elem.size()+IO_pagesize)+header.size()+elem.size());
	  XX.read(*x,IO_pagesize);
	}
      CGAL_DOGN_Cache(
		 test_buffer_items();
		 )
      return true;
    }
  }

  // the element at position pos is updated with x.
  bool update(difference_type pos, char **x){
    bool ret = true;
    if(pos <0 || pos > header.plast)
      return false;
    else{
      if(IO_pages>0)
	{
	  BufferItem *bi;
	  if((bi=find_item(pos))!=0){
	    memcpy(bi->raw_data,*x,IO_pagesize);
	    bi->elem.dirty=true;
	    bi->elem.deleted=false;
	    if(bi!=buffer_first) //internal update
	      {
		if(bi!=buffer_last)
		  bi->next->prev=bi->prev;
		else
		  buffer_last=bi->prev;
		bi->prev->next=bi->next;
		bi->next=buffer_first;
		bi->prev=0;
		bi->next->prev=bi;
		buffer_first=bi;
	      }
	  }
	  else  //external update
	    {
	      // write buffer_last to disk, read the new element in
	      // and update it
	      if(buffer_last->elem.dirty)
		{
		  buffer_last->elem.dirty=false;
		  XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
		           + header.size());
		  buffer_last->elem.write(elem_string);
		  XX.write(*elem_string,elem_size);
		  XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
		           + header.size()+ elem.size());
		  XX.write(buffer_last->raw_data,IO_pagesize);
		}
	      XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	      XX.read(*elem_string,elem_size);
	      buffer_last->elem.read(elem_string);
	      if(buffer_last->elem.deleted)
	      {
		CGAL_DOGN_Cache(
			   std::cerr << "Warning: element was declared deleted"
			             << std::endl;
			   )
		ret = false;
	      }
	      memcpy(buffer_last->raw_data,*x,IO_pagesize);
	      buffer_last->elem.dirty=true;
	      buffer_last->diskpos=pos;
	      buffer_last->prev->next=0;
	      buffer_last->next=buffer_first;
	      buffer_first->prev=buffer_last;
	      buffer_first=buffer_last;
	      buffer_last=buffer_last->prev;
	      buffer_first->prev=0;
	    }
	}     
      else
	{
	  XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	  XX.read(*elem_string,elem_size);
	  elem.read(elem_string);

	  if(elem.deleted)
	  {
	    CGAL_DOGN_Cache(
		       std::cerr << "Warning: element was declared deleted"
		                 << std::endl;
		       )
	    ret = false;
	  }
	  XX.seekp((pos)*(elem.size() + IO_pagesize) + header.size());
	  CGAL_DOGN_Cache(
		     std::cerr << "******Write elem at pos: "
		               << (pos)*(elem.size() + IO_pagesize)
			          + header.size()
			       << std::endl;
		     std::cerr << "******Write data at pos: "
		               << (pos)*(elem.size() + IO_pagesize)
			          + header.size()+ elem.size()
			       << std::endl;
		     )
	  elem.write(elem_string);
	  XX.write(*elem_string,elem_size); 
	  XX.seekp((pos)*(elem.size() + IO_pagesize) + header.size()
	           + elem.size());
	  XX.write(*x,IO_pagesize);
	}
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
      if(IO_pages >0)
	{
	  BufferItem *bi;
	  if((bi=find_item(pos))!=0) //element is in internal memory
	    { 
	      if(!bi->elem.deleted)
		header.number_of_elements--;
	      bi->elem.deleted=true;
	      bi->elem.dirty=true;
	      bi->elem.pnext=header.pfirst;
	      header.pfirst=pos;
	      XX.seekp(0);
	      header.write(header_string);
	      XX.write(*header_string,header_size);

	      // position bi on last position
	      if(bi!=buffer_last){
		if(bi!=buffer_first)
		  bi->prev->next=bi->next;
		else
		  buffer_first=bi->next;
		bi->next->prev=bi->prev;
		bi->prev=buffer_last;
		bi->next=0;
		buffer_last->next=bi;
		buffer_last=bi;
	      }
	    }
	  else
	    {
	      CGAL_DOGN_Cache(
			 std::cerr << "Element in external memory\n";
			 )
	      XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	      XX.read(*elem_string,elem_size);
	      elem.read(elem_string);

	      if(elem.deleted)
		return false;
	      
	      elem.deleted = true;
	      elem.pnext=header.pfirst;
	      header.pfirst=pos;
	      header.number_of_elements--;
	      XX.seekp((pos)*(elem.size() + IO_pagesize) + header.size());
	      elem.write(elem_string);
	      XX.write(*elem_string,elem_size);
	      XX.seekp(0);
	      header.write(header_string);
	      XX.write(*header_string,header_size);
	    }
	}
      else
	{  
	  XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	   XX.read(*elem_string,elem_size);
	   elem.read(elem_string);
		
	   if(elem.deleted)
	     return false;  
	   elem.deleted = true;  
	   elem.pnext=header.pfirst; 
	   header.pfirst=pos;
	   header.number_of_elements--;
	   XX.seekp((pos)*(elem.size() + IO_pagesize) + header.size());
	   elem.write(elem_string);
	   XX.write(*elem_string,elem_size);
	   XX.seekp(0);
	   header.write(header_string);
	   XX.write(*header_string,header_size);
	}
    }
    CGAL_DOGN_Cache(
	       test_buffer_items();
	       )
    return true;
  }
  
  //returns true if the element at position pos is deleted.
  bool deleted(difference_type pos)
  {
    if(pos <0 || pos > header.plast)
      return true;
    else{
      if(IO_pages>0)
	{ 
	  BufferItem *bi;
	  if((bi=find_item(pos))!=0) //element is in internal memory
	    if(bi->elem.deleted)
	      return true;
	    else
	      return false;
	  else{
	    //write buffer_last to disk, read the new element in and update it
	    if(buffer_last->elem.dirty)
	      {
		buffer_last->elem.dirty=false;
		XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
		         + header.size());
		buffer_last->elem.write(elem_string);
		XX.write(*elem_string,elem_size);
		XX.seekp((buffer_last->diskpos)*(elem.size() + IO_pagesize)
		         + header.size()+ elem.size());
		XX.write(buffer_last->raw_data,IO_pagesize);
	      }
	    XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	    XX.read(*elem_string,elem_size);
	    buffer_last->elem.read(elem_string);
	    if(buffer_last->elem.deleted)
	      {
		CGAL_DOGN_Cache(
			   std::cerr << "Warning: element was declared deleted"
			             << std::endl;
			   )
		ret = false;
	      }
	    buffer_last->diskpos=pos;
	    buffer_last->prev->next=0;
	    buffer_last->next=buffer_first;
	    buffer_first->prev=buffer_last;
	    buffer_first=buffer_last;
	    buffer_last=buffer_last->prev;
	    buffer_first->prev=0;  
	    if(buffer_last->elem.deleted)
	      return true;
	  }
	}
      else
	{
	  XX.seekg((pos)*(elem.size() + IO_pagesize) + header.size());
	  XX.read(*elem_string,elem_size);
	  elem.read(elem_string);
	  if(elem.deleted)
	    return true;
	}
    }
    CGAL_DOGN_Cache( 
	       test_buffer_items();
	       )
    return false;
  }
};

CGAL_END_NAMESPACE
#endif














































