// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/ch_map.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/Hash_map.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Hashing Map by chained hashing
// ============================================================================
#ifndef CGAL_GENERIC_CH_MAP_H
#define CGAL_GENERIC_CH_MAP_H

#if !(defined(__exportC) && defined(nil))
#define LEDANOTHERE
#endif

#ifdef LEDANOTHERE
#define __exportC
#define nil 0
#endif
#define GenPtr T

namespace leda {

template <typename T> class __exportC ch_map;
template <typename T> class __exportC ch_map_elem;

template <typename T>
class __exportC ch_map_elem 
{
  friend class ch_map<T>;

  unsigned long    k;
  GenPtr           i;
  ch_map_elem<T>*  succ;
};

//typedef ch_map_elem*  ch_map_item;

template <typename T>
class __exportC ch_map   
{
   const unsigned long NULLKEY; 
   const unsigned long NONNULLKEY;

   ch_map_elem<T> STOP;

   ch_map_elem<T>* table;
   ch_map_elem<T>* table_end;
   ch_map_elem<T>* free;
   int table_size;           
   int table_size_1;  

   ch_map_elem<T>* old_table;
   ch_map_elem<T>* old_table_end;
   ch_map_elem<T>* old_free;
   int old_table_size;           
   int old_table_size_1;  

   unsigned long old_index;

   // no virtual necessary            
   //virtual void clear_inf(GenPtr&)  const { }
   //virtual void copy_inf(GenPtr&)   const { }
   //virtual void init_inf(GenPtr&)   const { }
protected:
   GenPtr& xdef() { return STOP.i; }
   const GenPtr& cxdef() const { return STOP.i; }
private:
   void init_inf(GenPtr& x)   const { x = STOP.i; }

   
   ch_map_elem<T>*  HASH(unsigned long x)  const
   { return table + (x & table_size_1);  }
   
   void init_table(int t);
   void rehash();
   void del_old_table();

   inline void insert(unsigned long x, GenPtr y);

protected:
   typedef ch_map_elem<T>*  ch_map_item;
   typedef ch_map_item item;

   unsigned long index(ch_map_item it) const { return it->k; }
   GenPtr&       inf(ch_map_item it) const { return it->i; }

   ch_map(int n = 1); 
   ch_map(const ch_map<T>& D);
   ch_map& operator=(const ch_map<T>& D);
   

   void clear_entries();
   void clear();
   ~ch_map() { if (old_table) delete[] old_table; delete[] table; }

   GenPtr& access(ch_map_item p, unsigned long x);
   GenPtr& access(unsigned long x);
   ch_map_item lookup(unsigned long x) const;
   ch_map_item first_item() const;
   ch_map_item next_item(ch_map_item it) const; 

   void statistics() const;
};

template <typename T>
inline GenPtr& ch_map<T>::access(unsigned long x)
{ ch_map_item p = HASH(x);

  if (old_table) del_old_table();
  if ( p->k == x ) { 
     old_index = x; 
     return p->i;
  }
  else {
    if ( p->k == NULLKEY ) {
      p->k = x;
      init_inf(p->i);  // initializes p->i to xdef
      old_index = x;
      return p->i;
    } else 
      return access(p,x);
  }
}

template <typename T>
void ch_map<T>::init_table(int t)
{ 
  table_size = t;
  table_size_1 = t-1;
  table = new ch_map_elem<T>[t + t/2];
  free = table + t;
  table_end = table + t + t/2;      

  for (ch_map_item p = table; p < free; p++) 
  { p->succ = &STOP; 
    p->k = NULLKEY;
  }
  table->k = NONNULLKEY;
}


template <typename T>
inline void ch_map<T>::insert(unsigned long x, GenPtr y)
{ ch_map_item q = HASH(x);                                    
  if ( q->k == NULLKEY ) {      
    q->k = x;                                                  
    q->i = y; 
  } else { 
    free->k = x;                                                
    free->i = y;                                                
    free->succ = q->succ;                                       
    q->succ = free++; 
  }                                         
}

                                                                            
template <typename T>
void ch_map<T>::rehash()
{ 
  old_table = table;
  old_table_end = table_end;
  old_table_size = table_size;
  old_table_size_1 = table_size_1;
  old_free = free;

  ch_map_item old_table_mid = table + table_size;

  init_table(2*table_size);

  ch_map_item p;

  for(p = old_table + 1; p < old_table_mid; p++)
  { unsigned long x = p->k;
    if ( x != NULLKEY ) // list p is non-empty
    { ch_map_item q = HASH(x);  
      q->k = x;
      q->i = p->i;
    }
  }

  while (p < old_table_end)
  { unsigned long x = p->k;
    insert(x,p->i);
    p++;
  }
}


template <typename T>
void ch_map<T>::del_old_table()
{
  ch_map_item save_table = table;
  ch_map_item save_table_end = table_end;
  ch_map_item save_free = free;
  int save_table_size = table_size;
  int save_table_size_1 = table_size_1;

  table = old_table;
  table_end = old_table_end;
  table_size = old_table_size;
  table_size_1 = old_table_size_1;
  free = old_free;

  old_table = 0;

  GenPtr p = access(old_index);

  delete[] table;

  table = save_table;
  table_end = save_table_end;
  table_size = save_table_size;
  table_size_1 = save_table_size_1;
  free = save_free;
  access(old_index) = p;
}

template <typename T>
GenPtr& ch_map<T>::access(ch_map_item p, unsigned long x)
{
  STOP.k = x;
  ch_map_item q = p->succ; 
  while (q->k != x) q = q->succ;
  if (q != &STOP) 
  { old_index = x;
    return q->i;
  }

  // index x not present, insert it

  if (free == table_end)   // table full: rehash
  { rehash();
    p = HASH(x);
  }

  if (p->k == NULLKEY)
  { p->k = x;
    init_inf(p->i);  // initializes p->i to xdef
    return p->i;
  }

  q = free++;
  q->k = x;
  init_inf(q->i);    // initializes q->i to xdef
  q->succ = p->succ;
  p->succ = q;
  return q->i;
}


template <typename T>
ch_map<T>::ch_map(int n) : NULLKEY(0), NONNULLKEY(1), old_table(0)
{ 
  if (n < 512)
    init_table(512); 
  else {
    int ts = 1;
    while (ts < n) ts <<= 1;
    init_table(ts);
  }
}


template <typename T>
ch_map<T>::ch_map(const ch_map<T>& D) : 
  NULLKEY(0), NONNULLKEY(1), old_table(0)
{ 
  init_table(D.table_size);
  STOP.i = D.STOP.i; // xdef

  for(ch_map_item p = D.table + 1; p < D.free; p++) 
  { if (p->k != NULLKEY || p >= D.table + D.table_size)
    { insert(p->k,p->i);
      //D.copy_inf(p->i);  // see chapter Implementation
    }
  }
}

template <typename T>
ch_map<T>& ch_map<T>::operator=(const ch_map<T>& D)
{ 
  clear_entries();
  delete[] table;
  init_table(D.table_size);
  STOP.i = D.STOP.i; // xdef

  for(ch_map_item p = D.table + 1; p < D.free; p++) 
  { if (p->k != NULLKEY || p >= D.table + D.table_size)
    { insert(p->k,p->i);
      //copy_inf(p->i);    // see chapter Implementation
    }
  }
  return *this;
}

template <typename T>
void ch_map<T>::clear_entries() 
{ for(ch_map_item p = table + 1; p < free; p++)
    if (p->k != NULLKEY || p >= table + table_size) 
      p->i = T();  //clear_inf(p->i);  // see chapter Implementation
}

template <typename T>
void ch_map<T>::clear() 
{ clear_entries();
  delete[] table;
  init_table(512); 
}

template <typename T>
typename ch_map<T>::ch_map_item 
ch_map<T>::lookup(unsigned long x) const 
{ ch_map_item p = HASH(x);
  ((unsigned long &)STOP.k) = x;  // cast away const
  while (p->k != x) 
  { p = p->succ;
#if defined(LEDA_MULTI_THREAD)
    if (p == &STOP) break;
#endif
  }
  return (p == &STOP) ? nil : p;
}


template <typename T>
typename ch_map<T>::ch_map_item ch_map<T>::first_item() const
{ return next_item(table); }

template <typename T>
typename ch_map<T>::ch_map_item ch_map<T>::next_item(ch_map_item it) const 
{ if (it == nil) return nil;
  do it++; while (it < table + table_size && it->k == NULLKEY);
  return (it < free ? it : nil);
}

template <typename T>
void ch_map<T>::statistics() const
{ cout << "table_size: " << table_size <<"\n";
  int n = 0;
  for (ch_map_item p = table + 1; p < table + table_size; p++)
     if (p ->k != NULLKEY) n++;
  int used_in_overflow = free - (table + table_size );
  n += used_in_overflow;
  cout << "number of entries: " << n << "\n";
  cout << "fraction of entries in first position: " << 
          ((double) (n - used_in_overflow))/n <<"\n";
  cout << "fraction of empty lists: " << 
          ((double) (n - used_in_overflow))/table_size<<"\n";

}

} // namespace leda

#ifdef LEDANOTHERE
#undef __exportC
#undef nil
#endif
#undef GenPtr

#endif // CGAL_GENERIC_CH_MAP_H

