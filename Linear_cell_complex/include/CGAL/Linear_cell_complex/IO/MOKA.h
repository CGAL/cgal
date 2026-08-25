// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef IMPORT_MOKA_H
#define IMPORT_MOKA_H

#include <vector>
#include <stack>
#include <fstream>
#include <string>
#include <iostream>
#include <cassert>
#include <unordered_map>

namespace CGAL::IO
{
template<typename LCC>
struct GDart
{
  unsigned int alpha[4];
  typename LCC::Dart_descriptor dh;
  typename LCC::Vertex_attribute_descriptor vh;

  GDart() : dh(LCC::null_descriptor), vh(LCC::null_descriptor)
  {}

  GDart(const GDart& adart) : dh(adart.dh),
    vh(adart.vh)
  {
    for (unsigned int i=0; i<4; ++i)
    { alpha[i]=adart.alpha[i]; }
  }
};

template<typename LCC>
bool read_MOKA(LCC& lcc, const char* filename)
{
  typedef typename LCC::Point Point;

  std::ifstream ifile(filename);
  if (!ifile)
  {
    std::cout<<"Error opening file "<<filename<<"."<<std::endl;
    return false;
  }

  std::string line;
  std::getline(ifile, line);

  if ( line == "Moka file [binary]" )
  {
    std::cout<<"Binary file not (yet) considered.\n";
    return false;
  }
  else if ( line != "Moka file [ascii]" )
  {
    std::cout<<"File "<<filename<<" is not a moka file.\n";
    std::cout<< line;
    return false;
  }

  // To skip the masks mark (TODO read the marks ?)
  std::getline(ifile, line);

  std::vector<GDart<LCC>> gdarts;
  unsigned int nbLoaded = 0;
  unsigned int number;
  double x,y,z;

  // First load all the gdarts, and create vertex attributes
  while(ifile)
  {
    GDart<LCC> agdart;
    ifile>>agdart.alpha[0]>>agdart.alpha[1]
        >>agdart.alpha[2]>>agdart.alpha[3]; // the 4 alpha
    ifile>>number>>number>>number>>number; // to skip the 4*8 marks
    if ( agdart.alpha[0]==nbLoaded )
    {
      std::cout<<"Impossible to load a moka file with 0-free darts.\n";
      return false;
    }
    if ( ifile )
    {
      ifile>>number; // bool to know if dart has a vertex of not.
      if (number)
      {
        ifile>>x>>y>>z;
        agdart.vh = lcc.create_vertex_attribute(Point(x, y, z));
      }

      gdarts.push_back(agdart);
      ++nbLoaded;
    }
  }
  ifile.close();

  // Second orient the gmap, and create oriented darts.
  std::stack<unsigned int> totreat;
  for (unsigned int startingdart = 0; startingdart<nbLoaded; ++startingdart)
  {
    bool orient=(gdarts[startingdart].dh==lcc.null_descriptor);
    for (unsigned int dim=0; orient && dim<4; ++dim)
      if (gdarts[gdarts[startingdart].alpha[dim]].dh!=lcc.null_descriptor) orient=false;

    if ( orient )
    {
      totreat.push(startingdart);
      gdarts[startingdart].dh=lcc.create_dart();

      while ( !totreat.empty() )
      {
        unsigned int i=totreat.top();
        totreat.pop();

        assert(gdarts[i].dh!=lcc.null_descriptor);

        for (unsigned int dim=1; dim<4; ++dim)
        {
          if (gdarts[i].alpha[dim]!=i &&
              gdarts[gdarts[i].alpha[dim]].vh!=lcc.null_descriptor)
          {
            gdarts[i].vh = gdarts[gdarts[i].alpha[dim]].vh;
            gdarts[gdarts[i].alpha[dim]].vh=lcc.null_descriptor;
          }

          unsigned int alpha0 = gdarts[i].alpha[0];
          assert( alpha0!=i );

          if (gdarts[alpha0].alpha[dim]!=alpha0)
          {
            if ( gdarts[gdarts[alpha0].alpha[dim]].dh==lcc.null_descriptor )
            {
              totreat.push(gdarts[alpha0].alpha[dim]);
              gdarts[gdarts[alpha0].alpha[dim]].dh = lcc.create_dart();
              lcc.basic_link_beta(gdarts[i].dh,
                                  gdarts[gdarts[alpha0].alpha[dim]].dh,
                  dim);
            }
            else if (lcc.is_free(gdarts[i].dh, dim))
            {
              lcc.basic_link_beta(gdarts[i].dh,
                                  gdarts[gdarts[alpha0].alpha[dim]].dh,
                  dim);
            }
          }
        }
      }
    }
  }

  // Test that the gmap was orientable.
  bool orientable = true;
  for (unsigned int i = 0; i<nbLoaded; ++i)
  {
    if (gdarts[i].dh!=lcc.null_descriptor)
    {
      for (unsigned int dim=0; dim<4; ++dim)
      {
        if (orientable &&
            gdarts[i].alpha[dim]!=i &&
            gdarts[gdarts[i].alpha[dim]].dh!=lcc.null_descriptor)
        {
          std::cout<<"Pb, the gmap is NOT orientable."<<std::endl;
          orientable=false;
          // lcc.clear();
        }
      }

      /* if ( lcc.template attribute<3>(gdarts[i].dh) == lcc.null_descriptor )
      {
        lcc.template set_attribute<3>(gdarts[i].dh, lcc.template create_attribute<3>());
      } */
    }
    if (gdarts[i].vh!=lcc.null_descriptor)
    {
      lcc.set_vertex_attribute(gdarts[i].dh, gdarts[i].vh);
    }
  }

  return true;
}

template<typename LCC>
bool write_MOKA(LCC& lcc, const char* filename)
{
  std::ofstream os(filename);
  if (!os)
  {
    std::cout<<"Error opening file "<<filename<<"."<<std::endl;
    return false;
  }

  os<<"Moka file [ascii]"<<std::endl;
  os<<"0 0 0 0 0 0 0 0"<<std::endl; // For now, marks are not saved

  // Number darts
  std::unordered_map<typename LCC::Dart_const_descriptor,
                     typename LCC::size_type> ids;
  typename LCC::Dart_range::const_iterator it(lcc.darts().begin());
  typename LCC::size_type num=0;
  for(; it!=lcc.darts().end(); num+=2, ++it)
  { ids[it]=num; }

  for(num=0, it=lcc.darts().begin(); it!=lcc.darts().end(); num+=2, ++it)
  {
    assert(!lcc.is_free(it, 0));
    assert(!lcc.is_free(it, 1));

    // First g-dart
    os<<num+1<<" "<<ids[lcc.beta(it, 0)]+1<<" ";
    if(lcc.is_free(it, 2)) { os<<num<<" "; }
    else { os<<ids[lcc.beta(it, 2)]+1<<" "; }
    if(lcc.is_free(it, 3)) { os<<num<<" "; }
    else { os<<ids[lcc.beta(it, 3)]+1<<" "; }
    os<<"0 0 0 0 "; // 4 bytes => 32 Boolean marks (unsaved for now)

    if(lcc.template dart<0>(it)==it)
    { os<<"1 "<<lcc.point(it)<<std::endl; } // Point
    else
    { os<<"0"<<std::endl; } // No point

    // Second g-dart
    os<<num<<" "<<ids[lcc.beta(it, 1)]<<" ";
    if(lcc.is_free(it, 2)) { os<<num+1<<" "; }
    else { os<<ids[lcc.beta(it, 2)]<<" "; }
    if(lcc.is_free(it, 3)) { os<<num+1<<" "; }
    else { os<<ids[lcc.beta(it, 3)]<<" "; }
    os<<"0 0 0 0 "; // 4 bytes => 32 Boolean marks (unsaved for now)

    os<<"0"<<std::endl; // No point
  }

  os.close();
  return true;
}

}

#endif
