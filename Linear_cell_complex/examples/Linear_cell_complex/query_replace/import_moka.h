// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
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

namespace CGAL
{
template<typename LCC>
struct GDart
{
  unsigned int alpha[4];
  typename LCC::Dart_handle dh;
  typename LCC::Vertex_attribute_handle vh;

  GDart() : dh(nullptr), vh(nullptr)
  {}

  GDart(const GDart& adart) : dh(adart.dh),
    vh(adart.vh)
  {
    for (unsigned int i=0; i<4; ++i)
    { alpha[i]=adart.alpha[i]; }
  }
};

template<typename LCC>
bool import_from_moka(LCC& lcc, const char* filename)
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
    bool orient=(gdarts[startingdart].dh==nullptr);
    for (unsigned int dim=0; orient && dim<4; ++dim)
      if (gdarts[gdarts[startingdart].alpha[dim]].dh!=nullptr) orient=false;

    if ( orient )
    {
      totreat.push(startingdart);
      gdarts[startingdart].dh=lcc.create_dart();

      while ( !totreat.empty() )
      {
        unsigned int i=totreat.top();
        totreat.pop();

        assert(gdarts[i].dh!=nullptr);

        for (unsigned int dim=1; dim<4; ++dim)
        {
          if (gdarts[i].alpha[dim]!=i &&
              gdarts[gdarts[i].alpha[dim]].vh!=nullptr)
          {
            gdarts[i].vh = gdarts[gdarts[i].alpha[dim]].vh;
            gdarts[gdarts[i].alpha[dim]].vh = nullptr;
          }

          unsigned int alpha0 = gdarts[i].alpha[0];
          assert( alpha0!=i );

          if (gdarts[alpha0].alpha[dim]!=alpha0)
          {
            if ( gdarts[gdarts[alpha0].alpha[dim]].dh==nullptr )
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
    if (gdarts[i].dh!=nullptr)
    {
      for (unsigned int dim=0; dim<4; ++dim)
      {
        if (orientable &&
            gdarts[i].alpha[dim]!=i &&
            gdarts[gdarts[i].alpha[dim]].dh!=nullptr)
        {
          std::cout<<"Pb, the gmap is NOT orientable."<<std::endl;
          orientable=false;
          // lcc.clear();
        }
      }

      /* if ( lcc.template attribute<3>(gdarts[i].dh) == nullptr )
      {
        lcc.template set_attribute<3>(gdarts[i].dh, lcc.template create_attribute<3>());
      } */
    }
    if (gdarts[i].vh!=nullptr)
    {
      lcc.set_vertex_attribute(gdarts[i].dh, gdarts[i].vh);
    }
  }

  return true;
}

}

#endif