// Copyright (c) 2005, 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// author(s)     : Bastien Manuel
//                 Laurent Rincon
//                 Jerome Piovano
//                 Raphaelle Chaine (Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr)
//
// ======================================================================


#ifndef CGAL_PARSER_CNT_H
#define CGAL_PARSER_CNT_H

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>

#include <cmath>

#include <CGAL/Parser_SLC.h>

CGAL_BEGIN_NAMESPACE

template <class handler_CNT>
//handler_CNT class supporting insert_slice(numslice,plane) insert(Point,numslice) 
// with vertex supporting set_next
class parser_CNT : public parser_SLC<handler_CNT>
{
  /*==================================*/
  /*         Type Definitions         */
  /*==================================*/
  typedef parser_CNT<handler_CNT> CNT_p;
  typedef parser_SLC<handler_CNT> SLC_p;
  typedef typename handler_CNT::Plane Plane;
  typedef typename handler_CNT::Vertex_handle Vertex_handle;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using SLC_p::ifs;
  using SLC_p::ignore_blanks;
#endif
  
 public:
  
  /*==================================*/
  /*     Public Members Functions     */
  /*==================================*/
  
  parser_CNT(const char * filename,handler_CNT * slch);
  virtual ~parser_CNT(){;}
private :
    
  /*==================================*/
  /*     Private Members data         */
  /*==================================*/

  int current_group;
  Vertex_handle first_vertex_in_group;
  Vertex_handle previous_vertex_in_group;
    
  /*==================================*/
  /*     Private Members Functions    */
  /*==================================*/
  virtual void read_slice_body(int num_slice, int verticesNumber, const Plane & pl) throw (Parse_exception);
  virtual Vertex_handle read_point(int num_slice, const Plane & pl) throw (Parse_exception);
};


/*==================================================*/


template<class handler_CNT>
parser_CNT<handler_CNT>::parser_CNT(const char * filename, handler_CNT * slch)
  : SLC_p(filename,slch)
{
  current_group = 0;
  first_vertex_in_group=NULL;
  previous_vertex_in_group=NULL;
}

template<class handler_CNT>
void parser_CNT<handler_CNT>::read_slice_body(int num_slice, int verticesNumber, const Plane & pl) throw (Parse_exception)
{
  char c;
  int i;
  current_group = 0;
  for(i=0; (i<verticesNumber) && (!ifs.eof());)
    {
      ignore_blanks();
      if(ifs.peek() == '{')
	{
	  ifs.get(c);
#ifdef CGAL_DUMP
	  std::cout<<"\t\t\t\tParsing vertices set..."<<std::endl;
#endif
	  first_vertex_in_group = previous_vertex_in_group = NULL;
	  read_point_set(num_slice,i,verticesNumber,pl);
	  
	  previous_vertex_in_group->set_next(first_vertex_in_group);
	  current_group++;

#ifdef CGAL_DUMP
	  std::cout<<"\t\t\t\tPoint set parsed successfully"<<std::endl;
#endif
	}
    }
}

template<class handler_CNT>
typename parser_CNT<handler_CNT>::Vertex_handle parser_CNT<handler_CNT>::read_point(int num_slice, const Plane & pl) throw (Parse_exception)
{  
  Vertex_handle v=SLC_p::read_point(num_slice,pl);
  if(previous_vertex_in_group!=NULL)
    previous_vertex_in_group->set_next(v);
  else
    first_vertex_in_group=v;
  previous_vertex_in_group=v;	
  return v;
}
CGAL_END_NAMESPACE
#endif // CGALPARSER_CNT_H
