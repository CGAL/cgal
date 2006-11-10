// ======================================================================
//
// file          : include/CGAL/Parser_SLC.h
// package       : Reconstruction_from_slices
//
// author(s)     : Bastien Manuel
//                 Laurent Rincon
//                 Jerome Piovano
//                 Raphaelle Chaine (Raphaelle.Chaine@sophia.inria.fr,
//                                   raphaelle.chaine@liris.cnrs.fr)
//
// ======================================================================

#ifndef CGAL_PARSER_SLC_H
#define CGAL_PARSER_SLC_H

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>

#include <cmath>

CGAL_BEGIN_NAMESPACE

class Exception;
class Number_exception;
class Parse_exception;

template <class handler_SLC>
//handler_SLC class supporting insert_slice(numslice,plane) insert(Point,numslice)
class parser_SLC
{
  /*==================================*/
  /*         Type Definitions         */
  /*==================================*/
  typedef typename handler_SLC::Point Point;
  typedef typename handler_SLC::Vector Vector;
  typedef typename handler_SLC::Plane Plane;
  typedef typename handler_SLC::Vertex_handle Vertex_handle;
  
 public:
  
  /*==================================*/
  /*     Public Members Functions     */
  /*==================================*/
  
  parser_SLC(const char * filename,handler_SLC * slch);
  virtual ~parser_SLC();
  void parse() throw (Parse_exception);
  void setInput(const char * filename);
  
  template<class SLC_H>
    friend std::ostream & operator <<(std::ostream & os, const parser_SLC<SLC_H> & parser);
  
 enum
  {DEFAULT,PARSING,PARSED,FAILED};  
 
protected :
    
  /*==================================*/
  /*     Protected Members data         */
  /*==================================*/

  std::string input;
  std::ifstream ifs;
  int status;
  bool parallel;
  int nb_slices;
  int current_slice;
  handler_SLC * slch;

  /*==================================*/
  /*     Protected Members Functions    */
  /*==================================*/
    
  void read_header() throw (Parse_exception);
  void read_slice_style() throw (Parse_exception);
  void read_slices_number() throw (Parse_exception);

  void read_body() throw (Parse_exception);
  void read_slice(int num_slice) throw (Parse_exception);
  void read_slice_header(int num_slice, int & verticesNumber, Plane & pl) throw (Parse_exception);
  int  read_vertices_number() throw (Parse_exception);
  double read_z_value() throw (Parse_exception);
  Plane read_slice_plane() throw (Parse_exception);
  virtual void read_slice_body(int num_slice, int verticesNumber, const Plane & pl) throw (Parse_exception);
  virtual void read_point_set(int num_slice, int & currentVertexNumber, int verticesNumber, const Plane & pl) throw (Parse_exception);
  virtual Vertex_handle read_point(int num_slice, const Plane & pl) throw (Parse_exception);
  inline int read_integer() throw (Number_exception);
  inline double read_double() throw (Number_exception);
  inline void ignore_blanks();
};

class Exception
{
public:

  Exception();
  Exception(std::string &);
  Exception(const Exception &);
  std::string getMessage() const;
  friend std::ostream & operator <<(std::ostream &, const Exception &);
  
private:
  std::string message;
};

class Number_exception : public Exception
{
public:
  
  Number_exception(): Exception(){};
  Number_exception(std::string &msg) : Exception(msg){};
  
  friend std::ostream & operator <<(std::ostream &, const Number_exception &);
  
};

class Parse_exception : public Exception
{
public:

  Parse_exception():Exception(){};
  Parse_exception(std::string &msg):Exception(msg){};
  friend std::ostream & operator <<(std::ostream &, const Parse_exception &);
};

/*===========================================================================*/
/*         Implementation                                                    */
/*===========================================================================*/


std::ostream & operator <<(std::ostream & os, const Parse_exception & e)
{
  os<<"PARSER EXCEPTION : "<<e.getMessage()<<std::endl;
  return os;
}

std::ostream & operator <<(std::ostream & os, const Number_exception & e)
{
  os<<"NUMBER EXCEPTION : "<<e.getMessage()<<std::endl;
  return os;
}

Exception::Exception()
{
  message = "No description message";
}

Exception::Exception(std::string &msg)
{
  message = msg;
}

Exception::Exception(const Exception & e)
{
  message = e.message;
}

std::string Exception::getMessage() const
{
  return message;
}

std::ostream & operator <<(std::ostream & os, const Exception & e)
{
  os<<"EXCEPTION : "<<e.message<<std::endl;
  return os;
}

template<class handler_SLC>
parser_SLC<handler_SLC>::parser_SLC(const char * filename, handler_SLC * slch)
{
  setInput(filename);
  status = DEFAULT;
  current_slice = 0;
  this->slch = slch;
}

template<class handler_SLC>
parser_SLC<handler_SLC>::~parser_SLC()
{
  ifs.close();
}

template<class handler_SLC>
void parser_SLC<handler_SLC>::parse() throw (Parse_exception)
{
  status = PARSING;
#ifdef DUMP
  std::cout<<std::endl<<"\tParsing header..."<<std::endl;
#endif
  read_header();
#ifdef DUMP
  std::cout<<"\tHeader successfully parsed"<<std::endl; 
  std::cout<<"\t=>Slice type = ";
  if(parallel)
    std::cout<<"parallel";
  else
    std::cout<<"non parallel";
  std::cout<<std::endl;
  std::cout<<"\t=>Slices number = "<<nb_slices<<std::endl;
  std::cout<<std::endl;
  std::cout<<"\tParsing body..."<<std::endl;
#endif
  read_body();
#ifdef DUMP
  std::cout<<"\tBody successfully parsed"<<std::endl;
  std::cout<<std::endl;
#endif
  status = PARSED;
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::setInput(const char * filename)
{
  input = filename;
  ifs.open(input.c_str(), std::ios::in);
}


template<class handler_SLC>
std::ostream & operator <<(std::ostream & os, const parser_SLC<handler_SLC> & parser)
{
  os<<"PARSER : "<<std::endl;
  os<<"  FILE : "<<parser.input<<std::endl;
  os<<"  STATUS : ";
  switch(parser.status)
    {
    case parser_SLC<handler_SLC>::PARSING : os<<"PARSING..."; break;
    case parser_SLC<handler_SLC>::PARSED : os<<"PARSING DONE"; break;
    case parser_SLC<handler_SLC>::DEFAULT : os<<"PARSING NOT DONE"; break;
    case parser_SLC<handler_SLC>::FAILED : os<<"PARSING FAILED"; break;
    }
  return os;
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::read_header() throw (Parse_exception)
{
#ifdef DUMP
  std::cout<<"\t\tParsing slice style..."<<std::endl;
#endif
  read_slice_style();
#ifdef DUMP
  std::cout<<"\t\tSlice style parsed"<<std::endl;
  std::cout<<std::endl;
  std::cout<<"\t\tParsing slices number..."<<std::endl;
#endif
  read_slices_number();
#ifdef DUMP
  std::cout<<"\t\tSlices number parsed"<<std::endl;
  std::cout<<std::endl;
#endif
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::read_slice_style() throw (Parse_exception)
{
  bool exit = false;
  int k = 0;
  // the SLC file must have his first line like this
  //   SLC n 
  //   with n = 0  if slc file of slices parallel to xy plane
  //     or n = 1  else
  while(!exit)
    {
      char cc[5];
      char c;
      ignore_blanks();
      c = ifs.peek();
      if( ( ((k == 0) && ((c == 'S') || (c == 's')))
	    ||((k == 1) && ((c == 'L') || (c == 'l')))
	    ||((k == 2) && ((c == 'C') || (c == 'c')))
	    ||((k == 3) && ((c == '0') || (c == '1')))
	    )
	  && (!ifs.eof()) )
	{
	  ifs.get(cc[k]);
	  k++;
	}
      else
	{
	  status = FAILED;
	  std::string msg = "";
	  msg += ("Badly formatted header : \n");
	  msg += ("  SLC followed by type of slices was expected.");
	  throw (Parse_exception(msg));
	}
      if(k == 4) 
	{
	  parallel = (cc[3]=='0');
	  exit = true;
	}
    }
}



template<class handler_SLC>
void parser_SLC<handler_SLC>::read_slices_number() throw (Parse_exception)
{ 
  char c;
  ignore_blanks();
  c = ifs.peek();
  if((c == 'S') || (c == 's'))
    {
      ifs.get(c);
      ignore_blanks();
      try
	{
	  nb_slices = read_integer();
	}
      catch(Number_exception ne)
	{
	  status = FAILED;
	  std::string msg = "";
	  msg += "Badly formatted header : \n";
	  msg += "  " + ne.getMessage();
	  throw Parse_exception(msg);
	}
    }
  else
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted header : \n";
      msg += "  S followed by number of slices was expected.";
      throw Parse_exception(msg);   
    }
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::read_body() throw (Parse_exception)
{  
  while((current_slice < nb_slices)
	&& (!ifs.eof()))
    {
#ifdef DUMP
      std::cout<<"\t\tParsing slice n°"<<current_slice<<"..."<<std::endl;
#endif
      read_slice(current_slice);
      current_slice++;
#ifdef DUMP
      std::cout<<"\t\tSlice n°"<<current_slice<< "/" << nb_slices <<" parsed successfully"<<std::endl;
#endif
    }
  ignore_blanks();
  if(!ifs.eof())
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted file : \n";
      msg += "  Too many slices.";
      throw Parse_exception(msg);
    }
  if(ifs.eof() && (current_slice < nb_slices))
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted file : \n";
      msg += "  Not enough slices.";
      throw Parse_exception(msg);
    }
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::read_slice(int current_slice) throw (Parse_exception)
{
  int verticesNumber;
  Plane pl;
#ifdef DUMP
  std::cout<<"\t\t\tParsing slice header..."<<std::endl;
#endif
  read_slice_header(current_slice,verticesNumber,pl);
#ifdef DUMP
  std::cout<<"\t\t\tSlice header parsed successfully"<<std::endl;
  std::cout<<"\t\t=>Vertices number = "<<verticesNumber<<std::endl;
  std::cout<<"\t\t=>Plane = "<<pl<<std::endl;
#endif
  slch->insert_slice(current_slice,pl);
#ifdef DUMP
  std::cout<<"\t\t\tParsing slice body..."<<std::endl;
#endif
  read_slice_body(current_slice,verticesNumber,pl);
#ifdef DUMP
  std::cout<<"\t\t\tSlice body parsed successfully"<<std::endl;  
#endif
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::read_slice_header(int current_slice, int & verticesNumber, typename handler_SLC::Plane & pl) throw (Parse_exception)
{
#ifdef DUMP
  std::cout<<"\t\t\t\tParsing vertices number..."<<std::endl;
#endif
  verticesNumber = read_vertices_number();
#ifdef DUMP
  std::cout<<"\t\t\t\tVertices number parsed successfully"<<std::endl;
  std::cout<<"\t\t\t\tParsing slice plane"<<std::endl;
#endif
  if (parallel)
    {
      double zValue = read_z_value();
      pl=Plane(Point(0,0,zValue),Vector(0,0,1));
    }
  else
    pl=read_slice_plane();
#ifdef DUMP
  std::cout<<"\t\t\t\tSlice plane parsed successfully"<<std::endl;
#endif
}


template <class handler_SLC>
int parser_SLC<handler_SLC>::read_vertices_number() throw (Parse_exception)
{
  char c;
  ignore_blanks();
  c = ifs.peek();
  if(!ifs.eof() && ((c == 'v') || (c == 'V')))
    try
      {
	ifs.get(c);
	return read_integer();
      }
  catch(Number_exception ne)
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted slice :\n";
      msg += "  Number of vertices :\n";
      msg += "  " + ne.getMessage();
      throw Parse_exception(msg);
    }
  status = FAILED;
  std::string msg = "";
  msg += "Badly formatted slice :\n";
  msg += "  Expected v followed by a number of vertex.";
  throw Parse_exception(msg);
}


template <class handler_SLC>
double parser_SLC<handler_SLC>::read_z_value() throw (Parse_exception)
{
  char c;
  
  ignore_blanks();
  c = ifs.peek();
  if(!ifs.eof() && ((c == 'z') || (c == 'Z')))
    try
      {
	ifs.get(c);
	return read_double();
      }
  catch(Number_exception ne)
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted slice :\n";
      msg += "  Z Value :\n";
      msg += "  ";
      msg += ne.getMessage();
      throw Parse_exception(msg);
    }
  status = FAILED;
  std::string msg = "";
  msg += "Badly formatted slice :\n";
  msg += "  Expected z followed by a double value.";
  throw Parse_exception(msg);
}

template <class handler_SLC>
typename handler_SLC::Plane parser_SLC<handler_SLC>::read_slice_plane() throw (Parse_exception)
{ 
  Point origin;
  Vector v;
  double x, y, z;
  char c;  
 
  ignore_blanks();
  c = ifs.peek();
  if(!ifs.eof() && ((c == 'o') || (c == 'O')))
    try
      {
	ifs.get(c);
	ignore_blanks();
	x = read_double();
	ignore_blanks();
	y = read_double();
	ignore_blanks();
	z = read_double();
	origin=Point(x,y,z);
      }
  catch(Number_exception ne)
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted slice :\n";
      msg += "  O Value :\n";
      msg += "  ";
      msg += ne.getMessage();
      throw Parse_exception(msg);
    }
  
  ignore_blanks();
  c = ifs.peek();
  if(!ifs.eof() && ((c == 'n') || (c == 'N')))
    try
      {
	ifs.get(c);
	ignore_blanks();
	x = read_double();
	ignore_blanks();
	y = read_double();
	ignore_blanks();
	z = read_double();
	v=Vector(x,y,z);
	return Plane(origin,v);
      }
  catch(Number_exception ne)
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted slice :\n";
      msg += "  N Value :\n";
      msg += "  ";
      msg += ne.getMessage();
      throw Parse_exception(msg);
    }
  status = FAILED;
  std::string msg = "";
  msg += "Badly formatted slice :\n";
  msg += "Expected o followed by a point value + n followed by a vector value\n";
  throw Parse_exception(msg);
}

template<class handler_SLC>
void parser_SLC<handler_SLC>::read_slice_body(int num_slice, int verticesNumber, const Plane & pl) throw (Parse_exception)
{
  char c;
  int i;
  for(i=0; (i<verticesNumber) && (!ifs.eof());)
    {
      ignore_blanks();
      if(ifs.peek() == '{')
	{
	  ifs.get(c);
#ifdef DUMP
	  std::cout<<"\t\t\t\tParsing vertices set..."<<std::endl;
#endif
	  read_point_set(num_slice,i,verticesNumber,pl);
#ifdef DUMP
	  std::cout<<"\t\t\t\tPoint set parsed successfully"<<std::endl;
#endif
	}
    }
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::read_point_set(int num_slice, int & currentVertexNumber, int verticesNumber, const Plane & pl) throw (Parse_exception)
{
  char c;
  while(!ifs.eof() && (ifs.peek() != '}') && (currentVertexNumber < verticesNumber))
    {
#ifdef DUMP
      std::cout<<"\t\t\t\t\tParsing vertex n°"<< currentVertexNumber+1 << "/" << verticesNumber << "..."<<std::endl;
#endif
      read_point(num_slice,pl);
#ifdef DUMP
      std::cout<<"\t\t\t\t\tVertex parsed successfully"<<std::endl;
#endif
      currentVertexNumber++;
    }
  ignore_blanks();

if(!ifs.eof())
    {
      if(ifs.peek() == '}')
	{
	  ifs.get(c);
	  ignore_blanks();
	  if(currentVertexNumber<verticesNumber)
	    {
	      if(ifs.peek() != '{')
		{
		  status = FAILED;
		  std::string msg = "";
		  msg += "Badly formatted slice :\n";
		  msg += "  Not enough vertices for this slice.";
		  throw Parse_exception(msg);
		}
	    }
	  else
	    {
	      if(ifs.peek() == '{')
		{
		  status = FAILED;
		  std::string msg = "";
		  msg += "Badly formatted slice :\n";
		  msg += "  Too many vertices for this slice.";
		  throw Parse_exception(msg);
		}
	    }
	}
      else
	{
	  status = FAILED;
	  std::string msg = "";
	  msg += "Badly formatted slice :\n";
	  msg += "  Expected end of set : \'}\'.";
	  throw Parse_exception(msg);
	}
    }
  else
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted slice :\n";
      msg += "  Unexpected end of file : \'}\' expected.";
      throw Parse_exception(msg);
    }
}


template<class handler_SLC>
typename parser_SLC<handler_SLC>::Vertex_handle parser_SLC<handler_SLC>::read_point(int num_slice, const Plane & pl) throw (Parse_exception)
{
  double x, y, z=0;  
  try
    {
      ignore_blanks();
      x = read_double();
      ignore_blanks();
      y = read_double();
      ignore_blanks();
      if(!parallel)
	{
	  z = read_double();
	  ignore_blanks();
	}      
#ifdef DUMP
      std::cout<<"\t\t\t\t\tx = "<<x<<" y = "<<y<<" z = "<<z<<std::endl;
      std::cout<<"\t\t\t\t\t "<<pl.projection(Point(x,y,z))<<std::endl;
#endif
      if(!parallel)     
	return slch->insert(pl.projection(Point(x,y,z)),current_slice);
      else
	return slch->insert(Point(x,y,-pl.d()),current_slice);
    }
  catch(Number_exception ne)
    {
      status = FAILED;
      std::string msg = "";
      msg += "Badly formatted vertex :\n";
      msg += "  ";
      msg += ne.getMessage();
      throw Parse_exception(msg);
    }
}


template<class handler_SLC>
int parser_SLC<handler_SLC>::read_integer() throw (Number_exception)
{
  char number[50];
  char c;
  int i = 0;
  
  ignore_blanks();
  c = ifs.peek();
  if(c == '-'){
    ifs.get(c);
    number[i++] = c; 
  }
  
  c = ifs.peek();
  if(!isdigit(c)){
    std::string msg = "Badly formatted number.\n";
    throw Number_exception(msg);
  }
  
  while(isdigit(ifs.peek())){
    ifs.get(c);
    number[i++] = c;
  }
  
  return std::atoi(number);
}


template<class handler_SLC>
double parser_SLC<handler_SLC>::read_double() throw (Number_exception)
{
  char number[50];
  char c;
  int i = 0;
  
  ignore_blanks();
  c = ifs.peek();
  if(c == '-')
    {
      ifs.get(c);
      number[i++] = c; 
    }
  
  c = ifs.peek();
  if(!isdigit(c))
    {
      std::string msg = "Badly formatted number.\n";
      throw Number_exception(msg);
    }
  
  while(isdigit(ifs.peek()))
    {
      ifs.get(c);
      number[i++] = c;
    }
  
  c = ifs.peek();
  if(c == '.')
    {
      ifs.get(c);
      number[i++] = c; 
    }
  
  while(isdigit(ifs.peek()))
    {
      ifs.get(c);
      number[i++] = c;
    }
  number[i] = '\0';
  return std::atof(number);
}


template<class handler_SLC>
void parser_SLC<handler_SLC>::ignore_blanks()
{
  char c;
  while( ifs.peek()=='\n' || ifs.peek()=='\t' || ifs.peek()=='\r' || ifs.peek()==' ' )
    ifs.get(c);
}

CGAL_END_NAMESPACE

#endif // CGAL_PARSER_SLC_H
 
