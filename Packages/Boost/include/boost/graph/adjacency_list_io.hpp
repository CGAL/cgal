//=======================================================================
// Copyright 2001 Universite Joseph Fourier, Grenoble.
// Author: François Faure
//
// This file is part of the Boost Graph Library
//
// You should have received a copy of the License Agreement for the
// Boost Graph Library along with the software; see the file LICENSE.
// If not, contact Office of Research, University of Notre Dame, Notre
// Dame, IN 46556.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================


#ifndef ______adj_list_io_______
#define ______adj_list_io_______

#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <cctype>

// Method read to parse an adjacency list from an input stream. Examples:
// cin >> read( G );
// cin >> read( G, NodePropertySubset(), EdgepropertySubset() );
//
// Method write to print an adjacency list to an output stream. Examples:
// cout << write( G );
// cout << write( G, NodePropertySubset(), EdgepropertySubset() );

namespace boost {

/* outline
        - basic property input
        - get property subset
        - graph parser
        - property printer
        - graph printer
        - user methods
*/

//===========================================================================
// basic property input

template<class Tag, class Value, class Next>
std::istream& operator >> ( std::istream& in, property<Tag,Value,Next>& p )
{
        in >> p.m_value >> *(static_cast<Next*>(&p)); // houpla !!
        return in;
}

template<class Tag, class Value>
std::istream& operator >> ( std::istream& in, property<Tag,Value,no_property>& p )
{
        in >> p.m_value;
        return in;
}

inline std::istream& operator >> ( std::istream& in, no_property& )
{
        return in;
}

// basic property input
//===========================================================================
// get property subsets

// get a single property tagged Stag
template<class Tag, class Value, class Next, class V, class Stag>
void get
( property<Tag,Value,Next>& p, const V& v, Stag s )
{
        get( *(static_cast<Next*>(&p)),v,s );
}

template<class Value, class Next, class V, class Stag>
void get
( property<Stag,Value,Next>& p, const V& v, Stag )
{
        p.m_value = v;
}

// get a subset of properties tagged Stag
template<class Tag, class Value, class Next, 
        class Stag, class Svalue, class Snext>
void getSubset
( property<Tag,Value,Next>& p, const property<Stag,Svalue,Snext>& s )
{
        get( p, s.m_value, Stag() );
        getSubset( p, Snext(s) );
}

template<class Tag, class Value, class Next, 
        class Stag, class Svalue>
void getSubset
( property<Tag,Value,Next>& p, const property<Stag,Svalue,no_property>& s )
{
        get( p, s.m_value, Stag() );
}

inline void getSubset
( no_property& p, const no_property& s )
{
}

//  get property subset
//===========================================================================
// graph parser

template<class Graph_t, class VertexProperty, class EdgeProperty, class VertexPropertySubset,
class EdgePropertySubset>
struct GraphParser
{

        typedef Graph_t Graph;
        
        GraphParser( Graph* g ): graph(g)
        {}      
        
        GraphParser& operator () ( std::istream& in )
        {
                typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
                std::vector<Vertex> nodes;

                typedef enum{ PARSE_NUM_NODES, PARSE_VERTEX, PARSE_EDGE } State;
                State state = PARSE_VERTEX;

                unsigned int numLine = 1;
                char c;
                while ( in.get(c) )
                {
                        if( c== '#' ) skip(in);
                        else if( c== 'n' ) state = PARSE_NUM_NODES;
                        else if( c== 'v' ) state = PARSE_VERTEX;
                        else if( c== 'e' ) state = PARSE_EDGE;
                        else if( c== '\n' ) numLine++;
                        else if( !std::isspace(c) ){
                                in.putback(c);
                                if( state == PARSE_VERTEX ){
                                        VertexPropertySubset readProp;
                                        if( in >> readProp )
                                        {
                                                VertexProperty vp;
                                                getSubset( vp, readProp );
                                                nodes.push_back( add_vertex(vp, *graph) );
                                        }
                                        else
                                                std::cerr<<"read vertex, parse error at line"<<numLine<<std::endl;
                                }
                                else if( state == PARSE_EDGE ) {
                                        int source, target;
                                        EdgePropertySubset readProp;
                                        in >> source >> target;
                                        if( in >> readProp ) 
                                        {
                                                EdgeProperty ep;
                                                getSubset( ep, readProp );
                                                add_edge(nodes[source], nodes[target], ep, *graph);
                                        }
                                        else
                                                std::cerr<<"read edge, parse error at line"<<numLine<<std::endl;
                                }
                                else { // state == PARSE_NUM_NODES
                                        int n;
                                        if( in >> n ){
                                                for( int i=0; i<n; ++i )
                                                        nodes.push_back( add_vertex( *graph ));
                                        }
                                        else 
                                                std::cerr<<"read num_nodes, parse error at line "<< numLine << std::endl;
                                }
                        }
                }
        return (*this);
        }
        
        
protected:

        Graph* graph;
        
        void skip( std::istream& in )
        {
                char c = 0;
                while( c!='\n' && !in.eof() ) 
                       in.get(c);
                in.putback(c);
        }
};

// parser
//=======================================================================
// property printer

template<class Graph, class Property>
struct PropertyPrinter
{
        typedef typename Property::value_type Value;
        typedef typename Property::tag_type Tag;
        typedef typename Property::next_type Next;
        
        PropertyPrinter( Graph& g ):graph(&g){}
        
        template<class Iterator>
        PropertyPrinter& operator () ( std::ostream& out, Iterator it )
        {
                typename property_map<Graph,Tag>::type ps = get(Tag(), *graph);
                out << ps[ *it ] <<" ";
                PropertyPrinter<Graph,Next> print(*graph);
                print(out, it);
                return (*this);
        }
private:
        Graph* graph;
};
template<class Graph>
struct PropertyPrinter<Graph, no_property>
{
        PropertyPrinter( Graph& ){}

        template<class Iterator>
        PropertyPrinter& operator () ( std::ostream&, Iterator it ){ return *this; }
};

// property printer
//=========================================================================
// graph printer

template<class Graph_t, class EdgeProperty>
struct EdgePrinter
{

        typedef Graph_t Graph;
        typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
        
        EdgePrinter( Graph& g )
                : graph(g)
        {}      
        
        const EdgePrinter& operator () ( std::ostream& out ) const
        {
                // assign indices to vertices
                std::map<Vertex,int> indices;
                int num = 0;
                typename graph_traits<Graph>::vertex_iterator vi;
                for (vi = vertices(graph).first; vi != vertices(graph).second; ++vi){
                        indices[*vi] = num++;
                }

                // write edges
                PropertyPrinter<Graph, EdgeProperty> print_Edge(graph);
                out << "e" << std::endl;
                typename graph_traits<Graph>::edge_iterator ei;
                for (ei = edges(graph).first; ei != edges(graph).second; ++ei){
                        out << indices[source(*ei,graph)] <<  " " << indices[target(*ei,graph)] << "  "; 
                        print_Edge(out,ei); 
                        out << std::endl;
                }
                out << std::endl;            
                return (*this);
        }
        
protected:

        Graph& graph;
        
};

template<class Graph, class V, class E>
struct GraphPrinter: public EdgePrinter<Graph,E>
{
        GraphPrinter( Graph& g )
          : EdgePrinter<Graph,E>(g)
        {}
        
        const GraphPrinter& operator () ( std::ostream& out ) const
        {
                PropertyPrinter<Graph, V> printNode(this->graph);
                out << "v"<<std::endl;
                typename graph_traits<Graph>::vertex_iterator vi;
                for (vi = vertices(this->graph).first; vi != vertices(this->graph).second; ++vi){
                        printNode(out,vi); 
                        out << std::endl;
                }
                
                EdgePrinter<Graph,E>::operator ()( out );
                return (*this);
        }
};

template<class Graph, class E>
struct GraphPrinter<Graph,no_property,E> 
  : public EdgePrinter<Graph,E>
{
        GraphPrinter( Graph& g )
          : EdgePrinter<Graph,E>(g)
        {}
        
        const GraphPrinter& operator () ( std::ostream& out ) const
        {
                out << "n "<< num_vertices(this->graph) << std::endl;
                EdgePrinter<Graph,E>::operator ()( out );
                return (*this);
        }
};

// graph printer
//=========================================================================
// user methods

/// input stream for reading a graph
template<class Graph, class VP, class EP, class VPS, class EPS>
std::istream& operator >> ( std::istream& in, GraphParser<Graph,VP,EP,VPS,EPS> gp ) 
{ 
        gp(in); 
        return in; 
}

/// graph parser for given subsets of internal vertex and edge properties
template<class EL, class VL, class D, class VP, class EP, class GP, class VPS, class EPS>
GraphParser<adjacency_list<EL,VL,D,VP,EP,GP>,VP,EP,VPS,EPS> 
read( adjacency_list<EL,VL,D,VP,EP,GP>& g, VPS vps, EPS eps )
{
        return GraphParser<adjacency_list<EL,VL,D,VP,EP,GP>,VP,EP,VPS,EPS>(&g);
}

/// graph parser for all internal vertex and edge properties
template<class EL, class VL, class D, class VP, class EP, class GP>
GraphParser<adjacency_list<EL,VL,D,VP,EP,GP>,VP,EP,VP,EP> 
read( adjacency_list<EL,VL,D,VP,EP,GP>& g )
{
        return GraphParser<adjacency_list<EL,VL,D,VP,EP,GP>,VP,EP,VP,EP>(&g);
}


/// output stream for writing a graph
template<class Graph, class VP, class EP>
std::ostream& operator << ( std::ostream& out, const GraphPrinter<Graph,VP,EP>& gp ) 
{ 
        gp(out); 
        return out; 
}

/// write the graph with given property subsets
template<class EL, class VL, class D, class VP, class EP, class GP, class VPS, class EPS>
GraphPrinter<adjacency_list<EL,VL,D,VP,EP,GP>,VPS,EPS> 
write( adjacency_list<EL,VL,D,VP,EP,GP>& g, VPS, EPS )
{
        return GraphPrinter<adjacency_list<EL,VL,D,VP,EP,GP>,VPS,EPS>(g);
}

/// write the graph with all internal vertex and edge properties
template<class EL, class VL, class D, class VP, class EP, class GP>
GraphPrinter<adjacency_list<EL,VL,D,VP,EP,GP>,VP,EP> 
write( adjacency_list<EL,VL,D,VP,EP,GP>& g )
{
        return GraphPrinter<adjacency_list<EL,VL,D,VP,EP,GP>,VP,EP>(g);
}

// user methods
//=========================================================================
}// boost
#endif
