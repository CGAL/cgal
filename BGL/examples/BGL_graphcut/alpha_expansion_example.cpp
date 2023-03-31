#include <CGAL/boost/graph/alpha_expansion_graphcut.h>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>

struct Vertex_property
{
int label;
std::vector<double> cost;
};

struct Edge_property
{
double weight;
};

using Graph = boost::adjacency_list <boost::setS,
boost::vecS,
boost::undirectedS,
Vertex_property,
Edge_property>;
using GT = boost::graph_traits<Graph>;
using vertex_descriptor = GT::vertex_descriptor;
using edge_descriptor = GT::edge_descriptor;

int main()
{
std::array<char, 3> labels = { 'X', ' ', 'O' };

std::array<std::array<int, 6>, 5> input
= { { { 0, 2, 0, 1, 1, 1 },
{ 0, 0, 1, 0, 1, 2 },
{ 2, 0, 1, 1, 2, 2 },
{ 0, 1, 1, 2, 2, 0 },
{ 1, 1, 2, 0, 2, 2 } } };

std::array<std::array<vertex_descriptor, 6>, 5> vertices;

// Init vertices from values
Graph g;
for (std::size_t i = 0; i < input.size(); ++ i)
for (std::size_t j = 0; j < input[i].size(); ++ j)
{
vertices[i][j] = boost::add_vertex(g);
g[vertices[i][j]].label = input[i][j];
  // Cost of assigning this vertex to any label is positive except
  // for current label which is 0 (favor init solution)
  g[vertices[i][j]].cost.resize(3, 1);
  g[vertices[i][j]].cost[std::size_t(input[i][j])] = 0;
}
// Display input values
std::cerr << "Input:" << std::endl;
for (std::size_t i = 0; i < vertices.size(); ++ i)
{
for (std::size_t j = 0; j < vertices[i].size(); ++ j)
std::cerr << labels[std::size_t(g[vertices[i][j]].label)];
std::cerr << std::endl;
}

// Init adjacency
double weight = 0.5;
for (std::size_t i = 0; i < vertices.size(); ++ i)
for (std::size_t j = 0; j < vertices[i].size(); ++ j)
{
// Neighbor vertices are connected
if (i < vertices.size() - 1)
{
edge_descriptor ed = boost::add_edge (vertices[i][j], vertices[i+1][j], g).first;
g[ed].weight = weight;
}
if (j < vertices[i].size() - 1)
{
edge_descriptor ed = boost::add_edge (vertices[i][j], vertices[i][j+1], g).first;
g[ed].weight = weight;
}
}

std::cerr << std::endl << "Alpha expansion..." << std::endl << std::endl;

// Change alpha expansion parameter
double alpha = 0.9;
// Modify vertex labels to initialize the expansion
for (std::size_t i = 0; i < input.size(); ++ i)
for (std::size_t j = 0; j < input[i].size(); ++ j)
g[vertices[i][j]].label = (input[i][j] + 1) % 3;

// Display initial labels
std::cerr << "Initial labels:" << std::endl;
for (std::size_t i = 0; i < vertices.size(); ++ i)
{
for (std::size_t j = 0; j < vertices[i].size(); ++ j)
std::cerr << labels[std::size_t(g[vertices[i][j]].label)];
std::cerr << std::endl;
}

// Run alpha expansion graph cut with modified vertex labels
std::cerr << std::endl << "Alpha expansion with modified labels..." << std::endl << std::endl;
CGAL::alpha_expansion_graphcut (g,
get (&Edge_property::weight, g),
get (&Vertex_property::cost, g),
get (&Vertex_property::label, g),
CGAL::parameters::vertex_index_map (get (boost::vertex_index, g)));

// Display output graph
std::cerr << "Output:" << std::endl;
for (std::size_t i = 0; i < vertices.size(); ++ i)
{
for (std::size_t j = 0; j < vertices[i].size(); ++ j)
std::cerr << labels[std::size_t(g[vertices[i][j]].label)];
std::cerr << std::endl;
}

return 0;
}
