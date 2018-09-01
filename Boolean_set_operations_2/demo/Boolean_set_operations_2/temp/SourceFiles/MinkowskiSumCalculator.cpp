//#include <fstream>

#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>

#include "MinkowskiSumCalculator.h"
#include "MinkowskiSumModified.h"

using namespace std;

MinkowskiSumCalculator::MinkowskiSumCalculator()
{

}

MinkowskiSumCalculator::~MinkowskiSumCalculator()
{

}

void MinkowskiSumCalculator::setInput(const std::list<PolygonWithHoles>& input)
{
  input_ = input;
  
  cache_.clear();

  //Add to cache
  cache_.insert(make_pair(1, std::make_shared<MinkowskiSumResult>(input)));
}

std::shared_ptr<const MinkowskiSumResult> MinkowskiSumCalculator::getSum(int k)
{
  if (k <= 0 || cache_.size() == 0)
    return std::make_shared<MinkowskiSumResult>();

  //Look for result in cache
  if (cache_.find(k) != cache_.end())
  {
    cout << "Final result: " << endl;
    cache_[k]->print();

    return cache_[k];
  }

  auto inputResult = cache_[1];
  auto decomposedInput = inputResult->getDecomposedSum();

  //Take closest value to the requested one
  int initK = 0;
  std::shared_ptr<MinkowskiSumResult> initSum;
  for(auto it = cache_.begin(); it != cache_.end(); ++it)
  {
    if (it->first <= k)
    {
      initK = it->first;
      initSum = it->second;
    }
  }

  std::list<PolygonWithHoles> currentSum;
  for(auto it = initSum->getSum().begin(); it != initSum->getSum().end(); ++it)
  {
Kernel::FT temp;
temp=(double)initK;    
currentSum.push_back(scalePolygon(*it, temp));
  }

  //Calculate the sum
  for(int i = initK + 1; i <= k; ++i)
  {
    cout << "K = " << i << endl;

    //Decompose current sum
    std::list<Polygon> decomposedPolygons;
    for(auto polygon = currentSum.begin(); polygon != currentSum.end(); ++polygon)
    {
      decomposePolygon(*polygon, decomposedPolygons);
    }

    //Calculate minkowski sum
    CGAL::Minkowski_sum_by_decomposition_2_Modified<CGAL::Small_side_angle_bisector_decomposition_2<Kernel>, Container> minkowskiSum;
    currentSum = minkowskiSum(decomposedInput, decomposedPolygons);

    //Scale it and remove collinear points and add to cache
    std::list<PolygonWithHoles> scaledSum;
    for(auto it = currentSum.begin(); it != currentSum.end(); ++it)
    {
Kernel::FT temp;
temp=(double)1.0/i;      
scaledSum.push_back(removeCollinearPoints(scalePolygon(*it, temp)));
    }

    cache_[i] = std::make_shared<MinkowskiSumResult>(scaledSum);

    //Print the sum to screen
    if (i != k)
      cache_[i]->print();
  }

  //Print final result
  cout << "Final result: " << endl;
  cache_[k]->print();

  return cache_[k];
}


void MinkowskiSumCalculator::decomposePolygon(const PolygonWithHoles& polygon, std::list<Polygon>& result)
{
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
  typedef CGAL::Exact_predicates_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

  CDT triangulation;

  //Add outer boundary vertices and holes
  triangulation.insert(polygon.outer_boundary().vertices_begin(), polygon.outer_boundary().vertices_end());
  for(auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
  {
    triangulation.insert(hole->vertices_begin(), hole->vertices_end());
  }

  //Add constraints
  for(auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
  {
    for(auto it = hole->edges_begin(); it != hole->edges_begin(); ++it)
    {
      triangulation.insert_constraint(it->source(), it->target());
    }
  }

  //Now go over all triangles and check them
  for(auto it = triangulation.finite_faces_begin(); it != triangulation.finite_faces_end(); ++it)
  {
    //Calculate mid point of the triangle
    auto point = Point(
      (it->vertex(0)->point().x() + it->vertex(1)->point().x() + it->vertex(2)->point().x()) / 3,
      (it->vertex(0)->point().y() + it->vertex(1)->point().y() + it->vertex(2)->point().y()) / 3);

    //Check if it is inside the outer boundary
    if (CGAL::bounded_side_2(polygon.outer_boundary().vertices_begin(), polygon.outer_boundary().vertices_end(), point) == CGAL::ON_UNBOUNDED_SIDE)
    {
      //Don't add this triangle
      continue;
    }

    //Check if it is in one of the holes
    bool isInHole = false;
    for(auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
    {
      if (CGAL::bounded_side_2(hole->vertices_begin(), hole->vertices_end(), point) == CGAL::ON_BOUNDED_SIDE)
      {
        isInHole = true;
        break;
      }
    }

    if (!isInHole)
    {
      //cout << "---" << endl;
      //cout << it->vertex(0)->point() << "," << it->vertex(1)->point() << "," << it->vertex(2)->point() << endl;
      //cout << point << endl;
      Polygon newPolygon;
      newPolygon.push_back(it->vertex(0)->point());
      newPolygon.push_back(it->vertex(1)->point());
      newPolygon.push_back(it->vertex(2)->point());

      result.push_back(newPolygon);
    }
  }
}


PolygonWithHoles MinkowskiSumCalculator::removeCollinearPoints(const PolygonWithHoles& polygon)
{
  PolygonWithHoles newPolygon(removeCollinearPoints(polygon.outer_boundary()));

  for(auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
  {
    newPolygon.add_hole(removeCollinearPoints(*hole));
  }

  return newPolygon;
}

Polygon MinkowskiSumCalculator::removeCollinearPoints(const Polygon& polygon)
{
  if (polygon.size() <= 3)
    return polygon;

  Polygon newPolygon;

  //Add non collinear vertices
  for(auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end(); ++vertex)
  {
    auto nextVertex = vertex + 1;
    if (nextVertex == polygon.vertices_end())
      continue;

    auto nextNextVertex = vertex + 2;
    if (nextNextVertex == polygon.vertices_end())
      continue;

    //Check for colinearity
    if (!CGAL::collinear(*vertex, *nextVertex, *nextNextVertex))
    {
      newPolygon.push_back(*nextVertex);
    }
  }

  //The two edge cases that needs to be considered: (last-1,last,first) and (last,first,first+1)
  auto lastMinusOne = polygon.vertex(polygon.size() - 2);
  auto last = polygon.vertex(polygon.size() - 1);
  auto first = polygon.vertex(0);
  auto firstPlusOne = polygon.vertex(1);

  if (!CGAL::collinear(lastMinusOne, last, first))
  {
    newPolygon.push_back(last);
  }

  if (!CGAL::collinear(last, first, firstPlusOne))
  {
    newPolygon.push_back(first);
  }

  return newPolygon;
}

PolygonWithHoles MinkowskiSumCalculator::scalePolygon(const PolygonWithHoles& polygon, Kernel::FT scale)
{
  PolygonWithHoles newPolygon;
  for(auto it = polygon.outer_boundary().vertices_begin(); it != polygon.outer_boundary().vertices_end(); ++it)
  {
    newPolygon.outer_boundary().push_back(Point(it->x()*scale, it->y()*scale));
  }

  for(auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
  {
    Polygon newHole;
    for(auto it = hole->vertices_begin(); it != hole->vertices_end(); ++it)
    {
      newHole.push_back(Point(it->x()*scale, it->y()*scale));
    }

    newPolygon.add_hole(newHole);
  }

  return newPolygon;	
}

MinkowskiSumResult::MinkowskiSumResult()
{
}

MinkowskiSumResult::MinkowskiSumResult(const std::list<PolygonWithHoles>& sum)
  : sum_(sum)
{
  //Decompose input and add to cache
  for(auto polygon = sum.begin(); polygon != sum.end(); ++polygon)
  {
    MinkowskiSumCalculator::decomposePolygon(*polygon, decomposedSum_);
  }
}

double MinkowskiSumResult::getArea() const
{
  //Calculate the area
  Kernel::FT area = 0;
  for(auto it = decomposedSum_.begin(); it != decomposedSum_.end(); ++it)
  {
    area += it->area();
  }

  return CGAL::to_double(area);
}

int MinkowskiSumResult::getSize() const
{
  //Calculate the size
  int numVertices = 0;
  std::list<Polygon> decomposedScaledSum;
  for(auto it = sum_.begin(); it != sum_.end(); ++it)
  {
    numVertices += it->outer_boundary().size();
    std::for_each(it->holes_begin(), it->holes_end(), [&numVertices](const Polygon& hole) { numVertices += hole.size(); });
  }

  return numVertices;
}

void MinkowskiSumResult::print() const
{
  //Print the sum to screen
  cout << "Polygons:" << endl;
  for(auto it = sum_.begin(); it != sum_.end(); ++it)
  {
    cout << *it << endl;
  }

  //Calculate the area and size
  cout << "Size: " << getSize() << endl;
  cout << "Area: " << getArea() << endl;
}
