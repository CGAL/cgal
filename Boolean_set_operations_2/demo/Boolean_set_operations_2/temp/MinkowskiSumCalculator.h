#ifndef CGAL_MINKOWSKISUMCALCULATOR_H
#define CGAL_MINKOWSKISUMCALCULATOR_H

#include <list>
#include <map>
//#include <memory>
#include "Typedefs.h"

class MinkowskiSumResult
{
public:
  MinkowskiSumResult();
  MinkowskiSumResult(const std::list<PolygonWithHoles>& sum);

  const std::list<PolygonWithHoles>& getSum() const { return sum_; }
  const std::list<Polygon>& getDecomposedSum() const { return decomposedSum_; }
  double getArea() const;
  int getSize() const;
  void print() const;

private:
  std::list<PolygonWithHoles> sum_;
  std::list<Polygon> decomposedSum_;
};

class MinkowskiSumCalculator
{
public:
  MinkowskiSumCalculator();
  ~MinkowskiSumCalculator();

  void setInput(const std::list<PolygonWithHoles>& input);
  std::shared_ptr<const MinkowskiSumResult> getSum(int k);

  static void decomposePolygon(const PolygonWithHoles& polygon, std::list<Polygon>& result);

private:
  PolygonWithHoles scalePolygon(const PolygonWithHoles& polygon, Kernel::FT scale);
  PolygonWithHoles removeCollinearPoints(const PolygonWithHoles& polygon);
  Polygon removeCollinearPoints(const Polygon& polygon);

  std::map<int, std::shared_ptr<MinkowskiSumResult>> cache_;
  std::list<PolygonWithHoles> input_;
};

#endif // CGAL_MINKOWSKISUMCALCULATOR_H
