#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include<CGAL/Polygon_2.h>
#include<CGAL/create_offset_polygons_2.h>
#include<CGAL/draw_straight_skeleton_2.h>

#include <memory>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_2                                              Point;

typedef CGAL::Polygon_2<K>                                      Polygon_2;
typedef std::shared_ptr<Polygon_2>                            PolygonPtr;

void low_precision_run()
{
  Polygon_2 poly;
  poly.push_back(Point(3.4641, 16.9915));
  poly.push_back(Point(3.4641, 15.4955));
  poly.push_back(Point(3.4641, 13.0989));
  poly.push_back(Point(3.4641, 9.94113));
  poly.push_back(Point(3.4641, 6.20559));
  poly.push_back(Point(3.4641, 2.10939));
  poly.push_back(Point(3.4641, -2.10939));
  poly.push_back(Point(3.4641, -6.20559));
  poly.push_back(Point(3.4641, -9.94113));
  poly.push_back(Point(3.4641, -13.0989));
  poly.push_back(Point(3.4641, -15.4955));
  poly.push_back(Point(3.4641, -16.9915));
  poly.push_back(Point(3.4641, -17.5));
  poly.push_back(Point(-6.52835e-17, -19.5));
  poly.push_back(Point(-6.33865e-17, -18.9334));
  poly.push_back(Point(-5.78057e-17, -17.2664));
  poly.push_back(Point(-4.88654e-17, -14.596));
  poly.push_back(Point(-3.70853e-17, -11.0773));
  poly.push_back(Point(-2.31499e-17, -6.9148));
  poly.push_back(Point(-7.86906e-18, -2.35047));
  poly.push_back(Point(7.86906e-18, 2.35047));
  poly.push_back(Point(2.31499e-17, 6.9148));
  poly.push_back(Point(3.70853e-17, 11.0773));
  poly.push_back(Point(4.88654e-17, 14.596));
  poly.push_back(Point(5.78057e-17, 17.2664));
  poly.push_back(Point(6.33865e-17, 18.9334));
  poly.push_back(Point(6.52835e-17, 19.5));
  poly.push_back(Point(3.4641, 17.5));

  assert(poly.is_simple());

  if(poly.is_clockwise_oriented())
    poly.reverse_orientation();

  std::vector<PolygonPtr> exteriorSkeleton =
    CGAL::create_exterior_skeleton_and_offset_polygons_2(1e-5, poly);

  assert(exteriorSkeleton.size() == 2);
  assert(exteriorSkeleton[0]->size() == 4);
  assert(exteriorSkeleton[0]->is_simple());
  assert(exteriorSkeleton[1]->is_simple());
}

void high_precision_run()
{
  Polygon_2 poly;
  poly.push_back(Point(3.4641015529632568, 16.991481781005859));
  poly.push_back(Point(3.4641015529632568, 15.495480537414551));
  poly.push_back(Point(3.4641015529632568, 13.09893798828125));
  poly.push_back(Point(3.4641015529632568, 9.9411334991455078));
  poly.push_back(Point(3.4641015529632568, 6.2055854797363281));
  poly.push_back(Point(3.4641015529632568, 2.1093919277191162));
  poly.push_back(Point(3.4641015529632568, -2.1093919277191162));
  poly.push_back(Point(3.4641015529632568, -6.2055854797363281));
  poly.push_back(Point(3.4641015529632568, -9.9411334991455078));
  poly.push_back(Point(3.4641015529632568, -13.09893798828125));
  poly.push_back(Point(3.4641015529632568, -15.495480537414551));
  poly.push_back(Point(3.4641015529632568, -16.991481781005859));
  poly.push_back(Point(3.4641015529632568, -17.5));
  poly.push_back(Point(-6.5283512761760263e-17, -19.5));
  poly.push_back(Point(-6.3386490615069335e-17, -18.933364868164063));
  poly.push_back(Point(-5.7805677252904074e-17, -17.266391754150391));
  poly.push_back(Point(-4.8865411228468343e-17, -14.595959663391113));
  poly.push_back(Point(-3.7085263204542273e-17, -11.077262878417969));
  poly.push_back(Point(-2.314985300804045e-17, -6.9147953987121582));
  poly.push_back(Point(-7.8690580132517397e-18, -2.3504652976989746));
  poly.push_back(Point(7.8690580132513576e-18, 2.3504652976989746));
  poly.push_back(Point(2.3149853008040067e-17, 6.9147953987121582));
  poly.push_back(Point(3.7085263204541891e-17, 11.077262878417969));
  poly.push_back(Point(4.8865411228467961e-17, 14.595959663391113));
  poly.push_back(Point(5.7805677252903704e-17, 17.266391754150391));
  poly.push_back(Point(6.3386490615068965e-17, 18.933364868164063));
  poly.push_back(Point(6.5283512761759894e-17, 19.5));
  poly.push_back(Point(3.4641015529632568, 17.5));

  assert(poly.is_simple());

  if(poly.is_clockwise_oriented())
    poly.reverse_orientation();

  std::vector<PolygonPtr> exteriorSkeleton =
    CGAL::create_exterior_skeleton_and_offset_polygons_2(1e-5, poly, K());

  assert(exteriorSkeleton.size() == 2);
  assert(exteriorSkeleton[0]->is_simple());
  assert(exteriorSkeleton[0]->size() == 4);
  assert(exteriorSkeleton[1]->is_simple());
}

int main()
{
  std::cout << "------------------- low precision test ---------------------" << std::endl;
  low_precision_run();

  std::cout << "------------------- high precision test ---------------------" << std::endl;
  high_precision_run();

  return EXIT_SUCCESS;
}
