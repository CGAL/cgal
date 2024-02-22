#include "happly.h"
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <sys/time.h>
using namespace CGAL;
using namespace std;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<int, K> Cb;
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Triangulation_data_structure_3<Vb, Cb, Parallel_tag> Tds;
#else
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
#endif
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
struct node {
  float priority;
  Delaunay::Cell_handle ch;
  struct node *link;
} tempn;
bool Comparen(node n1, node n2) { return n1.priority < n2.priority; }
std::priority_queue<node, std::vector<node>, std::function<bool(node, node)>> pq(Comparen);
std::ofstream myfile;
clock_t start, _end;
int group = 1, gcount = 0, option, maxg = 0, maingroup, max1 = 0, secondgroup;
float minx = 9999, miny = 9999, minz = 9999, maxx = -9999, maxy = -9999, maxz = -9999;
double par, bbdiaglen, tlen;
Delaunay T;
double sq(double a) { return a * a; }
Point cc(Delaunay::Cell_handle v) {
  CGAL::Tetrahedron_3<K> t1(v->vertex(0)->point(), v->vertex(1)->point(),v->vertex(2)->point(), v->vertex(3)->point());
  return CGAL::circumcenter(t1);
}
double distance(Point p1, Point p2) {
  return sqrt(sq(p2.x() - p1.x()) + sq(p2.y() - p1.y()) + sq(p2.z() - p1.z()));
}
double ir(Delaunay::Cell_handle fh1, Delaunay::Cell_handle fh2) {
  Point p1 = cc(fh1);
  Point p2 = cc(fh2);
  double d = distance(p1, p2);
  double a = distance(fh1->vertex(0)->point(), p1);
  double b = distance(fh2->vertex(0)->point(), p2);
  if (std::isnan(d))
    return 0;
  if ((a + b - d) / a > (a + b - d) / b)
    return (a + b - d) / a;
  return (a + b - d) / b;
}

int main(int argc, char **argv) {
  const char *fname = argv[1];
  const char *inFilename = argv[1];
  std::ifstream inStream(inFilename);
  std::vector<std::array<unsigned char, 3>> meshVertexColors;
  Point p;
  float minx = 9999, miny = 9999, minz = 9999, maxx = -9999, maxy = -9999, maxz = -9999;
  double min[3], max[3];
  int vIndex = 0;
  std::vector<std::pair<Point, unsigned>> points;
  for (int i = 0; i < 3; i++) {
    min[i] = std::numeric_limits<double>::max();
    max[i] = -std::numeric_limits<double>::max();
  }
  std::cout<<"Reading the file\n";
  while (!inStream.eof()) {
    string line;
    getline(inStream, line);
    istringstream str(line);
    if (str >> p) {
      int r, g, b;
      if (str >> r) {
        str >> g;
        str >> b;
        meshVertexColors.push_back({(unsigned char)r, (unsigned char)g, (unsigned char)b});
      }
      if (p.x() < minx)
        minx = p.x();
      if (p.x() < miny)
        miny = p.y();
      if (p.x() < minz)
        minz = p.z();
      if (p.x() > maxx)
        maxx = p.x();
      if (p.x() > maxy)
        maxy = p.y();
      if (p.x() > maxz)
        maxz = p.z();
      points.push_back(std::make_pair(Point(p.x(), p.y(), p.z()), vIndex++));
      for (int i = 0; i < 3; i++) {
        if (p[i] < min[i])
          min[i] = p[i];

        if (p[i] > max[i])
          max[i] = p[i];
      }
    }
  }
  std::cout<<"Computing Delaunay\n";
#ifdef CGAL_LINKED_WITH_TBB
  Bbox_3 bbox = Bbox_3(min[0], min[1], min[2], max[0], max[1], max[2]);
  Delaunay::Lock_data_structure locking_ds(bbox, 50);
  Delaunay T(points.begin(), points.end(), &locking_ds);
#else
  Delaunay T(points.begin(), points.end());
#endif
std::cout<<"Running BallMerge\n";
  double maxr = -9999;
  Delaunay::Cell_handle maxcell;
  Delaunay::Finite_cells_iterator vit;
  for (vit = T.finite_cells_begin(); vit != T.finite_cells_end(); vit++) {//Loop to pick the seed tetrahedron to grow from
    CGAL::Tetrahedron_3<K> t1(vit->vertex(0)->point(), vit->vertex(1)->point(),vit->vertex(2)->point(), vit->vertex(3)->point());
    Point tp1 = CGAL::circumcenter(t1);
    double dis1 = distance(tp1, vit->vertex(0)->point());
    if (dis1 > maxr) {
      Delaunay::Cell_handle ch = T.locate(tp1);
      Delaunay::Cell_handle ch1 = vit;
      if (!T.is_infinite(ch)  && (T.is_infinite(vit->neighbor(0)) || T.is_infinite(vit->neighbor(2)) || T.is_infinite(vit->neighbor(1)) || T.is_infinite(vit->neighbor(3)))){//Tetrahedron with largest circumradius which is part of the convex hull
        maxcell = ch;
        maxr = dis1;
      }
    }
    vit->info() = 0;
  }
  maxcell->info() = 1;
  for (int i = 0; i < 4; i++)
    if (!T.is_infinite(maxcell->neighbor(i))) {
      tempn.priority = ir(maxcell, maxcell->neighbor(i));
      tempn.ch = maxcell->neighbor(i);
      pq.push(tempn);
    }
  double maxir = ir(maxcell, maxcell->neighbor(0));
  for (int i = 0; i < 4; i++)
    if (maxir < ir(maxcell, maxcell->neighbor(i)))
      maxir = ir(maxcell, maxcell->neighbor(i));
  int loop = 0;
  maxir = 2.0;
  while (maxir > 1.0) {//Loop for 200 iterations maximum since maxir decreases in 0.01 steps
    loop++;
    while (!pq.empty()) {//Create groups from the seed tetrahedron for a particular parameter
      node tch = pq.top();
      pq.pop();
      if (tch.ch->info() == 0) {
        if (tch.priority > maxir) {
          tch.ch->info() = loop;
          for (int i = 0; i < 4; i++)
            if ((!T.is_infinite(tch.ch->neighbor(i))) && (tch.ch->neighbor(i)->info() == 0)) {          
              tempn.priority = ir(tch.ch, tch.ch->neighbor(i));
              tempn.ch = tch.ch->neighbor(i);
              pq.push(tempn);
            }
        } else {
          pq.push(tch);
          break;
        }
      }
    }
    int good = 0, masked = 0, separate = 0;
    for (Delaunay::Finite_vertices_iterator vIter = T.finite_vertices_begin();
         vIter != T.finite_vertices_end(); vIter++) {
      std::vector<Delaunay::Cell_handle> vif;
      T.incident_cells(vIter, back_inserter(vif));
      int convfl = 0, naddedfl = 0, addedfl = 0;
      for (int i = 0; i < vif.size(); i++) {
        if (T.is_infinite(vif.at(i)))
          convfl = 1;
        else if (vif.at(i)->info() == 0)
          naddedfl = 1;
        else
          addedfl = 1;
      }
      if (addedfl == 1 && (naddedfl == 1 || convfl == 1))
        good++;
      if (addedfl == 1 && (naddedfl == 0 && convfl == 0))
        masked++;
      if (addedfl == 0)
        separate++;
    }
    if (good < points.size() * 0.9999999 && masked < points.size() * 0.05)//Check how many points got removed compared to the result generated using the previous parameter
      maxir -= 0.01;//Step size for parameter change
    else
      break;
  }
  std::cout << "Automatic Parameter is " << maxir + 0.01 << "\n";
  std::cout<<"Writing to the file\n";
  std::vector<std::array<double, 3>> meshVertexPositions(points.size());
  std::vector<std::vector<int>> meshFaceIndices;
  for (Delaunay::Finite_vertices_iterator vIter = T.finite_vertices_begin();
       vIter != T.finite_vertices_end(); vIter++) {
    Point point = vIter->point();
    int vIndex = vIter->info();
    meshVertexPositions[vIndex][0] = point.x();
    meshVertexPositions[vIndex][1] = point.y();
    meshVertexPositions[vIndex][2] = point.z();
  }
  std::string st = argv[1];
  st = st + "out.ply";
  char outname[100];
  strcpy(outname, st.c_str());
  bbdiaglen = distance(Point(minx, miny, minz), Point(maxx, maxy, maxz));
  for (vit = T.finite_cells_begin(); vit != T.finite_cells_end(); vit++)
    for (int i = 0; i < 4; i++)
      if ((vit->info() != 0) && (vit->info() != loop) &&  ((vit->neighbor(i)->info() == 0 || vit->neighbor(i)->info() == loop) || T.is_infinite(vit->neighbor(i)))) {
        std::vector<int> indices(3);
        for (int j = 0; j < 3; j++)
          indices[j] = vit->vertex((i + 1 + j) % 4)->info();
        meshFaceIndices.push_back(indices);
      }
  happly::PLYData plyOut;
  plyOut.addVertexPositions(meshVertexPositions);
  if (meshVertexColors.size() > 0)
    plyOut.addVertexColors(meshVertexColors);
  plyOut.addFaceIndices(meshFaceIndices);
  plyOut.write(outname, happly::DataFormat::Binary);
}
