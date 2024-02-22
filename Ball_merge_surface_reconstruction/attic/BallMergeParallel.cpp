#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/IO/Color.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <string>
#include <sys/time.h>
//#include <bits/stdc++.h>
#include "happly.h"
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
std::ofstream myfile;
clock_t start, _end;
int group = 1, gcount = 0, option, maxg = 0, maingroup, max1 = 0, secondgroup;
float minx = 9999, miny = 9999, minz = 9999, maxx = -9999, maxy = -9999, maxz = -9999;
double par, bbdiaglen, tlen;
Delaunay T;
double sq(double a) { return a * a; } //Function to compute the square
double distance(Point p1, Point p2) { return sqrt(sq(p2.x() - p1.x()) + sq(p2.y() - p1.y()) + sq(p2.z() - p1.z())); } //To compute Euclidean distance between two 3D points
Point cc(Delaunay::Cell_handle v){ //Function to compute the circumcenter of a Cell
    CGAL::Tetrahedron_3<K> t1(v->vertex(0)->point(), v->vertex(1)->point(), v->vertex(2)->point(), v->vertex(3)->point());
    return CGAL::circumcenter(t1);
}
int _function(Delaunay::Cell_handle fh1, Delaunay::Cell_handle fh2){ //Function that computes the Intersection Ratio (IR) given two cell handles
    Point p1 = cc(fh1); //p1 is the circumcenter of the first cell
    Point p2 = cc(fh2); //p2 is the circumcenter of the second cell 
    double d = distance(p1, p2); //distance between the circumcenters
    double a = distance(fh1->vertex(0)->point(), p1);//circumradius of circumsphere of the first cell
    double b = distance(fh2->vertex(0)->point(), p2);//circumradius of circumsphere of the second cell
    if ((a + b - d) / a > par || (a + b - d) / b > par)//Check whether the IR is less than a user given parameter
        return 1;
    return 0;
}
void Group(Delaunay &T, Delaunay::Cell_handle fh){ //Function to recursively group all the cells that is mergable from the cell fh with user specified parameter
    std::queue<Delaunay::Cell_handle> q; //A queue to do the grouping
    q.push(fh); //The cell itself is mergable
    while (!q.empty()){ //A traversal
        fh = q.front();
        q.pop();
        fh->info() = group; //A label
        gcount++;
        for (int i = 0; i < 4; i++) //For each neighbors
            if ((!T.is_infinite(fh->neighbor(i)))&&(fh->neighbor(i)->info() == 0 && _function(fh, fh->neighbor(i)) == 1)){//If it is unlabeled and can be merged - using the function that computes IR
                fh->neighbor(i)->info() = group;//If mergable, then change the label
                q.push(fh->neighbor(i));//Push it into the traversal list
            }
    }
}
int bblen(Point a, Point b, Point c){//Function to check whether a triangle made by three points has an edge greater than bbdiag/(user defined paramter - 200 by default) - used for filtering long triangles in local ballmerge
    if (distance(a, b) > bbdiaglen / tlen || distance(c, b) > bbdiaglen / tlen || distance(a, c) > bbdiaglen / tlen)
        return 0;
    return 1;
}
int main(int argc, char **argv){
    const char *fname = argv[1];//Filename
    const char *inFilename = argv[1];//Filename
    std::ifstream inStream(inFilename);//Read the file
    std::vector<std::array<unsigned char, 3>> meshVertexColors;//Variable to hold if the input has textures
    Point p;
    par = atof(argv[2]);//Parameter to check IR
    option = atoi(argv[3]);//Option to select global and local variants - 1 for global and 0 for local
    if (argc >= 5)//If local, there is an option to give extra optional parameter to filter long triangles
        tlen = atof(argv[4]);
    else
        tlen = 200.0;//If not provided, default 200 will be taken
    float minx = 9999, miny = 9999, minz = 9999, maxx = -9999, maxy = -9999, maxz = -9999;
    double min[3], max[3];
    int vIndex = 0;
    std::vector<std::pair<Point, unsigned>> points;
    for (int i = 0; i < 3; i++){
        min[i] = std::numeric_limits<double>::max();
        max[i] = -std::numeric_limits<double>::max();
    }
    while (!inStream.eof()){//Read the data from the input file - coordinates, texture information (if available)
        string line;
        getline(inStream, line);
        istringstream str(line);
        if (str >> p){
            int r, g, b;
            if (str >> r){
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
            for (int i = 0; i < 3; i++){
                if (p[i] < min[i])
                    min[i] = p[i];
                if (p[i] > max[i])
                    max[i] = p[i];
            }
        }
    }
#ifdef CGAL_LINKED_WITH_TBB //Parallel Delaunay computation
   Bbox_3 bbox = Bbox_3(min[0], min[1], min[2], max[0], max[1], max[2]);
   Delaunay::Lock_data_structure locking_ds(bbox, 50);
    Delaunay T(points.begin(), points.end(), &locking_ds);
 #else
    Delaunay T(points.begin(), points.end());//General Delaunay computation
 #endif
    Delaunay::Finite_cells_iterator vit;
    for (vit = T.finite_cells_begin(); vit != T.finite_cells_end(); vit++)//Initialize the labels of all tetrahedrons
        vit->info() = 0;
    if (option == 1){//If the user opted for global algorithm
        for (vit = T.finite_cells_begin(); vit != T.finite_cells_end(); vit++){//For each cell
            if (vit->info() == 0){//If the cell label is not altered
                vit->info() = group;//Assign a label
                Group(T,vit);//Group all the mergable cells
                if (maxg < gcount){//Remember the largest group - based on the number of tetrahedrons in the group
                    max1 = maxg;
                    secondgroup = maingroup;//To remember the second largest group
                    maxg = gcount;
                    maingroup = group;
                }
                else if (max1 < gcount && gcount != maxg){
                    max1 = gcount;
                    secondgroup = group;//To remember the second largest group
                }
                gcount = 0;
                group++;//Update the label for next cell
            }
        }
    }
    std::vector<std::array<double, 3>> meshVertexPositions(points.size());//Preparing the tetrahedrons for creating the PLY file & PLY file writing starts
    std::vector<std::vector<int>> meshFaceIndices;
    for (Delaunay::Finite_vertices_iterator vIter = T.finite_vertices_begin(); vIter != T.finite_vertices_end(); vIter++){
        Point point = vIter->point();
        int vIndex = vIter->info();
        meshVertexPositions[vIndex][0] = point.x();
        meshVertexPositions[vIndex][1] = point.y();
        meshVertexPositions[vIndex][2] = point.z();
    }
    std::string st = argv[1];
    st = st + std::to_string(par) + "out"+std::to_string(tlen)+".ply";
    char outname[100];
    strcpy(outname, st.c_str());
    bbdiaglen = distance(Point(minx, miny, minz), Point(maxx, maxy, maxz));
    for (vit = T.finite_cells_begin(); vit != T.finite_cells_end(); vit++)
        for (int i = 0; i < 4; i++)
            if (option == 1){//If global, write the cell details of the largest group to the PLY file
                if (vit->info() == maingroup && (vit->neighbor(i)->info() != maingroup || T.is_infinite(vit->neighbor(i)))){//Write the triangles between cells if they have have different labels and one of them is labeled as the same as the largest group
                    std::vector<int> indices(3);
                    for (int j = 0; j < 3; j++)
                        indices[j] = vit->vertex((i + 1 + j) % 4)->info();
                    meshFaceIndices.push_back(indices);
                }
            }
            else if (vit->neighbor(i)->info() != 9999)//If local
                if (bblen(vit->vertex((i + 1) % 4)->point(), vit->vertex((i + 2) % 4)->point(), vit->vertex((i + 3) % 4)->point()))//If the triangle crosses our bbdiagonal based criteria
                    if (!_function(vit, vit->neighbor(i)) == 1||T.is_infinite(vit->neighbor(i))){//If the cells cannot be merged, then wirte the triangle between these two cells to the PLY file
                        vit->info() = 9999;
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
    meshFaceIndices.clear();
    if (option == 1){//Sometimes, in the gloabl case, the largest group would be a mould created by the convex hull, just to avoid it, we will write the second largest group as well
        st = st + "1.ply";
        strcpy(outname, st.c_str());
        for (vit = T.finite_cells_begin(); vit != T.finite_cells_end(); vit++)
            for (int i = 0; i < 4; i++)
                if (vit->info() == secondgroup && (vit->neighbor(i)->info() != secondgroup || T.is_infinite(vit->neighbor(i)))){//Write the triangles between cells if they have have different labels and one of them is labeled as the same as the second largest group
                    std::vector<int> indices(3);
                    for (int j = 0; j < 3; j++)
                        indices[j] = vit->vertex((i + 1 + j) % 4)->info();
                    meshFaceIndices.push_back(indices);
                }
        happly::PLYData plyOut2;
        plyOut2.addVertexPositions(meshVertexPositions);
        if (meshVertexColors.size() > 0)
            plyOut2.addVertexColors(meshVertexColors);
        plyOut2.addFaceIndices(meshFaceIndices);
        plyOut2.write(outname, happly::DataFormat::Binary);
    }
}
