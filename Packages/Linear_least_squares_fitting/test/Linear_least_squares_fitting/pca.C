// Test program for the linear_least_square_fitting() functions.
// Pierre Alliez

#include <vector>
#include <cassert>
#include <stdlib.h>

#include <CGAL/Cartesian.h>

#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Point_2           Point_2;
typedef K::Line_2            Line_2;
typedef K::Point_3           Point_3;


// dump a point set and a line to a postscript file
void dump_ps(char *pFilename,
             std::vector<Point_2>& points,
             Line_2& line)
{
  // try opening file
  FILE *pFile = fopen(pFilename,"wt");
  if(pFile == NULL)
  {
    std::cerr << "unable to open file " << pFilename << std::endl;
    return;
  }

  // header
  fprintf(pFile,"%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(pFile,"%%%%BoundingBox: 0 0 500 500\n");
  fprintf(pFile,"%%%%EndComments\n");
  fprintf(pFile,"gsave\n");

  // options
  fprintf(pFile,"1 setlinewidth\n");
  fprintf(pFile,"0 0 0 setrgbcolor\n");

  // stroke macro -> S
  fprintf(pFile,"\n%% stroke - x1 y1 x2 y2 S\n");
  fprintf(pFile,"/S {moveto lineto stroke} bind def\n\n");

  // disc macro -> D (dot)
  fprintf(pFile,"\n%% disc macro - x y D\n");
  fprintf(pFile,"/D {1 0 360 arc closepath fill} bind def\n");

  for(std::vector<Point_2>::iterator it = points.begin();
      it != points.end();
      it++)
  {
    const Point_2& p = *it;
    fprintf(pFile,"%9.3f %9.3f D\n",500.0f*p.x(),500.0f*p.y());
  }

  // output line segment
  Point_2 a = line.point(-100);
  Point_2 b = line.point(100);
  fprintf(pFile,"%9.3f %9.3f %9.3f %9.3f S\n",
          500.0f*a.x(),500.0f*a.y(),
          500.0f*b.x(),500.0f*b.y());

  // emit EPS trailer
  fputs("grestore\n\n",pFile);
  fputs("showpage\n",pFile);
  fclose(pFile);
}


Point_2 random_point_2()
{
  FT x = rand() / (FT)RAND_MAX;
  FT y = rand() / (FT)RAND_MAX;
  return Point_2(x,y);
}

void test_2(const unsigned int nb_points,
            const FT epsilon)
{
  std::cout << "2D: fit a line to a point set" << std::endl;

  // create random points nearby a segment
  std::vector<Point_2> points;
  Point_2 p = random_point_2();
  Point_2 q = random_point_2();
  std::cout << "  generate " << nb_points << 
       " 2D random points on a segment...";
  points_on_segment_2(p,q,nb_points,std::back_inserter(points));
  perturb_points_2(points.begin(),points.end(),epsilon);
  std::cout << "done" << std::endl;

  // fit a line
  std::cout << "  fit a 2D line...";
  Line_2 line;
  FT quality = linear_least_squares_fitting_2(points.begin(),points.end(),line);
  std::cout << "done (quality: " << quality << ")" << std::endl;

  // dump to ps
  std::cout << "  dump to ps...";
  dump_ps("test.ps",points,line);
  std::cout << "done" << std::endl;
}

int main()
{
  test_2(1000,0.05);
  return 0;
}
