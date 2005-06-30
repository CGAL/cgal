// utils for least_squares_linear_fitting
// Pierre Alliez

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
  fprintf(pFile,"1 0 0 setrgbcolor\n");
  fprintf(pFile,"%9.3f %9.3f %9.3f %9.3f S\n",
          500.0f*a.x(),500.0f*a.y(),
          500.0f*b.x(),500.0f*b.y());

  // emit EPS trailer
  fputs("grestore\n\n",pFile);
  fputs("showpage\n",pFile);
  fclose(pFile);
}

// dump a 2D triangle set and a line to a postscript file
void dump_ps(char *pFilename,
             std::vector<Triangle_2>& triangles,
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

  // triangle macro -> T
  fprintf(pFile,"\n%% triangle macro - x1 y1 x2 y2 x3 y3 T\n");
  fprintf(pFile,"/T {moveto lineto lineto closepath stroke} bind def\n");

  for(std::vector<Triangle_2>::iterator it = triangles.begin();
      it != triangles.end();
      it++)
  {
    const Triangle_2& t = *it;
    const Point_2& a = t[0];
    const Point_2& b = t[1];
    const Point_2& c = t[2];
    fprintf(pFile,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f T\n",
            500.0f*a.x(),500.0f*a.y(),
            500.0f*b.x(),500.0f*b.y(),
	    500.0f*c.x(),500.0f*c.y());
  }

  // output line segment
  Point_2 a = line.point(-100);
  Point_2 b = line.point(100);
  fprintf(pFile,"1 0 0 setrgbcolor\n");
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

