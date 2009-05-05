#ifndef RENDER_POVRAY
#define RENDER_POVRAY

template <class Point>
struct Compare_lex_point {
  bool operator()(Point p1, Point p2) {
    if (p1.x() == p2.x()) {
      if (p1.y() == p2.y())
	return (p1.z() < p2.z());
      return (p1.y() < p2.y());
    }
    return (p1.x() < p2.x());
  }
};

template <class PDT, class ITYPE, class SegmentSet>
class Render_povray
{
  typedef typename PDT::Point Point;
  typedef typename PDT::Cell_handle                        Cell_handle;
  typedef typename PDT::Vertex_handle                      Vertex_handle;
  typedef typename PDT::Vertex_iterator                    Vertex_iterator;
  typedef typename PDT::Cell_iterator                      Cell_iterator;
  typedef typename PDT::Edge_iterator                      Edge_iterator;
  typedef typename PDT::Periodic_segment_iterator          Segment_iterator;

public:
  Render_povray() : it_type(0), ddomain(true), sset(), clipping(false) {}

  Render_povray(ITYPE _it_type, bool _ddomain, SegmentSet & _sset,
		bool _clipping) :
    it_type(_it_type), ddomain(_ddomain), sset(_sset),
    clipping(_clipping) {}

  ~Render_povray() {}

private:
  ITYPE it_type;
  bool ddomain;
  SegmentSet sset;
  bool clipping;

public:
  bool render(PDT& pdt,
              const char *pov_name,
	      const char *ini_name
	      )
  {
    printf("file: %s, ini: %s\n",pov_name, ini_name);
    // open files
    FILE* ini = fopen(ini_name,"wt");
    FILE* pov = fopen(pov_name,"wt");
    if(!ini || !pov)
      return false;
    
    // width & height
    GLint viewport[4];
    glGetIntegerv( GL_VIEWPORT, viewport );
    unsigned int w = viewport[2];
    unsigned int h = viewport[3];

    // write INI file
    fprintf(ini, "Input_File_Name = %s\n",pov_name);
    fprintf(ini, "Width  = %d\n", w);
    fprintf(ini, "Height = %d\n", h);
    fprintf(ini, "Antialias = on\n");
    fprintf(ini, "Verbose = off\n");
    fprintf(ini, "Display = on\n");

    // background color
	  GLfloat cc[4];
    glGetFloatv( GL_COLOR_CLEAR_VALUE, cc );
    fprintf(pov, "background\n{\n  color rgb <%f, %f, %f>\n}\n\n", cc[0], cc[1], cc[2]);

    // materials
    fprintf(pov, 
            "#declare CUBE_MATERIAL = texture {\n"
            "  pigment { rgbf <%f, %f, %f, %f> }\n"
            "  finish  { ambient 0.%02d  diffuse 0.%02d }\n"
            "}\n\n",0.9,0.925,0.99,0.7,30,60);

    fprintf(pov, 
            "#declare EDGE_MATERIAL = texture {\n"
            "  pigment { rgb <%f, %f, %f> }\n"
            "  finish  { ambient 0.%02d  diffuse 0.%02d phong 0.%02d phong_size %d }\n"
            "}\n\n",0.78,0.27,0.11,30,80,30,100);

    fprintf(pov, 
            "#declare SPHERE_MATERIAL = texture {\n"
            "  pigment { rgb <%f, %f, %f> }\n"
            "  finish  { ambient 0.%02d  diffuse 0.%02d phong 0.%02d phong_size %d }\n"
            "}\n\n",0.85,0.71,0.23,30,80,30,100);
    
    fprintf(pov, 
            "#declare SPHERE_B_MATERIAL = texture {\n"
            "  pigment { rgb <%f, %f, %f> }\n"
            "  finish  { ambient 0.%02d  diffuse 0.%02d phong 0.%02d phong_size %d }\n"
            "}\n\n",0.8,0.8,0.8,30,80,30,100);
    
    // camera settings
    {
		  // find near plane distance
		  GLdouble modelview[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
		  GLdouble projection[16];
		  glGetDoublev( GL_PROJECTION_MATRIX, projection );
      GLint viewport[4];
		  glGetIntegerv( GL_VIEWPORT, viewport );
		  GLdouble ox, oy, oz, wx, wy, wz;
		  wx = viewport[0] + 0.5 * viewport[2];
		  wy = viewport[1] + 0.5 * viewport[3];
		  wz = 0.0;
		  gluUnProject( wx, wy, wz, modelview, projection, viewport, &ox, &oy, &oz );
  		
		  double n            = -oz;
      double height       = h;
      double width        = (height * w) / h;
      double height_near  = 2.0 * n * tan(45.0*3.141592/360.0);
      double width_near   = height_near / height * width;
      fprintf(pov, 
	      "camera\n{\n"
	      "  location  <0, 0, 0>\n"
	      "  direction <0, 0, %f>\n"
	      "  up        <0, %f, 0>\n"
	      "  right     <%f, 0, 0>\n"
	      "}\n\n",
	      (float)n, (float)height_near, (float)width_near);
    }

	  // modelview matrix
    GLdouble modelview[16];
	  glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    #define MV(i,j) modelview[j*4+i]
    fprintf(pov, 
            "#declare MODELVIEW = transform {\n  matrix <\n"
            "    %f, %f, %f,\n"
            "    %f, %f, %f,\n"
            "    %f, %f, %f,\n"
            "    %f, %f, %f\n"
            "  >\n}\n\n",
            MV(0,0), MV(1,0), -MV(2,0),
            MV(0,1), MV(1,1), -MV(2,1),
            MV(0,2), MV(1,2), -MV(2,2),
            MV(0,3), MV(1,3), -MV(2,3));
    #undef MV

    // lights
    const float scene_radius = 5.;
    write_lights(pov,scene_radius);

    // write boundary mesh
    write_scene(pdt,pov);

    // close files
    fclose(pov);
    fclose(ini);
    return true;
  }

  void write_lights(FILE* pov,
                    const float scene_radius)
  {
    float r = scene_radius;
    Vector u(1.0,0.0,0.0);
    Vector v(0.0,0.0,1.0);
    Vector up(0.0,1.0,0.0);

    fprintf(pov, "global_settings { ambient_light rgb <1.0, 1.0, 1.0> }\n\n");

    fprintf(pov, 
					  "light_source\n{\n"
					  "  <%f, %f, %f>\n"
					  "  rgb <0.7, 0.7, 0.7> shadowless\n"
					  "}\n\n",
					  0.5*r, 0.5*r, -0.1*r);

    fprintf(pov, 
					  "light_source\n{\n"
					  "  <%f, %f, %f>\n"
					  "  rgb <0.7, 0.7, 0.7> shadowless\n"
					  "}\n\n",
					  -0.5*r, 0.5*r, -0.1*r);


	  // point light
   // if (dialog_->point_light_radiobutton->isChecked())
   // {
   //   fprintf(pov, 
   //           "light_source\n{\n"
   //           "  <%f, %f, %f>\n"
   //           "  rgb <1.0, 1.0, 1.0>\n"
   //           "}\n\n",
   //           10.0*r*up[0], 10.0*r*up[1], 10.0*r*-up[2] );
   // }

   // // area light
	  //else
    {
      fprintf(pov, 
              "light_source\n{\n"
              "  <%f, %f, %f>\n"
              "  rgb <1.0, 1.0, 1.0>\n"
              "  area_light <%f, %f, %f>, <%f, %f, %f>, 10, 10\n"
              "  adaptive 1\n"
              "  jitter\n"
              "}\n\n",
              4.0*r, 2.0*r, -2.0*r,
              0.0,3.5,3.5,
              1.5,-1.5,1.5 );
    }
  }

  void write_scene(PDT& pdt, FILE* pov)
  {
    // initializations
    fprintf(pov, "#declare sph_radius=0.005;\n");
    fprintf(pov, "#declare cyl_radius=0.0002;\n\n");

    fprintf(pov, "#declare offsetx=0;\n");
    fprintf(pov, "#declare offsety=0;\n");
    fprintf(pov, "#declare offsetz=0;\n");

    // write vertex positions
    //int nbv = pdt.number_of_vertices();
    //fprintf(pov, "  vertex_vectors\n  {\n    %d,\n",nbv);
    for(Point_iterator vit = pdt.periodic_points_begin(PDT::UNIQUE);
	vit != pdt.periodic_points_end(PDT::UNIQUE); vit++)
      {
	fprintf(pov, "sphere \n");
	fprintf(pov, "{\n");
	fprintf(pov, "    < %lf+offsetx, %lf+offsety, %lf+offsetz >, sph_radius\n",
		vit->first.x(),vit->first.y(),vit->first.z());
	fprintf(pov, "  texture { SPHERE_MATERIAL }\n");
	fprintf(pov, "  transform MODELVIEW\n");
	fprintf(pov, "}\n\n");
      }
    //fprintf(pov, "  }\n\n");
    

    // write edge
    //    fprintf(pov, "  edges\n  {\n    %d,\n",nbf);
    std::set<Point, Compare_lex_point<Point> > ptset;
    typename SegmentSet::iterator he;
    for(he = sset.begin(); he != sset.end(); he++)
    {
      Point s = he->source();
      Point t = he->target();
      if (   s.x() >= 1.0 || s.x() < 0.0
	  || s.y() >= 1.0 || s.y() < 0.0
	  || s.z() >= 1.0 || s.z() < 0.0 )
	ptset.insert(s);
      if (   t.x() >= 1.0 || t.x() < 0.0
	  || t.y() >= 1.0 || t.y() < 0.0
	  || t.z() >= 1.0 || t.z() < 0.0 )
	ptset.insert(t);
      fprintf(pov, "cylinder \n");
      fprintf(pov, "{\n");
      fprintf(pov, "    < %lf+offsetx, %lf+offsety, %lf+offsetz >, < %lf+offsetx, %lf+offsety, %lf+offsetz >, cyl_radius\n",
	      s[0],s[1],s[2],t[0],t[1],t[2]);
      if (clipping) {
      fprintf(pov, "  clipped_by { box { < 0 0 0 >, < 1 1 1 > } }\n");
      }
      fprintf(pov, "  texture { EDGE_MATERIAL }\n");
      fprintf(pov, "  transform MODELVIEW\n");
      fprintf(pov, "}\n\n");
    }

    typename std::set<Point, Compare_lex_point<Point> >::iterator pit;
    for(pit = ptset.begin(); pit != ptset.end(); pit++)
    {
      fprintf(pov, "sphere \n");
      fprintf(pov, "{\n");
      fprintf(pov, "    < %lf, %lf, %lf >, 0.003\n",
	      pit->x(),pit->y(),pit->z());
      fprintf(pov, "  texture { SPHERE_B_MATERIAL }\n");
      fprintf(pov, "  transform MODELVIEW\n");
      fprintf(pov, "}\n\n");
    }


    // write cube
    if (ddomain) {
      std::vector<const char *> cube_edges;
      cube_edges.push_back("    < 0., 0., 0. >, < 1., 0., 0. >, 0.01\n");
      cube_edges.push_back("    < 0., 0., 0. >, < 0., 1., 0. >, 0.01\n");
      cube_edges.push_back("    < 1., 0., 0. >, < 1., 1., 0. >, 0.01\n");
      cube_edges.push_back("    < 0., 1., 0. >, < 1., 1., 0. >, 0.01\n");
      cube_edges.push_back("    < 0., 0., 0. >, < 0., 0., 1. >, 0.01\n");
      cube_edges.push_back("    < 0., 1., 0. >, < 0., 1., 1. >, 0.01\n");
      cube_edges.push_back("    < 1., 0., 0. >, < 1., 0., 1. >, 0.01\n");
      cube_edges.push_back("    < 1., 1., 0. >, < 1., 1., 1. >, 0.01\n");
      cube_edges.push_back("    < 0., 0., 1. >, < 1., 0., 1. >, 0.01\n");
      cube_edges.push_back("    < 0., 0., 1. >, < 0., 1., 1. >, 0.01\n");
      cube_edges.push_back("    < 1., 0., 1. >, < 1., 1., 1. >, 0.01\n");
      cube_edges.push_back("    < 0., 1., 1. >, < 1., 1., 1. >, 0.01\n");
   
      for (unsigned int i=0 ; i<cube_edges.size() ; i++) {
	fprintf(pov, "cylinder \n");
	fprintf(pov, "{\n");
	fprintf(pov, cube_edges[i]);
	fprintf(pov, "  texture { CUBE_MATERIAL }\n");
	fprintf(pov, "  transform MODELVIEW\n");
	fprintf(pov, "}\n\n");
      }

      std::vector<const char *> cube_vertices;
      cube_vertices.push_back("    < 0., 0., 0. >, 0.01\n");
      cube_vertices.push_back("    < 0., 0., 1. >, 0.01\n");
      cube_vertices.push_back("    < 0., 1., 0. >, 0.01\n");
      cube_vertices.push_back("    < 0., 1., 1. >, 0.01\n");
      cube_vertices.push_back("    < 1., 0., 0. >, 0.01\n");
      cube_vertices.push_back("    < 1., 0., 1. >, 0.01\n");
      cube_vertices.push_back("    < 1., 1., 0. >, 0.01\n");
      cube_vertices.push_back("    < 1., 1., 1. >, 0.01\n");

      for (unsigned int i=0 ; i<cube_vertices.size() ; i++) {
	fprintf(pov, "sphere \n");
	fprintf(pov, "{\n");
	fprintf(pov, cube_vertices[i]);
	fprintf(pov, "  texture { CUBE_MATERIAL }\n");
	fprintf(pov, "  transform MODELVIEW\n");
	fprintf(pov, "}\n\n");
      }
    }
    // material & modelview
    //    fprintf(pov, "  texture { MESH_MATERIAL }\n\n");
    //fprintf(pov, "  transform MODELVIEW\n\n");

  }

}; // end class Render_povray

#endif // RENDER_POVRAY
