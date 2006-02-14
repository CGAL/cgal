/***************************************************************************
dump_eps.h  -  postscript output for CGAL polyhedron
----------------------------------------------------------------------------
begin                : june 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/
#ifndef DUMP_EPS
#define DUMP_EPS

#include "Enriched_polyhedron.h"
#include <vector>
#include <GL/glu.h>

template <class Polyhedron,class kernel>
class Dumper_eps
{
	typedef typename kernel::Point_3 Point;
	typedef typename kernel::Vector_3 Vector;
	typedef typename Polyhedron::Vertex                             Vertex;
	typedef typename Polyhedron::Vertex_iterator                    Vertex_iterator;
	typedef typename Polyhedron::Halfedge_handle                    Halfedge_handle;
	typedef typename Polyhedron::Edge_iterator                      Edge_iterator;
	typedef typename Polyhedron::Facet_iterator                     Facet_iterator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator  HV_circulator;
	typedef typename Polyhedron::Halfedge_around_facet_circulator   HF_circulator;

	class Projected_facet : public std::vector<double>
	{
		double m_z; // z-key
		bool m_is_valid;
	public:

		Projected_facet()
		{
			m_is_valid = true;
			m_z = 0.0;
		}

		const double z() const { return m_z; }
		void z(double z) { m_z = z; }
		bool is_valid() { return m_is_valid; }
		void invalid() { m_is_valid = false; }

		static bool closer(const Projected_facet& facet1,
											const Projected_facet& facet2)
		{
			return facet1.z() < facet2.z();
		}
	};

public:
	// life cycle
	Dumper_eps()	 {}
	~Dumper_eps() {}

	// dump	
  void dump(Polyhedron *pMesh,
		        char *pFileName,
		        double *modelMatrix,
						double *projMatrix,
            int *viewport,
						int FlagBack)
	{
	  CGAL_assertion(pMesh != NULL);
		double min_z = 1e308;
		double max_z = -1e308;
		// fill array
		std::vector<Projected_facet> pfacets;
		Facet_iterator pFacet;
		for(pFacet = pMesh->facets_begin();
				pFacet != pMesh->facets_end();
				pFacet++)
		{
			// add one facet
			Projected_facet projected_facet;

			// project center	
			Point center;
			pMesh->compute_facet_center(pFacet,center);

			// project center
			double px,py,pz;
			gluProject(center.x(),
								center.y(),
								center.z(),
								modelMatrix,
								projMatrix,
								viewport,
								&px,&py,&pz);
			projected_facet.z(-pz);

			if(-pz < min_z) min_z = -pz;
			if(-pz > max_z) max_z = -pz;

			// and vertices (xyxyxyxy)
 			HF_circulator h = pFacet->facet_begin();
			do
			{
				double x = h->vertex()->point().x();
				double y = h->vertex()->point().y();
				double z = h->vertex()->point().z();
    		gluProject(x,y,z,
    								modelMatrix,
    								projMatrix,
    								viewport,
    								&px,&py,&pz);
    		projected_facet.push_back(px);
    		projected_facet.push_back(py);
			}
			while(++h != pFacet->facet_begin());

			// add projected facet to sort
			pfacets.push_back(projected_facet);
		}
		// double range_z = max_z - min_z;

		// sort
		std::sort(pfacets.begin(),pfacets.end(),Projected_facet::closer);

		// crop
		int nbf = pfacets.size();
		int i;
		for(i=0;i<nbf;i++)
		{
			Projected_facet& f = pfacets[i];
			int degree = f.size()/2;
			for(int j=0;j<degree;j++)
			{
				if(f[2*j]   < viewport[0] || // xmin
					f[2*j+1] < viewport[1] || // ymin
					f[2*j]   > viewport[2] || // xmax
					f[2*j+1] > viewport[3])   // ymax
					f.invalid();
			}
		}
		
		// measure bounding box
		double minx = 1e308;
		double miny = 1e308;
		double maxx = -1e308;
		double maxy = -1e308;
		for(i=0;i<nbf;i++)
		{
			Projected_facet& f = pfacets[i];
			if(!f.is_valid())
				continue;

  		int degree = f.size()/2;
			for(int j=0;j<degree;j++)
  		{
				minx = std::min(minx,f[2*j]);
  			miny = std::min(miny,f[2*j+1]);
  			maxx = std::max(maxx,f[2*j]);
  			maxy = std::max(maxy,f[2*j+1]);
  		}
		}
		double rangex = maxx-minx;
		double rangey = maxy-miny;
		bool max_range_in_x = (rangex > rangey);
		double ratio = rangex/rangey;
		double range = std::max(rangex,rangey);
		double tmp = 500/range;
		double xoffest = 0.0;
		double yoffest = 0.0;
		
		double bbox_x = 500;
		if(!max_range_in_x)
			bbox_x = 500*ratio;
		double bbox_y = 500;
		if(max_range_in_x)
			bbox_y = 500/ratio;

		FILE *pFile = fopen(pFileName,"wt");
		CGAL_assertion(pFile != NULL);

		// Emit EPS header
		fprintf(pFile,"%%!PS-Adobe-2.0 EPSF-2.0\n");
		fprintf(pFile,"%%%%BoundingBox: 0 0 %g %g\n",bbox_x,bbox_y);
		fprintf(pFile,"%%%%EndComments\n");
		fprintf(pFile,"gsave\n");
		fprintf(pFile, "\n%g setlinewidth\n", 0.1);

		// color macros
		fprintf(pFile,"\n%% RGB color command - r g b C\n");
		fprintf(pFile,"/C { setrgbcolor } bind def\n");
		fprintf(pFile,"/white { 1 1 1 C } bind def\n");
		fprintf(pFile,"/black { 0 0 0 C } bind def\n");
		fprintf(pFile,"/min_width { 0.1 } bind def\n");
		fprintf(pFile,"/range_width { 0.3 } bind def\n");

		// triangle
		fprintf(pFile,"\n%% Outlined triangle - x3 y3 x2 y2 x1 y1 T\n");
		fprintf(pFile,"/T {6 copy white newpath moveto lineto lineto closepath fill \n \
		black newpath moveto lineto lineto closepath stroke} \n \
		bind def\n");

		fprintf(pFile,"\n%% Quad - x4 y4 x3 y3 x2 y2 x1 y1 P4\n");
		fprintf(pFile,"/P4 {8 copy white newpath moveto lineto lineto lineto closepath fill \n \
			black newpath moveto lineto lineto lineto closepath stroke} \n \
			bind def\n");

		fprintf(pFile,"\n%% Pentagon - x5 y5 x4 y4 x3 y3 x2 y2 x1 y1 P5\n");
		fprintf(pFile,"/P5 {10 copy white newpath moveto lineto lineto lineto lineto closepath fill \n \
			black newpath moveto lineto lineto lineto lineto closepath stroke} \n \
			bind def\n");

		fprintf(pFile,"\n%% Hexagon - x6 y6 x5 y5 x4 y4 x3 y3 x2 y2 x1 y1 P6\n");
		fprintf(pFile,"/P6 {12 copy white newpath moveto lineto lineto lineto lineto lineto closepath fill \n \
			black newpath moveto lineto lineto lineto lineto lineto closepath stroke} \n \
			bind def\n");

		for(i=0;i<nbf;i++)
		{
			Projected_facet& f = pfacets[i];
			if(!f.is_valid())
				continue;
			int degree = f.size()/2;
			switch(degree)
			{
				case 3:
				{
					double x1 = xoffest+(f[0]-minx)*tmp;
					double y1 = yoffest+(f[1]-miny)*tmp;
					double x2 = xoffest+(f[2]-minx)*tmp;
					double y2 = yoffest+(f[3]-miny)*tmp;
					double x3 = xoffest+(f[4]-minx)*tmp;
					double y3 = yoffest+(f[5]-miny)*tmp;
					fprintf(pFile,"%g %g %g %g %g %g T\n",
						x1,y1,x2,y2,x3,y3);
					break;
				}
				case 4:
				{
					double x1 = xoffest+(f[0]-minx)*tmp;
					double y1 = yoffest+(f[1]-miny)*tmp;
					double x2 = xoffest+(f[2]-minx)*tmp;
					double y2 = yoffest+(f[3]-miny)*tmp;
					double x3 = xoffest+(f[4]-minx)*tmp;
					double y3 = yoffest+(f[5]-miny)*tmp;
					double x4 = xoffest+(f[6]-minx)*tmp;
					double y4 = yoffest+(f[7]-miny)*tmp;
					fprintf(pFile,"%g %g %g %g %g %g %g %g P4\n",
						x1,y1,x2,y2,x3,y3,x4,y4);
					break;
				}
				case 5:
				{
					double x1 = xoffest+(f[0]-minx)*tmp;
					double y1 = yoffest+(f[1]-miny)*tmp;
					double x2 = xoffest+(f[2]-minx)*tmp;
					double y2 = yoffest+(f[3]-miny)*tmp;
					double x3 = xoffest+(f[4]-minx)*tmp;
					double y3 = yoffest+(f[5]-miny)*tmp;
					double x4 = xoffest+(f[6]-minx)*tmp;
					double y4 = yoffest+(f[7]-miny)*tmp;
					double x5 = xoffest+(f[8]-minx)*tmp;
					double y5 = yoffest+(f[9]-miny)*tmp;
					fprintf(pFile,"%g %g %g %g %g %g %g %g %g %g P5\n",
						x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
					break;
				}
				case 6:
				{
					double x1 = xoffest+(f[0]-minx)*tmp;
					double y1 = yoffest+(f[1]-miny)*tmp;
					double x2 = xoffest+(f[2]-minx)*tmp;
					double y2 = yoffest+(f[3]-miny)*tmp;
					double x3 = xoffest+(f[4]-minx)*tmp;
					double y3 = yoffest+(f[5]-miny)*tmp;
					double x4 = xoffest+(f[6]-minx)*tmp;
					double y4 = yoffest+(f[7]-miny)*tmp;
					double x5 = xoffest+(f[8]-minx)*tmp;
					double y5 = yoffest+(f[9]-miny)*tmp;
					double x6 = xoffest+(f[10]-minx)*tmp;
					double y6 = yoffest+(f[11]-miny)*tmp;
					fprintf(pFile,"%g %g %g %g %g %g %g %g %g %g %g %g P6\n",
						x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6);
					break;
				}
				default: // > degree 6 (hope rare)
				{
					// Draw a filled polygon
					fprintf(pFile,"newpath\n");
					fprintf(pFile,"%f %f %f C\n",1.0f,1.0f,1.0f);
					double x = xoffest+(f[0]-minx)*tmp;
					double y = yoffest+(f[1]-miny)*tmp;
					fprintf(pFile,"%g %g moveto\n",x,y);
					int j;
					for(j=1;j<degree;j++)
					{
  					double x = xoffest+(f[2*j]-minx)*tmp;
	  				double y = yoffest+(f[2*j+1]-miny)*tmp;
						fprintf(pFile,"%g %g lineto\n",x,y);
					}
					fprintf(pFile, "closepath fill\n");

					// Outline polygon

					fprintf(pFile,"newpath\n");
					fprintf(pFile,"%f %f %f C\n",0.0f,0.0f,0.0f);
					x = xoffest+(f[0]-minx)*tmp;
					y = yoffest+(f[1]-miny)*tmp;
					fprintf(pFile,"%g %g moveto\n",x,y);
					for(j=1;j<degree;j++)
					{
  					double x = xoffest+(f[2*j]-minx)*tmp;
	  				double y = yoffest+(f[2*j+1]-miny)*tmp;
						fprintf(pFile,"%g %g lineto\n",x,y);
					}
					fprintf(pFile, "closepath stroke\n");
				}
			}
		}
			  
		// Emit EPS trailer
		fputs("grestore\n\n",pFile);
		fputs("showpage\n",pFile);
		fclose(pFile);
	}
};

#endif // DUMP_EPS


