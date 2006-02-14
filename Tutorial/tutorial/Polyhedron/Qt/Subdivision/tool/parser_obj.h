#include <CGAL/Polyhedron_incremental_builder_3.h>
#include "lib/Enriched_polyhedron.h"

template <class HDS>
class Builder_obj : public CGAL::Modifier_base<HDS>
{
private:
	typedef typename HDS::Vertex::Point Point;
  typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
  FILE *m_pFile;

public:
  Builder_obj(FILE *pFile)
  {
    m_pFile = pFile;
  }
	~Builder_obj() {}

  void operator()(HDS& hds)
  {
    Builder builder(hds,true);
    builder.begin_surface(3,1,6);
			read_vertices(builder);
      read_facets(builder);
    builder.end_surface();
  }

private:
  // read vertex coordinates
  void read_vertices(Builder &builder)
  {
    fseek(m_pFile,0,SEEK_SET);
  	char pLine[512];
  	while(fgets(pLine,512,m_pFile))
  		if(pLine[0] == 'v')
  		{
      	float x,y,z;
  			if(sscanf(pLine,"v %f %f %f",&x,&y,&z) == 3)
					builder.add_vertex(Point(x,y,z));
  		}
  }

  // read facets and uv coordinates per halfedge
  void read_facets(Builder &builder)
  {
    fseek(m_pFile,0,SEEK_SET);
  	char pLine[512];
  	while(fgets(pLine,512,m_pFile))
  	{
  		char *pTmp = pLine;
  		if(pTmp[0] == 'f')
  		{
				int index,n;
				char index_ascii[512],n_ascii[512];

        // create facet        
        builder.begin_facet();

				pTmp += 2; // jump after 'f '
				if(strstr(pTmp,"//"))
					while(sscanf(pTmp,"%d//%d",&index,&n))
					{
						itoa(index,index_ascii,10);
						itoa(n,n_ascii,10);
		        builder.add_vertex_to_facet(index-1);
						pTmp += (2 + strlen(index_ascii) + strlen(n_ascii));
						if(strlen(pTmp) < 3)
							break;
						else
							pTmp += 1;
					}
				else
					while(sscanf(pTmp,"%d",&index))
					{
						itoa(index,index_ascii,10);
						pTmp += strlen(index_ascii);
		        builder.add_vertex_to_facet(index-1);
						if(strlen(pTmp) < 3)
							break;
						else
							pTmp += 1;
					}
        builder.end_facet();
			}
		}
  }
};

template <class	kernel,	class	items>
class Parser_obj
{
public:
    typedef typename Enriched_polyhedron<kernel,items>::HalfedgeDS HalfedgeDS;
    Parser_obj() {}
    ~Parser_obj() {}

public:
    bool read(const char*pFilename,
			        Enriched_polyhedron<kernel,items> *pMesh)
    {
      CGAL_assertion(pMesh != NULL);
      FILE *pFile = fopen(pFilename,"rt");
      if(pFile == NULL)
        return false;
      Builder_obj<HalfedgeDS> builder(pFile);
      pMesh->delegate(builder);
      fclose(pFile);
      return true;
    }
};

