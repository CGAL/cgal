// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Bruno Levy, Pierre Alliez

#include <CGAL/basic.h>

#include "CGAL/parameterization.h"
#include "CGAL/Circular_border_parametizer_3.h"
#include "CGAL/Square_border_parametizer_3.h"
#include "CGAL/Barycentric_mapping_parametizer_3.h"
#include "CGAL/Discrete_conformal_map_parametizer_3.h"
#include "CGAL/Discrete_authalic_parametizer_3.h"
#include "CGAL/Mean_value_coordinates_parametizer_3.h"

#include "options.h"
#include "cgal_types.h"
#include "Mesh_adaptor_polyhedron_ex.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <fstream>


// USAGE EXAMPLES

//----------------------------------------------------------
// conformal parameterization
// circle boundary
// output is a ps map
// input file is mesh.off (the latter must be at the very end)
//----------------------------------------------------------
// ./polyhedron_ex_parameterization -t conformal -b circle -o map mesh.off > mesh.eps

//----------------------------------------------------------
// floater parameterization
// square boundary
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// ./polyhedron_ex_parameterization -t floater -b square -o map mesh.off > mesh.eps

//----------------------------------------------------------
// natural parameterization
// no explicitly pinned vertices
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// ./polyhedron_ex_parameterization -t natural -o map mesh.off > mesh.eps

//----------------------------------------------------------
// natural parameterization
// explicitly pinned vertices (index 5 and 30)
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// ./polyhedron_ex_parameterization -t natural -p -i 5 -j 30 -o map mesh.off > mesh.eps
// equivalent to
// ./polyhedron_ex_parameterization --type natural --pin -i 5 -j 30 --output map mesh.off > mesh.eps

//----------------------------------------------------------
// natural parameterization
// no explicitly pinned vertices
// output is a .obj
// input file is mesh.off
//----------------------------------------------------------
// ./polyhedron_ex_parameterization -t natural -o obj mesh.off > mesh.obj

extern "C" {
   void exit(int);
   long atol(const char *);
}

static const char *  optv[] = {	

   //     '|' -- indicates that the option takes NO argument;
   //     '?' -- indicates that the option takes an OPTIONAL argument;
   //     ':' -- indicates that the option takes a REQUIRED argument;
   //     '*' -- indicates that the option takes 0 or more arguments;
   //     '+' -- indicates that the option takes 1 or more arguments;

   "t:type <string>", // -t or --type
   // -t    conformal   -> default
   //       natural     -> free boundaries
   //       floater     -> mean coordinate values
   //       uniform     -> weight = 1
   //       inverse     -> weight = 1/len
   //       centripetal -> weight = 1/sqrt(len)
   //       david       -> david's stuff
   //       authalic    -> weak area-preserving

   "o:output <string>", // -o or --output
   // -o map   -> eps map (-o map file.off > file.eps)
   //    obj -> .obj) (-o obj file.off > file.obj)

   "p|pin", // -p or --pin (for natural parameterization only)
   // -p

   "i:indexi <number>", // -i or --indexi (for natural parameterization only)
   // -i [0;#boundary vertices-1]

   "j:indexj <number>", // -j or --indexj (for natural parameterization only)
   // -j [0;#boundary vertices-1]

   "b:boundary <string>", // -b or --boundary
   // -b circle
   //    square
   NULL
} ;

// Border parameterization methods
enum boundary_type {BOUNDARY_CIRCLE, BOUNDARY_SQUARE};

int main(int argc,char * argv[])
{
   std::cerr << "\nPARAMETERIZATION" << std::endl;

   // options
   const char *type = "conformal";  // default: conformal param
   const char *boundary = "circle"; // default: circle boundary
   boundary_type boundary_type = BOUNDARY_CIRCLE;
   const char *output = "";         // default:  output nothing
   bool pin = false;          // let the system choose two pinned vertices
   int index_pin[2] = {0,0};  // indices of the two pinned vertices

   // misc
   char  optchar;
   const char *optarg;
   int  errors = 0;
   int  npos = 0;

   Options opts(*argv, optv);
   OptArgvIter iter(--argc,++argv);
   while ((optchar = opts(iter,optarg)) != 0)
   {
      switch(optchar)
      {
        // type of param.
        case 't' :
          type = optarg;
          std::cerr << "  " << type << " parameterization" << std::endl;
          break;

        // output: "map" or "param" or "obj"
        case 'o' :
          output = optarg;
          std::cerr << "  " << "output: " << output << std::endl;
          break;

        // pinned vertices
        case 'p' :
          pin = true; // you decided to choose the pinned vertices
          std::cerr << "  vertices are pinned" << std::endl;
          break;

        // first pinned vertex
        case 'i' :
          index_pin[0] = (int)::atol(optarg);
          std::cerr << "  first pinned vertex: " << index_pin[0] << std::endl;
          break;

        // first pinned vertex
        case 'j' :
          index_pin[1] = (int)::atol(optarg);
          std::cerr << "  second pinned vertex: " << index_pin[1] << std::endl;
          break;

        // boundary
        case 'b' :
          boundary = optarg;
          std::cerr << "  " << boundary << " boundary" << std::endl;
          if(strcmp(boundary,"circle") == 0)
            boundary_type = BOUNDARY_CIRCLE;
          else
            boundary_type = BOUNDARY_SQUARE;
          break;


        // help
        case '?' :
        case 'H' :
           opts.usage(cerr, "files ...");
           ::exit(0);
           break;

        default :
           ++errors;
           break;
      } // end switch
   }

   // check
   int index = iter.index();
   if(errors || ((index == argc) && (npos == 0)))
   {
      if(! errors)
         cerr << opts.name() << ": no filenames given." << endl ;
      opts.usage(cerr, "files ...");
      ::exit(1);
   }

	//***************************************     		
	// for each file...
	//***************************************     		
	if((npos > 0) || (index < argc))
   {
      int first = (npos > 0) ? 0 : index;
      int limit = (npos > 0) ? npos : argc;
      for(int i=first;i<limit;++i)
      {
     		//***************************************     		
     		// check file
     		//***************************************     		
     		fprintf(stderr,"\n  read file...%s...",argv[i]);
     		std::ifstream stream(argv[i]);
     		if(!stream) {
     			fprintf(stderr,"failed: cannot open file!\n");
     			return 1;
     		}
      		
     		//***************************************     		
     		// reads the mesh
     		//***************************************     		
      		Polyhedron_ex mesh;
     		fprintf(stderr,"ok\n  fill mesh...");
     		stream >> mesh;
      		
     		// print mesh info
     		fprintf(stderr,"(%d faces, ",mesh.size_of_facets());
     		fprintf(stderr,"%d vertices)\n",mesh.size_of_vertices());
      		
     		//// compute misc.
     		//mesh.compute_normals();
     		//mesh.compute_facet_centers();
     		//mesh.compute_mean_curvature_normal();

     		//***************************************     		
     		// switch parameterization
     		//***************************************     	
			/*
     		CParameterization parameterization(&mesh);
     		if(strcmp(type,"conformal") == 0)
	     		parameterization.conformal(boundary_type);
			else
    		if(strcmp(type,"authalic")  == 0)
       			parameterization.authalic(boundary_type);
			else
    		if(strcmp(type,"natural") == 0)
       			parameterization.natural(pin,index_pin);
			else
    		if(strcmp(type,"floater") == 0)
       			parameterization.floater(boundary_type);
			else
    		if(strcmp(type,"uniform") == 0)
       			parameterization.hooke(0,boundary_type);
			else
    		if(strcmp(type,"inverse") == 0)
       			parameterization.hooke(-1,boundary_type);
			else
    		if(strcmp(type,"centripetal") == 0)
     			parameterization.hooke(-0.5,boundary_type);
			else
			{
     			fprintf(stderr,"invalid choice\n");
     			return 1;
     		}
			*/
			CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::ErrorCode err;
			if ( (strcmp(type,"uniform") == 0) && (strcmp(boundary,"circle") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Barycentric_mapping_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else if ( (strcmp(type,"uniform") == 0) && (strcmp(boundary,"square") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Barycentric_mapping_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else if ( (strcmp(type,"conformal") == 0) && (strcmp(boundary,"circle") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Discrete_conformal_map_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else if ( (strcmp(type,"conformal") == 0) && (strcmp(boundary,"square") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Discrete_conformal_map_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else if ( (strcmp(type,"authalic") == 0) && (strcmp(boundary,"circle") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Discrete_authalic_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else if ( (strcmp(type,"authalic") == 0) && (strcmp(boundary,"square") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Discrete_authalic_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else if ( (strcmp(type,"floater") == 0) && (strcmp(boundary,"circle") == 0) )
			{
				// Floater/circle is the default parameterization algorithm
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh));
			}
			else if ( (strcmp(type,"floater") == 0) && (strcmp(boundary,"square") == 0) )
			{
				err = CGAL::parameterize(&Mesh_adaptor_polyhedron_ex(&mesh),
										 CGAL::Mean_value_coordinates_parametizer_3<Mesh_adaptor_polyhedron_ex, 
										 CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
			}
			else 
			{
     			fprintf(stderr,"invalid choice\n");
     			return -1;
     		}
			if (err != CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::OK)
			{
				fprintf(stderr,"parameterization failure: error = %d\n", (int)err);
     			return err;
     		}
     		
     		//***************************************     		
     		// output
     		//***************************************     		
     		if(strcmp(output,"map") == 0)
       			mesh.dump_param(); // .eps stream
     		else
     		if(strcmp(output,"obj") == 0)
       			mesh.write_file_obj();			 // write .obj file
		}

		// Flush cout and cerr
		cout << endl ;
		cerr << endl ;
	}

	return EXIT_SUCCESS;
}

