        //std::list<Cell_handle> incidents_c;

        //incident_cells(static_cast<Vertex_handle>(fv),
         // back_inserter(incidents_c));

       // typename std::list<Cell_handle>::iterator cit = incidents_c.begin();

        /*if (is_infinite(*chit)) as->set_is_on_chull(true);*/

        //while (is_infinite(*cit)) ++cit; //skip infinte cells

        //for( ; cit != incidents_c.end(); ++cit) 
        //{

        //}


        /*alpha = (*chit)->get_alpha();
        as->set_alpha_mid(alpha);
        as->set_alpha_max(alpha);
        for( ; chit != incidents.end(); ++chit) {
          if (is_infinite(*chit)) as->set_is_on_chull(true);
          else {
            alpha = (*chit)->get_alpha();
            if (alpha < as->alpha_mid()) as->set_alpha_mid(alpha);
            if (alpha > as->alpha_max()) as->set_alpha_max(alpha);
          }*/
        
//        Finite_cell_iterator fc;
//        this->incident_cells(fv->vertex(),fc);
        
//        for(; fc != this->finite_facets_end(); ++ff)
//        {
        
//        Facet_circulator fcir_begin = ff->facet_begin();
//        Facet_circulator fcir = fcir_begin;





  Cell_handle ray_intersection(Vertex_handle vh, Ray_3 ray) 
  {
    // create a list of all the incident cells to the concerned vertex handle
    std::list<Cell_handle> incident_c;
    incident_cells(static_cast<Vertex_handle>(vh),back_inserter(incident_c));

    // create an iterator for these cells
    typename std::list<Cell_handle>::iterator cell_it = incident_c.begin();

    Cell_handle result = NULL;

    //while (is_infinite(*cell_it)) ++cell_it; //skip infinite cells

    for( ; cell_it != incident_c.end(); ++cell_it) 
    {
      for (int k = 0; k<4; ++k) 
      {
        Triangle_3 cell_facet = construct_triangle((*cell_it)->vertex( (k+1)&3 )->point(),
          (*cell_it)->vertex( (k+2)&3 )->point(),
          (*cell_it)->vertex( (k+3)&3 )->point());     
        if(do_intersect(cell_facet,ray))
        {
          (*cell_it)->info() = OUTSIDE;
          result = (*cell_it);
        }
      }
    }
    return result;
  }



void inside_outside() 
  {
    Finite_vertices_iterator fv = this->finite_vertices_begin();

    int counter = 1;

    for(;fv != this->finite_vertices_end(); ++fv)
    {
      counter++;

      std::vector<Point_3> list_of_cams = fv->info().get_list_of_cameras();

      // construct rays emanating from the vertex to cameras.
      for(int j = 0; j<fv->info().get_list_of_cameras().size(); ++j)
      {
        Ray_3  ray (fv->point(), list_of_cams[j]);

        // the tetrahedron behind the vertex is labelled as inside
        std::list<Cell_handle> incident_c;
        incident_cells(static_cast<Vertex_handle>(fv),back_inserter(incident_c));

        // create an iterator for these cells
        typename std::list<Cell_handle>::iterator cell_it = incident_c.begin();

        while (is_infinite(*cell_it)) ++cell_it; //skip infinite cells

        for( ; cell_it != incident_c.end(); ++cell_it) 
        {
          for (int k = 0; k<4; ++k) 
          {
            Triangle_3 cell_facet = construct_triangle((*cell_it)->vertex( (k+1)&3 )->point(),
              (*cell_it)->vertex( (k+2)&3 )->point(),
              (*cell_it)->vertex( (k+3)&3 )->point());     
            if(do_intersect(cell_facet,ray.opposite()))
            {
              (*cell_it)->info() = INSIDE;
            }
          }
        }
        
        

        Vertex_handle v1,v2,v3;
        v1 = fv;
        v2 = v1;      
        v3 = v1;

        Cell_handle c1,c2,c3;
        c1 = ray_intersection(v1,ray);
        c2 = c1;
        c3 = c1;

        while(!this->is_infinite(c1) || !this->is_infinite(c2) || !this->is_infinite(c3))
        {

          if(this->is_cell(c1))
          {
            int k = c1->index(v1);
            v1 = c1->vertex((k+1)&3);
            v2 = c1->vertex((k+2)&3);
            v3 = c1->vertex((k+3)&3);
            c1 = ray_intersection(v1,ray);
            c2 = ray_intersection(v2,ray);
            c3 = ray_intersection(v3,ray);
          }

          else if(this->is_cell(c2))
          {
            int k = c2->index(v2);
            v1 = c2->vertex((k+1)&3);
            v2 = c2->vertex((k+2)&3);
            v3 = c2->vertex((k+3)&3);
            c1 = ray_intersection(v1,ray);
            c2 = ray_intersection(v2,ray);
            c3 = ray_intersection(v3,ray);
          }

          else if(this->is_cell(c3))
          {
            int k = c3->index(v3);
            v1 = c3->vertex((k+1)&3);
            v2 = c3->vertex((k+2)&3);
            v3 = c3->vertex((k+3)&3);
            c1 = ray_intersection(v1,ray);
            c2 = ray_intersection(v2,ray);
            c3 = ray_intersection(v3,ray);
          }
        }
      } 
    }
  }





  void inside_outside() 
  {
    int counter_vertex = 0;
    int counter_boucle = 0;
    int counter_ray = 0;
    int counter_behind = 0;
    int counter_rayintersect = 0;
    int counter_k = 0;

    Finite_vertices_iterator fv = this->finite_vertices_begin();
    int counter_opposite = 0;

    for(;fv != this->finite_vertices_end(); ++fv)
    {
      counter_vertex++;

      std::vector<Point_3> list_of_cams = fv->info().get_list_of_cameras();

      // construct rays emanating from the vertex to cameras.
      for(int j = 0; j<fv->info().get_list_of_cameras().size(); ++j)
      {
        counter_ray++;

        Ray_3  ray (fv->point(), list_of_cams[j]);

        // Hack to make it work : take the point on ray.
        // point(0) is the source. point(1) is the one after the source.
        double epsilon = 0.01;
        Vector_3 ray_direction = ray.to_vector();

        double norm = sqrt(ray_direction.x()*ray_direction.x() + 
          ray_direction.y()*ray_direction.y() + ray_direction.z()*ray_direction.z());

        // a changer
        Point_3 next_on_ray(fv->point().x()+epsilon*ray_direction.x()/norm,
          fv->point().y()+epsilon*ray_direction.y()/norm,
          fv->point().z()+epsilon*ray_direction.z()/norm);

        Ray_3  adapted_ray (next_on_ray, list_of_cams[j]);

        std::list<Cell_handle> incident_c;
        incident_cells(static_cast<Vertex_handle>(fv),back_inserter(incident_c));

        // create an iterator for these cells
        typename std::list<Cell_handle>::iterator cell_it = incident_c.begin();

        // while (is_infinite(*cell_it)) ++cell_it; //skip infinite cells
        Cell_handle next_cell_h;
        int facet_number;

        for( ; cell_it != incident_c.end(); ++cell_it) 
        {
          // Hack: to be able to label the tetrahedron behind the vertex as inside
          counter_opposite = 0;

          for (int k = 0; k<4; ++k) 
          {

            // skip infinite facet
            if (is_infinite(*cell_it, k))
              continue;

            Triangle_3 cell_facet = construct_triangle((*cell_it)->vertex( (k+1)&3 )->point(),
              (*cell_it)->vertex( (k+2)&3 )->point(),
              (*cell_it)->vertex( (k+3)&3 )->point());     

            if(do_intersect(cell_facet,ray.opposite()))
              counter_opposite++;

            if(counter_opposite == 4)
            {
              counter_behind++;
              TRACE("INSIDE\n");
              (*cell_it)->info() = INSIDE;
            }

            if(do_intersect(cell_facet,adapted_ray))
            {
              TRACE("OUTSIDE\n");
              (*cell_it)->info() = OUTSIDE;
              next_cell_h = (*cell_it)->neighbor(k); // A VOIRRRRRR EDGE BON COTE
              facet_number = this->mirror_index((*cell_it),k);
            }           
          }
        }

//        while(next_cell_h != NULL && !is_infinite(next_cell_h->vertex(facet_number)))
        while(next_cell_h != NULL && !is_infinite(next_cell_h))
        {
          counter_boucle++;
          TRACE("while: next_cell_h=%x, facet_number=%d\n", next_cell_h, facet_number);

          //Facet_handle old_facet_h = next_facet_h;
          Cell_handle old_cell_h = next_cell_h;//old_facet_h->first;

          next_cell_h = NULL;

          for (int k = 0; k<4; ++k) 
          {
            // skip previous facet
            if (k != facet_number)
            {   

              Triangle_3 cell_facet = construct_triangle(old_cell_h->vertex( (k+1)&3 )->point(),
                old_cell_h->vertex( (k+2)&3 )->point(),
                old_cell_h->vertex( (k+3)&3 )->point());     

              if(do_intersect(cell_facet,ray))
              {
                counter_rayintersect++;
                TRACE("OUTSIDE\n");
                old_cell_h->info() = OUTSIDE;
                next_cell_h = old_cell_h->neighbor(k); // A VOIRRRRRR EDGE BON COTE
                facet_number = this->mirror_index(old_cell_h,k);
              }           
            }
          }
        }
      }
    }
  }



//  Finite_vertices_iterator fv = this->finite_vertices_begin();

// for(;fv != this->finite_vertices_end(); ++fv)
// {
// vertex per vertex
//Cell_handle c = fv->cell();

// corresponding cell is labelled as inside.
//c->info()     = INSIDE;

//      std::vector<Point_3> list_of_cams = fv->info().get_list_of_cameras();

// construct rays emanating from the vertex to cameras.
// for(int j = 0; j<fv->info().get_list_of_cameras().size(); ++j)
// {
// Ray_3  ray (fv->point(), list_of_cams[j]);
//Finite_facets_iterator ff = this->finite_facets_begin();;
//
//        // All the tetrahedra intersected by a ray emanating from the vertex to
//        // the camera center of one of these views should be labbelled as outside.
//        for(; ff != this->finite_facets_end(); ++ff)
//        {
//          int k = ff->second;
//          Triangle_3 facet_tri = 
//            construct_triangle(ff->first->vertex( (k+1)&3 )->point(),
//            ff->first->vertex( (k+2)&3 )->point(),
//            ff->first->vertex( (k+3)&3 )->point());
//
//          if(do_intersect(facet_tri, ray))
//          {
//            // both tetrahedrons sharing this face are labelled as outside.
//            //Cell_handle c1 = ff->first;
//            Cell_handle c2 = ff->first->neighbor(ff->second);
//
//            //c1->info() = OUTSIDE;
//            c2->info() = OUTSIDE;


//}
//}

        std::vector<Point_3> list_of_cameras;
        for (int i = 0 ; i< sizeof(rest_of_line); ++i) rest_of_line[i] = '\0'; // initialize to EOL 

        if(sscanf(pLine,"%lg\t%lg\t%lg\t%512c",&x,&y,&z,&rest_of_line)== 4)
        {
        
        TRACE("Line number : %d, List_cameras : %s\n", lineNumber,rest_of_line);
        
          int value;
          std::istringstream iss(rest_of_line);
          while (iss >> value)
          
          
          
          
          
          
          
          while(next_cell_h != NULL && !is_infinite(next_cell_h))
        {

          //TRACE("while: next_cell_h=%x, facet_number=%d\n", next_cell_h, facet_number);

          Cell_handle old_cell_h = next_cell_h;
          int old_facet_number = current_facet_number;
          next_cell_h = NULL;
          
          std::vector<int> intersected_facets;

          for (int k = 0; k<4; ++k) 
          {
            // skip previous facet
            if (k != facet_number)
            {   
              Triangle_3 cell_facet = this->triangle(old_cell_h,k);     

              if(do_intersect(cell_facet,segment))
              {
                intersected_facets.push_back(k);
              }
            }
          }

          assert(intersected_facets.size() >=0);
          //TRACE("intersected_facets.size()=%d\n", intersected_facets.size());

          if(intersected_facets.size()== 1) //intersects single facet
          {
            //TRACE("Case SINGLE FACET\n");
            old_cell_h->info() = OUTSIDE;
            next_cell_h = old_cell_h->neighbor(intersected_facets[0]);
            facet_number = this->mirror_index(old_cell_h,intersected_facets[0]);
          }

          else if(intersected_facets.size()== 2) //intersects an edge
          {
            TRACE("Case EDGE\n");
            old_cell_h->info() = OUTSIDE;
            Edge edge(old_cell_h,intersected_facets[0],intersected_facets[1]);
            Cell_circulator ccir = this->incident_cells(edge,old_cell_h);
            Cell_circulator cdone = ccir;

            // these are the indexes of the vertices in the ccir cell
            int index_vertex_0, index_vertex_1; 

            do {

              if(is_infinite(ccir)) ccir++;

              if((*ccir).has_vertex(old_cell_h->vertex(intersected_facets[0]),index_vertex_0)){}

              if((*ccir).has_vertex(old_cell_h->vertex(intersected_facets[1]),index_vertex_1)){}

              Triangle_3 cell_facet_0 = this->triangle(ccir,index_vertex_0);
              Triangle_3 cell_facet_1 = this->triangle(ccir,index_vertex_1);

              if(do_intersect(cell_facet_0,segment) && !do_intersect(cell_facet_1,segment))
              {
                next_cell_h = old_cell_h->neighbor(intersected_facets[0]);
                facet_number = this->mirror_index(old_cell_h,intersected_facets[0]);
              }

              else if(do_intersect(cell_facet_1,segment)&& !do_intersect(cell_facet_0,segment))
              {
                next_cell_h = old_cell_h->neighbor(intersected_facets[1]);
                facet_number = this->mirror_index(old_cell_h,intersected_facets[1]);
              }

              else
              {// go to next cell (ATTENTION : CASE intersects both faces = Intersects common edge)
                ++ccir;
              }
            } while ( ccir != cdone );
          }

          else //intersects a vertex
          {
            TRACE("Case VERTEX\n");
            std::list<Cell_handle> in_c;
            incident_cells(old_cell_h->vertex(facet_number),back_inserter(in_c));

            // create an iterator for these cells
            typename std::list<Cell_handle>::iterator cell_it = in_c.begin();

            for(; cell_it != in_c.end(); ++cell_it) 
            {

              TRACE("for: old_cell_h=%x, cell_it=%x, next_cell_it=%x\n", old_cell_h, *cell_it, next_cell_h);

              // skip infinite cell
              if (is_infinite(*cell_it) || next_cell_h == old_cell_h)
                continue;

              // this is the index of the vertex in the cell_it cell
              int index_in_cell_it;

              if((*cell_it)->has_vertex(old_cell_h->vertex(facet_number),index_in_cell_it))
              {
                Triangle_3 cell_facet = this->triangle((*cell_it),index_in_cell_it); 

                if(do_intersect(cell_facet,segment))
                {
                  TRACE("OUTSIDE\n");
                  (*cell_it)->info() = OUTSIDE;
                  next_cell_h = (*cell_it)->neighbor(index_in_cell_it);
                  facet_number = this->mirror_index((*cell_it),index_in_cell_it);
                }
              }
            }
          }
        }
        
        
        
         void inside_outside() 
  {
    int counter_vertex = 0;
    int counter_boucle = 0;
    int counter_ray = 0;
    int counter_behind = 0;
    int intersected_facets = 0;
    int counter_k = 0;

    Finite_vertices_iterator fv = this->finite_vertices_begin();
    int counter_opposite = 0;

    for(;fv != this->finite_vertices_end(); ++fv)
    {
      counter_vertex++;

      TRACE("Vertex : %x\n",fv);

      //if(counter_vertex == 33){TRACE("Vertex : %x, List_of_cameras_size : %d\n",fv,fv->info().get_list_of_cameras().size());}

      std::vector<Point_3> list_of_cams = fv->info().get_list_of_cameras();

      TRACE("Vertex : %d, List_of_cameras_size : %d\n",counter_vertex,fv->info().get_list_of_cameras().size());

      // construct rays emanating from the vertex to cameras.
      for(int j = 0; j<fv->info().get_list_of_cameras().size(); ++j)
      {

        Segment_3  segment (fv->point(), list_of_cams[j]);

        std::list<Cell_handle> incident_c;
        incident_cells(static_cast<Vertex_handle>(fv),back_inserter(incident_c));

        // create an iterator for these cells
        typename std::list<Cell_handle>::iterator cell_it = incident_c.begin();

        Cell_handle next_cell_h;
        int facet_number;

        for( ; cell_it != incident_c.end(); ++cell_it) 
        {

          // skip infinite facet
          if (is_infinite(*cell_it))
            continue;

          // this is the index of the vertex in the cell_it cell
          int index_in_cell_it;
          if((*cell_it)->has_vertex(fv,index_in_cell_it))
          {
            Triangle_3 cell_facet = this->triangle((*cell_it),index_in_cell_it); 

            if(do_intersect(cell_facet,segment))
            {
              //TRACE("OUTSIDE\n");
              (*cell_it)->info() = OUTSIDE;
              next_cell_h = (*cell_it)->neighbor(index_in_cell_it);
              facet_number = this->mirror_index((*cell_it),index_in_cell_it);
            }           
          }          
        }

        while(next_cell_h != NULL && !is_infinite(next_cell_h))
        {

          //TRACE("while: next_cell_h=%x, facet_number=%d\n", next_cell_h, facet_number);

          Cell_handle old_cell_h = next_cell_h;
          next_cell_h = NULL;

          std::vector<int> intersected_facets;

          for (int k = 0; k<4; ++k) 
          {
            // skip previous facet
            if (k != facet_number)
            {   
              Triangle_3 cell_facet = this->triangle(old_cell_h,k);     

              if(do_intersect(cell_facet,segment))
              {
                intersected_facets.push_back(k);
              }
            }
          }

          assert(intersected_facets.size() >=0);
          //TRACE("intersected_facets.size()=%d\n", intersected_facets.size());

          if(intersected_facets.size()== 1) //intersects single facet
          {
            //TRACE("Case SINGLE FACET\n");
            old_cell_h->info() = OUTSIDE;
            next_cell_h = old_cell_h->neighbor(intersected_facets[0]);
            facet_number = this->mirror_index(old_cell_h,intersected_facets[0]);
          }



          else if(intersected_facets.size()== 2) //intersects an edge
          {
            TRACE("Case EDGE\n");

            old_cell_h->info() = OUTSIDE;

            int ind_0,ind_1;

            other_two_indices(intersected_facets[0],intersected_facets[1],&ind_0, &ind_1);

            Edge edge(old_cell_h, ind_0, ind_1);

            Cell_circulator ccir = this->incident_cells(edge);
            Cell_circulator cdone = ccir;

            // these are the indexes of the vertices in the ccir cell
            int index_vertex_0, index_vertex_1; 

            do {

              if(!is_infinite(ccir) && (Cell_handle)ccir != old_cell_h)
              {

                if((*ccir).has_vertex(old_cell_h->vertex(ind_0),index_vertex_0)){}

                if((*ccir).has_vertex(old_cell_h->vertex(ind_1),index_vertex_1)){}

                Triangle_3 cell_facet_0 = this->triangle(ccir,index_vertex_0);
                Triangle_3 cell_facet_1 = this->triangle(ccir,index_vertex_1);

                if(do_intersect(cell_facet_0,segment) && !do_intersect(cell_facet_1,segment))
                {
                  TRACE("Case do_intersect(cell_facet_0,segment)\n");
                  // A verifier que je sois pas en train de sauter une etape
                  (*ccir).info() = OUTSIDE;
                  next_cell_h = (*ccir).neighbor(index_vertex_0);
                  facet_number = index_vertex_0;
                  TRACE("Next_cell_h :%x, Facet_number:%d\n", next_cell_h, facet_number);
                }

                else if(do_intersect(cell_facet_1,segment)&& !do_intersect(cell_facet_0,segment))
                {
                  TRACE("Case do_intersect(cell_facet_1,segment)\n");
                  (*ccir).info() = OUTSIDE;
                  next_cell_h = (*ccir).neighbor(index_vertex_1);
                  facet_number = index_vertex_1;
                }
              }

              // go to next cell (ATTENTION : CASE intersects both faces = Intersects common edge)
              TRACE("ELSE\n");
              ++ccir;

            } while ( ccir != cdone );
          }

          else //intersects a vertex
          {
            TRACE("Case VERTEX\n");

            old_cell_h->info() = OUTSIDE;

            std::list<Cell_handle> in_c;
            incident_cells(old_cell_h->vertex(facet_number),back_inserter(in_c));

            // create an iterator for these cells
            typename std::list<Cell_handle>::iterator cell_it = in_c.begin();

            for(; cell_it != in_c.end(); ++cell_it) 
            {

              TRACE("for: old_cell_h=%x, cell_it=%x, next_cell_it=%x\n", old_cell_h, *cell_it, next_cell_h);

              // skip infinite cell
              if (is_infinite(*cell_it) || *cell_it == old_cell_h)
                continue;

              // this is the index of the vertex in the cell_it cell
              int index_in_cell_it;

              if((*cell_it)->has_vertex(old_cell_h->vertex(facet_number),index_in_cell_it))
              {
                Triangle_3 cell_facet = this->triangle((*cell_it),index_in_cell_it); 

                if(do_intersect(cell_facet,segment))
                {
                  TRACE("OUTSIDE\n");
                  (*cell_it)->info() = OUTSIDE;
                  next_cell_h = (*cell_it)->neighbor(index_in_cell_it);
                  facet_number = this->mirror_index((*cell_it),index_in_cell_it);
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  //void inside_outside() 
  //{
  //  Finite_facets_iterator ff = this->finite_facets_begin();

  //  for(;ff != this->finite_facets_end(); ++ff)
  //  {
  //    // vertex per vertex
  //    for(int i=0; i<3; ++i)
  //    {
  //      int o;

  //      if (ff->second%2==1) 
  //        o = (ff->second+i+1)%4;

  //      else 
  //      {
  //        o= (ff->second+(2-i)+1)%4;
  //      }

  //      Vertex_handle v =  ff->first->vertex(o);

  //      // corresponding cell is labelled as inside.
  //      ff->first->info()     = INSIDE;

  //      std::vector<Point_3> list_of_cams = v->info().get_list_of_cameras();

  //      // construct rays emanating from the vertex to cameras.
  //      for(int j = 0; j<v->info().get_list_of_cameras().size(); ++j)
  //      {
  //        Ray_3  ray (v->point(), list_of_cams[j]);

  //        Finite_facets_iterator ffi = ff;

  //        // All the tetrahedra intersected by a ray emanating from the vertex to
  //        // the camera center of one of these views should be labbelled as outside.
  //        for(++ffi ; ffi != this->finite_facets_end(); ++ffi)
  //        {
  //          int k = ffi->second;
  //          Triangle_3 facet_tri = 
  //            construct_triangle(ffi->first->vertex( (k+1)&3 )->point(),
  //            ffi->first->vertex( (k+2)&3 )->point(),
  //            ffi->first->vertex( (k+3)&3 )->point());


  //          if(do_intersect(facet_tri, ray))
  //          {
  //            // both tetrahedrons sharing this face are labelled as outside.
  //            Cell_handle c1 = ffi->first;
  //            //Cell_handle c2 = ffi->first->neighbor(ffi->second);

  //            c1->info() = OUTSIDE;
  //            //c2->info() = OUTSIDE;


  //Finite_facets_iterator ff = this->finite_facets_begin();;

  //      // All the tetrahedra intersected by a ray emanating from the vertex to
  //      // the camera center of one of these views should be labbelled as outside.
  //      for(; ff != this->finite_facets_end(); ++ff)
  //      {
  //        int k = ff->second;
  //        Triangle_3 facet_tri = 
  //          construct_triangle(ff->first->vertex( (k+1)&3 )->point(),
  //          ff->first->vertex( (k+2)&3 )->point(),
  //          ff->first->vertex( (k+3)&3 )->point());

  //        if(do_intersect(facet_tri, ray))
  //        {
  //          // both tetrahedrons sharing this face are labelled as outside.
  //          //Cell_handle c1 = ff->first;
  //          Cell_handle c2 = ff->first->neighbor(ff->second);

  //          //c1->info() = OUTSIDE;
  //          c2->info() = OUTSIDE;
  //          }
  //        } 
  //      }
  //    }
  //  }
  //}
  
  
  //void inside_outside() {

//  Finite_vertices_iterator fv = this->finite_vertices_begin();

//  for(; fv != this->finite_vertices_end(); ++fv)
//  {
//    // the tetrahedron behind the vertex should be labelled as inside.
//    Cell_handle ch = this->locate(fv->point());
//    ch->info()     = INSIDE;


//    for(int i = 0; i<fv->info().get_vector().size(); ++i)
//    {
//      Ray_3  ray = construct_ray(fv->point(), fv->info().get_vector()[i]);

//      // All the tetrahedra intersected by a ray emanating from the vertex to
//      // the camera center of one of these views should be labbelled as outside
//      for ( ++f_it ; f_it != this->finite_facets_end(); f_it++)
//      {
//        Triangle_3 facet_tri = 
//          construct_triangle(f_it->halfedge()->opposite()->vertex()->point(),
//          f_it->halfedge()->vertex()->point(),
//          f_it->halfedge()->next()->vertex()->point());

//        if (do_intersect(facet_tri, ray)){
//          return false;
//        }
//     
//     
//     
//     
//     
//      Cell_handle ch = ff->first;



//      iter++;
//    }              

//      // tetrahedron linked to vertex is INSIDE

//      // construct rays from the vertex to the different cameras

//      // see if one of these rays intersects any other cells

//      // if true put in the info of the cell OUTSIDE



//  }
//}