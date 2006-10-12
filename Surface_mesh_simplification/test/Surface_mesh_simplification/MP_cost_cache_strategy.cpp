typedef SMS::Edge_length_cost<Surface> ActualCost ;

ActualCost actual_cost ;

SMS::Set_cost_cache<Surface,ActualCost> cache(actual_cost) ; 

SMS::Cached_cost       <Surface> cost ;
SMS::Midpoint_placement<Surface> placement ;
                  
