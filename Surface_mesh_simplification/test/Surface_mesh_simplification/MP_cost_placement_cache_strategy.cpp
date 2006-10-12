typedef SMS::Edge_length_cost  <Surface> ActualCost ;
typedef SMS::Midpoint_placement<Surface> ActualPlacement ;

ActualCost      actual_cost ;
ActualPlacement actual_placement ;

SMS::Set_cost_and_placement_cache<Surface,ActualCost,ActualPlacement> cache(actual_cost,actual_placement) ; 

SMS::Cached_cost     <Surface> cost ;
SMS::Cached_placement<Surface> placement ;
                  
