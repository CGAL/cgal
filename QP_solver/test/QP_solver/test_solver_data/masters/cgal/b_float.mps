* Number-type: floating-point
* Description: optimally solvable lp for which glpk claims it's unbounded
* Generated-by: http://members.jcom.home.ne.jp/masashi777/exlp.html 
NAME          from_lp_
ROWS
 N  r_0
 G  r_1
 L  r_2
COLUMNS
    x         r_0                 -1   r_1                  1
    x         r_2             0.0001
    y         r_1                  1   r_2              50001
    z         r_1                 -1
RHS
    RHS       r_2                  1
ENDATA
