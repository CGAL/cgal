* There is a random behavior for the test_polygons_ac.off case. Sometimes it works and sometimes not.
* Polygon_splitter bugs for the case test_1_polygon_b.off and for the case test_2_polygons_ad.off.
* When we create a new face, should we use the initial k for it, e.g. 2 or should we copy the last counted k, e.g. 2 or 1?
* Should we decrease k in the polygon splitter when two polygons are intersected at the very beginning?
* Should we count the number of intersections per polygon or per mesh face?
* Should we keep facei-facej number of intersections?
* When we do a test e.g. 6 polygons intersect k = 6 times all intersections between polygons are already inserted for k = 5 and do not change for k = 6. For 3 input polygons, they are inserted only for k = 3. Is it correct behaviour?