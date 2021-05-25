To discuss:
* When computing weights one by one, we recompute certain areas/distances multiple times and so using them e.g. in barycentric coordinates is not efficient. What should we do about that?
* Why do we have in the orbifold_Tutte_parameterizer "+tangent_weight" while in the MVC_post_processor it is "-tangent_weight"? I have fixed that.
* Do I need to cut-zero the dh weights in the orbifold_Tutte_parameterizer? (line 837)? I do cut it for the moment. Otherwise, the results are not equal to the original version.
* Skeletonization uses the weird secure version for the cotangent weights.
* In skeletonization, the final example results are not determenistic.
* Should I remove the positive area from the Tangent_weight and substitute it by computing tan(alpha/2)? In this case, I will keep the correct sign in any configuration.
* What about using a solution that we currently use in the triangulate_hole_with_cdt() instead of flattening?

Later:
* Cleanup tests.
* Comment the code.
* Add more pictures/figures.
* Add a concept test as in the heat_method.
* Try to combine cotangent wrappers from the tools.h. Not sure if this is necessary.
* Try to combine mvc and dhc in the orbifold parameterization. Not sure if this is necessary.
* Mention that tangent_weight_3 uses positive areas (no distortions) and can be used only for PMP, while mean_value_weight_2/3 e.g. can have different signs/distortions for 2D and 3D versions due to the flattening of the 3D region.
* Fix all other packages with respect to the new interface.
* Check memory leaking.

To do now:
* What happens with WP/MV/DH weights on the polygon boundary?
* Check if this code works with the Projection_traits_xy/xz/yz classes.
* Add missing functions in Triangulation_2_projection_3 traits.
