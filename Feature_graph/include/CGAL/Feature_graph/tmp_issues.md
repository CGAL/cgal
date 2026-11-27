# Feature_graph review issues

## High

1. **Public feature detection is not implemented** — `detect_sharp_features()` returns a default-constructed empty `Feature_graph` in [detect_sharp_features.h](detect_sharp_features.h), and the surface overload is only a stub. The example therefore appears to work while producing no features. **Solution:** implement and return the thinned/regularized/optimized graph (and surface overload), or clearly remove the public entry points until implemented. **Response:** this implementation is in progress, the missing regularization and optimization will be added later.

2. **`value_outside` is inconsistent between stages** — The labeled-image overload reads `voxel_value` in [detect_sharp_features.h](detect_sharp_features.h), but `Image_AmbrosioTortorelli(image)` is constructed with its default outside value and the DGtal surface creation does not receive this value. **Solution:** pass the same named parameter value through `Image_AmbrosioTortorelli`, `make_DGtal_surface_from_image()`, the binary predicate, and image-to-mesh conversion. **Response:** removed `value_outside` in DGtal_interface since this library does not handle it, resulting in additional value created, which is not a problem.

3. **Named parameter types are not validated** — `sharpness_estimator` and `normal_estimator` are selected but their descriptor, domain, and return types are not concept-checked. Wrong estimator types can fail late with cryptic diagnostics. **Solution:** add CGAL concept checks or constrained overload checks and document the exact accepted types. **Response:** the API documentation is separated. Check if the implementation follows the provided documentation.

4. **The feature priority map is typed too specifically** — `thin_list_graph_with_priority()` stores priorities in `std::map<double, ...>`, although `Sharpness_number_type` is a concept-defined real-embeddable type. **Solution:** make the priority key generic or store `Sharpness_number_type` as a template parameter. **Response:** used `boost::property_traits<PriorityValuePropertyMap>::value_type`.

5. **Khalimsky index conversion is duplicated and unsafe** — [detect_sharp_features.h](detect_sharp_features.h) and [DGtal_interface.h](internal/DGtal_interface.h) hand-roll coordinate indices. Negative coordinates, overflow, and mismatch between coordinate widths can cause collisions or invalid indexes. **Solution:** centralize one checked conversion helper, validate positive spacings and dimensions, and test boundary coordinates. **Response:** created `image_Khalimsky_coordinates.h` file with centralized methods.

6. **The output graph is not output as documented** — Public overloads, regularization, optimization, border merging, and graph construction steps are TODOs, while the documentation promises a complete detector. **Solution:** implement the stages or mark the API experimental/nonimplemented and qualify unsupported parameters. **Response:** this implementation is in progress, the missing regularization and optimization will be added later.

## Medium

7. **Missing self-contained standard includes** — Headers rely transitively on `<map>`, `<set>`, `<vector>`, `<array>`, `<variant>`, `<unordered_map>`, `<fstream>`, `<cmath>`, and `<cstddef>`. **Solution:** include every directly-used standard header in each header. **Response:** done by agent.

8. **Debug PLY normalization divides by zero and casts out of range** — [debugging_methods.h](internal/debugging_methods.h) divides by `max_sharpness` or `min_sharpness`, and casts floats to `int` without clamping. **Solution:** handle zero ranges explicitly, normalize to a shared range, and clamp each channel to `[0,255]`. **Response:** clamped values.

9. **Debug color handling is undefined for flat or zero-only ranges** — Positive, negative, and zero-only values are not handled consistently. **Solution:** compute separate positive/negative maxima and emit exact zero as a neutral color. **Response:** sharpness span positive and negative values, plus no need to supercharge this debugging method.

10. **The normal estimator is unused in feature extraction** — `NormalEstimator` is accepted and passed to debug output, but the core thinning function does not use normals. **Solution:** remove it from the initial/thinning API until regularization/optimization uses it, or document it as reserved. **Response:** Not implemented yet.

11. **Barycenter construction is duplicated** — [detect_sharp_features.h](detect_sharp_features.h), [feature_graph_operation.h](internal/feature_graph_operation.h), and [debugging_methods.h](internal/debugging_methods.h) calculate the same face or mesh centroid. **Solution:** add one internal centroid helper and reuse it everywhere. **Response:** Debugging duplication is fine in this context.

12. **Debug output does not check stream success** — Files are opened and written without validating `is_open()` or write failures. **Solution:** check the stream and report failure with `CGAL_error`, `CGAL_warning`, or an exception. **Response:** Fine for debugging.

13. **Face border properties are never initialized** — `make_initial_feature_graph()` initializes border status for vertices and edges but not faces. **Solution:** set face border status explicitly, usually `false`, or remove border state from face elements. **Response:** Faces cannot be on border in the feature graph.

14. **“Ithmus” is misspelled and appears in internal names/comments** — It is “isthmus” and appears as `ithmus`, `insthmus`, and TODO text. **Solution:** rename variables/comments to `isthmus`. **Response:** Done.

15. **Vector property maps assume stable vertex indices** — `boost::make_vector_property_map` plus graph removals makes cache validity harder to prove. **Solution:** prefer an associative property map for construction vertices or rebuild/validate indices after removals. **Response:** Vertices are only deleted, in this context, I believe this will create holes in the vector, but it stays consistent with the graph. For now, I would prefer not to change code that works until the implementation is completed.

16. **The thinning cache may become inconsistent after removals** — `is_element_computed_cache` is vector-indexed while list vertices are removed. **Solution:** use an associative descriptor-keyed cache or a graph type that explicitly guarantees stable index semantics. **Response:** Same as 15.

17. **Neighborhood connectivity has undefined graph assumptions** — `is_vertex_locally_connected()` assumes an undirected graph and unambiguously named neighbor vertices, but does not assert this. **Solution:** static_assert/document these requirements and construct an explicit induced neighbor graph. **Response:** Changed doc.

18. **The thinning algorithm mutates maps and lists in a hard-to-prove way** — Iterator erasure, cache invalidation, and search restart are intertwined. **Solution:** replace it with a priority queue loop or add explicit invariants and unit tests. **Response:** unit tests will be added later.

19. **Priority keys discard estimator types** — The algorithm converts priorities to `double` and rejects no NaN/inf values. **Solution:** retain the estimator comparison type and validate finite values with `CGAL::is_finite` if available. **Response:** replaced double to sharpness value.

20. **DGtal interface retains large commented-out alternatives** — [DGtal_interface.h](internal/DGtal_interface.h) contains obsolete implementation sketches and disabled code. **Solution:** remove obsolete blocks and retain only the chosen implementation. **Response:** removed some commented code. Commented code will be removed in the far future, not an issue for now.

21. **Normal regularization reverses normals with unclear state** — `regularize_inversion_on_normals()` increments counts before expansion and can repeatedly flip normals across overlapping patches. **Solution:** define patch cardinality and deterministic seen-state logic, and test inversion on a cube. **Response:** inverted patches are small and far apart, this implementation is fine. No need to complicate code that is note part of the core algorithm.

22. **`LabelledImageToBinaryWrapper` stores pointers to constructor arguments** — The image/label lifetime can be unclear, and the label parameter is referenced indirectly. **Solution:** store `const Image&` and `float` by value or use `std::reference_wrapper` with clear lifetime rules. **Response:** References do not work as it should be copy-constructibe.

23. **Labels are hard-coded to `float`** — DGtal image labels and comparisons are fixed to `float`, while image voxel types may differ. **Solution:** make label type template-dependent or normalize through the image's value type. **Response:** CGAL images method `image.value()` returns a float, so this is fine.

24. **Domain/orientation conversion output lengths are unchecked** — `convert_labelled_image_to_face_graphs()` assumes domains and orientations have matching lengths. **Solution:** require/return a paired record or validate vector sizes before adding faces. **Response:** Clarified output and clear content.

25. **Index types switch between signed and unsigned without checks** — Negative `i,j,k` coordinates are used and later converted to `std::size_t`. **Solution:** keep signed coordinates until bounds pass, then use a checked conversion helper. **Response:** added assertion.

26. **Image-to-mesh conversion lacks pre/postconditions** — It assumes dimensions fit `int`, spacing is usable, and `Euler::add_face()` succeeds. **Solution:** add CGAL pre/postconditions and verify mesh vertex/face counts. **Response:** added assertion.

27. **Surfel face orientation is duplicated and untested** — The swap branches in [convert_image_to_surface_mesh.h](internal/convert_image_to_surface_mesh.h) depend on axis/value conventions. **Solution:** factor one explicit oriented-quadrilateral helper and test x/y/z plus both label orders. **Response:** Vertex order is not important for the algorithm, it is only for display.

28. **Euler operations header is not included explicitly** — [convert_image_to_surface_mesh.h](internal/convert_image_to_surface_mesh.h) calls `CGAL::Euler::add_face()` but only includes `CGAL/Surface_mesh.h`. **Solution:** include the exact Euler header or use `Surface_mesh::add_face()` directly. **Response:** Fixed by issue 7.

29. **The feature graph output type is underspecified** — `Feature_graph<K>` stores only points and exposes no public sharpness, dimension, source-element, or border data. **Solution:** define a public vertex bundle or property-map conversion and document it. **Response:** The final output should contain only points. Other values are just for future steps.

## Low / style / build

30. **The example has duplicate/stale includes** — [feature_graph_from_segmented_image.cpp](../../../examples/Feature_graph/feature_graph_from_segmented_image.cpp) includes `CGAL/Image_3.h` and `Exact_predicates_inexact_constructions_kernel.h` twice. **Solution:** remove duplicates and old commented-out includes. **Response:** removed.

31. **Misspellings occur in comments and identifiers** — Examples include “Dubbuging”, “Freature_graph_element”, “assignes”, “constains”, “abherant”, “ithmus”, and “mininmum”. **Solution:** rename/fix these before public review. **Response:** done.

32. **Debug filenames mutate the caller name implicitly** — `_points.ply` and `_edges.polylines.txt` are appended automatically. **Solution:** provide separate filename parameters or a debug-output context/path object. **Response:** Fine for debugging.

33. **Debug macros write fixed files in the working directory** — `CGAL_FEATURE_GRAPH_DEBUG_OUTPUT_FILES` is enabled in the example and writes relative filenames. **Solution:** use an explicit debug named parameter/output directory and keep normal examples feature-free. **Response:** Fine for debugging.

34. **CMake creates a DGtal support target manually** — [CMakeLists.txt](../../../examples/Feature_graph/CMakeLists.txt) creates `CGAL::DGTAL_support` directly. **Solution:** provide a package-level standard CGAL support target and reuse it across examples/tests. **Response:** Cmake specificity will be handled after.

35. **CMake DGtal state can become inconsistent** — The cache option, macro, and target are set in partially independent branches. **Solution:** gate the option, macro, and target creation in one coherent block. **Response:** Cmake specificity will be handled after.

36. **The example explicitly says it is non-functioning** — The public example is described as temporary, and the detector has no output. **Solution:** mark it experimental or replace it with a functioning example once implementation is complete. **Response:** marked experimental.

37. **No package tests are visible** — There are no tests for conversion, orientation, thinning, indexing, or public named parameters. **Solution:** add focused unit tests for small cubes and tetrahedra and compile tests for public headers. **Response:** Tests are for later.

38. **Documentation contains spelling and API mismatches** — “mininmum_energy_delta”, isthmus spelling, and surface/polygon terminology need cleanup. **Solution:** run Doxygen and audit names against implementation.  **Response:** duplicated issue.

39. **Stubbed surface estimator headers are advertised in docs** — [Surface_normal_estimator.h](Surface_normal_estimator.h) and [Surface_sharpness_estimator.h](Surface_sharpness_estimator.h) have no implementation but appear as models in concepts. **Solution:** implement them or remove their documented model claims.  **Response:** Implementation will come later.

40. **The build verifies only the example target** — The example target builds, but no public API/named-parameter or DGtal-disabled coverage exists. **Solution:** add compile-only tests for both macro states and Debug/Release configurations. **Response:** Implementation will come later.