set(Eigen_HEADERS AdolcForward BVH IterativeSolvers MatrixFunctions MoreVectorization AutoDiff AlignedVector3 Polynomials
                  CholmodSupport FFT NonLinearOptimization SparseExtra SuperLUSupport UmfPackSupport IterativeSolvers
                  NumericalDiff Skyline MPRealSupport OpenGLSupport
   )

install(FILES
  ${Eigen_HEADERS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/unsupported/Eigen COMPONENT Devel
  )

add_subdirectory(src)
