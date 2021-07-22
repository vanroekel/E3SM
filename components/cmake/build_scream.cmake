# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_scream)

  # Include machine file here
  message("Found scream component")
  include(${CMAKE_SOURCE_DIR}/scream/cmake/machine-files/${MACH}.cmake)
  add_subdirectory("scream")
# set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)

endfunction(build_scream)
