# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_scream)

# if (USE_SAMXX  OR COMP_NAMES MATCHES ".*scream.*") 
    # Include machine file here
    message("Found scream component")
    include(${CMAKE_SOURCE_DIR}/scream/cmake/machine-files/${MACH}.cmake)
    add_subdirectory("scream")
    set(SCREAM_SOURCE_DIR ${CMAKE_SOURCE_DIR}/scream/src CACHE STRING "Scream Source Dir")
    set(SCREAM_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/scream CACHE STRING "Scream Binary Dir")
# endif()

endfunction(build_scream)
