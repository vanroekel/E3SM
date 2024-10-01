# For this file, fixes non-BFB behavior of stealth feature on pm-cpu with -O2
list(APPEND NOOPT_FILES cice/src/mpi/ice_boundary.F90 cice/src/serial/ice_boundary.F90 cice/src/source/ice_domain.F90 eam/src/physics/cam/gw_common.F90)

set(NOOPT
  eam/src/physics/cam/zm_conv.F90)

set(O0MODELSRC
  eam/src/chemistry/aerosol/dust_sediment_mod.F90
  eam/src/chemistry/modal_aero/modal_aero_convproc.F90
  eam/src/chemistry/utils/modal_aero_calcsize.F90
  eam/src/physics/cam/zm_conv.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS O0MODELSRC)
    e3sm_add_flags("${ITEM}" "-O0 -vector0")
  endforeach()

endif()

if (NOT DEBUG)
  foreach(ITEM IN LISTS NOOPT)
    e3sm_deoptimize_file("${ITEM}")
  endforeach()
endif()




