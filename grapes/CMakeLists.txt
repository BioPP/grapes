# CMake script for Grapes
# Author: Julien Dutheil
# Created: 30/05/2018

add_executable(grapes Grapes.cpp)
add_executable(multi_grapes MultiGrapes.cpp)
set (grapes-targets grapes multi_grapes)

foreach (target ${grapes-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC})
    target_link_libraries (${target} ${LIBS})
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED})
    target_link_libraries (${target} ${LIBS})
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${grapes-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
