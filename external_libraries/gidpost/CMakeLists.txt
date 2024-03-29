set(CMAKE_INCLUDE_CURRENT_DIR ON)

## generate variables with the sources
set( GIDPOST_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/source/gidpost.c
 	${CMAKE_CURRENT_SOURCE_DIR}/source/gidpostInt.c
 	${CMAKE_CURRENT_SOURCE_DIR}/source/gidpostHash.c
 	${CMAKE_CURRENT_SOURCE_DIR}/source/hashtab.c
 	${CMAKE_CURRENT_SOURCE_DIR}/source/recycle.c
 	${CMAKE_CURRENT_SOURCE_DIR}/source/lookupa.c
     )
##gidpostfor.c gidpostforAPI.c

add_definitions( -DUSE_CONST )
add_definitions( -w )

###############################################################
add_library(gidpost STATIC ${GIDPOST_SOURCES})

# install(TARGETS gidpost DESTINATION libs )

if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_DEPRECATE)
endif (MSVC)

if (HDF5)
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(gidpost ${HDF5_LIBRARIES})
endif (HDF5)

target_link_libraries(gidpost zlibstatic)

find_package( Threads )
if ( Threads_FOUND )
  target_link_libraries(gidpost ${CMAKE_THREAD_LIBS_INIT})
endif (Threads_FOUND)

if (UNIX)
  target_link_libraries(gidpost m)
endif (UNIX)

