# message("SUBFOLDER:")

# message(${CMAKE_SOURCE_DIR})
# message(${CMAKE_CURRENT_SOURCE_DIR})
# message(${CMAKE_BINARY_DIR})
# message(${CMAKE_CURRENT_BINARY_DIR})
# message(${PROJECT_NAME})
# message(${CMAKE_PROJECT_NAME})

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Triangulation.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Triangulation1.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Triangulation2.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Polyhedra.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Dimension.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshExport.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Dual.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/CamminoMinimo.cpp)

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_headers ${src_headers} PARENT_SCOPE)
set(src_sources ${src_sources} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)
