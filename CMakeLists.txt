cmake_minimum_required (VERSION 2.8)
project (Turbine)

FIND_PACKAGE(VTK REQUIRED)
include(${VTK_USE_FILE})

set(OSUFLOW_DIR "/usr/local/osuflow" CACHE PATH "OSUFlow dir")
set(OSUFLOW_INCLUDE_DIR "${OSUFLOW_DIR}/include")
set(OSUFLOW_LIB_DIR "${OSUFLOW_DIR}/lib")

# OpenMP
#find_package(OpenMP)
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


# convertVTK

add_executable(convertVTK convertVTK.cpp)
target_link_libraries( convertVTK ${VTK_LIBRARIES} )

## streamline render

include_directories(${OSUFLOW_INCLUDE_DIR} ${OSUFLOW_INCLUDE_DIR}/VTK)
link_directories(${OSUFLOW_LIB_DIR})

#add_executable(StreamlineRenderPlot3D StreamlineRenderPlot3D.cpp)
#target_link_libraries( StreamlineRenderPlot3D  OSUFlowVTK OSUFlow_r ${VTK_LIBRARIES})

add_executable(StreamlineRenderUnstructured StreamlineRenderUnstructured.cpp)
target_link_libraries( StreamlineRenderUnstructured  OSUFlowVTK OSUFlow_r ${VTK_LIBRARIES})

add_executable(streamline streamline.cpp)
target_link_libraries( streamline  OSUFlowVTK OSUFlow_r ${VTK_LIBRARIES})

add_executable(computeFTLE computeFTLE.cpp)
target_link_libraries( computeFTLE   OSUFlow_r)

add_executable(regular_sampling regular_sampling.cpp)
target_link_libraries( regular_sampling ${VTK_LIBRARIES} )

add_executable(pmi_importance_2var pmi_importance_2var.cpp)
target_link_libraries( pmi_importance_2var ${VTK_LIBRARIES} )

#add_definitions(-std=gnu++0x)

add_executable(stat stat.cpp)
target_link_libraries( stat ${VTK_LIBRARIES} )
