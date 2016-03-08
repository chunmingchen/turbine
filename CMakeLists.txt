cmake_minimum_required (VERSION 2.8)
project (Turbine)

option (USE_PARAVIEW "use paraview")

if (USE_PARAVIEW)
  find_package(Paraview)
  include (${PARAVIEW_USE_FILE})
else ()
  FIND_PACKAGE(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif ()

find_package(osuflow)

#set(OSUFLOW_DIR "/usr/local/osuflow" CACHE PATH "OSUFlow dir")
set(JCLIB_INSTALL_DIR "${CMAKE_SOURCE_DIR}/lib/jclib/" CACHE PATH "JCLib installed dir")

#set(OSUFLOW_INCLUDE_DIR "${OSUFLOW_DIR}/include/osuflow")
#set(OSUFLOW_LIB_DIR "${OSUFLOW_DIR}/lib/osuflow")

include_directories(${JCLIB_INSTALL_DIR}/include)
link_directories(${JCLIB_INSTALL_DIR}/lib/jclib)

set(DATA_DIR "/data" CACHE PATH "Directory to PLOT3D data")
add_definitions(-DDATA_DIR=\"${DATA_DIR}/\")

# OpenMP
find_package(OpenMP)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

## GLM
set(GLM_DIR "/usr/lib" CACHE PATH "Directory to GLM")
include_directories(${GLM_DIR})


# convertVTK

add_executable(convertVTK convertVTK.cpp)
target_link_libraries( convertVTK ${VTK_LIBRARIES} )

## streamline render

#include_directories(${OSUFLOW_INCLUDE_DIR} ${OSUFLOW_INCLUDE_DIR}/VTK
#    ${JCLIB_DIR})
#link_directories(${OSUFLOW_LIB_DIR} ${JCLIB_DIR})

#add_executable(StreamlineRenderPlot3D StreamlineRenderPlot3D.cpp)
#target_link_libraries( StreamlineRenderPlot3D  OSUFlowVTK OSUFlow ${VTK_LIBRARIES})

add_executable(StreamlineRenderUnstructured StreamlineRenderUnstructured.cpp)
target_link_libraries( StreamlineRenderUnstructured ${osuflow_LIBRARIES} ${VTK_LIBRARIES})

#add_executable(streamline streamline.cpp)
#target_link_libraries( streamline  ${osuflow_LIBRARIES} ${VTK_LIBRARIES})

add_executable(computeFTLE computeFTLE.cpp)
target_link_libraries( computeFTLE   ${osuflow_LIBRARIES})

add_executable(regular_sampling regular_sampling.cpp vtkImageProbeFilter.cxx vtkImageCompositeDataProbeFilter.cxx )
target_link_libraries( regular_sampling ${VTK_LIBRARIES} libjclib.a)

add_executable(regular_sampling_1passage regular_sampling_1passage.cpp vtkImageProbeFilter.cxx)
target_link_libraries( regular_sampling_1passage ${VTK_LIBRARIES} )

add_executable(pmi_importance_2var pmi_importance_2var.cpp)
target_link_libraries( pmi_importance_2var ${VTK_LIBRARIES} )

#add_definitions(-std=gnu++0x)

add_executable(stat stat.cpp)
target_link_libraries( stat ${VTK_LIBRARIES} )

add_executable(stat_raw stat_raw.cpp)
target_link_libraries( stat_raw ${VTK_LIBRARIES} jclib ${osuflow_LIBRARIES})

add_executable(surfaces2volume surfaces2volume.cpp)
target_link_libraries( surfaces2volume ${VTK_LIBRARIES} jclib ${osuflow_LIBRARIES})

add_executable(flowsepcount flowsepcount.cpp)
target_link_libraries( flowsepcount ${VTK_LIBRARIES} jclib )

add_executable(VortexDetector VortexDetector.C)
target_link_libraries( VortexDetector  ${osuflow_LIBRARIES} ${VTK_LIBRARIES} jclib)

#add_executable(test_vtkopengl test_vtkopengl.cpp)
#target_link_libraries( test_vtkopengl  ${osuflow_LIBRARIES} ${VTK_LIBRARIES} )

add_executable(vts_vortex_merge vts_vortex_merge.cpp)
target_link_libraries( vts_vortex_merge  ${VTK_LIBRARIES} )

#add_executable(pressure_probe pressure_probe.cpp)
#target_link_libraries( pressure_probe ${VTK_LIBRARIES} )

add_executable(extractVariables extractvariables.cpp)
target_link_libraries( extractVariables  ${VTK_LIBRARIES} jclib)

add_executable(extractVariables_1passage extractvariables_1passage.cpp)
target_link_libraries( extractVariables_1passage  ${VTK_LIBRARIES} jclib)

add_executable(convertUnstructured convertunstructured.cpp)
target_link_libraries( convertUnstructured  ${VTK_LIBRARIES} )

add_executable(gmm_emd gmm_emd.cpp)
target_link_libraries( gmm_emd  ${VTK_LIBRARIES} )

################################################################33
# Paraview
#ADD_PARAVIEW_PLUGIN(LoadTurbine "1.0"
#  SERVER_MANAGER_XML CellDerivatives.xml)


