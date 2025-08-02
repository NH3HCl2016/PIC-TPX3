
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(TPX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Find OpenMP (required package for parallel processing)
#
option(USE_MT "Build the project using Geant4's built-in multi-threaded simulation" OFF)
if (USE_MT_THREAD)
  find_package(OpenMP REQUIRED)
endif()

# if(OpenMP_CXX_FOUND)
#   set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OPENMP_C_FLAGS})
#   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}, ${OPENMP_CXX_FLAGS})
# endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
                    ${OpenMP_CXX_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the related libraries (Geant4, ROOT and OpenMP if necessary)
#
add_executable(TPX3 TPX3.cc ${sources} ${headers})
target_link_libraries(TPX3 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

if (USE_MT)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(OpenMP::OpenMP_CXX)
  endif()
  add_definitions(-DTPXMT)
endif()

option(USE_MT_THREAD "Build the project using OpenMP-based multi-thread processing in each thread" OFF)
if (USE_MT_THREAD)
  add_definitions(-DTPXInThreadParallel)
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TPX. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEEC_SCRIPTS
  init_vis.mac
  vis.mac
  )

foreach(_script ${EXAMPLEEC_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(TPX DEPENDS TPX)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TPX3 DESTINATION bin)


