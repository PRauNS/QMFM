cmake_minimum_required(VERSION 3.26)
project(QMFM)

# --------------------------------------------------
# Compiler settings
# --------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Optimization flags
set(CMAKE_BUILD_TYPE Debug) #when debugging
#add_compile_options(-O3 -O2) #for production runs

# Preprocessor definitions
add_compile_definitions(
        OLD_HEADER_FILENAME
        HDF_NO_NAMESPACE
        NO_STATIC_CAST
)

# --------------------------------------------------
# Find required packages
# --------------------------------------------------

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

find_library(FFTW3_LIB fftw3 REQUIRED
		HINTS /cm/shared/apps/fftw/openmpi/gcc/64/3.3.10/lib)
find_library(FFTW3_OMP_LIB fftw3_omp)

if(NOT FFTW3_LIB)
    message(FATAL_ERROR "FFTW3 not found")
endif()

find_library(SHTNS_LIB shtns_omp
        HINTS /usr/local/lib
)

# --------------------------------------------------
# Source files
# --------------------------------------------------

set(SOURCES
        main.cpp
        common.cpp
        initial_conditions.cpp
        boundary_conditions.cpp
        conservation_checks.cpp
        microphysics.cpp
        microphysics_Bdep.cpp
        field_evolution.cpp
)

add_executable(main ${SOURCES})

# --------------------------------------------------
# Include directories
# --------------------------------------------------

target_include_directories(main PRIVATE
        ${MPI_CXX_INCLUDE_PATH}
        ${HDF5_INCLUDE_DIRS}
)

# --------------------------------------------------
# Link libraries
# --------------------------------------------------

target_link_libraries(main PRIVATE
        MPI::MPI_CXX
        ${HDF5_LIBRARIES}
        GSL::gsl
        GSL::gslcblas
        #FFTW3::fftw3
        #FFTW3::fftw3_omp
        OpenMP::OpenMP_CXX
        ${FFTW3_LIB}
)

if(FFTW3_OMP_LIB)
    target_link_libraries(main PRIVATE ${FFTW3_OMP_LIB})
endif()

# Link SHTns if found
if(SHTNS_LIB)
    target_link_libraries(main PRIVATE ${SHTNS_LIB})
endif()
