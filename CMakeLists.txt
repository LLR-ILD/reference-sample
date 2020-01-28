###############################################################################
# cmake file for building the reference-sample analysis package.
# @author Remi Ete, DESY (original)
# @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
###############################################################################

# Project name.
PROJECT( reference-sample )

# Project version.
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ##############################################################
FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# Load default settings from ILCSOFT_CMAKE_MODULES.
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # Minimum required Marlin version.
INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( DD4hep )
INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} )
ADD_DEFINITIONS( ${DD4hep_DEFINITIONS} )

FIND_PACKAGE( MarlinUtil 1.0 REQUIRED ) # Minimum required MarlinUtil version.
INCLUDE_DIRECTORIES( SYSTEM ${MarlinUtil_INCLUDE_DIRS} )
LINK_LIBRARIES( ${MarlinUtil_LIBRARIES} )
ADD_DEFINITIONS( ${MarlinUtil_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED COMPONENTS GenVector)
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

### LIBRARY ###################################################################
# Register all processors/ packages with cmake.
MACRO( ADD_PROCESSOR _input_dir )
    MESSAGE( STATUS "Including your processor: ${_input_dir}" )
    INCLUDE_DIRECTORIES( BEFORE ${_input_dir}/include )
    AUX_SOURCE_DIRECTORY( ${_input_dir}/src project_cxx_srcs )
ENDMACRO()

# -----------------------------------------------------------------------------
# Add the packages you want to compile here.
ADD_PROCESSOR( ./processors/shared_functionality )
ADD_PROCESSOR( ./processors/check_isolated_lepton_tagging )
ADD_PROCESSOR( ./processors/check_relation )
ADD_PROCESSOR( ./processors/example )
ADD_PROCESSOR( ./processors/draft )
ADD_PROCESSOR( ./processors/mc_as_reconstructed )
ADD_PROCESSOR( ./processors/skip_empty_events )
ADD_PROCESSOR( ./processors/slcio_to_root )
ADD_PROCESSOR( ./processors/split_off_z )
ADD_PROCESSOR( ./processors/tau_cones )

# Build the project library.
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${project_cxx_srcs} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# Display some variables and write them to cache.
DISPLAY_STD_VARIABLES()