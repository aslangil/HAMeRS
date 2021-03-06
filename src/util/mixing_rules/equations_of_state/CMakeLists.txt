# Define a variable equations_of_state_source_files containing
# a list of the source files for the equations_of_state library
set(equations_of_state_source_files
    EquationOfStateMixingRulesManager.cpp
    ./EquationOfStateMixingRules.cpp
    ./ideal_gas/EquationOfStateIdealGas.cpp
    ./ideal_gas/EquationOfStateMixingRulesIdealGas.cpp
    ./stiffened_gas/EquationOfStateStiffenedGas.cpp
    ./stiffened_gas/EquationOfStateMixingRulesStiffenedGas.cpp
)

# Create a library called equations_of_state which includes the 
# source files defined in equations_of_state_source_files
add_library(equations_of_state ${equations_of_state_source_files})

TARGET_LINK_LIBRARIES(equations_of_state
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  ${HDF5_C_LIBRARIES} gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(equations_of_state z)
endif()
