
include_directories(${SEAL_INCLUDE_DIRS})
include_directories(${SEAL_INCLUDE_DIRS2})
include_directories(${SEAL_INCLUDE_DIRS3})
# include_directories(${SEAL_INCLUDE_DIRS4})

add_library(SEAL STATIC IMPORTED)
set_target_properties(SEAL PROPERTIES IMPORTED_LOCATION "${SEAL_LIBRARIES}")

set(CMAKE_CXX_STANDARD 17)

add_executable (runSEAL runSEAL.cpp)
target_link_libraries(runSEAL SEAL -pthread -lz)

