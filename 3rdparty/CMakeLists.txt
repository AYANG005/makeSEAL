INCLUDE(ExternalProject)

ExternalProject_Add(
        build-seal
        URL ${CMAKE_CURRENT_SOURCE_DIR}/v4.1.1.tar.gz
        # GIT_REPOSITORY https://github.com/microsoft/SEAL.git
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

ExternalProject_Get_Property(build-seal INSTALL_DIR)
set(SEAL_INCLUDE_DIRS ${INSTALL_DIR}/include CACHE PATH "SEAL include directory")
set(SEAL_INCLUDE_DIRS2 ${INSTALL_DIR}/include/SEAL-4.1 CACHE PATH "SEAL include directory")
set(SEAL_INCLUDE_DIRS3 ${INSTALL_DIR}/include/SEAL-4.1/seal CACHE PATH "SEAL include directory")
# set(SEAL_INCLUDE_DIRS4 ${INSTALL_DIR}/include/SEAL-4.1/gsl CACHE PATH "SEAL include directory")
set(SEAL_LIBRARIES "${INSTALL_DIR}/lib/libseal-4.1.a" CACHE PATH "SEAL library path")

