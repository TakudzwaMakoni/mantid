set(SRC_FILES
    src/CatalogPublish.cpp
    src/CatalogAlgorithmHelper.cpp
    src/CatalogDownloadDataFiles.cpp
    src/CatalogGetDataFiles.cpp
    src/CatalogGetDataSets.cpp
    src/CatalogKeepAlive.cpp
    src/CatalogListInstruments.cpp
    src/CatalogListInvestigationTypes.cpp
    src/CatalogLogin.cpp
    src/CatalogLogout.cpp
    src/CatalogMyDataSearch.cpp
    src/CatalogSearch.cpp
    src/CatalogSearchParam.cpp
    src/GSoap.cpp
    src/ICat3/ICat3GSoapGenerated.cpp
    src/ICat3/ICat3ErrorHandling.cpp
    src/ICat3/ICat3Catalog.cpp
    src/ICat3/ICat3Helper.cpp
    src/ICat4/ICat4GSoapGenerated.cpp
    src/ICat4/ICat4Catalog.cpp)

set(SRC_UNITY_IGNORE_FILES)

set(INC_FILES
    inc/MantidICat/CatalogPublish.h
    inc/MantidICat/CatalogAlgorithmHelper.h
    inc/MantidICat/CatalogDownloadDataFiles.h
    inc/MantidICat/CatalogGetDataFiles.h
    inc/MantidICat/CatalogGetDataSets.h
    inc/MantidICat/CatalogKeepAlive.h
    inc/MantidICat/CatalogListInstruments.h
    inc/MantidICat/CatalogListInvestigationTypes.h
    inc/MantidICat/CatalogLogin.h
    inc/MantidICat/CatalogLogout.h
    inc/MantidICat/CatalogMyDataSearch.h
    inc/MantidICat/CatalogSearch.h
    inc/MantidICat/CatalogSearchParam.h
    inc/MantidICat/GSoap/soapserializersH.h
    inc/MantidICat/GSoap/soapserializersStub.h
    inc/MantidICat/GSoap/stdsoap2.h
    inc/MantidICat/ICat3/GSoapGenerated/ICat3H.h
    inc/MantidICat/ICat3/GSoapGenerated/ICat3ICATPortBindingProxy.h
    inc/MantidICat/ICat3/GSoapGenerated/ICat3Stub.h
    inc/MantidICat/ICat3/ICat3ErrorHandling.h
    inc/MantidICat/ICat3/ICat3Catalog.h
    inc/MantidICat/ICat3/ICat3Helper.h
    inc/MantidICat/ICat4/GSoapGenerated/ICat4H.h
    inc/MantidICat/ICat4/GSoapGenerated/ICat4ICATPortBindingProxy.h
    inc/MantidICat/ICat4/GSoapGenerated/ICat4Stub.h
    inc/MantidICat/ICat4/ICat4Catalog.h)

set(TEST_FILES
    CatalogDownloadDataFilesTest.h
    CatalogGetDataFilesTest.h
    CatalogGetDataSetsTest.h
    CatalogListInstrumentsTest.h
    CatalogListInvestigationTypesTest.h
    CatalogLoginTest.h
    CatalogMyDataSearchTest.h
    CatalogSearchTest.h
    CompositeCatalogTest.h)

if(COVERALLS)
  foreach(loop_var ${SRC_FILES} ${INC_FILES})
    set_property(GLOBAL APPEND
                 PROPERTY COVERAGE_SRCS
                          "${CMAKE_CURRENT_SOURCE_DIR}/${loop_var}")
  endforeach(loop_var)
endif()

if(UNITY_BUILD)
  include(UnityBuild)
  enable_unity_build(ICat SRC_FILES SRC_UNITY_IGNORE_FILES 10)
endif(UNITY_BUILD)

# Add ssl dependency
include_directories(${OPENSSL_INCLUDE_DIR})
add_definitions(-DWITH_OPENSSL -DWITH_NONAMESPACES)

# Add a precompiled header where they are supported
enable_precompiled_headers(inc/MantidICat/PrecompiledHeader.h SRC_FILES)
# Add the target for this directory
add_library(ICat ${SRC_FILES} ${INC_FILES})
# Set the name of the generated library
set_target_properties(ICat
                      PROPERTIES OUTPUT_NAME
                                 MantidICat
                                 COMPILE_DEFINITIONS
                                 IN_MANTID_ICAT)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(ICat
                        PROPERTIES INSTALL_RPATH
                                   "@loader_path/../Contents/MacOS")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  set_target_properties(ICat PROPERTIES INSTALL_RPATH "\$ORIGIN/../${LIB_DIR}")
endif()

# Add to the 'Framework' group in VS
set_property(TARGET ICat PROPERTY FOLDER "MantidFramework")

include_directories(inc)

target_link_libraries(ICat
                      LINK_PRIVATE
                      ${MANTIDLIBS}
                      ${OPENSSL_LIBRARIES}
                      ${JSONCPP_LIBRARIES})

# Add the unit tests directory
add_subdirectory(test)

# Auto-generate exports header
target_include_directories(ICat PUBLIC ${CMAKE_BINARY_DIR}/Framework/ICat)
generate_mantid_export_header(ICat TRUE)

# Installation settings
install(TARGETS ICat ${SYSTEM_PACKAGE_TARGET} DESTINATION ${PLUGINS_DIR})
