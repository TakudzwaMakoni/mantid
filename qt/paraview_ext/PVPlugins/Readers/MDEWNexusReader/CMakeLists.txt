project(MantidParaViewMDEWNexusReader)
include_directories(SYSTEM ${PARAVIEW_INCLUDE_DIRS})
add_paraview_plugin(MantidParaViewMDEWNexusReaderSMPlugin
                    "1.0"
                    SERVER_MANAGER_XML
                    MDEWNexusReader.xml
                    SERVER_MANAGER_SOURCES
                    vtkMDEWNexusReader.cxx)
set_pvplugin_properties(MantidParaViewMDEWNexusReaderSMPlugin QT_VERSION 4)

include_directories(SYSTEM ${QWT5_INCLUDE_DIR})

target_link_libraries(MantidParaViewMDEWNexusReaderSMPlugin
                      LINK_PRIVATE
                      ${vtkjsoncpp_LIBRARIES}
                      ${CORE_MANTIDLIBS}
                      DataObjects
                      VatesAPI
                      ${POCO_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${QWT5_LIBRARIES}
                      Qt4::QtCore)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(
    MantidParaViewMDEWNexusReaderSMPlugin
    PROPERTIES
      INSTALL_RPATH
      "@loader_path/../../../Contents/Frameworks;@loader_path/../../../Contents/MacOS"
    )
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  set_target_properties(MantidParaViewMDEWNexusReaderSMPlugin
                        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../../${LIB_DIR}")
endif()

install_pvplugin(MantidParaViewMDEWNexusReaderSMPlugin QT_VERSION 4)
