project(MantidParaViewPeaksReader)
include_directories(SYSTEM ${PARAVIEW_INCLUDE_DIRS})
add_paraview_plugin(MantidParaViewPeaksReaderSMPlugin
                    "1.0"
                    SERVER_MANAGER_XML
                    PeaksReader.xml
                    SERVER_MANAGER_SOURCES
                    vtkPeaksReader.cxx)
set_pvplugin_properties(MantidParaViewPeaksReaderSMPlugin QT_VERSION 4)

# Add to the 'VatesParaViewPlugins' group in VS
set_property(TARGET MantidParaViewPeaksReaderSMPlugin
             PROPERTY FOLDER "MantidVatesParaViewPlugins")

target_link_libraries(MantidParaViewPeaksReaderSMPlugin
                      LINK_PRIVATE
                      ${CORE_MANTIDLIBS}
                      DataObjects
                      VatesAPI
                      ${POCO_LIBRARIES}
                      ${Boost_LIBRARIES})

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(
    MantidParaViewPeaksReaderSMPlugin
    PROPERTIES
      INSTALL_RPATH
      "@loader_path/../../../Contents/Frameworks;@loader_path/../../../Contents/MacOS"
    )
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  set_target_properties(MantidParaViewPeaksReaderSMPlugin
                        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../../${LIB_DIR}")
endif()

install_pvplugin(MantidParaViewPeaksReaderSMPlugin QT_VERSION 4)
