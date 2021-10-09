set(UGENE_PLUGIN_NAME biostruct3d_view)

file(GLOB SRCS src/*.cpp src/*.c src/*.h src/gl2ps/*)

include(../../Plugin.cmake)

if(APPLE)
    target_link_libraries(biostruct3d_view "-framework OpenGL")
elseif(WIN32)
    find_package(OpenGL REQUIRED)
    include_directories( ${OPENGL_INCLUDE_DIRS} )
    target_link_libraries(biostruct3d_view ${OPENGL_LIBRARIES})
endif()
