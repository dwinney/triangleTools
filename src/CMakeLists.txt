# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB_RECURSE SRC "*.cpp" "*.c")

include_directories(.)
include_directories(old)
include_directories(cubature)
include_directories(plotting)
include_directories(${LOOPTOOLS_INC})

# Build LoopTools
add_library(LOOPTOOLS SHARED ${STATIC_LOOPTOOLS})
set_target_properties(LOOPTOOLS PROPERTIES LINKER_LANGUAGE Fortran )
target_link_libraries(LOOPTOOLS PUBLIC
        -Wl,--whole-archive
        ${STATIC_LOOPTOOLS}
        -Wl,--no-whole-archive
        )

# Build hadMolee and ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__TRIANGLETOOLS
        plotter.hpp
        plot.hpp
        LINKDEF TRIANGLETOOLSLinkDef.h
        )
add_library( TRIANGLETOOLS SHARED ${SRC} G__TRIANGLETOOLS.cxx )
set_target_properties( TRIANGLETOOLS PROPERTIES LINKER_LANGUAGE CXX )
target_link_libraries( TRIANGLETOOLS ROOT::Core ROOT::MathMore )
target_link_libraries( TRIANGLETOOLS Boost::boost)
set_target_properties( TRIANGLETOOLS PROPERTIES LINKER_LANGUAGE Fortran )
target_link_libraries( TRIANGLETOOLS LOOPTOOLS )

# INSTALLATION 
install(TARGETS LOOPTOOLS
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

install(TARGETS TRIANGLETOOLS
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libTRIANGLETOOLS_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable( triangleTools "cling/triangleTools.cpp")
set_target_properties( triangleTools PROPERTIES LINKER_LANGUAGE CXX )
target_link_libraries( triangleTools ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( triangleTools Boost::boost)
target_link_libraries( triangleTools TRIANGLETOOLS)
set_target_properties( triangleTools PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries( triangleTools LOOPTOOLS )
