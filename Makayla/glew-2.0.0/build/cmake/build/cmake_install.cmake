# Install script for directory: C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/glew")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/lib/libglew32.dll.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/bin/glew32.dll")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/glew32.dll" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/glew32.dll")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "C:/MinGW/bin/strip.exe" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/glew32.dll")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/lib/libglew32.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/bin/glewinfo.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/glewinfo.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/glewinfo.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "C:/MinGW/bin/strip.exe" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/glewinfo.exe")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/bin/visualinfo.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/visualinfo.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/visualinfo.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "C:/MinGW/bin/strip.exe" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/visualinfo.exe")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/../../glew.pc")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES
      )
  endif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg]|[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/GL" TYPE FILE FILES
    "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/../../include/GL/wglew.h"
    "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/../../include/GL/glew.h"
    "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/../../include/GL/glxew.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/glew/glew-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/glew/glew-targets.cmake"
         "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/CMakeFiles/Export/lib/cmake/glew/glew-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/glew/glew-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/glew/glew-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/glew" TYPE FILE FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/CMakeFiles/Export/lib/cmake/glew/glew-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/glew" TYPE FILE FILES "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/CMakeFiles/Export/lib/cmake/glew/glew-targets-release.cmake")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/glew" TYPE FILE FILES
    "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/glew-config.cmake"
    "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/CopyImportedTargetProperties.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/MinGW/msys/1.0/home/Makayla/glew-2.0.0/build/cmake/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
