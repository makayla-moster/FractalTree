# CMAKE generated file: DO NOT EDIT!
# Generated by "MSYS Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "/C/Program Files/CMake/bin/cmake.exe"

# The command to remove a file.
RM = "/C/Program Files/CMake/bin/cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/splitview.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/splitview.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/splitview.dir/flags.make

examples/CMakeFiles/splitview.dir/splitview.c.obj: examples/CMakeFiles/splitview.dir/flags.make
examples/CMakeFiles/splitview.dir/splitview.c.obj: examples/CMakeFiles/splitview.dir/includes_C.rsp
examples/CMakeFiles/splitview.dir/splitview.c.obj: ../examples/splitview.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/CMakeFiles/splitview.dir/splitview.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/splitview.dir/splitview.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/examples/splitview.c

examples/CMakeFiles/splitview.dir/splitview.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/splitview.dir/splitview.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/examples/splitview.c > CMakeFiles/splitview.dir/splitview.c.i

examples/CMakeFiles/splitview.dir/splitview.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/splitview.dir/splitview.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/examples/splitview.c -o CMakeFiles/splitview.dir/splitview.c.s

examples/CMakeFiles/splitview.dir/splitview.c.obj.requires:

.PHONY : examples/CMakeFiles/splitview.dir/splitview.c.obj.requires

examples/CMakeFiles/splitview.dir/splitview.c.obj.provides: examples/CMakeFiles/splitview.dir/splitview.c.obj.requires
	$(MAKE) -f examples/CMakeFiles/splitview.dir/build.make examples/CMakeFiles/splitview.dir/splitview.c.obj.provides.build
.PHONY : examples/CMakeFiles/splitview.dir/splitview.c.obj.provides

examples/CMakeFiles/splitview.dir/splitview.c.obj.provides.build: examples/CMakeFiles/splitview.dir/splitview.c.obj


examples/CMakeFiles/splitview.dir/glfw.rc.obj: examples/CMakeFiles/splitview.dir/flags.make
examples/CMakeFiles/splitview.dir/glfw.rc.obj: ../examples/glfw.rc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building RC object examples/CMakeFiles/splitview.dir/glfw.rc.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/windres.exe -O coff $(RC_DEFINES) $(RC_INCLUDES) $(RC_FLAGS) /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/examples/glfw.rc CMakeFiles/splitview.dir/glfw.rc.obj

examples/CMakeFiles/splitview.dir/glfw.rc.obj.requires:

.PHONY : examples/CMakeFiles/splitview.dir/glfw.rc.obj.requires

examples/CMakeFiles/splitview.dir/glfw.rc.obj.provides: examples/CMakeFiles/splitview.dir/glfw.rc.obj.requires
	$(MAKE) -f examples/CMakeFiles/splitview.dir/build.make examples/CMakeFiles/splitview.dir/glfw.rc.obj.provides.build
.PHONY : examples/CMakeFiles/splitview.dir/glfw.rc.obj.provides

examples/CMakeFiles/splitview.dir/glfw.rc.obj.provides.build: examples/CMakeFiles/splitview.dir/glfw.rc.obj


examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj: examples/CMakeFiles/splitview.dir/flags.make
examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj: examples/CMakeFiles/splitview.dir/includes_C.rsp
examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj: ../deps/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/splitview.dir/__/deps/glad.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c

examples/CMakeFiles/splitview.dir/__/deps/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/splitview.dir/__/deps/glad.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c > CMakeFiles/splitview.dir/__/deps/glad.c.i

examples/CMakeFiles/splitview.dir/__/deps/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/splitview.dir/__/deps/glad.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c -o CMakeFiles/splitview.dir/__/deps/glad.c.s

examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.requires:

.PHONY : examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.requires

examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.provides: examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.requires
	$(MAKE) -f examples/CMakeFiles/splitview.dir/build.make examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.provides.build
.PHONY : examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.provides

examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.provides.build: examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj


# Object files for target splitview
splitview_OBJECTS = \
"CMakeFiles/splitview.dir/splitview.c.obj" \
"CMakeFiles/splitview.dir/glfw.rc.obj" \
"CMakeFiles/splitview.dir/__/deps/glad.c.obj"

# External object files for target splitview
splitview_EXTERNAL_OBJECTS =

examples/splitview.exe: examples/CMakeFiles/splitview.dir/splitview.c.obj
examples/splitview.exe: examples/CMakeFiles/splitview.dir/glfw.rc.obj
examples/splitview.exe: examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj
examples/splitview.exe: examples/CMakeFiles/splitview.dir/build.make
examples/splitview.exe: src/libglfw3dll.a
examples/splitview.exe: examples/CMakeFiles/splitview.dir/linklibs.rsp
examples/splitview.exe: examples/CMakeFiles/splitview.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable splitview.exe"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && "/C/Program Files/CMake/bin/cmake.exe" -E remove -f CMakeFiles/splitview.dir/objects.a
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/ar.exe cr CMakeFiles/splitview.dir/objects.a @CMakeFiles/splitview.dir/objects1.rsp
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && /C/MinGW/bin/gcc.exe   -mwindows -Wl,--whole-archive CMakeFiles/splitview.dir/objects.a -Wl,--no-whole-archive  -o splitview.exe -Wl,--out-implib,libsplitview.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/splitview.dir/linklibs.rsp

# Rule to build all files generated by this target.
examples/CMakeFiles/splitview.dir/build: examples/splitview.exe

.PHONY : examples/CMakeFiles/splitview.dir/build

examples/CMakeFiles/splitview.dir/requires: examples/CMakeFiles/splitview.dir/splitview.c.obj.requires
examples/CMakeFiles/splitview.dir/requires: examples/CMakeFiles/splitview.dir/glfw.rc.obj.requires
examples/CMakeFiles/splitview.dir/requires: examples/CMakeFiles/splitview.dir/__/deps/glad.c.obj.requires

.PHONY : examples/CMakeFiles/splitview.dir/requires

examples/CMakeFiles/splitview.dir/clean:
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/splitview.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/splitview.dir/clean

examples/CMakeFiles/splitview.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MSYS Makefiles" /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1 /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/examples /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/examples/CMakeFiles/splitview.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/splitview.dir/depend

