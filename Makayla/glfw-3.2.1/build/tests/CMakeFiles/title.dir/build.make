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
include tests/CMakeFiles/title.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/title.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/title.dir/flags.make

tests/CMakeFiles/title.dir/title.c.obj: tests/CMakeFiles/title.dir/flags.make
tests/CMakeFiles/title.dir/title.c.obj: tests/CMakeFiles/title.dir/includes_C.rsp
tests/CMakeFiles/title.dir/title.c.obj: ../tests/title.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object tests/CMakeFiles/title.dir/title.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/title.dir/title.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/title.c

tests/CMakeFiles/title.dir/title.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/title.dir/title.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/title.c > CMakeFiles/title.dir/title.c.i

tests/CMakeFiles/title.dir/title.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/title.dir/title.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/title.c -o CMakeFiles/title.dir/title.c.s

tests/CMakeFiles/title.dir/title.c.obj.requires:

.PHONY : tests/CMakeFiles/title.dir/title.c.obj.requires

tests/CMakeFiles/title.dir/title.c.obj.provides: tests/CMakeFiles/title.dir/title.c.obj.requires
	$(MAKE) -f tests/CMakeFiles/title.dir/build.make tests/CMakeFiles/title.dir/title.c.obj.provides.build
.PHONY : tests/CMakeFiles/title.dir/title.c.obj.provides

tests/CMakeFiles/title.dir/title.c.obj.provides.build: tests/CMakeFiles/title.dir/title.c.obj


tests/CMakeFiles/title.dir/__/deps/glad.c.obj: tests/CMakeFiles/title.dir/flags.make
tests/CMakeFiles/title.dir/__/deps/glad.c.obj: tests/CMakeFiles/title.dir/includes_C.rsp
tests/CMakeFiles/title.dir/__/deps/glad.c.obj: ../deps/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object tests/CMakeFiles/title.dir/__/deps/glad.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/title.dir/__/deps/glad.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c

tests/CMakeFiles/title.dir/__/deps/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/title.dir/__/deps/glad.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c > CMakeFiles/title.dir/__/deps/glad.c.i

tests/CMakeFiles/title.dir/__/deps/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/title.dir/__/deps/glad.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c -o CMakeFiles/title.dir/__/deps/glad.c.s

tests/CMakeFiles/title.dir/__/deps/glad.c.obj.requires:

.PHONY : tests/CMakeFiles/title.dir/__/deps/glad.c.obj.requires

tests/CMakeFiles/title.dir/__/deps/glad.c.obj.provides: tests/CMakeFiles/title.dir/__/deps/glad.c.obj.requires
	$(MAKE) -f tests/CMakeFiles/title.dir/build.make tests/CMakeFiles/title.dir/__/deps/glad.c.obj.provides.build
.PHONY : tests/CMakeFiles/title.dir/__/deps/glad.c.obj.provides

tests/CMakeFiles/title.dir/__/deps/glad.c.obj.provides.build: tests/CMakeFiles/title.dir/__/deps/glad.c.obj


# Object files for target title
title_OBJECTS = \
"CMakeFiles/title.dir/title.c.obj" \
"CMakeFiles/title.dir/__/deps/glad.c.obj"

# External object files for target title
title_EXTERNAL_OBJECTS =

tests/title.exe: tests/CMakeFiles/title.dir/title.c.obj
tests/title.exe: tests/CMakeFiles/title.dir/__/deps/glad.c.obj
tests/title.exe: tests/CMakeFiles/title.dir/build.make
tests/title.exe: src/libglfw3dll.a
tests/title.exe: tests/CMakeFiles/title.dir/linklibs.rsp
tests/title.exe: tests/CMakeFiles/title.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable title.exe"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && "/C/Program Files/CMake/bin/cmake.exe" -E remove -f CMakeFiles/title.dir/objects.a
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/ar.exe cr CMakeFiles/title.dir/objects.a @CMakeFiles/title.dir/objects1.rsp
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe   -mwindows -Wl,--whole-archive CMakeFiles/title.dir/objects.a -Wl,--no-whole-archive  -o title.exe -Wl,--out-implib,libtitle.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/title.dir/linklibs.rsp

# Rule to build all files generated by this target.
tests/CMakeFiles/title.dir/build: tests/title.exe

.PHONY : tests/CMakeFiles/title.dir/build

tests/CMakeFiles/title.dir/requires: tests/CMakeFiles/title.dir/title.c.obj.requires
tests/CMakeFiles/title.dir/requires: tests/CMakeFiles/title.dir/__/deps/glad.c.obj.requires

.PHONY : tests/CMakeFiles/title.dir/requires

tests/CMakeFiles/title.dir/clean:
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/title.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/title.dir/clean

tests/CMakeFiles/title.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MSYS Makefiles" /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1 /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests/CMakeFiles/title.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/title.dir/depend

