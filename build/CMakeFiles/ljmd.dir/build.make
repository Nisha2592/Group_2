# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build"

# Include any dependencies generated for this target.
include CMakeFiles/ljmd.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ljmd.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ljmd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ljmd.dir/flags.make

CMakeFiles/ljmd.dir/src/main.c.o: CMakeFiles/ljmd.dir/flags.make
CMakeFiles/ljmd.dir/src/main.c.o: /home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/src/main.c
CMakeFiles/ljmd.dir/src/main.c.o: CMakeFiles/ljmd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/ljmd.dir/src/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/ljmd.dir/src/main.c.o -MF CMakeFiles/ljmd.dir/src/main.c.o.d -o CMakeFiles/ljmd.dir/src/main.c.o -c "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/src/main.c"

CMakeFiles/ljmd.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/ljmd.dir/src/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/src/main.c" > CMakeFiles/ljmd.dir/src/main.c.i

CMakeFiles/ljmd.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/ljmd.dir/src/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/src/main.c" -o CMakeFiles/ljmd.dir/src/main.c.s

# Object files for target ljmd
ljmd_OBJECTS = \
"CMakeFiles/ljmd.dir/src/main.c.o"

# External object files for target ljmd
ljmd_EXTERNAL_OBJECTS =

ljmd: CMakeFiles/ljmd.dir/src/main.c.o
ljmd: CMakeFiles/ljmd.dir/build.make
ljmd: libmdlib.so
ljmd: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
ljmd: CMakeFiles/ljmd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ljmd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ljmd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ljmd.dir/build: ljmd
.PHONY : CMakeFiles/ljmd.dir/build

CMakeFiles/ljmd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ljmd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ljmd.dir/clean

CMakeFiles/ljmd.dir/depend:
	cd "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2" "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2" "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build" "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build" "/home/christian/Desktop/MHPC_LECTURES/11_week*/Group_2/build/CMakeFiles/ljmd.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : CMakeFiles/ljmd.dir/depend

