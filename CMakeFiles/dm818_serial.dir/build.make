# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/latho12/git/DM818-Mandatory-Assignment-2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/latho12/git/DM818-Mandatory-Assignment-2

# Include any dependencies generated for this target.
include CMakeFiles/dm818_serial.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dm818_serial.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dm818_serial.dir/flags.make

CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o: CMakeFiles/dm818_serial.dir/flags.make
CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o: serial-other/source/serial.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/latho12/git/DM818-Mandatory-Assignment-2/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o -c /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial.cpp

CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial.cpp > CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.i

CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial.cpp -o CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.s

CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.requires:
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.requires

CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.provides: CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.requires
	$(MAKE) -f CMakeFiles/dm818_serial.dir/build.make CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.provides.build
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.provides

CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.provides.build: CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o

CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o: CMakeFiles/dm818_serial.dir/flags.make
CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o: serial-other/source/serial_main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/latho12/git/DM818-Mandatory-Assignment-2/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o -c /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial_main.cpp

CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial_main.cpp > CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.i

CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial_main.cpp -o CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.s

CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.requires:
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.requires

CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.provides: CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.requires
	$(MAKE) -f CMakeFiles/dm818_serial.dir/build.make CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.provides.build
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.provides

CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.provides.build: CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o

CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o: CMakeFiles/dm818_serial.dir/flags.make
CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o: serial-other/source/common.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/latho12/git/DM818-Mandatory-Assignment-2/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o -c /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/common.cpp

CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/common.cpp > CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.i

CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/common.cpp -o CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.s

CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.requires:
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.requires

CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.provides: CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.requires
	$(MAKE) -f CMakeFiles/dm818_serial.dir/build.make CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.provides.build
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.provides

CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.provides.build: CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o

CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o: CMakeFiles/dm818_serial.dir/flags.make
CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o: serial-other/source/grid.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/latho12/git/DM818-Mandatory-Assignment-2/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o -c /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/grid.cpp

CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/grid.cpp > CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.i

CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/grid.cpp -o CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.s

CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.requires:
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.requires

CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.provides: CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.requires
	$(MAKE) -f CMakeFiles/dm818_serial.dir/build.make CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.provides.build
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.provides

CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.provides.build: CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o

CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o: CMakeFiles/dm818_serial.dir/flags.make
CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o: serial-other/source/serial_base.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/latho12/git/DM818-Mandatory-Assignment-2/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o -c /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial_base.cpp

CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial_base.cpp > CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.i

CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/latho12/git/DM818-Mandatory-Assignment-2/serial-other/source/serial_base.cpp -o CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.s

CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.requires:
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.requires

CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.provides: CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.requires
	$(MAKE) -f CMakeFiles/dm818_serial.dir/build.make CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.provides.build
.PHONY : CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.provides

CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.provides.build: CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o

# Object files for target dm818_serial
dm818_serial_OBJECTS = \
"CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o" \
"CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o" \
"CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o" \
"CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o" \
"CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o"

# External object files for target dm818_serial
dm818_serial_EXTERNAL_OBJECTS =

dm818_serial: CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o
dm818_serial: CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o
dm818_serial: CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o
dm818_serial: CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o
dm818_serial: CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o
dm818_serial: CMakeFiles/dm818_serial.dir/build.make
dm818_serial: CMakeFiles/dm818_serial.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable dm818_serial"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dm818_serial.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dm818_serial.dir/build: dm818_serial
.PHONY : CMakeFiles/dm818_serial.dir/build

CMakeFiles/dm818_serial.dir/requires: CMakeFiles/dm818_serial.dir/serial-other/source/serial.cpp.o.requires
CMakeFiles/dm818_serial.dir/requires: CMakeFiles/dm818_serial.dir/serial-other/source/serial_main.cpp.o.requires
CMakeFiles/dm818_serial.dir/requires: CMakeFiles/dm818_serial.dir/serial-other/source/common.cpp.o.requires
CMakeFiles/dm818_serial.dir/requires: CMakeFiles/dm818_serial.dir/serial-other/source/grid.cpp.o.requires
CMakeFiles/dm818_serial.dir/requires: CMakeFiles/dm818_serial.dir/serial-other/source/serial_base.cpp.o.requires
.PHONY : CMakeFiles/dm818_serial.dir/requires

CMakeFiles/dm818_serial.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dm818_serial.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dm818_serial.dir/clean

CMakeFiles/dm818_serial.dir/depend:
	cd /home/latho12/git/DM818-Mandatory-Assignment-2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/latho12/git/DM818-Mandatory-Assignment-2 /home/latho12/git/DM818-Mandatory-Assignment-2 /home/latho12/git/DM818-Mandatory-Assignment-2 /home/latho12/git/DM818-Mandatory-Assignment-2 /home/latho12/git/DM818-Mandatory-Assignment-2/CMakeFiles/dm818_serial.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dm818_serial.dir/depend
