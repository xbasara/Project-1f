# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.0

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files (x86)\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files (x86)\CMake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\Tristan Jones\Desktop\Cis441\Project1F"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\Tristan Jones\Desktop\Cis441\Project1F"

# Include any dependencies generated for this target.
include CMakeFiles/project1F.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/project1F.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/project1F.dir/flags.make

CMakeFiles/project1F.dir/project1F.cxx.obj: CMakeFiles/project1F.dir/flags.make
CMakeFiles/project1F.dir/project1F.cxx.obj: CMakeFiles/project1F.dir/includes_CXX.rsp
CMakeFiles/project1F.dir/project1F.cxx.obj: project1F.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report "C:\Users\Tristan Jones\Desktop\Cis441\Project1F\CMakeFiles" $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/project1F.dir/project1F.cxx.obj"
	C:\MinGW\bin\g++.exe   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles\project1F.dir\project1F.cxx.obj -c "C:\Users\Tristan Jones\Desktop\Cis441\Project1F\project1F.cxx"

CMakeFiles/project1F.dir/project1F.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project1F.dir/project1F.cxx.i"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_FLAGS) -E "C:\Users\Tristan Jones\Desktop\Cis441\Project1F\project1F.cxx" > CMakeFiles\project1F.dir\project1F.cxx.i

CMakeFiles/project1F.dir/project1F.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project1F.dir/project1F.cxx.s"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_FLAGS) -S "C:\Users\Tristan Jones\Desktop\Cis441\Project1F\project1F.cxx" -o CMakeFiles\project1F.dir\project1F.cxx.s

CMakeFiles/project1F.dir/project1F.cxx.obj.requires:
.PHONY : CMakeFiles/project1F.dir/project1F.cxx.obj.requires

CMakeFiles/project1F.dir/project1F.cxx.obj.provides: CMakeFiles/project1F.dir/project1F.cxx.obj.requires
	$(MAKE) -f CMakeFiles\project1F.dir\build.make CMakeFiles/project1F.dir/project1F.cxx.obj.provides.build
.PHONY : CMakeFiles/project1F.dir/project1F.cxx.obj.provides

CMakeFiles/project1F.dir/project1F.cxx.obj.provides.build: CMakeFiles/project1F.dir/project1F.cxx.obj

# Object files for target project1F
project1F_OBJECTS = \
"CMakeFiles/project1F.dir/project1F.cxx.obj"

# External object files for target project1F
project1F_EXTERNAL_OBJECTS =

project1F.exe: CMakeFiles/project1F.dir/project1F.cxx.obj
project1F.exe: CMakeFiles/project1F.dir/build.make
project1F.exe: C:/VTKbin/lib/libvtkalglib-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkChartsCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonColor-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonDataModel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonMath-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtksys-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonMisc-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonSystem-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonTransforms-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInfovisCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersExtraction-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonExecutionModel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersGeneral-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonComputationalGeometry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersStatistics-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingFourier-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingContext2D-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersGeometry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersSources-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingFreeType-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkfreetype-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkzlib-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkftgl-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingOpenGL-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingHybrid-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOImage-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkDICOMParser-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkmetaio-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkjpeg-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkpng-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtktiff-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkDomainsChemistry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOXML-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOGeometry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkjsoncpp-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOXMLParser-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkexpat-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkexoIIc-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkNetCDF-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkNetCDF_cxx-6.1.dll.a
project1F.exe: C:/VTKbin/lib/vtkhdf5_hl-6.1.lib
project1F.exe: C:/VTKbin/lib/vtkhdf5-6.1.lib
project1F.exe: C:/VTKbin/lib/libvtkFiltersAMR-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkParallelCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOLegacy-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersFlowPaths-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersGeneric-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersHybrid-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingSources-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersHyperTree-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersImaging-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingGeneral-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersModeling-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersParallel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersParallelImaging-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersProgrammable-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersSelection-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersSMP-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersTexture-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersVerdict-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkverdict-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkGeovisCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInfovisLayout-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInteractionStyle-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInteractionWidgets-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingAnnotation-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingColor-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingVolume-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkViewsCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkproj4-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkgl2ps-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingMath-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingMorphological-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingStatistics-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingStencil-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInteractionImage-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOAMR-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOEnSight-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOExodus-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOExport-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingGL2PS-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingLabel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOImport-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOInfovis-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtklibxml2-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOLSDyna-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOMINC-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOMovie-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkoggtheora-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIONetCDF-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOParallel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOPLY-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOSQL-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtksqlite-6.1.a
project1F.exe: C:/VTKbin/lib/libvtkIOVideo-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingFreeTypeOpenGL-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingImage-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingLIC-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingLOD-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingVolumeAMR-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingVolumeOpenGL-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkViewsContext2D-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkViewsGeovis-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkViewsInfovis-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkgl2ps-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkexoIIc-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersParallel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIONetCDF-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkNetCDF_cxx-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkNetCDF-6.1.dll.a
project1F.exe: C:/VTKbin/lib/vtkhdf5_hl-6.1.lib
project1F.exe: C:/VTKbin/lib/vtkhdf5-6.1.lib
project1F.exe: C:/VTKbin/lib/libvtkFiltersAMR-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkParallelCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOLegacy-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkGeovisCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOXML-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOGeometry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkjsoncpp-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOXMLParser-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkexpat-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkproj4-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkChartsCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonColor-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingContext2D-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingOpenGL-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersImaging-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInfovisLayout-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInfovisCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkViewsCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInteractionWidgets-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingHybrid-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOImage-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkDICOMParser-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkIOCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkmetaio-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkpng-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtktiff-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkjpeg-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersHybrid-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingGeneral-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingSources-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersModeling-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkInteractionStyle-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingAnnotation-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingColor-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingVolume-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingLabel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingFreeType-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkRenderingCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersExtraction-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersStatistics-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkalglib-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingFourier-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkImagingCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersGeometry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersSources-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersGeneral-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkFiltersCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonExecutionModel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonComputationalGeometry-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonDataModel-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonMisc-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonTransforms-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonMath-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonSystem-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkCommonCore-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtksys-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkftgl-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkfreetype-6.1.dll.a
project1F.exe: C:/VTKbin/lib/libvtkzlib-6.1.dll.a
project1F.exe: CMakeFiles/project1F.dir/objects1.rsp
project1F.exe: CMakeFiles/project1F.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable project1F.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\project1F.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/project1F.dir/build: project1F.exe
.PHONY : CMakeFiles/project1F.dir/build

CMakeFiles/project1F.dir/requires: CMakeFiles/project1F.dir/project1F.cxx.obj.requires
.PHONY : CMakeFiles/project1F.dir/requires

CMakeFiles/project1F.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\project1F.dir\cmake_clean.cmake
.PHONY : CMakeFiles/project1F.dir/clean

CMakeFiles/project1F.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\Tristan Jones\Desktop\Cis441\Project1F" "C:\Users\Tristan Jones\Desktop\Cis441\Project1F" "C:\Users\Tristan Jones\Desktop\Cis441\Project1F" "C:\Users\Tristan Jones\Desktop\Cis441\Project1F" "C:\Users\Tristan Jones\Desktop\Cis441\Project1F\CMakeFiles\project1F.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/project1F.dir/depend

