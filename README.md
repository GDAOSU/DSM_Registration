# DSM_registration
Large-scale DSM registration using two-stage ICP 

![Picture1](https://github.com/Ggs1mida/DSM_registration/assets/32317924/66703956-d501-4873-b82c-c7c5ae0ba753)
![Clipboard Image (6-15-2023, 3 56 52 PM)](https://github.com/Ggs1mida/DSM_registration/assets/32317924/dedefeab-5ac2-418f-937c-104ecdbaf1ee)

# Requirement & Install
It is tested on windows 10, visual studio 2019, c++ 17 (required for argparse.hpp)  
Library: eigen, GDAL

# Usage
reg.exe -src source.tif -dst reference.tif  
Support rigid transformation (6DoF) and translation (3DoF).

![usage](https://github.com/Ggs1mida/DSM_registration/assets/32317924/07c12b67-e692-4e91-bc76-893eeca6c1ba)


# Requirement & Build

Tested on Windows 11 with Visual Studio 2019, c++17, Eigen3 3.4.0, GDAL 3.8.4, and nanoflann 1.5.1.

1. Open the project folder in CMake GUI.
2. Select Visual Studio 16 2019, x64, and the vcpkg toolchain file.
3. Click Configure and then Generate.
4. Open the generated `.sln` file in Visual Studio and build it in Release mode.

The main executable for DSM registration will be generated as `reg.exe`.
