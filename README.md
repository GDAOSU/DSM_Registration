# DSM_registration
Large-scale DSM registration using two-stage ICP 

![Picture1](https://github.com/Ggs1mida/DSM_registration/assets/32317924/66703956-d501-4873-b82c-c7c5ae0ba753)
![Clipboard Image (6-15-2023, 3 56 52 PM)](https://github.com/Ggs1mida/DSM_registration/assets/32317924/dedefeab-5ac2-418f-937c-104ecdbaf1ee)

# Usage
reg.exe -src source.tif -dst reference.tif  
Support rigid transformation (6DoF) and translation (3DoF).

![usage](https://github.com/Ggs1mida/DSM_registration/assets/32317924/07c12b67-e692-4e91-bc76-893eeca6c1ba)


# Requirement & Build

This project was tested on Windows 11 with Visual Studio 2019, c++17, and vcpkg.

vcpkg environment:

```text
vcpkg tool version: 2026-04-08-e0612b42ce44e55a0e630f2ee9d3c533a63d8bc1
vcpkg commit: 2b65c20fc66eda893aa15a15a453c3cf09500b19
vcpkg triplet: x64-windows
````

dependency versions:
```text
Eigen3: 5.0.1
GDAL: 3.12.4#1
nanoflann: 1.9.0
```

Build steps:
1. Install Visual Studio 2019 with c++ development tools.
2. Install and bootstrap vcpkg.
3. Install the required packages with vcpkg:
```bat
vcpkg install eigen3:x64-windows gdal:x64-windows nanoflann:x64-windows
```
4. Open the project folder in CMake GUI.
5. Select the generator:
```text
Visual Studio 16 2019
```
6. Select the platform:

```text
x64
```
7. Set the vcpkg toolchain file in CMake:
```text
CMAKE_TOOLCHAIN_FILE=path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
```
8. Click `Configure`, then click `Generate`.
9. Open the generated `.sln` file in Visual Studio.
10. Build the project in `Release` mode.
The main executable for DSM registration will be generated as:

```text
reg.exe
```
