# FluidSimulation
This is the fluid simulation we have implemented mini-project of the UU GMT 2016-2017 Game Physics course.
The fluid is simulated using Smoothed Particle Hydrodynamics. We have implemented the technique described in the paper *“[Particle-Based Fluid Simulation for Interactive Applications](http://matthias-mueller-fischer.ch/publications/sca03.pdf)”* by Matthias Müller, David Charypar and Markus Gross.

## Compilation
All libraries have been included, so the compilation should be _out of the box_. Simply open the Visual Studio 2015 solution in VS2015 and compile for 64-bit. The software was developed on Windows.

## Usage
The camera can be moved around using the mouse. Hold the left mouse button to rotate the camera or the right mouse button to pan the camera.

Initially, the simulation starts as paused, with a square grid of particles in the center. Untick the "Paused" setting in the AntTweakBar to start the simulation.

During the simulation, the grid in which the fluid resides can be moved around. The WASD keys move the grid in the XZ plane, Space moves the grid upwards and Ctrl moves the grid downwards. While Shift is held down, the movement speed is doubled.

The program initially launches in an empty scene: no mesh objects are visible. One can change the currently active scene using the function keys:

- F1: Empty scene
- F2: Scene with a horizontal plane
- F3: Scene with a vertical plane
- F4: Scene with a horizontal plane, where fluid can flow downwards past the plane
- F5: Scene with two horizontal planes where fluid can flow past
- F6: Scene featuring the Stanford Bunny

Note that when the scene is switched, the fluid is not reset (allowing the user to change the objects in the scene on-the-fly). To reset the particles of the fluid to a cube, press the "Reset particles" button in the AntTweakBar.

Several simulation parameters can be changed on-the-fly using the AntTweakBar. One can for example change the mass m. By default, m is set to 0.3 . The particles always spawn in a fixed grid; with a mass 3.0, they spawn "too close together", so the particles will push away from each other when the simulation starts. To prevent this from happening, one could for example set the particle mass to 0.15 .

## Feature overview
- We simulate a fluid using an SPH-based technique described in the paper *“[Particle-Based Fluid Simulation for Interactive Applications](http://matthias-mueller-fischer.ch/publications/sca03.pdf)”* by Matthias Müller, David Charypar and Markus Gross.
- Parameters of the simulation can be changed on-the-fly using AntTweakBar.
- A mesh can be constructed representing the fluid (this can be changed using the "Construct mesh" setting in the AntTweakBar). This is done by using a grid of `61 * 51 * 51 = 158661` voxels. In every voxel, a density is calculated. If this density is above a threshold, the voxel is considered solid. Marching Cubes (implemented in the PolyVox library) is then used to construct a mesh. A fast enough performance (the density is computed for 158661 voxels in every frame) is achieved by using OpenMP multithreading, SIMD instructions and a lookup table for the poly6 kernel.
- The grid in which the fluid resides can be moved around.
- Collisions between particles and meshes are handled. For this, we use the "Intersection Method" described in Section 3.1 of the paper *“[Smoothed Particle Hydrodynamics in Flood Simulations](http://www.sccg.sk/~durikovic/publications/Pub09_11_files/SCCG2010_SPH_Flood.pdf)”* by Michal Chládel amd Roman Durikovic. No library was used for the collision detection code.
- One can switch between various scenes with meshes in them using the function keys. .obj files can be loaded using tinyobjloader.

## Used libraries
- The simulation is drawn using OpenGL. For this, we use GLEW and GLFW.
- For our vectors, matrices etc. GLM is used.
- AntTweakBar allows for changing parameters during the simulation.
- We perform Marching Cubes to construct a mesh representing the water. For this, PolyVox is used.
- We use tinyobjloader to load .obj files.

## Screenshots
Note that these screenshots are in a 1280 x 1024 resolution, but the software runs at a 1280 x 800 resolution by default (mainly to make it fit on laptop screens with a low resolution). The resolution of the software can be changed in main.cpp.
