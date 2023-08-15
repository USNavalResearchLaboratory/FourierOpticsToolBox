# fourierOptics ToolBox 'fOptics' version 1.0


Object-Oriented Matlab approach to do fourier optics based simulation of light propagation. Used to simulate propagation of light and images for applications in Astronomy, Microscropy, Imaging, etc. 

An optical system is created using this codebase by creating two files:
1. Parameters file
2. main_script

The Parameters file contains the optical components and simulation parameters (such as # of pixels) in one place. This gives an overview of all the components and parameters in a given simulation, with a central point of choosing simulation parameters. 
The main_script imports the parameters and calls the fOptics classes and methods to do the simulation. The fOptics package makes it easy to define common optics and perform propagation, where the user does not need to worry about issues with array sizes, sampling, and other small issues (ideally). 

### Motivation

There are many texts and sample codes out there for fourier optics propagation. However, many suffer from problems in readability, organization, and good coding practices for future extendability. This open source software builds a clear framework using object-oriented programming practices and modern Matlab features for ease of use and extendability. Matlab is used because of its popularity, support, and ease of use. However, the software here can also be used within Python using the matlab-engine python package. 
My main focus for this packaage is wavefront sensors and atmospheric propagation, but additional features can be added over time by contributors. 


## Getting started

Check out \Examples\ptsrc_aberrated.m for an example script 

For your own project, create a new parameterSetup file and start defining your optics and propagation parameters. 
Then create a new script (ex. ptsrc_aberrated) to perform your simulation. Feel free to modify/add methods and classes to existing code. 

The fOptics tool box is organized as a Matlab Project. You can integrate fourierOpticsToolBox.prj into your own matlab project. Assuming you have a matlab Project, click on your project in the matlab tab, and click +References. Add FourieropticsToolBox.prj 
This way, keep your project specific files separate from this Toolbox and stay up to date with changes to this toolbox. 

## File Organization:
(fOptics/optics) - Includes classes and functions for optical components such as Lenses and Sensors. 
(fOptics/environment) - Atmospheric turbulence and propagation mediums 
(fOptics/lightSource) - Sources of light, such as lasers and LEDs 
(fOptics/Examples) - Example scripts 
(fOptics/field) - Mainly the Efield class and functions for propagate and handle the Efield. 
(fOptics/utils) - Uncategorized utilities 
(fOptics/analysis) - Code for analyzing and visualizing 

## How to Contribute:
If this code is on Github, create issues and pull requests to add to the existing repository for others to use. 
Email for additional questions: harshil.dave@nrl.navy.mil

## Roadmap

As of Aug 09 2023
- Write unit tests to verify propagation, optics, and test optical systems. 
- Fix issues with optimization and sampling criteria 
- Create scripts for Python integration
- Add class for Spatial Light Modulator 
- Add additional propagation methods to include additional atmospheric methods

## Authors and acknowledgment
Created by Harshil Dave - Naval Research Lab 
The authors acknowledge contributions from past and present NRL colleagues, open-source code, and academic publications that have made code available. Sources and licenses are included within specific functions and scripts where appropriate. 

## Distribution Statement A. Approved for public release: distribution is unlimited

## License

The source code is in the public domain and not licensed or under copyright. The information and software may be used freely by the public. As required by 17 U.S.C. 403, third parties producing copyrighted works consisting predominantly of the material produced by U.S. government agencies must provide notice with such work(s) identifying the U.S. Government material incorporated and stating that such material is not subject to copyright protection.

Derived works shall not identify themselves in a manner that implies an endorsement by or an affiliation with the Naval Research Laboratory.

RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.



