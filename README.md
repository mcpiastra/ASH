ASH: an Automatic pipeline to generate realistic and individualized chronic Stroke volume conduction Head models.
=========================

<p align="center">
  <img src="https://github.com/mcpiastra/ASH/blob/master/ash_repo.png" width="500">
</p>

Description:
-------
ASH is a Matlab-pipeline that automatically generates realistic and individualized volume conduction head models of chronic stroke patients, by combining the already existing software SimNIBS, for the mesh generation, and LINDA, for the lesion identification.

Table of Contents
-------
- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [License](#license)

Installation:
-------
ASH is based on several external software packages, namely, Matlab, SimNIBS (including SPM12 and CAT12), LINDA (including R) and FieldTrip. 
The four software packages should be then installed in your computer. We refer to the documentation given by the software packages: [Matlab](http://matlab.com), [SimNIBS](https://simnibs.github.io/simnibs/build/html/installation/simnibs_installer.html), [LINDA](http://dorianps.github.io/LINDA/), [FieldTrip](http://www.fieldtriptoolbox.org/download/).

Usage: 
-------
In the Matlab script Pipeline_SimNIBS_LINDA.m the main steps necessary to generate the individual and realistic volume conduction head model of a chronic stroke patients (with LINDA and SimNIBS) are implemented:
### Segmentation and meshing of the whole head: 
The T1w MRI is processed by the SimNIBS software  which gives as output a tetrahedral volumetric mesh with 6 homogeneous and isotropic compartments: scalp, skull, eyes, cerebrospinal fluid (CSF), gray matter and white matter. In particular, we utilize the function headreco with the option cat that leads to the use of SPM12. 
### Segmentation of the lesion: 
Since the segmentation and meshing of the lesion compartment are not performed by SimNIBS, we use the LINDA software. LINDA (Lesion Identification with Neighborhood Data Analysis) is a neuroimaging toolkit for automatic segmentation of chronic stroke lesions based on machine learning techniques. It requires a T1w MRI as input and gives as output a volumetric mask of the lesion.
### Generation of the final mesh: 
The volumetric mesh generated in the previous step is modified to incorporate the lesion compartment. First, the coordinates of the lesion mask voxels are identified. Second, the mesh elements whose centroids are within the lesion mask are labeled as “lesion”. Finally, we make sure that the resulting lesion compartment does not contain elements of the scalp, skull or eye compartments.

Furthermore, SimNIBS simulations are performed, once the head model is generated. Two transcranial Direct Current Simulation (tDCS) electrodes were modelled as 3 mm thick, 1 cm x 1 cm elliptical patches. The electrode pairs at C3-Fp2 and the one at C4-Fp1 were selected for the ipsi- and contra-lesional stimulation of the primary motor cortex, respectively.
  

Credits: 
-------
Contributors are:
- Maria Carla Piastra, *Department of Cognitive Neuroscience, Donders Institute for Brain, Cognition and Behaviour, Radboud University Nijmegen Medical Center, Nijmegen, The Netherlands*
- Joris van der Cruijsen, *Department of Rehabilitation, Erasmus MC- University Medical Center Rotterdam, Rotterdam, the Netherlands*
- Floor EM Jeukens, *Department of Biomechanical Engineering, Delft University of Technology, Delft, the Netherlands*
- Mana Manoochehri, *Department of Biomechanical Engineering, Delft University of Technology, Delft, the Netherlands*
- Alfred C Schouten, *Department of Biomechanical Engineering, Delft University of Technology, Delft, the Netherlands and Department of Biomechanical Engineering, University of Twente, Enschede, The Netherlands*
- Ruud W Selles, *Department of Rehabilitation, Erasmus MC- University Medical Center Rotterdam, Rotterdam, the Netherlands and Department of Plastic and Reconstructive Surgery, Erasmus MC- University Medical Center Rotterdam, Rotterdam, the Netherlands*
- Thom Oostendorp, *Department of Cognitive Neuroscience, Donders Institute for Brain, Cognition and Behaviour, Radboud University Nijmegen Medical Center, Nijmegen, The Netherlands*

License: 
-------
The ASH pipeline is copyrighted free software. You can use, modify and/or redistribute it under the terms of GNU General Public License v3.0.

