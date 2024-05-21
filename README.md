# Paper-Geometric-Mechanics-of-Contact-Switching-Systems

This is the code repository for the paper "Geometric Mechanics of Contact-Switching Systems": <10.1109/LRA.2023.3327930>. This repository can generate the plots provided in the paper.

## Information on using this repository

Each file within the "Code" folder is intended to be run sectionwise. For more details on the terminology and techniques, please refer to the paper linked above.

1. The stratified panels shown in the paper (in Fig.1(c), Fig.3(d), Fig.3(e), Fig.4(b), and Fig.6) are generated in the "EfficiencyAndPanels.m" code.

2. The contact interpolation Fig.3(a) is generated in the "ContactInterpolationFunction.m" code.

3. The local connection vector fields and the constraint curvature functions shown in figures 3(b) and 3(c) are generated using the GeometricSystemPlotter repository developed by collaborators at OSU. A link to this repository can be found here: <https://github.com/OSU-LRAM/GeometricSystemPlotter>. The code for the class of systems explored in this paper is provided in the "GeometricSystemPlotter" folder in the exact same structure of the GeometricSystemPlotter repository. So, if you're interested in testing things out on the systems used in this paper, please clone the repository from that link and then copy the contents of the "GeometricSystemPlotter" folder into their corresponding location.

4. The shape change or gait trajectory in Fig.4(a) is generated in the "ContactInterpolationFunction.m" code.

5. All animations provided with this paper can be recreated using the "Animations_HybridContact.m" code.

6. The main video for this paper can be found in the "Video" folder. This video is intended to be a condensed illustrative aid to the paper (with less jargon).
