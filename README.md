# **Numerical 2D-PES of Reactions with PTSB**

* Using several programs/scripts to general a numerical two-dimensional potential energy surface (2D-PES) of reaction with post-transition-state bifurcation (PTSB).
* Supporting information of [Construction of Two-Dimensional Potential Energy Surfaces of Reactions with Post-Transition-State Bifurcations.][1]
* **Author** : Hsiao-Han (Grace) Chuang - *Initial work* - 2018 May.
  
[1]: https://pubs.acs.org/doi/10.1021/acs.jctc.0c00172

## Table of Content
- [**Numerical 2D-PES of Reactions with PTSB**](#numerical-2d-pes-of-reactions-with-ptsb)
  - [Table of Content](#table-of-content)
  - [Processes to generate 2D-PESs](#processes-to-generate-2d-pess)
    - [Step 1: Calculate stationary points and IRCs.](#step-1-calculate-stationary-points-and-ircs)
    - [Step 2: Asymmetric cases: generate artificial reaction coordinate](#step-2-asymmetric-cases-generate-artificial-reaction-coordinate)
    - [Step 3: Construct x- and y-axes (optional)](#step-3-construct-x--and-y-axes-optional)
      - [*Symmetric cases*](#symmetric-cases)
      - [*Asymmetric cases*](#asymmetric-cases)
    - [Step 4: Select 1D grid points for all IRC paths](#step-4-select-1d-grid-points-for-all-irc-paths)
    - [Step 5: Scan a 2D-PES](#step-5-scan-a-2d-pes)
    - [Step 6: Plot 2D/3D figures](#step-6-plot-2d3d-figures)

## Processes to generate 2D-PESs


### Step 1: Calculate stationary points and IRCs.
- Code: 
    - `checkGau`, `runIRC.sh` 
    <!-- - FIXME: check runIRC.sh and checkGau -->
1. Quickly run the standard procedure; guess the TSS between two minimums. 
   -  Use `checkGau` to check all important info. from gaussian output files. 
        <div style='float: center'>
            <img style='width: 500px' src="./aux/Fig/checkGau.png"></img>
        </div> 
2. Re-run the whole process again, because 
   1. double check the accuracy of assigned mechanism
   2. double check the atomic index and order for all systems are consistent
      - check the order of atoms between TS1 and TS2 systems. Reorder them via gaussview if they are different. 
3. Gain enough grid points of an IRC path
   1. use the homemade script `runIRC.sh` to test different combinations 
   2. the direction of IRC path in the output files may be wrong, thus the name of output file may need to be modified manually. 
4. Double check the orientation for all systems.
    - Use *gaussview* or *jmol* to check the molecular orientation.
    - Show axes is easier to check.
    
### Step 2: Asymmetric cases: generate artificial reaction coordinate
*For asymmetric cases, build the artificial reaction coordinate in the TSS1 forward direction.* 

- Code: 
    - `runArticIRC.sh`, `genArticStruc.f90` 
- Input:
  - $(selectCoord).log, TS2.log
- Output: 
  - Artic1D.xyz, Artic1D_PEC.dat
  
1. Compare the energy difference between TSS1 and TSS2:
    1. TSS1 > TSS2
        - go to the following steps.
    2. TSS1 < TSS2 and the energy difference is small 
       - beyond this project, talk to Grace : ) 
2. Select a point from TSS1 IRC forward direction.
   - near the shoulder of its energy profile.
   - near the gradient which is close to zero, or has a turn over point.
3. Execute `runArticIRC.sh` to generate a serious of structures which use program `genArticStruc.f90` to build artificial reaction coordinate. 
   - After it fulfills the criteria (i.e. energy difference between the last point and TS2 is less than 0.0001 hartree; less than 1 kcal/mol), the program will stop automatically. 
   - Use `jmol` to check the structures via file *Artic1D.xyz*.
   - Use `gnuplot` to plot energy profile via file *Artic1D_PEC.dat*. 

### Step 3: Construct x- and y-axes (optional) 
*The following steps is part of the detail in the macro script, `run1Dgrid.sh`. Step 3 can be skipped if using above script, or, if something wrong that this step can be a reference to debug.*

- Code: 
    - `getIRCcurve`, `getIRCstruc`, `rev1Dstruc`,`checkGau`,`getCoord`

#### *Symmetric cases*
- Output:
  - x.xyz, y_F.xyz, y_R.xyz 
1. Copy 4 IRC output files from /OPT to /1D as the input files for the following steps.
2. Generate x.xyz 
   1. Use `getIRCcurve` to extract energy profiles, and then use *gnuplot* to plot the potential energy curves.   
       - The sign of coordinate may need to be modified; use *awk*.
            <div style='float: center'>
            <img style='width: 400px' src="./aux/Fig/getIRCcurve.png"></img>
            </div> 
    1. Use `getIRCstruc` to extract all the structures.
       - In the IRC file, it starts from TS structure by default. So use `rev1Dstruc` to reverse the direction for reverse direction. 
           <div style='float: center'>
           <img style='width: 400px' src="./aux/Fig/getIRCstruc.png"></img>
           </div> 
    2. Combine the forward and reverse direction to form a whole path as x.xyz.
          - x.xyz =  R $\to$  TSS1 $\to$  TSS2.
3. Generate y_F.xyz and y_R.xyz.
   - y_F.xyz = TSS2 $\to$ P1.
   - y_R.xyz = TSS2 $\to$ P2.
4.  Check the orientation via *jmol* and then open the axes function; sometimes the TS structure and IRC structures has shift, and cannot properly overlap to each other. 

#### *Asymmetric cases* 
- Input:
  - TS1_F.xyz, Artic1D.xyz
- Output:  
  - TS1_F_Artic1D.xyz, x.xyz, y_F.xyz and y_R.xyz
1. Generate x.xyz
   1. For TS1 forward direction, combine the TS1_F.xyz from reactant to the selected point and Artic1D.xyz. After that, name the new file as TS1_F_Artic1D.xyz .
   2. Use *jmol* to double check the orientation of TS1_F_Artic1D.xyz .     
   3. Combine TS1_F_Artic1D.xyz and TS1_R.xyz as x.xyz
      - x.xyz = R $\to$ TSS1 $\to$ selected point $\to$ TSS2
2. Generate y_F.xyz and y_R.xyz: same process as the symmetric case.

### Step 4: Select 1D grid points for all IRC paths

- Code: `run1Dgrid.sh`, `get1Dgrid.sh`
  - `run1Dgrid.sh`: `getIRCcurve`, `getIRCstruc`, `selectIRCstruc.sh`, `writeGauInpV`
  - `get1Dgrid.sh`: `checkGau`, `rev1Dstruc`, `rot.py`
- Input:
  - Symmetric cases (5 files): 
    1. header.dat (route section)
    2. TS1_F/R.log 
    3. TS2_F/R.log. 
  - Asymmetric cases (7 files):
    1. header.dat
    2. TS1_F/R.log
    3. TS2_F/R.log 
    4. \$(coordinate of select point).log 
    5. Artic1D.xyz 
- Output:  
  - Symmetric cases
    - $(NGrid)_TS1_F.xyz
  - Asymmetric cases
    - $(NGrid)_TS1_F_Artic1D.xyz
  - $(NGrid)_TS1_R.xyz
  - $(NGrid)_TS2_F.xyz
  - $(NGrid)_TS2_R.xyz
  - E_scan.*.dat 

1. Execute script `run1Dgrid.sh` to generate selected 1D grid points.
   - If the file, Artic1D.xyz, exists, then it is an asymmetric case. Or, it is a symmetric case. Do not remove it durning the calculation. 
   - The default amount of selected grid point is 30. If one wants to change it, go to this part of `run1Dgrid.sh`, and then uncomment *line 211* and comment *line 212*. 
        <div style='float: center'>
            <img style='width: 300px' src="./aux/Fig/gridP.png"></img>
        </div> 
    - This script uses `selectIRCstruc.sh` to select points along the new reaction coordinate in order to reduce the computational cost. 
        <div style='float: center'>
            <img style='width: 400px' src="./aux/Fig/selectIRCstruc.png"></img>
        <div>
    - Use gnuplot to double check the geometries and the selected geometries. 
2. Execute script `get1Dgrid` to collect data.
   - Plot the potential energy curves of selected grid points to see if the selected points are reasonable.
   - **Double check the geometries! Really important!** Use gnuplot to check the geometries of x.xyz; it should start from R to TSS2. Recall that, 
     - x.xyz = R $\to$ TSS1 $\to$ selected point $\to$ TSS2.

### Step 5: Scan a 2D-PES
- Code: 
    - Generate numerical 2D-PES
      - `getIRCvec`, `Vsca1D_IRCvec.f90`, `writeGauInpV`, `qsubGau`
    - Modify the topology 
      - `getGrad.sh`, `genGradStruc.f90`, `genNewPES.sh`
    - Collect rawdata 
      - `getPESwStruc.sh`, `checkGau`
    - Plot contour and 3D version figure 
      - `plot2DPES.py`
- Input:
  - header.dat, x.xyz, y_F.xyz and y_R.xyz 
- Output: 
  - $(Pts)_E.dat, $(Pts)_Struc.xyz
  
1. Generate numerical 2D-PES 
   1. Create sub-directories and named as *Forward* and *Reverse*, which is the directions along y-axis. 
   2. In /1D directory, rename $(NGrid)_TS2_F.xyz as y_F.xyz, and $(NGrid)_TS2_R.xyz as y_R.xyz. After that, copy x.xyz, y_F.xyz and y_R.xyz from /1D to /2DPES. 
   3. Execute `getIRCvec` to calculate translational vectors from y_F.xyz and y_R.xyz.
       <div style='float: center'>
           <img style='width: 500px' src="./aux/Fig/getIRCvec.png"></img>
       </div> 
   5. Generate a serious of structures via script `Vscan1D_IRCvec`, which has an output file, *pes.txt*. 
      <div style='float: center'>
           <img style='width: 400px' src="./aux/Fig/Vscan1D_IRCvec.png"></img>
       </div> 
   6. Generate gaussian input files via `writeGauInpV`.
       <div style='float: center'>
           <img style='width: 500px' src="./aux/Fig/writeGauInpV.png"></img>
       </div> 
   7. Submit all the gaussian input files via `qsubGau`. 
       <div style='float: center'>
           <img style='width: 500px' src="./aux/Fig/qsubGau.png"></img>
       </div> 
2. Collect rawdata 
   1. Make sure all the single points are successfully calculated. 
       - Count the amount of successful output files as the primarily checking. If it is different from the expect amount of jobs, go to the next step. 
           ```
           > cat *.log | grep -c Normal
           ```
       - Execute `checkGau` to extract the list of fail jobs (i.e. red box) and the name is recorded in *Fail.txt*. 
           <div style='float: center'>
               <img style='width: 500px' src="./aux/Fig/listFailJobs.png"></img>
           </div> 

         - Usually, recalculate the fail jobs can solve the problem.
            ````
            > for name in `cat Fail.txt`
            > do
            >   g16sub $name.com $name.log 
            > done
            ````
        - If above strategy cannot solve the problem, add other SCF keywords in the fail gaussian input files, eg. SCF=xqc.
   2. Use `getPESwStruc.sh` to extract electronic energy and structures.
        <div style='float: center'>
           <img style='width: 600px' src="./aux/Fig/getPESwStruc.png"></img>
        </div> 
3. Plot contour and 3D version figure; execute step 1, 2 and 4 in `plotPESandTraj.py`, which is *line 215* , *line 218* and *line 224*. 
    <div style='float: center'>
           <img style='width: 300px' src="./aux/Fig/plotPES.png"></img>
    </div> 
4. Modify this potential if it is not 'reasonable'.
    - Code:
       - `getGrad.sh`, `genNewPES.sh`, `GenGradStruc.f90`
    - Add *force* in the header.dat file, and then recalculate all the grid points.
    - Use `getGrad.sh` to extract gradient and named as $name.grad.
    - Use `genNewPES.sh` to execute `GenGradStruc.f90`, and then generate several new potentials which are named as *G_\*.xyz*.
      - Modify the range of *ds* in `GenGradStruc.f90`.
    - Follow the original processes to calculate those potentials. 

### Step 6: Plot 2D/3D figures  
5. Plot 2D and 3D figures with the projection of trajectories; execute step 1, 3 and 5 in `plotPESandTraj.py`, which is *line 215* , *line 221* and *line 227*. 
    <div style='float: center'>
        <img style='width: 300px' src="./aux/Fig/plotPESwTraj.png"></img>
    </div>  
