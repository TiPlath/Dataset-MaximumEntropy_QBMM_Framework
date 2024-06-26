*Version: 1.0, Date: 25.06.2024*
 
### __Data to reproduce the paper: "Population balance modelling and reconstruction by quadrature method of moments for wet granulation"__

__Authors: Plath, Timo¹*; Luding, Stefan¹; Weinhart, Thomas¹__

¹ Thermal and Fluid Engineering (ET) | University of Twente | P.O. Box 217, 7500 AE Enschede, The Netherlands |  
\* Corresponding author: Plath, Timo; t.plath@utwente.nl

***General Introduction***

The software consists of a fully working MATLAB (Version R2023a) script for a non-dimensional and dimensional direct quadrature method of moments coupled to a maximum entropy reconstruction. Additional methods that are needed to run the quadrature method of moments (reading in data, adaptive Wheeler algorithm, kernels, etc.) and the maximum entropy reconstruction (Gaussian qudrature, Cholesky inversion, etc.) including input data consisting of different density distributions are available. The input data was derived from a different dataset, see [Plath et al. 2021](https://doi.org/10.4121/14248433.v1).
Furthermore, this set of software allows to reproduce all figures and verifies integrity of the title paper. It makes the software reusable and accessible according to the FAIR principles.  
This research was funded by NWO VIDI grant 16604, *“Virtual Prototyping of Particulate Processes”*.

***Abstract***

Population balance methods utilised in multiphase flow simulations mark a significant advancement in computational fluid dynamics. However, existing approaches exhibit shortcomings, such as being prone to inaccuracies or being computationally prohibitive. Addressing these challenges, a recent innovation in closure for the method of moments is the introduction of quadrature based moments methods (QBMM). Discretising a distribution by a number of discrete elements, QBMM facilitate efficient and accurate tracking of density distributions, particularly for particle size distributions (PSD). However, obtaining the full particle size distribution information using these methods requires reconstructing the distribution from a finite set of moments, which is not a trivial step.
This study introduces a novel integration of the maximum entropy reconstruction (MER) into QBMM, establishing a robust and rapid framework for the time evolution and reconstruction of PSDs. As proof of concept for this framework, we focus on the direct quadrature method of moments (DQMOM) with spatially homogeneous and monovariate distributions. We show that coupling of MER with DQMOM has numerous advantages. To verify the framework, special cases of constant growth, aggregation, and breakage are considered for which analytical solutions can be found.
Furthermore, we show the advantage of using DQMOM with volume-based over length-based distributions, and address numerical as well as theoretical issues. Validation of the framework is successfully conducted on the evolution of the PSD from a twin-screw wet granulation dataset, considering all primary physical mechanisms inherent in a wet granulation process, namely growth, aggregation, and breakage. This showcases the consistency of the proposed framework and underscores its applicability to real-world scenarios.

***Description of the data in this software set***

Wherever possible, all data is provided in established open file formats (.csv, .md, .eps, .gnu and .txt). Non-open formats include MATLAB files (.m, .fig) and a file used for the Curve Fitting Toolbox (.sfit) where curve fitting can be done in a graphical user interface. We refer to the paper for variable definition and names in generated figures from the software. Please do not move files or folders, as the scripts are fully functional in the current state of this dataset.
The MATLAB version used was *R2023a*. MATLAB scripts are commented to make them self explanatory.


***Detailed description of the data and folder structure***

Subfolders are named intuitively and the data inside subfolders should be described by the subfolders name. See folder structure of this dataset below

```
Dataset/  
├── Chapter3_InitialMER.m						: MATLAB script for the initial MER reconstruction tested on a log-normal density distribution
├── Chapter4_ComputePSDError.m						: MATLAB script which compares the length-based and volume-based volume conservation error
├── Chapter4_QBMMNonDimensionalised.m					: MATLAB script which runs aggregation, breakage and growth special cases and compares it to analytical solutions
├── Chapter5_ConsistencyCheck.m						: MATLAB script which applies the DQMOM-MER framework to a twin-screw wet granulation dataset to check consistency
├── Chapter5_DatasetAnalysis.m						: MATLAB script which performs a data analysis of the twin-screw wet granulation data
├── Chapter5_MRTd32Match.m						: MATLAB functions which is used in the data analysis to run DQMOM for different parameter pairs of SFL and L/S
├── Chapter5_ModelExtrapolation.m					: MATLAB script that extrapolates to parameter pairs outside the data range of the twin-screw wet granulation dataset
├── CholeskyInversion.m							: MATLAB function that performs a Cholesky inversion to check if a matrix is positive definite
├── ComputeMoments.m							: MATLAB function that computes moments from particle size distribution data
├── ComputePSDError.m							: MATLAB function for the initial reconstruction bisection method to compute the error for a certain c-value 
├── CutHighSizeRatio.m							: MATLAB function that cuts the PSD at a certain minimum and maximum quantile
├── DensityDistributionAdvectionEquation.m				: MATLAB script that compute the solution of a constant advection equation applied to a density distribution
├── DirectQuadratureMethodOfMoments.m					: MATLAB function that performs a volume-based DQMOM timestep
├── DirectQuadratureMethodOfMomentsLengthBased.m			: MATLAB function that performs a length-based DQMOM timestep
├── DirectQuadratureMethodOfMomentsNonDimensionalized.m			: MATLAB function that performs a non-dimensional volume-based DQMOM timestep
├── DirectQuadratureMethodOfMomentsNonDimensionalizedLengthBased.m	: MATLAB function that performs a non-dimensional length-based DQMOM timestep
├── Gauss.m								: MATLAB function that performs a Gauss-Legendre quadrature for the maxium entropy reconstruction
├── MaximumEntropyUnitTest.m						: MATLAB script that tests the functionality of the MER
├── ReadData.m								: MATLAB function that reads in a .csv file consisting of particle size distribution data
├── Wheeler.m								: MATLAB function for the adaptive Wheeler algorithm (Marchisio and Fox 2013)
├── convertCDFtoPDF.m							: MATLAB function to convert a cumulative density function (CDF) into a (probability) density function (PDF)	
├── convertPDFtoCDF.m							: MATLAB function to convert a probability density function (PDF) into a cumulative density function (CDF)
├── getMomenta.m							: MATLAB function to compute the moments from discrete quadrature distribution nodes and weights
├── getPSD.m								: MATLAB function that reconstructs a PSD by MER using a Newton-Raphson optimisation with Lagrangian multipliers
├── jensenShannonDivergence.m						: MATLAB function that computes the Jenson-Shannon divergence of two density distributions
├── plotInitLeastErrorPSD.m						: MATLAB function that finds a c-value that reconstructs the PSD with the least error to the initial distribution
├── plotMaximumEntropyReconstruction.m					: MATLAB function to plot the reconstructed continuous PSD by MER with a certain resolution
├── validateCDF.m							: MATLAB function that validates a CDF by checking if it starts with a zero
├── Data/
│   ├── InitialLactoseMCCPVP_CDF-v_L.csv				: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data
│   ├── MCCLactosePVP_CenterPointLongScrew_N22_Validation.csv		: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N22
│   ├── MCCLactosePVP_CenterPointShortScrew_N28_Validation.csv		: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N28
│   ├── MCCLactosePVP_HighSFL_HighLS_N16.csv				: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N16
│   ├── MCCLactosePVP_HighSFL_LowLS_N15.csv				: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N15
│   ├── MCCLactosePVP_LowSFL_HighLS_N14.csv				: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N14
│   ├── MCCLactosePVP_LowSFL_LowLS_N13.csv				: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of experiment N13
│   └── MCCLactosePVP_MediumSFL_MediumLS_N26.csv			: Spreadsheet file with length-based volumetric cumulative particle size distribution (Q_3, v_L) data of Experiment N26
├── Figures/
│   ├── AlgorithmGeneral.eps						: Figure showing the general algorithm flowchart of the DQMOM-MER framework
│   ├── CenterPointConsistencyCheckMoments.eps				: Figure showing the moment evolution from DQMOM simulations of experiment N22 and N28
│   ├── CenterPointConsistencyCheckReconstruction.eps			: Figure showing the MER from DQMOM simualations of experiment N22 and N28
│   ├── ExtrapolationHighParameters.eps					: Figure showing the extrapolated MER from DQMOM simlations for very high SFL and L/S values
│   ├── ExtrapolationLowParameters.eps					: Figure showing the extrapolated MER from DQMOM simlations for very low SFL and L/S values
│   ├── MEInitialReconstructionBimodal.eps				: Figure showing the initial MER for an optimal c-value on a bimodal log-normal distribution
│   ├── MEInitialReconstructionUnimodal.eps				: Figure showing the initial MER for an optimal c-value on a unimodal log-normal distribution
│   ├── PureAggregationMomentsNondimensional.eps			: Figure showing the non-dimensional moment evolution from DQMOM simulations of pure constant aggregation
│   ├── PureAggregationNondimensional.eps				: Figure showing the MER from non-dimensional DQMOM simulations of pure constant aggregation
│   ├── PureBreakageMomentsNondimensional.eps				: Figure showing the non-dimensional moment evolution from DQMOM simulations of pure constant breakage
│   ├── PureBreakageNondimensional.eps					: Figure showing the MER from non-dimensional DQMOM simulations of pure constant breakage
│   ├── PureGrowthMomentsNondimensional.eps				: Figure showing the non-dimensional moment evolution from DQMOM simulations of pure constant growth
│   ├── PureGrowthNondimensional.eps					: Figure showing the MER from non-dimensional DQMOM simulations of pure constant growth
│   ├── VolumeConservationError.eps					: Figure showing the volume conservation for volume-based and length-based DQMOM using different timestep magnitudes
│   ├── aggregationRateContourPlot.eps					: Figure showing the contour plot for the miminum aggregation rate a_s against MRT and d_32
│   ├── aggregationRateContourPlot.fig					: MATLAB figure to explore the contour plot for the miminum aggregation rate a_s against MRT and d_32
│   ├── aggregationRateMeasuredVsPredicted.eps				: Figure showing the correlation of measured against predicted minimum aggregation rate
│   ├── breakageRateFit.eps						: Figure showing the non-dimensional breakage rate fit against d_32
│   ├── d32ContourPlot.eps						: Figure showing the contour plot for the Sauter mean diameter d_32 against SFL and L/S
│   ├── d32ContourPlot.fig						: MATLAB figure to explore the contour plot for the Sauter mean diameter d_32 against SFL and L/S
│   ├── d32MeasuredVsPredicted.eps					: Figure showing the correlation of measured against predicted Sauter mean diameter
│   └── testFit.sfit							: MATLAB CurveFitter file to explore the conducted fits in a GUI of the curve fitting toolbox
├── Fits/
│   ├── aggregationRateVsD32VsMRT.gnu					: Gnuplot file to plot the contour plot for the minimum aggregation rate a_s against MRT and d_32
│   ├── aggregationRateVsD32VsMRT.txt					: Text file containing data for the contour plot of minimum aggregation a_s against MRT and d_22
│   ├── aggregationRateVsMRT.gnu					: Gnuplot file to plot the minmum aggregation rate fit against MRT
│   ├── aggregationRateVsMRT.txt					: Text file containing data for the minimum aggregation rate fit against MRT
│   ├── breakageRateVsD32.gnu						: Gnuplot file to plot the non-dimensional breakage rate fit against d_32
│   ├── breakageRateVsD32.txt						: Text file containing data for the non-dimensional breakage rate fit against d_32
│   ├── d32vsLS.gnu							: Gnuplot file to plot the Sauter mean diameter fit against L/S
│   ├── d32vsLS.txt							: Text file containing data for the Sauter mean diameter fit against L/S
│   ├── d32vsLSvsSFL.gnu						: Gnuplot file to plot the contour plot for the Sauter mean diameter d_32 against SFL and L/S
│   └── d32vsLSvsSFL.txt						: Text file containing data for the contour plot of Sauter mean diameter d_32 against SFL and L/S
└── Kernels/
    ├── GrowthAggregationBreakage.m					: MATLAB function to compute the volume-based source terms for constant growth, aggregation and breakage
    ├── GrowthAggregationBreakageLengthBased.m				: MATLAB function to compute the length-based source terms for constant growth, aggregation and breakage
    ├── GrowthAggregationBreakageNonDimensionalized.m			: MATLAB function to compute the non-dimensional volume-based source terms for constant growth, aggregation and breakage
    ├── GrowthAggregationBreakageNonDimensionalizedLengthBased.m	: MATLAB function to compute the non-dimensional length-based source terms for constant growth, aggregation and breakage
    ├── GrowthHydrodynamicAggregationPowerLawBreakage.m			: MATLAB function to compute the non-dimensional volume-based source terms for growth, hydrodynamic aggregation and power-law breakage
    └── GrowthHydrodynamicAggregationPowerLawBreakageLengthBased.m	: MATLAB function to compute the non-dimensional length-based source terms for growth, hydrodynamic aggregation and power-law breakage
```