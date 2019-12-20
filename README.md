# Orts-Del-Immagine_CurrentBio_2020
MATLAB scripts used for calcium imaging analysis, in Orts-Del'Immagine et al. 2020, Current Biology:


## AutoCalciumvsContraction_FlashPlusRepeatedStimulus_V180814.m

This code has been used to generate analysis of CSF-cNs responses upon passive tail curvature in Figure 1B and Figure 3B.
It calculates dF/F for several ROIs in one image after repeated stimulations happening at known frequency. 


## AutoCalciumvsContraction_RepeatedStimulus_V180802.m

This code has been used to generate analysis of CSF-cNs responses after active stimulation in Figure 1A and Figure 3A.
It calculates dF/F for several ROIs in one image after repeated stimulations happening at known frequency, after a user-defined flash light pulse to detect the first stimulation.


## M2avi_color_ROIs

Generates a ".avi" file within the scripts AutoCalciumvsContraction_FlashPlusRepeatedStimulus_V180814.m and AutoCalciumvsContraction_RepeatedStimulus_V180802.m
It allows to visualize the registered image.

## getROIcell.m

Allows to define regions of interest to be analyzed within the scripts AutoCalciumvsContraction_FlashPlusRepeatedStimulus_V180814.m and AutoCalciumvsContraction_RepeatedStimulus_V180802.m

## linspecer.m

Returns a color map for the analysis in the scripts AutoCalciumvsContraction_FlashPlusRepeatedStimulus_V180814.m and AutoCalciumvsContraction_RepeatedStimulus_V180802.m

## multitiff2M.m

Converts .tiff files informations to be analyzed in the scripts AutoCalciumvsContraction_FlashPlusRepeatedStimulus_V180814.m and AutoCalciumvsContraction_RepeatedStimulus_V180802.m

