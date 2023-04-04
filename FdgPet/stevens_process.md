```matlab
%-----------------------------------------------------------------------
% Job saved on 04-Dec-2022 21:45:41 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
subjects=[001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041 042 043 044 045 046 047 048 049 050 051 052 053 054 055 056 057 058 59 060 061 062 063 064 065 066 067 068 069 070 071 072 073 074 075];
```

```matlab
for subject=subjects

subject1=num2str(subject, '%03d');

matlabbatch{1}.spm.spatial.realign.estwrite.data = {
                            {['.../OneDrive-VUMC/HN_images/Carboplatin_Taxol/patient_' subject1 '/crop/pre_crop.nii,1']}
                            {['.../OneDrive-VUMC/HN_images/Carboplatin_Taxol/patient_' subject1 '/crop/post_crop.nii,1']}
                            }';
```
Loop over `subjects` array, convert Numeric to String, concatenate to make path.

Realign a time-series of images acquired from the same subject using a least squares approach and a 6 parameter (rigid body) spatial transformation [43]. The first image in the list specified by the user is used as a reference to which all subsequent scans are realigned. The reference scan does not have to the the first chronologically and it may be wise to chose a "representative scan" in this role. The aim is primarily to remove movement artefact in fMRI and PET time-series (or more generally longitudinal studies). The headers are modified for each of the input images, such that. they reflect the relative orientations of the data. The details of the transformation are displayed in the results window as plots of translation and rotation. A set of realignment parameters are saved for each session, named rp_*.txt. These can be modelled as confounds within the general linear model [43].

**Data**

Add new sessions for this subject. In the coregistration step, the sessions are first realigned to each other, by aligning the first scan from each session to the first scan of the first session. Then the images within each session are aligned to the first image of the session. The parameter estimation is performed this way because it is assumed (rightly or not) that there may be systematic differences in the images between sessions.


```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
```
Set the quality of realignment. This parameter is involved in selecting number of voxels used in estimation.
```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
```
The separation (in mm) between the points sampled in the reference image. Smaller sampling distances gives more accurate results, but will be slower.

```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 8;
```
The FWHM of the Gaussian smoothing kernel (mm) applied to the images before estimating the realignment parameters. * PET images typically use a 7 mm kernel. * MRI images typically use a 5 mm kernel.

```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
```
Register to first: Images are registered to the first image in the series. Register to mean: A two pass procedure is used in order to register the images to the mean of the images after the first realignment. PET images are typically registered to the mean. This is because PET data are more noisy than fMRI and there are fewer of them, so time is less of an issue. MRI images are typically registered to the first image. The more accurate way would be to use a two pass procedure, but this probably wouldn’t improve the results so much and would take twice as long to run.
```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 1;
```
Register to first: Images are registered to the first image in the series. Register to mean: A two pass procedure is used in order to register the images to the mean of the images after the first realignment. PET images are typically registered to the mean. This is because PET data are more noisy than fMRI and there are fewer of them, so time is less of an issue. MRI images are typically registered to the first image. The more accurate way would be to use a two pass procedure, but this probably wouldn’t improve the results so much and would take twice as long to run.
```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
```
Directions in the volumes the values should wrap around in. For example, in MRI scans, the images wrap around in the phase encode direction, so (e.g.) the subject’s nose may poke into the back of the subject’s head. These are typically: No wrapping - for PET or images that have already been spatially transformed. Also the recommended option if you are not really sure. Wrap in Y - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).
```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
```
Optional weighting image to weight each voxel of the reference image differently when estimating the realignment parameters. The weights are proportional to the inverses of the standard deviations. This would be used, for example, when there is a lot of extra-brain motion - e.g., during speech, or when there are serious artifacts in a particular region of the images.
```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];


%for context, [2 1] corresponds to the definition below. from [1] 
which.labels = {
                ' All Images (1..n)'
                'Images 2..n'
                ' All Images + Mean Image'
                ' Mean Image Only'
}';
which.values = {[2 0] [1 0] [2 1] [0 1]};
```
[1](https://github.com/neurodebian/spm12/blob/master/config/spm_cfg_realign.m)

Specify the images to reslice. All Images (1..n) : This reslices all the images - including the first image selected - which will remain in its original position. Images 2..n : Reslices images 2..n only. Useful for if you wish to reslice (for example) a PET image to fit a structural MRI, without creating a second identical MRI volume. All Images + Mean Image : In addition to reslicing the images, it also creates a mean of the resliced image. Mean Image Only : Creates the mean resliced image only.

```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 1;
```
The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not recommended for image realignment. Trilinear Interpolation is probably OK for PET, but not so suitable for fMRI because higher degree interpolation generally gives better results [111, 112, 113]. Although higher degree methods provide better interpolation, but they are slower because they use more neighbouring voxels. Fourier Interpolation [35, 29] is another option, but note that it is only implemented for purely rigid body transformations. Voxel sizes must all be identical and isotropic.

```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
```
This indicates which directions in the volumes the values should wrap around in. For example, in MRI scans, the images wrap around in the phase encode direction, so (e.g.) the subject’s nose may poke into the back of the subject’s head. These are typically: No wrapping - for PET or images that have already been spatially transformed. Wrap in Y - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).

```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
```
Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).

```matlab
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
```
Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ’r’.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).source(1) = cfg_dep(
    'Realign: Estimate & Reslice: Realigned Images (Sess 1)',
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),
    substruct('.','sess', '()',{1}, '.','cfiles'));

matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';

matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).resample(1) = cfg_dep(
    'Realign: Estimate & Reslice: Realigned Images (Sess 1)', 
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), 
    substruct('.','sess', '()',{1}, '.','cfiles'));

matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(2).source(1) = cfg_dep(
    'Realign: Estimate & Reslice: Realigned Images (Sess 2)', 
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), 
    substruct('.','sess', '()',{2}, '.','cfiles'));

matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(2).wtsrc = '';

matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(2).resample(1) = cfg_dep(
    'Realign: Estimate & Reslice: Realigned Images (Sess 2)', 
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),
    substruct('.','sess', '()',{2}, '.','cfiles'));
```
Dependencies for operations under `spm.tools.oldnorm.etswrite`. These are filled by batch operation at run time.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.template = {'/Users/stevenbishay/Downloads/spm8/templates/PET.nii,1'};
```
Specify a template image to match the source image with. The contrast in the template must be similar to that of the source image in order to achieve a good registration. It is also possible to select more than one template, in which case the registration algorithm will try to find the best linear combination of these images in order to best model the intensities in the source image.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
```
Applies a weighting mask to the template(s) during the parameter estimation. With the default brain mask, weights in and around the brain have values of one whereas those clearly outside the brain are zero. This is an attempt to base the normalisation purely upon the shape of the brain, rather than the shape of the head (since low frequency basis functions can not really cope with variations in skull thickness). The option is now available for a user specified weighting image. This should have the same dimensions and mat file as the template images, with values in the range of zero to one.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
```
Smoothing to apply to a copy of the source image. The template and source images should have approximately the same smoothness. Remember that the templates supplied with SPM have been smoothed by 8mm, and that smoothnesses combine by Pythagoras’ rule.
```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
```
Smoothing to apply to a copy of the template image. The template and source images should have approximately the same smoothness. Remember that the templates supplied with SPM have been smoothed by 8mm, and that smoothnesses combine by Pythagoras’ rule.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
```
Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking). The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). If registering to an image in ICBM/MNI space, then choose the first option. If registering to a template that is close in size, then select the second option. If you do not want to regularise, then choose the third.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
```
Cutoff of DCT bases. Only DCT bases of periods longer than the cutoff are used to describe the warps. The number used will depend on the cutoff and the field of view of the template image(s).

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
```
Number of iterations of nonlinear warping performed.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
```
The amount of regularisation for the nonlinear part of the spatial normalisation. Pick a value around one. However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude) - or even just use an affine normalisation. The regularisation influences the smoothness of the deformation fields.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
```
Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images. Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity.

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.bb = [-85 -120 -90 85 80 95];
```
The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.vox = [1.5 1.5 1.5];
```
The voxel sizes (x, y & z, in mm) of the written normalised images.


```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
```
The method by which the images are sampled when being written in a different space. (Note that Inf or NaN values are treated as zero, rather than as missing data) Nearest Neighbour: - Fastest, but not normally recommended. Trilinear Interpolation: - OK for PET, realigned fMRI, or segmentations B-spline Interpolation: - Better quality (but slower) interpolation [111], especially with higher degree splines. Can produce values outside the original range (e.g. small negative values from an originally all positive image).

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];

% context [1]
wrap.labels  = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z'
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrap.values  = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
```
[1](https://github.com/neurodebian/spm12/blob/master/toolbox/OldNorm/spm_cfg_normalise.m)

These are typically: No wrapping: for PET or images that have already been spatially transformed. Wrap in Y: for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).

```matlab
matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
```
Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ’w’.

```matlab
matlabbatch{3}.spm.spatial.smooth.data(1) = cfg_dep(
    'Old Normalise: Estimate & Write: Normalised Images (Subj 1)',
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),
    substruct('()',{1}, '.','files'));
```
Dependency for operations under `spm.spatial.smooth`

```matlab
matlabbatch{3}.spm.spatial.smooth.fwhm = [8 8 8];
```
Full width at half maximum (FWHM) of the Gaussian smoothing kernel in mm. Three values should be entered, denoting the FWHM in the x, y and z directions.

```matlab
matlabbatch{3}.spm.spatial.smooth.dtype = 0;
```
Data type of the output images. ’SAME’ indicates the same data type as the original images.

```matlab
matlabbatch{3}.spm.spatial.smooth.im = 0;
```
An "implicit mask" is a mask implied by a particular voxel value (0 for images with integer type, NaN for float images). If set to ’Yes’, the implicit masking of the input image is preserved in the smoothed image.

```matlab
matlabbatch{3}.spm.spatial.smooth.prefix = 's';
```
String to be prepended to the filenames of the smoothed image file(s). Default prefix is ’s’.

```matlab
matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep(
    'Old Normalise: Estimate & Write: Normalised Images (Subj 2)', 
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), 
    substruct('()',{2}, '.','files'));
```
Dependencies for operations under `spm.spatial.smooth`

```matlab
matlabbatch{4}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{4}.spm.spatial.smooth.dtype = 0;
matlabbatch{4}.spm.spatial.smooth.im = 0;
matlabbatch{4}.spm.spatial.smooth.prefix = 's';
```
Same as previous block.

```matlab
spm_jobman('run', matlabbatch);
end
```
