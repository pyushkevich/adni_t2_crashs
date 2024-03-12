#!/bin/bash

# Project paths
ROOT=/data/pauly2/t2_crashs

# CRASHS location
CRASHS_HOME=/data/pauly2/tk/crashs

# Add to path
PATH=$ROOT/bin:$PATH

# Set up experiment for particular fold
function set_atlas_subject_vars()
{
    local tracestate id args
    tracestate=$(shopt -po xtrace); set +x

    # Read the fold information
    id=${1?}

    # Restore tracestate
    set +vx; eval $tracestate
}

function set_atlas_subject_side_vars()
{
    local tracestate id side args
    # tracestate=$(shopt -po xtrace); set +x

    # Read the fold information
    id=${1?}
    side=${2?}

    # Read the subject-wide stuff
    set_atlas_subject_vars $id

    # Source data 
    T2_CHUNK_MRI=$ROOT/input/pmc_t2/${id}_tse_native_chunk_${side}.nii.gz
    T2_CHUNK_SEG=$ROOT/input/pmc_t2/${id}_tse_native_chunk_${side}_seg.nii.gz
    T1SR_CHUNK_MRI=$ROOT/input/pmc_t1/${id}/${id}_${side}_T1_SR_reg.nii.gz
    T1SR_CHUNK_SEG=$ROOT/input/pmc_t1/${id}/${id}_${side}_sf_full_SR_manualseg_amyg.nii.gz

    # Directories for this subject
    CRASHS_INPUT_DIR=$ROOT/work/crashs_build/input/$id/$side

    # Upsample inputs
    T2_CHUNK_SEG_GMDG=$CRASHS_INPUT_DIR/${id}_${side}_chunk_gm_dg_seg.nii.gz

    # Upsample model inference on nnunet folder
    # UPSAMPLE_NNUNET_INFER_DIR=$FOLDDIR/upsample/nnunet_infer00

    # Restore tracestate
    # set +vx; eval $tracestate
}

# Set up experiment for a particular subject in a fold. The second parameter is the 
# subdirectory, since for testing we may take the subject's manual segmentation or
# their automatic nnU-net segmentation
function set_fold_subject_vars()
{
    local tracestate id subdir args
    # tracestate=$(shopt -po xtrace); set +x

    id=${1?}
    subdir=${2?}
    WDIR=$FOLDDIR/${subdir}/${id}
    IDIR=$ROOT/input/${id}
    FAKE_ASHS_DIR=$WDIR/fake_ashs

    side=$(echo $id | sed -e "s/.*L$/left/" -e "s/.*R$/right/")

    # set +vx; eval $tracestate
}

# Import training data, ground truth, etc
function import_inputs()
{
    # Copy the data for building CRASHS template 
    mkdir -p $ROOT/input/pmc_t1/
    rsync -av bscsub1.pmacs.upenn.edu:/project/hippogang_2/pauly/wolk/anterior_mtl/atlas_set/work/ $ROOT/input/pmc_t1/

    # Copy the T2 images and segmentations
    rsync -avL bscsub1.pmacs.upenn.edu:/project/hippogang_2/pauly/wolk/upenn_pmc_2021/chunks/ $ROOT/input/pmc_t2/
}


# Command to post-process one upsampling result, creating a CRASHS input folder
function postprocess_upsample()
{
    # The ID of the file
    local id side sid
    id=${1?}
    side=${2?}

    # Combined id
    sid=${id}_${side}

    # ASHS-T2 segmentation in native space
    SEG_T2_ASHS=${3?}

    # Segmentation of the DG and GM in native space
    SEG_DG_GM=${4?}

    # Segmentation of the WM in T1 space
    SEG_WM_T1=${5?}

    # The work directory
    WDIR=${6?}

    # Go to the work directory
    pushd $WDIR

    # Define "fake ASHS dir"
    FAKE_ASHS_DIR=$WDIR/fake_ashs

    # At what level do we trim the GM probability?
    L=0.7

    # Now further group into just CA, Sub+ERC and Cortex (for finding fused CA/PHC)
    c3d $SEG_T2_ASHS -replace 1 5 2 5 3 5 4 5 7 0 8 7 10 7 11 1 12 1 13 1 14 0 \
        -type uchar -o ${sid}_ivseg_hc.nii.gz

    # Split the upsampled image into GM and DG components
    c3d -mcs ${sid}_ivseg_unet_upsample.nii.gz -oo ${sid}_ivseg_unet_upsample_gm.nii.gz ${sid}_ivseg_unet_upsample_dg.nii.gz

    # Register the GM component to the in vivo GM
    c3d $SEG_DG_GM -thresh 1 1 1 0 -o ${sid}_ivseg_gm.nii.gz
    greedy -d 3 -a -dof 6 -n 100x100 -m SSD -ia-identity -threads 1 \
        -i ${sid}_ivseg_gm.nii.gz ${sid}_ivseg_unet_upsample_gm.nii.gz \
        -o ${sid}_ivseg_unet_upsample_rigid.mat
    greedy -d 3  -threads 1 \
        -rf ${sid}_ivseg_unet_upsample_gm.nii.gz \
        -rm ${sid}_ivseg_unet_upsample_gm.nii.gz ${sid}_ivseg_unet_upsample_gm_shift.nii.gz \
        -rm ${sid}_ivseg_unet_upsample_dg.nii.gz ${sid}_ivseg_unet_upsample_dg_shift.nii.gz \
        -r ${sid}_ivseg_unet_upsample_rigid.mat    

    # Generate the band to add to white matter where hippocampus and cortex touch
    c3d ${sid}_ivseg_unet_upsample_gm_shift.nii.gz -thresh 0.7 inf 1 0 -as X ${sid}_ivseg_hc.nii.gz \
        -int 0 -reslice-identity -as Y -thresh 1 1 1 0 -sdt \
        -push Y -thresh 5 5 1 0 -sdt \
        -push Y -thresh 7 7 1 0 -sdt \
        -foreach -scale -1 -endfor \
        -vote -shift 1 -push X -times -as Z -thresh 1 1 1 0 -as A -push Z -thresh 2 2 1 0 -as B \
        -dilate 1 1x1x1 -times -push A -dilate 1 1x1x1 -push B -times -add \
        -o ${sid}_hc_overlap.nii.gz

    # Generate probability maps for GM, WM and CSF, zeroing out the hippocampus labels
    mkdir -p ./mtl_only

    # Reslice the WM segmentation into an upsampled T2 segmentation space and clean up,
    # getting rid of WM pixels that are too far from any GM, as well as any CSF stuck 
    # between GM and WM
    set -x -e
    c3d ${sid}_ivseg_unet_upsample_gm_shift.nii.gz -thresh ${L} inf 1 0 -as G \
        $SEG_WM_T1 -thresh 20 20 1 0 -smooth 0.2mm -reslice-identity -thresh 0.5 inf 1 0 -popas WG \
        -push G -stretch 0 1 1 0 -push WG -fm 5.0 -thresh 0 3 1 0 -o test_wm_pad.nii.gz \
        -push G -stretch 0 1 1 0 -times -as W -push G -fm 20.0 -thresh 0 5 1 0 -push W -times \
        -o ${sid}_wm_upsample_shift_pad_trim.nii.gz

    # Generate four probability maps: GM, WM and CSF. 
    c3d \
        -verbose ${sid}_ivseg_unet_upsample_gm_shift.nii.gz -pad 10x10x10 10x10x10 0 -as Q \
        -thresh 0 ${L} 1 0 -popas M -push Q -stretch 0 ${L} 0 0.5 -clip 0 0.5 -push M -times \
        -push Q -stretch ${L} 1 0.5 1.0 -clip 0.5 1.0 -push M -stretch 0 1 1 0 -times -max -info \
        -as B -dup ${sid}_hc_overlap.nii.gz -int 0 -reslice-identity -as O -stretch 0 1 1 0 -times -popas C \
        -push B -push C -push B -scale -1 -add -smooth 0.5mm -add -clip 0 1 \
        -as PG -o ./mtl_only/${sid}_prob_gm.nii.gz -stretch 0 1 1 0 -as PNG \
        -push Q ${sid}_wm_upsample_shift_pad_trim.nii.gz -reslice-identity -push O -max -as W \
        -push PNG -times -o ./mtl_only/${sid}_prob_wm.nii.gz \
        -push PNG -push W -stretch 0 1 1 0 -times -o ./mtl_only/${sid}_prob_csf.nii.gz    

    # For CRASHS we only care about assigning the gray matter labels to the gray matter probability
    # map, all the other labels 
    mkdir -p ./mtl_only/fake_ashs/bootstrap/fusion
    c3d -verbose ./mtl_only/${sid}_prob_gm.nii.gz -popas GM ./mtl_only/${sid}_prob_wm.nii.gz -popas WM ./mtl_only/${sid}_prob_csf.nii.gz -popas CSF \
        -push GM $SEG_T2_ASHS -int 0 -reslice-identity \
        -push WM -thresh 0.5 inf 15 0 -max -popas S \
        $(for ((i=0;i<=15;i++)); do echo "-push S -thresh $i $i 1 0"; done) \
        -push CSF -push GM -push GM -push CSF -push GM \
        -push CSF -push CSF -push CSF -push GM -push CSF \
        -push GM -push GM -push GM -push GM -push CSF -push WM \
        -foreach-comp 16 -as P -thresh 0.5 inf 1 0 -times -insert P 1 -fast-marching 20 -reciprocal -endfor \
        -foreach -scale 10 -endfor -softmax -omc ./mtl_only/${sid}_ivseg_ashs_upsample_posteriors.nii.gz \
        -oo ./mtl_only/fake_ashs/bootstrap/fusion/posterior_corr_usegray_${side}_%03d.nii.gz \
        -vote -o ./mtl_only/${sid}_ivseg_ashs_upsample.nii.gz 


<<'SKIP_BAND'
    # Find the interface between CA and DG and define this as the artificial WM label
    c3d \
        ${id}_ivseg_unet_upsample_gm_shift.nii.gz -as G -stretch 1 ${L} -1 0 \
        ${id}_ivseg_unet_upsample_dg_shift.nii.gz -as D -stretch 0 1 1 -1 \
        -levelset-curvature 0.2 -levelset 20 -stretch -1 1 1 0 \
        -thresh 0.5 inf 1 0 -o ${id}_ivseg_unet_upsample_dg_expand.nii.gz -as DX -dilate 0 2x2x2 -stretch 0 1 1 0 -push DX -times \
        -push G -thresh 0.5 inf 1 0 -dilate 1 2x2x2 -times -o ${id}_ivseg_unet_upsample_dg_band.nii.gz 

    # Generate four probability maps: GM, DG, WM and CSF. 
    c3d \
        -verbose ${id}_ivseg_unet_upsample_gm_shift.nii.gz -pad 10x10x10 10x10x10 0 -as Q \
        -thresh 0 ${L} 1 0 -popas M -push Q -stretch 0 ${L} 0 0.5 -clip 0 0.5 -push M -times \
        -push Q -stretch ${L} 1 0.5 1.0 -clip 0.5 1.0 -push M -stretch 0 1 1 0 -times -max -info \
        -as B -dup ${id}_hc_overlap.nii.gz -int 0 -reslice-identity -as O -stretch 0 1 1 0 -times -popas C \
        -push B -push C -push B -scale -1 -add -smooth 0.5mm -add -clip 0 1 \
        -as PG -o ${id}_prob_gm.nii.gz -stretch 0 1 1 0 -as PNG \
        ${id}_ivseg_unet_upsample_dg_shift.nii.gz ${id}_ivseg_unet_upsample_dg_expand.nii.gz -times \
        ${id}_ivseg_unet_upsample_dg_band.nii.gz -stretch 0 1 1 0 -times \
        -reslice-identity -push PNG -min -o ${id}_prob_dg.nii.gz -push PG -max -stretch 0 1 1 0 -as PNG \
        -shift -0.25 \
        -dup $SEG_WM_T1 -thresh 20 20 1 0 -dilate 1 1x1x0 -reslice-identity -stretch 0 1 1 -1 \
        -as A -levelset-curvature 0.8 -levelset 100 \
        -thresh -inf 0 1 0 -dilate 1 2x2x2 -push A -thresh -inf 0 1 0 -dilate 1 3x3x3 -times -as W \
        -push Q ${id}_ivseg_unet_upsample_dg_band.nii.gz -reslice-identity -max -as W \
        -push PNG -times -o ${id}_prob_wm.nii.gz \
        -push PNG -push W -stretch 0 1 1 0 -times -o ${id}_prob_csf.nii.gz

    # Within each probability map, vote between corresponding labels 
    mkdir -p $FAKE_ASHS_DIR/bootstrap/fusion
    c3d -verbose ${id}_prob_gm.nii.gz -popas GM ${id}_prob_dg.nii.gz -popas DG ${id}_prob_wm.nii.gz -popas WM ${id}_prob_csf.nii.gz -popas CSF \
        -push GM $SEG_T2_ASHS -int 0 -reslice-identity \
        -push GM $SEG_WM_T1 -thresh 20 20 15 0 -int 0 -reslice-identity -add -clip 0 15 \
        -push GM ${id}_ivseg_unet_upsample_dg_band.nii.gz -replace 1 15 -reslice-identity -max -popas S \
        $(for ((i=0;i<=15;i++)); do echo "-push S -thresh $i $i 1 0"; done) \
        -push CSF -push GM -push GM -push DG -push GM \
        -push CSF -push CSF -push CSF -push GM -push CSF \
        -push GM -push GM -push GM -push GM -push CSF -push WM \
        -foreach-comp 16 -as P -thresh 0.5 inf 1 0 -times -insert P 1 -fast-marching 20 -reciprocal -endfor \
        -softmax -omc ${id}_ivseg_ashs_upsample_posteriors.nii.gz \
        -oo $FAKE_ASHS_DIR/bootstrap/fusion/posterior_corr_usegray_${side}_%03d.nii.gz \
        -vote -o ${id}_ivseg_ashs_upsample.nii.gz 

    # Copy the affine matrix to ASHS dir
    # mkdir -p $FAKE_ASHS_DIR/affine_t1_to_template/
    # cp $ROOT/input/${id}/t1_to_template_affine.mat $FAKE_ASHS_DIR/affine_t1_to_template/

SKIP_BAND
   popd
}


# Upsample one subject/side
function upsample_atlas_qsub()
{
    id=${1?}
    side=${2?}
    set_atlas_subject_side_vars $id $side

    mkdir -p $CRASHS_INPUT_DIR

<<'SKIP'
    # Combine labels for upsampling
    c3d $T2_CHUNK_SEG \
        -replace 7 0 14 0 2 1 3 2 4 1 8 1 10 1 11 1 12 1 13 1 \
        -o $T2_CHUNK_SEG_GMDG

    # Perform upsampling
    python /data/pauly2/ashs_xv/scripts/upsample_net.py apply \
        -t /data/pauly2/ashs_xv/work/fold_5/upsample/train00 \
        -o $CRASHS_INPUT_DIR \
        -g $T2_CHUNK_MRI \
        -s $T2_CHUNK_SEG_GMDG \
        -i ${id}_${side}
SKIP

    # Run postprocessing script
    postprocess_upsample ${id} ${side} $T2_CHUNK_SEG $T2_CHUNK_SEG_GMDG $T1SR_CHUNK_SEG $CRASHS_INPUT_DIR

    # Copy the affine matrix to ASHS dir
    mkdir -p $CRASHS_INPUT_DIR/fake_ashs/affine_t1_to_template/
    cp $ROOT/input/pmc_t2/${id}_t1_to_template_affine.mat $FAKE_ASHS_DIR/affine_t1_to_template/t1_to_template_affine.mat

    mkdir -p $CRASHS_INPUT_DIR/mtl_only/fake_ashs/affine_t1_to_template/
    cp $ROOT/input/pmc_t2/${id}_t1_to_template_affine.mat $CRASHS_INPUT_DIR/mtl_only/fake_ashs/affine_t1_to_template/t1_to_template_affine.mat
}


# Perform inference and post-processing for a fold
function upsample_atlas()
{
    pushd $ROOT

    # For each subject and side, run upsampling
    for id in $(awk -F, '{print $1}' $ROOT/manifest/atlas_ids.csv); do
        for side in left right; do
            upsample_atlas_qsub $id $side
        done
    done

    popd 
}


function apply_crashs_qsub()
{
    expid=${1?}
    id=${2?}
    side=${3?}
    set_atlas_subject_side_vars $id $side

    pushd $ROOT

    rm -rf $ROOT/tmp/crashs_test_${expid}_${id}_${side}
    mkdir -p $ROOT/tmp/crashs_test_${expid}_${id}_${side} 

    # Perform CRASHS inference on this case
    python \
        $CRASHS_HOME/crashs.py \
        $CRASHS_INPUT_DIR/$expid/fake_ashs \
        $ROOT/manifest/template_init_dir/$expid \
        $ROOT/tmp/crashs_test_${expid}_${id}_${side} \
        -i $id -s $side --keops --device cuda

    popd
}

function build_template()
{
    expid=${1?}
    pushd $ROOT
    
    # Create a build directory for this template
    CDIR=$ROOT/work/crashs_build/$expid
    rm -rf $CDIR
    mkdir -p $CDIR/manifest $CDIR/work

    # Create a manifest for this template
    python $ROOT/scripts/select_random_training_set.py --seed 42 -n 24 $ROOT/manifest/atlas_ids.csv $CDIR/manifest/crashs_atlas_ids.csv

    # Generate a JSON of crashs build inputs
    while IFS=, read -r id side; do
        set_atlas_subject_side_vars $id $side
        # jq -n --arg id ${id}_${side} --arg side ${side} --arg path="$CRASHS_INPUT_DIR" '{id: $id, side: $side}'
        jq -n --arg id ${id}_${side} --arg side $side \
            --arg ashs "$CRASHS_INPUT_DIR/$expid/fake_ashs" \
            '{id: $id, side: $side, path: $ashs}'
    done < $CDIR/manifest/crashs_atlas_ids.csv | jq -n '. |= [inputs]' > $CDIR/manifest/crashs_inputs.json 
    cat $CDIR/manifest/crashs_inputs.json

    # Clear out the CRASHS build folder
    rm -rf $CDIR/work

    # Run the template build code using the specified JSON
    python $CRASHS_HOME/build_template.py \
        $ROOT/manifest/template_init_dir/${expid} \
        $CDIR/manifest/crashs_inputs.json \
        $CDIR/work $CDIR/crashs_template

    popd 
}

function usage()
{
  echo "runme.sh : T2-CRASHS scripts for ADNI3"
  echo "Usage:"
  echo "  runme.sh [options] <function> [args]"
  echo "Options:"
  echo "  -d            : Turn on command echoing (for debugging)"
  echo "  -s <level>    : Skip expensive operations if results exist"
  echo "Primary functions:"
  echo "  import_inputs                           : Copy inputs"
  echo "  apply_upsample                          : Perform inference with upsampling network"
  echo "  build_template                          : Build the CRASHS template"
}

# Read the command-line options
while getopts "dhs:" opt; do
  case $opt in
    d) set -x -e;;
    h) usage; exit 0;;
    s) SKIPLEVEL=$OPTARG;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;
  esac
done

# Get remaining args
shift $((OPTIND - 1))

# No parameters? Show usage
if [[ "$#" -lt 1 ]]; then
  usage
  exit 255
fi

# Main entrypoint into script
COMMAND=$1
shift
$COMMAND "$@"