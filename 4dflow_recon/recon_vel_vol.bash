
#!/usr/bin/env bash



# RECON VEL VOLUME

# e.g., bash recon_vel_vol.bash ~/path/to/top/level/recon/directory/ vel_vol


# Input 

RECONDIR=$1
VOLDESC=$2


# Check that Recon Directory Exists

if [ ! -d "$RECONDIR" ]; then
  echo directory $RECONDIR does not exist
  exit 1
else


# Manage Paths to Allow Queueing of Jobs using TaskSpooler

ORIG_PATH=$(pwd)
SCRIPT_PATH=$(dirname `which $0`)

RECONVOLDIR=$RECONDIR/$VOLDESC
mkdir -p $RECONVOLDIR
cd $RECONVOLDIR

echo RECON CINE VOLUME
echo $RECONVOLDIR


# Variables 

RECON=$VOLDESC.nii.gz
STACKS="../data/s*_rlt_ph_corr_uterus.nii.gz"
GRADMOMVALSFILE="../data/grad_moment_vals.nii.gz"
GRADMOMDIRSFILE="../data/grad_moment_dirs.nii.gz"
THICKNESSFILE="../data/slice_thickness.txt"
EXCLUDECINEVOLFILE="../data/force_exclude_cine_vol.txt"
RESOLUTION=1.25
NMC=0
NSR=40
NSRLAST=20
NUMCARDPHASE=25
ALPHA=3
STACKDOFDIR="../dc_vol/stack_transformations"
CINETRANSDIR="../cine_vol/transformations"
MASKCINEVOL="../cine_vol/mask_cine_vol.nii.gz"
MEANRRFILE="../cardsync/mean_rrinterval.txt"
RRINTERVALSFILE="../cardsync/rrintervals.txt"
CARDPHASESFILE="../cardsync/cardphases_interslice_cardsync.txt"


# Setup

ITER=$(($NMC+1))
THICKNESS=$(cat $THICKNESSFILE)
MEANRR=$(cat $MEANRRFILE)
RRINTERVALS=$(cat $RRINTERVALSFILE)
CARDPHASES=$(cat $CARDPHASESFILE)
NUMSTACK=$(ls -1 ../data/s*_dc_ab.nii.gz | wc -l);
NUMSLICE=$(eval "wc -w $RRINTERVALSFILE | awk -F'[ ]' '{print \$1}'" )
NUMFRAME=$(eval "wc -w $CARDPHASESFILE | awk -F'[ ]' '{print \$1}'" )

EXCLUDECINEVOL=$(cat $EXCLUDECINEVOLFILE)
NUMEXCLUDECINEVOL=$(eval "wc -w $EXCLUDECINEVOL | awk -F'[ ]' '{print \$1}'" )
EXCLUDESTACK=$(cat $EXCLUDESTACKFILE)
NUMEXCLUDESTACK=$(eval "wc -w $EXCLUDESTACKFILE | awk -F'[ ]' '{print \$1}'" )
EXCLUDESLICE=$(cat $EXCLUDESLICEFILE)
NUMEXCLUDESLICE=$(eval "wc -w $EXCLUDESLICEFILE | awk -F'[ ]' '{print \$1}'" )
EXCLUDEFRAME=$(cat $EXCLUDEFRAMEFILE)
NUMEXCLUDEFRAME=$(eval "wc -w $EXCLUDEFRAMEFILE | awk -F'[ ]' '{print \$1}'" )


# Recon Velocity Cine Volume

echo reconstructing cine volume: $RECON
CMD="reconstructCardiacVelocity $RECON $NUMSTACK $STACKS -thickness $THICKNESS -dofin $STACKDOFDIR/stack-transformation*.dof -transformations $CINETRANSDIR -mask $MASKCINEVOL -no_stack_intensity_matching -no_robust_statistics -alpha $ALPHA -limit_intensities -iterations $ITER -rec_iterations $NSR -rec_iterations_last $NSRLAST -resolution $RESOLUTION -force_exclude_stack $NUMEXCLUDESTACK $EXCLUDESTACK -force_exclude_sliceloc $NUMEXCLUDESLICE $EXCLUDESLICE -force_exclude $NUMEXCLUDECINEVOL $EXCLUDECINEVOL -numcardphase $NUMCARDPHASE -rrinterval $MEANRR -rrintervals $NUMSLICE $RRINTERVALS -cardphase $NUMFRAME $CARDPHASES -debug > log-main.txt"
echo $CMD > recon.bash
eval $CMD


# Clean Up

echo "" >> recon.bash
CMD="mkdir transformations; mv transformation*.dof transformations;"
echo $CMD >> recon.bash
eval $CMD
echo "" >> recon.bash
CMD="mkdir sr_iterations; mv *_mc*sr* sr_iterations;"
echo $CMD >> recon.bash
eval $CMD


# Finish

echo "volume reconstruction complete"

cd $ORIG_PATH


fi

