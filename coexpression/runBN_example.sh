
    /hpc/users/akersn01/software/multiscale-run_bn-kipp/BN_driver_extra.sh \
    -d /sc/orga/scratch/akersn01/starnet2/bayesian/BigNet/BigNet.disc.MAM \
    -b /hpc/users/akersn01/software/multiscale-run_bn-kipp/ \
    -o /sc/orga/scratch/akersn01/starnet2/bayesian/BigNet/MAMnet \
    -e /sc/orga/scratch/akersn01/starnet2/bayesian/eQTLS_MAMorAOR.list \
    -E /sc/orga/scratch/akersn01/starnet2/bayesian/eQTLS_MAMorAOR.list \
    -q alloc \
    -A acc_apollo \
    -W 12:00 \
    -m 5000 \
    -C TRUE  \
    -w /sc/orga/scratch/akersn01/starnet2/bayesian/BigNet/BigNet.cont.MAM


# -d: discrete gene expr data
# -w continuous gene expr data
# -o output
# -e eQTLs (just a list of genes with eQTLs I beleive.  
# -E same as above, I don't remember why both are needed

# Joe and Ariella had everything set up to run automatically.  everything submits, then when things finish, the next step submits.  
# I personally preferred to have things a bit more controlled, so it's in two steps.  
#  0) run this script
#  1) run ${output_path}/submit_jobs.sh
#  2) run ${output_path}/final_step.sh
