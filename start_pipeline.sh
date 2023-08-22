###############################################################################
### COLORS ###
##############
red="\033[1;31m"   # red
green="\033[1;32m" # green
ylo="\033[1;33m"   # yellow
blue="\033[1;34m"  # blue
nc="\033[0m"       # no color

###############################################################################
### ABOUT ###
#############
workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd) # Get working directory
workflow_base_version="2023.08"                        # Workflow base version

###############################################################################
### OPERATING SYSTEM ###
########################
echo -e "
${green}------------------------------------------------------------------------${nc}
${green}#####${nc} ${red}OPERATING SYSTEM${nc} ${green}#####${nc}
${green}----------------------------${nc}
"

# Get and print operating system 
case "$OSTYPE" in
  linux*)   os="linux" ;;
  bsd*)     os="bsd" ;;
  darwin*)  os="osx" ;;
  solaris*) os="solaris" ;;
  msys*)    os="windows" ;;
  cygwin*)  os="windows" ;;
  *)        os="unknown (${OSTYPE})" ;;
esac
echo -e "${blue}Operating system${nc} _______ ${red}${os}${nc}"

# Get and print shell
shell=$SHELL
echo -e "${blue}Shell${nc} __________________ ${ylo}${shell}${nc}"





###############################################################################
### WORKFLOW-BASE INSTALLATION ###
##################################

# Delete previous session
find results -type f -delete
find resources/reads -name "*.index" -delete
find resources/reads -name "*.index.fai" -delete
find resources/reads -name "*.index.gzi" -delete
find resources/reads -name "*.index.readdb" -delete

echo -e "
${green}------------------------------------------------------------------------${nc}
${green}#####${nc} ${red}WORKFLOW-BASE INSTALLATION${nc} ${green}#####${nc}
${green}---------------------------------------${nc}
"


# Test if latest 'workflow-base' environment exist
if [[ $(conda info --envs | grep -o -E "^workflow-base_v.${workflow_base_version}") ]]
then
    echo -e "
Conda environment ${ylo}workflow-base_v.${workflow_base_version}${nc} it's already created!
"
else
    echo -e "
Conda environment ${red}workflow-base_v.${workflow_base_version}${nc} will be now created!
"
    conda env create -f ${workdir}/workflow/environment/${os}/workflow-base_v.${workflow_base_version}.yaml
fi


###############################################################################
### CONDA ACTIVATION ###
########################
echo -e "
${green}------------------------------------------------------------------------${nc}
${green}#####${nc} ${red}CONDA ACTIVATION${nc} ${green}#####${nc}
${green}----------------------------${nc}
"

echo -e "conda activate ${ylo}workflow-base_v.${workflow_base_version}${nc}"
# intern shell source conda
source ~/miniconda3/etc/profile.d/conda.sh 2> /dev/null          # local user
source /usr/local/miniconda3/etc/profile.d/conda.sh 2> /dev/null # HPC server
conda activate workflow-base_v.${workflow_base_version}          # conda activate workflow-base

# ###############################################################################
# ### SNAKEMAKE ACTIVATION ###
# ########################
# echo -e "
# ${green}------------------------------------------------------------------------${nc}
# ${green}#####${nc} ${red}SNAKEMAKE ACTIVATION${nc} ${green}#####${nc}
# ${green}----------------------------${nc}
# "
# # Test if snakemake environment exist
# if [[ $(conda info --envs | grep -o -E "^snakemake") ]]
# then
#     echo -e "
# Snakemake it's already created!
# "
# else
#     echo -e "
# Snakemake will be now created!
# "
#     conda env create -n snakemake
# fi

# echo -e "conda activate snakemake"
# conda activate snakemake        # conda activate snakemake

###############################################################################
### SETTINGS ###
################
echo -e "
${green}------------------------------------------------------------------------${nc}
${green}#####${nc} ${red}SETTINGS${nc} ${green}#####${nc}
${green}--------------------${nc}
"

conda_version=$(conda --version | sed 's/conda //')                   # Conda version (ver. >= 23.3.1)
snakemake_version=$(snakemake --version)                              # Snakemake version (ver. 7.25.0)
mamba_version=$(mamba --version | sed 's/mamba //' | head -n 1)       # Mamba version (ver. 1.4.2)

fastq=$(expr $(ls -l ${workdir}/resources/reads/*.fastq.gz 2> /dev/null | wc -l)) # Get fastq.gz files count
if [[ "${fastq}" == "0" ]]                                                        # If no sample,
then                                                                               # start GeVarLi with at least 1 sample
    echo -e "${red}¡${nc} No fastq file detected in ${ylo}resources/reads/${nc} ${red}!${nc}
${red}¡${nc} Please, add at least 1 sample in ${ylo}resources/reads/${nc} ${red}!${nc}
"
    exit 1
fi

config_file="${workdir}/config.yaml"       # Get configuration file
# conda_frontend=$(yq e '.conda.frontend' ${config_file} | sed 's/^- //' | sed 's/\[\"//' | sed 's/\"\]//') # Get user config: conda frontend
conda_frontend=$(yq -r '.conda.frontend' ${config_file})
max_threads=$(yq -r '.resources.cpus' ${config_file})    # Get user config: max threads






##############################################################################
## SNAKEMAKE PIPELINES ###
##########################
echo -e "
${green}------------------------------------------------------------------------${nc}
${green}#####${nc} ${red}SNAKEMAKE PIPELINES${nc} ${green}#####${nc}
${green}-------------------------------${nc}
"

snakefiles_list="quality_control pipeline"

echo -e "
${blue}## Unlocking Working Directory ##${nc}
${blue}---------------------------------${nc}
"
conda activate snakemake
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# Remove a lock on the working directory.

for snakefile in ${snakefiles_list} ; do
    echo -e "${blue}-- ${snakefile} --${nc}" ;
    snakemake \
	--directory ${workdir}/ \
        --snakefile ${workdir}/workflow/snakefiles/${snakefile}.smk \
        --config os=${os} \
        --rerun-incomplete \
        --unlock ;
done

# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# List all conda environments and their location on disk.

for snakefile in ${snakefiles_list} ; do
    echo -e "${blue}-- ${snakefile} --${nc}" ;
    snakemake \
        --directory ${workdir}/ \
        --snakefile ${workdir}/workflow/snakefiles/${snakefile}.smk \
        --cores ${max_threads} \
        --config os=${os} \
        --rerun-incomplete \
        --list-conda-envs ;
done

echo -e "
${blue}## Conda Environments Setup ##${nc}
${blue}------------------------------${nc}
"
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# If defined in the rule, run job in a conda environment.
# If mamba package manager is not available, or if you still prefer to use conda, you can enforce that with this setting.
## Default "mamba", recommended because much faster !
# If specified, only creates the job-specific conda environments then exits. The –use-conda flag must also be set.

for snakefile in ${snakefiles_list} ; do
    echo -e "${blue}-- ${snakefile} --${nc}" ;
    snakemake \
        --directory ${workdir}/ \
        --snakefile ${workdir}/workflow/snakefiles/${snakefile}.smk \
        --cores ${max_threads} \
        --config os=${os} \
        --rerun-incomplete \
        --use-conda \
        --conda-frontend ${conda_frontend} \
        --conda-create-envs-only ;
done

echo -e "
${blue}## Dry Run ##${nc}
${blue}-------------${nc}
"
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# If defined in the rule, run job in a conda environment.
# If mamba package manager is not available, or if you still prefer to use conda, you can enforce that with this setting.
## Default "mamba", recommended because much faster !
# Tell the scheduler to assign creation of given targets (and all their dependencies) highest priority.
# Do not execute anything, and display what would be done. If very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.
# Do not output any progress or rule information.

for snakefile in ${snakefiles_list} ; do
    echo -e "${blue}-- ${snakefile} --${nc}" ;
    snakemake \
        --directory ${workdir}/ \
        --snakefile ${workdir}/workflow/snakefiles/${snakefile}.smk \
        --cores ${max_threads}\
        --config os=${os} \
        --rerun-incomplete \
        --use-conda \
        --conda-frontend ${conda_frontend} \
        --dry-run \
        --quiet ;
done

echo -e "
${blue}## Let's Run! ##${nc}
${blue}----------------${nc}
"
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Define a global maximum number of threads available to any rule. Snakefiles requesting more threads will have their values reduced to the maximum. 
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# Go on with independent jobs if a job fails.
# If defined in the rule, run job in a conda environment.
# If mamba package manager is not available, or if you still prefer to use conda, you can enforce that with this setting.
## Default "mamba", recommended because much faster !
# Tell the scheduler to assign creation of given targets (and all their dependencies) highest priority.
# Print out the shell commands that will be executed.

for snakefile in ${snakefiles_list} ; do
    echo -e "${blue}-- ${snakefile} --${nc}" ;
    snakemake \
        --directory ${workdir}/ \
        --snakefile ${workdir}/workflow/snakefiles/${snakefile}.smk \
        --cores ${max_threads} \
        --max-threads ${max_threads} \
        --config os=${os} \
        --rerun-incomplete \
        --keep-going \
        --use-conda \
        --conda-frontend ${conda_frontend} \
        --printshellcmds ;
done
