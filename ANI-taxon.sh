#!/bin/bash

set -e
# Parameters: $1 - taxon ID or name to download type strain sequences;
#                  "off" if you want no sequences to be downloaded and want to compare local files only;
#             $2 - output directory
#             $3 - A file containing paths to input fasta files: genomic (.fna) sequences, one per line
#             $4 - path to fastANI executable
#             $5 - Number of threads to use. Default: 1
#             Dependencies:
#                  ncbi datasets
#                  fastANI

VERSION='1.0.a'
COLOR='\x1B[0;33m'
RESET_COLOR='\x1B[0m'

echo " "
echo "   _____            _____  ______  _____             _             _     "
echo "  / ____|    /\    / ____||  ____||  __ \           | |           | |    "
echo " | |        /  \  | |  __ | |__   | |__) |  ______  | |      __ _ | |__  "
echo " | |       / /\ \ | | |_ ||  __|  |  _  /  |______| | |     / _\` || '_ \ "
echo " | |____  / ____ \| |__| || |____ | | \ \           | |____| (_| || |_) |"
echo "  \_____|/_/    \_\\\\_____||______||_|  \_\          |______|\__,_||_.__/ "
echo "                                                                         "
echo " ------------------------------------------------------------------------"
echo " - presents - the program ANI-taxon! Version ${VERSION}"
echo " - $(date)"
echo " ------------------------------------------------------------------------"

if [[ -n "$1" ]]; then
    TAXON_TO_DOWNLOAD="${1}"
    echo " "
    echo -e " ${COLOR} Welcome to our crimson submarine! ${RESET_COLOR}"
    echo -e " ${COLOR} Please take a deep breath: we are about to begin! ${RESET_COLOR}"
    echo "----------------------------------------------------"
    echo " "
    echo " You have chosen Taxon ID / Name: '${TAXON_TO_DOWNLOAD}'"
else
    echo "Error: Taxon ID / Name (\$1) is not specified!"
    exit 1
fi


if [[ -n "$2" ]]; then
    WORKDIR_ROOT=$(realpath $2)
    echo " Output directory: '${WORKDIR_ROOT}'"
else
    echo " Error: output directory (\$2) if not specified!"
    exit 1
fi

if [[ -n "$3" ]]; then
    QUERY_LIST_FILE="${3}"
    echo " Query fna/faa file: '${QUERY_LIST_FILE}'"
else
    echo " Error: query sequence file in fna or faa format (\$3) is not specified!"
    exit 1
fi
while read query_file; do
    if [[ "${query_file}" != *.fna ]]; then
        echo " Error: query file extension is inappropriate: '${query_file}'"
        echo " The extension must be '.fna'. It is indeed important."
        echo ' Please provide the correct file.'
        echo ' Or, if you are sure that the file is correct, just rename it.'
        exit 1
    fi
done < "${QUERY_LIST_FILE}"

if [[ -n "$4" ]]; then
    FASTANI="${4}"
else
    echo "Error: path to fastANI executable (\$4) is not specified!"
    echo 'Please download it here: https://github.com/ParBLiSS/FastANI'
    exit 1
fi

THREADS_NUM=1
if [[ -n "$5" ]]; then
    THREADS_NUM="${5}"
    if [[ ! ${THREADS_NUM} =~ ^[0-9]+$ ]]; then
        echo "Error: invalid number of threads provided: '${THREADS_NUM}'"
        echo 'It must be a positive integer number'
        exit 1
    fi
    max_threads_num=$(nproc)
    if [[ ${THREADS_NUM} > ${max_threads_num} ]]; then
        echo "Error: number of threads provided (${THREADS_NUM}) is greater than the computer has: ${max_threads_num}"
        exit 1
    fi
fi

set -eu

echo 'Run parameters:'
echo " - Taxon ID / Name: '${TAXON_TO_DOWNLOAD}'"
echo " - Output directory: '${WORKDIR_ROOT}'"
echo " - Query list file path: '${QUERY_LIST_FILE}'"
echo " - fastANI path: '${FASTANI}'"
echo " - Number of CPU threads to use: '${THREADS_NUM}'"
echo '--------------------------'



# Check dependencies
dependencies=( datasets )
for exe_file_name in "${dependencies[@]}"; do
    if [[ -z $(which "${exe_file_name}" ) ]]; then
        echo "Error: cannot find program '${exe_file_name}'"
        exit 1
    fi
done


# Make useful variables
TAXON_NO_SPACES=${TAXON_TO_DOWNLOAD/ /_}
datadir="${WORKDIR_ROOT}/${TAXON_NO_SPACES}_db"
fastANI_work_dir="${datadir}/fastANI_workdir"
fastANI_input_list_file="${fastANI_work_dir}/input_file_list.txt"
raw_result_tsv="${WORKDIR_ROOT}/raw_ani_result.tsv"
result_tsv="${WORKDIR_ROOT}/ani_result.tsv"
log_file="${WORKDIR_ROOT}/result.log"
seq_data_dir="${datadir}/ncbi_dataset/data"
taxonomy_file="${WORKDIR_ROOT}/taxonomy.tsv"

# Create necessary dirs and check paths
echo " "
if [[ -d "${datadir}" ]]; then
    echo "Directory ${datadir} already exists"
else
    mkdir -pv "${datadir}"
fi
if [ -d "${fastANI_work_dir}" ]; then
    echo "Directory ${fastANI_work_dir} already exists"
else
    mkdir -pv "${fastANI_work_dir}"
fi

# Create or empty log file
echo -n '' > "${log_file}"

{
    if [[ "${TAXON_TO_DOWNLOAD}" != 'off' ]]; then
        echo "----------------------------------------------"
        echo -e " ${COLOR} $(date) -- Downloading genomes of type strains of ${TAXON_TO_DOWNLOAD} ${RESET_COLOR}"
        echo "----------------------------------------------"

        zip_file="${WORKDIR_ROOT}/${TAXON_NO_SPACES}_RefSeq_${MODE}.zip"

        # useful options:
        # --assembly-level
        # --assembly-source
        # --from-type
        datasets download genome taxon "${TAXON_TO_DOWNLOAD}" \
            --include genome \
            --from-type \
            --dehydrated \
            --assembly-source RefSeq \
            --filename "${zip_file}"

        unzip -d ${fastANI_work_dir} "${zip_file}"

        datasets rehydrate --directory ${fastANI_work_dir}

        echo " $(date) -- Creating taxonomy file: '${taxonomy_file}'..."

        # Create the file and add header
        echo -e "genome_id\tspecies_name" \
            > "${taxonomy_file}"

        assembly_data_report="${seq_data_dir}/assembly_data_report.jsonl"
        # Select assembly IDs
        cat "${assembly_data_report}" \
            | grep -Eo '"accession":"GCF_[0-9\.]+"' \
            | grep -Eo 'GCF_[0-9\.]+' \
            > "${TMPDIR}/genome_ids.txt"
        # Select species names
        cat "${assembly_data_report}" \
            | grep -Eo '"submittedSpecies":"[^"]+"' \
            | sed 's/"submittedSpecies":"//' \
            | sed 's/"//' \
            > "${TMPDIR}/species_names.txt"
        # Combine assembly IDs and species names into single TSV file
        #   of two rows
        paste -d '\t' \
            "${TMPDIR}/genome_ids.txt" "${TMPDIR}/species_names.txt" \
            >> "${taxonomy_file}"
        rm -v "${TMPDIR}/genome_ids.txt" "${TMPDIR}/species_names.txt"
        echo -e "query\tNA NA" \
            >> "${taxonomy_file}"
        echo ' Taxonomy file is created.'
    fi

    # Create input sequence list file
    cat "${QUERY_LIST_FILE}" > "${fastANI_input_list_file}"
    if [[ -d "${seq_data_dir}" ]]; then
        find "${seq_data_dir}" -type f -name *.fna >> "${fastANI_input_list_file}"
    fi

    echo "----------------------------------------------------------------------"
    echo -e "${COLOR} $(date) Starting FastANI - calculate ANI value from genomic sequences ${RESET_COLOR}"
    echo "----------------------------------------------------------------------"

    "${FASTANI}" \
        --refList "${fastANI_input_list_file}" \
        --queryList "${fastANI_input_list_file}" \
        --threads "${THREADS_NUM}" \
        --output "${raw_result_tsv}"

    echo "----------------------------------------------------------------------"
    echo -e "${COLOR} $(date) Writing final output file ${RESET_COLOR}"
    echo "----------------------------------------------------------------------"
    header_col_names=(
        "query_genome\t"
        "reference_genome\t"
        "ANI\t"
        "bidirectional_frag_map_num\t"
        "total_query_fragments\n"
    )
    echo -n '' > "${result_tsv}"
    for col_name in "${header_col_names[@]}"; do
        echo -en "${col_name}"
    done >> "${result_tsv}"
    while read -r query_genome reference_genome col3 col4 col5; do
        query_basename="$(basename ${query_genome})"
        query_label="${query_basename%.fna}"
        ref_basename="$(basename ${reference_genome})"
        ref_label="${ref_basename%.fna}"
        echo -en "${query_label}\t${ref_label}\t${col3}\t${col4}\t${col5}\n"
    done < "${raw_result_tsv}" >> "${result_tsv}"

    echo "created '${result_tsv}'"
    rm -v "${raw_result_tsv}"

} |& tee "${log_file}"

echo ""
echo "---------------------------------------------------------------------------------"
echo -e "${COLOR} $(date) -- Completed! ${RESET_COLOR}"
echo " Result table: '${result_tsv}'"
echo " Taxonomy file: '${taxonomy_file}'"
echo " Have fun and please come again!"

exit 0
