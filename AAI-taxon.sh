#!/bin/bash

set -eu
# Параметры: $1 - taxon ID or name
#            $2 - output directory
#            $3 - Input multi-fasta file: genome (.fna) or protein sequences (.faa)
#            $4 - Mode: 'genome' or 'protein'
#            $5 - Path to EzAAI-1.2.3_masikol.0.2-jar-with-dependencies.jar
#            $6 - Number of threads to use. Default: 1
#            $7 - temporary directory for EzAAI. Default: /tmp/ezaai
#            Dependencies:
#               ncbi datasets
#               ezaai
#               A jar file of ezaai version 1.2.3_masikol.0.2: EzAAI-1.2.3_masikol.0.2-jar-with-dependencies.jar

VERSION='1.0.b'
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
echo " - presents - the program AAI-taxon! Version ${VERSION}"
echo " - `date`"
echo " ------------------------------------------------------------------------"

if [[ -n "$1" ]]; then
    TAXON="${1}"
    echo " "
    echo -e " ${COLOR} Welcome to our crimson submarine! ${RESET_COLOR}"
    echo -e " ${COLOR} Please take a deep breath: we are about to begin! ${RESET_COLOR}"
    echo "----------------------------------------------------"
    echo " "
    echo " You have chosen Taxon ID / Name: '${TAXON}'"
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
    QUERY_FILE="${3}"
    echo " Query fna/faa file: '${QUERY_FILE}'"
else
    echo " Error: query sequence file in fna or faa format (\$3) is not specified!"
    exit 1
fi

if [[ -n "$4" ]]; then
    MODE="${4}"
    echo " Mode: '${MODE}'"
    allowed_modes=( 'genome' 'protein' )
    if [[ ! "${allowed_modes[@]}" =~ "${MODE}" ]]; then
        echo " Invalid mode: '${MODE}'. Allowed modes: ${allowed_modes[@]}"
        exit 1
    fi
    if [[ "${MODE}" == 'genome' ]]; then
        QUERY_EXTENSION='fna'
    elif [[ "${MODE}" == 'protein' ]]; then
        QUERY_EXTENSION='faa'
    fi

    if [[ "${QUERY_FILE}" != *.${QUERY_EXTENSION} ]]; then
        echo " Error: query file extension is inappropriate: '${QUERY_FILE}'"
        echo " The extension must be '.${QUERY_EXTENSION}'. It is indeed important."
        echo ' Please provide the correct file.'
        echo ' Or, if you are sure that the file is correct, just rename it.'
        exit 1
    fi
else
    echo " Error: sequence mode (genome or protein, \$4) is not specified!"
    exit 1
fi

if [[ -n "$5" ]]; then
    EZAAI_MASIKOL="${5}"
else
    echo "Error: path to EzAAI-1.2.3_masikol.0.2-jar-with-dependencies.jar (\$5) is not specified!"
    echo 'Please download it here: https://github.com/masikol/ezaai'
    exit 1
fi

THREADS_NUM=1
if [[ -n "$6" ]]; then
    THREADS_NUM="${6}"
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

TMPDIR='/tmp/ezaai'
if [[ -n "$7" ]]; then
    TMPDIR="$7"
fi
if [[ ! -d "${TMPDIR}" ]]; then
    mkdir -pv "${TMPDIR}"
fi
if [[ ! -d "${TMPDIR}" ]]; then
    echo "Error: temporary directory does not exist and we cannot craete it: '${TMPDIR}'"
    exit 1
fi

echo 'Run parameters:'
echo " - Taxon ID / Name: '${TAXON}'"
echo " - Output directory: '${WORKDIR_ROOT}'"
echo " - Query file path: '${QUERY_FILE}'"
echo " - Sequence mode: '${MODE}'"
echo " - EzAAI masikol version path: '${EZAAI_MASIKOL}'"
echo " - Number of CPU threads to use: '${THREADS_NUM}'"
echo " - Temporary directory: '${TMPDIR}'"
echo '--------------------------'



# Check dependencies
dependencies=( datasets ezaai )
for exe_file_name in "${dependencies[@]}"; do
    if [[ -z $(which "${exe_file_name}" ) ]]; then
        echo "Error: cannot find program '${exe_file_name}'"
        exit 1
    fi
done


# Make useful variables
query_basename=$(basename "${QUERY_FILE}")
query_name="${query_basename}%${QUERY_EXTENSION}"
query_label="${query_basename}"
datadir="${WORKDIR_ROOT}/${TAXON}_db"
ezaai_db_dir="${datadir}/ezaai-db"
log_file="${WORKDIR_ROOT}/result.log"
seq_data_dir="${datadir}/ncbi_dataset/data"
taxonomy_file="${WORKDIR_ROOT}/taxonomy.tsv"
result_tsv="${WORKDIR_ROOT}/aai_result.tsv"
result_tree="${WORKDIR_ROOT}/aai_result.nwk"

if [[ "${MODE}" == 'genome' ]]; then
    seq_mode_plural='genomes'
elif [[ "${MODE}" == 'protein' ]]; then
    seq_mode_plural='protein sequences'
fi


# Create necessary dirs and check paths
echo " "
if [[ -d "${datadir}" ]]; then
    echo "Directory ${datadir} уже существует"
else
    mkdir -pv "${datadir}"
fi
if [ -d "${ezaai_db_dir}" ]; then
    echo "Directory ${ezaai_db_dir} уже существует"
else
    mkdir -pv "${datadir}/ezaai-db"
fi

# Create or empty log file
echo -n '' > "${log_file}"

{
    echo "----------------------------------------------"
    echo -e " ${COLOR} Downloading ${seq_mode_plural} of type strains of ${TAXON} ${RESET_COLOR}"
    echo "----------------------------------------------"

    if [[ "${MODE}" == 'genome' ]]; then
        include_arg='genome'
    elif [[ "${MODE}" == 'protein' ]]; then
        include_arg='protein'
    fi

    zip_file="${WORKDIR_ROOT}/${TAXON}_RefSeq_${MODE}.zip"

    datasets download genome taxon "${TAXON}" \
        --include "${include_arg}" \
        --from-type \
        --dehydrated \
        --assembly-source RefSeq \
        --filename "${zip_file}"

    unzip -d ${datadir} "${zip_file}"

    datasets rehydrate --directory ${datadir}

    # Rename 'protein.faa' files to [ASSEMBLY_ID].faa
    if [[ "${MODE}" == 'protein' ]]; then
        for prot_dir in ${seq_data_dir}/GCF_*; do
            dir_basename=$(basename "${prot_dir}")
            prot_file="${prot_dir}/protein.faa"
            if [[ -f "${prot_file}" ]]; then
                mv -v "${prot_dir}/protein.faa" "${prot_dir}/${dir_basename}.faa"
            fi
        done
    fi

    echo " Creating taxonomy file: '${taxonomy_file}'..."

    # Craete the file and add header
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


    echo "-------------------------------------------------------------------"
    echo -e " ${COLOR} Starting EzAAI - make profile DB from sequences ${RESET_COLOR}"
    echo "-------------------------------------------------------------------"

    # Check if seq_data_dir exists
    if [ ! -d "${seq_data_dir}" ]; then
        echo "Error: sequence data directory does not exist: '${seq_data_dir}'"
        exit 1
    fi

    if [[ "${MODE}" == 'genome' ]]; then
        ezaai_subprogram='extract'
        ezaai_convert_s_option=''
    elif [[ "${MODE}" == 'protein' ]]; then
        ezaai_subprogram='convert'
        ezaai_convert_s_option='-s prot'
    fi

    # Проходим по всем подкаталогам в указанной директории
    find "${seq_data_dir}" -type f -name "*.${QUERY_EXTENSION}" | while read -r seq_file; do

        # Извлекаем имя файла без расширения
        file_dir_name=$(dirname "${seq_file}")
        file_basename=$(basename "${seq_file}")
        genome_id=$(basename "${file_dir_name}")

        db_file="${ezaai_db_dir}/${genome_id}.db"

        echo "----------------------------"
        echo "Путь к файлу: '${seq_file}'"
        echo "Genome ID: '${genome_id}'"

        # Example for genome:
        # ezaai extract \
        #     -i $file_path \
        #     -o ${datadir}/ezaai-db/$file_name.db \
        #     -l "${genome_id}" 
        # Example for protein:
        # ezaai convert \
        #     -i $seq_file \
        #     -s prot \
        #     -o "${ezaai_db_dir}/${file_name}.db" \
        #     -l "${genome_id}"

        ezaai "${ezaai_subprogram}" ${ezaai_convert_s_option} \
            -i "${seq_file}" \
            -o "${db_file}" \
            -l "${genome_id}"

    done
    echo "----------------------------"

    # Agg query genome to the database
    genome_id='query'
    file_basename=$(basename "${QUERY_FILE}")
    db_file="${ezaai_db_dir}/${genome_id}.db"
    ezaai "${ezaai_subprogram}" ${ezaai_convert_s_option} \
        -i "${QUERY_FILE}" \
        -o "${db_file}" \
        -l "${genome_id}"

    echo "----------------------------------------------------------------------"
    echo -e " ${COLOR} Starting EzAAI - calculate AAI value from profile DBs ${RESET_COLOR}"
    echo "----------------------------------------------------------------------"

    java -jar "${EZAAI_MASIKOL}" \
        calculate \
        -i "${ezaai_db_dir}" \
        -j "${ezaai_db_dir}" \
        -t "${THREADS_NUM}" \
        -self 1 \
        -tmp "${TMPDIR}" \
        -o "${result_tsv}"

    echo "---------------------------------------------------------------------------------"
    echo -e " ${COLOR} Starting EzAAI - hierarchical clustering of taxa with AAI values ${RESET_COLOR}"
    echo "---------------------------------------------------------------------------------"

    ezaai cluster \
        -i "${result_tsv}" \
        -o "${result_tree}"

} |& tee "${log_file}"

echo ""
echo "---------------------------------------------------------------------------------"
echo -e " ${COLOR} Completed! ${RESET_COLOR}"
echo " Result table: '${result_tsv}'"
echo " Result tree: '${result_tree}'"
echo " Taxonomy file: '${taxonomy_file}'"
echo " Have fun and please come again!"

exit 0
