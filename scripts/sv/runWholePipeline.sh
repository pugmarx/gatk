#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [5] absolute path to the reference fasta on the cluster (2bit is assumed accompanying with same basename and extension \".2bit\", skiplist with extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [6] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "  [*] extra command-line arguments to StructuralVariationDiscoveryPipelineSpark"
    echo -e "Example:"
    echo -e " bash svDiscover.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-test-cluster \\"
    echo -e "      /test-sample \\"
    echo -e "      /data/NA12878_test.bam \\"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta"
    echo -e "      /mnt/1/reference/Homo_sapiens_assembly38.fasta.img (the /mnt/1/ prefix is the default mounting location of the ssd if the cluster is created with a local ssd for the worker nodes)"
    exit 1
fi


GATK_DIR="$1"
CLUSTER_NAME="$2"
MASTER_NODE="hdfs://${CLUSTER_NAME}-m:8020"
PROJECT_OUTPUT_DIR="${MASTER_NODE}$3"
INPUT_BAM="${MASTER_NODE}$4"
REF_FASTA="${MASTER_NODE}$5"
REF_INDEX_IMAGE="$6"
INTERVAL_KILL_LIST=$(echo "${REF_FASTA}" | sed 's/.fasta$/.kill.intervals/')
KMER_KILL_LIST=$(echo "${REF_FASTA}" | sed 's/.fasta$/.kill.kmers/')
ALTS_KILL_LIST=$(echo "${REF_FASTA}" | sed 's/.fasta$/.kill.alts/')
REF_TWOBIT=$(echo "${REF_FASTA}" | sed 's/.fasta$/.2bit/')

# extract any extra arguments to StructuralVariationDiscoveryPipelineSpark
shift $(($# < 6 ? $# : 6))
SV_ARGS=${*:-${SV_ARGS:-""}}
# expand any local variables passed as strings (e.g. PROJECT_OUTPUT_DIR)
eval "SV_ARGS=\"${SV_ARGS}\""

# Choose 3 executors per worker, one on master
# NOTE: this would find preemptible workers, but it produces
# (erroneous?) deprecation warnings
#NUM_WORKERS=$(gcloud compute instances list --filter="name ~ ${CLUSTER_NAME}-[sw].*" | grep RUNNING | wc -l)
# this works but does not see preemptible workers
#NUM_WORKERS=$(gcloud dataproc clusters list --filter "clusterName = ${CLUSTER_NAME}" | tail -n 1 | awk '{print $2}')
# This gets all. the -1 is to subtract away master
NUM_WORKERS=$(gcloud dataproc clusters describe ${CLUSTER_NAME} | grep numInstances | awk '{sum += $2} END {print sum - 1}')
if [ -z "${NUM_WORKERS}" ]; then
    echo "Cluster \"${CLUSTER_NAME}\" not found"
    exit 1
fi
# was: 3 * NUM_WORKERS + 1
#      --driver-memory 31G
#      --executor-memory 31G
#      --executor-cores 5
#      --conf spark.yarn.executor.memoryOverhead=3300
# tried: 5 * NUM_WORKERS + 2
#      --driver-memory 19G
#      --executor-memory 19G
#      --executor-cores 3
#      --conf spark.yarn.executor.memoryOverhead=1600
#  got a crash the time I tried it...
# next time, try first way, but with 7 executor cores instead of 5
NUM_EXECUTORS=$((3 * ${NUM_WORKERS} + 1))

"${GATK_DIR}/gatk-launch" StructuralVariationDiscoveryPipelineSpark \
    -I "${INPUT_BAM}" \
    -O "${PROJECT_OUTPUT_DIR}/variants/inv_del_ins.vcf" \
    -R "${REF_TWOBIT}" \
    --fastaReference "${REF_FASTA}" \
    --alignerIndexImage "${REF_INDEX_IMAGE}" \
    --exclusionIntervals "${INTERVAL_KILL_LIST}" \
    --kmersToIgnore "${KMER_KILL_LIST}" \
    --crossContigsToIgnore "${ALTS_KILL_LIST}" \
    --breakpointIntervals "${PROJECT_OUTPUT_DIR}/intervals" \
    --fastqDir "${PROJECT_OUTPUT_DIR}/fastq" \
    --contigSAMFile "${PROJECT_OUTPUT_DIR}/assemblies.sam" \
    --targetLinkFile ${PROJECT_OUTPUT_DIR}/target_links.bedpe \
    ${SV_ARGS} \
    -- \
    --sparkRunner GCS \
    --cluster "${CLUSTER_NAME}" \
    --num-executors ${NUM_EXECUTORS} \
    --driver-memory 31G \
    --executor-memory 31G \
    --executor-cores 7 \
    --conf spark.yarn.executor.memoryOverhead=3300 \
    --conf spark.network.timeout=600 \
    --conf spark.executor.heartbeatInterval=120
