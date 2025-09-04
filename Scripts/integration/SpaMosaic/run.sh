#!/bin/bash

# SpaMosaic Integration Execution Script
# This script provides easy access to run SpaMosaic integration tasks
# Compatible with SMOBench project structure

set -e  # Exit on any error

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/run_spamosaic.py"

# Default environment
CONDA_ENV="smobench"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Help function
show_help() {
    cat << EOF
SpaMosaic Integration Runner

USAGE:
    $0 [OPTIONS] DATASET TASK MODALITIES SAMPLES...

ARGUMENTS:
    DATASET     Dataset name (e.g., HLN, MISAR, Mouse_Brain)
    TASK        Integration task: vertical, horizontal, or mosaic
    MODALITIES  Space-separated modalities (rna, adt, atac)
    SAMPLES     Space-separated sample names

OPTIONS:
    -d, --data-dir DIR      Base data directory (default: \$PROJECT_ROOT/Dataset/withGT)
    -o, --output-dir DIR    Output directory (default: \$PROJECT_ROOT/Results/SpaMosaic/DATASET)
    -g, --gpu ID            GPU device ID (default: 0)
    -e, --env NAME          Conda environment name (default: smobench)
    -c, --clusters NUM      Number of clusters (default: auto-detect)
    --epochs NUM            Training epochs (default: 100)
    --lr FLOAT              Learning rate (default: 0.01)
    --seed NUM              Random seed (default: 2024)
    --dry-run               Show command without executing
    -h, --help              Show this help message

EXAMPLES:
    # Vertical integration (RNA+ADT) for HLN dataset
    $0 HLN vertical "rna adt" A1 D1

    # Horizontal integration (RNA only) for Mouse_Brain
    $0 Mouse_Brain horizontal rna sample1 sample2 sample3

    # Mosaic integration (RNA+ATAC) for MISAR
    $0 MISAR mosaic "rna atac" E11 E13 E15 E18

    # Custom data directory and output
    $0 -d /path/to/data -o /path/to/output HLN vertical "rna adt" A1

    # Specify GPU and clusters
    $0 -g 1 -c 10 HLN vertical "rna adt" A1 D1

DATASET STRUCTURE:
    Data should be organized as:
    DATA_DIR/DATASET/SAMPLE/adata_RNA.h5ad
    DATA_DIR/DATASET/SAMPLE/adata_ADT.h5ad
    DATA_DIR/DATASET/SAMPLE/adata_ATAC.h5ad

EOF
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    # Check conda
    if ! command -v conda &> /dev/null; then
        log_error "Conda not found. Please install Miniconda or Anaconda."
        exit 1
    fi
    
    # Check conda environment
    if ! conda env list | grep -q "^${CONDA_ENV}\s"; then
        log_error "Conda environment '${CONDA_ENV}' not found."
        log_info "Please create the environment with: conda create -n ${CONDA_ENV} python=3.9"
        exit 1
    fi
    
    # Check Python script
    if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
        log_error "Python script not found: ${PYTHON_SCRIPT}"
        exit 1
    fi
    
    log_success "All dependencies checked"
}

# Validate arguments
validate_arguments() {
    local dataset="$1"
    local task="$2"
    local modalities="$3"
    shift 3
    local samples=("$@")
    
    # Check required arguments
    if [[ -z "$dataset" || -z "$task" || -z "$modalities" || ${#samples[@]} -eq 0 ]]; then
        log_error "Missing required arguments"
        show_help
        exit 1
    fi
    
    # Validate task
    case "$task" in
        vertical|horizontal|mosaic)
            ;;
        *)
            log_error "Invalid task: $task. Must be vertical, horizontal, or mosaic"
            exit 1
            ;;
    esac
    
    # Validate modalities
    for mod in $modalities; do
        case "$mod" in
            rna|adt|atac)
                ;;
            *)
                log_error "Invalid modality: $mod. Must be rna, adt, or atac"
                exit 1
                ;;
        esac
    done
    
    log_success "Arguments validated"
}

# Check data availability
check_data() {
    local data_dir="$1"
    local dataset="$2"
    local modalities="$3"
    shift 3
    local samples=("$@")
    
    log_info "Checking data availability..."
    
    local dataset_dir="${data_dir}/${dataset}"
    if [[ ! -d "$dataset_dir" ]]; then
        log_error "Dataset directory not found: $dataset_dir"
        exit 1
    fi
    
    local missing_files=()
    for sample in "${samples[@]}"; do
        local sample_dir="${dataset_dir}/${sample}"
        if [[ ! -d "$sample_dir" ]]; then
            log_warning "Sample directory not found: $sample_dir"
            continue
        fi
        
        for mod in $modalities; do
            local mod_upper=$(echo "$mod" | tr '[:lower:]' '[:upper:]')
            local file_path="${sample_dir}/adata_${mod_upper}.h5ad"
            if [[ ! -f "$file_path" ]]; then
                missing_files+=("$file_path")
            fi
        done
    done
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        log_warning "Some data files are missing:"
        for file in "${missing_files[@]}"; do
            log_warning "  - $file"
        done
        log_warning "Integration will continue with available data"
    else
        log_success "All data files found"
    fi
}

# Parse command line arguments
parse_arguments() {
    # Default values
    DATA_DIR="${PROJECT_ROOT}/Dataset/withGT"
    OUTPUT_DIR=""
    GPU_ID=0
    N_CLUSTERS=""
    EPOCHS=100
    LEARNING_RATE=0.01
    SEED=2024
    DRY_RUN=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--data-dir)
                DATA_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -g|--gpu)
                GPU_ID="$2"
                shift 2
                ;;
            -e|--env)
                CONDA_ENV="$2"
                shift 2
                ;;
            -c|--clusters)
                N_CLUSTERS="$2"
                shift 2
                ;;
            --epochs)
                EPOCHS="$2"
                shift 2
                ;;
            --lr)
                LEARNING_RATE="$2"
                shift 2
                ;;
            --seed)
                SEED="$2"
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
            *)
                break
                ;;
        esac
    done
    
    # Remaining arguments
    DATASET="$1"
    TASK="$2"
    MODALITIES="$3"
    shift 3
    SAMPLES=("$@")
    
    # Set default output directory if not specified
    if [[ -z "$OUTPUT_DIR" ]]; then
        OUTPUT_DIR="${PROJECT_ROOT}/Results/adata/SpaMosaic/${DATASET}"
    fi
}

# Build command
build_command() {
    local cmd="conda run -n ${CONDA_ENV} python ${PYTHON_SCRIPT}"
    
    # Required arguments
    cmd+=" --dataset '${DATASET}'"
    cmd+=" --task '${TASK}'"
    cmd+=" --modalities ${MODALITIES}"
    cmd+=" --samples ${SAMPLES[*]}"
    cmd+=" --data_dir '${DATA_DIR}'"
    cmd+=" --output_dir '${OUTPUT_DIR}'"
    
    # Optional arguments
    cmd+=" --gpu ${GPU_ID}"
    cmd+=" --epochs ${EPOCHS}"
    cmd+=" --lr ${LEARNING_RATE}"
    cmd+=" --seed ${SEED}"
    
    if [[ -n "$N_CLUSTERS" ]]; then
        cmd+=" --n_clusters ${N_CLUSTERS}"
    fi
    
    echo "$cmd"
}

# Main execution
main() {
    log_info "Starting SpaMosaic Integration Runner"
    log_info "Project root: ${PROJECT_ROOT}"
    
    # Parse arguments
    parse_arguments "$@"
    
    # Validate inputs
    validate_arguments "$DATASET" "$TASK" "$MODALITIES" "${SAMPLES[@]}"
    
    # Check dependencies
    check_dependencies
    
    # Check data
    check_data "$DATA_DIR" "$DATASET" "$MODALITIES" "${SAMPLES[@]}"
    
    # Build command
    local cmd
    cmd=$(build_command)
    
    # Display execution plan
    log_info "Execution Plan:"
    log_info "  Dataset: ${DATASET}"
    log_info "  Task: ${TASK}"
    log_info "  Modalities: ${MODALITIES}"
    log_info "  Samples: ${SAMPLES[*]}"
    log_info "  Data Directory: ${DATA_DIR}"
    log_info "  Output Directory: ${OUTPUT_DIR}"
    log_info "  GPU: ${GPU_ID}"
    log_info "  Environment: ${CONDA_ENV}"
    log_info "  Command: ${cmd}"
    
    # Execute or dry run
    if [[ "$DRY_RUN" == true ]]; then
        log_info "Dry run mode - command would be:"
        echo "$cmd"
        exit 0
    fi
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Execute command
    log_info "Executing SpaMosaic integration..."
    if eval "$cmd"; then
        log_success "SpaMosaic integration completed successfully!"
        log_success "Results saved to: ${OUTPUT_DIR}"
    else
        log_error "SpaMosaic integration failed!"
        exit 1
    fi
}

# Run main function with all arguments
main "$@"