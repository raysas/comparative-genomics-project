#!/bin/bash
################################################################################
# Parallel Processing Helper Functions
# 
# Utility functions for parallel execution across both pipelines
# 
# CG Project - Refactor Branch
# Date: 2025-11-23
################################################################################

# Detect available CPU cores
detect_cores() {
    local cores=1
    if command -v nproc &> /dev/null; then
        cores=$(nproc)
    elif [ -f /proc/cpuinfo ]; then
        cores=$(grep -c ^processor /proc/cpuinfo)
    elif command -v sysctl &> /dev/null; then
        cores=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
    fi
    echo "$cores"
}

# Detect available memory (in GB)
detect_memory() {
    local mem_gb=4  # default
    if [ -f /proc/meminfo ]; then
        local mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        mem_gb=$((mem_kb / 1024 / 1024))
    elif command -v sysctl &> /dev/null; then
        local mem_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 4294967296)
        mem_gb=$((mem_bytes / 1024 / 1024 / 1024))
    fi
    echo "$mem_gb"
}

# Calculate optimal job count based on available resources
calculate_optimal_jobs() {
    local task_type="${1:-default}"  # blast, align, ks, default
    local cores=$(detect_cores)
    local mem_gb=$(detect_memory)
    local optimal_jobs=$cores
    
    case "$task_type" in
        blast)
            # BLAST is CPU and memory intensive
            # Use 75% of cores, ensure at least 2GB per job
            optimal_jobs=$((cores * 3 / 4))
            local max_by_memory=$((mem_gb / 2))
            [ $optimal_jobs -gt $max_by_memory ] && optimal_jobs=$max_by_memory
            ;;
        align)
            # Alignment (MUSCLE/MAFFT) is moderately intensive
            # Use 80% of cores, ensure at least 1GB per job
            optimal_jobs=$((cores * 4 / 5))
            local max_by_memory=$mem_gb
            [ $optimal_jobs -gt $max_by_memory ] && optimal_jobs=$max_by_memory
            ;;
        ks)
            # Ks calculation (yn00) is lightweight
            # Can use all cores
            optimal_jobs=$cores
            ;;
        *)
            # Conservative default: 75% of cores
            optimal_jobs=$((cores * 3 / 4))
            ;;
    esac
    
    # Ensure at least 1 job
    [ $optimal_jobs -lt 1 ] && optimal_jobs=1
    
    echo "$optimal_jobs"
}

# Print system resource summary
print_resource_summary() {
    local cores=$(detect_cores)
    local mem_gb=$(detect_memory)
    
    echo "System Resources:"
    echo "  CPU Cores: $cores"
    echo "  Memory: ${mem_gb} GB"
    echo ""
    echo "Recommended Parallel Jobs:"
    echo "  BLAST:     $(calculate_optimal_jobs blast) jobs"
    echo "  Alignment: $(calculate_optimal_jobs align) jobs"
    echo "  Ks Calc:   $(calculate_optimal_jobs ks) jobs"
}

# Check if GNU parallel is available
check_parallel_available() {
    if command -v parallel &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Run command with GNU parallel if available, otherwise sequentially
run_parallel_or_sequential() {
    local jobs="$1"
    local input_list="$2"
    local command_template="$3"
    
    if check_parallel_available; then
        echo "Running with GNU parallel ($jobs jobs)..."
        cat "$input_list" | parallel -j "$jobs" "$command_template"
    else
        echo "GNU parallel not found - running sequentially..."
        while IFS= read -r item; do
            eval "$command_template" "$item"
        done < "$input_list"
    fi
}

# Split FASTA file into chunks for parallel processing
split_fasta_for_parallel() {
    local input_fasta="$1"
    local num_chunks="$2"
    local output_dir="$3"
    
    mkdir -p "$output_dir"
    
    # Count total sequences
    local total_seqs=$(grep -c '^>' "$input_fasta")
    local seqs_per_chunk=$(( (total_seqs + num_chunks - 1) / num_chunks ))
    
    echo "Splitting $total_seqs sequences into $num_chunks chunks..."
    echo "Approximately $seqs_per_chunk sequences per chunk"
    
    # Use awk to split by sequence count
    awk -v chunk_size="$seqs_per_chunk" -v outdir="$output_dir" '
        BEGIN { n=0; chunk=1; file=outdir"/chunk_"chunk".fa" }
        /^>/ { 
            if (n >= chunk_size && n > 0) { 
                chunk++; 
                n=0; 
                file=outdir"/chunk_"chunk".fa" 
            }
            n++
        }
        { print > file }
    ' "$input_fasta"
    
    local actual_chunks=$(ls -1 "$output_dir"/chunk_*.fa 2>/dev/null | wc -l)
    echo "Created $actual_chunks chunk files in $output_dir"
}

# Merge results from parallel processing
merge_parallel_results() {
    local input_pattern="$1"
    local output_file="$2"
    local has_header="${3:-no}"
    
    echo "Merging results from: $input_pattern"
    
    if [ "$has_header" = "yes" ]; then
        # Take header from first file, then all data
        local first_file=$(ls $input_pattern 2>/dev/null | head -n1)
        if [ -n "$first_file" ]; then
            head -n1 "$first_file" > "$output_file"
            for file in $input_pattern; do
                tail -n+2 "$file" >> "$output_file"
            done
        fi
    else
        # Simple concatenation
        cat $input_pattern > "$output_file" 2>/dev/null || true
    fi
    
    echo "Merged to: $output_file"
}

# Monitor parallel job progress
monitor_parallel_jobs() {
    local total_jobs="$1"
    local output_dir="$2"
    local pattern="${3:-.done}"
    
    local completed=0
    local last_completed=0
    
    echo "Monitoring progress in: $output_dir"
    
    while [ $completed -lt $total_jobs ]; do
        completed=$(find "$output_dir" -name "*$pattern" 2>/dev/null | wc -l)
        
        if [ $completed -ne $last_completed ]; then
            local percent=$((completed * 100 / total_jobs))
            echo "Progress: $completed/$total_jobs ($percent%)"
            last_completed=$completed
        fi
        
        sleep 5
    done
    
    echo "All jobs completed!"
}

# Export functions for GNU parallel
export -f detect_cores
export -f detect_memory
export -f calculate_optimal_jobs
export -f check_parallel_available

################################################################################
# Example Usage
################################################################################

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    # Script is being run directly
    echo "Parallel Processing Helper Functions"
    echo "======================================"
    echo ""
    print_resource_summary
    echo ""
    
    if check_parallel_available; then
        echo "✓ GNU parallel is installed"
        parallel --version | head -n1
    else
        echo "✗ GNU parallel not found"
        echo "  Install with: sudo apt-get install parallel"
        echo "  Or: conda install -c conda-forge parallel"
    fi
fi
