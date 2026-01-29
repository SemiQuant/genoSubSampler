#!/bin/bash
#
# genoSubSampler Installation Script
# Author: SemiQuant (JasonLimberis@ucsf.edu)
#

set -e

echo "=============================================="
echo "   genoSubSampler Dependency Installer"
echo "=============================================="
echo ""

# Detect package manager
detect_package_manager() {
    if command -v micromamba &> /dev/null; then
        echo "micromamba"
    elif command -v mamba &> /dev/null; then
        echo "mamba"
    elif command -v conda &> /dev/null; then
        echo "conda"
    else
        echo "none"
    fi
}

PKG_MANAGER=$(detect_package_manager)
INSTALL_MTB_FILTER=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --with-filter)
            INSTALL_MTB_FILTER=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --with-filter    Install additional dependencies for M.tb filtering"
            echo "                   (bowtie2, samtools, bedtools)"
            echo "  -h, --help       Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "Detected package manager: $PKG_MANAGER"
echo ""

# Core dependencies
CORE_DEPS="seqtk p7zip bc"

# Additional dependencies for M.tb filtering
FILTER_DEPS="bowtie2 samtools bedtools htslib"

install_with_conda() {
    local manager="$1"
    shift
    local deps="$@"
    
    echo "Installing with $manager: $deps"
    $manager install -y -c bioconda -c conda-forge $deps
}

install_with_brew() {
    local deps="$@"
    echo "Installing with Homebrew: $deps"
    
    for dep in $deps; do
        case "$dep" in
            p7zip)
                brew install p7zip
                ;;
            seqtk)
                brew install seqtk
                ;;
            bc)
                # bc is usually pre-installed on macOS, but install if needed
                brew install bc 2>/dev/null || true
                ;;
            bowtie2)
                brew install bowtie2
                ;;
            samtools)
                brew install samtools
                ;;
            bedtools)
                brew install bedtools
                ;;
            htslib)
                brew install htslib
                ;;
            *)
                echo "Unknown package: $dep"
                ;;
        esac
    done
}

install_with_apt() {
    local deps="$@"
    echo "Installing with apt: $deps"
    
    sudo apt-get update
    for dep in $deps; do
        case "$dep" in
            p7zip)
                sudo apt-get install -y p7zip-full
                ;;
            seqtk)
                sudo apt-get install -y seqtk
                ;;
            bc)
                sudo apt-get install -y bc
                ;;
            bowtie2)
                sudo apt-get install -y bowtie2
                ;;
            samtools)
                sudo apt-get install -y samtools
                ;;
            bedtools)
                sudo apt-get install -y bedtools
                ;;
            htslib)
                sudo apt-get install -y tabix
                ;;
            *)
                echo "Unknown package: $dep"
                ;;
        esac
    done
}

# Main installation logic
case "$PKG_MANAGER" in
    micromamba|mamba|conda)
        echo "Using $PKG_MANAGER for installation..."
        install_with_conda "$PKG_MANAGER" $CORE_DEPS
        
        if [[ "$INSTALL_MTB_FILTER" == true ]]; then
            echo ""
            echo "Installing M.tb filtering dependencies..."
            install_with_conda "$PKG_MANAGER" $FILTER_DEPS
        fi
        ;;
    none)
        echo "No conda/mamba/micromamba found."
        echo ""
        
        if command -v brew &> /dev/null; then
            echo "Using Homebrew for installation..."
            install_with_brew $CORE_DEPS
            
            if [[ "$INSTALL_MTB_FILTER" == true ]]; then
                echo ""
                echo "Installing M.tb filtering dependencies..."
                install_with_brew $FILTER_DEPS
            fi
        elif command -v apt-get &> /dev/null; then
            echo "Using apt for installation..."
            install_with_apt $CORE_DEPS
            
            if [[ "$INSTALL_MTB_FILTER" == true ]]; then
                echo ""
                echo "Installing M.tb filtering dependencies..."
                install_with_apt $FILTER_DEPS
            fi
        else
            echo "ERROR: No supported package manager found."
            echo ""
            echo "Please install the following manually:"
            echo "  Core: seqtk, 7zip (p7zip), bc"
            if [[ "$INSTALL_MTB_FILTER" == true ]]; then
                echo "  Filtering: bowtie2, samtools, bedtools, bgzip (htslib)"
            fi
            exit 1
        fi
        ;;
esac

echo ""
echo "=============================================="
echo "   Verifying Installation"
echo "=============================================="

check_command() {
    if command -v "$1" &> /dev/null; then
        echo "  ✓ $1"
        return 0
    else
        echo "  ✗ $1 (NOT FOUND)"
        return 1
    fi
}

echo ""
echo "Core dependencies:"
check_command seqtk
check_command 7z
check_command bc

if [[ "$INSTALL_MTB_FILTER" == true ]]; then
    echo ""
    echo "Filtering dependencies:"
    check_command bowtie2
    check_command samtools
    check_command bedtools
    check_command bgzip
fi

# Build bowtie2 index if filtering is enabled
if [[ "$INSTALL_MTB_FILTER" == true ]]; then
    echo ""
    echo "=============================================="
    echo "   Building Bowtie2 Index"
    echo "=============================================="
    
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REF_FASTA="${SCRIPT_DIR}/Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta"
    REF_INDEX="${SCRIPT_DIR}/Mycobacterium_tuberculosis_H37Rv_genome_v4"
    
    if [[ -f "$REF_FASTA" ]]; then
        # Check if index already exists
        if [[ -f "${REF_INDEX}.1.bt2" ]]; then
            echo "Bowtie2 index already exists, skipping build."
        else
            echo "Building bowtie2 index from: $REF_FASTA"
            echo "This may take a few minutes..."
            bowtie2-build "$REF_FASTA" "$REF_INDEX"
            echo "Bowtie2 index built successfully."
        fi
    else
        echo "WARNING: Reference FASTA not found at: $REF_FASTA"
        echo "You will need to provide your own reference with --ref"
    fi
fi

echo ""
echo "=============================================="
echo "   Installation Complete!"
echo "=============================================="
echo ""
echo "You can now run: ./subsample.sh -h"
echo ""
