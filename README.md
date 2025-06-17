# hpc_genomics_assignment

# Computational Assignment:  Integrated Human Genome Analysis Pipeline

--

## **Overview**

This comprehensive assignment integrates three fundamental bioinformatics skill sets you practiced in this module into a unified pipeline for analyzing human genomic data. You will work with high-performance computing environments, process large-scale genomic datasets, extract genomic features, and perform comparative sequence alignment analyses. This assignment mirrors real-world computational genomics workflows used in medical research and precision medicine applications.

## **Learning Objectives**

By completing this assignment, you will:

* Master high-performance computing environments for genomic data processing
* Develop proficiency in handling large-scale genomic datasets and file formats
* Apply command-line tools for genomic feature extraction and coordinate manipulation
* Compare and evaluate different sequence alignment algorithms
* Integrate multiple bioinformatics tools into a cohesive analytical pipeline
* Interpret results in the context of medical genomics and precision medicine

---

## **Part 1: High-Performance Computing with Human Genome Data (25 points)**

### **Objectives**
* Access and configure HPC resources for genomic analysis
* Download and process large genomic datasets
* Implement efficient file management strategies for big data

### **Instructions**

#### **Step 1: HPC Environment Setup**
```bash
# Access Palmetto2 cluster via OnDemand interface
# Navigate to: https://ondemand.rcd.clemson.edu/

# Request interactive session with:
# - 8 CPU cores
# - 64 GB memory  
# - 24-hour walltime
# - Appropriate partition for genomic analysis
```

#### **Step 2: Workspace Configuration**
```bash
# Navigate to scratch space
cd /scratch/$USER

# Create project directory structure
mkdir -p human_genome_pipeline/{data,results,scripts,logs}
cd human_genome_pipeline

# Set up environment variables
export GENOME_DIR=$PWD/data
export RESULTS_DIR=$PWD/results
export SCRIPTS_DIR=$PWD/scripts
```

#### **Step 3: Genomic Data Acquisition**
```bash
# Download human genome reference (GRCh38)
cd $GENOME_DIR
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download human cDNA sequences
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Download GTF annotation file
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Verify downloads and check file integrity
ls -lh *.gz
md5sum *.gz > download_checksums.txt
```

#### **Step 4: Data Processing and Quality Assessment**
```bash
# Decompress files efficiently
parallel gunzip ::: *.gz

# Generate comprehensive dataset statistics
echo "=== GENOME ASSEMBLY STATISTICS ===" > $RESULTS_DIR/dataset_summary.txt
echo "Genome file: $(basename $GENOME_DIR/*.fa)" >> $RESULTS_DIR/dataset_summary.txt
echo "Number of chromosomes/contigs: $(grep -c '^>' $GENOME_DIR/*.fa)" >> $RESULTS_DIR/dataset_summary.txt
echo "Total genome size: $(grep -v '^>' $GENOME_DIR/*.fa | wc -c) bp" >> $RESULTS_DIR/dataset_summary.txt

echo -e "\n=== cDNA STATISTICS ===" >> $RESULTS_DIR/dataset_summary.txt
echo "Number of cDNA sequences: $(grep -c '^>' $GENOME_DIR/*cdna*.fa)" >> $RESULTS_DIR/dataset_summary.txt

# Calculate N50 and other assembly metrics
seqkit stats $GENOME_DIR/*.fa >> $RESULTS_DIR/assembly_metrics.txt
```

### **Deliverables for Part 1:**
1. Screenshot of successful HPC session allocation
2. Complete dataset summary report
3. File size comparison (compressed vs. uncompressed)
4. Resource utilization log

---

## **Part 2: Genomic Feature Extraction and Coordinate Analysis (30 points)**

### **Objectives**
* Parse GTF files to extract genomic coordinates
* Manipulate genomic intervals and coordinate systems
* Generate feature-specific datasets for downstream analysis

### **Instructions**

#### **Step 1: GTF File Analysis**
```bash
cd $GENOME_DIR

# Examine GTF structure and content
head -20 Homo_sapiens.GRCh38.109.gtf
grep -v '^#' Homo_sapiens.GRCh38.109.gtf | cut -f3 | sort | uniq -c

# Extract comprehensive exon coordinates
awk '$3 == "exon" {print $1"\t"$4"\t"$5"\t"$7"\t"$9}' Homo_sapiens.GRCh38.109.gtf > $RESULTS_DIR/human_exons_detailed.bed

# Add header to BED file
echo -e "chromosome\tstart\tend\tstrand\tattributes" | cat - $RESULTS_DIR/human_exons_detailed.bed > temp && mv temp $RESULTS_DIR/human_exons_detailed.bed
```

#### **Step 2: Advanced Coordinate Manipulation**
```bash
# Extract gene coordinates with metadata
awk '$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, gene_id);
    match($9, /gene_name "([^"]+)"/, gene_name);
    match($9, /gene_biotype "([^"]+)"/, biotype);
    print $1"\t"$4"\t"$5"\t"gene_id[1]"\t"gene_name[1]"\t"biotype[1]"\t"$7
}' Homo_sapiens.GRCh38.109.gtf > $RESULTS_DIR/human_genes_annotated.bed

# Calculate feature statistics per chromosome
awk '$3 == "exon" {exons[$1]++} 
     $3 == "gene" {genes[$1]++} 
     END {
         print "Chromosome\tGenes\tExons\tExons_per_Gene"
         for (chr in genes) {
             printf "%s\t%d\t%d\t%.2f\n", chr, genes[chr], exons[chr], exons[chr]/genes[chr]
         }
     }' Homo_sapiens.GRCh38.109.gtf | sort -k1,1V > $RESULTS_DIR/chromosome_feature_stats.txt
```

#### **Step 3: Specialized Feature Extraction**
```bash
# Extract protein-coding genes only
awk '$3 == "gene" && $9 ~ /gene_biotype "protein_coding"/ {
    match($9, /gene_id "([^"]+)"/, gene_id);
    match($9, /gene_name "([^"]+)"/, gene_name);
    print $1"\t"$4"\t"$5"\t"gene_id[1]"\t"gene_name[1]"\t"$7
}' Homo_sapiens.GRCh38.109.gtf > $RESULTS_DIR/protein_coding_genes.bed

# Calculate gene length distribution
awk '$3 == "gene" {
    length = $5 - $4 + 1;
    lengths[NR] = length;
    sum += length;
    if (length > max) max = length;
    if (length < min || min == 0) min = length;
}
END {
    n = asort(lengths);
    median = (n % 2) ? lengths[(n+1)/2] : (lengths[n/2] + lengths[n/2+1])/2;
    printf "Gene Length Statistics:\n";
    printf "Total genes: %d\n", n;
    printf "Mean length: %.0f bp\n", sum/n;
    printf "Median length: %.0f bp\n", median;
    printf "Min length: %d bp\n", min;
    printf "Max length: %d bp\n", max;
}' Homo_sapiens.GRCh38.109.gtf > $RESULTS_DIR/gene_length_statistics.txt
```
 
### **Deliverables for Part 2:**
1. Annotated BED files for exons and genes
2. Chromosome-wise feature statistics
3. Gene length distribution analysis
4. Protein-coding gene subset

---

## **Part 3: Comparative Sequence Alignment Analysis (35 points)**

### **Objectives**
* Implement and compare BLAST and Smith-Waterman algorithms
* Analyze alignment sensitivity and specificity
* Evaluate computational performance trade-offs

### **Instructions**

#### **Step 1: Alignment Tool Installation and Setup**
```bash
cd $SCRIPTS_DIR

# Install BLAST+ suite
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.14.0+-x64-linux.tar.gz
export PATH=$PWD/ncbi-blast-2.14.0+/bin:$PATH

# Install EMBOSS for Smith-Waterman
module load emboss  # or install locally if needed

# Verify installations
blastn -version
water -help
```

#### **Step 2: Database Preparation**
```bash
# Create BLAST database from genome
makeblastdb -in $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            -dbtype nucl \
            -out $RESULTS_DIR/human_genome_db \
            -title "Human Genome GRCh38"

# Extract subset of cDNA sequences for testing (first 100 sequences)
head -200 $GENOME_DIR/Homo_sapiens.GRCh38.cdna.all.fa > $RESULTS_DIR/test_cdna_subset.fa

# Create smaller test datasets for Smith-Waterman
seqkit sample -n 10 $RESULTS_DIR/test_cdna_subset.fa > $RESULTS_DIR/sw_test_sequences.fa
```

#### **Step 3: BLAST Analysis**
```bash
# Run BLAST alignment with detailed output
blastn -query $RESULTS_DIR/test_cdna_subset.fa \
       -db $RESULTS_DIR/human_genome_db \
       -out $RESULTS_DIR/blast_results.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -num_threads 8

# Generate BLAST statistics
echo "BLAST Alignment Results Summary:" > $RESULTS_DIR/blast_summary.txt
echo "Total alignments: $(wc -l < $RESULTS_DIR/blast_results.txt)" >> $RESULTS_DIR/blast_summary.txt
echo "Average percent identity: $(awk '{sum+=$3; count++} END {print sum/count}' $RESULTS_DIR/blast_results.txt)" >> $RESULTS_DIR/blast_summary.txt
echo "Average alignment length: $(awk '{sum+=$4; count++} END {print sum/count}' $RESULTS_DIR/blast_results.txt)" >> $RESULTS_DIR/blast_summary.txt
```

#### **Step 4: Smith-Waterman Analysis**
```bash
# Run Smith-Waterman alignments for comparison
mkdir -p $RESULTS_DIR/smith_waterman

# Extract individual sequences for pairwise alignment
seqkit split -i $RESULTS_DIR/sw_test_sequences.fa -O $RESULTS_DIR/smith_waterman/

# Run Smith-Waterman alignments
for seq_file in $RESULTS_DIR/smith_waterman/*.fa; do
    base_name=$(basename $seq_file .fa)
    water -asequence $seq_file \
          -bsequence $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
          -gapopen 10 \
          -gapextend 0.5 \
          -outfile $RESULTS_DIR/smith_waterman/${base_name}_sw.txt \
          -aformat pair
done
```

#### **Step 5: Comparative Analysis**
```bash
# Create comprehensive comparison script
cat > $SCRIPTS_DIR/alignment_comparison.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import time

def parse_blast_results(blast_file):
    """Parse BLAST output file"""
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = pd.read_csv(blast_file, sep='\t', names=columns)
    return df

def analyze_alignment_performance():
    """Compare BLAST and Smith-Waterman performance"""
    results_dir = Path('../results')
    
    # Load BLAST results
    blast_df = parse_blast_results(results_dir / 'blast_results.txt')
    
    # Generate performance comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Percent Identity Distribution
    axes[0,0].hist(blast_df['pident'], bins=30, alpha=0.7, edgecolor='black')
    axes[0,0].set_xlabel('Percent Identity (%)')
    axes[0,0].set_ylabel('Frequency')
    axes[0,0].set_title('BLAST Percent Identity Distribution')
    
    # Plot 2: Alignment Length vs E-value
    axes[0,1].scatter(blast_df['length'], -np.log10(blast_df['evalue']), alpha=0.6)
    axes[0,1].set_xlabel('Alignment Length (bp)')
    axes[0,1].set_ylabel('-log10(E-value)')
    axes[0,1].set_title('Alignment Length vs Statistical Significance')
    
    # Plot 3: E-value distribution
    axes[1,0].hist(np.log10(blast_df['evalue']), bins=30, alpha=0.7, edgecolor='black')
    axes[1,0].set_xlabel('log10(E-value)')
    axes[1,0].set_ylabel('Frequency')
    axes[1,0].set_title('E-value Distribution')
    
    # Plot 4: Bit Score vs Percent Identity
    axes[1,1].scatter(blast_df['pident'], blast_df['bitscore'], alpha=0.6)
    axes[1,1].set_xlabel('Percent Identity (%)')
    axes[1,1].set_ylabel('Bit Score')
    axes[1,1].set_title('Bit Score vs Percent Identity')
    
    plt.tight_layout()
    plt.savefig(results_dir / 'alignment_analysis.png', dpi=300, bbox_inches='tight')
    
    # Generate summary statistics
    summary_stats = {
        'Total_Alignments': len(blast_df),
        'Mean_Percent_Identity': blast_df['pident'].mean(),
        'Median_Percent_Identity': blast_df['pident'].median(),
        'Mean_Alignment_Length': blast_df['length'].mean(),
        'High_Quality_Alignments': len(blast_df[blast_df['pident'] > 95]),
        'Significant_Alignments': len(blast_df[blast_df['evalue'] < 1e-10])
    }
    
    return summary_stats

if __name__ == "__main__":
    stats = analyze_alignment_performance()
    print("Alignment Analysis Summary:")
    for key, value in stats.items():
        print(f"{key}: {value:.2f}" if isinstance(value, float) else f"{key}: {value}")
EOF

# Run the analysis
cd $SCRIPTS_DIR
python3 alignment_comparison.py > $RESULTS_DIR/alignment_comparison_summary.txt
```

### **Deliverables for Part 3:**
1. BLAST and Smith-Waterman alignment results
2. Performance comparison analysis
3. Statistical significance evaluation
4. Computational efficiency assessment

---

## **Part 4: Integrated Pipeline and Medical Applications (10 points)**

### **Objectives**
* Integrate all modules into a cohesive analysis pipeline
* Apply results to medical genomics contexts
* Develop reproducible workflows

### **Instructions**

#### **Step 1: Pipeline Integration**
```bash
# Create master pipeline script
cat > $SCRIPTS_DIR/integrated_pipeline.sh << 'EOF'
#!/bin/bash
set -e

echo "Starting Integrated Human Genome Analysis Pipeline..."
echo "Timestamp: $(date)"

# Module 1: Data acquisition and processing
echo "=== MODULE 1: Data Processing ==="
bash module1_data_processing.sh

# Module 2: Feature extraction
echo "=== MODULE 2: Feature Extraction ==="
bash module2_feature_extraction.sh

# Module 3: Sequence alignment
echo "=== MODULE 3: Alignment Analysis ==="
bash module3_alignment_analysis.sh

# Generate final report
echo "=== GENERATING FINAL REPORT ==="
python3 generate_final_report.py

echo "Pipeline completed successfully!"
echo "Results available in: $RESULTS_DIR"
EOF

chmod +x $SCRIPTS_DIR/integrated_pipeline.sh
```

#### **Step 2: Medical Genomics Application**
Create a focused analysis on clinically relevant genes:

```bash
# Extract clinically relevant gene sets
grep -E "(BRCA1|BRCA2|TP53|EGFR|KRAS|PIK3CA|APC|MLH1|MSH2|MSH6)" \
     $RESULTS_DIR/protein_coding_genes.bed > $RESULTS_DIR/clinical_genes.bed

# Analyze these genes specifically
awk 'NR>1 {print $4}' $RESULTS_DIR/clinical_genes.bed | \
while read gene_id; do
    grep "$gene_id" $RESULTS_DIR/human_exons_detailed.bed | wc -l
done > $RESULTS_DIR/clinical_gene_exon_counts.txt
```
### **Deliverables for Part 4:**
1. Pipeline script as a plain text file.
2. Clinically relevant disease genes results.
---

## **Assessment and Deliverables**

### **Final Submission Requirements:**

1. **Technical Implementation (40%)**
   - Complete, executable scripts for all modules
   - Proper error handling and logging
   - Efficient resource utilization documentation

2. **Scientific Analysis (35%)**
   - Comprehensive results interpretation
   - Statistical analysis of alignment performance
   - Medical relevance discussion

3. **Documentation and Reproducibility (25%)**
   - Detailed README with pipeline instructions
   - Well-commented code
   - Complete results summary report

### **Submission Instructions:**
1. **Submit your Parts 1-4 deliverables files through Canvas in a single double compressed (tar + gzip = tarball). Please put each Part in a seperate directory.  
2 **Late submissions will be subject to the policy in the syllabus.

---

## **Resources and References**

### **Primary Resources:**
- [ENSEMBL Database](https://www.ensembl.org/)
- [NCBI BLAST Documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [EMBOSS Documentation](http://emboss.sourceforge.net/docs/)
- [Palmetto2 User Guide](https://docs.rcd.clemson.edu/palmetto/)


---

*This assignment is designed to provide hands-on experience with real-world computational genomics workflows. Success requires careful attention to detail, systematic approach to problem-solving, and integration of multiple bioinformatics concepts and tools.*
