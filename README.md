# pumilioAnalysis
Analysis of Bocas del Toro exome sequecing Oophaga pumilio samples

## 1 Process reads:

- Trimmomatic:
  - Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
  - Remove leading low quality or N bases (below quality 3) (LEADING:3)
  - Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
  - Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
  - Drop reads below the 75 bases long (MINLEN:75)
- Prinseq:
  - Remove low complexity reads
  - Remove overrepresented reads
- Fastqc
  - Asses quality of reads

## 2 Mapping:

- Bwa mem
- Samtools:
  - Samfile flags to exclude:
    - read unmapped (0x4)
    - mate unmapped (0x8)
    - not primary alignment (0x100)
    - read fails platform/vendor quality checks (0x200)
    - read is PCR or optical duplicate (0x400)
    - 1804
