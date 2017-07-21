# Heterogeneous BLAST (H-BLAST) 
Heterogeneous BLAST (H-BLAST) is a fast parallel search tool for a 
heterogeneous computer that couples CPUs and GPUs, to accelerate 
BLASTX and BLASTP - basic modules of NCBI-BLAST. H-BLAST employs a 
locally decoupled seed-extension algorithm to take advantages of 
GPUs, and offers a performance tuning mechanism for better efficiency 
among various CPUs and GPUs combinations. H-BLAST produces identical 
alignment results as NCBI-BLAST and its computational speed is much 
faster than that of NCBI-BLAST.

H-BLAST has been tested on different software and hardware settings.
For example, CentOS Linux (v5.4, v6.0 and v7.1) and Ubuntu Linux 
(v11.04) were used as H-BLAST's benchmark software platforms;
NVIDIA Tesla C2050 GPU card, NVIDIA Tesla K20/K40 GPU card and NVIDIA 
Geforce GTX 560 GPU card were used as H-BLAST's benchmark hardware 
platforms.

H-BLAST is free for academic and non-commercial use. The current implementation 
is based on NCBI-BLAST, and immigrates some code segments from GPU-BLAST. 

Please cite the authors in any work or product based on this material:
Weicai Ye, Ying Chen, Yongdong Zhang, Yuesheng Xu; H-BLAST: a fast protein sequence 
alignment toolkit on heterogeneous computers with GPUs. Bioinformatics 2017; 33 (8): 
1130-1138. doi: 10.1093/bioinformatics/btw769.

For any questions and feedback about H-BLAST, contact cai_rcy@163.com 
or lnszyd@mail.sysu.edu.cn.

Features:
* H-BLAST accelerates BLASTX and BLASTP - basic modules of NCBI-BLAST.
* H-BLAST can take advantage of multiply CPU cores and GPUs in a computer.
* H-BLAST offers a performance tuning mechanism for better efficiency 
  among various CPUs and GPUs combinations, including manual and automatic 
  tuning.
* The GPU memory demand for H-BLAST is low, that for on-the-fly results and 
  alignment results are as less as 1/8 and 1/2 of those in GPU-BLAST, 
  respectively. Even a Geforce GTX 560 GPU card with 1GB ram can be used 
  for H-BLAST to accelerate BLASTX.
* H-BLAST v1.1 does not support PSI BLAST. The substitution matrix is 
  fixed to the Blosum62 matrix. 

Requirement:
* Nvidia GPU card with compute capability >= 2.0
* CUDA version >= 5.5.
* NCBI-BLAST v2.2.28+
