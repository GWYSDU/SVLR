SVLR is a tool designed to detect structural variation with long reads sequenced by the third generation sequencing (PacBio or Oxford Nanopore). 
It detects seven types of SVs (10bp+) which include deletion, insertion, inversion, tandem duplication, interspersed duplication, cut & paste isertion, block replacement, block interchange, and translocation by following four steps: SV signatures marking, SV signatures Clustering, SV clusters combination, SV clusters optimization. 
The current version of SVLR utilizes the alignments from NGMLR(.sam) or LAST(.maf).
If you have any susgetions or problems when you apply SVLR to detect structure variants, please contact to us by: gwy@mail.sdu.edu.cn.






How to build SVLR?
==================
==```wget https://codeload.github.com/GWYSDU/SVLR/zip/master     
unzip SVLR-master.zip  
cd /SVLR-master  
javac ComputeOverlap.java DataInput.java getCandidatemetrics.java Global.java MaximalCliquesWithPivot.java Parameter_Setting.java SV_Candidate.java SV_Candidate_Filter.java SVsignature.java SV_Signature_Cluster.java SV_Signature_Distinguish.java T_I_C_SVClusters_Distinguish.java```===
