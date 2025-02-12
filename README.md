# Genomic insights into sustained mammalian transmission of avian influenza A viruses
<img src="/img/intro.png" width="300">

Very few lineages of influenza A virus (IAV) have evolved sustained transmission in mammals. Why this is so remains largely unknown and the possibility of avian IAVs achieving sustained mammalian transmission is a pervasive concern. All available whole genomes of IAVs were assessed for their GC content and GC containing dinucleotide frequencies. Persistent mammalian lineages showed declining trends and could be accurately separated from avian lineages. Early viruses of the persistent mammalian lineages were also distinguishable, by reduced GC related content, from IAVs in birds and sporadic mammalian infections. Of note, the genotype of recent highly pathogenic 2.3.4.4b H5 viruses mostly detected in gulls also showed lower GC related content, and subsequently infected mink and foxes with evident mammalian transmissibility. The findings here suggest that the pattern of reduced GC related content in avian IAVs is one of predictors of enhanced ability to persist in mammals and should be included in risk assessment tools for pandemic influenza.

# Content
The repo contains all analyzed data in this study including the phylogenetic trees of eight protein coding regions and H5 gene for recently emerging H5Nx, the summary of statistics for genomic GC content, GC dinucleotide frequency and SVM model files. The influenza A virus genomic data are downloaded from GISAID (https://www.gisaid.org/). 

# Risk assessment tool
+ The tool is publicly available at https://iav-transmission.org/.
+ It is designed to assess the risk of sustained transmission in mammals for avian or recently zoonotic influenza A viruses.
+ It classifies mammalian IAVs as from either persistent (positive class) or sporadic (negative class) mammalian infections that achieved high accuracies of 99.61% for training (5-fold cross validation balanced accuracy = 99.13%) and 99.16% for testing.  
+ It was developed by a support vector machine (SVM) model implemented in LIBSVM with linear kernel (https://github.com/cjlin1/libsvm).
+ It is easy to use by uploading eight protein coding regions (HA, NA, NP, PA, PB1, PB2, M1 and NS1) or eight nucleotide segments (HA, NA, NP, PA, PB1, PB2, MP and NS) for this risk assessment.
 
