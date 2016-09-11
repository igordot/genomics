# RNA-Seq Strand

| | forward (transcript) | reverse (rev comp of transcript) |
|:----------|:----------|:----------|
| TopHat/Cufflinks | fr-secondstrand | fr-firststrand |
| STAR | 1st read strand | 2nd read strand |
| Picard CollectRnaSeqMetrics | FIRST_READ_TRANSCRIPTION_STRAND | SECOND_READ_TRANSCRIPTION_STRAND |
| htseq-count | yes | reverse |
| subread featureCounts | 1 | 2 |
| RSEM | forward-prob 1 | forward-prob 0 |
| Library Kit | Illumina ScriptSeq | Illumina TruSeq Stranded Total RNA |

***

Strand-specific protocols
([ref](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html)):
![](http://3.bp.blogspot.com/-BkupUsIrnXk/UBbmmmx6T8I/AAAAAAAAAUU/_rcrd_ahT48/s1600/strand.png)

***

Three widely used protocols for strand-specific RNA sequencing
([ref](http://www.nature.com/neuro/journal/v17/n11/full/nn.3814.html)):
![](http://www.nature.com/neuro/journal/v17/n11/images/nn.3814-F3.jpg)

***

Illumina TruSeq Stranded Total RNA Kit
([ref](https://www.abmgood.com/marketing/knowledge_base/next_generation_sequencing_experimental_design.php)):
![](https://www.abmgood.com/marketing/knowledge_base/img/NGS/Next_Generation_Sequencing_NGS_TruSeq_Stranded_Total_RNA.png)
