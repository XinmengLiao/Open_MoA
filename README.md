# Open MoA


## 1. Introduction

Open MoA is a potent tool for identifying the underlying drug mechanism of actions (MoAs), core biological pathways, and key targets.

With the transcriptomic data, by inputting the known starting point and endpoints, Open MoA is able to give out the significant confidence score for each interaction in the context-specific subnetworks, thus leading to the identification of the most possible pathway of drug action and key candidates. 


## 2. Input and output data

(1) Input data:

1. Context-specific transcriptomic data

2. The starting point of the expected biological reaction.

3. The endpoint of the expected biological reaction.

(2) Output data:

1. Context-specific weighted subnetwork with reliable confidence scores.

2. The predicted core pathway (the shortest path) from the provided starting point to the endpoint.

3. Significant targets and interactions of the expected biological context.


## 3. Datasets

(1) DrugBank

Drug-target data of homo sapiens from DrugBank v5.1.9 is used in Open MoA. The screened Drug-targets dataset can be found in `Datasets/DrugBank.csv`

(2) STRING

Experimental and transferred experimental protein-protein interactions (PPIs) of homo sapiens from STRING database v11.5 are used in Open MoA. Only the PPIs with high confidence combined scores are retained. The screened PPIs dataset can be found in `Datasets/STRING.csv`

(3) RegNetowrk

Regulatory interactions of homo sapiens from the RegNetwork database are used in Open MoA. The ones with miRNA are filtered out. The screened dataset can be found in `Datasets/RegNetwork.csv`

(4) The reference integrated network (IN)

Interactions from DrugBank, STRING, and RegNetwork databases are combined together as `Networks/reference_IN.csv`. The network contains the names and types of two vertexes, interaction type, combined scores, and penalty score. The reference IN is an unweighted network, thus, all penalty scores are set to be 1.00. Since only the STRING database provides the combined scores for each PPI, only PPIs have exact combined scores, while the others are set to be 0.

(5) HepG2 expressed genes

The gene transcriptomic data of the HepG2 cell line from the Human Protein Atlas and the Cancer Cell Line Encyclopedia can be found in `Datasets/HepG2_genes.csv`

(6) The HepG2 cell-specific subnetwork

The reference integrated network is screened with genes expressed in the HepG2 cell line to generate a HepG2 cell-specific subnetwork. It can be found in `Datasets/HepG2_specific.csv`


## 4. The algorithms of Open MoA

(1) Construction of the reference IN and the HepG2-specific IN

The reference IN is the simple combination of three interaction datasets by `reference_IN <- rbind(DrugBank, STRING, RegNetwork)` in R. The reference IN is screened under Only the HepG2 expressed genes are kept for establishing the HepG2-specific IN. Based on the aim of study, the transcriptomic data of other cell lines could be used for cell line-specific IN.

(2) Construction of the context-specific weighted subnetworks

The transcriptomic data of the context target will be used in this step. The starting point should be the target, while the endpoints should be the perturbed genes of the target. The FDR values of the perturbation will be the weight for calculation. For instance, as we constructed a TGFβ1-specific weighted subnetwork in our study, the transcriptomic data of TGFβ1 were first downloaded from the Connectivity Map and the FDR values of each perturbed genes were calculated. Subsequently, the TGFβ1 was set as the starting point, while all the perturbed genes were set as the endpoints. The FDR values were used for weight calculation. Eventually, each edge in the context-specific weighted subnetwork will have a computed penalty score.

The test data of TGFβ1 transcriptomic data are located in `Test/TGFB1.csv`

The  Detailed process is shown in `Algorithms/Weights_Calculation.R`

(3) Identification of the potential core pathway (the shortest path)

With the context-specific weighted subnetworks, the potential core pathways could then be predicted by Open MoA. The shortest path function from `igaph` package in R is used for discovering the shortest path. The penalty scores of edges are used as the weight parameter.

The  Detailed process is shown in `Algorithms/Shortest_Path.R`


## 5. Following analyses for the shortest path and the context-specific centric subnetworks.

The following analyses were conducted by these tools and packages:

(1) KEGG enrichment analysis: The Database for Annotation, Visualization, and Integrated Discovery (DAVID) Website.

(2) OMIM analysis: `Enrichr v3.1` in R.

(3) Pathway reconstruction: `rWikipathways package v1.14.0` in R.

(4) Centric subnetwork visualisation: Cytoscape v3.9.1.
