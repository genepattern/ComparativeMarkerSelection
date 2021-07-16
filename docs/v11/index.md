# ComparativeMarkerSelection (v11)

**Description**: Identify differentially expressed genes that can discriminate between distinct classes of samples.

**Authors**: Joshua Gould, Gad Getz, Stefano Monti - Cancer Program, Broad Institute; Barbara Hill - Mesirov Lab, Broad Institute

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)

## Summary

When analyzing genome-wide transcription profiles derived from microarray or RNA-seq experiments, the first step is often to identify differentially expressed genes that can discriminate between distinct classes of samples (usually defined by a phenotype, such as tumor or normal).  This process is commonly referred to as marker (or feature) selection.  Marker genes are identified by calculating, for each profiled gene, a test statistic (e.g., t-test) which assesses correlation of the gene's expression profile with a class template.  If the value of the test statistic for a specific gene, and thus the degree of differential expression presented by that gene, is significantly greater than what one would expect to see under the null hypothesis (gene is not differentially expressed between classes), that gene is identified as a statistically significant marker gene.  

The ComparativeMarkerSelection module takes as input a dataset of expression profiles from samples belonging to two classes and, implementing the statistical tests described above, identifies marker genes which discriminate between the classes.

The ComparativeMarkerSelection module includes several approaches to determine the features that are most closely correlated with a class template and the significance of that correlation.  The module computes significance values for features using several metrics, including FDR(BH), Q-Value, maxT, FWER, Feature-Specific P-Value, and Bonferroni. The results from the ComparativeMarkerSelection algorithm can be viewed with the ComparativeMarkerSelectionViewer.   ExtractComparativeMarkerResults creates a derived dataset and feature list file from the results of ComparativeMarkerSelection.

By default, ComparativeMarkerSelection expects the data in the input file to **not** be log transformed. Some of the calculations such as the fold change are not accurate when log transformed data is provided and not indicated. To indicate that your data is log transformed, be sure to set the “log transformed data” parameter to “yes”. Also, ComparativeMarkerSelection requires **at least three** samples per class to run successfully.

## Algorithm

The analytic module takes as input a dataset of expression profiles from samples belonging to two phenotypes. If a dataset contains more than two phenotypes, then there is the option to perform all pairwise comparisons or all one-versus-all comparisons. A test statistic (e.g. t-test) is chosen to assess the differential expression between the two classes of samples. Note that technical and biological replicates are handled the same way as independent samples. The significance (nominal P-value) of marker genes is computed using a permutation test, which is a commonly used method for assessing the significance of marker genes; see (4) for details.

Selecting class markers is a particular instance of the general multiple hypothesis testing problem. Since several thousand hypotheses are usually tested at once (one per gene), the nominal P-values have to be corrected to account for the increased number of potential false positives. For example, if we test 20,000 genes for differential expression, a nominal P-value threshold of 0.01 would only ensure that the expected number of false positives is <200 (0.01 x 20,000).  ComparativeMarkerSelection includes several methods of correcting for multiple hypothesis testing, including FDR(BH), Q-Value, maxT, FWER, Feature-Specific P-Value, and Bonferroni;  (4) describes their applicability.

## References

1. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. _Journal of the Royal Statistical Society. Series B (Methodological)._ 1995;57(1):289-300. 
2. Golub T, Slonim D, et al. Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression. _Science_. 1999;286:531-537. 
3. Good P. _Permutation Tests: A Practical Guide for Testing Hypotheses_, 2nd Ed. New York: Springer-Verlag. 2000. 
4. Gould J, Getz G, Monti S, Reich M, Mesirov JP. Comparative gene marker selection suite. _Bioinformatics_. 2006;22;1924-1925, doi:10.1093/bioinformatics/btl196. 
5. Lu J, Getz G, Miska E, et al. MicroRNA Expression Profiles Classify Human Cancers. _Nature_. 2005;435:834-838. 
6. Storey JD, Tibshirani R. Statistical significance for genomewide studies. _PNAS_. 2003;100(16):9440-9445. 
7. Westfall PH, Young SS. Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment, in _Wiley Series in Probability and Statistics_. New York: Wiley, 1993.

## Source Links
* [ComparativeMarkerSelection v2 source repository](https://github.com/genepattern/ComparativeMarkerSelection/tree/v11)
* ComparativeMarkerSelection v2 uses the [genepattern/comparativemarkerselection:v10.5 Docker image](https://hub.docker.com/layers/genepattern/comparativemarkerselection/v10.5/images/sha256-5c66f127fd7a32b8de3fe03fb2cf517b1f0306b8f3257cc5056e37a97741876b?context=explore)
* [The Dockerfile used to build that image is here.](https://hub.docker.com/r/genepattern/comparativemarkerselection/dockerfile)

## Parameters
<!-- short description of the module parameters and their default values, as well as whether they are required -->

| Name | Description <!--short description--> | Default Value |
---------|--------------|----------------
| input file * |  The input file.  [GCT](https://genepattern.org/file-formats-guide#GCT), [RES](https://genepattern.org/file-formats-guide#RES) <br><br> Note the following constraints: <br><br> - If the expression data contains duplicate identifiers, ComparativeMarkerSelection generates the error message: "An error occurred while running the algorithm." The UniquifyLabels module provides one way of handling duplicate identifiers. <br> - If the expression data contains fewer than three samples per class, ComparativeMarkerSelection appears to complete successfully but test statistic scores are not shown in the results. <br> - If the expression data contains missing values, ComparativeMarkerSelection completes successfully but does not compute test statistic scores for rows that contain missing values. <br> - If your data is log transformed, you will need to set the "log transformed data" parameter above to "yes".Note that if your data is log transformed, you will need to set the "log transformed data" parameter below to "yes". |
| cls file * | The class file. [CLS](https://genepattern.org/file-formats-guide#CLS) <br><br> ComparativeMarkerSelection analyzes two phenotype classes at a time. If the expression data set includes samples from more than two classes, use the phenotype test parameter to analyze each class against all others (one-versus-all) or all class pairs (all pairs). |
| confounding variable cls file   | The class file containing the confounding variable.  [CLS](https://genepattern.org/file-formats-guide#CLS) <br><br> If you are studying two variables and your data set contains a third variable that might distort the association between the variables of interest, you can use a confounding variable class file to correct for the affect of the third variable. For example, the data set in Lu, Getz, et. al. (2005) contains tumor and normal samples from different tissue types. When studying the association between the tumor and normal samples, the authors use a confounding variable class file to correct for the effect of the different tissue types. <br><br> The phenotype class file identifies the tumor and normal samples: <br><br>75 2 1 <br> \# Normal Tumor <br> 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 <br><br> The confounding variable class file identifies the tissue type of each sample: <br><br> 75 6 1 <br> \# colon kidney prostate uterus human-lung breast <br> 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 <br><br> Given these two class files, when performing permutations, ComparativeMarkerSelection shuffles the tumor/normal labels only among samples with the same tissue type.  |
| test direction * | The test to perform.  By default, ComparativeMarkerSelection performs a two-sided test; that is, the test statistic score is calculated assuming that the differentially expressed gene can be up-regulated in either phenotype class. Optionally, use the test direction parameter to specify a one-sided test, where the differentially expressed gene must be up-regulated for class 0 or for class 1.  | 2 Sided |
| test statistic * | The statistic to use: <br> 

\*  required

## Input Files
<!-- longer descriptions of the module input files. Include information about format and/or preprocessing...etc -->

1. filename  
    A long form explanation of the parameter. For example: This is the file which will be read in by the python script and to which text will be added, if add_custom_message is set to true. The parameter expects a text file with a .txt extension (e.g. file.txt)
    
## Output Files
<!-- list and describe any files output by the module -->

1. \<output_filename\>.txt  
    The input file plus any text you added, if you chose to add text.
2. stdout.txt
    This is standard output from the Python script. Sometimes helpful for debugging.

## Example Data
<!-- provide links to example data so that users can see what input & output should look like and so that they and we can use it to test -->

Input:  
[data_placeholder.txt](https://github.com/genepattern/ExampleModule/blob/v2/data/data_placeholder.txt)

Output:  
[created_file_ground_truth.txt](https://github.com/genepattern/ExampleModule/blob/v2/gpunit/output/basic_test/created_file_ground_truth.txt)


## Requirements
<!--This section is typically used to list any special requirements for running the module, such as, language/operating system requirements and Docker images. -->

Requires the [genepattern/example-module:2 Docker image](https://hub.docker.com/layers/150060459/genepattern/example-module/2/images/sha256-ae4fffff67672e46b251f954ad226b7ad99403c456c1c19911b6ac82f1a27f2f?context=explore).

## License

`ExampleModule` is distributed under a modified BSD license available at [https://github.com/genepattern/ExampleModule/blob/v2/LICENSE.](https://github.com/genepattern/ExampleModule/blob/v2/LICENSE)

## Version Comments
<!--For each version of a module, provide a short comment about what was changed in the new version of a module. Version comments consist of 3 parts: a date, a version number, and a short description. The date should be the release date of that version of the module, and the version number should match the version of the module for which it corresponds to. The description can be short, but should be informative (e.g. "added support for log transformed data", or "fixed bug with out of memory exception"). When a user views the documentation, all version comments up to and including the current version will be displayed, and act as a short version history for the module. -->

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
|  1.4  | May 17, 2021 | Added all GenePattern Team module release requirements and renamed as ExampleModule, from ABasicModule. |
| 1 | May 1, 2018 | Initial version for team use. |
