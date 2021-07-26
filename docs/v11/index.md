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
| test statistic * | The statistic to use: <br><br> **t-test** <br> This is the standardized mean difference between the two classes.  It is the difference between the mean expression of class 1 and class 2 divided by the variability of expression, which is the square root of the sum of the standard deviation for each class divided by the number of samples in each class. <br><br> <center>μ A − μ B σ A n A + σ B n B</center> <br> where <br> μ is the average <br> σ is the standard deviation <br> n is the number of samples <br> ___ <br><br> **t-test (median)** <br> Same as t-test, but uses median rather than average. <br> ___ <br><br> **t-test (min std)** <br> Same as t-test, but enforces a minium value for sigma (minimal standard deviation). <br> ___ <br><br> **t-test (median, min std)** <br> Same as t-test, but uses median rather than average and enforces a minimum value for sigma (minimal standard deviation). <br> ___ <br><br> **SNR** <br> The signal-to-noise ratio is computed by dividing the difference of class means by the sum of their standard deviation. <br><br> <center> μ A − μ B σ A + σ B </center> <br> where μ is the average and σ is the standard deviation <br> ___ <br><br> **SNR (median)** <br> Same as SNR, but uses the median rather than average. <br> ___ <br><br> **SNR (min std)** <br> Same as SNR, but enforces a minimum value for sigma (minimal standard deviation). <br> ___ <br><br> **SNR (median, min std)** <br> Same as SNR, but uses median rather than average and enforces a minimum value for sigma (minimal standard deviation). <br> ___ <br><br> **Paired t-test** <br> The Paired T-Test can be used to analyze paired samples; for example, samples taken from patients before and after treatment. This test is used when the cross-class differences (e.g. the difference before and after treatment) are expected to be smaller than the within-class differences (e.g. the difference between two patients). For example if you are measuring weight gain in a population of people, the weights may be distributed from 90 lbs. to say 300 lbs. and the weight gain/loss (the paired variable) may be on the order of 0-30 lbs. So the cross-class difference ("before" and "after") is less than the within-class difference (person 1 and person 2). <br> The standard T-Test takes the mean of the difference between classes, the Paired T-Test takes the mean of the differences between pairs: <br><br> <center> X ¯ D − μ 0 s D / N </center> <br> where the differences between all pairs are calculated and XD is the average of the differences and SD the standard deviation.  μ0 is the mean difference between paired samples under the null hypothesis, typically 0. <br> _**Note**_: For the Paired T-Test, paired samples in the expression data file must be arranged by class, where the first samples in each class are paired, the second samples are paired, and so on. For example, sample pairs A1/B1, A2/B2 and A3/B3 would be ordered in an expression data file as A1, A2, A3, B1, B2, B3. Note that your data must contain the same number of samples in each class in order to use this statistic. | t-test |
| min std | The minimum standard deviation if test statistic includes min std option.  If σ is less than _min std_, σ is set to _min std_ . |
| number of permutations * | The number of permutations to perform (use 0 to calculate asymptotic p-values using the standard independent two-sample t-test).  ComparativeMarkerSelection uses a permutation test to estimate the significance (p-value) of the test statistic score. The number of permutations you specify depends on the number of hypotheses being tested and the significance level that you want to achieve (3). If the data set includes at least 10 samples per class, use the default value of 10000 permutations to ensure sufficiently accurate p-values. <br> If the data set includes fewer than 10 samples in any class, permuting the samples cannot give an accurate p-value. Specify a value of 0 permutations to use asymptotic p-values instead. In this case, ComparativeMarkerSelection computes p-values assuming the test statistic scores follow Student's t-distribution (rather than using the test statistic to create an empirical distribution of the scores). Asymptotic p-values are calculated using the p-value obtained from the standard independent two-sample t-test. | 10000 |
| log transformed data * | Whether the input data has been log transformed.  By default ComparativeMarkerSelection expects the data in the input file to not be log transformed. Some of the calculations such as the fold are not accurate when log transformed data is provided and not indicated. To indicate that your data is log transformed, set this parameter to “yes”. | no |
| complete * | Whether to perform all possible permutations.  When the complete parameter is set to yes, ComparativeMarkerSelection ignores the number of permutations parameter and computes the p-value based on all possible sample permutations. Use this option only with small data sets, where the number of all possible permutations is less than 1000. | no |
| balanced * | Whether to perform balanced permutations.  When the balanced parameter is set to yes, ComparativeMarkerSelection requires an equal and even number of samples in each class (e.g. 10 samples in each class, not 11 in each class or 10 in one class and 12 in the other). | no |
| random seed * | he seed of the random number generator used to produce permutations. | 779948241 |
| smooth p values * | Whether to smooth p-values by using the Laplace’s Rule of Succession. By default, smooth p values is set to yes, which means p-values are always less than 1.0 and greater than 0.0. | yes |
| phenotype test | Tests to perform when cls file has more than two classes: one-versus-all, all pairs. (Note: The p-values obtained from the one-versus-all comparison are not fully corrected for multiple hypothesis testing.) | one versus all |
| output filename * | The name of the output file. | <input.file_basename>.comp.marker.odf |

\*  required

## Input Files

1. _input file_: <br>
   [GCT](https://genepattern.org/file-formats-guide#GCT) - or [RES](https://genepattern.org/file-formats-guide#RES) -formatted file containing the expression dataset.
2. _cls_ file: <br>
    [CLS](https://genepattern.org/file-formats-guide#CLS) -formatted class file.
3.  _confounding variable cls file_:<br> 
    [CLS](https://genepattern.org/file-formats-guide#CLS) -formatted file containing the confounding variable.   
    
## Output Files

1. _\<input.file_basename\>.comp.marker.odf_ <br>
   [ODF](https://genepattern.org/file-formats-guide#ODF) -formatted file containing the following columns:<br> <br>
   * _**Rank**_: The rank of the feature within the dataset based on the value of the test statistic. If a two-sided p-value is computed, the rank is with respect to the absolute value of the statistic. 
   * _**Feature**_: The feature name. 
   * _**Description**_: The description of the feature. 
   * _**Score**_: The value of the test statistic. 
   * _**Feature P**_: The feature-specific p-value based on permutation testing. 
   * _**Feature P Low**_: The estimated lower bound for the feature p-value. 
   * _**Feature P High**_: The estimated upper bound for the feature p-value. 
   * _**FDR (BH)**_: An estimate of the false discovery rate by the Benjamini and Hochberg procedure (1). The FDR is the expected proportion of erroneous rejections among all rejections. 
   * _**Q Value**_: An estimate of the FDR using the procedure developed by Storey and Tibshirani (6). 
   * _**Bonferroni**_: The value of the Bonferroni correction applied to the feature specific p-value. 
   * _**maxT**_: The adjusted p-values for the maxT multiple testing procedure described in Westfall (7), which provides strong control of the FWER. 
   * _**FWER (Family Wise Error Rate)**_: The probability of at least one null hypothesis/feature having a score better than or equal to the observed one. This measure is not feature-specific. 
   * **Fold Change**: The class zero mean divided by the class one mean. 
   * **Class Zero Mean**: The class zero mean. 
   * **Class Zero Standard Deviation**: The class zero standard deviation. 
   * **Class One Mean**: The class one mean. 
   * **Class One Standard Deviation**: The class one standard deviation. 
   * **k**: If performing a two-sided test or a one-sided test for markers of class zero, the number of permuted scores greater than or equal to the observed score. If testing for markers of class one, then the number of permuted scores less than or equal to the observed score.  

## Example Data
<!-- provide links to example data so that users can see what input & output should look like and so that they and we can use it to test -->

Input:  
[all_aml_test.gct](https://datasets.genepattern.org/data/all_aml/all_aml_test.gct) <br>
[all_aml_test.cls](https://datasets.genepattern.org/data/all_aml/all_aml_test.cls)

Output:  
[created_file_ground_truth.txt](https://datasets.genepattern.org/data/all_aml/all_aml_test.comp.marker.odf)


## Requirements

Requires the [genepattern/comparativemarkerselection:v10.5 Docker image](https://hub.docker.com/layers/genepattern/comparativemarkerselection/v10.5/images/sha256-5c66f127fd7a32b8de3fe03fb2cf517b1f0306b8f3257cc5056e37a97741876b?context=explore).

## License

`ComparativeMarkerSelection` is distributed under a modified BSD license available at [https://github.com/genepattern/ComparativeMarkerSelection/blob/v11/LICENSE](https://github.com/genepattern/ComparativeMarkerSelection/blob/v11/LICENSE)

## Version Comments

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 11 | 2021-07-20 | Updated to run in its own Docker container - allowing q values to be computed |
| 10 | 2013-12-04 | Updated to html doc |
| 9 | 2012-03-26 | Changed default number of permutations to 10000 |
| 8 | 2011-08-30 | added parameter to specify whether data is log transformed |
| 7 | 2010-05-28 | Made improvements to error messages |
| 6 | 2009-12-30 | Fixed bug with using res file with paired t-test |
| 5	| 2008-10-24 | Added Paired T-Test |
| 4 | 2008-02-19 | Added Paired T-Test |
| 3	 | 2006-03-03 | Added additional metrics |
| 2 | 2005-06-08 | Added restricted permutations option and maxT p-value |
