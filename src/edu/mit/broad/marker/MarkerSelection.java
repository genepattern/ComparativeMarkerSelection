package edu.mit.broad.marker;
import java.io.*;
import java.io.FileReader;

import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;
import org.genepattern.data.expr.*;

import org.genepattern.data.matrix.*;
import org.genepattern.io.expr.*;

import edu.mit.broad.marker.permutation.*;
import org.genepattern.module.AnalysisUtil;
import org.genepattern.stats.*;

/**
 * @author     Joshua Gould
 * @created    September 17, 2004
 */
public class MarkerSelection {

   /**  input expression values */
   DoubleMatrix2D dataset;

   /**  input classes */
   ClassVector classVector;

   /**  number of permutations to perform */
   int numPermutations;
   final static int CLASS_ZERO_GREATER_THAN_CLASS_ONE = 0;
   final static int CLASS_ZERO_LESS_THAN_CLASS_ONE = 1;
   final static int TWO_SIDED = 2;

   /**  whether to do one-sided >, one-sided <, or two-sided test */
   int testDirection;

   int metric;
   final static int T_TEST = 0;
   final static int SNR = 1;
   final static int T_TEST_MEDIAN = 2;
   final static int SNR_MEDIAN = 3;
   final static int T_TEST_MIN_STD = 4;
   /**  minimum standard deviation when metric==T_TEST_MIN_STD */
   double minStd;
   /**  whether permutations are balanced */
   boolean balanced;

   /**  whether to compute all possible permutations */
   boolean complete;
   /**  whether to read permutations from a file */
   boolean readPermutations = false;
   /**  reader when reading permutations from file */
   BufferedReader testClassPermutationsReader;

   /**  name of output file */
   String outputFileName;

   /**  unpermuted scores */
   double[] scores;

   /**  number of rows in input data */
   int N;

   /**  Whether to fix the standard deviation, as is done in GeneCluster */
   boolean fixStdev;
   ITestStatistic statisticalMeasure;

   /**  Input dataset file name */
   String datasetFile;
   /**  Input cls file name */
   String clsFile;

   /**  Covariate assignments */
   ClassVector covariate;

   /** Seed for generating permutations */
   int seed;
   static Debugger debugger;
   boolean seedUsed = false;

   public MarkerSelection(String datasetFile, String clsFile,
         int _numPermutations, int _side,
         String _outputFileName, boolean _balanced,
         boolean complete, boolean fixStdev, int metric, double minStd, int seed, String confoundingClsFile) {
      this.datasetFile = datasetFile;
      this.clsFile = clsFile;
      this.minStd = minStd;
      this.seed = seed;
      IExpressionDataReader reader = AnalysisUtil.getExpressionReader(datasetFile);
      ExpressionData expressionData = AnalysisUtil.readExpressionData(reader, datasetFile);

      this.dataset = expressionData.getExpressionMatrix();
      if(confoundingClsFile != null) {
         covariate = AnalysisUtil.readClassVector(confoundingClsFile);
          /*  ClassVector[] cv = new edu.mit.broad.internal.MultiClassReader().read(clsFile, expressionData);
            this.classVector = cv[0];
            if(cv.length > 1) {
               covariate = cv[1];
               for(int i = 2; i < cv.length; i++) {
                  covariate = covariate.union(cv[i]);
               }
            }
            */
      } 
		
		this.classVector = AnalysisUtil.readClassVector(clsFile);
      AnalysisUtil.checkDimensions(dataset, classVector);
      if(classVector.getClassCount() != 2) {
         AnalysisUtil.exit("Class file must contain 2 classes.");
      }
      this.numPermutations = _numPermutations;
      this.testDirection = _side;
      this.balanced = _balanced;
      this.outputFileName = _outputFileName;
      this.complete = complete;
      this.N = dataset.getRowCount();
      this.fixStdev = fixStdev;
      this.metric = metric;
      String permutationsFile = System.getProperty("edu.mit.broad.marker.perm");
      if(permutationsFile != null) {
         System.out.println("reading permutations from file " + permutationsFile);
         readPermutations = true;
         initFile(permutationsFile);
      }
      if(metric == T_TEST) {
         statisticalMeasure = new TestStatistics.TTest(fixStdev);
      } else if(metric == T_TEST_MEDIAN) {
         statisticalMeasure = new TestStatistics.TTestMedian(fixStdev);
      } else if(metric == SNR) {
         statisticalMeasure = new TestStatistics.SNR(fixStdev);
      } else if(metric == SNR_MEDIAN) {
         statisticalMeasure = new TestStatistics.SNRMedian(fixStdev);
      } else if(metric == T_TEST_MIN_STD) {
         if(minStd <= 0) {
            AnalysisUtil.exit("Minimum standard deviation must be greater than zero.");
         }
         statisticalMeasure = new TestStatistics.TTestMinStd(minStd, fixStdev);
      } else {
         AnalysisUtil.exit("Unknown test statistic");
      }

      if(testDirection != CLASS_ZERO_GREATER_THAN_CLASS_ONE && testDirection != CLASS_ZERO_LESS_THAN_CLASS_ONE && testDirection != TWO_SIDED) {
         AnalysisUtil.exit("Unknown test direction.");
      }

      computePValues();

      if(readPermutations) {

         if(testClassPermutationsReader != null) {

            try {
               testClassPermutationsReader.close();
            } catch(Exception e) {
            }
         }
      }
   }


   public static void main(String[] args) {
      try {
         run(args);  
      } catch(Throwable t) {
         AnalysisUtil.exit("An error occurred while running the algorithm.");
      }
   }
   
   
   static void run(String[] args) {
      String datasetFile = args[0];
      String clsFile = args[1];
      int _numPermutations = 0;
      try {
         _numPermutations = Integer.parseInt(args[2]);
      } catch(NumberFormatException nfe) {
         AnalysisUtil.exit("number of permutations is not an integer");
      }
      int testDirection = Integer.parseInt(args[3]);
      String outputFileName = args[4];
      boolean balanced = Boolean.valueOf(args[5]).booleanValue();
      boolean complete = Boolean.valueOf(args[6]).booleanValue();
      boolean fixStdev = Boolean.valueOf(args[7]).booleanValue();
      int metric = Integer.parseInt(args[8]);
      int seed = 0;
      try {
         seed = Integer.parseInt(args[9]);
      } catch(NumberFormatException nfe) {
         AnalysisUtil.exit("random seed is not an integer");
      }
      double minStd = -1;
      String confoundingClsFile = null;
      for(int i = 10; i < args.length; i++) {
         String arg = args[i].substring(0, 2);
         String value = args[i].substring(2, args[i].length());
         if(value.equals("")) {
            continue;
         }
         if(arg.equals("-m")) {
            try {
               minStd = Double.parseDouble(value); 
            } catch(NumberFormatException nfe) {
               AnalysisUtil.exit("minimum standard deviation is not a number");  
            }
         } else if(arg.equals("-c")) {
            confoundingClsFile = value;
         }
      }
     
   
      new MarkerSelection(datasetFile,
            clsFile,
            _numPermutations, testDirection, outputFileName, balanced,
            complete, fixStdev, metric, minStd, seed, confoundingClsFile);
   }


   void computePValues() {
      int[] classZeroIndices = classVector.getIndices(0);
      int[] classOneIndices = classVector.getIndices(1);
      if(balanced) {
         if(classZeroIndices.length != classOneIndices.length) {
            AnalysisUtil.exit(
                  "The number of items in each class must be equal for balanced permutations.");
         }
         if((classZeroIndices.length % 2) != 0) {
            AnalysisUtil.exit(
                  "The number of items in class 0 must be an even number for balanced permutations.");
         }
         if((classOneIndices.length % 2) != 0) {
            AnalysisUtil.exit(
                  "The number of items in class 1 must be an even number for balanced permutations.");
         }
      }
      scores = new double[N];
		double[] classZeroMeans = new double[N];
		double[] classOneMeans = new double[N];
		double[] classZeroStds = new double[N];
		double[] classOneStds = new double[N];
	
		for(int i = 0; i < N; i++) {
			double classZeroMean = Util.mean(dataset, classZeroIndices, i);
			double classZeroStd = Util.standardDeviation(dataset, classZeroIndices,
         	i, classZeroMean);
		
			double classOneMean = Util.mean(dataset, classOneIndices, i);
			double classOneStd = Util.standardDeviation(dataset, classOneIndices,
         	i, classOneMean);
			classZeroMeans[i] = classZeroMean;
			classOneMeans[i] = classOneMean;
			classZeroStds[i] = classZeroStd;
			classOneStds[i] = classOneStd;
		}
      statisticalMeasure.compute(dataset,
            classZeroIndices,
            classOneIndices, scores);

      Permuter permuter = null;
      if(covariate != null) {
         permuter = new UnbalancedRandomCovariatePermuter(classVector, covariate, seed);
         if(complete || balanced) {
            AnalysisUtil.exit("Covariate permutations not yet implemented for complete or balanced permutations.");
         }
      } else if(!complete && balanced) {
         permuter = new BalancedRandomPermuter(classZeroIndices, classOneIndices, seed);
         seedUsed = true;
      } else if(!complete && !balanced) {
         permuter = new UnbalancedRandomPermuter(classVector.size(), classVector.getIndices(1).length, seed);
         seedUsed = true;
      } else if(complete && !balanced) {
         permuter = new UnbalancedCompletePermuter(classVector.size(),
               classZeroIndices.length);
         java.math.BigInteger totalPermutations = ((UnbalancedCompletePermuter) (permuter)).getTotal();
         if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
            AnalysisUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
         }
         numPermutations = totalPermutations.intValue();
      } else if(complete && balanced) {
         permuter = new BalancedCompletePermuter(classZeroIndices, classOneIndices);
         
         java.math.BigInteger totalPermutations = ((BalancedCompletePermuter) (permuter)).getTotal();
         if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
            AnalysisUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
         }
      }

      int[] descendingIndices = Sorting.index(scores, Sorting.DESCENDING);
      double[] absoluteScores = (double[]) scores.clone();
      for(int i = 0; i < N; i++) {
         absoluteScores[i] = Math.abs(absoluteScores[i]);
      }
      int[] descendingAbsIndices = Sorting.index(absoluteScores, Sorting.DESCENDING);

      double[] fwer = new double[N];
      double[] fpr = new double[N];
      double[] featureSpecificPValues = new double[N];
      double[] permutedScores = new double[N];
		double[] maxT = new double[N];
			
      for(int perm = 0; perm < numPermutations; perm++) {
         int[] permutedClassZeroIndices = null;
         int[] permutedClassOneIndices = null;
         int[] permutedAssignments = null;
         if(readPermutations) {
            permutedAssignments = nextPermutation();
         } else {
            permutedAssignments = permuter.next();
         }
        
         if(debugger != null) {
            debugger.addAssignment(permutedAssignments);
         }
         java.util.List zeroIndices = new java.util.ArrayList();
         java.util.List oneIndices = new java.util.ArrayList();
         for(int i = 0, length = permutedAssignments.length; i < length; i++) {
            if(permutedAssignments[i] == 0) {
               zeroIndices.add(new Integer(i));
            } else {
               oneIndices.add(new Integer(i));
            }
         }
         permutedClassZeroIndices = new int[zeroIndices.size()];
         for(int i = 0, length = permutedClassZeroIndices.length; i < length; i++) {
            permutedClassZeroIndices[i] = ((Integer) zeroIndices.get(i)).intValue();
         }

         permutedClassOneIndices = new int[oneIndices.size()];
         for(int i = 0, length = permutedClassOneIndices.length; i < length; i++) {
            permutedClassOneIndices[i] = ((Integer) oneIndices.get(i)).intValue();
         }

         statisticalMeasure.compute(dataset,
               permutedClassZeroIndices,
               permutedClassOneIndices, permutedScores);// compute scores using permuted class labels

         for(int i = 0; i < N; i++) {
            double score = scores[i];
            if(testDirection == TWO_SIDED || testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
               if(permutedScores[i] >= score) {
                  featureSpecificPValues[i] += 1.0;
               }
            } else {
               if(permutedScores[i] <= score) {
                  featureSpecificPValues[i] += 1.0;
               }
            }
         }

			for(int i = permutedScores.length-2; i >= 0; i--) {
				double score_i = permutedScores[i];
				double score_i_plus_1 = permutedScores[i+1];
				permutedScores[i] = Math.max(score_i, score_i_plus_1);	
			}
			
			for(int i = descendingAbsIndices.length-1; i >= 0; i--) {
				int index = descendingAbsIndices[i];
				double score = absoluteScores[index];
				if(permutedScores[i] >= score) {
					maxT[index] += 1.0;
            }
			}	
				
         for(int i = 0; i < N; i++) {
            permutedScores[i] = Math.abs(permutedScores[i]);
         }
         Sorting.sort(permutedScores, Sorting.DESCENDING);

         int j = 0;
         int count = 0;
         for(int i = 0; i < N; i++) {
            double score = absoluteScores[descendingAbsIndices[i]];
            while(j < N && score < permutedScores[j]) {
               count++;
               j++;
            }
            fpr[descendingAbsIndices[i]] += count;

            if(count > 0) {
               fwer[descendingAbsIndices[i]]++;
            }
         }
      }

      for(int i = 0; i < N; i++) {
         fpr[i] /= N;
         fpr[i] /= numPermutations;
         featureSpecificPValues[i] /= numPermutations;
			maxT[i] /= numPermutations;
         if(testDirection == TWO_SIDED) {
            featureSpecificPValues[i] = 2.0 * Math.min(featureSpecificPValues[i], 1.0 - featureSpecificPValues[i]);
         }
         fwer[i] /= numPermutations;
      }

      double[] fdr = new double[N];
      int[] pValueIndices = Sorting.index(featureSpecificPValues, Sorting.ASCENDING);
      int[] ranks = Sorting.rank(pValueIndices);

      // check for ties
      for(int i = pValueIndices.length - 1; i > 0; i--) {
         double bigPValue = featureSpecificPValues[pValueIndices[i]];
         double smallPValue = featureSpecificPValues[pValueIndices[i - 1]];
         if(bigPValue == smallPValue) {
            ranks[pValueIndices[i - 1]] = ranks[pValueIndices[i]];
         }
      }

      for(int i = 0; i < N; i++) {
         int rank = ranks[i];
         double p = featureSpecificPValues[i];
         fdr[i] = (p * N) / rank;
      }
		
		int[] featurePIndices = Sorting.index(featureSpecificPValues, Sorting.ASCENDING);

      // ensure fdr is monotonically decreasing
      for(int i = 1; i < featurePIndices.length; i++) { 
			int index = featurePIndices[i];
			int indexMinusOne = featurePIndices[i-1];
			fdr[index] = Math.min(fdr[index], fdr[indexMinusOne]);
		}
        
      int[] _ranks = null;

      if(testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE || testDirection == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
         _ranks = Sorting.rank(descendingIndices);
      } else if(testDirection == TWO_SIDED) {
         _ranks = Sorting.rank(descendingAbsIndices);
      }

      // write p values to temporary file
      PrintWriter tempFileWriter = null;
      String tempFileName = "pvalues.txt";
      new File(tempFileName).deleteOnExit();
      try {
         tempFileWriter = new PrintWriter(new FileWriter(tempFileName));
         for(int i = 0; i < N; i++) {
            tempFileWriter.println(featureSpecificPValues[descendingIndices[i]]);
         }
      } catch(IOException ioe) {
         AnalysisUtil.exit("An error occurred while saving the output file.");
      } finally {
         if(tempFileWriter != null) {
            tempFileWriter.close();
         }
      }

      edu.mit.broad.marker.qvalue.QValue.qvalue(tempFileName);
      String[] qvalues = new String[N];
      int qvalueHeaderLines = 4;
      String[] qvalueHeaders = new String[qvalueHeaderLines];
      BufferedReader br = null;
      String qvalueOutputFileName = "qvalues.txt";
      new File(qvalueOutputFileName).deleteOnExit();
      try {// read in results from qvalue
         br = new BufferedReader(new FileReader(qvalueOutputFileName));
         for(int i = 0; i < qvalueHeaderLines; i++) {
            qvalueHeaders[i] = br.readLine();
         }
         for(int i = 0; i < N; i++) {
            qvalues[i] = br.readLine();
         }
      } catch(IOException ioe) {
         AnalysisUtil.exit("An error occurred while saving the output file.");
      } finally {
         if(br != null) {
            try {
               br.close();
            } catch(IOException x) {
            }
         }
      }

      PrintWriter pw = null;

      try {
         if(!outputFileName.toLowerCase().endsWith(".odf")) {
            outputFileName += ".odf";
         }
         pw = new PrintWriter(new FileWriter(outputFileName));
         pw.println("ODF 1.0");
         String numHeaderLines = seedUsed?"19":"18";
         pw.println("HeaderLines="+numHeaderLines);
         pw.println("COLUMN_NAMES:Rank\tFeature\tScore\tFeature Specific P Value\tFPR\tFWER\tmaxT P Value Value\tFDR(BH)\tBonferroni\tQ Value\tClass Zero Mean\tClass One Mean\tClass Zero Std\tClass One Std");
         pw.println("COLUMN_TYPES:int\tString\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat\tfloat");
         pw.println("Model=Comparative Marker Selection");
         pw.println("Dataset File=" + AnalysisUtil.getFileName(datasetFile));
         pw.println("Class File=" + AnalysisUtil.getFileName(clsFile));
         pw.println("Permutations=" + numPermutations);
         pw.println("Balanced=" + balanced);
         pw.println("Complete=" + complete);
         if(testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
            pw.println("Test Direction=Class 0");
         } else if(testDirection == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
            pw.println("Test Direction=Class 1");
         } else {
            pw.println("Test Direction=2 Sided");
         }
         pw.println("Class 0=" + classVector.getClassName(0));
         pw.println("Class 1=" + classVector.getClassName(1));
         if(metric == SNR) {
            pw.println("Test Statistic=SNR");
         } else if(metric == SNR_MEDIAN) {
            pw.println("Test Statistic=SNR (median)");
         } else if(metric == T_TEST) {
            pw.println("Test Statistic=T-Test");
         } else if(metric == T_TEST_MEDIAN) {
            pw.println("Test Statistic=T-Test (median)");
         } else if(metric == T_TEST_MIN_STD) {
            pw.println("Test Statistic=T-Test (min std=" + minStd + ")");
         }

         pw.println("Fix Standard Deviation=" + fixStdev);
         for(int i = 0; i < qvalueHeaderLines; i++) {
            pw.println(qvalueHeaders[i]);
         }
         if(seedUsed) {
            pw.println("Random Seed=" + seed);  
         }
         pw.println("DataLines=" + N);
         for(int i = 0; i < N; i++) {
            int index = descendingIndices[i];

            int rank = _ranks[index];
            if(testDirection == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
               rank = N - rank + 1;
            }
				
            double bonferroni = Math.min(featureSpecificPValues[index] * N, 1.0);
				pw.print(rank);
				pw.print("\t");
				pw.print(dataset.getRowName(index));
				pw.print("\t");
				pw.print(scores[index]);
				pw.print("\t");
				pw.print(featureSpecificPValues[index]);
				pw.print("\t");
				pw.print(fpr[index]);
				pw.print("\t");
				pw.print(fwer[index]);
				pw.print("\t");
				pw.print(maxT[index]);
				pw.print("\t");
				pw.print(fdr[index]);
				pw.print("\t");
				pw.print(bonferroni);
				pw.print("\t");
				pw.print(qvalues[i]);
				pw.print("\t");
				pw.print(classZeroMeans[index]);
				pw.print("\t");
				pw.print(classOneMeans[index]);
				pw.print("\t");
				pw.print(classZeroStds[index]);
				pw.print("\t");
				pw.print(classOneStds[index]);
				pw.println();
         }

      } catch(Exception e) {
         AnalysisUtil.exit("An error occurred while saving the output file.");
      } finally {

         if(pw != null) {

            try {
               pw.close();
            } catch(Exception x) {
            }
         }
      }

      if(debugger != null) {
         debugger.print();
      }
   }



   /**
    *  Reads the next permutation from a file
    *
    * @return    the next permutation
    */
   int[] nextPermutation() {
      try {
         String s = testClassPermutationsReader.readLine();
         StringTokenizer st = new StringTokenizer(s, " \t");
         int cols = dataset.getColumnCount();
         int[] a = new int[cols];

         for(int j = 0; j < cols; j++) {
            a[j] = Integer.parseInt(st.nextToken());
         }

         return a;
      } catch(Exception e) {
         AnalysisUtil.exit("An error occurred while reading the permutations file.");
         return null;
      }

   }


   /**
    *  Opens the permutations file for reading
    *
    * @param  fileName  the file name
    */
   void initFile(String fileName) {

      try {
         testClassPermutationsReader = new BufferedReader(new FileReader(
               fileName));
      } catch(Exception e) {
         AnalysisUtil.exit(
               "An error occurred while reading the permutations file.");
      }
   }



   static class Debugger {
      java.util.Map assignment2Occurences = new java.util.HashMap();


      public static String toString(int[] a) {
         StringBuffer buf = new StringBuffer();
         for(int i = 0; i < a.length; i++) {
            buf.append(a[i]);
         }
         return buf.toString();
      }


      public void addAssignment(int[] a) {
         String s = toString(a);
         Integer occurences = (Integer) assignment2Occurences.get(s);
         if(occurences == null) {
            occurences = new Integer(0);
         }
         assignment2Occurences.put(s, new Integer(occurences.intValue() + 1));
      }


      public void print() {
         System.out.println("Unique permutations " + assignment2Occurences.size());
         System.out.println("permutations");
         System.out.println(assignment2Occurences);
      }
   }

   static {
      if("true".equalsIgnoreCase(System.getProperty("edu.mit.broad.marker.debug"))) {
         debugger = new Debugger();
      }
   }

}

