package edu.mit.broad.marker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.StringTokenizer;
import java.math.BigInteger;

import org.genepattern.data.expr.ExpressionData;
import org.genepattern.data.matrix.ClassVector;
import org.genepattern.data.matrix.DoubleMatrix2D;
import org.genepattern.io.expr.ExpressionDataCreator;
import org.genepattern.io.expr.IExpressionDataReader;
import org.genepattern.module.AnalysisUtil;
import org.genepattern.stats.Function;
import org.genepattern.stats.ITestStatistic;
import org.genepattern.stats.Sorting;
import org.genepattern.stats.TestStatistics;
import org.genepattern.stats.ZeroFinder;
import org.genepattern.ioutil.Util;

import edu.mit.broad.marker.permutation.BalancedCompletePermuter;
import edu.mit.broad.marker.permutation.BalancedRandomPermuter;
import edu.mit.broad.marker.permutation.Permuter;
import edu.mit.broad.marker.permutation.UnbalancedCompletePermuter;
import edu.mit.broad.marker.permutation.UnbalancedRandomCovariatePermuter;
import edu.mit.broad.marker.permutation.UnbalancedRandomPermuter;

/**
 * @author     Joshua Gould
 */
public class MarkerSelection {
	private static boolean debug;
	private static Debugger debugger;

	/** input expression values */
	private DoubleMatrix2D dataset;

	/** data array stored in <tt>dataset</tt> */
	private double[][] dataArray;
	
	/** input classes */
	private ClassVector classVector;

	/** total number of permutations to perform */
	private int numPermutations;

	/** whether to do one-sided >, one-sided <, or two-sided test */
	private int testDirection;

	private int metric;

	/**
	 * minimum standard deviation when metric==T_TEST_MIN_STD or SNR_MIN_STD or
	 * SNR_MEDIAN_MIN_STD
	 */
	private double minStd;

	/** whether permutations are balanced */
	private boolean balanced;

	/** whether to compute all possible permutations */
	private boolean complete;

	/** whether to read permutations from a file */
	private boolean readPermutations = false;

	/** reader when reading permutations from file */
	private BufferedReader testClassPermutationsReader;

	/** name of output file */
	private String outputFileName;

	/** unpermuted scores */
	private double[] scores;

	/** number of rows in input data */
	private int numFeatures;

	private ITestStatistic statisticalMeasure;

	/** Input dataset file name */
	private String datasetFile;

	/** Input cls file name */
	private String clsFile;
	
	/** Input confounding variable cls file name */
	private String confoundingClsFile;

	/** Covariate assignments */
	private ClassVector covariate;

	/** Seed for generating permutations */
	private int seed;

	/** Whether a seed was used for generating permutations */
	private boolean seedUsed = false;

	/** Family-wise error rate, size=numFeatures */
	private double[] fwer;

	/** size=numFeatures */
	private double[] featureSpecificPValues;

	/** Temporary array used in calculation of permuted scores */
	private double[] permutedScores;

	/** Temporary array used in calculation of permuted score */
	private double[] monotonicPermutedScores;

	/** size=numFeatures */
	private double[] maxT;

	/** Absolute value of the test statistic, size=numFeatures */
	private double[] absoluteScores;

	/** descending sorted indices into absoluteScores */
	private int[] descendingAbsIndices;

	private Permuter permuter;

	/** Lower confidence interval for feature p value */
	private double[] lowerBound;
	
	/** Upper confidence interval for feature p value */
	private double[] upperBound;
	
	/** Whether to try to trim the number of features to permute */
	private boolean removeFeatures;

	/** Keeps track of number of permutations performed per feature */
	private int[] permutationsPerFeature;
	
	/** Threshold for removing features when <tt>removeFeatures</tt> is <tt>true</tt> */
	private double theta = 0.25;
	
	/** Whether to smooth p values */
	private boolean smoothPValues;

	public MarkerSelection(DoubleMatrix2D dataset, String datasetFile, 
			ClassVector classVector, String clsFile,
			int _numPermutations, int _side, String _outputFileName,
			boolean _balanced, boolean complete, int metric,
			double minStd, int seed, ClassVector confoundingClassVector,
			String confoundingClsFile, boolean removeFeatures, double theta, boolean smoothPValues) {
		this.dataset = dataset;
		this.dataArray = dataset.getArray();
		this.datasetFile = datasetFile;
		this.classVector = AnalysisUtil.readClassVector(clsFile);
		AnalysisUtil.checkDimensions(dataset, classVector);
		if (classVector.getClassCount() != 2) {
			AnalysisUtil.exit("Class file must contain 2 classes.");
		}
		this.clsFile = clsFile;
		this.numPermutations = _numPermutations;
		this.testDirection = _side;
		this.outputFileName = _outputFileName;
		this.balanced = _balanced;
		this.complete = complete;
		this.metric = metric;
		this.minStd = minStd;
		this.seed = seed;
		this.covariate = confoundingClassVector;
		this.confoundingClsFile = confoundingClsFile;
		this.removeFeatures = removeFeatures;
		this.theta = theta;
		this.smoothPValues = smoothPValues;
		if(complete && removeFeatures) {
			AnalysisUtil.exit("Speedup option can only be used when performing random permutations.");
		}
		
		this.numFeatures = dataset.getRowCount();
		
		String permutationsFile = System
				.getProperty("edu.mit.broad.marker.perm");
		if (permutationsFile != null) {
			System.out.println("reading permutations from file "
					+ permutationsFile);
			readPermutations = true;
			initFile(permutationsFile);
		}
		if (metric == Constants.T_TEST) {
			statisticalMeasure = new TestStatistics.TTest();
		} else if (metric == Constants.T_TEST_MEDIAN) {
			statisticalMeasure = new TestStatistics.TTestMedian();
		} else if (metric == Constants.T_TEST_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.TTestMinStd(minStd);
		}  else if (metric == Constants.T_TEST_MEDIAN_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.TTestMedianMinStd(minStd);
					
		} else if (metric == Constants.SNR) {
			statisticalMeasure = new TestStatistics.SNR();
		} else if (metric == Constants.SNR_MEDIAN) {
			statisticalMeasure = new TestStatistics.SNRMedian();
		
		} else if (metric == Constants.SNR_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.SNRMinStd(minStd);
		} else if (metric == Constants.SNR_MEDIAN_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.SNRMedianMinStd(minStd);
			
		} else {
			AnalysisUtil.exit("Unknown test statistic");
		}

		if (testDirection != Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE
				&& testDirection != Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE
				&& testDirection != Constants.TWO_SIDED) {
			AnalysisUtil.exit("Unknown test direction.");
		}

		computePValues();

		if (readPermutations) {

			if (testClassPermutationsReader != null) {

				try {
					testClassPermutationsReader.close();
				} catch (Exception e) {
				}
			}
		}
	}

	public static void main(String[] args) {
		debug = true;
		try {
			run(args);
		} catch (Throwable t) {
			if(debug) {
				t.printStackTrace();
			}
			AnalysisUtil.exit("An error occurred while running the algorithm.");
		}
	}

	static void run(String[] args) {
		String datasetFile = args[0];
		String clsFile = args[1];
		int _numPermutations = 0;
		try {
			_numPermutations = Integer.parseInt(args[2]);
		} catch (NumberFormatException nfe) {
			AnalysisUtil.exit("Number of permutations is not an integer.");
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
		} catch (NumberFormatException nfe) {
			AnalysisUtil.exit("Random seed is not an integer.");
		}
		boolean removeFeatures = Boolean.valueOf(args[10]).booleanValue();
		boolean smoothPValues = Boolean.valueOf(args[11]).booleanValue();
		
		double minStd = -1;
		String confoundingClsFile = null;
		double theta = -1;
		for (int i = 12; i < args.length; i++) {
			String arg = args[i].substring(0, 2);
			String value = args[i].substring(2, args[i].length());
			if (value.equals("")) {
				continue;
			}
			if (arg.equals("-m")) {
				try {
					minStd = Double.parseDouble(value);
				} catch (NumberFormatException nfe) {
					AnalysisUtil
							.exit("Minimum standard deviation is not a number.");
				}
			} else if (arg.equals("-c")) {
				confoundingClsFile = value;
			} else if(arg.equals("-t")) {
				try {
					theta = Double.parseDouble(value);
				} catch (NumberFormatException nfe) {
					AnalysisUtil
							.exit("Theta is not a number.");
				}
			}
		}

		ClassVector classVector = AnalysisUtil.readClassVector(clsFile);
		IExpressionDataReader reader = AnalysisUtil
				.getExpressionReader(datasetFile);
		ExpressionData expressionData = AnalysisUtil
				.readExpressionData(reader, datasetFile);

		DoubleMatrix2D dataset = expressionData.getExpressionMatrix();
		
		ClassVector confoundingClassVector = null;
		if (confoundingClsFile != null) {
			confoundingClassVector = AnalysisUtil.readClassVector(confoundingClsFile);
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


		if(classVector.getClassCount() > 2) {
			ClassVector[] oneVersusAll = classVector.getOneVersusAll();
			String baseOutputFileName = Util.getBaseFileName(outputFileName);
			for(int i = 0; i < classVector.getClassCount(); i++) {
				String tempOutputFileName = baseOutputFileName + "." + 
					classVector.getClassName(i) + ".vs.All.odf";
				new MarkerSelection(dataset, datasetFile, oneVersusAll[i], 
				clsFile, _numPermutations, testDirection, tempOutputFileName, 
				balanced, complete, metric, minStd, seed, 
				confoundingClassVector, confoundingClsFile, removeFeatures, theta, smoothPValues);
			}
		} else {
			new MarkerSelection(dataset, datasetFile, classVector, clsFile,
				_numPermutations, testDirection, outputFileName, balanced, 
				complete, metric, minStd, seed, 
				confoundingClassVector, confoundingClsFile, removeFeatures, theta, smoothPValues);
		}
	}

	final void computePValues() {
		int[] classZeroIndices = classVector.getIndices(0);
		int[] classOneIndices = classVector.getIndices(1);
		if (balanced) {
			if (classZeroIndices.length != classOneIndices.length) {
				AnalysisUtil
						.exit("The number of items in each class must be equal for balanced permutations.");
			}
			if ((classZeroIndices.length % 2) != 0) {
				AnalysisUtil
						.exit("The number of items in class 0 must be an even number for balanced permutations.");
			}
			if ((classOneIndices.length % 2) != 0) {
				AnalysisUtil
						.exit("The number of items in class 1 must be an even number for balanced permutations.");
			}
		}
		scores = new double[numFeatures];
		lowerBound = new double[numFeatures];
		upperBound = new double[numFeatures];
		
		statisticalMeasure.compute(dataArray, classZeroIndices, classOneIndices,
				scores, null, -1); 
		
		// index of observed scores, greatest to smallest
		int[] descendingIndices = Sorting.index(scores, Sorting.DESCENDING);
		
		if (covariate != null) {
			seedUsed = true;
			permuter = new UnbalancedRandomCovariatePermuter(classVector,
					covariate, seed);
			if (complete || balanced) {
				AnalysisUtil
						.exit("Covariate permutations not yet implemented for complete or balanced permutations.");
			}
		} else if (!complete && balanced) {
			seedUsed = true;
			permuter = new BalancedRandomPermuter(classZeroIndices,
					classOneIndices, seed);

		} else if (!complete && !balanced) {
			seedUsed = true;
			permuter = new UnbalancedRandomPermuter(classVector.size(),
					classVector.getIndices(1).length, seed);
		} else if (complete && !balanced) {
			seedUsed = false;
			permuter = new UnbalancedCompletePermuter(classVector.size(),
					classZeroIndices.length);
			java.math.BigInteger totalPermutations = ((UnbalancedCompletePermuter) (permuter))
					.getTotal();
			if ((totalPermutations.compareTo(new java.math.BigInteger(""
					+ Integer.MAX_VALUE))) == 1) {
				AnalysisUtil.exit("Number of permutations exceeds maximum of "
						+ Integer.MAX_VALUE + ".");
			}
			numPermutations = totalPermutations.intValue();
		} else if (complete && balanced) {
			seedUsed = false;
			permuter = new BalancedCompletePermuter(classZeroIndices,
					classOneIndices);

			java.math.BigInteger totalPermutations = ((BalancedCompletePermuter) (permuter))
					.getTotal();
			if ((totalPermutations.compareTo(new java.math.BigInteger(""
					+ Integer.MAX_VALUE))) == 1) {
				AnalysisUtil.exit("Number of permutations exceeds maximum of "
						+ Integer.MAX_VALUE + ".");
			}
			numPermutations = totalPermutations.intValue();
		}

		featureSpecificPValues = new double[numFeatures];
		permutedScores = new double[numFeatures];
		absoluteScores = (double[]) scores.clone();
		for (int i = 0; i < numFeatures; i++) {
				absoluteScores[i] = Math.abs(absoluteScores[i]);
		}
		descendingAbsIndices = Sorting
					.index(absoluteScores, Sorting.DESCENDING);
					
		if (!removeFeatures) {
			monotonicPermutedScores = new double[numFeatures];
			maxT = new double[numFeatures];
			fwer = new double[numFeatures];
		}
			
		long startTime = System.currentTimeMillis();
		permute();
		long endTime = System.currentTimeMillis();
		long elapsed = endTime-startTime;
		System.out.println("elapsed time " + elapsed/1000.0 + " seconds");
			
		// free temporary storage
		monotonicPermutedScores = null;
		absoluteScores = null;
		permutedScores = null; 
		
		for (int i = 0; i < numFeatures; i++) {
			// fpr[i] /= numFeatures;
			// fpr[i] /= numPermutations;
			int N = permutationsPerFeature[i];
			double k = featureSpecificPValues[i];
			double p;
			if (smoothPValues) {
				p = (k + 1)/(N + 2);
			} else {
				p = k/N;
			}
			if (testDirection == Constants.TWO_SIDED) {
				double oneMinusP = 1.0 - p;
				if(oneMinusP < p) {
					 p = oneMinusP;
				}
				p *= 2;
				
				if(p==0) {
					// ensure not degenerate case where profile is completely flat
					// TODO handle cases where profile is flat (but not completely)
					double[] row = dataArray[i];
					double val = row[0];
					boolean flat = true;
					for(int j = 1, cols = row.length; j < cols && flat; j++) {
						if(row[j]!=val) {
							flat = false;
						}
					}
					if(flat) {
						p = 1;
					}
				} 
			}
			featureSpecificPValues[i] = p;
			
	
			if (testDirection == Constants.TWO_SIDED) {
				 k = 2*Math.min(k,N-k); 
			}
			
			double shape1 = k + 1;
			double shape2 = N - k + 1;
			final JSci.maths.statistics.BetaDistribution betaDist = 
				new JSci.maths.statistics.BetaDistribution(shape1, shape2); 
			double plow = betaDist.inverse(0.025);
			double phigh = betaDist.inverse(0.975);
			
			Function function = new Function() {
				public double evaluate(double x) {
					double d = Math.min(1, 0.95+betaDist.cumulative(x));
					return  betaDist.probability(x)-
						betaDist.probability(
							betaDist.inverse(d));
				}
			};
			
			if(k==0) {
				plow = 0;
				phigh = betaDist.inverse(0.05);
			} else if(k==N) {
				plow = betaDist.inverse(0.95);
				phigh = 1;
			} else if(k < (N/2)) {
				double accuracy = Math.min(p/1000, 0.000001);
				// search between 0 and plow
				plow = ZeroFinder.bisection(0, plow, accuracy,
         function);
				phigh = betaDist.inverse(0.95+betaDist.cumulative(plow));
			} else if(k > (N/2)) {
				double accuracy = Math.min(p/1000, 0.000001);
				 // search between plow and beta.inverse(0.05)
				plow = ZeroFinder.bisection(plow, betaDist.inverse(0.05), accuracy,
         function);
				phigh = betaDist.inverse(0.95+betaDist.cumulative(plow));
			}
			
			
			lowerBound[i] = plow;
			upperBound[i] = phigh;
			
			if (!removeFeatures) {
				fwer[i] /= N;
				// rankBasedPValues[i] /= N;
				maxT[i] /= N;
			}
		}
		

		double[] fdr = new double[numFeatures];
		int[] pValueIndices = Sorting.index(featureSpecificPValues,
				Sorting.ASCENDING);
		int[] ranks = Sorting.rank(pValueIndices);

		// check for ties
		for (int i = pValueIndices.length - 1; i > 0; i--) {
			double bigPValue = featureSpecificPValues[pValueIndices[i]];
			double smallPValue = featureSpecificPValues[pValueIndices[i - 1]];
			if (bigPValue == smallPValue) {
				ranks[pValueIndices[i - 1]] = ranks[pValueIndices[i]];
			}
		}

		for (int i = 0; i < numFeatures; i++) {
			int rank = ranks[i];
			double p = featureSpecificPValues[i];
			fdr[i] = (p * numFeatures) / rank;
		}

		// ensure fdr is monotonically decreasing
		int[] pIndices = Sorting.index(featureSpecificPValues,
				Sorting.DESCENDING);
		for (int i = 0; i < pIndices.length - 1; i++) {
			int highIndex = pIndices[i];
			int lowIndex = pIndices[i + 1];
			fdr[lowIndex] = Math.min(fdr[lowIndex], fdr[highIndex]);
		}

		for (int i = 0; i < numFeatures; i++) {
			fdr[i] = Math.min(fdr[i], 1);
		}

		if (!removeFeatures) {
			for (int i = 1; i < numFeatures; i++) { // ensure maxT is monotonic
				maxT[descendingAbsIndices[i]] = Math.max(
						maxT[descendingAbsIndices[i]],
						maxT[descendingAbsIndices[i - 1]]);
			}
		}
		int[] _ranks = null;

		if (testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE
				|| testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
			_ranks = Sorting.rank(descendingIndices);
		} else if (testDirection == Constants.TWO_SIDED) {
			_ranks = Sorting.rank(descendingAbsIndices);
		}

		// write p values to temporary file
		PrintWriter tempFileWriter = null;
		String tempFileName = "pvalues.txt";
		new File(tempFileName).deleteOnExit();
		try {
			tempFileWriter = new PrintWriter(new FileWriter(tempFileName));
			for (int i = 0; i < numFeatures; i++) {
				tempFileWriter
						.println(featureSpecificPValues[descendingIndices[i]]);
			}
		} catch (IOException ioe) {
			if(debug) {
				ioe.printStackTrace();
			}
			AnalysisUtil.exit(
					"An error occurred while saving the output file.");
		} finally {
			if (tempFileWriter != null) {
				tempFileWriter.close();
			}
		}

		edu.mit.broad.marker.qvalue.QValue.qvalue(tempFileName);
		String[] qvalues = new String[numFeatures];
		int qvalueHeaderLines = 4;
		String[] qvalueHeaders = new String[qvalueHeaderLines];
		BufferedReader br = null;
		String qvalueOutputFileName = "qvalues.txt"; // produced by R code above
		new File(qvalueOutputFileName).deleteOnExit();
		
		try {// read in results from qvalue
			br = new BufferedReader(new FileReader(qvalueOutputFileName));
			for (int i = 0; i < qvalueHeaderLines; i++) {
				qvalueHeaders[i] = br.readLine();
			}
			for (int i = 0; i < numFeatures; i++) {
				qvalues[i] = br.readLine();
			}
		} catch (IOException ioe) {
			if(debug) {
				ioe.printStackTrace();
			}
			if(!debug) {
				AnalysisUtil.exit(
					"An error occurred while saving the output file.");
			}
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException x) {
				}
			}
		}

		PrintWriter pw = null;

		try {
			if (!outputFileName.toLowerCase().endsWith(".odf")) {
				outputFileName += ".odf";
			}
			pw = new PrintWriter(new FileWriter(outputFileName));
			pw.println("ODF 1.0");
			int numHeaderLines = 18;
			if(seedUsed) {
				numHeaderLines++;
			}
			if(covariate!=null) {
				numHeaderLines++;
			}
			
			pw.println("HeaderLines=" + numHeaderLines);
			String[] columnNames = null;
			if(!removeFeatures) {
					columnNames = new String[]{"Rank", "Feature", "Score", "Feature P", "Feature P Low", "Feature P High", "FWER", "FDR(BH)", "Bonferroni", "Q Value", "maxT"};
			
			} else {
				columnNames = new String[]{"Rank", "Feature", "Score", "Feature P", "Feature P Low", "Feature P High", "FDR(BH)", "Bonferroni", "Q Value", "Permutations"};
			}
			pw.print("COLUMN_NAMES:");
			for(int j = 0; j < columnNames.length; j++) {
				if(j > 0) {
					pw.print("\t");
				}
				pw.print(columnNames[j]);
			}
			pw.println();
			pw.println("Model=Comparative Marker Selection");
			pw.println("Dataset File=" + AnalysisUtil.getFileName(datasetFile));
			pw.println("Class File=" + AnalysisUtil.getFileName(clsFile));
			if(covariate!=null) {
				pw.println("Confounding Class File=" + 
					AnalysisUtil.getFileName(confoundingClsFile));
			}
			if(!removeFeatures) {
				pw.println("Permutations=" + numPermutations);
			}
			pw.println("Balanced=" + balanced);
			pw.println("Complete=" + complete);
			if (testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
				pw.println("Test Direction=Class 0");
			} else if (testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
				pw.println("Test Direction=Class 1");
			} else {
				pw.println("Test Direction=2 Sided");
			}
			pw.println("Class 0=" + classVector.getClassName(0));
			pw.println("Class 1=" + classVector.getClassName(1));
			pw.println("Test Statistic=" + statisticalMeasure.toString());
			
			for (int i = 0; i < qvalueHeaderLines; i++) {
				pw.println(qvalueHeaders[i]);
			}
			if (seedUsed) {
				pw.println("Random Seed=" + seed);
			}
			
			pw.println("Remove Features=" + removeFeatures);
			if(removeFeatures) {
				pw.println("Theta=" + theta);
			}
			pw.println("Smooth p-values=" + smoothPValues);
			pw.println("DataLines=" + numFeatures);
	
			for (int i = 0; i < numFeatures; i++) {
				int index = descendingIndices[i];

				int rank = _ranks[index];
				if (testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
					rank = numFeatures - rank + 1;
				}
				double bonferroni = Math.min(featureSpecificPValues[index] * numFeatures,
						1.0);
				pw.print(rank);
				pw.print("\t");
				pw.print(dataset.getRowName(index));
				pw.print("\t");
				pw.print(scores[index]);
				pw.print("\t");
				pw.print(featureSpecificPValues[index]);
				pw.print("\t");
				pw.print(lowerBound[index]);
				pw.print("\t");
				pw.print(upperBound[index]);
				pw.print("\t");
				
				if(!removeFeatures) {
					pw.print(fwer[index]);
					pw.print("\t");
				} 
				
				pw.print(fdr[index]);
				pw.print("\t");
				pw.print(bonferroni);
				pw.print("\t");
				pw.print(qvalues[i]);
				if(!removeFeatures) {
					pw.print("\t");
					pw.print(maxT[index]);
				} 
				if(removeFeatures) {
					pw.print("\t");
					pw.print(permutationsPerFeature[index]);
				}
				pw.println();
			}

		} catch (Exception e) {
			if(debug) {
				e.printStackTrace();
			}
			AnalysisUtil.exit(
					"An error occurred while saving the output file.");
		} finally {

			if (pw != null) {

				try {
					pw.close();
				} catch (Exception x) {
				}
			}
		}

		if (debugger != null) {
			debugger.print();
		}
	}

	private final void permute() {
		BigInteger maxLong = BigInteger.valueOf(Long.MAX_VALUE);
		BigInteger bigIntBank = BigInteger.valueOf(numFeatures).multiply(
			BigInteger.valueOf(numPermutations));
		if(bigIntBank.compareTo(maxLong) > 0) {
			AnalysisUtil.exit(
				"Arithmetic overflow. Try decreasing the number of permutations.");
		}
		long bank = bigIntBank.longValue();
		int[] featureIndicesToPermute = null;
		permutationsPerFeature = new int[numFeatures];
		if(removeFeatures) {
			featureIndicesToPermute = new int[numFeatures];
			for(int i = 0; i < numFeatures; i++) {
				featureIndicesToPermute[i] = i;
			}
		} else {
			Arrays.fill(permutationsPerFeature, numPermutations);
		}
	
		int lengthOfIndicesToPermute = numFeatures;
		
		while(bank > 0) {
			int[] permutedClassZeroIndices = null;
			int[] permutedClassOneIndices = null;
			int[] permutedAssignments = null;
			if (readPermutations) {
				permutedAssignments = nextPermutation();
			} else {
				permutedAssignments = permuter.next();
			}

			if (debugger != null) {
				debugger.addAssignment(permutedAssignments);
			}
			java.util.List zeroIndices = new java.util.ArrayList();
			java.util.List oneIndices = new java.util.ArrayList();
			for (int i = 0, length = permutedAssignments.length; i < length; i++) {
				if (permutedAssignments[i] == 0) {
					zeroIndices.add(new Integer(i));
				} else {
					oneIndices.add(new Integer(i));
				}
			}
			permutedClassZeroIndices = new int[zeroIndices.size()];
			for (int i = 0, length = permutedClassZeroIndices.length; i < length; i++) {
				permutedClassZeroIndices[i] = ((Integer) zeroIndices.get(i))
						.intValue();
			}

			permutedClassOneIndices = new int[oneIndices.size()];
			for (int i = 0, length = permutedClassOneIndices.length; i < length; i++) {
				permutedClassOneIndices[i] = ((Integer) oneIndices.get(i))
						.intValue();
			}

			statisticalMeasure.compute(dataArray, permutedClassZeroIndices,
					permutedClassOneIndices, permutedScores, featureIndicesToPermute, 
					lengthOfIndicesToPermute);// compute scores using permuted class labels
			
			int length = featureIndicesToPermute!=null?lengthOfIndicesToPermute:numFeatures;
			lengthOfIndicesToPermute = 0;
			if(length==0) {
				bank = 0;
			}
			for (int i = 0; i < length && bank > 0; i++) {
				bank--;
				int index = featureIndicesToPermute!=null?featureIndicesToPermute[i]:i;
				double score = scores[index];
				if (testDirection == Constants.TWO_SIDED
						|| testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
					if (permutedScores[index] >= score) {
						featureSpecificPValues[index] += 1.0;
					}
				} else {
					if (permutedScores[index] <= score) {
						featureSpecificPValues[index] += 1.0;
					}
				}
				
				if(removeFeatures) {
					permutationsPerFeature[index]++;
					double k = featureSpecificPValues[index];
					int N = permutationsPerFeature[index];
					double d = 2*1.96*Math.sqrt((N+1-k)/(N*(k+1)));
					double d2 = -Double.MAX_VALUE;
					
					if(testDirection == Constants.TWO_SIDED) {
						k = N-k;
						d2 = 2*1.96*Math.sqrt((N+1-k)/(N*(k+1)));
					}
					
					if (d >= theta || d2 >= theta) { // include marker
						featureIndicesToPermute[lengthOfIndicesToPermute] = index;
						lengthOfIndicesToPermute++;
					} 
					
				}
			}
			
			
			/*Sorting.sort(permutedScores, Sorting.DESCENDING);

			 for(int i = 0; i < numFeatures; i++) {
			 double score = scores[descendingIndices[i]];

			 if(testDirection == TWO_SIDED) {
			 if(score >= 0) {
			 if(permutedScores[i] >= score) {
			 rankBasedPValues[i] += 1.0;
			 }
			 } else {
			 if(permutedScores[i] <= score) {
			 rankBasedPValues[i] += 1.0;
			 }
			 }
			 } else if(testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
			 if(permutedScores[i] >= score) {
			 rankBasedPValues[i] += 1.0;
			 }
			 } else {
			 if(permutedScores[i] <= score) {
			 rankBasedPValues[i] += 1.0;
			 }
			 }
			 }*/

			if (!removeFeatures) {
				for (int i = 0; i < numFeatures; i++) {
					permutedScores[i] = Math.abs(permutedScores[i]);
					monotonicPermutedScores[i] = permutedScores[i];
				}

				for (int i = numFeatures - 2; i >= 0; i--) { // make permuted
																// scores
																// monotonicically
																// non-decreasing
																// for maxT
					int index = descendingAbsIndices[i];
					int indexPlusOne = descendingAbsIndices[i + 1];

					if (monotonicPermutedScores[indexPlusOne] > monotonicPermutedScores[index]) {
						monotonicPermutedScores[index] = monotonicPermutedScores[indexPlusOne];
					}

				}

				for (int i = 0; i < numFeatures; i++) {
					if (monotonicPermutedScores[i] >= absoluteScores[i]) {
						maxT[i] += 1.0;
					}
				}

				Sorting.sort(permutedScores, Sorting.DESCENDING);
				int j = 0;
				int count = 0;
				for (int i = 0; i < numFeatures; i++) {
					double score = absoluteScores[descendingAbsIndices[i]];

					while (j < numFeatures && score < permutedScores[j]) {
						count++;
						j++;
					}
					// fpr[descendingAbsIndices[i]] += count;

					if (count > 0) {
						fwer[descendingAbsIndices[i]]++;
					}
				}
			}

		}
	}

	/**
	 *  Reads the next permutation from a file
	 *
	 * @return    the next permutation
	 */
	final int[] nextPermutation() {
		try {
			String s = testClassPermutationsReader.readLine();
			StringTokenizer st = new StringTokenizer(s, " \t");
			int cols = dataset.getColumnCount();
			int[] a = new int[cols];

			for (int j = 0; j < cols; j++) {
				a[j] = Integer.parseInt(st.nextToken());
			}

			return a;
		} catch (Exception e) {
			AnalysisUtil
					.exit(
							"An error occurred while reading the permutations file.",
							e);
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
		} catch (Exception e) {
			AnalysisUtil
					.exit(
							"An error occurred while reading the permutations file.",
							e);
		}
	}

	static class Debugger {
		java.util.Map assignment2Occurences = new java.util.HashMap();

		public static String toString(int[] a) {
			StringBuffer buf = new StringBuffer();
			for (int i = 0; i < a.length; i++) {
				buf.append(a[i]);
			}
			return buf.toString();
		}

		public void addAssignment(int[] a) {
			String s = toString(a);
			Integer occurences = (Integer) assignment2Occurences.get(s);
			if (occurences == null) {
				occurences = new Integer(0);
			}
			assignment2Occurences
					.put(s, new Integer(occurences.intValue() + 1));
		}

		public void print() {
			System.out.println("Unique permutations "
					+ assignment2Occurences.size());
			System.out.println("permutations");
			System.out.println(assignment2Occurences);
		}
	}

	static {
		if ("true".equalsIgnoreCase(System
				.getProperty("edu.mit.broad.marker.debug"))) {
			debugger = new Debugger();
		}
	}

}
