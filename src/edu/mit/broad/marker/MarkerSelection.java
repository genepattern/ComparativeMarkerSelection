package edu.mit.broad.marker;
import java.io.*;
import java.io.FileReader;

import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;

import edu.mit.broad.dataobj.*;
import edu.mit.broad.dataobj.microarray.*;
import edu.mit.broad.gp.*;
import edu.mit.broad.io.microarray.*;

import edu.mit.broad.marker.permutation.*;

/**
 *@author     Joshua Gould
 *@created    September 17, 2004
 */
public class MarkerSelection {

	/**  input expression values */
	Dataset dataset;

	/**  input classes */
	ClassVector classVector;
	/**  number of features to select */
	int numFeatures;
	/**  number of permutations to perform */
	int numPermutations;
	final static int CLASS_ZERO_GREATER_THAN_CLASS_ONE = 0;
	final static int CLASS_ZERO_LESS_THAN_CLASS_ONE = 1;
	final static int TWO_SIDED = 2;

	/**  whether to do one-sided >, one-sided <, or two-sided test */
	int side;

	int metric;
	final static int T_TEST = 0;
	final static int SNR = 1;

	/**  whether permutations are balanced */
	boolean balanced;

	/**  whether to compute all possible permutations */
	boolean complete;
	/**  whether to read permutations from a file */
	boolean readPermutations = false;
	/** reader when reading permutations from file */
	BufferedReader testClassPermutationsReader;
	
	/** name of output file */
	String outputFileName;

	double[][] rankBasedScores;
	double[][] geneSpecificScores;

	/** unpermuted scores */
	double[] unpermutedScores;

	/** sorted indices in the input dataset of the top numFeatures */
	int[] topFeaturesIndices;

	/** number of rows in input data */
	int N;
	public MarkerSelection(Dataset _dataset, ClassVector _classVector,
			int _numFeatures, int _numPermutations, int _side,
			String _outputFileName, boolean _balanced,
			boolean complete, String permutationsFile) {
		this.dataset = _dataset;
		this.classVector = _classVector;
		
		GPUtil.checkDimensions(dataset, classVector);
		
		this.numFeatures = _numFeatures;
		this.numPermutations = _numPermutations;
		this.side = _side;
		this.balanced = _balanced;
		this.outputFileName = _outputFileName;
		this.complete = complete;
		this.N = dataset.getRowDimension();
		
		if(permutationsFile != null) {
			readPermutations = true;
			initFile(permutationsFile);
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


	


	boolean containsAtLeastOneGreater(double score, double[] values) {
		for(int i = 0; i < values.length; i++) {
			if(Math.abs(values[i]) >= score) {
				return true;
			}
		}
		return false;
	}


	int countColumnsGreaterThan(double score, double[][] values) {
		int count = 0;
		for(int j = numPermutations - 1; j >= 0; j--) {
			boolean found = false;
			for(int i = 0; i < values.length; i++) {
				if(Math.abs(values[i][j]) >= score) {
					count++;
					break;
				}
			}
		}
		return count;
	}


	public static void main(String[] args)
			 throws Exception {

		String datasetFile = args[0];
		String clsFile = args[1];
		int _numFeatures = Integer.parseInt(args[2]);
		int _numPermutations = Integer.parseInt(args[3]);
		int side = Integer.parseInt(args[4]);
		String outputFileName = args[5];
		boolean balanced = Boolean.valueOf(args[6]).booleanValue();
		boolean complete = Boolean.valueOf(args[7]).booleanValue();
		String permutationsFile = null;
		if(args.length==9) {
			permutationsFile = args[8];
		}
		ExpressionDataReader reader = GPUtil.getExpressionReader(datasetFile);
		ExpressionDataImpl e = (ExpressionDataImpl) reader.read(datasetFile);
		new MarkerSelection(e.getExpressionMatrix(),
				GPUtil.getClassVector(clsFile), _numFeatures,
				_numPermutations, side, outputFileName, balanced,
				complete, permutationsFile);
	}
	
	static int countNumberGreater(double score, double[] permutedScores) {
		int count = 0;
		for(int i = 0, length = permutedScores.length; i < length; i++) {
			if(Math.abs(permutedScores[i]) >= Math.abs(score)) {
				count++;
			}
		}
		return count;
	}

	void computePValues() {
		int[] classZeroIndices = classVector.getIndices(0);
		int[] classOneIndices = classVector.getIndices(1);

		if(balanced) {

			if(classZeroIndices.length != classOneIndices.length) {
				GPUtil.exit(
						"The number of items in each class must be equal for balanced permutations.");
			}

			if((classZeroIndices.length % 2) != 0) {
				GPUtil.exit(
						"The number of items in each class 0 must be an even number.");
			}

			if((classOneIndices.length % 2) != 0) {
				GPUtil.exit(
						"The number of items in each class 1 must be an even number.");
			}
		}
						
		
		GeneSpecificFeatureSelector geneSpecificFeatureSelector = new GeneSpecificFeatureSelector();
		
		unpermutedScores = new double[N];
		
		geneSpecificFeatureSelector.compute(dataset,
				classZeroIndices,
				classOneIndices, unpermutedScores); 

		int sortOrder = -1;
		if(side==CLASS_ZERO_LESS_THAN_CLASS_ONE) {
			sortOrder = Util.ASCENDING;
		} else if(side==CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
			sortOrder = Util.DESCENDING;
		} else {
			sortOrder = Util.ABSOLUTE;	
		}
		topFeaturesIndices = Util.getIndices(unpermutedScores, sortOrder); 
	
		Permuter permuter = null;

		if(!complete && balanced) {
			permuter = new BalancedRandomPermuter(classZeroIndices, classOneIndices);
		} else if(!complete && !balanced) {
			permuter = new UnbalancedRandomPermuter(classVector.size(), classVector.getIndices(1).length);
		} else if(complete && !balanced) {
			permuter = new UnbalancedCompletePermuter(classVector.size(),
					classZeroIndices.length);
			java.math.BigInteger totalPermutations = ((UnbalancedCompletePermuter) (permuter)).getTotal();
			if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
				GPUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
			}
			numPermutations = totalPermutations.intValue();
		} else if(complete && balanced) {
			permuter = new BalancedCompletePermuter(classZeroIndices, classOneIndices);
			java.math.BigInteger totalPermutations = ((BalancedCompletePermuter) (permuter)).getTotal();
			if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
				GPUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
			}
		}
		
		rankBasedScores = new double[numFeatures][numPermutations]; // FIXME
		
		geneSpecificScores = new double[N][numPermutations]; // FIXME
		double[] fwer = new double[numFeatures];
		double[] fpr = new double[N];
		
		double[] permutedScores = new double[N];
		int[] levels = {0, 1};
		
		for(int perm = 0; perm < numPermutations; perm++) {
			int[] permutedAssignments = null;
			int[] permutedClassZeroIndices = null;
			int[] permutedClassOneIndices = null;

			if(readPermutations) {
				permutedAssignments = nextPermutation();
				ClassVector temp = new ClassVector(permutedAssignments,
						levels);
				permutedClassZeroIndices = temp.getIndices(0);
				permutedClassOneIndices = temp.getIndices(1);

			} else {
				permutedAssignments = permuter.next();
				ClassVector temp = new ClassVector(permutedAssignments,
						levels);

				permutedClassZeroIndices = temp.getIndices(0);
				permutedClassOneIndices = temp.getIndices(1);

			}

			geneSpecificFeatureSelector.compute(dataset,
					permutedClassZeroIndices,
					permutedClassOneIndices, permutedScores); // compute scores using permuted class labels
			
			for(int i = 0; i < N; i++) {
				geneSpecificScores[i][perm] = permutedScores[i];
			}

			if(side == TWO_SIDED) {
				for(int i = 0; i < permutedScores.length; i++) {
					permutedScores[i] = Math.abs(permutedScores[i]);
				}

			}
			
			Arrays.sort(permutedScores); // ascending sort
			if(side != CLASS_ZERO_LESS_THAN_CLASS_ONE) {
				for(int feature = 0; feature < numFeatures; feature++) {
					rankBasedScores[feature][perm] = permutedScores[permutedScores.length - 1 - feature];
				}
			} else {
				for(int feature = 0; feature < numFeatures; feature++) {
					rankBasedScores[feature][perm] = permutedScores[feature];
				}
			}
			for(int i = 0; i < numFeatures; i++) {
				if(containsAtLeastOneGreater(unpermutedScores[topFeaturesIndices[i]], permutedScores)) {
					fwer[i] = fwer[i] + 1.0;
				}
			}
			for(int i =0; i < N; i++) {
				fpr[i] += countNumberGreater(unpermutedScores[i], permutedScores);
			}
		}
		
		for(int i = 0; i < numFeatures; i++) {
			fwer[i] /= numPermutations;
		}
		
		for(int i = 0; i < N; i++) {
			fpr[i] /= N;
			fpr[i] /= numPermutations;
		}
		
		double[] rankBasedPValues = computeRankBasedPValues();
		double[] geneSpecificPValues = computeGeneSpecificPValues();

		double[] fdr = new double[numFeatures];
		int[] pValueIndices = Util.getIndices(fpr, Util.ASCENDING);
		int[] ranks = Util.rank(pValueIndices);
		
		
		// check for ties
		for(int i = pValueIndices.length-1; i > 0; i--) {
			double bigPValue = fpr[pValueIndices[i]];
			double smallPValue = fpr[pValueIndices[i-1]];
			if(bigPValue==smallPValue) {
				ranks[pValueIndices[i-1]] = ranks[pValueIndices[i]];
			}
		}
		
		
		for(int i = 0; i < numFeatures; i++) {
			int index = topFeaturesIndices[i];
			int rank = ranks[index];
			double p = fpr[index];
			fdr[i] = (p * N)/rank;
		}
		
		PrintWriter pw = null;
		java.util.Vector features = new java.util.Vector(numFeatures);

		try {
			pw = new PrintWriter(new FileWriter(outputFileName));

			for(int feature = 0; feature < numFeatures; feature++) {
				int originalIndex = topFeaturesIndices[feature];
				features.add(dataset.getRowName(originalIndex));
				
				pw.println(features.get(feature) + "\t" +
						rankBasedPValues[feature] + "\t" +
						geneSpecificPValues[originalIndex] + "\t " +
						fwer[feature] + "\t" + fpr[originalIndex] + "\t" + fdr[feature]);
			}
		} catch(Exception e) {
			e.printStackTrace();
			GPUtil.exit("An error occurred while saving the output file.", e);
			
		}

		if(pw != null) {

			try {
				pw.close();
			} catch(Exception x) {
			}
		}
		
		double[] topFeaturesGeneSpecificPValues = new double[numFeatures];
		double[] topFeatures_fpr = new double[numFeatures];
		
			
		for(int i = 0; i < numFeatures; i++) {
			int originalIndex = topFeaturesIndices[i];
			topFeatures_fpr[i] = fpr[originalIndex];
			topFeaturesGeneSpecificPValues[i] = geneSpecificPValues[originalIndex];
			
		}
		new MarkerSelectionFrame(features, topFeaturesGeneSpecificPValues,
				rankBasedPValues, fwer, fdr, topFeatures_fpr);
	}

	double[] computeRankBasedPValues() {
		double[] rankBasedPValues = new double[numFeatures];
			for(int i = 0; i < numFeatures; i++) {
				double score = unpermutedScores[topFeaturesIndices[i]];
				for(int j = 0; j < numPermutations; j++) {
					if(side==CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
						if(rankBasedScores[i][j] >= score) {
							rankBasedPValues[i] += 1.0;
						}
					} else if(side==CLASS_ZERO_LESS_THAN_CLASS_ONE) {
						if(rankBasedScores[i][j] <= score) {
							rankBasedPValues[i] += 1.0;
						}
					} else {
						if(Math.abs(rankBasedScores[i][j]) >= Math.abs(score)) {
							rankBasedPValues[i] += 1.0;
						}
					}
				}
				rankBasedPValues[i] /=numPermutations;
			}	
			return rankBasedPValues;
	}
	
	double[] computeGeneSpecificPValues() {
		double[] geneSpecificPValues = new double[N];
			for(int i = 0; i < N; i++) {
				double score = unpermutedScores[i];
				for(int j = 0; j < numPermutations; j++) {
					if(side==CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
						if(geneSpecificScores[i][j] >= score) {
							geneSpecificPValues[i] += 1.0;
						}
					} else if(side==CLASS_ZERO_LESS_THAN_CLASS_ONE) {
						if(geneSpecificScores[i][j] <= score) {
							geneSpecificPValues[i] += 1.0;
						}
					} else {
						if(Math.abs(geneSpecificScores[i][j]) >= Math.abs(score)) {
							geneSpecificPValues[i] += 1.0;
						}
					}
					
				}
				geneSpecificPValues[i] /=numPermutations;
			}	
		return geneSpecificPValues;
	}

	/** 
	* Reads the next permutation from a file 
	*/
	int[] nextPermutation() {
		try {
			String s = testClassPermutationsReader.readLine();
			StringTokenizer st = new StringTokenizer(s, " \t");
			int cols = dataset.getColumnDimension();
			int[] a = new int[cols];
	
			for(int j = 0; j < cols; j++) {
				a[j] = Integer.parseInt(st.nextToken());
			}
	
			return a;	
		} catch(Exception e) {
			GPUtil.exit("An error occurred while reading the permutations file.", e);
			return null;
		}
		
	}
	
	/** Opens the permutations file for reading */
	void initFile(String fileName) {

		try {
			testClassPermutationsReader = new BufferedReader(new FileReader(
					fileName));
		} catch(Exception e) {
			GPUtil.exit(
					"An error occurred while reading the permutations file.",
					e);
		}
	}


	/**
	 *@author     Joshua Gould
	 */
	static class GeneSpecificFeatureSelector {

		public void compute(Dataset dataset,
				int[] classOneIndices,
				int[] class2Indices,
				double[] scores) {

			int rows = dataset.getRowDimension();

			for(int i = 0; i < rows; i++) {

				double class1Mean = Util.mean(dataset, classOneIndices, i);
				double class2Mean = Util.mean(dataset, class2Indices, i);
				double class1Std = Util.standardDeviation(dataset, classOneIndices,
						i, class1Mean);
				double class2Std = Util.standardDeviation(dataset, class2Indices, i,
						class2Mean);
				double Sxi = (class1Mean - class2Mean) / (class1Std +
						class2Std);
				scores[i] = Sxi;
			}
		}
	}
}

