package edu.mit.broad.marker.permutation;
import java.util.Arrays;
import java.util.HashMap;

import edu.mit.broad.data.matrix.ClassVector;

/**
 *  Permuter using covariates
 *
 * @author    Joshua Gould
 */
public class UnbalancedRandomCovariatePermuter implements Permuter {
   cern.jet.random.engine.MersenneTwister random;
   long[] values;
   int[] classOneIndices;
   ClassVector covariate;
   int[] choose;
   int size;
   
   public UnbalancedRandomCovariatePermuter(ClassVector classVector, ClassVector covariate) {
      this.random = new cern.jet.random.engine.MersenneTwister();
      this.covariate = covariate;
      this.classOneIndices = classVector.getIndices(1);
      this.choose = new int[covariate.getClassCount()];
      this.size = classVector.size();
      int max = 0;
      for(int i = 0, levels = covariate.getClassCount(); i < levels; i++) {
         int[] covariateIndices = covariate.getIndices(i);
         // find number of class one that belong to the ith covariate
         int numClassOne = 0;
         for(int j = 0; j < classOneIndices.length; j++) {
            if(Arrays.binarySearch(covariateIndices, classOneIndices[j]) >= 0) {
               numClassOne++;
            }
         }
         choose[i] = numClassOne;
         max = Math.max(max, numClassOne);
      }
      values = new long[max];
   }


   public int[] next() {
      int[] assignments = new int[size];
      for(int i = 0, levels = covariate.getClassCount(); i < levels; i++) {
         int[] covariateIndices = covariate.getIndices(i);
         int n = choose[i];
         cern.jet.random.sampling.RandomSampler.sample(n, covariateIndices.length, n, 0, values, 0, random);
         for(int j = 0; j < n; j++) {
            int index = covariateIndices[(int) values[j]];
            assignments[index] = 1;
         }
      }
      return assignments;
   }
}
