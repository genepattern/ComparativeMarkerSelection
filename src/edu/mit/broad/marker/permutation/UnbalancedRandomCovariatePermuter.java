package edu.mit.broad.marker.permutation;

import edu.mit.broad.dataobj.ClassVector;
import java.util.*;

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
      this.choose = new int[covariate.levels()];
      this.size = classVector.size();
      int max = 0;
      for(int i = 0, levels = covariate.levels(); i < levels; i++) {
         int[] covariateIndices = covariate.getIndices(covariate.getLevel(i));
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
      for(int i = 0, levels = covariate.levels(); i < levels; i++) {
         int[] covariateIndices = covariate.getIndices(covariate.getLevel(i));
         int n = choose[i];
         cern.jet.random.sampling.RandomSampler.sample(n, covariateIndices.length, n, 0, values, 0, random);
         for(int j = 0; j < n; j++) {
            int index = covariateIndices[(int)values[j]];
            assignments[index] = 1;
         }
      }

      return assignments;
   }


   public static void main(String[] args) {
      int[] levels = {0,1};
      ClassVector assignments = new ClassVector(new int[]{1,1,0,0,0,0}, levels);
      ClassVector covariate = new ClassVector(new int[]{0,1,1,1,1,0}, levels);
      UnbalancedRandomCovariatePermuter perm = new UnbalancedRandomCovariatePermuter(assignments, covariate);
      HashMap map = new HashMap();
      for(int i = 0; i < 100; i++) {
         String p = toString(perm.next());
         Integer count = (Integer) map.get(p);
         if(count == null) {
            count = new Integer(1);
         } else {
            count = new Integer(count.intValue() + 1);
         }
         map.put(p, count);
      }
      System.out.println(map);
   }

   static String toString(int[] a) {
      String s = "";
      for(int i = 0; i < a.length; i++) {
         s += a[i];
      }
      return s;
   }
}
