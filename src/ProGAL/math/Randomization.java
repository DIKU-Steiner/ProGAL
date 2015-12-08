package ProGAL.math;

import java.util.Random;

/**
 * A wrapper for static randomization functions.
 * @author rfonseca
 */
public class Randomization {
	private static java.util.Random rand = new java.util.Random(); 

	
	public static Random getGenerator() {
		return rand;
	}
	
	/**
	 * Return a uniform random number between (including) i1 and (not including) i2.
	 * @param i1 lower bound of random number
	 * @param i2 upper bound of random number 
	 */
	public static int randBetween(int i1, int i2){
		int max = Math.max(i1, i2);
		int min = Math.min(i1, i2);
		return (int)(rand.nextDouble()*(max-min)+min);
	}
	/**
	 * Return a uniform random double between (including) d1 and (not including) d2
	 * @param d1 lower bound of random number
	 * @param d2 upper bound of random number
	 */
	public static double randBetween(double d1, double d2){
		double max = Math.max(d1, d2);
		double min = Math.min(d1, d2);
		return rand.nextDouble()*(max-min)+min;
	}


	/**
	 * Generate a random permutation of integers between (including) 0 and (not including) max
	 * @param max the length of the permutation
	 */
	public static int[] randomPermutation(int max){
		int[] ret = new int[max];
		for(int i=0;i<max;i++){ ret[i] = i;}
		randomizeInPlace(ret);
		return ret;
	}
	
	/**
	 * Randomize an array.
	 * @param as the array to be randomized
	 */
	public static void randomizeInPlace(int[] as){
		for(int i=0;i<as.length;i++){
			int r = randBetween(i,as.length);
			int t = as[i];
			as[i] = as[r];
			as[r] = t;
		}
	}


	/**
	 * Seeds the random generater used by this class.
	 * @param s seed
	 */
	public static void seed(long s){
		rand = new java.util.Random(s);
	}
	
	public static Object[] shuffle(Object[] arr){
		for (int i=arr.length; i>1; i--){
			int r = randBetween(0,i);
			Object tmp = arr[i-1];
			arr[i-1] = arr[r];
			arr[r] = tmp;
		}
		return arr;
	}
	
}
