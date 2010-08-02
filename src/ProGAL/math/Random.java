package ProGAL.math;

import java.util.ArrayList; 
import java.util.List;

public class Random {
	private static java.util.Random rand = new java.util.Random(); 

	/**
	 * Return a uniform random number between (including) i1 and (not including) i2. 
	 */
	public static int randBetween(int i1, int i2){
		int max = Math.max(i1, i2);
		int min = Math.min(i1, i2);
		return (int)(rand.nextDouble()*(max-min)+min);
	}
	/**
	 * Return a uniform random double between (including) d1 and (not including) d2
	 */
	public static double randBetween(double d1, double d2){
		double max = Math.max(d1, d2);
		double min = Math.min(d1, d2);
		return rand.nextDouble()*(max-min)+min;
	}


	/**
	 * Generate a random permutation of integers between (including) 0 and (not including) max
	 */
	public static int[] randomPermutation(int max){
		int[] ret = new int[max];
		for(int i=0;i<max;i++){ ret[i] = i;}
		randomizeInPlace(ret);
		return ret;
	}
	
	/**
	 * Randomize the array as.
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
	 * Generates all permutations of integers from 0 (inclusive) to max (exclusive).
	 * TODO: Move this method
	 * @param max the length of each permutation
	 * @return a list of all permutations of length max
	 */
	public static List<int[]> getAllPermutations(int max){
		int[] begin = {};

		int[] end = new int[max];
		for(int i=0;i<max;i++) end[i] = i;
		
		return permute(begin, end);
	}

	/** 
	 * @see <a href="http://www.java2s.com/Tutorial/Java/0100__Class-Definition/RecursivemethodtofindallpermutationsofaString.htm">link</a>
	 */
	private static List<int[]> permute(int[] beginning, int[] ending) {
		List<int[]> ret = new ArrayList<int[]>();
		if (ending.length <= 1){
			int[] perm = new int[beginning.length+ending.length];
			System.arraycopy(beginning, 0, perm, 0, beginning.length);
			System.arraycopy(ending, 0, perm, beginning.length, perm.length-beginning.length);
			ret.add(perm);
		}else{
			for (int i = 0; i < ending.length; i++) {
				try {
					int[] newEnd = new int[ending.length-1];
					System.arraycopy(ending, 0, newEnd, 0, i);
					System.arraycopy(ending, i+1, newEnd, i, ending.length-1-i);

					int[] newBegin = new int[beginning.length+1];
					System.arraycopy(beginning, 0, newBegin, 0, beginning.length);
					newBegin[beginning.length] = ending[i];
					ret.addAll(permute(newBegin, newEnd));
					
				} catch (StringIndexOutOfBoundsException exception) {
					exception.printStackTrace();
				}
			}
		}
		return ret;
	}

	/**
	 * Calculate the binomial coefficient.
	 * TODO: Move this method
	 * @param n a number larger than or equal 0
	 * @param k a number larger than or equal 0 but smaller than or equal n
	 * @return n choose k
	 * @see <a href="http://www.brpreiss.com/books/opus5/html/page460.html">link</a>
	 */
	public static int binom(int n, int m){
		int[] b = new int[n+1];
		b[0] = 1;
		for(int i=1;i<=n;++i){
			b[i] = 1;
			for(int j=i-1;j>0;--j)
				b[j]+=b[j-1];
		}
		return b[m];
	}

	/**
	 * Seeds the random generater used by the randPermutation and randBetween methods.
	 * @param s
	 */
	public static void seed(long s){
		rand = new java.util.Random(s);
	}
}
