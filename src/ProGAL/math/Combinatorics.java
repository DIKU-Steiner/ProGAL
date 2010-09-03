package ProGAL.math;

import java.util.ArrayList;
import java.util.List;

/** 
 * A wrapper for static combinatorial helper functions.
 * @author rfonseca
 */
public class Combinatorics {

	/**
	 * Generates all permutations of integers from 0 (inclusive) to max (exclusive).
	 * @param max the length of each permutation
	 * @return a list of all permutations of length max
	 */
	public static List<int[]> getAllPermutations(int max){
		if(max<0) throw new RuntimeException("max must be nonnegative");
		
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
	 * @param n a number larger than or equal 0
	 * @param k a number larger than or equal 0 but smaller than or equal n
	 * @return n choose k
	 * @see <a href="http://www.brpreiss.com/books/opus5/html/page460.html">link</a>
	 */
	public static int binom(int n, int m){
		if(n<0) return binom(-n,m);
		if(m<0 || m>n) return 0;
		
		int[] b = new int[n+1];
		b[0] = 1;
		for(int i=1;i<=n;++i){
			b[i] = 1;
			for(int j=i-1;j>0;--j)
				b[j]+=b[j-1];
		}
		return b[m];
	}

}
