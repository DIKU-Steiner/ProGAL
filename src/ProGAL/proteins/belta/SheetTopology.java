package ProGAL.proteins.belta;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/** 
 * Representation of the topology of a single sheet. Where the BetaTopology is represented 
 * using a pairing matrix, the sheet topology uses a list of strand pairs that can be accessed 
 * via the public <code>strandPairs</code> field. The list of strand pairs is ordered such that 
 * the pairs follow the order along the sheet, and the first strand from the N-terminal is in 
 * the first half of the list.
 * 
 * <code>getStrandOrder</code> and </code>getStrandOrientation</code> are alternative methods 
 * to access the topology. The following example shows the relationship between the pairing-
 * matrix in the beta-topology and the order and orientation-arrays.
 * <pre>
 *   SecondaryStructure ss = new SecondaryStructure(" EE EE EE EE ");
 *   //Beta-topology with a single sheet
 *   BetaTopology bTop = new BetaTopology(ss);
 *   bTop.setPaired(1,0);
 *   bTop.setPaired(2,0);
 *   bTop.setPaired(2,3);
 *   SheetTopology sTop = new SheetTopology(bTop,0);
 *   System.out.println(sTop.getStrandOrderString());       //Prints "1 0 2 3"
 *   System.out.println(sTop.getStrandOrientationString()); //Prints "1110"
 * </pre>
 * @author R.Fonseca
 */
public class SheetTopology{
	
	/** The secondary structure that the sheet-topology is related to */
	public final SecondaryStructure secondaryStructure;
	
	/** An ordered list of strand-pairs. The order follows the order of pairs when 
	 * "climbing" the sheet-ladder from one end to the other. The first strand from 
	 * the N-terminal is in the first half of this list. If it is in the middle, then 
	 * the second strand is in the first half. This corresponds to the convention in 
	 * Ruczinskis paper.*/
	public final List<StrandPair> strandPairs = new ArrayList<StrandPair>();
	
	/** A sorted list of strand-indices. A strand-index is used in the 
	 * <code>secondaryStructure.getStrands()</code>-array.*/
	public final List<Integer> strands = new ArrayList<Integer>();

	/** 
	 * Construct the sheet containing the specified strand. The order of the 
	 * strandPairs list will be such that the lowest index strand will be in the first 
	 * half of the list. If it is in the center the second-lowest index strand will 
	 * be in the first part of the list. 
	 */
	public SheetTopology(BetaTopology bTop, int strand){
		this.secondaryStructure = bTop.secondaryStructure;
		
		ArrayList<StrandPair> tmp = new ArrayList<StrandPair>();
		Set<Integer> fringe = new TreeSet<Integer>();
		fringe.add(strand);

		while(!fringe.isEmpty()){
			Integer s0 = fringe.iterator().next();
			fringe.remove(s0);

			for(int s=0;s<bTop.N;s++){
				if(bTop.pair(s0,s) && !pairExists(tmp, s0,s)) {
					tmp.add(new StrandPair(s0,s,s<s0));
					fringe.add(s);
				}
				if(bTop.pair(s,s0) && !pairExists(tmp, s0,s)){
					tmp.add(new StrandPair(s0,s,s>s0));
					fringe.add(s);
				}
			}
		}
		if(tmp.isEmpty()){
			System.out.println("Sheet(...) tmp empty .. ");
		}
		//Store all strands
		for(StrandPair bp: tmp){
			if(!strands.contains(bp.strand1))	strands.add(bp.strand1);
			if(!strands.contains(bp.strand2))	strands.add(bp.strand2);
			Collections.sort(strands);
		}

		//Sort the sheets
		//  Find an edge strand
		int[] pairs = new int[bTop.N]; for(int s=0;s<bTop.N;s++) pairs[s] = 0;
		for(StrandPair bp: tmp){
			pairs[bp.strand1]++;
			pairs[bp.strand2]++;
		}
		int edgeStrand = -1;
		for(int s=0;s<bTop.N;s++) if(pairs[s]==1){ edgeStrand = s; break; }
		
		//  Sort
		int nextStrand = edgeStrand;
		if(nextStrand<0) {//Its a barrel
			nextStrand = tmp.get(0).strand1;
			while(!tmp.isEmpty()){
				StrandPair nextPair = null;
				for(StrandPair bp: tmp) if(bp.contains(nextStrand)) { nextPair = bp; break;}
				strandPairs.add(nextPair); 
				tmp.remove(nextPair); 

				if(nextPair.strand2==nextStrand){ 
					int t = nextPair.strand1;
					nextPair.strand1 = nextPair.strand2;
					nextPair.strand2 = t;
				}		

				nextStrand = strandPairs.get(strandPairs.size()-1).strand2;
			}
		}else{
			while(!tmp.isEmpty()){
				StrandPair nextPair = null;
				for(StrandPair bp: tmp) if(bp.contains(nextStrand)) { nextPair = bp; break;}
				strandPairs.add(nextPair); 
				tmp.remove(nextPair); 

				if(nextPair.strand2==nextStrand){ 
					int t = nextPair.strand1;
					nextPair.strand1 = nextPair.strand2;
					nextPair.strand2 = t;
				}

				nextStrand = strandPairs.get(strandPairs.size()-1).strand2;
			}
		}

		//Find the min-index strand and ensure that it occurs early in the list
		int minStrand = Integer.MAX_VALUE;
		int minIdx = -1, c = 0;
		for(StrandPair bp: strandPairs){
			if(minStrand>bp.strand1) {
				minStrand = bp.strand1;
				minIdx = c;
			}
			c++;
		}
		if(strands.size()%2==1 && minIdx==(strands.size()/2)){
			//Repeat for the second lowest strand;
			c=0;
			int minStrand2 = Integer.MAX_VALUE;
			for(StrandPair bp: strandPairs){
				if(bp.strand1<minStrand2 && bp.strand1!=minStrand){
					minStrand2 = bp.strand1;
					minIdx = c;
				}
				c++;
			}
		}

		if(minIdx>=strands.size()/2){ //Reverse list
			Collections.reverse(strandPairs);
			for(StrandPair bp: strandPairs) {
				int t=bp.strand1;
				bp.strand1=bp.strand2;
				bp.strand2 = t;
			}
		}
	}


	/**
	 * Return true if the specified strand is in this sheet.
	 * @param strand A strand-index (i.e. an index in the 
	 * <code>secondaryStructure.getStrands()</code>-array). 
	 */
	public boolean containsStrand(int strand){
		return strands.contains(strand);
	}
 
	/**
	 * Return an array indicating the strand order. This method can be more 
	 * convenient than traversing the strand-pairs in order. 
	 */
	public int[] getStrandOrder(){
		int[] ret = new int[strandPairs.size()+1];
		for(int i=0;i<strandPairs.size();i++)
			ret[i] = strandPairs.get(i).strand1;
		ret[ret.length-1] = strandPairs.get(strandPairs.size()-1).strand2;
		return ret;
	}

	/**
	 * Return an array indicating the normalized strand order. Members of the normalized 
	 * order refers to indexes in the <code>SheetTopology.this.strands</code>-array and 
	 * not to <code>SheetTopology.this.secondaryStructure.getStrands()</code> as strand-
	 * indices usually do.
	 */
	public int[] getNormalizedStrandOrder(){
		int[] realOrder = getStrandOrder();
		int[] sortedOrder = getStrandOrder();
		Arrays.sort(sortedOrder);
		int[] normOrder = new int[sortedOrder.length];
		for(int i=0;i<sortedOrder.length;i++){
			normOrder[i] = indexOf(realOrder[i], sortedOrder);
		}
		return normOrder;
	}
	
	/**
	 * Return an array indicating orientation of strands. A 1-entry indicates the strand 
	 * is pointing up and a 0-entry down.  
	 */
	public int[] getStrandOrientation(){
		boolean lastUp = true;
		boolean[] ups = new boolean[strands.size()];
		int c=0;
		ups[c++] = lastUp;
		for(StrandPair bp: strandPairs){
			if(!lastUp && !bp.parallel) lastUp = true;
			else if(lastUp && !bp.parallel) lastUp = false;
			ups[c++] = lastUp;
		}
		int minStrand = getMinStrand();
		c=0;
		for(StrandPair bp: strandPairs){
			if(bp.strand1==minStrand){
				if(!ups[c]) for(c=0;c<ups.length;c++) ups[c] = !ups[c];
				break;
			}
			c++;
		}
		int[] ret = new int[strands.size()];
		c=0;
		for(Boolean up: ups) ret[c++] = up?1:0;
		return ret;
	}

	public String toString(){
		return "Sheet< "+getStrandOrderString()+", "+getStrandOrientationString()+" >";
	}
	
	/** String representation of strand order. */
	public String getStrandOrderString(){
		StringBuilder sb = new StringBuilder();
		for(StrandPair bp: strandPairs){
			sb.append(bp.strand1+" ");
		}
		sb.append(strandPairs.get(strandPairs.size()-1).strand2);
		return sb.toString();
	}
	
	/** String representation of strand orientation. */
	public String getStrandOrientationString(){
		StringBuilder ret = new StringBuilder();
		boolean lastUp = true;
		boolean[] ups = new boolean[strands.size()];
		int c=0;
		ups[c++] = lastUp;
		for(StrandPair bp: strandPairs){
			if(!lastUp && !bp.parallel) lastUp = true;
			else if(lastUp && !bp.parallel) lastUp = false;
			if(c==strands.size()) break;
			ups[c++] = lastUp;
		}
		int minStrand = getMinStrand();
		c=0;
		for(StrandPair bp: strandPairs){
			if(bp.strand1==minStrand){
				if(!ups[c]) for(c=0;c<ups.length;c++) ups[c] = !ups[c];
				break;
			}
			c++;
		}
		for(Boolean up: ups) ret.append(up?"1":"0");
		return ret.toString();
	}
	
	private int getMinStrand(){
		int minStrand = strandPairs.get(0).strand1;
		for(StrandPair bp: strandPairs){
			if(bp.strand1<minStrand) minStrand = bp.strand1;
		}
		int t= strandPairs.get(strandPairs.size()-1).strand2;
		if(t<minStrand) minStrand = t;
		return minStrand;
	}
	
	private static int indexOf(int c, int[] array){
		for(int i=0;i<array.length;i++) if(array[i]==c) return i;
		return -1;
	}

	private static boolean pairExists(List<StrandPair> bPairs, int s1, int s2){
		for(StrandPair bp: bPairs){
			if( (bp.strand1==s1 && bp.strand2==s2)||(bp.strand1==s2&&bp.strand2==s1) )
				return true;
		}
		return false;
	}

	
	/**
	 * A pair of strands specified by two strand-indices (i.e. indices in the 
	 * <code>secondaryStructure.getStrands()</code>-array).
	 */
	public static class StrandPair{
		public int strand1, strand2;
		public final boolean parallel; 
		private StrandPair(int strand1, int strand2, boolean parallel){
			this.strand1 = Math.min(strand1,strand2);
			this.strand2 = Math.max(strand1,strand2);
			this.parallel = parallel;
		}
		boolean contains(int strand){
			return strand1==strand || strand2==strand;
		}
	}
}


