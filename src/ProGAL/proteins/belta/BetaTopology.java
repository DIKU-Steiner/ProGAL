package ProGAL.proteins.belta;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import ProGAL.math.Matrix;

/**
 * Representation of a beta-topology. The pairingMatrix is initialized to a zero-matrix, 
 * but can freely be edited using the setPaired-method. The <code>isValid</code> method 
 * will tell if the pairing matrix is valid (i.e. no strand paired with itself, no two 
 * strands paired both parallel and anti-parallel and all strands paired at least once 
 * and at most twice with other strands).
 * 
 * To find out which strands are paired either read the pairing matrix or use 
 * <code>getSheets</code> to retrieve a list of sheets that can be individually managed.  
 * @author R.Fonseca
 */
public class BetaTopology {

	/** The secondary structure that the beta-topology is related to */
	public final SecondaryStructure secondaryStructure;

	/** The pairing matrix that indicates which strands are paired. Columns and rows in the 
	 * pairing matrix correspond to indices in the 
	 * <code>secondaryStructure.getStrands()</code>-array. If an entry is nonzero in the 
	 * upper-triangle it indicates an antiparallel pairing and in the lower-triangle it 
	 * indicates a parallel pairing.  */
	public final Matrix pairingMatrix;
	
	/** The number of strands in the secondary structure */
	public final int N;
	
	

	/** 
	 * Constructs a beta-topology associated with the specified secondary 
	 * structure but with no pairings. All fields are initialized but the 
	 * pairing matrix has only zero-entries.
	 */
	public BetaTopology(SecondaryStructure ss){
		this.secondaryStructure = ss;
		this.N = ss.getStrands().length;
		this.pairingMatrix = new Matrix(N,N);
	}

	/** 
	 * Constructs a beta-topology associated with the specified secondary 
	 * structure and pairings. All fields are initialized.
	 */
	public BetaTopology(SecondaryStructure ss, boolean[][] pairs){
		this(ss);
		
		for(int r=0;r<N;r++)
			for(int c=0;c<N;c++)
				if(pairs[r][c]) pairingMatrix.set(r, c, 1);
	}

	/** Returns true if i and j are paired. */
	public boolean pair(int i, int j){
		return pairingMatrix.get(i, j)!=0;
	}

	/** Set strands i and j to be paired. */
	public void setPaired(int i, int j){
		pairingMatrix.set(i, j, 1);
	}
	
	/** Set strands i and j to be unpaired. */
	public void setNotPaired(int i, int j){
		pairingMatrix.set(i, j, 0);
	}
	
	
	/** 
	 * Checks if this beta-topology is valid. A valid topology has
	 * <ul><li>no strand paired with itself</li>
	 * <li>no two strands paired both parallel and anti-parallel</li>
	 * <li>all strands paired at least once and at most twice with other strands</li></ul> 
	 */
	public boolean isValid(){
		//Check for double pairings (both anti-parallel and parallel)
		//No self-pairings (diagonal)
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				if(pair(i,j) && pair(j,i)) return false;
			}
		}

		//At least one pairing per strand
		//At most two pairings per strand
		for(int i=0;i<N;i++) {
			int sum = 0;
			for(int j=0;j<N;j++) {
				if(pair(i,j) || pair(j,i)) sum++;
			}
			if(sum<1 || sum>2) return false;
		}
		return true;
	}
	
	/** 
	 * Construct a list of sheet-topologies from the beta-topology.
	 */
	public List<SheetTopology> getSheets(){
		if(!isValid()) throw new RuntimeException("Topology must be valid");
		List<SheetTopology> sheets = new LinkedList<SheetTopology>();
		for(int s=0;s<N;s++){
			//Check if strand s is in a sheet already
			boolean strandInSheets = false;
			for(SheetTopology sh: sheets){ if(sh.containsStrand(s)) {strandInSheets = true; break;} }

			//If not construct a sheet
			if(!strandInSheets) sheets.add(new SheetTopology(this,s));
		}
		return sheets;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("BetaTopology:\n");
		for(int i=0;i<N;i++){
			sb.append("> ");
			for(int j=0;j<N;j++){
				if(pair(i,j)) 	sb.append("1");
				else			sb.append("0");
			}
			if(i!=N-1)	sb.append('\n');
		}
		return sb.toString();
	}

	/**
	 * Determine if this beta-topology matches some other topology. Two topologies, t1 and t2, 
	 * match each other iff t1 respects t2 and t2 respects t1. 
	 */
	public boolean matches(BetaTopology bTop){
		if(bTop.N!=N) return false;
		if(!secondaryStructure.matches(bTop.secondaryStructure)) return false;
		
		for(int r=0;r<N;r++){
			for(int c=0;c<N;c++){
				if(pair(r,c)!=bTop.pair(r,c)) return false;
			}
		}
		return true;
	}

	/**
	 * Determine if this beta-topology respects some other topology. A topology, t1, respects 
	 * another topology, t2, iff <code>t1.secondaryStructure.respects(t2.secondaryStructure)</code> 
	 * and if a pairing of two strands in t1 implies a pairing of two overlapping strands in t2.  
	 */
	public boolean respects(BetaTopology bTop){
		if(bTop.N<N) return false;
		List<int[]> strandAlignment = SecondaryStructure.getStrandPairing(secondaryStructure.getStrands(), bTop.secondaryStructure.getStrands());
		if(strandAlignment==null) return false;//Also implies that secondaryStructures doesnt respect bTop.secondaryStructure
		
		List<Integer> strandsInSecond = new ArrayList<Integer>();
		for(int[] p: strandAlignment) strandsInSecond.add(p[1]);
		Collections.sort(strandsInSecond);

		Matrix minor = bTop.pairingMatrix;
		for(int i=bTop.N-1;i>=0;i--){
			if(!strandsInSecond.contains(i)) //minor = minor(minor, i,i);
				minor = minor.minor(i, i);
		}
		
		for(int r=0;r<N;r++){
			for(int c=0;c<N;c++){
				if(pair(r,c) && minor.get(r,c)==0) return false;
			}
		}
		return true;
	}


}
