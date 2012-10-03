package ProGAL.proteins.belta;

import java.util.ArrayList;
import java.util.List;

/** 
 * Representation of a secondary structure assignment. The constructor reads
 * a string representation of an assignment and parses it into segments. The 
 * start and end-index of a segment, as well as the segments type can be 
 * accessed in the following way
 * <pre>
 * SecondaryStructure ss = new SecondaryStructure("   HHHHHHHH  EEE  EEE   ");
 * Segment seg = ss.segments[1];
 * System.out.printf("Segment %d: start %d, end %d, type: %s\n", seg.segmentIndex, seg.start, seg.end, seg.type );  
 * </pre>
 * Three convenient methods are also supplied to return all helix, strand or 
 * coil segments separate from each other. 
 * @author R.Fonseca
 */
public class SecondaryStructure {
	/** 
	 * A secondary structure segment of either helix, coil or strand. The only function 
	 * is to hold indices of start and end residues, the length, the index of the segment 
	 * in the secondary structure and the type.
	 */
	public static class SSSegment{
		
		/** The residue-index of the first amino acid in this segment */ 
		public final int start;
		
		/** The residue-index of the last amino acid in this segment PLUS ONE. */
		public final int end;
		
		/**  */
		public final int length;
		public final int segmentIndex;
		public final SSType type;
	
		private SSSegment(int start, int end, int segmentIndex, SSType type){
			this.start = start;
			this.end = end;
			this.length = end-start;
			this.segmentIndex = segmentIndex;
			this.type = type;
		}
		
		public boolean overlaps(SSSegment seg ){
			if(seg.start>=end) return false;
			if(seg.end<=start) return false;
			return true;
		}
		public String toString(){
			return String.format("SSSegment[%c,%d,%d]", type.toChar(), start, end);
		}

		public boolean contains(int res) {
			return start<=res && res<end;
		}

		public int midpoint() {
			return (end+start-1)/2;
		}
	}
	
	/** Secondary structure segments (helix, coil or strand) */
	public final SSSegment[] segments;

	/** The primary structure related to this secondary structure. Might be null. */
	public PrimaryStructure primaryStructure;
	
	/** 
	 * Construct a secondary structure object from a string representation. 
	 * The primary structure is stored as well.  
	 */
	public SecondaryStructure(PrimaryStructure ps, String ssString){
		this(ssString);
		primaryStructure = ps;
	}
	
	/** 
	 * Construct a secondary structure object from a string representation. 
	 * The string representation may be both the DSSP (helix: H,G; strand: E;
	 * coil: B,I,T,S) or the normal three-class representation (helix: H; 
	 * strand: E; coil: C or space). 
	 */
	public SecondaryStructure(String ssString){
		//Construct segments
		int tmpStart = 0, count = 0;
		char[] tmpSS = ssString.toCharArray();
		ArrayList<SSSegment> tmpSegments = new ArrayList<SSSegment>();
		for(int i=1;i<tmpSS.length;i++){
			if( type(tmpSS[i]) != type(tmpSS[i-1]) ){ 
				tmpSegments.add(  new SSSegment(tmpStart, i, count, type(tmpSS[i-1]))  );
				tmpStart = i;
				count++;
			}
		}
		if(tmpSS.length>0) 
			tmpSegments.add(  new SSSegment(tmpStart, tmpSS.length, count, type(tmpSS[tmpSS.length-1]))  );
		
		segments = new SSSegment[tmpSegments.size()];
		tmpSegments.toArray(segments);

		//Collect strands
		tmpSegments.clear();
		for(SSSegment seg: segments) if(seg.type==SSType.STRAND) tmpSegments.add(seg);
		strands = new SSSegment[tmpSegments.size()];
		tmpSegments.toArray(strands);
		
		//Collect helices
		tmpSegments.clear();
		for(SSSegment seg: segments) if(seg.type==SSType.HELIX) tmpSegments.add(seg);
		helices = new SSSegment[tmpSegments.size()];
		tmpSegments.toArray(helices);
		
		//Collect coils
		tmpSegments.clear();
		for(SSSegment seg: segments) if(seg.type==SSType.COIL) tmpSegments.add(seg);
		coils = new SSSegment[tmpSegments.size()];
		tmpSegments.toArray(coils);
		
	}
	
	private SSType type(char c){
		switch(Character.toUpperCase(c)){
		case ' ':
		case 'I':
		case 'T':
		case 'S':
		case 'C':
		case 'B':
		case '.': return SSType.COIL;
		case 'E': return SSType.STRAND;
		case 'G':
		case 'H': return SSType.HELIX;
		}
		throw new RuntimeException("Unknown SS-type: "+c);
	}
	
	private final SSSegment[] strands, helices, coils;
	
	/** Return array containing all segments, s, where <code>s.type==SSType.STRAND</code>. These are in sequential order.*/
	public SSSegment[] getStrands(){ return strands; }
	
	/** Return array containing all segments, s, where <code>s.type==SSType.HELIX</code>. These are in sequential order. */
	public SSSegment[] getHelices(){ return helices; }
	
	/** Return array containing all segments, s, where <code>s.type==SSType.COIL</code>. These are in sequential order. */
	public SSSegment[] getCoils(){   return coils;	}

	/** 
	 * Determine if this secondary structure matches ss. Two secondary structures, s1 and s2,
	 * match iff s1 respects s2 and s2 respects s1. 
	 */
	public boolean matches(SecondaryStructure ss){
		return this.respects(ss) && ss.respects(this);
	}
	
	/**
	 * Determine if this secondary structure respects ss. One secondary structure, s1, respects 
	 * another, s2, iff there exists a one-to-one pairing of every strand in s1 to a subset of 
	 * strands in s2 such that each pair of strands overlap.
	 */
	public boolean respects(SecondaryStructure ss){
		if(strands.length==0) return true; 
		List<int[]> strandPairing = getStrandPairing(strands, ss.strands);
		return strandPairing!=null;
	}
	
	/** 
	 * If possible, find a list of index-pairs such that for each member, (i1, i2), 
	 * <code>strands1[i1].overlaps(strands2[i2])</code>. If every member of strands1  
	 * can not be paired with a member of strands2 then null is returned. 
	 */
	static List<int[]> getStrandPairing(SSSegment[] strands1, SSSegment[] strands2){
		if(strands1.length>strands2.length) return null;
		List<int[]> ret = new ArrayList<int[]>();
		if(strands1.length==0) return ret;
		boolean[][] overlaps = new boolean[strands1.length][strands2.length];
		for(int r=0;r<strands1.length;r++)
			for(int c=0;c<strands2.length;c++)
				overlaps[r][c] = strands1[r].overlaps(strands2[c]);
		boolean[][] pairing = findPairing(overlaps,new boolean[strands1.length][strands2.length], 0,0);
		if(pairing==null) return null;

		for(int r=0;r<strands1.length;r++)
			for(int c=0;c<strands2.length;c++)
				if(pairing[r][c]) ret.add(new int[]{r,c});
		return ret;
	}
	
	/** 
	 * Recursive brute-force method to determine a pairing of every row with a subset of columns given a 
	 * matrix of potential pairs. Assumes that <code>overlaps</code> has no more rows than columns. 
	 */
	private static boolean[][] findPairing(boolean[][] overlaps, boolean[][] curPairing, int r, int c){

		if(c==0 && r>0){
			boolean lastRowPaired = false;
			for(int i=0;i<overlaps[0].length;i++) if(curPairing[r-1][i]){ lastRowPaired = true; break; }
			if(!lastRowPaired) return null;
			if(r==overlaps.length) return curPairing;
		}
		
		int nextR = r;
		int nextC = (c+1)%overlaps[0].length;
		if(nextC==0) nextR++;
		
		if(!overlaps[r][c]) {
			return findPairing(overlaps, curPairing, nextR, nextC);
		}else{
			//Is either already chosen
			boolean taken = false; 
			for(int i=0;i<c;i++) if(curPairing[r][i]) {taken = true; break;}
			if(taken) return findPairing(overlaps,curPairing,r+1, 0);
			for(int i=0;i<r;i++) if(curPairing[i][c]) {taken = true; break;}
			if(taken) return findPairing(overlaps,curPairing,nextR, nextC);
			
			//Its available .. try to take it and try not to
			curPairing[r][c] = true;
			boolean[][] res = findPairing(overlaps, curPairing, r+1,0);
			if(res!=null) return res;
			curPairing[r][c] = false;
			return findPairing(overlaps,curPairing, nextR, nextC);
			
		}
	}
	
	/**
	 * TODO: Comment and test
	 */
	public SSSegment getSegmentContainingResidue(int res){
		for(SSSegment seg: segments) 
			if(seg.start<=res && seg.end>res) return seg;
		return null;
	}

	/**
	 * TODO: Comment and test
	 */
	public SSType getType(int r) {
		return getSegmentContainingResidue(r).type;
	}
	/**
	 * TODO: Comment
	 */
	public String toString(){
		StringBuilder ret = new StringBuilder();
		for(SSSegment seg: segments) for(int i=0;i<seg.length;i++) ret.append(seg.type.toChar());
		return ret.toString();
	}
	
	
}
