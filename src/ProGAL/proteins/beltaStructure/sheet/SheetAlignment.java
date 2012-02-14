package ProGAL.proteins.beltaStructure.sheet;

import java.util.Arrays; 

import ProGAL.proteins.PDBFile.AtomRecord;
import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.PDBFile;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.belta.SheetTopology.StrandPair;

public class SheetAlignment {
	public SheetTopology sTop;
	
	/**
	 * The entire alignment is specified by <code>sTop.strandPairs.size()</code> residue indices. There is an entry 
	 * for each StrandPair in <code>sTop.strandPairs</code>.
	 * Each entry is a residue-index indicating which residue in the latter strand is paired to the first residue
	 * in the first strand. Notice that this residue might be outside the actual strand. TODO: Explain better.
	 */
	public int[] alignmentPairs;
	
	/** Construct an alignment of strands in the given topology. Strands are 
	 * aligned such that their centers align.*/
	public SheetAlignment(SheetTopology sTop){
		this.sTop = sTop;
		this.alignmentPairs = new int[sTop.strandPairs.size()];
		
		SSSegment[] strands = sTop.secondaryStructure.getStrands();
		for(StrandPair sp: sTop.strandPairs){
			SSSegment strand1 = strands[sp.strand1];
			SSSegment strand2 = strands[sp.strand2];
			setAligned( (strand1.end+strand1.start)/2,(strand2.end+strand2.start)/2 );
		}
	}
	
	public SheetAlignment(SheetTopology sTop, PDBFile f){
		this(sTop);
		for(StrandPair bp: sTop.strandPairs){
			SSSegment s1 = sTop.secondaryStructure.getStrands()[bp.strand1];
			SSSegment s2 = sTop.secondaryStructure.getStrands()[bp.strand2];
			double minDist = Double.POSITIVE_INFINITY;
			int minRes1 = 0, minRes2 = 0;
			for(int r1=s1.start;r1<s1.end;r1++){
				AtomRecord n1 = f.getAtom(r1, "N");
				AtomRecord o1 = f.getAtom(r1, "O"); 
				for(int r2=s2.start;r2<s2.end;r2++){
					AtomRecord n2 = f.getAtom(r2, "N");
					AtomRecord o2 = f.getAtom(r2, "O");
					double dist = Math.min(n1.coords.distance(o2.coords), o1.coords.distance(n2.coords));
					if(dist<minDist){
						minDist = dist;
						minRes1 = r1;
						minRes2 = r2;
					}
				}
			}
			this.setAligned(minRes1, minRes2);
		}
		
	}
	
	public static void main(String[] args){
		SecondaryStructure ss = new SecondaryStructure("  EEE  EEEE  EEEEE");
		BetaTopology bt = new BetaTopology(ss);
		bt.setPaired(0, 1);
		bt.setPaired(2, 1);
		SheetTopology st = bt.getSheets().get(0);
		SheetAlignment sa = new SheetAlignment(st);
		System.out.println(Arrays.toString(sa.alignmentPairs));
		//Expects 9(10) and 13(14)
		sa.setAligned(10, 17);
		System.out.println(Arrays.toString(sa.alignmentPairs));
	}
	
	/** 
	 * Aligns the strand with residue res1 with another strand containing 
	 * res2 such that res1 and res2 are aligned. 
	 */
	public void setAligned(int res1, int res2){
		if(res1>res2){//Swap
			int tmp = res1; 	res1 = res2;	res2 = tmp;
		}
		
		//Find strand-indices corresponding to the two residues
		SSSegment[] strands = sTop.secondaryStructure.getStrands();
		int strand1 = indexOf(strands, sTop.secondaryStructure.getSegmentContainingResidue(res1));
		int strand2 = indexOf(strands, sTop.secondaryStructure.getSegmentContainingResidue(res2));

		
		if(strand1<0 || strand2<0)
			throw new RuntimeException(String.format("Residue %d or %d is not in a strand", res1, res2));
		
		//Find index of strand-pair
		int pairIdx = -1;
		StrandPair pair = null;
		for(StrandPair sp: sTop.strandPairs){
			if( 	(sp.strand1==strand1 && sp.strand2==strand2) || 
					(sp.strand2==strand1 && sp.strand1==strand2)	){
				pairIdx = sTop.strandPairs.indexOf(sp);
				pair = sp;
				break;
			}
		}
		
		//Fail if the residues are not paired
		if(pair==null) throw new RuntimeException(String.format("Strands %d and %d are not paired in this topology",strand1, strand2));
		
//		System.out.printf("res1 %d, res2 %d, strand1 %d, strand2 %d, strand1.start %d, strand2.start %d (%s)\n",res1,res2,strand1, strand2, strands[strand1].start,strands[strand2].start,pair.strand2==strand2?"no reverse":"reverse");
		
		//Update alignmentPairs
		if(pair.strand2==strand2){
			if(pair.parallel)	alignmentPairs[pairIdx] = res2 - (res1-strands[strand1].start);
			else				alignmentPairs[pairIdx] = res2 + (res1-strands[strand1].start);
		}else{
			if(pair.parallel)	alignmentPairs[pairIdx] = res1 - (res2-strands[strand2].start);
			else				alignmentPairs[pairIdx] = res1 + (res2-strands[strand2].start);
		}
		
	}
	
	private int indexOf(Object[] arr, Object entry){
		for(int i=0;i<arr.length;i++){
			if(entry==arr[i]) return i;
		}
		return -1;
	}
}
