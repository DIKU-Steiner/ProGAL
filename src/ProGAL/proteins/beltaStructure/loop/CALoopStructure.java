package ProGAL.proteins.beltaStructure.loop;

import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.beltaStructure.sheetLoop.PartialStructure;
import ProGAL.proteins.chainTree.ChainTree;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.Atom;

public class CALoopStructure implements PartialStructure{
	public final SecondaryStructure secondaryStructure;
	public final SSSegment segment1, segment2;
	public SegmentStructure[] segments; 


	public CALoopStructure(SecondaryStructure ss, int seg1, int seg2){
		this.secondaryStructure = ss;
		this.segment1 = ss.segments[Math.min(seg1, seg2)];
		this.segment2 = ss.segments[Math.max(seg1, seg2)];
	}

	
	@Override
	public void updateAtoms(AminoAcidChain chain) {
		// TODO Auto-generated method stub
		
	}

	
	public class SegmentStructure implements PartialStructure{
		private SSSegment seg;
		
		SegmentStructure(SSSegment seg){
			this.seg = seg;
		}
		
		public void updateAtoms(AminoAcidChain chain) {
			
		}
		
	}
	
	public static void main(String[] args){

	}
}
