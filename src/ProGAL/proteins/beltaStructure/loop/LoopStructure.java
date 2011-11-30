package ProGAL.proteins.beltaStructure.loop;

import ProGAL.geom3d.Point;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.beltaStructure.sheetLoop.PartialStructure;
import ProGAL.proteins.chainTree.ChainTree;
import ProGAL.proteins.dunbrack.RamachandranDistribution;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.Atom;

public class LoopStructure implements PartialStructure{
	public final SecondaryStructure secondaryStructure;
	public final SSSegment segment1, segment2;
	public Atom[] targetAtoms;
	private ChainTree chaintree;


	public LoopStructure(SecondaryStructure ss, int seg1, int seg2, Atom[] targetAtoms){
		this.secondaryStructure = ss;
		this.segment1 = ss.segments[Math.min(seg1, seg2)];
		this.segment2 = ss.segments[Math.max(seg1, seg2)];
		this.targetAtoms = targetAtoms;
		this.chaintree = new ChainTree(ss.primaryStructure, segment1.start, segment2.end);

		//Lock helices
		for(int s=seg1+1;s<seg2;s++){
			if(ss.segments[s].type==SSType.HELIX){
				for(int r=ss.segments[s].start+1;r<ss.segments[s].end-1;r++){
					chaintree.setLocked(r-segment1.start, 0);
					chaintree.setLocked(r-segment1.start, 1);
					chaintree.setLocked(r-segment1.start, 2);
					chaintree.setTorsionAngle(r-segment1.start, 0, -60*Math.PI/180);
					chaintree.setTorsionAngle(r-segment1.start, 1, -30*Math.PI/180);
				}
			}
		}
	}

	public LoopStructure(SecondaryStructure ss, int seg1, int seg2){
		this(ss,seg1,seg2, new Atom[]{});
	}

	public void setFirstTransform(Matrix m){
		chaintree.setFirstTransformation(m);
	}

	public String toString(){
		return String.format("LoopStructure[%d-%d]",segment1.start, segment2.end-1);
	}

	/** 
	 * Indicates if the loop is closed, ie. if the end of the loop matches up with the target atoms.
	 */
	public boolean isClosed(){
		int loopLength = segment2.end-segment1.start;
		double sqDistSum = 0;
		for(int i=0;i<targetAtoms.length;i++){
			sqDistSum += chaintree.getBackboneAtom(loopLength, i).distanceSquared(targetAtoms[i]);
		}
		double rms = Math.sqrt(sqDistSum/(targetAtoms.length>0?targetAtoms.length:1)); 
		return rms<0.4;
	}

	public void enforceClosureCCD(){
		chaintree.closeCCD(targetAtoms, 20);
	}

	public void enforceClosureAnalytically(){
		//TODO
	}

	public void enforceClosureJacobian(){
		//TODO
	}

	public void rebuildCCD(){
		resampleFromRama();
		enforceClosureCCD();
	}

	public void rebuildAnalytically(){
		resampleFromRama();
		//TODO
	}

	public void rebuildJacobian(){
		resampleFromRama();
		//TODO
	}

	public void rebuildACO(){
		//TODO
	}

	private void resampleFromRama(){
		RamachandranDistribution distr = RamachandranDistribution.getDistribution();
		
		int length = segment2.end-segment1.start;
		for(int r=0;r<length;r++){
			if(chaintree.isLocked(r, 0) || chaintree.isLocked(r, 1)) continue;

//			double rand = Randomization.randBetween(0.0, 1.0), phi,psi;
//			if(rand<0.4){//Beta area
//				phi = Randomization.randBetween(-170.0, -50.0)*Math.PI/180;
//				psi = Randomization.randBetween( 100.0, 160.0)*Math.PI/180;
//			}else if(rand<0.8){//Alpha area
//				phi = Randomization.randBetween(-100.0, -50.0)*Math.PI/180;
//				psi = Randomization.randBetween( -40.0,  10.0)*Math.PI/180;
//			}else{//Left helix area
//				phi = Randomization.randBetween(50.0, -70.0)*Math.PI/180;
//				psi = Randomization.randBetween(20.0,  60.0)*Math.PI/180;
//			}
//			chaintree.setTorsionAngle(r, 0, phi);
//			chaintree.setTorsionAngle(r, 1, psi);
			int type = distr.getTypes().indexOf(secondaryStructure.primaryStructure.getThreeLetterType(r+segment1.start).toUpperCase());
			double[] phiPsi = distr.samplePhiPsi(type);
			chaintree.setTorsionAngle(r, 0, phiPsi[0]);
			chaintree.setTorsionAngle(r, 1, phiPsi[1]);

		}
	}

	@Override
	public void updateAtoms(AminoAcidChain chain) {
		Point[] atomCoords = chaintree.getAllBackboneAtoms();
		int c = 0;
		for(int r=segment1.start;r<segment2.end;r++){
			for(int a=0;a<3;a++){
				chain.atom(r, a).set(atomCoords[c++]);
			}
		}
	}


}
