package ProGAL.proteins.beltaStructure;

import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;

public class ChainTreeLoop {
	public final SecondaryStructure secondaryStructure;
	public final SSSegment segment1, segment2;
	
	public ChainTreeLoop(SecondaryStructure ss, int seg1, int seg2){
		this.secondaryStructure = ss;
		this.segment1 = ss.segments[Math.min(seg1, seg2)];
		this.segment2 = ss.segments[Math.max(seg1, seg2)];
	}

	public boolean isClosed(){
		//TODO
		return false;
	}

	public void setPsi(int res, double newValue){ 	setTorsion(res, 0, newValue);	}
	public void setPhi(int res, double newValue){ 	setTorsion(res, 1, newValue);	}
	public void setOmega(int res, double newValue){	setTorsion(res, 2, newValue);	}
	
	public void setTorsion(int res, int torsion, double newValue){
		//TODO
	}
	
	public void enforceClosureCCD(){
		//TODO
	}
	
	public void enforceClosureAnalytically(){
		//TODO
	}
	
	public void enforceClosureJacobian(){
		//TODO
	}
	
	public void rebuildCCD(){
		//TODO
	}
	
	public void rebuildAnalytically(){
		//TODO
	}
	
	public void rebuildJacobian(){
		//TODO
	}
	
	public void rebuildACO(){
		//TODO
	}
}
