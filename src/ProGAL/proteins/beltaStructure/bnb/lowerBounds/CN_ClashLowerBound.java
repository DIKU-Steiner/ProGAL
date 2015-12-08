package ProGAL.proteins.beltaStructure.bnb.lowerBounds;

import ProGAL.proteins.beltaStructure.bnb.BnBNode;
import ProGAL.proteins.beltaStructure.bnb.BnBSolver;
import ProGAL.proteins.structure.AminoAcidChain;

public class CN_ClashLowerBound implements LowerBound {
	private double[] contactNumbers;
	
	public CN_ClashLowerBound(double[] contactNumbers){
		this.contactNumbers = contactNumbers;
	}
	@Override
	public double lowerBound(BnBNode n, BnBSolver solver) {
		AminoAcidChain chain = solver.getChain();
		
		
		
		// TODO Auto-generated method stub
		return 0;
	}

}
