package ProGAL.proteins.beltaStructure.bnb;

import ProGAL.proteins.beltaStructure.bnb.lowerBounds.LowerBound;
import ProGAL.proteins.structure.AminoAcidChain;

public class CN_ClashLowerBound implements LowerBound {
	private double[] contactNumbers;
	
	public CN_ClashLowerBound(double[] contactNumbers){
		this.contactNumbers = contactNumbers;
	}
	@Override
	public double lowerBound(BnBNode n, BnBSolver solver) {
		AminoAcidChain chain = solver.chain;
		
		
		
		// TODO Auto-generated method stub
		return 0;
	}

}
