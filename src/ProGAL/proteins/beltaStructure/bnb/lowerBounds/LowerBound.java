package ProGAL.proteins.beltaStructure.bnb.lowerBounds;

import ProGAL.proteins.beltaStructure.bnb.BnBNode;
import ProGAL.proteins.beltaStructure.bnb.BnBSolver;

public interface LowerBound {
	double lowerBound(BnBNode n, BnBSolver solver);
}
