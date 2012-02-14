package ProGAL.proteins.beltaStructure.bnb;

import java.util.List;

/** A part of a protein chain. This could be either a beta-sheet or a coil/helix segment */
public interface Branchable {
	public void setStructure(int s);
	public int getStructures();
	public List<Integer> definedResidues();
}
