package ProGAL.proteins.structure.generators;

import ProGAL.proteins.structure.AminoAcid;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.CBond;

public interface AtomGenerator {
	public Atom[] generateAtoms(char type);
	public CBond[] generateBonds(AminoAcid[] aminoAcids);
}
