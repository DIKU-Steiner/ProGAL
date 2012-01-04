package ProGAL.proteins.structure;

import ProGAL.proteins.structure.generators.AtomGenerator;
import ProGAL.proteins.structure.generators.HeavyAtomAminoAcidGenerator;

/**
 * A chain of amino acids. This class together with AminoAcid, Atom and CBond are meant as 
 * convenient methods to maintain and traverse a protein structure. There are no methods 
 * for changing a protein structure (such classes could appropriately extend AminoAcidChain).
 * Furthermore the pointers between amino acids and atoms are not meant to be changed. 
 * The types of atoms and covalent bonds are defined by an AtomGenerator.  
 * 
 * In a Model-View-Control framework this class together with AminoAcid, Atom and CBond are 
 * meant as the model. 
 *  
 * @author R.Fonseca
 */
public class AminoAcidChain {
	protected AminoAcid[] aminoAcids;
	protected CBond[] covalentBonds;

	public AminoAcidChain(String sequence){
		this(sequence, new HeavyAtomAminoAcidGenerator());
	}
	
	public AminoAcidChain(String sequence, AtomGenerator ag){
		this.aminoAcids = new AminoAcid[sequence.length()];
		for(int r=0;r<sequence.length();r++){
			aminoAcids[r] = new AminoAcid(sequence.charAt(r), ag.generateAtoms(sequence.charAt(r)));
			aminoAcids[r].chain = this;
			aminoAcids[r].index = r;
		}
		this.covalentBonds = ag.generateBonds(aminoAcids);
	}

	public AminoAcid[] aminoAcids(){ return aminoAcids; }

	public AminoAcid aminoAcid(int res){
		return aminoAcids[res];
	}

	public Atom atom(int res, String name){			return aminoAcids[res].atom(name);			}//Very convenient
	public Atom atom(int res, int atomNumber){		return aminoAcids[res].atom(atomNumber);	}//Very fast. Requires no array-traversal

	public Atom[] atoms(){
		int count = 0;
		for(AminoAcid aa: aminoAcids) count+=aa.atoms.length;
		Atom[] ret = new Atom[count];

		count=0;
		for(AminoAcid aa: aminoAcids) 
			for(Atom a: aa.atoms) ret[count++] = a;

		return ret;
	}

	public CBond[] covalentBonds(){
		return covalentBonds;

	}
}
