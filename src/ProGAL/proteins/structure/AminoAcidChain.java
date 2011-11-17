package ProGAL.proteins.structure;

/**
 * A chain of amino acids. This class together with AminoAcid, Atom and CBond are meant as 
 * convenient methods to maintain and traverse a protein structure. There are no methods 
 * for changing a protein structure (such classes could appropriately extend AminoAcidChain).
 * Furthermore the pointers between amino acids and atoms are not meant to be changed. 
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
		this.aminoAcids = new AminoAcid[sequence.length()];
		for(int r=0;r<sequence.length();r++){
			aminoAcids[r] = new AminoAcid(sequence.charAt(r));
			aminoAcids[r].chain = this;
			aminoAcids[r].index = r;
		}
		this.covalentBonds = HeavyAtomAminoAcidGenerator.generateBonds(aminoAcids);
	}

	public AminoAcidChain(AminoAcid[] aminoAcids){
		this.aminoAcids = aminoAcids;
		for(int r=0;r<aminoAcids.length;r++){
			aminoAcids[r].chain = this;
			aminoAcids[r].index = r;
		}
		this.covalentBonds = HeavyAtomAminoAcidGenerator.generateBonds(aminoAcids);
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
