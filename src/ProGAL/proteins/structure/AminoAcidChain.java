package ProGAL.proteins.structure;

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
