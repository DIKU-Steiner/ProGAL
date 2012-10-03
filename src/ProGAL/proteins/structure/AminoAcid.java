package ProGAL.proteins.structure;

/**
 * An amino acid containing the following properties:
 * <ul>
 * <li>Pointers to all atoms within this amino acid</li>
 * <li>a pointer to the <code>AminoAcidChain</code> that contains this amino acid</li>
 * <li>the type of the amino acid (Alanine, arginine etc.)</li>
 * <li>an index specifying this amino acids placement in the 
 * <code>aminoAcidChain().aminoAcids()</code>-array (for convenience)</li>
 * </ul>
 * @author R.Fonseca
 */
public class AminoAcid {
	protected Atom[] atoms;
	private final char type;
	
	protected AminoAcidChain chain;
	protected int index; 
	
//	public AminoAcid(char type){
//		this(type,HeavyAtomAminoAcidGenerator.generateAtoms(type));
//	}

	public AminoAcid(char type, Atom[] atoms){
		this.type = type;
		this.atoms = atoms;
		for(int a=0;a<atoms.length;a++){
			atoms[a].aminoAcid = this;
			atoms[a].index = a;
		}
	}

	public Atom[] atoms(){ return atoms; }

	public AminoAcid clone(){
		Atom[] retAtoms = new Atom[atoms.length];
		for(int a=0;a<atoms.length;a++) retAtoms[a] = atoms[a].clone();
		AminoAcid ret = new AminoAcid(type, retAtoms);
		return ret;
	}
	
	public AminoAcidChain aminoAcidChain(){ return chain; }
	
	public int index(){ return index; }

	static long stop = 10000000;
	static long count = 0;
	
	public Atom atom(String name){
//		if((count++%stop)==0){
//			Thread.dumpStack();
//		}
		String ucName = name.toUpperCase();
		int id = ucName.hashCode();
		for(Atom a: atoms)
			if(a.id==id) return a;// && a.name().equalsIgnoreCase(ucName)) return a;
		throw new RuntimeException(name+" not found in "+typeThreeLetter());
	}
	
	public Atom atom(int atomNumber){
		return atoms[atomNumber];
	}
	
	public char type(){
		return type;
	}
	
	public String typeName(){
		switch(type){
		case 'A': return "Alanine";
		case 'R': return "Arginine";
		case 'N': return "Asparagine";
		case 'D': return "Aspartic acid";
		case 'C': return "Cysteine";
		case 'E': return "Glutamic acid";
		case 'Q': return "Glutamine";
		case 'G': return "Glycine";
		case 'H': return "Histidine";
		case 'I': return "Isoleucine";
		case 'L': return "Leucine";
		case 'K': return "Lysine";
		case 'M': return "Methionine";
		case 'F': return "Phenylalanine";
		case 'P': return "Proline";
		case 'S': return "Serine";
		case 'T': return "Threonine";
		case 'W': return "Tryptophan";
		case 'Y': return "Tyrosine";
		case 'V': return "Valine";
		case 'U': return "Selenocysteine";
		}
		return String.format("Unknown (%c)",type);
	}
	public String typeThreeLetter(){
		switch(type){
		case 'A': return "Ala";
		case 'R': return "Arg";
		case 'N': return "Asn";
		case 'D': return "Asp";
		case 'C': return "Cys";
		case 'E': return "Glu";
		case 'Q': return "Gln";
		case 'G': return "Gly";
		case 'H': return "His";
		case 'I': return "Ile";
		case 'L': return "Leu";
		case 'K': return "Lys";
		case 'M': return "Met";
		case 'F': return "Phe";
		case 'P': return "Pro";
		case 'S': return "Ser";
		case 'T': return "Thr";
		case 'W': return "Trp";
		case 'Y': return "Tyr";
		case 'V': return "Val";
		case 'U': return "Sec";
		}
		return "Unk";

	}
}
