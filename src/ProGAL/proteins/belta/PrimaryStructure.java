package ProGAL.proteins.belta;

public class PrimaryStructure {
	/** The sequence of a protein chain */
	public final String sequence;

	public PrimaryStructure(String seq){
		this.sequence = seq;
	}
	
	public String toString(){
		return sequence;
	}
	
	public String getThreeLetterType(int res){
		switch(sequence.charAt(res)){
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
