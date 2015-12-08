package ProGAL.proteins.belta;

/** Types of secondary structure */
public enum SSType{ COIL, HELIX, STRAND;

	public char toChar(){
		switch(this){
		case COIL: return 'C';
		case HELIX: return 'H';
		case STRAND: return 'E';
		}
		return '?';
	}
}

