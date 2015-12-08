package ProGAL.dataStructures;

/**
 * Adapted from http://stackoverflow.com/questions/521171/a-java-collection-of-value-pairs-tuples
 */ 
public class Pair<L,R> {

	public final L fst;
	public final R snd;

	public Pair(L first, R second) {
		this.fst = first;
		this.snd = second;
	}

	@Override
	public int hashCode() { return fst.hashCode() ^ snd.hashCode(); }

	@Override
	public boolean equals(Object o) {
		if (o == null) return false;
		if (!(o instanceof Pair)) return false;

		@SuppressWarnings("rawtypes")
		Pair pairo = (Pair) o;
		return this.fst.equals(pairo.fst) &&
				this.snd.equals(pairo.snd);
	}
	
	@Override
	public String toString(){
		return String.format("(%s,%s)",fst.toString(), snd.toString());
	}

}