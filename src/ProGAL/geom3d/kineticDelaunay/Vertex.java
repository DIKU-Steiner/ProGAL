package ProGAL.geom3d.kineticDelaunay;

import ProGAL.geom3d.Point;

public class Vertex extends Point implements Comparable<Vertex>{
	private static final long serialVersionUID = 1L;

	final int index;
	private static int indexCounter = -1;
	
	public Vertex(Point p) {
		super(p);
		this.index = indexCounter++;
	}

	public int compareTo(Vertex arg0) {
		return index-arg0.index;
	}
	
	public String toString(){
		return Integer.toString(index);
	}
}
