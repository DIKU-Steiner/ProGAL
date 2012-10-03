package ProGAL.geomNd.complex;

import ProGAL.geomNd.Point;

public class Vertex extends Point {
	private static final long serialVersionUID = 1L;
	
	private boolean degenerate=false;	
	public  enum DegenerateCase { ONFACE, ONEDGE };	
	private final boolean bigPoint;
	
	public boolean isBigpoint() {
		return bigPoint;
	}	

	public boolean isDegenerate() {
		return degenerate;
	}

	public void setDegenerate(boolean degenerate) {
		this.degenerate = degenerate;
	}
	
	public Vertex(Point p){
		this(p, false);
	}
	
	public Vertex(Point p, boolean bigpoint){ 	
		super(p);
		setDegenerate(false);
		this.bigPoint = bigpoint;
	}

	public String toString(){ return toString(2); }
	
	public String toString(int dec) {
		String tstr = "Vertex[";
		for (int i = 0; i < coords.length; ++i) {
			tstr += String.format("%.4f", coords[i]);
			if (i != coords.length-1) {
				tstr += ", ";
			}
		}
		tstr += "]";
		return tstr;
	}
}