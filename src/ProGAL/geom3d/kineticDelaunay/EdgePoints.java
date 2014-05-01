package ProGAL.geom3d.kineticDelaunay;

public class EdgePoints {
	Vertex p0;
	Vertex p1;
	
	public EdgePoints(Vertex p0, Vertex p1) {
		if (p0.getId()<p1.getId()) {
			this.p0 = p0;
			this.p1 = p1;
		} else {
			this.p1 = p0;
			this.p0 = p1;
		}
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof EdgePoints)) return false;
		return ((EdgePoints)o).p0.equals(p0)&&((EdgePoints)o).p1.equals(p1);
	}
	
	public int hashCode() {
		return p0.getId()+p1.getId()*2;
	}
}
