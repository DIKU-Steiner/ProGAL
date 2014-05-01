package ProGAL.geom3d.kineticDelaunay;

public class TrianglePoints {
	Vertex p0;
	Vertex p1;
	Vertex p2;
	
	public TrianglePoints(Vertex p0, Vertex p1, Vertex p2) {
		this.p0 = p0;
		this.p1 = p1;
		this.p2 = p2;
		sortCorners();
	}
	
	private void sortCorners() {
		int id0 = p0.getId();
		int id1 = p1.getId();
		int id2 = p2.getId();
		Vertex tmp;
		if (id0>id1) {
			if (id2>id1) {
				tmp = p0;
				p0 = p1;
				if (id0>id2) {
					p1 = p2;
					p2 = tmp;
				} else p1 = tmp;
			} else {
				tmp = p0;
				p0 = p2;
				p2 = tmp;
			}
		}
		if (id1>id2) {
			if (id0>id2){
				tmp = p0;
				p0 = p2;
				p2 = p1;
				p1 = tmp;
			} else {
				tmp = p1;
				p1 = p2;
				p2 = tmp;
			}
		}
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof TrianglePoints)) return false;
		return (((TrianglePoints)o).p0.equals(p0) && 
				((TrianglePoints)o).p1.equals(p1) && 
				((TrianglePoints)o).p2.equals(p2));
	}
	
	public int hashCode() {
		return p0.getId()+p1.getId()+p2.getId();
	}
}
