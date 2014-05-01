package ProGAL.geom3d.kineticDelaunay;

public class TetPoints {
	Vertex p0;
	Vertex p1;
	Vertex p2;
	Vertex p3;
	
	public TetPoints(Vertex p0, Vertex p1, Vertex p2, Vertex p3) {
		this.p0 = p0;
		this.p1 = p1;
		this.p2 = p2;
		this.p3 = p3;
		sortCorners();
	}
	
	protected final void sortCorners(){
		Vertex[] vs = new Vertex[]{p0, p1, p2, p3};
		Vertex tmp;
		if(p0.compareTo(p1)>0) {
//			swap(0,1);//First two elements are sorted
			tmp = vs[0];
			vs[0] = vs[1];
			vs[1] = tmp;
		}
		if(p2.compareTo(p3)>0) {
//			swap(2,3);//Last two elements are sorted
			tmp = vs[2];
			vs[2] = vs[3];
			vs[3] = tmp;
		}
		if(p0.compareTo(p2)>0) {
//			swap(0,2);//Smallest element is at 0
			tmp = vs[0];
			vs[0] = vs[2];
			vs[2] = tmp;
		}
		if(p1.compareTo(p3)>0) {
//			swap(1,3);//Largest element is at 3
			tmp = vs[1];
			vs[1] = vs[3];
			vs[3] = tmp;
		}
		if(p1.compareTo(p2)>0) {
//			swap(1,2);//Check middle elements
			tmp = vs[1];
			vs[1] = vs[2];
			vs[2] = tmp;
		}
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof TetPoints)) return false;
		return (((TetPoints)o).p0.equals(p0) && ((TetPoints)o).p1.equals(p1) && ((TetPoints)o).p2.equals(p2) && ((TetPoints)o).p3.equals(p3));
	}
}
