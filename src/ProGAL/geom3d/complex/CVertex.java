package ProGAL.geom3d.complex;

import java.util.ArrayList; 
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointWeighted;

public class CVertex extends PointWeighted {
	private static final long serialVersionUID = 1L;
	
	private boolean degenerate=false;	
	public  enum DegenerateCase { ONFACE, ONEDGE };	
	private DegenerateCase degCase;
	private CVertex degPointOpposite;
	private CVertex degPointA;
	private CVertex degPointB;
	private final boolean bigPoint;
	
	private final List<CEdge> adjacentEdges = new ArrayList<CEdge>();

	
	public boolean isBigpoint() {
		return bigPoint;
	}
	
	public CVertex getDegPointA() {
		return degPointA;
	}

	public void setDegPointA(CVertex degPointA) {
		this.degPointA = degPointA;
	}

	public CVertex getDegPointB() {
		return degPointB;
	}

	public void setDegPointB(CVertex degPointB) {
		this.degPointB = degPointB;
	}	

	public DegenerateCase getDegCase() {
		return this.degCase;
	}

	public void setDegCase(DegenerateCase degCase) {
		this.degCase = degCase;
	}

	public boolean isDegenerate() {
		return degenerate;
	}

	public void setDegenerate(boolean degenerate) {
		this.degenerate = degenerate;
	}

	
	public CVertex(Point p){
		this(p, false);
	}
	public CVertex(Point p, boolean bigpoint){ 	
		super(p);
		setDegenerate(false);
		this.bigPoint = bigpoint;
	}

	public CVertex getDegPointOpposite() {
		return degPointOpposite;
	}

	public void setDegPointOpposite(CVertex degPointOpposite) {
		this.degPointOpposite = degPointOpposite;
	}
	
	public void addAdjacentEdge(CEdge e){
		adjacentEdges.add(e);
	}
	public List<CEdge> getAdjacentEdges(){
		return adjacentEdges;
	}

	public String toString(){ return toString(2); }
	
	public String toString(int dec) {
		return String.format("CVertex[x=%."+dec+"f,y=%."+dec+"f,z=%."+dec+"f%s]",x(),y(),z(),bigPoint?",big point":"");
	}
}
