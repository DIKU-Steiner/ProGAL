package ProGAL.geomNd.complex;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geomNd.Point;

public class Tessel {
	private Tessel[] neighbours;
	private Vertex[] corners;
	private int dimension;
	private boolean modified;
	
	public Tessel(Vertex[] points) {
		int len = points.length;
		dimension = len-1;
		corners = points;//new Point[len];
		neighbours = new Tessel[len];
		/*for (int i=0;i<len;i++) {
			corners[i] = points[i];
		}*/
	}
	
	public void setModified(boolean modi) {
		this.modified = modi;
	}
	
	public boolean isModified() {
		return modified;
	}
	
	public void setNeighbour(int index, Tessel t) {	
		neighbours[index] = t;
	}
	
	public Vertex getVertex(int i) {						
		return corners[i];
	}
	
	public Point getPoint(int i) {
		return corners[i];
	}
	
	public Vertex[] getPoints() {
		return corners;
	}
	
	public Tessel getNeighbour(int index) {		
		return neighbours[index];
	}
	
	public boolean containsBigPoint() {
		for (int i=0;i<corners.length;i++) {
			if (getVertex(i).isBigpoint()){
				return true;
			}
		}
		return false;
	}
	
	//find id of point
	public int findpoint(Vertex p){
		for(int i = 0; i<=dimension; i++){
			if(getVertex(i)==p) {
				return i;
			}
		}
		return -1;
	}
	
	public Vertex findVertex(Tessel tess) {
		Vertex p;
		for (int i = 0; i <= dimension; i++) {
			p = getVertex(i);
			if (!tess.containsPoint(p)) return p; 
		}
		return null;
	}
	
	public int getID(Vertex v) {
		for (int i=0;i<=dimension;i++) {
			if (v == getVertex(i)) return i;
		}
		return -1;
	}
	
	public List<Vertex> oppositeVertices(int i) {
		List<Vertex> points = new ArrayList<Vertex>();
		for (int j=0;j<=dimension;j++) {
			if (findpoint(getVertex(j))!=i) {
				points.add(getVertex(j));
			}
		}
		return points;
	}
	
	public boolean containsPoint(Vertex p) {
		for (int i = 0; i <= dimension; i++) {
			if (getVertex(i) == p) return true;
		}
		return false;
	}
	
	public int apexid(int index){
		Tessel apex_tess= getNeighbour(index);
		if(apex_tess!= null){
			for(int i=0;i<=dimension;i++){
				if(apex_tess.getNeighbour(i)== this){
					return i;
				}
			}
		}
		//never happens:
		return -1;
	}
	
	public List<Vertex> getCommonVertices(Tessel tess) {
		List<Vertex> points = new ArrayList<Vertex>();
		if (tess!=null) {
			for (int i = 0; i < dimension+1; i++) {
				if (tess.containsPoint(this.getVertex(i))) {
					points.add(this.getVertex(i));
				}
			}
		}
		return points;
	}
	
	public String toString() {
		String vertices = new String("");
		for (int i=0;i<dimension+1;i++) {
			vertices = vertices.concat(this.getVertex(i).toString());
			if (i != dimension) {
				vertices += ", ";
			}
		}
		return String.format("Tessel[%s]",vertices);
	}
}

