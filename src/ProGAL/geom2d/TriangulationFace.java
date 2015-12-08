package ProGAL.geom2d;

import java.awt.Color;
import java.util.List;

import ProGAL.dataStructures.DisjointSet.DSNode;
import ProGAL.geom2d.viewer.J2DScene;

public class TriangulationFace {
	
	protected int id;
	protected TriangulationVertex[] corners = new TriangulationVertex[3];
	protected TriangulationFace[] neighbors = new TriangulationFace[3];
	protected boolean shorty = false;
	protected Shape[] edgeShape = new Shape[3];
	protected Shape[] triangleShape = new Shape[3];
	protected Circle circumCircle = null;
	protected Double circumRadius = null;
	public boolean displayed = true;
//	protected Shape circumCircleShape;
	protected Line oppLine = null; 
	protected boolean alive = true;
	protected Integer count = null;
	protected Integer rotCount = null;
	protected ProGAL.geom3d.Triangle liftedTriangle = null;
	protected ProGAL.geom3d.Triangle groundTriangle = null;
	protected ProGAL.geom3d.Plane plane = null;
	protected int delCount = 0;
	protected int uIndx;
	
	public TriangulationFace(TriangulationVertex v0, TriangulationVertex v1, TriangulationVertex v2) {
		corners[0] = v0;
		corners[1] = v1;
		corners[2] = v2;
	}
	public TriangulationFace(TriangulationVertex v0, TriangulationVertex v1, TriangulationVertex v2, 
							 J2DScene scene, boolean testing) {
		corners[0] = v0;
		corners[1] = v1;
		corners[2] = v2;
		if (scene != null) draw(scene, Color.black);
	}
	
	public void killShapes(J2DScene scene) {
		for (int i = 0; i < 3; i++) {
			if (edgeShape[i] != null) scene.removeShape(edgeShape[i]);
			scene.removeShape(circumCircle);
		}
	}
	
	public TriangulationVertex getCorner(int i) { return corners[i]; }
	public TriangulationFace getNeighbor(int i) { return neighbors[i]; }	
	public boolean isShort() { return shorty; }
	public void setShort(boolean shorty) { this.shorty = shorty; }
	public Shape getEdgeShape(int i) { return edgeShape[i]; }
	public void setEdgeShape(int i, Shape shape) { edgeShape[i] = shape; }
	
	/* computes the circumscribing circle */
	public void setCircumCircle() { 
		circumCircle = new Circle(corners[0], corners[1], corners[2]);
		circumRadius = circumCircle.radius;
	}

	/* sets the circumscribing circle to c - not safe */
	public void setCircumCircle(Circle c) { 
		if (c == null) setCircumCircle();
		else {
			circumCircle = c; 
			circumRadius = circumCircle.radius;
		}
	}
	
	public void setCircumCircleCenterX(double x) { circumCircle.center.setCoord(0, x); }
	public void setCircumCircleCenterY(double y) { circumCircle.center.setCoord(1, y); }
	
	public void setCircumCircleRadius(Double r) { circumCircle.radius = r; }
	
	public void setCircumCircleRadius() { setCircumCircle(); }
//	public Shape getCircumCircleShape() { return circumCircleShape; } 
//	public void setCircumCircleShape(Circle c) { circumCircleShape = c; }
	public boolean isAlive() { return alive; }
	public boolean hasShape() { return circumCircle != null; }

	/** returns the bit info about rotating vertices of this face - as a number between 0 and 7 */
	public int getCount() { 
		if (count != null) return count; 
		else {
			count = 0;
			for (int i = 0; i < 3; i++) 
				if (corners[i].getType() == TriangulationVertex.VertexType.R) count = count + (int)Math.pow(2,i);
			return count;
		}
	}

	public int getId() { return id; }
	public void setId(int id) { this.id = id; }
	public void setNeighbor(int indx, TriangulationFace face) { neighbors[indx] = face; }
	public void setAlive(boolean alive) { 
		this.alive = alive; 
	}
	public boolean isBigFace() {
		return corners[0].isBigPoint() || corners[1].isBigPoint() || corners[2].isBigPoint();
	}
	public boolean isFlat() { return Point.collinear(corners[0], corners[1], corners[2]); }
	
	// return TRUE if the rotating vertex v has the smallest index among the rotating vertices of this face 
	public boolean hasLowestRotIndex(TriangulationVertex v, List<Integer> rotIndx, int size) {
		int vIndx = v.getId();
		int indx = size;
		for (int i = 0; i < 3; i++) {
			int cIndx = corners[i].getId();
			if ((cIndx < indx) && (corners[i].getType() == TriangulationVertex.VertexType.R)) indx = cIndx;
		}
		return vIndx == indx;
	}
	
	// returns the number of rotating points of this face
	public int nrRotatingCorners() {
		if (rotCount == null) {
			rotCount = 0;
			for (int i = 0; i < 3; i++)
				if (corners[i].getType() == TriangulationVertex.VertexType.R) rotCount++;
		}
		return rotCount;
	}
	
	public boolean hasRotatingCorners(List<Integer> rotIndx) {
		if (rotIndx.contains(getCorner(0).getId())) return true;
		if (rotIndx.contains(getCorner(1).getId())) return true;
		if (rotIndx.contains(getCorner(2).getId())) return true;
		return false;
	}
	
	public boolean circumCircleContains(List<TriangulationVertex> vertices, double eps) {
		Circle circle = getCircumCircle();
		for (Point p : vertices) if (circle.contains(p, eps)) return true;
		return false;
	}
	
	public Triangle getTriangle() { return new Triangle(corners[0], corners[1], corners[2]); }
	
	/** returns the circumcircle of the face. Requires that the circle is deleted when any vertex of the face is moved */
	public Circle getCircumCircle() {
		if (circumCircle == null) setCircumCircle();
		return circumCircle;
	}
	
	/** returns the radius of the circumcircle of the face. Requires that the radius is deleted when 1 or 2 vertices move */
	public Double getCircumRadius() {
		if (circumRadius == null) getCircumCircle();
		return circumCircle.radius;
	}
	
//	public void destroyCircumCircle(J2DScene scene) { 
//		if (circumCircleShape != null) scene.removeShape(circumCircleShape);
//	}
	
	public Line getOppLine(TriangulationVertex v) {
		if (oppLine == null) {
			int indx = getIndex(v);
			oppLine = new Line(corners[(indx+2)%2], corners[(indx+1)%3]);
		}
		return oppLine;
	}
	
	/** returns the face opposite the vertex v in this face, null if such face does not exist */
	public TriangulationFace getOppFace(TriangulationVertex v) { return getNeighbor(getIndex(v)); }
	
	/** returns TRUE if the face has vertex v as one of its corners */
	public boolean hasVertex(TriangulationVertex v) {
		return ((corners[0] == v) || (corners[1] == v) || (corners[2] == v));
	}
		
	/** returns the index of the specified vertex in the face */
	public int getIndex(TriangulationVertex v) {
		for (int i = 0; i < 3; i++) if (corners[i] == v) return i;
		return -1;
	}
		
	/** returns the index of the specified neighbor face */
	public int getIndex(TriangulationFace t) {
		for (int i = 0; i < 3; i++) if (neighbors[i] == t) return i;
		return -1;
	}
		
	/** returns third vertex of the face */
	public TriangulationVertex getThirdVertex(TriangulationVertex u, TriangulationVertex v) {
		if ((corners[0] != u) && (corners[0] != v)) return corners[0];
		if ((corners[1] != u) && (corners[1] != v)) return corners[1];
		return corners[2];
	}
	
	/** returns vertex of the face not in this.face */
	public TriangulationVertex getThirdVertex(TriangulationFace face) {
		if (!hasVertex(face.getCorner(0))) return face.getCorner(0);
		if (!hasVertex(face.getCorner(1))) return face.getCorner(1);
		return face.getCorner(2);
	}
	
	// removes "red" edges (if any) and draws "black" edges
	public void draw(J2DScene scene, Color clr, double offset) {
		Vector offVector;
		TriangulationVertex a, b, c;
		TriangulationFace oppT;
		int oppIndx;
		b = corners[0];
		c = corners[1];
		if (!isBigFace()) {
			for (int i = 0; i < 3; i++) {
				scene.removeShape(triangleShape[i]);
				scene.removeShape(edgeShape[i]);
				a = b;
				b = c;
				c = corners[(i+2)%3];
				oppT = getOppFace(c);
				oppIndx = oppT.getIndex(b);
				scene.removeShape(oppT.edgeShape[oppIndx]);
				offVector = new Vector(a.y()-b.y(),b.x()-a.x()).scaleToLength(offset);
				triangleShape[i] = new LineSegment(corners[i].add(offVector), corners[(i+1)%3].add(offVector));
				scene.addShape(triangleShape[i], clr, 0.01);
			}
		}
	}
	public void draw(J2DScene scene, Color clr) {
		if (!isBigFace()) {
			TriangulationFace oppT;
			int indxOppT;
			for (int i = 0; i < 3; i++) {
				scene.removeShape(edgeShape[i]);
				edgeShape[i] = new LineSegment(corners[i], corners[(i+1)%3]);
				scene.addShape(edgeShape[i], clr, 0.002);
				if (!isShort()) {
					oppT = getOppFace(getCorner((i+2)%3));
					if ((oppT != null) && oppT.isShort()) {
						indxOppT = oppT.getIndex(getCorner((i+1)%3));
						scene.removeShape(oppT.edgeShape[indxOppT]);
						scene.addShape(oppT.edgeShape[indxOppT]);
					}
				}
				else {
					oppT = getOppFace(getCorner((i+2)%3));
					if (!oppT.isShort()) {
//						indxOppT = oppT.getIndex(getCorner((i+1)%3));
//						scene.removeShape(oppT.edgeShape[indxOppT]);
//						scene.addShape(oppT.edgeShape[indxOppT]);
//						scene.removeShape(edgeShape[i]);
//						scene.addShape(edgeShape[i]);
					}
					
				}
			}
		}
	}
	
	public void draw(J2DScene scene) { draw(scene, Color.black); }	
	
	public void draw(J2DScene scene, boolean testing) {
		if (testing && (scene != null) && !isBigFace()) {
			for (int i = 0; i < 3; i++) {
				scene.removeShape(edgeShape[i]);
				edgeShape[i] = new LineSegment(corners[i], corners[(i+1)%3]);
				scene.addShape(edgeShape[i], Color.black, 0.001);
			}
		}
	}

	public void draw(J2DScene scene, double alpha, boolean testing) {
		shorty = getCircumRadius() <= alpha;
		if (!isShort()) {
			for (int i = 0; i < 3; i++) {
				if (corners[i].distance(corners[(i+1)%3]) < 2*alpha) {
					scene.removeShape(edgeShape[i]);
					edgeShape[i] = new LineSegment(corners[i], corners[(i+1)%3]);
					scene.addShape(edgeShape[i], Color.red, 0.01);
					TriangulationFace oppT = getOppFace(getCorner((i+2)%3));
					if (oppT.isShort()) {
						int indxOppT = oppT.getIndex(getCorner((i+1)%3));
						scene.removeShape(oppT.edgeShape[indxOppT]);
						scene.addShape(oppT.edgeShape[indxOppT]);
					}
				}
			}
		}
		else draw(scene, testing);
	}

	
	public void redraw(J2DScene scene, Color clr) {
		if (!isBigFace()) {
			for (int i = 0; i < 3; i++) {
				scene.removeShape(edgeShape[i]);
				scene.addShape(edgeShape[i], clr, 0.001);
			}
		}
	}
	public void redraw(J2DScene scene) { redraw(scene, Color.black); }	
	
	public void redraw(J2DScene scene, boolean testing) {
		if ((scene != null) && !isBigFace()) {
			Color clr;
			for (int i = 0; i < 3; i++) {
//				clr = Color.black;
				if (corners[i].getType() == corners[(i+1)%3].getType()) clr = Color.black; else clr = Color.black;
				scene.addShape(edgeShape[i], clr, 0.001);
			}
		}
	}

	public Shape drawCircumCircle(J2DScene scene, Color clr) {
		scene.removeShape(circumCircle);
		circumCircle = new Circle(corners[0], corners[1], corners[2]);
		scene.addShape(circumCircle, clr);
		return circumCircle;
	}

	public Shape drawRotatedCircumCircle(J2DScene scene, Color clr, double angle, List<Integer> rotIndx) {
		Point[] p = new Point[3]; 
		for (int i = 0; i < 3; i++) {
			if (rotIndx.contains(corners[i].getId())) p[i] = corners[i].rotationClone(angle); else p[i] = corners[i].clone();
		}
		Shape circle = new Circle(p[0], p[1], p[2]);
		scene.addShape(circle, clr);
		return circle;
	}

	public void hide(J2DScene scene) {
		if (scene != null) 
			for (int i = 0; i < 3; i++) scene.removeShape(edgeShape[i]);
//		scene.repaint();
	}
	
	public void hideCircumCircle(J2DScene scene) {
		scene.removeShape(circumCircle);
		scene.repaint();
	}
	
	public void reshape(J2DScene scene, boolean testing) {
		if (testing && hasShape()) {
			Circle c = new Circle(getCorner(0), getCorner(1), getCorner(2));
			circumCircle.center = c.center;
			circumCircle.radius = c.radius;
		}
	}
	
	public String toString() {
		return "[" + corners[0].id + "," + corners[1].id + "," + corners[2].id + "]"; 
	}
	public TriangulationVertex[] getCorners() {
		return corners;
	}
}
