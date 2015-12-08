package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.viewer.J3DScene;



public class Hole {
	public List<Vertex> vertices = new ArrayList<Vertex>();
	public List<Face> faces = new ArrayList<Face>();
	J3DScene scene;

	class Face {
		Vertex[] vertices = new Vertex[3];
		Face[] neighbors = new Face[3];
		Tet tet;
		Tet oppTet;
		Shape shape;
		boolean processed = false;
		
		public Face(Vertex a, Vertex b, Vertex c, Vertex u, Tet tet, Tet oppTet) {
			vertices[0] = a;
			if (Point.orientation(a, b, c, u) > 0) {
				vertices[1] = c;
				vertices[2] = b;
			}
			else {
				vertices[1] = b;
				vertices[2] = c;
			}
			this.tet = tet;
			this.oppTet = oppTet;
		}
		
		public Tet getTet() { return tet; }
		
		public int indexOf(Vertex v) {
			for (int i = 0; i < 3; i++) {
				if (vertices[i] == v) return i;
			}
			return -1;
		}
		
		public Vertex getFreeVertex(Face f) {
			for (int i = 0; i < 3; i++) {
				Vertex v = vertices[i]; 
				boolean found = false;
				int j = 0;
				while (!found && (j < 3)) found = v == f.vertices[j++];
				if (!found) return v;	
			}
			return null;
		}

		public int IndexOfFreeVertex(Face f) {
			for (int i = 0; i < 3; i++) {
				Vertex v = vertices[i]; 
				boolean found = false;
				int j = 0;
				while (!found && (j < 3)) found = v == f.vertices[j++];
				if (!found) return i;	
			}
			return -1;
		}

		public boolean containsVertex(Vertex v) {
			return ((v == vertices[0]) || (v == vertices[1]) || (v == vertices[2]));
		}
		
		public Face commonNeighbor(Face face) {
			for (int i = 0; i < 3; i++) 
				for (int j = 0; j < 3; j++) 
					if (neighbors[i] == face.neighbors[j]) return neighbors[i];
			return null;
		}
		
		/* returns common vertex of three tetrahedra */
		public Vertex getCommonVertex(Face f1, Face f2) {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					if (vertices[i] == f1.vertices[j]) {
						for (int k = 0; k < 3; k++) 
							if (vertices[i] == f2.vertices[k]) return vertices[i];
					}
				}
			return null;
		}
		
		public String toString() {
			return "[" + vertices[0] + " " + vertices[1] + " " + vertices[2] + "]";
		}
	}


	
	public Hole(KineticAlphaComplex  triangulation, Vertex u, J3DScene scene, boolean testing) {
		// processes tetrahedra incident with u
		this.scene = scene;
		Stack<Tet> stack = new Stack<Tet>();
		for (Tet tet : triangulation.getTetrahedra()) {
			tet.onStack = false;
			tet.toConsole();
		}
		Tet tet = triangulation.getTetrahedron(u);
		Tet nTet, oppTet;
		int indx;
		stack.push(tet);
		tet.onStack = true;
		while (!stack.isEmpty()) {
			tet = stack.pop();
			indx = tet.indexOf(u);
			oppTet = tet.neighbors[indx];
			if (oppTet == null) 
				System.out.println("oppTet does not exist");
			Face face = new Face(tet.corners[(indx+1)%4], tet.corners[(indx+2)%4], tet.corners[(indx+3)%4], oppTet.corners[oppTet.apex(tet)], tet, oppTet);
//			Face face = new Face(tet.getCorner((indx+1)%4), tet.getCorner((indx+2)%4), tet.getCorner((indx+3)%4), oppTet.getCorner(oppTet.apex(tet)), tet, oppTet);
			tet.selectedFace = face;
			faces.add(face);
			if (testing) face.shape = tet.toSceneFace(scene, indx, Color.red);
			for (int i = 1; i < 4; i++) {
				nTet = tet.neighbors[(indx+i)%4];
				if ((nTet != null) && !nTet.onStack) {
					stack.push(nTet);
					nTet.onStack = true;
				}
			}
		}
		
		// updates neighbor face information
		for (Face f : faces) {
			for (int i = 0; i < 3; i++) {
				Vertex v = f.vertices[i];
				indx = f.tet.indexOf(v);
				nTet = f.tet.neighbors[indx];
				f.neighbors[i] = nTet.selectedFace;
			}
		}
		
		// removes tetrahedra incident with u
		for (Face f : faces) {
			nTet = f.tet.neighbors[f.tet.indexOf(u)];
			for (int i = 0; i < 4; i++) if (nTet.neighbors[i] == f.tet) nTet.neighbors[i] = null;
//			triangulation.removeTetrahedron(f.tet);
			
			if (testing) 
				for (int i = 0; i < 6; i++) scene.removeShape(f.tet.LSSs[i]);
			
			throw new RuntimeException("Sorry .. we removed the removeTetrahedron method");
		}
	}
}
