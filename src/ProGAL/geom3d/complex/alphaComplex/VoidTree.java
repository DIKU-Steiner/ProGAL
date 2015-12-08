package ProGAL.geom3d.complex.alphaComplex;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Simplex;
import ProGAL.math.Randomization;
import ProGAL.proteins.PDBFile;
import ProGAL.proteins.PDBWebReader;

public class VoidTree {
	private AlphaFiltration alphaFil;
	public Node root = null;
//	J3DScene scene = J3DScene.createJ3DSceneInFrame();
	
	public VoidTree(List<Point> points, double interval) {
		this.alphaFil = new AlphaFiltration(points);
		createTree(interval);
	}
	
	// Collects tetrahedra on each side of marked triangle
	private ArrayList<LinkedList<CTetrahedron>> getVoid(CTriangle m, LinkedList<CTriangle> tris){
		ArrayList<LinkedList<CTetrahedron>> v = new ArrayList<LinkedList<CTetrahedron>>();
		LinkedList<CTetrahedron> vOne = new LinkedList<CTetrahedron>();
		LinkedList<CTetrahedron> vTwo = new LinkedList<CTetrahedron>();
		v.add(vOne);
		v.add(vTwo);
		beginning:
			for (int i=0; i<2; i++){
				LinkedList<CTetrahedron> tempVoid = v.get(i);// Contains the tetrahedra set
				LinkedList<CTriangle> tempTris = new LinkedList<CTriangle>();// Contains triangle of tetrahedra
				CTetrahedron inTet = m.getAdjacentTetrahedron(i);// inTet = first tetrahedron
				if (inTet.containsBigPoint()){
					continue beginning;
				}
				tempVoid.addFirst(inTet);
				// Store all but the marked triangle in tempTris
				for (int j=0; j<4; j++){
					tempTris.addFirst(inTet.getTriangle(j));
				}
				tempTris.remove(m);
				while (tempTris.size()>0){
					CTriangle tri = tempTris.removeFirst();
					if (tris.contains(tri)){// If triangle is part of the complex look at next triangle in tempTris
						continue;
					}
					for (int k=0; k<2; k++){
						CTetrahedron adTetra = tri.getAdjacentTetrahedron(k);
						if (tempVoid.contains(adTetra)){// This tetrahedron has been added to tempVoid
							continue;
						}
						if (adTetra.containsBigPoint()){// There is no bounded void on this side of marked triangle
							tempVoid.clear();
							continue beginning;
						}
						tempVoid.addFirst(adTetra);
						// Store all triangles bounding the tetrahedron in tempTris, except the current triangle
						for (int l=0; l<4; l++){
							tempTris.addFirst(adTetra.getTriangle(l));
						}
						tempTris.remove(tri);
					}
				}
			}
		return v;
	}
	
	// Create the root
	public Node setRoot(){
		if (root == null){
			LinkedList<CTetrahedron> orgTetra = new LinkedList<CTetrahedron>(alphaFil.getTetrahedra());
			root = new Node(0, orgTetra);
		}
		return root;
	}
	
	
	// Returns the list of all leaf nodes
	public LinkedList<Node> getLeaves(Node n){
		LinkedList<Node> nodes = new LinkedList<Node>();
		if (n == null){
			return null;
		}
		if (n.left == null){
			nodes.add(n);
			return nodes;
		}
		else { 
			nodes.addAll(getLeaves(n.right));
			nodes.addAll(getLeaves(n.left));
			return nodes;
		}
	}
	
	// Returns the "rest" leaf
	public Node getRest(){
		Node node = root;
		while (node.left != null){
			node = node.left;
		}
		return node;
	}
	
	// Identifies and returns the leaf in the tree with a tetrahedra set containing the argument
	public Node find(LinkedList<CTetrahedron> list){
		LinkedList<Node> nodes = getLeaves(root);
		while (!nodes.isEmpty()){
			Node node = nodes.removeFirst();
			if (node.getTetra().contains(list.getFirst())){
				return node;
			}
		}
		return null;
	}
	
	public void createTree(double interval){
		root = setRoot();
		LinkedList<CTriangle> tris = new LinkedList<CTriangle>();// Set of triangles in complex
		int[][] table = alphaFil.getBettiNumbers();// Table holding information about the dimension of a simplex and if it is marked
		List<Simplex> simplices = alphaFil.getSimplices();
		for (int i=0; i<table[0].length;i++){
			if(table[5][i]==2){// All triangles (dim=2) met in table is added to tris
				tris.add((CTriangle)simplices.get(i));
			}
			// Vertices are added to the scene
//			if(table[5][i]==0) 
//				scene.addShape(
//						new Sphere( (Point)simplices.get(i),0.3 ), 
//						java.awt.Color.BLUE );
			// When a tetrahedron is met a void is destroyed and a leaf obtains a death time different than default
			if (table[5][i]==3){
				LinkedList<CTetrahedron> tets = new LinkedList<CTetrahedron>();
				CTetrahedron tetra = (CTetrahedron) simplices.get(i);
				double d = alphaFil.getInAlpha(simplices.get(i));
				tets.add(tetra);
				Node m = find(tets);
				m.changeDeath(d);
			}
			// A marked triangle is met:
			if (table[4][i]==1 && table[5][i]==2){
				double alpha = alphaFil.getInAlpha(simplices.get(i));// Alpha-value of triangle
				CTriangle marked = (CTriangle) simplices.get(i);// The marked triangle
				ArrayList<LinkedList<CTetrahedron>> newVoids = getVoid(marked, tris);// Collect the voids on each side of marked triangle
				if (newVoids.get(0).isEmpty() && newVoids.get(1).isEmpty()){
					return; //Error!
				}
				// If only one bounded void is created the rest leaf gains children
				if (newVoids.get(0).isEmpty()){
					Node rest = getRest();
					Node newNode = new Node(alpha, newVoids.get(1));
					rest.getTetra().removeAll(newNode.getTetra());
					Node newRest = new Node(0, rest.getTetra());
					rest.setChild(newRest, 0);
					rest.setChild(newNode, 1);
				} else {
					if (newVoids.get(1).isEmpty()){
						Node rest = getRest();
						Node newNode = new Node(alpha, newVoids.get(0));
						rest.getTetra().removeAll(newNode.getTetra());
						Node newRest = new Node(0, rest.getTetra());
						rest.setChild(newRest, 0);
						rest.setChild(newNode, 1);
					// If two bounded voids are created the "parent" void gains children 
					} else {					
						LinkedList<CTetrahedron> v1 = newVoids.get(0);
						LinkedList<CTetrahedron> v2 = newVoids.get(1);
						Node match = find(v1);
						Node newNode1 = new Node(alpha, v1);
						Node newNode2 = new Node(alpha, v2);
						match.setChild(newNode1, 0);
						match.setChild(newNode2, 1);
					}
				}
			}
		}
	}
	
	public static void main(String[] args){
		Randomization.seed(3);
		//ArrayList<Point> points = new ArrayList<Point>();
		// disjoint spheres:
		//List<Point> points = new DisjointSpheresPointList(2, 40);
		
		// joint spheres:
		//List<Point> points = new OverlappingSpheresPointList(2, 40);
		
		// protein:
//		List<Point> points = new ProteinPointList("3SQF");
	
		// Eiffel:
//		List<Point> points = new EiffelPointList();
		
//		//cube: 
		/*points.add(new Point(0,1,0));
		points.add(new Point(1,1,0));
		points.add(new Point(1,2,0));
		points.add(new Point(0,2,0));
		points.add(new Point(0,1,1));
		points.add(new Point(1,1,1));
		points.add(new Point(1,2,1));
		points.add(new Point(0,2,1));*/
		
//		points.add(new Point());
		//2tetrahedra:
//		points.add(new Point(1,3,0));
//		points.add(new Point(1,6,1));
//		points.add(new Point(1,0,1));
//		points.add(new Point(2,3,2));
//		points.add(new Point(0,3,2));
		
		//4tetrahedra:
		/*points.add(new Point(0,0,3));
		points.add(new Point(2,3,6));
		points.add(new Point(1,6,3));
		points.add(new Point(6,5,3));
		points.add(new Point(5,4,4));
		points.add(new Point(2,3,0));*/
		
//		List<Point> points = ProGAL.geom3d.PointList.generatePointsInCube(300);
//		List<Point> points2 = new ArrayList<Point>();
//		for(Point p: points){
//			if(p.x()>0.9 || p.x()<-0.9) points2.add(p);
//			else if(p.y()>0.9 || p.y()<-0.9) points2.add(p);
//			else if(p.z()>0.9 || p.z()<-0.9) points2.add(p);
//		}
//		List<Point> points = ProGAL.geom3d.PointList.generatePointsOnSphere(9);
		
//		Special point set:
		//List<Point> points = new SphereVoidsPointList(300);

//		System.out.println(points.size());
//		VoidTree vt = new VoidTree(points, 0);
//		new ProGAL.dataStructures.viewer.BinaryTreePainter(vt.root);
		
		PDBFile pdb = new PDBFile(PDBWebReader.downloadPDBFile("2CRO"));
		List<Point> points = pdb.getAtomCoords();
		AlphaFiltration af = new AlphaFiltration(points);
		List<CTriangle> triangles = af.getAlphaShape(2.8);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for(CTriangle tri: triangles)
			scene.addShape(tri, java.awt.Color.BLUE);
		scene.centerCamera();
	}
}
