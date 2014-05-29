package ProGAL.steiner.bnb;

//import java.awt.Color;

import java.util.List;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.viewer.TextShape;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geomNd.Point;
import ProGAL.geomNd.Vector;
import ProGAL.math.Randomization;
import ProGAL.proteins.PDBWebReader;
import ProGAL.proteins.belta.PDBFile;

public class PointPlacementWS {
	public final Point[] points;

	protected final int D;
	protected final int N;

	private final double[][] EL;
	private final double[][] B;
	private final double[][] C;
	private final int[] eqnstack;
	private final int[] leafQ;
	private final int[] val;
	
	public PointPlacementWS(Point[] sites){
		N = sites.length;
		D = sites[0].getDimensions();
		
		Point min = new Point(D); min.fill(Double.POSITIVE_INFINITY);
		Point max = new Point(D); max.fill(Double.NEGATIVE_INFINITY);
		points = new Point[N+N-2];
		for(int i=0;i<N;i++){
			points[i] = sites[i];
			for(int d=0;d<D;d++){
				if(points[i].get(d)<min.get(d)) min.set(d, points[i].get(d));
				if(points[i].get(d)>max.get(d)) max.set(d, points[i].get(d));
			}
		}
		for(int i=N;i<N+N-2;i++){
			points[i] = new Point(D);
			for(int d=0;d<D;d++) points[i].set( d, Randomization.randBetween(min.get(d), max.get(d)) );
		}

		EL = new double[N-2][3];
		B = new double[N][3];
		C = new double[N][D];
		eqnstack = new int[N];
		leafQ = new int[N];
		val = new int[N];

	}


	public double updateSteinerPoints(Topology top, double tol){
//		System.out.println("updateSteinerPoints:");
//		System.out.println(top);
//		J2DScene scene = null; 
//		ProGAL.geom2d.Point[] scenePoints = new ProGAL.geom2d.Point[top.N*2-2]; 
		
		for(int i=N;i<2*N-2;i++)
			points[i].addThis(Vector.randomVector(D, 0.1));
		
		double q,r;
		int it=0;
		do{ 
//			if(it++==1000){
//				scene = J2DScene.createJ2DSceneInFrame();
//				for(int i=0;i<top.N*2-2;i++){
//					scenePoints[i] = new ProGAL.geom2d.Point();
//					scenePoints[i].set(points[i]);
//					scene.addShape(new ProGAL.geom2d.Circle(scenePoints[i], 0.1));
//					scene.addShape(new ProGAL.geom2d.viewer.TextShape(""+i, scenePoints[i], 0.2));
//				}
//				for(int e=0;e<top.edges.length;e++){
//					ProGAL.geom2d.Point p1 = scenePoints[top.edges[e][0]];
//					ProGAL.geom2d.Point p2 = scenePoints[top.edges[e][1]];
//					scene.addShape(new ProGAL.geom2d.LineSegment(p1,p2));
//				}
//			}
			q = length(top); r = error(top);
			optimize(top, tol*r/N);
			
//			if(scene!=null){
//				for(int i=0;i<top.N*2-2;i++)
//					scenePoints[i].set(points[i]);
//				scene.repaint();
//				try {
//					Thread.sleep(50);
//				} catch (InterruptedException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//			}
//		}while(r>0.002*q);
		}while(r>tol*q);
        
        return q;
	}
	
	

	/** Stores edge lengths of tree T in array EL[1..kl][0..2] and returns total length. */
	double length(Topology top)
	{
		int[][] adj = top.steinerAdjacencies;
		
		//#define dist(a,b) t=0.0;for(m=0; m<D; m++){ r=XX[a][m]-XX[b][m]; t +=r*r; }t=SQROOT(t);
		int m,i2,i,j; 
		int n0,n1,n2,kl; 
		double leng,t,r; 
		leng = 0.0; kl = N-2;
		for(i=0; i<top.sites-2; i++){
			i2 = i + N;
			n0 = adj[i][0]; n1 = adj[i][1]; n2 = adj[i][2]; 
			if(n0<i2){
//				t=0.0;for(m=0; m<D; m++){ r=points[n0].get(m)-points[i2].get(m); t +=r*r; } t=Math.sqrt(t);
				t = points[n0].distance(points[i2]);
				leng += t; EL[i][0] = t; n0 -= N;
				if(n0>=0) for(j=0; j<3; j++) if(adj[n0][j]==i2){ EL[n0][j] = t; break; } 
			}
			if(n1<i2){
//				t=0.0;for(m=0; m<D; m++){ r=points[n1].get(m)-points[i2].get(m); t +=r*r; } t=Math.sqrt(t);
				t = points[n1].distance(points[i2]);
				leng += t; EL[i][1] = t; n1 -= N;
				if(n1>=0) for(j=0; j<3; j++) if(adj[n1][j]==i2){ EL[n1][j] = t; break; }
			}
			if(n2<i2){
//				t=0.0;for(m=0; m<D; m++){ r=points[n2].get(m)-points[i2].get(m); t +=r*r; } t=Math.sqrt(t);
				t = points[n2].distance(points[i2]);
				leng += t; EL[i][2] = t; n2 -= N;
				if(n2>=0) for(j=0; j<3; j++) if(adj[n2][j]==i2){ EL[n2][j] = t; break; } 
			}
		}/* Have now figured out distance EL[i][0..2] from Steiner pt. i to neighbors. */
		return leng; 
	}


	/** Finds better coordinates XX[NUMSITES + 1..NUMSITES +kl][] for the kl Steiner points
	 * of tree T by: doing a relaxation iteration. Assumes edge lengths of old tree
	 * have been pre-stored in array EL[][] */
	public void optimize(Topology top, double tol){
		int[][] adj = top.steinerAdjacencies;
		
//		double[][] B = new double[N][3], C = new double[N][D];
//		int[] eqnstack = new int[N], leafQ = new int[N], val = new int[N];
		int i,m,j,i2; int n0,n1,n2,lqp,eqp,kl; double q0,q1,q2,t;
		lqp = eqp = 0; kl = N-2;

		// First: compute B array, C array, and valences. Set up leafQ. 
//		for(i=kl-1; i>=0; i--){
		for(i=top.sites-3; i>=0; i--){
			n0 = adj[i][0]; n1 = adj[i][1]; n2 = adj[i][2];
			q0 = 1.0/(EL[i][0]+tol); 
			q1 = 1.0/(EL[i][1]+tol); 
			q2 = 1.0/(EL[i][2]+tol); 
			//Have now figured out reciprocal distances q0,q1,q2 from
			//Steiner pt. i to neighbors n0,n1,n2 **/
			t = q0+q1+q2; q0/=t; q1/=t; q2/=t;
			val[i] = 0; B[i][0] = B[i][1] = B[i][2] = 0.0;
			for(m=0; m<D; m++){ C[i][m] = 0.0; }
			//if(n0>N){ val[i]++; B[i][0] = q0; } else for(m=0; m<D; m++) C[i][m] += XX[n0][m]*q0;
			if(n0>=N){ val[i]++; B[i][0] = q0; } else for(m=0; m<D; m++) C[i][m] += points[n0].get(m)*q0;
			if(n1>=N){ val[i]++; B[i][1] = q1; } else for(m=0; m<D; m++) C[i][m] += points[n1].get(m)*q1;
			if(n2>=N){ val[i]++; B[i][2] = q2; } else for(m=0; m<D; m++) C[i][m] += points[n2].get(m)*q2;
			//Now: Steiner point i has Steiner valence val[i];
			//coords obey eqns XX[i+N][] = sum(j)of B[i][j]*XX[nj][] + C[i][] */
			if(val[i] <= 1){ leafQ[lqp] = i; lqp++; }/* put leafs on leafQ */ 

		}
		//Have set up equations - now-to solve them. */
		//Second: eliminate leaves */
		while(lqp>1){
			lqp--; i = leafQ[lqp]; val[i]--; i2 = i+N;
			//Now to eliminate leaf i
			eqnstack[eqp] = i; eqp++; //Push i onto stack
			for(j=0; j<3; j++) if(B[i][j] != 0.0) break;// neighbor is 4+j 
			q0 = B[i][j];
			j = adj[i][j]-N;/* neighbor is j */
			val[j]--; if(val[j]==1){ leafQ[lqp] = j; lqp++; }/* new leaf? */
			for(m=0; m<3; m++) if(adj[j][m]==i2) break;
			q1 = B[j][m]; B[j][m] = 0.0;
			t = 1.0-q1*q0; t = 1.0/t;
			for(m=0; m<3; m++) B[j][m] *= t;
			for(m=0; m<D; m++){ C[j][m] += q1*C[i][m]; C[j][m] *= t; }
		}
		//Third: Solve 1-vertex tree!
		i = leafQ[0]; i2=i+N;
		for(m=0; m<D; m++) points[i2].set(m, C[i][m]); 
		//Fourth: backsolve
		while(eqp>0){
			eqp--; i=eqnstack[eqp]; i2=i+N;
			for(j=0; j<3; j++) if(B[i][j] != 0.0) break; //neighbor is #j
			q0 = B[i][j];
			j = adj[i][j]; //Neighbor is j
			for(m = 0; m < D; m++ ) points[i2].set(m,  C[i][m] + q0*points[j].get(m) );
		}
		
		return;
	}

	/** Returns the error figure of tree T with Steiner coords in XX[][].
	 * Assumes edge lengths have been pre-stored in array EL[][]. */ 
	double error(Topology top){
		int[][] adj = top.steinerAdjacencies;

		int i,m,i2,n0,n1,n2; 
		int kl; 
		double r,s,t,efig,d01,d12,d02; 
		kl=N-2;efig=0.0;
//		for(i=0; i<kl; i++){
		for(i=0; i<top.sites-2; i++){
			i2 = i+N;
			n0 = adj[i][0]; n1 = adj[i][1]; n2 = adj[i][2];
			d12 = d01 = d02 = 0.0;
			for(m=0; m<D; m++){
				t = points[i2].get(m);
				r = points[n0].get(m)-t; s = points[n1].get(m)-t; t = points[n2].get(m)-t; 
				d12 += s*t; d01 += r*s; d02 += r*t;
			}
			// only angles < 120 cause error 
			t = d12 + d12 + EL[i][1]*EL[i][2]; if(t>0.0) efig += t; 
			t = d01 + d01 + EL[i][0]*EL[i][1]; if(t>0.0) efig += t; 
			t = d02 + d02 + EL[i][0]*EL[i][2]; if(t>0.0) efig += t;
		}
		efig = Math.sqrt(efig); 
		return(efig);
	}

	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("PointPlacementWS:\n");
		for(int i=0;i<N;i++){
			sb.append("> Site ");
			sb.append(i);
			sb.append(": ");
			sb.append(points[i].toString(3));
			sb.append('\n');
		}
		for(int i=N;i<N+N-2;i++){
			sb.append("> Steiner point ");
			sb.append(i);
			sb.append(": ");
			sb.append(points[i].toString(3));
			sb.append('\n');
		}
		return sb.toString();
	}

	public static void main(String[] args){
		//		Node n=new Node();
		//		n = new Node(n,2,3);
		//		n = new Node(n,6,3);
		//		n = new Node(n,4,3);
		//		Topology top = new Topology(6);
		//		top.updateFromNode(n);
		//		
		//		Point[] sites = new Point[6];
		//		double theta = Math.PI/2, delta = Math.PI*5./180.;
		//		sites[0] = new ProGAL.geom2d.Point(Math.cos(theta-delta),Math.sin(theta-delta));
		//		sites[1] = new ProGAL.geom2d.Point(Math.cos(theta+delta),Math.sin(theta+delta));
		//		theta+=Math.PI*2./3.;
		//		sites[2] = new ProGAL.geom2d.Point(Math.cos(theta-delta),Math.sin(theta-delta));
		//		sites[3] = new ProGAL.geom2d.Point(Math.cos(theta+delta),Math.sin(theta+delta));
		//		theta+=Math.PI*2./3.;
		//		sites[4] = new ProGAL.geom2d.Point(Math.cos(theta-delta),Math.sin(theta-delta));
		//		sites[5] = new ProGAL.geom2d.Point(Math.cos(theta+delta),Math.sin(theta+delta));

//		Node n = new Node();
//		n = new Node(n, 0, 3);
//		Topology top = new Topology(5);
//		top.updateFromNode(n);
//		Point[] sites = new Point[5];
//		sites[0] = new ProGAL.geom2d.Point(0,0);
//		sites[1] = new ProGAL.geom2d.Point(0,1);
//		sites[2] = new ProGAL.geom2d.Point(6,1);
//		sites[3] = new ProGAL.geom2d.Point(6,1);
//		sites[4] = new ProGAL.geom2d.Point(5,0);
//
//		PointPlacementWS pp = new PointPlacementWS(sites);
//		double len = pp.updateSteinerPoints(top,0.0001);
//		System.out.println(top);
//		System.out.println(pp);
//		System.out.println("Len: "+len);
		
		
		

		PDBFile f = new PDBFile(PDBWebReader.downloadPDBFile("2CRO"));
		
		List<ProGAL.geom3d.Point> atoms = f.getAtomCoords();
		Point[] sites = new Point[atoms.size()];
		for(int i=0;i<sites.length;i++){
			
			sites[i] = new Point( atoms.get(i) );

		}
		System.out.println(atoms.size());
		Topology top = new Topology(sites.length);
//		MinimumSpanningTree mst = new MinimumSpanningTree(sites);
//		int[][] edges = mst.getEdges();
//		for(int e1 = 0;e1<edges.length-1;e1++){
//			for(int e2 = e1+1;e2<edges.length;e2++){
//				if(edges[e2][0]==edges[e1][0])
//			}
//		}
		
		Node n = new Node();
		for(int i=3;i<sites.length;i++)
			n = new Node(n, Randomization.randBetween(0, i), i);
		top.updateFromNode(n);
		
		PointPlacementWS pp = new PointPlacementWS(sites);
		long start = System.nanoTime();
		pp.updateSteinerPoints(top, 0.001);
		long end = System.nanoTime();
		System.out.printf("Took %.4f ms\n", (end-start)/1000000.0);
		
		ProGAL.geom3d.viewer.J3DScene scene = ProGAL.geom3d.viewer.J3DScene.createJ3DSceneInFrame();
		int i=0;
		for(Point p: pp.points){
			System.out.println(i++);
			scene.addShape(new ProGAL.geom3d.volumes.Sphere( new ProGAL.geom3d.Point(p), 0.1) , java.awt.Color.GRAY, 4);
		}
		for(int e=0;e<top.edges.length;e++){
			int s1 = top.edges[e][0];
			int s2 = top.edges[e][1];
			ProGAL.geom3d.Point p1 = new ProGAL.geom3d.Point(pp.points[s1]);
			ProGAL.geom3d.Point p2 = new ProGAL.geom3d.Point(pp.points[s2]);
			scene.addShape(new ProGAL.geom3d.volumes.LSS(p1,p2,0.05),java.awt.Color.GRAY,  3);
		}
		scene.centerCamera();
	}
}
