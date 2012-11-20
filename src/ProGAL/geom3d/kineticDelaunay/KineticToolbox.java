package ProGAL.geom3d.kineticDelaunay;

import java.util.Arrays;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;
import static java.lang.Math.atan;

import ProGAL.geom3d.Line;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

public class KineticToolbox {


	/** 
	 * Given a configuration of 5 vertices spanning two neighboring tetrahedra and an indication 
	 * of which vertices are moving (<code>moving</code> specifies indices in <code>vertices</code>),
	 * find the smallest rotation angle that causes the five vertices to be cospherical. The rotation 
	 * is assumed to be a right-hand-rotation around the z-axis.    
	 * @param vertices An array of 5 vertices spanning two neighboring tetrahedra
	 * @param moving An array indicating the indices of the points that are rotating. 
	 * @return 
	 */
	static double nextEvent(Vertex[] vertices, int[] moving){
		return 0;
	}
	
	
	/** 
	 * Given a configuration of 5 vertices and the index of one rotating vertex this method returns the 
	 * smallest absolute rotation angle that causes the five vertices to be cospherical. The rotation 
	 * is assumed to be a right-hand-rotation around the z-axis and the rotating point has already 
	 * rotated an angle of <code>curAlpha</code> even though this is not reflected in its coordinates.
	 *    
	 * @param vertices An array of 5 vertices
	 * @param moving The index of the point that is rotating.
	 * @param curAlpha Rotating points will not have their coordinates changed between events, so 
	 * <code>curAlpha</code> indicates how far they have rotated.   
	 * @return The absolute angle rotation angle of <code>vertices[moving]</code> around the z-axis 
	 * where the five vertices will be co-spherical. 
	 */
	static double nextEvent(Vertex[] vertices, int moving, double curAlpha){
		//Create the inSphere matrix (except the last row)
		Matrix m = new Matrix(5,5);
		for(int r=0;r<4;r++){
			int vIdx = r+(r<moving?0:1);
			for(int c=0;c<3;c++){
				m.set(r, c, vertices[vIdx].get(c));
			}
			m.set(r, 3, vertices[vIdx].dot(vertices[vIdx]));
			m.set(r, 4, 1);
		}

		Vertex v = vertices[moving];
		double r = Math.sqrt(v.x()*v.x()+v.y()*v.y());
		
		//Entering [rcos(alpha), rsin(alpha), v.z, v.dot(v), 1] as last row, expanding the determinant by 
		//minors on the last row and rewriting to the expression K1*cos(alpha) + K2*sin(alpha) + K3 = 0 gives:  
		double K1 = r*m.minor(4, 0).determinant();
		double K2 = -r*m.minor(4, 1).determinant();
		double K3 = v.z()*m.minor(4, 2).determinant()-v.dot(v)*m.minor(4, 3).determinant()+m.minor(4, 4).determinant();
		
		//The solutions to the above equation are of the form alpha = 2*atan( (-K2±Ã(K1^2+K2^2-K3^2))/(K3-K1) )
		double D = K1*K1 + K2*K2 - K3*K3;
		if(D<0) return Double.POSITIVE_INFINITY;

		double t1 = (-K2+Math.sqrt(D))/(K3-K1);
		double t2 = (-K2-Math.sqrt(D))/(K3-K1);
		if(t1>t2){ double tmp = t1; t1 = t2; t2=tmp; }
		
		double alpha1 = 2*Math.atan( t1 );
		double alpha2 = 2*Math.atan( t2 );
		
		//Find the current position of the rotating point and determine the correct earliest absolute rotation angle
		double vAlpha = Math.atan2(v.y(),v.x())+curAlpha;//TODO: Should be atan
		if(vAlpha>alpha2) return alpha1+2*Math.PI;
		if(vAlpha>alpha1) return alpha2;
		return alpha1;
	}
	

	static void testEvent1(){
		Vertex[] vs = {
				new Vertex(new Point(0,1,0)),
				new Vertex(new Point(1,1,0)),
				new Vertex(new Point(0,2,0)),
				new Vertex(new Point(0,1,1)),
				new Vertex(new Point(Math.sqrt(2),0,0))
		};
		double a = nextEvent(vs, 4, 0*Math.PI/180);
		System.out.printf("Vertex 4 will intersect circumsphere of 0-3 at the angle %.2f¡. Expects 45¡\n", a*180/Math.PI);
	}
	
	/** 
	 * Given a configuration of 5 vertices and the index of two rotating vertices this method returns the 
	 * smallest absolute rotation angle that causes the five vertices to be cospherical. The rotation 
	 * is assumed to be a right-hand-rotation around the z-axis and the rotating points have already 
	 * rotated an angle of <code>curAlpha</code> even though this is not reflected in its coordinates.
	 *    
	 * @param vertices An array of 5 vertices
	 * @param moving The index of the point that is rotating.
	 * @param curAlpha Rotating points will not have their coordinates changed between events, so 
	 * <code>curAlpha</code> indicates how far they have rotated.   
	 * @return The absolute angle rotation angle of <code>vertices[moving]</code> around the z-axis 
	 * where the five vertices will be co-spherical. 
	 */
	static double nextEvent(Vertex[] vertices, int moving1, int moving2, double curAlpha){
		//Create the inSphere matrix (except the last two rows)
		Matrix m = new Matrix(5,5);
		for(int r=0;r<3;r++){
			int vIdx = r+(r<moving1?0:1)+(r<moving2?0:1);
			for(int c=0;c<3;c++){
				m.set(r, c, vertices[vIdx].get(c));
			}
			m.set(r, 3, vertices[vIdx].dot(vertices[vIdx]));
			m.set(r, 4, 1);
		}

		Vertex v1 = vertices[moving1];
		double r1 = Math.sqrt(v1.x()*v1.x()+v1.y()*v1.y());
		Vertex v2 = vertices[moving2];
		double r2 = Math.sqrt(v2.x()*v2.x()+v2.y()*v2.y());
		
		//Entering [rcos(alpha), rsin(alpha), v.z, v.dot(v), 1] as last row, expanding the determinant by 
		//minors on the last row and rewriting to the expression K1*cos(alpha) + K2*sin(alpha) + K3 = 0 gives:
		double A12 = m.minor(4, 1).minor(3, 0).determinant();
		double A13 = m.minor(4, 2).minor(3, 0).determinant();
		double A14 = m.minor(4, 3).minor(3, 0).determinant();
		double A15 = m.minor(4, 4).minor(3, 0).determinant();
		double A23 = m.minor(4, 2).minor(3, 1).determinant();
		double A24 = m.minor(4, 3).minor(3, 1).determinant();
		double A25 = m.minor(4, 4).minor(3, 1).determinant();
		double A34 = m.minor(4, 3).minor(3, 2).determinant();
		double A35 = m.minor(4, 4).minor(3, 2).determinant();
		double A45 = m.minor(4, 4).minor(3, 3).determinant();
		double k1 = -r1*r2*A12;
		double k2 =  r2*(v1.z()*A13-v1.dot(v1)*A14+A15);
//		double k3 = -k1;
		double k4 = -r2*(v1.z()*A23-v1.dot(v1)*A24+A25);
		double k5 = -r1*(v2.z()*A13-v2.dot(v2)*A14+A15);
		double k6 =  r1*(v2.z()*A23-v2.dot(v2)*A24+A25);
		double k7 = A34*(v2.dot(v2)*v1.z()-v1.dot(v1)*v2.z()) + A35*(v2.z()-v1.z()) + A45*(v1.dot(v1)-v2.dot(v2));

		double a1 = Math.atan2(v1.y(),v1.x())+curAlpha;
		double a2 = Math.atan2(v2.y(),v2.x())+curAlpha;
		
		double C1 = k2*cos(a2)-k4*sin(a2)+k5*cos(a1)-k6*sin(a1);
		double C2 = k2*sin(a2)+k4*cos(a2)+k5*sin(a1)+k6*cos(a1);
		double C3 = k1*sin(a1-a2)+k7;

		//The solutions to the above equation are of the form alpha = 2*atan( (-K2±Ã(K1^2+K2^2-K3^2))/(K3-K1) )
		double D = C1*C1 + C2*C2 - C3*C3;
		if(D<-Constants.EPSILON) return Double.POSITIVE_INFINITY;
		if(D<Constants.EPSILON){
			double t = -C2/(C3-C1);
			double alpha = -2*atan( t );
			if(curAlpha>alpha) return alpha+2*Math.PI;
			return alpha;
		}

		double t1 = (-C2+Math.sqrt(D))/(C3-C1);
		double t2 = (-C2-Math.sqrt(D))/(C3-C1);
		if(t1>t2){ double tmp = t1; t1 = t2; t2=tmp; }
		
		double alpha1 = -2*atan( t1 );
		double alpha2 = -2*atan( t2 );
//		System.out.println(alpha1*180/Math.PI);
//		System.out.println(alpha2*180/Math.PI);
		
		//Find the current position of the rotating point and determine the correct earliest absolute rotation angle
		if(curAlpha>alpha2) return alpha1+2*Math.PI;
		if(curAlpha>alpha1) return alpha2;
		return alpha1;
	}
	
	
	static void testEvent2(){
		Vertex[] vs = {
				new Vertex(new Point( 0,-0.9, 0)),
				new Vertex(new Point( 0, 0.1, 1)),
				new Vertex(new Point( 0, 1.1, 0)),
				new Vertex(new Point( 1,   0, 0)),
				new Vertex(new Point(-1,   0, 0))
		};
		Line l = new Line(new Triangle(vs[0],vs[1],vs[2]).circumcenter(), new Vector(1,0,1).normalizeThis());
		vs[4] = new Vertex(new Sphere(new Triangle(vs[0],vs[1],vs[2]).circumcenter(), new Triangle(vs[0],vs[1],vs[2]).circumradius()).getIntersection(l).getA());
		l = new Line(new Triangle(vs[0],vs[1],vs[2]).circumcenter(), new Vector(1,0,0).normalizeThis());
		vs[3] = new Vertex(new Sphere(new Triangle(vs[0],vs[1],vs[2]).circumcenter(), new Triangle(vs[0],vs[1],vs[2]).circumradius()).getIntersection(l).getB());
		l = new Line(new Point(0,0,0), new Vector(0,0,1));
		double angle = 10*Math.PI/180;
		l.rotateIn(vs[3], -angle);
		l.rotateIn(vs[4], -angle);
		System.out.println(Arrays.toString(vs[3].getCoords()));
		System.out.println(Arrays.toString(vs[4].getCoords()));
		double a = nextEvent(vs, 3,4, 0*Math.PI/180);
		System.out.printf("Common sphere at angle %.10f¡\n", a*180/Math.PI);
//		Line l = new Line(new Point(0,0,0), new Vector(0,0,1));
//		l.rotateIn(vs[3], a);
//		l.rotateIn(vs[4], a);

		Sphere s0 = Sphere.getMinSphere(vs[0],vs[1],vs[2],vs[3]);
		Sphere s1 = Sphere.getMinSphere(vs[1],vs[2],vs[3],vs[4]);
		s0.toConsole();
		s1.toConsole();
	}
	
	
	static double nextEvent(Vertex[] vertices, int moving1, int moving2, int moving3, double curAlpha){
		//Create the inSphere matrix (except the last two rows)
		Matrix m = new Matrix(5,5);
		for(int r=0;r<3;r++){
			int vIdx = r+(r<moving1?0:1)+(r<moving2?0:1);
			for(int c=0;c<3;c++){
				m.set(r, c, vertices[vIdx].get(c));
			}
			m.set(r, 3, vertices[vIdx].dot(vertices[vIdx]));
			m.set(r, 4, 1);
		}

		Vertex v1 = vertices[moving1];
		double r1 = Math.sqrt(v1.x()*v1.x()+v1.y()*v1.y());
		Vertex v2 = vertices[moving2];
		double r2 = Math.sqrt(v2.x()*v2.x()+v2.y()*v2.y());

		//Entering [rcos(alpha), rsin(alpha), v.z, v.dot(v), 1] as last row, expanding the determinant by 
		//minors on the last row and rewriting to the expression K1*cos(alpha) + K2*sin(alpha) + K3 = 0 gives:
		double A12 = m.minor(4, 1).minor(3, 0).determinant();
		double A13 = m.minor(4, 2).minor(3, 0).determinant();
		double A14 = m.minor(4, 3).minor(3, 0).determinant();
		double A15 = m.minor(4, 4).minor(3, 0).determinant();
		double A23 = m.minor(4, 2).minor(3, 1).determinant();
		double A24 = m.minor(4, 3).minor(3, 1).determinant();
		double A25 = m.minor(4, 4).minor(3, 1).determinant();
		double A34 = m.minor(4, 3).minor(3, 2).determinant();
		double A35 = m.minor(4, 4).minor(3, 2).determinant();
		double A45 = m.minor(4, 4).minor(3, 3).determinant();
		double k1 = -r1*r2*A12;
		double k2 =  r2*(v1.z()*A13-v1.dot(v1)*A14+A15);
		//				double k3 = -k1;
		double k4 = -r2*(v1.z()*A23-v1.dot(v1)*A24+A25);
		double k5 = -r1*(v2.z()*A13-v2.dot(v2)*A14+A15);
		double k6 =  r1*(v2.z()*A23-v2.dot(v2)*A24+A25);
		double k7 = A34*(v2.dot(v2)*v1.z()-v1.dot(v1)*v2.z()) + A35*(v2.z()-v1.z()) + A45*(v1.dot(v1)-v2.dot(v2));

		double a1 = Math.atan2(v1.y(),v1.x())+curAlpha;
		double a2 = Math.atan2(v2.y(),v2.x())+curAlpha;

		double C1 = k2*cos(a2)-k4*sin(a2)+k5*cos(a1)-k6*sin(a1);
		double C2 = k2*sin(a2)+k4*cos(a2)+k5*sin(a1)+k6*cos(a1);
		double C3 = k1*sin(a1-a2)+k7;

		//The solutions to the above equation are of the form alpha = 2*atan( (-K2±Ã(K1^2+K2^2-K3^2))/(K3-K1) )
		double D = C1*C1 + C2*C2 - C3*C3;
		if(D<-Constants.EPSILON) return Double.POSITIVE_INFINITY;
		if(D<Constants.EPSILON){
			double t = -C2/(C3-C1);
			double alpha = 2*atan( t );
			if(curAlpha>alpha) return alpha+2*Math.PI;
			return alpha;
		}

		double t1 = (-C2+Math.sqrt(D))/(C3-C1);
		double t2 = (-C2-Math.sqrt(D))/(C3-C1);
		if(t1>t2){ double tmp = t1; t1 = t2; t2=tmp; }

		double alpha1 = 2*atan( t1 );
		double alpha2 = 2*atan( t2 );
		//				System.out.println(alpha1*180/Math.PI);
		//				System.out.println(alpha2*180/Math.PI);

		//Find the current position of the rotating point and determine the correct earliest absolute rotation angle
		if(curAlpha>alpha2) return alpha1+2*Math.PI;
		if(curAlpha>alpha1) return alpha2;
		return alpha1;
	}
	
	

	static void flip(Tet t, int face){
		Tet n = t.neighbors[face];
		for(int i=0;i<4;i++){
			if(i==face) continue;
			for(int j=0;j<4;j++) 
				if(n.neighbors[j]==t.neighbors[i]){
					flip32(t, n, t.neighbors[i]);
					return;
				}
		}
		flip23(t,n);
	}
	
	/**
	 * Performs a 3-2-flip of the three specified tetrahedra. Two of these tetrahedra are preserved 
	 * and the third (deleted) is returned. 
	 * @requires 	Arrays.equals(t0.corners, Arrays.sort(t0.corners) ) &&
	 * 				Arrays.equals(t1.corners, Arrays.sort(t1.corners) ) &&
	 * 				Arrays.equals(t1.corners, Arrays.sort(t1.corners) ) &&
	 * 				_two vertices shared by all three tetrahedra_ 
	 * @return The deleted tetrahedron.
	 */
	static Tet flip32(Tet t0, Tet t1, Tet t2){
		
		//Locate three non-shared vertices
		Vertex[] vs = new Vertex[3];
		vs[0] = t0.corners[t0.apex(t1)];
		vs[1] = t1.corners[t1.apex(t2)];
		vs[2] = t2.corners[t2.apex(t0)];

		//Locate two shared vertices
		Vertex v0 = t0.corners[0];
		for(int i=0;i<4&&(vs[0]==v0 || vs[1]==v0);i++)
			v0 = t0.corners[i];
		Vertex v1 = t0.corners[0];
		for(int i=0;i<4&&(vs[0]==v1 || vs[1]==v1 || v0==v1);i++)
			v1 = t0.corners[i];
		
		int a0 = t0.indexOf(v0); int b0 = t0.indexOf(v1);
		int a1 = t1.indexOf(v1); int b1 = t1.indexOf(v0);
		
		//Name neighbors
		Tet[][] ns = {
				{t0.neighbors[     a0       ], t0.neighbors[     b0       ]},
				{t1.neighbors[     b1       ], t1.neighbors[     a1       ]},
				{t2.neighbors[t2.indexOf(v0)], t2.neighbors[t2.indexOf(v1)]}
		};
		
		//Change corners of t0 and t1
		t0.corners[b0] = vs[2];
		t1.corners[b1] = vs[0];
		
		//Change neighbors
		t0.neighbors[a0] = t1; 							t1.neighbors[a1] = t0;
		
		t0.neighbors[t0.indexOf_slow(vs[0])] = ns[1][1]; 	if(ns[1][1]!=null) ns[1][1].neighbors[ns[1][1].apex(t1)] = t0;
		t0.neighbors[t0.indexOf_slow(vs[1])] = ns[2][1]; 	if(ns[2][1]!=null) ns[2][1].neighbors[ns[2][1].apex(t2)] = t0;
		///t0.neighbors[t0.indexOf(vs[2])] = ns[0][1]; 	if(ns[0][1]!=null) ns[0][1].neighbors[ns[0][1].apex(t0)] = t0;//Already there
		
		//t1.neighbors[t1.indexOf(vs[0])] = ns[1][0];		if(ns[1][0]!=null) ns[1][0].neighbors[ns[1][0].apex(t1)] = t1;//Already there
		t1.neighbors[t1.indexOf_slow(vs[1])] = ns[2][0];		if(ns[2][0]!=null) ns[2][0].neighbors[ns[2][0].apex(t2)] = t1;
		t1.neighbors[t1.indexOf_slow(vs[2])] = ns[0][0];		if(ns[0][0]!=null) ns[0][0].neighbors[ns[0][0].apex(t0)] = t1;

		//Resort
		t0.sortCorners();
		t1.sortCorners();
		
		return t2;
	}

	/**
	 * Performs a 2-3-flip of the two specified tetrahedra. The two existing tetrahedra are preserved and 
	 * a new is created and returned. 
	 * @requires 	Arrays.equals(t0.corners, Arrays.sort(t0.corners) ) &&
	 * 				Arrays.equals(t1.corners, Arrays.sort(t1.corners) ) &&
	 * 				convex(t0,t1)
	 * @return The newly created tetrahedron.
	 */
	static Tet flip23(Tet t0, Tet t1){
		int a0 = t0.apex(t1);
		int a1 = t1.apex(t0);
		
		//Make a local copy of the entire configuration
		Vertex[] vs = new Vertex[3];
		Tet[][] ns = new Tet[2][3];
		int[][] faceIds = new int[2][3];//Corner-indices of shared triangle
		for(int i=0;i<3;i++){
			faceIds[0][i] = i+(i<a0?0:1);
			faceIds[1][i] = i+(i<a1?0:1);
			vs[i] = t0.corners[faceIds[0][i]];
			ns[0][i] = t0.neighbors[faceIds[0][i]];
			ns[1][i] = t1.neighbors[faceIds[1][i]];
		}
		Vertex v0 = t0.corners[a0];
		Vertex v1 = t1.corners[a1];
		
		//Change corners
		Tet t2 = new Tet(new Vertex[]{ v0, vs[0], vs[1], v1 });
		t0.corners[faceIds[0][0]] = v1;
		t1.corners[faceIds[1][1]] = v0;

		//Change neighbors
		t0.neighbors[     a0      ] = ns[1][0];		if(ns[1][0]!=null) ns[1][0].neighbors[ns[1][0].apex(t1)] = t0; 
//		t0.neighbors[faceIds[0][0]] = ns[0][0]; 	if(ns[0][0]!=null) ns[0][0].neighbors[ns[0][0].apex(t0)] = t0; //Redundant
		t0.neighbors[faceIds[0][1]] = t1;
		t0.neighbors[faceIds[0][2]] = t2;

		t1.neighbors[      a1     ] = ns[0][1];		if(ns[0][1]!=null) ns[0][1].neighbors[ns[0][1].apex(t0)] = t1;
		t1.neighbors[faceIds[1][0]] = t0; 
//		t1.neighbors[faceIds[1][1]] = ns[1][1];		if(ns[1][1]!=null) ns[1][1].neighbors[ns[1][1].apex(t1)] = t1; //Redundant
		t1.neighbors[faceIds[1][2]] = t2;
		
		t2.neighbors[ 0 ] = ns[1][2];				if(ns[1][2]!=null) ns[1][2].neighbors[ns[1][2].apex(t1)] = t2;
		t2.neighbors[ 1 ] = t0;
		t2.neighbors[ 2 ] = ns[0][2];				if(ns[0][2]!=null) ns[0][2].neighbors[ns[0][2].apex(t0)] = t2;
		t2.neighbors[ 3 ] = t1;
		
		//Resort 
		t0.sortCorners();
		t1.sortCorners();
		t2.sortCorners();
		
		return t2;
	}
	
	public static void main(String[] args){
		testEvent2();
	}
	
	static void testSort(){
		Vertex v0 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v1 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v2 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v3 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Tet t = new Tet(new Vertex[]{v1,v0,v3,v2});

		System.out.println(t);
		t.sortCorners();
		System.out.println(t);
	}
	static void testFlip(){
		Vertex dummy = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v0 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v1 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v2 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v3 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		Vertex v4 = new Vertex(new ProGAL.geom3d.Point(0,0,0));
		
		Tet t0 = new Tet(new Vertex[]{ v0,v1,v3,v4 });
		Tet t1 = new Tet(new Vertex[]{ v1,v2,v3,v4 });
		
		Tet n1 = new Tet(new Vertex[]{dummy,v0,v3,v4}){ public String toString(){ return "N1";}};
		Tet n2 = new Tet(new Vertex[]{dummy,v0,v1,v4}){ public String toString(){ return "N2";}};
		Tet n3 = new Tet(new Vertex[]{dummy,v0,v1,v3}){ public String toString(){ return "N3";}};
		Tet n4 = new Tet(new Vertex[]{dummy,v2,v3,v4}){ public String toString(){ return "N4";}};
		Tet n5 = new Tet(new Vertex[]{dummy,v1,v2,v4}){ public String toString(){ return "N5";}};
		Tet n6 = new Tet(new Vertex[]{dummy,v1,v2,v3}){ public String toString(){ return "N6";}};

		t0.neighbors[0] = t1; 
		t0.neighbors[1] = n1;
		t0.neighbors[2] = n2;
		t0.neighbors[3] = n3;
		t1.neighbors[0] = n4;
		t1.neighbors[1] = t0;
		t1.neighbors[2] = n5;
		t1.neighbors[3] = n6;
		
		n1.neighbors[0] = t0;
		n2.neighbors[0] = t0;
		n3.neighbors[0] = t0;
		n4.neighbors[0] = t1;
		n5.neighbors[0] = t1;
		n6.neighbors[0] = t1;

		System.out.println(t0+Arrays.toString(t0.neighbors));
		System.out.println(t1+Arrays.toString(t1.neighbors));
		
		Tet t2 = flip23(t0, t1);
		
		System.out.println();
		System.out.println(t0+Arrays.toString(t0.neighbors));
		System.out.println(t1+Arrays.toString(t1.neighbors));
		System.out.println(t2+Arrays.toString(t2.neighbors));
//		if(true) return;
//		System.out.println();
		
//		t0 = new Tet(new Vertex[]{ v0,v1,v2,v4 });
//		t1 = new Tet(new Vertex[]{ v0,v1,v2,v3 });
//		t0.neighbors[0] = n3; 
//		t0.neighbors[1] = n2;
//		t0.neighbors[2] = n1;
//		t0.neighbors[3] = t1;
//		t1.neighbors[0] = n6;
//		t1.neighbors[1] = n5;
//		t1.neighbors[2] = n4;
//		t1.neighbors[3] = t0;
//
//		System.out.println(t0+Arrays.toString(t0.neighbors));
//		System.out.println(t1+Arrays.toString(t1.neighbors));
//		
//		t2 = flip23(t0, t1);
		
//		System.out.println();
//		System.out.println(t0+Arrays.toString(t0.neighbors));
//		System.out.println(t1+Arrays.toString(t1.neighbors));
//		System.out.println(t2+Arrays.toString(t2.neighbors));
		
		flip32(t0,t1,t2);
		
		System.out.println();
		System.out.println(t0+Arrays.toString(t0.neighbors));
		System.out.println(t1+Arrays.toString(t1.neighbors));

		System.out.println(n1+".neighbors: "+Arrays.toString(n1.neighbors));
		System.out.println(n2+".neighbors: "+Arrays.toString(n2.neighbors));
		System.out.println(n3+".neighbors: "+Arrays.toString(n3.neighbors));
		System.out.println(n4+".neighbors: "+Arrays.toString(n4.neighbors));
		System.out.println(n5+".neighbors: "+Arrays.toString(n5.neighbors));
		System.out.println(n6+".neighbors: "+Arrays.toString(n6.neighbors));
	}
	
}
