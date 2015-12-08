package ProGAL.geom3d;

/**
 * A planar rectangle in 3d. The only real functionality so far is the Rectangle-Rectangle distance method that 
 * is used for RSS intersection checks. 
 * @author P.Winter and R.Fonseca
 */
public class Rectangle {
	public Point center;
	public final Vector[] bases;
	private final Vector[] normBases;
	private final double[] extents; 

	public Rectangle(Point center, Vector[] bases){
		this.center = center;
		this.bases = new Vector[]{bases[0], bases[1]};
		this.normBases = new Vector[]{bases[0].normalize(), bases[1].normalize(),bases[0].cross(bases[1]).normalize()};
		this.extents = new double[]{bases[0].length()*2,bases[1].length()*2};
	}


	public double distance(Rectangle rect){
		return distance_optimized(rect);
	}

	/* 196..766 ~*/
	public double distance_optimized(Rectangle rect){
		Point[] thisCorners = getCorners();		//counter-clockwise order starting with bases[0]+bases[1]
		Point[] rectCorners = rect.getCorners();
		Vector[] thisNormals = {bases[1],bases[0].multiply(-1), bases[1].multiply(-1),bases[0]};//6HOps
		Vector[] rectNormals = {rect.bases[1],rect.bases[0].multiply(-1), rect.bases[1].multiply(-1),rect.bases[0]};//6HOps

		boolean[][] inVoronoi1 = geninVoronoi(thisCorners, thisNormals, rectCorners, rectNormals);//48HOps
		boolean[][] inVoronoi2 = geninVoronoi(rectCorners, rectNormals, thisCorners, thisNormals);//48HOps

		int[] perm = {0,2,1,3};
		for(int i: perm){
			for(int j: perm){
				if(!inVoronoi1[i][j] && !inVoronoi1[i][(j+1)%4]) continue;
				if(!inVoronoi2[j][i] && !inVoronoi2[j][(i+1)%4]) continue;
				double c = checkEdgePair(//38HOps
						thisCorners[i],thisCorners[(i+1)%4],thisNormals[i], 
						rectCorners[j],rectCorners[(j+1)%4],rectNormals[j]);
//				System.out.println("c = " + c);
				if(c>=0) {
//					System.out.println("Stopped when i = " + i + " and j = " + j);
					return c;
				}
			}
		}//38 -> 16*38 = 38 -> 608HOps

		double sep1 = axisSeparation(thisCorners,rectCorners);//25HOps
		double sep2 = axisSeparation(rectCorners,thisCorners);//25HOps
		return Math.max(sep1,sep2);
	}

	/** 16*3 = 48HOps */
	private static boolean[][] geninVoronoi(Point[] corners1, Vector[] normals1, Point[] corners2, Vector[] normals2){
		boolean[][] ret = new boolean[4][4];

		//e_0
		ret[0][0] = corners1[0].vectorTo(corners2[0]).dot(normals1[0])>=0;
		ret[0][1] = corners1[0].vectorTo(corners2[1]).dot(normals1[0])>=0;
		ret[0][2] = corners1[0].vectorTo(corners2[2]).dot(normals1[0])>=0;
		ret[0][3] = corners1[0].vectorTo(corners2[3]).dot(normals1[0])>=0;

		//e_1
		ret[1][0] = corners1[2].vectorTo(corners2[0]).dot(normals1[1])>=0;
		ret[1][1] = corners1[2].vectorTo(corners2[1]).dot(normals1[1])>=0;
		ret[1][2] = corners1[2].vectorTo(corners2[2]).dot(normals1[1])>=0;
		ret[1][3] = corners1[2].vectorTo(corners2[3]).dot(normals1[1])>=0;

		//e_2
		ret[2][0] = corners1[2].vectorTo(corners2[0]).dot(normals1[2])>=0;
		ret[2][1] = corners1[2].vectorTo(corners2[1]).dot(normals1[2])>=0;
		ret[2][2] = corners1[2].vectorTo(corners2[2]).dot(normals1[2])>=0;
		ret[2][3] = corners1[2].vectorTo(corners2[3]).dot(normals1[2])>=0;

		//e_3
		ret[3][0] = corners1[0].vectorTo(corners2[0]).dot(normals1[3])>=0;
		ret[3][1] = corners1[0].vectorTo(corners2[1]).dot(normals1[3])>=0;
		ret[3][2] = corners1[0].vectorTo(corners2[2]).dot(normals1[3])>=0;
		ret[3][3] = corners1[0].vectorTo(corners2[3]).dot(normals1[3])>=0;
		return ret;
	}

	/** 670HOps */
	public double distance_nonoptimized(Rectangle rect){
		Point[] thisCorners = getCorners();		//counter-clockwise order starting with bases[0]+bases[1]
		Point[] rectCorners = rect.getCorners();
		Vector[] thisNormals = {bases[1],bases[0].multiply(-1), bases[1].multiply(-1),bases[0]};//6HOps
		Vector[] rectNormals = {rect.bases[1],rect.bases[0].multiply(-1), rect.bases[1].multiply(-1),rect.bases[0]};//6HOps

		int[] perm = {0,2,1,3};
		for(int i: perm){
			for(int j: perm){
				double c = checkEdgePair(//38HOps
						thisCorners[i],thisCorners[(i+1)%4],thisNormals[i], 
						rectCorners[j],rectCorners[(j+1)%4],rectNormals[j]);
				if(c>=0) {
					return c;
				}
			}
		}//16*38 = 608HOps

		double sep1 = axisSeparation(thisCorners,rectCorners);//25HOps
		double sep2 = axisSeparation(rectCorners,thisCorners);//25HOps
		return Math.max(sep1,sep2);
	}

	/** 25HOps */
	private static double axisSeparation(Point[] corners1, Point[] corners2){
		Vector n = corners1[1].vectorTo(corners1[0]).cross(corners1[1].vectorTo(corners1[2])).normalizeThis();//6+7=13HOps
		boolean negatives = false;
		boolean positives = false;
		double min = Double.POSITIVE_INFINITY;
		for(Point c: corners2){
			Vector v = corners1[0].vectorTo(c);
			double dot = v.dot(n);//3HOps
			min = Math.min(min, Math.abs(dot));
			if(dot>0) 	positives = true;
			else		negatives = true;
		}//4*3=12HOps
		if(positives&&negatives) return 0;
		else return min;
	}

	/** 38HOps */
	private double checkEdgePair(Point p1, Point p2, Vector n1, Point q1, Point q2, Vector n2){
		Point[] minDist = closestSegmentPoint(p1,p2,q1,q2);//28HOps
		Vector v = minDist[0].vectorTo(minDist[1]);
		if(v.dot(n1)>0 && v.dot(n2)<0) return v.length();//3+3+4=10HOps
		return -1;
	}

	/** 
	 * 28HOps at most. 
	 */
	public static Point[] closestSegmentPoint(Point p1, Point p2, Point q1, Point q2){
		Point startPoint1 = p1;
		Point startPoint2 = q1;

		Vector dir1 = p1.vectorTo(p2);
		Vector dir2 = q1.vectorTo(q2);

		if(dir1.length()<0.000001 && dir2.length()<0.00001 )
			return new Point[]{startPoint1,startPoint2};
		if(dir1.length()<0.000001) return new Point[]{startPoint1,closestSegmentPoint(startPoint2, q2, startPoint1)};
		if(dir2.length()<0.000001) return new Point[]{closestSegmentPoint(startPoint1, p2, startPoint2),startPoint2};
		//System.out.println("len1 "+d1.length()+" .. len2 "+d2.length());

		Vector r = startPoint2.vectorTo(startPoint1);
		double a = dir1.dot(dir1);//|S1| squared       .. 3HOp
		double e = dir2.dot(dir2);//|S2| squared       .. 3HOp
		double f = dir2.dot(r);//                    .. 3HOp
		double c = dir1.dot(r);//                    .. 3HOp
		double b = dir1.dot(dir2);//                   .. 3HOp
		double denom = a*e-b*b;//                  .. 2HOp

		//If segments not parallel, compute closest point on L1 and L2
		//and clamp to S1
		double s, t;
		if(denom!=0.0f)  s = clamp( (b*f-c*e)/denom );//      .. 3HOp
		else             s = 0.0f;

		//Compute point on L2 closest to S1(S)
		double tnom = b*s+f;//                     .. 1HOp

		//If t in [0,1] done. Else clamp t and recompute and clamp s
		//.. 1 HOp
		if(tnom<0.0f){
			t = 0.0f;
			s = clamp(-c/a);
		}else if(tnom>e){
			t = 1.0f;
			s = clamp( (b-c)/a );
		}else{
			t = tnom/e;
		}

		Point c1 = startPoint1.add(dir1.multiplyThis(s));//      vec-scalar mult  .. 3HOp
		Point c2 = startPoint2.add(dir2.multiplyThis(t));//      vec-scalar mult  .. 3HOp
		return new Point[]{c1,c2};	    
	}

	public static Point closestSegmentPoint(Point p11, Point p12, Point p2){
		Line l = new Line(p11, p11.vectorTo(p12));
		double t = l.orthogonalProjectionParameter(p2);
		t = clamp(t)*p11.distance(p12);
		return l.getPoint(t);
	}

	private static double clamp(double s){
		if(s<0) return 0;
		if(s>1) return 1;
		return s;
	}

	/** Return corners of rectangle in counter-clockwise order. */
	public Point[] getCorners() {
		return new Point[]{
				center.add( bases[0]).addThis( bases[1]),
				center.subtract( bases[0]).addThis(bases[1]),
				center.subtract(bases[0]).subtractThis(bases[1]),
				center.add(bases[0]).subtractThis( bases[1])
		};
	}

	/*
	 * returns the plane through the rectangle
	 */
	public Plane getPlane() {
		return new Plane(center,bases[0].cross(bases[1]).normalizeThis());
	}


	/* ************************************************************************\

	  Copyright 1999 The University of North Carolina at Chapel Hill.
	  All Rights Reserved.

	  Permission to use, copy, modify and distribute this software and its
	  documentation for educational, research and non-profit purposes, without
	  fee, and without a written agreement is hereby granted, provided that the
	  above copyright notice and the following three paragraphs appear in all
	  copies.

	  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
	  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
	  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
	  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
	  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
	  DAMAGES.

	  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
	  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
	  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
	  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
	  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
	  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

	  The authors may be contacted via:

	  US Mail:             E. Larsen
	                       Department of Computer Science
	                       Sitterson Hall, CB #3175
	                       University of N. Carolina
	                       Chapel Hill, NC 27599-3175

	  Phone:               (919)962-1749

	  EMail:               geom@cs.unc.edu


	\**************************************************************************/

	private static double clipToRange(double val, double a, double b){
		if (val < a) return a;
		else if (val > b) return b;
		return val;
	}

	/** Finds the parameters t & u corresponding to the two closest points 
	 * on a pair of line segments 
	 *
	 * The first segment is defined as 
	 *
	 * Pa + A*t, 0 <= t <= a, 
	 * 
	 * where "Pa" is one endpoint of the segment, "A" is a unit vector 
	 * pointing to the other endpoint, and t is a scalar that produces
	 * all the points between the two endpoints. Since "A" is a unit
	 * vector, "a" is the segment's length.
	 *
	 * The second segment is 
	 *
	 * Pb + B*u, 0 <= u <= b
	 *
	 * In my application, many of the terms needed by the algorithm
	 * are already computed for other purposes, so I pass these terms to 
	 * the function instead of complete specifications of each segment. 
	 * "T" in the dot products is the vector between Pa and Pb.
	 *
	 * The algorithm is from
	 *
	 * Vladimir J. Lumelsky,
	 * On fast computation of distance between line segments.
	 * In Information Processing Letters, no. 21, pages 55-61, 1985.
	 */
	private static double[] segCoords(double a, double b, double A_dot_B, double A_dot_T, double B_dot_T){  
		double denom = 1 - (A_dot_B)*(A_dot_B);
		double t,u;

		if (denom == 0) t = 0;
		else{
			t = (A_dot_T - B_dot_T*A_dot_B)/denom; 
			t = clipToRange(t,0,a);
		}

		u = t*A_dot_B - B_dot_T;
		if (u < 0) {
			u = 0;
			t = A_dot_T;
			t = clipToRange(t,0,a);
		}else if (u > b){
			u = b;
			t = u*A_dot_B + A_dot_T;
			t = clipToRange(t,0,a);
		}
		return new double[]{t,u};
	}

	/** 
	 * Returns whether the nearest point on rectangle edge 
	 * Pb + B*u, 0 <= u <= b, to the rectangle edge,
	 * Pa + A*t, 0 <= t <= a, is within the half space 
	 * determined by the point Pa and the direction Anorm.
	 *
	 * A,B, and Anorm are unit vectors.
	 * T is the vector between Pa and Pb.
	 */
	private static boolean inVoronoi( double a, double b, double Anorm_dot_B, double Anorm_dot_T, double A_dot_B, double A_dot_T, double B_dot_T){ 
		if( Math.abs(Anorm_dot_B) < 1e-7 ) return false;

		double t, u, v;

		u = -Anorm_dot_T / Anorm_dot_B; 
		u = clipToRange(u,0,b); 

		t = u*A_dot_B + A_dot_T; 
		t = clipToRange(t,0,a); 

		v = t*A_dot_B - B_dot_T; 

		if (Anorm_dot_B > 0) {
			if (v > (u + 1e-7)) return true;
		} else {
			if (v < (u - 1e-7)) return true;
		}
		return false; 
	} 

	private static double[] MTxV(double[][] M1, double[] V1){
		double[] Vr = new double[3];
		Vr[0] = (M1[0][0] * V1[0] +
				M1[1][0] * V1[1] +
				M1[2][0] * V1[2]);
		Vr[1] = (M1[0][1] * V1[0] +
				M1[1][1] * V1[1] +
				M1[2][1] * V1[2]);
		Vr[2] = (M1[0][2] * V1[0] +
				M1[1][2] * V1[1] +
				M1[2][2] * V1[2]);
		return Vr;
	}

	//	private static double[] MTxV(double[][] M1, double[] V1){
	//		double[] Vr = new double[3];
	//		Vr[0] = (M1[0][0] * V1[0] +
	//				 M1[0][1] * V1[1] +
	//				 M1[0][2] * V1[2]);
	//		Vr[1] = (M1[1][0] * V1[0] +
	//				 M1[1][1] * V1[1] +
	//				 M1[1][2] * V1[2]);
	//		Vr[2] = (M1[2][0] * V1[0] +
	//				 M1[2][1] * V1[1] +
	//				 M1[2][2] * V1[2]);
	//		return Vr;
	//	}
	public double distance_Gottschalk(Rectangle b){
		double[][] Rab = new double[3][3];

		// Compute rotation matrix expressing b in aÕs coordinate frame
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				Rab[i][j] = normBases[i].dot(b.normBases[j]);	
		//Total: 3*3*3 = 27HOps

		// Compute translation vector t
		Vector tmp = center.subtract(bases[0]).subtractThis(bases[1]).vectorTo(
				b.center.subtract(b.bases[0]).subtractThis(b.bases[1])
		);

		// Bring translation into aÕs coordinate frame
		tmp = new Vector(tmp.dot(normBases[0]), tmp.dot(normBases[1]), tmp.dot(normBases[2]));//9HOps
		double[] Tab = new double[]{tmp.x(), tmp.y(), tmp.z()};
		return rectDist(Rab, Tab, extents, b.extents); 
	}

	/** 
	 * Finds the distance between two rectangles A and B.  A is assumed
	 * to have its corner on the origin, one side aligned with
	 * x, the other side aligned with y, and its normal aligned with z.
	 * 
	 * [Rab,Tab] gives the orientation and corner position of rectangle B
	 * 
	 * a[2] are the side lengths of A, b[2] are the side lengths of B
	 */
	private static double rectDist(double[][] Rab, double[] Tab, 
			double[] a, double[] b)	{
		double A0_dot_B0, A0_dot_B1, A1_dot_B0, A1_dot_B1;

		A0_dot_B0 = Rab[0][0];
		A0_dot_B1 = Rab[0][1];
		A1_dot_B0 = Rab[1][0];
		A1_dot_B1 = Rab[1][1];

		double aA0_dot_B0, aA0_dot_B1, aA1_dot_B0, aA1_dot_B1;
		double bA0_dot_B0, bA0_dot_B1, bA1_dot_B0, bA1_dot_B1; 

		aA0_dot_B0 = a[0]*A0_dot_B0;
		aA0_dot_B1 = a[0]*A0_dot_B1;
		aA1_dot_B0 = a[1]*A1_dot_B0;
		aA1_dot_B1 = a[1]*A1_dot_B1;
		bA0_dot_B0 = b[0]*A0_dot_B0;
		bA1_dot_B0 = b[0]*A1_dot_B0;
		bA0_dot_B1 = b[1]*A0_dot_B1;
		bA1_dot_B1 = b[1]*A1_dot_B1;

		double[] Tba = MTxV(Rab,Tab);

		Vector S;
		double t, u;

		// determine if any edge pair contains the closest points

		double ALL_x, ALU_x, AUL_x, AUU_x;
		double BLL_x, BLU_x, BUL_x, BUU_x;
		double LA1_lx, LA1_ux, UA1_lx, UA1_ux, LB1_lx, LB1_ux, UB1_lx, UB1_ux;

		ALL_x = -Tba[0];
		ALU_x = ALL_x + aA1_dot_B0;
		AUL_x = ALL_x + aA0_dot_B0;
		AUU_x = ALU_x + aA0_dot_B0;

		if (ALL_x < ALU_x)
		{ 
			LA1_lx = ALL_x;
			LA1_ux = ALU_x;
			UA1_lx = AUL_x;    
			UA1_ux = AUU_x;
		}
		else
		{ 
			LA1_lx = ALU_x;
			LA1_ux = ALL_x;
			UA1_lx = AUU_x;    
			UA1_ux = AUL_x;
		}

		BLL_x = Tab[0];
		BLU_x = BLL_x + bA0_dot_B1;
		BUL_x = BLL_x + bA0_dot_B0;
		BUU_x = BLU_x + bA0_dot_B0;

		if (BLL_x < BLU_x)
		{ 
			LB1_lx = BLL_x;
			LB1_ux = BLU_x;
			UB1_lx = BUL_x;    
			UB1_ux = BUU_x;
		}
		else
		{ 
			LB1_lx = BLU_x;
			LB1_ux = BLL_x;
			UB1_lx = BUU_x;    
			UB1_ux = BUL_x;
		}

		// UA1, UB1

		if ((UA1_ux > b[0]) && (UB1_ux > a[0]))
		{
			if (((UA1_lx > b[0]) || 
					inVoronoi(b[1],a[1],A1_dot_B0,aA0_dot_B0 - b[0] - Tba[0],
							A1_dot_B1, aA0_dot_B1 - Tba[1], 
							-Tab[1] - bA1_dot_B0))
							&&

							((UB1_lx > a[0]) || 
									inVoronoi(a[1],b[1],A0_dot_B1,Tab[0] + bA0_dot_B0 - a[0],
											A1_dot_B1,Tab[1] + bA1_dot_B0,Tba[1] - aA0_dot_B1)))
			{            
				double[] tu = segCoords(a[1],b[1],A1_dot_B1,Tab[1] + bA1_dot_B0,
						Tba[1] - aA0_dot_B1);
				t = tu[0];
				u = tu[1];

				S = new Vector(
						Tab[0] + Rab[0][0]*b[0] + Rab[0][1]*u - a[0],
						Tab[1] + Rab[1][0]*b[0] + Rab[1][1]*u - t,
						Tab[2] + Rab[2][0]*b[0] + Rab[2][1]*u );
				return S.length();
			}    
		}


		// UA1, LB1

		if ((UA1_lx < 0) && (LB1_ux > a[0]))
		{
			if (((UA1_ux < 0) ||
					inVoronoi(b[1],a[1],-A1_dot_B0,Tba[0] - aA0_dot_B0,
							A1_dot_B1, aA0_dot_B1 - Tba[1], -Tab[1]))
							&&

							((LB1_lx > a[0]) ||
									inVoronoi(a[1],b[1],A0_dot_B1,Tab[0] - a[0],
											A1_dot_B1,Tab[1],Tba[1] - aA0_dot_B1)))
			{
				double[] tu = segCoords(a[1],b[1],A1_dot_B1,Tab[1],Tba[1] - aA0_dot_B1);
				t = tu[0];
				u = tu[1];

				S = new Vector(
						Tab[0] + Rab[0][1]*u - a[0],
						Tab[1] + Rab[1][1]*u - t,
						Tab[2] + Rab[2][1]*u );
				return S.length();
			}
		}

		// LA1, UB1

		if ((LA1_ux > b[0]) && (UB1_lx < 0))
		{
			if (((LA1_lx > b[0]) || 
					inVoronoi(b[1],a[1],A1_dot_B0,-Tba[0] - b[0],
							A1_dot_B1,-Tba[1], -Tab[1] - bA1_dot_B0))
							&&

							((UB1_ux < 0) || 
									inVoronoi(a[1],b[1],-A0_dot_B1, -Tab[0] - bA0_dot_B0,
											A1_dot_B1, Tab[1] + bA1_dot_B0,Tba[1])))
			{

				double[] tu = segCoords(a[1],b[1],A1_dot_B1,Tab[1] + bA1_dot_B0,Tba[1]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*b[0] + Rab[0][1]*u,
						Tab[1] + Rab[1][0]*b[0] + Rab[1][1]*u - t,
						Tab[2] + Rab[2][0]*b[0] + Rab[2][1]*u );
				return S.length();
			}
		}

		// LA1, LB1

		if ((LA1_lx < 0) && (LB1_lx < 0))
		{   
			if (((LA1_ux < 0) || 
					inVoronoi(b[1],a[1],-A1_dot_B0,Tba[0],A1_dot_B1,
							-Tba[1],-Tab[1]))
							&&

							((LB1_ux < 0) || 
									inVoronoi(a[1],b[1],-A0_dot_B1,-Tab[0],A1_dot_B1,
											Tab[1], Tba[1])))
			{
				double[] tu = segCoords(a[1],b[1],A1_dot_B1,Tab[1],Tba[1]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*b[0] + Rab[0][1]*u,
						Tab[1] + Rab[1][0]*b[0] + Rab[1][1]*u - t,
						Tab[2] + Rab[2][0]*b[0] + Rab[2][1]*u );
				return S.length();
			}
		}

		double ALL_y, ALU_y, AUL_y, AUU_y;

		ALL_y = -Tba[1];
		ALU_y = ALL_y + aA1_dot_B1;
		AUL_y = ALL_y + aA0_dot_B1;
		AUU_y = ALU_y + aA0_dot_B1;

		double LA1_ly, LA1_uy, UA1_ly, UA1_uy, LB0_lx, LB0_ux, UB0_lx, UB0_ux;

		if (ALL_y < ALU_y)
		{ 
			LA1_ly = ALL_y;
			LA1_uy = ALU_y;
			UA1_ly = AUL_y;    
			UA1_uy = AUU_y;
		}
		else
		{ 
			LA1_ly = ALU_y;
			LA1_uy = ALL_y;
			UA1_ly = AUU_y;    
			UA1_uy = AUL_y;
		}

		if (BLL_x < BUL_x)
		{
			LB0_lx = BLL_x;
			LB0_ux = BUL_x;
			UB0_lx = BLU_x;
			UB0_ux = BUU_x;
		}
		else
		{
			LB0_lx = BUL_x;
			LB0_ux = BLL_x;
			UB0_lx = BUU_x;
			UB0_ux = BLU_x;
		}

		// UA1, UB0

		if ((UA1_uy > b[1]) && (UB0_ux > a[0]))
		{   
			if (((UA1_ly > b[1]) || 
					inVoronoi(b[0],a[1],A1_dot_B1, aA0_dot_B1 - Tba[1] - b[1],
							A1_dot_B0, aA0_dot_B0 - Tba[0], -Tab[1] - bA1_dot_B1))
							&&

							((UB0_lx > a[0]) || 
									inVoronoi(a[1],b[0],A0_dot_B0, Tab[0] - a[0] + bA0_dot_B1,
											A1_dot_B0, Tab[1] + bA1_dot_B1, Tba[0] - aA0_dot_B0)))
			{
				double[] tu = segCoords(a[1],b[0],A1_dot_B0,Tab[1] + bA1_dot_B1,
						Tba[0] - aA0_dot_B0);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][1]*b[1] + Rab[0][0]*u - a[0] ,
						Tab[1] + Rab[1][1]*b[1] + Rab[1][0]*u - t,
						Tab[2] + Rab[2][1]*b[1] + Rab[2][0]*u );
				return S.length();
			}
		}

		// UA1, LB0

		if ((UA1_ly < 0) && (LB0_ux > a[0]))
		{
			if (((UA1_uy < 0) || 
					inVoronoi(b[0],a[1],-A1_dot_B1, Tba[1] - aA0_dot_B1,A1_dot_B0,
							aA0_dot_B0 - Tba[0], -Tab[1]))
							&&

							((LB0_lx > a[0]) || 
									inVoronoi(a[1],b[0],A0_dot_B0,Tab[0] - a[0],
											A1_dot_B0,Tab[1],Tba[0] - aA0_dot_B0)))
			{
				double[] tu = segCoords(a[1],b[0],A1_dot_B0,Tab[1],Tba[0] - aA0_dot_B0);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*u - a[0],
						Tab[1] + Rab[1][0]*u - t,
						Tab[2] + Rab[2][0]*u );
				return S.length();
			}
		}

		// LA1, UB0

		if ((LA1_uy > b[1]) && (UB0_lx < 0))
		{
			if (((LA1_ly > b[1]) || 
					inVoronoi(b[0],a[1],A1_dot_B1,-Tba[1] - b[1],
							A1_dot_B0, -Tba[0], -Tab[1] - bA1_dot_B1))     
							&&

							((UB0_ux < 0) ||             
									inVoronoi(a[1],b[0],-A0_dot_B0, -Tab[0] - bA0_dot_B1,A1_dot_B0,
											Tab[1] + bA1_dot_B1,Tba[0])))
			{
				double[] tu = segCoords(a[1],b[0],A1_dot_B0,Tab[1] + bA1_dot_B1,Tba[0]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][1]*b[1] + Rab[0][0]*u,
						Tab[1] + Rab[1][1]*b[1] + Rab[1][0]*u - t,
						Tab[2] + Rab[2][1]*b[1] + Rab[2][0]*u );
				return S.length();
			}
		}

		// LA1, LB0

		if ((LA1_ly < 0) && (LB0_lx < 0))
		{
			if (((LA1_uy < 0) || 
					inVoronoi(b[0],a[1],-A1_dot_B1,Tba[1],A1_dot_B0,
							-Tba[0],-Tab[1]))
							&& 

							((LB0_ux < 0) || 
									inVoronoi(a[1],b[0],-A0_dot_B0,-Tab[0],A1_dot_B0,
											Tab[1],Tba[0])))
			{
				double[] tu = segCoords(a[1],b[0],A1_dot_B0,Tab[1],Tba[0]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*u,
						Tab[1] + Rab[1][0]*u - t,
						Tab[2] + Rab[2][0]*u );
				return S.length();
			}
		}

		double BLL_y, BLU_y, BUL_y, BUU_y;

		BLL_y = Tab[1];
		BLU_y = BLL_y + bA1_dot_B1;
		BUL_y = BLL_y + bA1_dot_B0;
		BUU_y = BLU_y + bA1_dot_B0;

		double LA0_lx, LA0_ux, UA0_lx, UA0_ux, LB1_ly, LB1_uy, UB1_ly, UB1_uy;

		if (ALL_x < AUL_x)
		{
			LA0_lx = ALL_x;
			LA0_ux = AUL_x;
			UA0_lx = ALU_x;
			UA0_ux = AUU_x;
		}
		else
		{
			LA0_lx = AUL_x;
			LA0_ux = ALL_x;
			UA0_lx = AUU_x;
			UA0_ux = ALU_x;
		}

		if (BLL_y < BLU_y)
		{
			LB1_ly = BLL_y;
			LB1_uy = BLU_y;
			UB1_ly = BUL_y;
			UB1_uy = BUU_y;
		}
		else
		{
			LB1_ly = BLU_y;
			LB1_uy = BLL_y;
			UB1_ly = BUU_y;
			UB1_uy = BUL_y;
		}

		// UA0, UB1

		if ((UA0_ux > b[0]) && (UB1_uy > a[1]))
		{
			if (((UA0_lx > b[0]) || 
					inVoronoi(b[1],a[0],A0_dot_B0, aA1_dot_B0 - Tba[0] - b[0],
							A0_dot_B1,aA1_dot_B1 - Tba[1], -Tab[0] - bA0_dot_B0))
							&&

							((UB1_ly > a[1]) || 
									inVoronoi(a[0],b[1],A1_dot_B1, Tab[1] - a[1] + bA1_dot_B0,
											A0_dot_B1,Tab[0] + bA0_dot_B0, Tba[1] - aA1_dot_B1)))
			{
				double[] tu = segCoords(a[0],b[1],A0_dot_B1,Tab[0] + bA0_dot_B0,
						Tba[1] - aA1_dot_B1);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*b[0] + Rab[0][1]*u - t,
						Tab[1] + Rab[1][0]*b[0] + Rab[1][1]*u - a[1],
						Tab[2] + Rab[2][0]*b[0] + Rab[2][1]*u );
				return S.length();
			}
		}

		// UA0, LB1

		if ((UA0_lx < 0) && (LB1_uy > a[1]))
		{
			if (((UA0_ux < 0) || 
					inVoronoi(b[1],a[0],-A0_dot_B0, Tba[0] - aA1_dot_B0,A0_dot_B1,
							aA1_dot_B1 - Tba[1],-Tab[0]))
							&&

							((LB1_ly > a[1]) || 
									inVoronoi(a[0],b[1],A1_dot_B1,Tab[1] - a[1],A0_dot_B1,Tab[0],
											Tba[1] - aA1_dot_B1)))
			{
				double[] tu = segCoords(a[0],b[1],A0_dot_B1,Tab[0],Tba[1] - aA1_dot_B1);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][1]*u - t,
						Tab[1] + Rab[1][1]*u - a[1],
						Tab[2] + Rab[2][1]*u );
				return S.length();
			}
		}

		// LA0, UB1

		if ((LA0_ux > b[0]) && (UB1_ly < 0))
		{
			if (((LA0_lx > b[0]) || 
					inVoronoi(b[1],a[0],A0_dot_B0,-b[0] - Tba[0],A0_dot_B1,-Tba[1],
							-bA0_dot_B0 - Tab[0]))
							&&

							((UB1_uy < 0) || 
									inVoronoi(a[0],b[1],-A1_dot_B1, -Tab[1] - bA1_dot_B0,A0_dot_B1,
											Tab[0] + bA0_dot_B0,Tba[1])))
			{
				double[] tu = segCoords(a[0],b[1],A0_dot_B1,Tab[0] + bA0_dot_B0,Tba[1]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*b[0] + Rab[0][1]*u - t,
						Tab[1] + Rab[1][0]*b[0] + Rab[1][1]*u,
						Tab[2] + Rab[2][0]*b[0] + Rab[2][1]*u );
				return S.length();
			}
		}

		// LA0, LB1

		if ((LA0_lx < 0) && (LB1_ly < 0))
		{
			if (((LA0_ux < 0) || 
					inVoronoi(b[1],a[0],-A0_dot_B0,Tba[0],A0_dot_B1,-Tba[1],
							-Tab[0]))
							&&

							((LB1_uy < 0) || 
									inVoronoi(a[0],b[1],-A1_dot_B1,-Tab[1],A0_dot_B1,
											Tab[0],Tba[1])))
			{
				double[] tu = segCoords(a[0],b[1],A0_dot_B1,Tab[0],Tba[1]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][1]*u - t,
						Tab[1] + Rab[1][1]*u,
						Tab[2] + Rab[2][1]*u );
				return S.length();
			}
		}

		double LA0_ly, LA0_uy, UA0_ly, UA0_uy, LB0_ly, LB0_uy, UB0_ly, UB0_uy;

		if (ALL_y < AUL_y)
		{
			LA0_ly = ALL_y;
			LA0_uy = AUL_y;
			UA0_ly = ALU_y;
			UA0_uy = AUU_y;
		}
		else
		{
			LA0_ly = AUL_y;
			LA0_uy = ALL_y;
			UA0_ly = AUU_y;
			UA0_uy = ALU_y;
		}

		if (BLL_y < BUL_y)
		{
			LB0_ly = BLL_y;
			LB0_uy = BUL_y;
			UB0_ly = BLU_y;
			UB0_uy = BUU_y;
		}
		else
		{
			LB0_ly = BUL_y;
			LB0_uy = BLL_y;
			UB0_ly = BUU_y;
			UB0_uy = BLU_y;
		}

		// UA0, UB0

		if ((UA0_uy > b[1]) && (UB0_uy > a[1]))
		{
			if (((UA0_ly > b[1]) || 
					inVoronoi(b[0],a[0],A0_dot_B1, aA1_dot_B1 - Tba[1] - b[1],
							A0_dot_B0, aA1_dot_B0 - Tba[0], -Tab[0] - bA0_dot_B1))
							&&

							((UB0_ly > a[1]) || 
									inVoronoi(a[0],b[0],A1_dot_B0,Tab[1] - a[1] + bA1_dot_B1,A0_dot_B0,
											Tab[0] + bA0_dot_B1, Tba[0] - aA1_dot_B0)))
			{
				double[] tu = segCoords(a[0],b[0],A0_dot_B0,Tab[0] + bA0_dot_B1,
						Tba[0] - aA1_dot_B0);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][1]*b[1] + Rab[0][0]*u - t,
						Tab[1] + Rab[1][1]*b[1] + Rab[1][0]*u - a[1],
						Tab[2] + Rab[2][1]*b[1] + Rab[2][0]*u );
				return S.length();
			}
		}

		// UA0, LB0

		if ((UA0_ly < 0) && (LB0_uy > a[1]))
		{
			if (((UA0_uy < 0) || 
					inVoronoi(b[0],a[0],-A0_dot_B1,Tba[1] - aA1_dot_B1,A0_dot_B0,
							aA1_dot_B0 - Tba[0],-Tab[0]))
							&&      

							((LB0_ly > a[1]) || 
									inVoronoi(a[0],b[0],A1_dot_B0,Tab[1] - a[1],
											A0_dot_B0,Tab[0],Tba[0] - aA1_dot_B0)))
			{
				double[] tu = segCoords(a[0],b[0],A0_dot_B0,Tab[0],Tba[0] - aA1_dot_B0);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*u - t,
						Tab[1] + Rab[1][0]*u - a[1],
						Tab[2] + Rab[2][0]*u );
				return S.length();
			}
		}

		// LA0, UB0

		if ((LA0_uy > b[1]) && (UB0_ly < 0))
		{  
			if (((LA0_ly > b[1]) ||
					inVoronoi(b[0],a[0],A0_dot_B1,-Tba[1] - b[1], A0_dot_B0,-Tba[0],
							-Tab[0] - bA0_dot_B1))
							&&

							((UB0_uy < 0) ||
									inVoronoi(a[0],b[0],-A1_dot_B0, -Tab[1] - bA1_dot_B1, A0_dot_B0,
											Tab[0] + bA0_dot_B1,Tba[0])))
			{
				double[] tu = segCoords(a[0],b[0],A0_dot_B0,Tab[0] + bA0_dot_B1,Tba[0]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][1]*b[1] + Rab[0][0]*u - t,
						Tab[1] + Rab[1][1]*b[1] + Rab[1][0]*u,
						Tab[2] + Rab[2][1]*b[1] + Rab[2][0]*u );
				return S.length();
			}
		}

		// LA0, LB0

		if ((LA0_ly < 0) && (LB0_ly < 0))
		{   
			if (((LA0_uy < 0) || 
					inVoronoi(b[0],a[0],-A0_dot_B1,Tba[1],A0_dot_B0,
							-Tba[0],-Tab[0]))
							&&

							((LB0_uy < 0) || 
									inVoronoi(a[0],b[0],-A1_dot_B0,-Tab[1],A0_dot_B0,
											Tab[0],Tba[0])))
			{
				double[] tu = segCoords(a[0],b[0],A0_dot_B0,Tab[0],Tba[0]);
				t = tu[0];
				u = tu[1];
				S = new Vector(
						Tab[0] + Rab[0][0]*u - t,
						Tab[1] + Rab[1][0]*u,
						Tab[2] + Rab[2][0]*u );
				return S.length();
			}
		}

		// no edges passed, take max separation along face normals

		double sep1, sep2;

		if (Tab[2] > 0.0)
		{
			sep1 = Tab[2];
			if (Rab[2][0] < 0.0) sep1 += b[0]*Rab[2][0];
			if (Rab[2][1] < 0.0) sep1 += b[1]*Rab[2][1];
		}
		else
		{
			sep1 = -Tab[2];
			if (Rab[2][0] > 0.0) sep1 -= b[0]*Rab[2][0];
			if (Rab[2][1] > 0.0) sep1 -= b[1]*Rab[2][1];
		}

		if (Tba[2] < 0)
		{
			sep2 = -Tba[2];
			if (Rab[0][2] < 0.0) sep2 += a[0]*Rab[0][2];
			if (Rab[1][2] < 0.0) sep2 += a[1]*Rab[1][2];
		}
		else
		{
			sep2 = Tba[2];
			if (Rab[0][2] > 0.0) sep2 -= a[0]*Rab[0][2];
			if (Rab[1][2] > 0.0) sep2 -= a[1]*Rab[1][2];
		}

		double sep = (sep1 > sep2? sep1 : sep2);
		return (sep > 0? sep : 0);
	}
}
