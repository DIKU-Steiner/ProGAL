package ProGAL.geom3d.tess;

import ProGAL.math.Randomization;

public class OldTESS {

//	private static final int NVERT = 15; // number of random vertices 
//
//	private static final int MAXVERT   =   25; // n+4 
//	private static final int MAXTET    =  250; // at least 8*n 
//	private static final int MAXCORNER = 1000; // 4*MAXTET 
//
//	private static double RANDCOORD() { return Randomization.randBetween(0.0, 100.0); }//(double)(random(0x100))
//
//	private static int TET(int corner) { return corner/4; }//((corner)>>2)
//	private static int mod4(int i) { return i&3; }//((i)&3)
//	private static int INDEX(int corner) { return corner&3; }//((corner)&3)
//	private static int BASECORNER(int corner) { return corner&0xFFFC; }//((corner)&0xfffc)
//	private static int CORNER(int tet, int index) { return tet*4+index; }//(((tet)<<2)+(index))
//
//	private void subtractV(TPoint v) {
//		if (!infiniteV(v)) {
//			v.xv = v.x - pv.x; v.yv = v.y - pv.y;
//			v.zv = v.z - pv.z; 
//		} 
//	}
//	private static double dot(TPlane pl,TPoint pv) {
//		return ((pl).x*(pv).x+(pl).y*(pv).y+(pl).z*(pv).z+(pl).w);
//	}
//	private static double dotInf(TPlane pl, TPoint pv) {
//		return ((pl).x*(pv).x+(pl).y*(pv).y+(pl).z*(pv).z); // AUDIT only 
//	}
//
//	//Cant figure out these datatypes yet
//	//	private static double DET2(p,q,i,j)((p)->i*(q)->j-(p)->j*(q)->i)
//
//	private void initQueue() { flink[0] = blink[0] = 0; }
//	private boolean emptyQueue() { return (flink[0] == 0); }
//	// enqueue at the back; remove from front
//	private void enqueue(int c) {  
//		//		assertNotInQueue(c);
//		//		ASSERT(blink[0] >= 0, "blink[0] negative"); 
//		flink[c] = 0; flink[blink[c] = blink[0]] = c; blink[0] = c; 
//	}
//	private void deletefromQueue(int c) {
//		//		assertInQueue(c);
//		//		ASSERT((blink[c]>=0 && flink[c] >= 0 && blink[flink[c]] >= 0), "bad link in delete");
//		//		ASSERT(c, "delete sentinel!"); 
//		blink[flink[blink[c]] = flink[c]] = blink[c]; 
//	}
//	private void dequeue(int c) { 
//		c = s[flink[0]].opp; blink[flink[0] = flink[flink[0]]] = 0; 
//		//		ASSERT((flink[0]>=0 && flink[flink[0]] >= 0 && blink[flink[flink[0]]] >= 0), "bad link in dequeue");\
//		/*printf(" DeQ(%dv%d,%dv%d)", s[c].opp, s[s[c].opp].v, c, s[c].v); /**/
//	}
//
//
//
//
//	// ============= From del3.c ==============
//	TPoint pv; // vertex index in outermost loop 
//	double[] od = new double[5]; // orientation determinants  
//
//	int nvert; // number of vertices 
//	TPoint[] vert = new TPoint[MAXVERT];  // vertex coordinates 
//	// Table-based tetrahedra 
//	Corner[] s = new Corner[MAXCORNER]; // per corner entries: v, opp, plane 
//	double[] orientDet = new double[MAXTET]; // per tetrahedron value: orientDet == v*plane > 0 for any tet corner 
//
//	int[] flink = new int[MAXCORNER]; // link for free list, doubly-linked queue  
//	int[] blink = new int[MAXCORNER]; // probably per tet, but per corner for now 
//	int psp;
//	int[] pst = new int[MAXTET]; // stack for tets incident on new point.   
//
//	int[] flag = new int[MAXTET]; // flag for AUDIT 
//
//
//	int freeTet; // head for free list for tetrahedra kept in flink 
//	int maxCorner = 0; // AUDIT only
//
//	private void startTet(int c, double od) {//WARNING: no write-back to c
//		c = freeTet; 
//		freeTet = flink[c]; 
//		orientDet[TET(c)] = od;
//		if (maxCorner < c+4) maxCorner = c+4; // AUDIT
//	}
//	private int freeBase(int b) {
//		//		printf("%%freeBase %d(%d,%d) ", b, TET(b), INDEX(b)); /**/ASSERT(INDEX(b)==0, "nonzero index in freeBase");
//		//		assertNotInQueue(b+3); 
//		blink[b] = -1; 
//		setCornerVC(b,0,MAXCORNER-4); 
//		setCornerVC(b+1,0,MAXCORNER-3); 
//		setCornerVC(b+2,0,MAXCORNER-2); 
//		setCornerVC(b+3,0,MAXCORNER-1); 
//		flink[b] = freeTet; freeTet = b;
//	}
//
//	/* Tetrahedra are groups of five corners in order of increasing 
//	vertex index, except that the first two may be swapped to ensure 
//	that the orientation determinant is positive.  
//	I.e., take the lex smallest alternating permutation with positive sign. 
//	There must always be an odd number of swaps between two permutations; 
//	we swap the first two if necessary to achieve this. 
//
//	The table indoff has the index offset for where 
//	the vertex of index i will be found in the tetrahedron opposite c.  
//	Vertex s[BASECORNER(c)+i].v will be at s[s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]].v 
//	(except that s[c].v and s[s[c].opp]].v are different, and opposite sides of plane s[c].plane).
//	Note that indoff[*][j] uses each offset -j:4-j exactly once, 
//	that indoff[i][j][i] = 0 for all possible i and j (so s[c].v is opposite of s[s[c].opp].v), 
//	and that the tet will never change, TET(s[c].opp) == TET(CORNERINOPP(i,c)).
//	 */
//
//
//	// for reference: this is the table of where the indices go
//	//  \ With V falling at position:
//	//Replacing:\0ABCD 1ABCD 2ABCD 3ABCD  
//	//A    0: bvcd  BVCD  CBVD  BCDV 
//	//B    1: VACD  vacd  ACVD  CADV  
//	//C    2: VBAD  AVBD  BAVD  ABDV  
//	//D    3: VABC  BVAC  ABVC  BACV 
//	//Is this a problem? Lowercase combinations are not allowed 
//	//because they have negative orientation determinants.
//	//2,0 3,1 & 4,0 require a 0,1 swap first, then move V to position.
//	//Since we always add points at the last position, I think this works.
//	//
//	//Table where letter i goes to entry[i]: e.g ABCD->CBVD=2103
//	//Replacing \0ABCD 1ABCD 2ABCD 3ABCD 
//	//A    0: vbcd  1023  2103  3012
//	//B    1: 1023  avcd  0213  1302
//	//C    2: 2103  0213  1023  0132
//	//D    3: 1230  2031  0132  1023
//	//
//	//Tables with offsets I=-1, Z=-2, B = -3
//	//Replacing \0ABCD 1ABCD 2ACBD 3ABCD 
//	//A    0: bvcd  0I12  0IZ1  0BZI  
//	//B    1: 1023  vacd  Z0I1  Z0BI  
//	//C    2: 2103  I102  IZ01  BZ0I  
//	//D    3: 1230  1I20  ZI10  ZBI0
//
//	private int CORNERINOPP(int i,int c) {
//		return (s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]);
//	}
//	private final int[][][] indoff = 
//	{/*        0ABCD        1ABCD         2ACBD          3ABCD    */
//			{/* 0: */ {5,5,5,5},   { 0,-1, 1, 2},  {  0,-1,-2, 1},   { 0,-3,-2,-1}},
//			{/* 1: */ {1,0,2,3},   { 5, 5, 5, 5},  { -2, 0,-1, 1},   {-2, 0,-3,-1}},
//			{/* 2: */ {2,1,0,3},   {-1, 1, 0, 2},  { -1,-2, 0, 1},   {-3,-2, 0,-1}},
//			{/* 3: */ {1,2,3,0},   { 1,-1, 2, 0},  { -2,-1, 1, 0},   {-2,-3,-1, 0}} 
//	};
//
//	// For computing plane equation dropping i, use b+drop[i][0:2].
//	// Order makes pv*plane >= 0, when pv is placed last
//	private final int[][] drop = {{2,1,3}, {0,2,3}, {1,0,3}, {0,1,2}};
//
//	private boolean EQUALV(int i,int j) {
//		return (vert[i].x == vert[j].x && vert[i].y == vert[j].y && vert[i].z == vert[j].z);
//	}
//
//	void setVert(int i, double xx, double yy, double zz, double weight, double pp) {
//		vert[i].x=xx; vert[i].y=yy; vert[i].z=zz;
//		vert[i].p=pp; vert[i].sq = (double) (xx*xx + yy*yy + zz*zz-weight);
//	}
//
//	void setInfVert(int i, double xx, double yy, double zz, double qq) {
//		vert[i].x=vert[i].xv=xx; vert[i].y=vert[i].yv=yy;
//		vert[i].z=vert[i].zv=zz; vert[i].p=i; vert[i].sq = qq; // set sq as well
//	}
//
//	//WARNING: wtf?
//	private final boolean infiniteV(int pv) { return ((pv) < vert+4); } // first four points are at infinity
//
//	private static final boolean STUB(int c) { return ((c) < 4); } // true if corner in stub tetrahedron (never change stub)
//
//	// set corner's vertex, opposite, and plane in tetrahedron structure
//	//void setCornerVC(int c, int vv, int op) {
//	private void setCornerVC(int c, TPoint vv, int op) { s[c].v = vv; s[c].opp = op; }
//	private void setOnePlane(int c, double ww, double xx, double yy, double zz) { 
//		s[c].plane.w = ww; s[c].plane.x = xx; s[c].plane.y = yy; s[c].plane.z = zz; 
//	}
//	private void setPlanes(int c, double ww,double xx,double yy,double zz) { 
//		setOnePlane(c, ww,xx,yy,zz); 
//		setOnePlane(s[c].opp, -ww, -xx, -yy, -zz); 
//	}
//	private void setCornerPairV(int c, int op, TPoint vv, TPoint ov) { 
//		setCornerVC(c, vv, op); 
//		setCornerVC(op, ov, c); 
//	}
//	//void setCornerVCP(int c, int v, int opp, double w,double x,double y,double z) { 
//	private void setCornerVCP(int c, TPoint v, int opp, double w, double x,double y,double z) { 
//		setCornerVC(c, v, opp); 
//		setOnePlane(c, w,x,y,z); 
//	}
//	//void setCornerVCNP(int c, int vv, int op) { // use negative of opp's plane; also set opp.opp=c
//	private void setCornerVCNP(int c, TPoint vv, int op) { 
//		setCornerVC(c, vv, op); s[op].opp = c; 
//		setOnePlane(c, -s[op].plane.w, -s[op].plane.x, -s[op].plane.y, -s[op].plane.z); 
//	}
//	// put them in suspect queue, as well as above
//	private void setCornerVCNPQ(int c, TPoint vv, int op) {
//		setCornerVCNP(c, vv, op); 
//		enqueue(c);
//	}
//
//
//
//	void initTets(){
//		int j;
//		int last = -1; /* set up free list of corners, kept in blocks of five */
//		int stub[] = {1,0,2,3}; // Stub will have vertices reversed to negate orientation determinant.
//
//		for (freeTet = CORNER(MAXTET,0); freeTet > CORNER(2,0); ) {
//			freeTet -= 4;
//			flink[freeTet] = last;
//			last = freeTet;
//		}
//
//		/* Tetrahedra 0 and 1 are special */
//		orientDet[0] = -0.0; //  tet 0 is stub used as header; should never contain points, never change
//		orientDet[1] =  0.0; //  tet 1 is infinite; contains all points
//
//		for (j=0; j<4; j++) {
//			setCornerVCP(CORNER(1,j), vert[j], CORNER(0,stub[j]), 1, 0, 0, 0); // make infinite plane
//			setCornerVCNP(CORNER(0,stub[j]), vert[j], CORNER(1,j));
//		}
//
//		maxCorner = 8;
//		initQueue();
//	}
//
//
//	void initVerts() {
//		int j;
//		double IEEEBITSMAX =  8.0*1024*1024*1024*1024; /* 2^53 */
//		double[] exponentShift = new double[4];
//		/* Here is a trick that uses the IEEE floating point to correctly
//		     compute the sign of a dot product where one of the terms can take
//		     values with different infinities, and the sign is determined by
//		     the largest non-zero infinity, or by the finite value.  If we use
//		     powers of two for the infinities that are powers of 2^53, then we
//		     get the result in different bit positions, and the first non-zero
//		     infinity determines the sign.  (Important that no term is more
//		     than an IEEE double, and that each infinity used is different;
//		     can't have cancellation from the infinities, or we'd lose bits.
//		     We can have up to floor(1024/53)= 19 different infinities... */
//
//		exponentShift[3] = IEEEBITSMAX;
//		for (j = 3; j > 0; j--)
//			exponentShift[j-1] = exponentShift[j]*IEEEBITSMAX;
//		setInfVert(0, 1,0,0, exponentShift[0]); // Five vertices at infinity
//		setInfVert(1, 0,1,0, exponentShift[1]);
//		setInfVert(2, 0,0,1, exponentShift[2]);
//		setInfVert(3, -1,-1,-1, exponentShift[3]);
//		nvert = 4;
//		setVert(nvert++, 0, 0, 0, 0, 6);
//		setVert(nvert++, 4, 9,17, 0, 10);
//		setVert(nvert++,19,34, 7, 0, 17);
//		setVert(nvert++, 1, 5, 8, 0, 7);
//		setVert(nvert++,26,11,20, 0, 19);
//		setVert(nvert++, 3, 2, 1, 0, 9);
//		setVert(nvert++,13, 8,22, 0, 14);
//		setVert(nvert++, 6, 7, 3, 0, 11);
//		setVert(nvert++,10, 1, 6, 0, 13);
//		setVert(nvert++,15,12,16, 0, 15);
//		setVert(nvert++,36,16, 4, 0, 20);
//		setVert(nvert++,17,39,29, 0, 16);
//		setVert(nvert++, 2, 4, 5, 0, 8);
//		setVert(nvert++,22,15, 2, 0, 18);
//		setVert(nvert++,43,35,12, 0, 21);
//		setVert(nvert++, 8,41,43, 0, 12);/**/
//
//	}
//
//	void cornerPrint(int c) {
//		System.out.printf(
//				"%%%3d(%2d,%1d)%2d=(%d %5.0f %5.0f %5.0f) opp:%4d(%3d,%2d) <%5.0f %5.0f %5.0f %5.0f> %7.0f=%7.0f sq:%g\n",
//				c, TET(c), INDEX(c), s[c].v-vert, (infiniteV(s[c].v)?0:1),
//				s[c].v->x, s[c].v->y, s[c].v->z,
//				s[c].opp,TET(s[c].opp),INDEX(s[c].opp),
//				s[c].plane.w, s[c].plane.x, s[c].plane.y, s[c].plane.z,
//				orientDet[TET(c)], (infiniteV(s[c].v)? dotInf(s[c].plane, s[c].v) : dot(s[c].plane, s[c].v)),
//				s[c].v->sq);
//	}
//
//	void cornerPrint4(int c) {
//		int b,j;
//		b = BASECORNER(c);
//		for (j=0;j<4;j++) cornerPrint(b+j);
//	}


}
