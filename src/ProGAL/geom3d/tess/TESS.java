package ProGAL.geom3d.tess;

import java.util.Locale;

import ProGAL.io.IOToolbox;
import ProGAL.math.Randomization;

public class TESS {


	// =========== From d3.h ===========


	private static final int MAXCOORD = 0x800;
	private static final int COORDMASK = 0x7ff;

	private static final int RANDBIT(){ return (Randomization.randBetween(0,2)); }             // one random bit
	private static final double RANDCOORD(){ return Randomization.randBetween(0.0, MAXCOORD); }
	private static final double RANDOM(int k) { return Randomization.randBetween(0,k ); }            // random number 0..k-1

	private static final int NVERT = 65000; // number of random vertices

	private static final int MAXVERT   =  65010; // n+4
	private static final int MAXTETRA  = 130200; // at least 20*n
	private static final int MAXCORNER = 520800; // 4*MAXTETRA

	private static final int MOD4(int a) { return (a & 3); }
	private static final int TETRA(int corner) { return ((corner) >> 2); }
	private static final int INDEX(int corner) { return (MOD4(corner)); }
	private static final int BASECORNER(int corner) { return ((corner) & 0xFFFFFFFC); }
	private static final int CORNER(int tetra,int index) { return (((tetra)<<2)+(index)); }

	private static final void ASSERT(boolean b, String msg){
		if(!b) System.err.println(msg);
	}

	int nvert; // number of vertices
	TPoint[] vert = new TPoint[MAXVERT];  // vertex coordinates

	// Table-based tetratopes
	TCorner[] s = new TCorner[MAXCORNER]; // per corner entries: v, opp
	TSphere[] sph = new TSphere[MAXTETRA]; // per tetra entries: sphere

	private boolean EQUALV(int i,int j) {
		return (vert[i].x == vert[j].x && vert[i].y == vert[j].y && vert[i].z == vert[j].z);
	}

	private void subtractV(TPoint v) {
		v.xv = v.x - vert[pv].x; v.yv = v.y - vert[pv].y; 
		v.zv = v.z - vert[pv].z; v.sqv = v.sq - vert[pv].sq;
	}

	private double spdot(TSphere sp, TPoint pv){
		return ((sp).x*(pv).x+(sp).y*(pv).y+(sp).z*(pv).z+(sp).sq*(pv).sq+(infiniteV(pv)?0:(sp).w));
	}


	// =========== From d3.c ===========
	//TPoint pv; // vertex index in outermost loop
	int pv;
	int freeTetra; // head for free list for tetrahedra kept in opp[CORNER(tetra,0)]
	int liveTetra; // latest tetra; known to be live. 
	int maxTetra; // AUDIT only

	private void PUSH(int value, int stack) { 
		//		stack##st[++stack##sp] = value; 
		//		if (stack##max < stack##sp) { stack##max = stack##sp; 
		//		if (stack##max >= STACKMAX) { 
		//			printf("ERROR: overflow stack %s pushing %d",  stack##n, value); exit(1); } 
		//		} 
		stackst[stack][++stacksp[stack]] = value;
		if(stackmax[stack]>=STACKMAX){
			System.err.printf("ERROR: overflow stack %s pushing %d",  stack, value);
			System.exit(-1);
		}
	}


	private static final int STACKMAX = 512;
	private final int POP(int stack){ return (stackst[stack][stacksp[stack]--]); }
	private final boolean isEMPTY(int stack) { return (stacksp[stack] < 0); }
	private final void stkINIT(int stack) { stacksp[stack] = -1; }
	//	private final void stkDECLARE(int stack) { stackmax[stack] = -1, stackst[stack][STACKMAX]; } 

	// stacks used in inserting pv
	private final int dfs  = 0;//	stkDECLARE(dfs);  // DFS stack to find dead tetras
	private final int idfs = 1;//	stkDECLARE(idfs); // DFS stack for tetras adj to infinite vertex (>4*30)
	private final int nbhr = 2;//	stkDECLARE(nbhr); // stack for dead corners with live neighbors
	private final int kill = 3;//	stkDECLARE(kill); // stack for base corners of tetras to recycle
	private final int[] stacksp = new int[4];
	private final int[] stackmax = {-1,-1,-1,-1};
	private final int[][] stackst = new int[4][STACKMAX];
	private final String[] stackn = {"dfs", "idfs", "nbhr", "kill" };


	// used during incremental insertion, and for printing/auditing
	private final int[] active = new int[MAXTETRA]; // flag -1 unused, 0 dead, 1 alive
	private boolean DEAD(int p){ return (active[p] <= 0); }// is this a dead or killed tetrahedron?
	private void KILL(int p) { active[p] = -1; } // kill tetrahedron

	// when keeping statistics...
	private int maxLocate =0, locateSideCnt = 0, sphereCnt = 0, startTetraCnt=0, freeTetraCnt=0, inSphereCnt=0;
	private int  locfail = 0, lochard = 0, locsteps = 0, rb = 0;

	//const int mod4[] = {0,1,2,3,0,1,2,3};
	private final int[][] offset = { 
			/*0*/{0,0,0,0}, /*1*/{1,1,1,-3}, /*2*/{2,2,-2,-2}, /*3*/{3,-1,-1,-1}};
	private int INCREMENT(int c) { return (c+offset[1][INDEX(c)]); }
	// drop[i] contains new vertex order after vertex i is dropped and replaced by pv on same side.
	// offdr[i] contains drop(i)-index(i)
	// invdrop[i][k] = j whenever drop[i][j] = k.  4s signal i=k; bad because i is dropped. 
	private final int[][] drop = {{2,1,3}, {0,2,3}, {1,0,3}, {0,1,2}};
	private final int[][] offdr = {{2,1,3}, {-1,1,2}, {-1,-2,1}, {-3,-2,-1}};
	private final int[][] invdrop = {{4,1,0,2}, {0,4,1,2}, {1,0,4,2}, {0,1,2,4}};

	private void STARTTETRA(){ 
		liveTetra = freeTetra; 
		freeTetra = s[CORNER(liveTetra,0)].opp; 
		startTetraCnt++; 
		active[liveTetra] = 1;
		if (maxTetra <= liveTetra) { 
			maxTetra = liveTetra+1; 
			if (liveTetra >= MAXTETRA) { System.err.printf("AUDIT: %d > MAXTETRA\n", liveTetra); System.exit(1); }
		}
	} // AUDIT /**/

	private void FREETETRA(int p) 
	{ 
		freeTetraCnt++;
		s[CORNER(p,0)].opp = freeTetra; freeTetra = p; /**/ active[p]=0; 
		setCornerVC(CORNER(p,1),null,MAXCORNER); 
		setCornerVC(CORNER(p,2),null,MAXCORNER); 
		setCornerVC(CORNER(p,3),null,MAXCORNER); 
		sph[p].w = 1.0; sph[p].x = sph[p].y = sph[p].z = sph[p].sq = 0.0; /* make sphere empty */
	}

	/* Tetrahedra are groups of four corners in order of increasing 
	   vertex index, except that the first two may be swapped to ensure 
	   that the orientation determinant is positive.  
	   I.e., take the lex smallest alternating permutation with positive sign. 
	   There must always be an odd number of swaps between two permutations; 
	   we swap the first two if necessary to achieve this. 

	   The table indoff has the index offset for where 
	   the vertex of index i will be found in the tetrahedron opposite c.  
	   Vertex s[BASECORNER(c)+i].v will be at s[s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]].v 
	   (except that s[c].v and s[s[c].opp]].v are different, and opposite sides of common pl).
	   Note that indoff[*][j] uses each offset -j:4-j exactly once, 
	   that indoff[i][j][i] = 0 for all possible i and j (so s[c].v is opposite of s[s[c].opp].v), 
	   and that the tetra will never change, TETRA(s[c].opp) == TETRA(CORNERINOPP(i,c)).
	 */

	// Lowercase combinations occur when b<a; xxxx configurations don't occur.
	// (Is this a problem?)
	//
	// Since we always add points at the last position, we create tetras with V falling at 4, 
	// which has no ambiguity---we never need to decide if A<B or b<a for construction.
	// Thereafter, the indices are right, and we can just follow.
	//
	// This is the table of where the indices go
	// Replacing  \ With V falling at position c.opp:
	//   corner c  \0ABCD 1ABCD 2ABCD 3ABCD 
	//	       A    0: xxxx  BVCD  CBVD  BCDV
	//	       B    1: VACD  xxxx  ACVD  CADV 
	//	       C    2: vBAD  AVBD  BAVD  ABDV 
	//	       D    3: Vabc  bvac  ABVC  BACV
	//
	// Offsets from c.opp  I=-1, Z=-2, B = -3, H = -4
	// Replacing \0ABCD 1ABCD 2ACBD 3ABCD 
	//	     A    0: xxxx  0I12  0IZ1  0BZI 
	//	     B    1: 1023  xxxx  Z0I1  Z0BI 
	//	     C    2: 1203  I102  IZ01  BZ0I 
	//	     D    3: 1230  1I20  ZI10  ZBI0 
	// Get these by subtracting 01234 from columns in order


	private int CORNERINOPP(int i, int c) { return (s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]); }
	private static final int[][][] indoff =
	{/*        0ABCD       1ABCD         2ACBD          3ABCD     */
		/* 0: */ {{5,5,5,5},  { 0,-1, 1, 2},  { 0,-1,-2, 1},  {  0,-3,-2,-1}},
		/* 1: */ {{1,0,2,3},  { 5, 5, 5, 5},  {-2, 0,-1, 1},  { -2, 0,-3,-1}},
		/* 2: */ {{2,1,0,3},  {-1, 1, 0, 2},  {-1,-2, 0, 1},  { -3,-2, 0,-1}},
		/* 3: */ {{1,2,3,0},  { 1,-1, 2, 0},  {-2,-1, 1, 0},  { -2,-3,-1, 0}}
	};

	// set vert coordinates
	private void setVert(int i, double xx, double yy, double zz, double pp) {
		if (i >= MAXVERT) { System.err.printf("AUDIT: %d > MAXVERT\n", i); System.exit(1); }
		if(vert[i]==null) vert[i] = new TPoint();
		vert[i].x=xx; vert[i].y=yy; vert[i].z=zz;
		vert[i].p=pp; vert[i].sq = (double) (xx*xx + yy*yy + zz*zz);
	}

	private void setInfVert(int i, double xx, double yy, double zz, double qq) {
		if(vert[i]==null) vert[i] = new TPoint();
		vert[i].x=vert[i].xv=xx; vert[i].y=vert[i].yv=yy; vert[i].z=vert[i].zv=zz;
		vert[i].p=0; vert[i].sq = vert[i].sqv = qq;
	}

	private boolean infiniteV(TPoint pv) { return ((pv) == vert[0]); } // first point is at infinity
	private boolean infiniteP(int p) { return (sph[p].sq == 0.0); } // if tetra uses infinite point

	// set corner's vertex and opposite in tetrahedron structure
	//void setCornerVC(int c, int vv, int op) {
	private void setCornerVC(int c, TPoint vv, int op) { s[c].v = vv; s[c].opp = op; }
	private void setCornerPairV(int c, int op, TPoint vv, TPoint ov) { setCornerVC(c, vv, op); setCornerVC(op, ov, c); }
	private void setCornerVCN(int c, TPoint vv, int op) { setCornerVC(c, vv, op); s[op].opp = c; } // set corner & adjust nbhr opp
	// set four corners vertices; pv will be fifth

	static class AtomType {
		int boxno;
		int atomno, atype;
		String aname;
		int x, y, z, sq;  // coordinates: 3 spatial 
	}

	void initVerts() {
		int j;
		int[] hashtb = new int[MAXCOORD];
		int[] nexttb = new int[NVERT];

		for (j = 0; j<MAXCOORD; j++) { hashtb[j] = 0; }
		nvert = 0;
		nexttb[nvert] = 0;
		setInfVert(nvert++, 0,0,0, 1); // point at +infinity
//		while (nvert < NVERT) {
//			setVert(nvert, RANDCOORD(), RANDCOORD(), RANDCOORD(), (double) nvert);
//			j = hashtb[(int) vert[nvert].z];
//			while ((j>0) && !EQUALV(nvert,j))
//				j = nexttb[j];
//			if (j<=0) {
//				nexttb[nvert] = hashtb[(int) vert[nvert].z];
//				hashtb[(int) vert[nvert].z] = nvert;
//				nvert++;
//			}
//		}

		setVert(nvert, 0, 1, 2, (double)nvert);
		j = hashtb[(int)vert[nvert].z];
		while ((j>0) && !EQUALV(nvert,j))
			j = nexttb[j];
		if (j<=0) {
			nexttb[nvert] = hashtb[(int) vert[nvert].z];
			hashtb[(int) vert[nvert].z] = nvert;
			nvert++;
		}
		setVert(nvert, 2, 2, 2, (double)nvert);
		j = hashtb[(int)vert[nvert].z];
		while ((j>0) && !EQUALV(nvert,j))
			j = nexttb[j];
		if (j<=0) {
			nexttb[nvert] = hashtb[(int) vert[nvert].z];
			hashtb[(int) vert[nvert].z] = nvert;
			nvert++;
		}
		setVert(nvert, 1, 1, 2, (double)nvert);
		j = hashtb[(int)vert[nvert].z];
		while ((j>0) && !EQUALV(nvert,j))
			j = nexttb[j];
		if (j<=0) {
			nexttb[nvert] = hashtb[(int) vert[nvert].z];
			hashtb[(int) vert[nvert].z] = nvert;
			nvert++;
		}
		setVert(nvert, 0, 1, 0, (double)nvert);
		j = hashtb[(int)vert[nvert].z];
		while ((j>0) && !EQUALV(nvert,j))
			j = nexttb[j];
		if (j<=0) {
			nexttb[nvert] = hashtb[(int) vert[nvert].z];
			hashtb[(int) vert[nvert].z] = nvert;
			nvert++;
		}
		setVert(nvert, 2, 0, 2, (double)nvert);
		j = hashtb[(int)vert[nvert].z];
		while ((j>0) && !EQUALV(nvert,j))
			j = nexttb[j];
		if (j<=0) {
			nexttb[nvert] = hashtb[(int) vert[nvert].z];
			hashtb[(int) vert[nvert].z] = nvert;
			nvert++;
		}
		setVert(nvert, 1, 2, 1, (double)nvert);
		j = hashtb[(int)vert[nvert].z];
		while ((j>0) && !EQUALV(nvert,j))
			j = nexttb[j];
		if (j<=0) {
			nexttb[nvert] = hashtb[(int) vert[nvert].z];
			hashtb[(int) vert[nvert].z] = nvert;
			nvert++;
		}
//		while (nvert < 5) {
//			setVert(nvert, RANDCOORD(), RANDCOORD(), RANDCOORD(), (double) nvert);
//			j = hashtb[(int) vert[nvert].z];
//			while ((j>0) && !EQUALV(nvert,j))
//				j = nexttb[j];
//			if (j<=0) {
//				nexttb[nvert] = hashtb[(int) vert[nvert].z];
//				hashtb[(int) vert[nvert].z] = nvert;
//				nvert++;
//			}
//		}
	}

	void cornerPrint(int c) {
		System.out.printf("%%%3d(%2d,%1d)=", c, TETRA(c), INDEX(c));//, s[c].v-vert);
		//if ((s[c].v >= vert) && (s[c].v < vert[MAXVERT))
		if(!(s[c].v==null))
			System.out.printf("(%d %5.0f %5.0f %5.0f) opp:%4d(%3d,%2d) \n", (infiniteV(s[c].v)?0:1),
					s[c].v.x, s[c].v.y, s[c].v.z, s[c].opp,TETRA(s[c].opp),INDEX(s[c].opp));
		System.out.flush();
	}

	void cornerPrint4(int c) {
		int b,j;
		TSphere sp;
		b = BASECORNER(c);
		sp = sph[TETRA(b)];
		System.out.printf("disp('Sphere(%d) = <%5.0f %5.0f %5.0f %5.0f %5.0f>')\n", TETRA(b),
				sp.w,sp.x,sp.y,sp.z,sp.sq);
		for (j=0;j<4;j++) cornerPrint(b+j);
	}

	double InSpherev(TSphere sp, TPoint pv) {
		double d = spdot(sp, pv); // Return true if inside (==negative)
		inSphereCnt++;
		System.out.printf("%% %f=InSphere(<%5.0f %5.0f %5.0f %5.0f %5.0f>*( %5.0f %5.0f %5.0f %5.0f))\n", d,
				     sp.w,sp.x,sp.y,sp.z,sp.sq, pv.x,pv.y,pv.z,pv.sq);

		return d; // perturb those on sphere to inside
	}

	double detdiv(double a, double b, double c, double d, double e) {
		// return integer 52 bits of det[a b; c d] / e, even if ad-bc has more than 52 bits.
		// e is assumed to divide evenly, and a=>0 and c <0.
		double na, nc, n; // n reused as temporary (so don't declare int).

		/*  if (!(a>=0 && c <0))
		      printf("AUDIT: a:%18.0f, b:%18.0f, c:%18.0f, d:%18.0f, e:%18.0f\n", a,b,c,d,e);/**/
		while (a > 0 && ((b > 0 && d < 0) || (b < 0 && d > 0)) ) { // while b,d opposite sign
			na = Math.floor(b/a); nc = Math.floor(d/c);
			n = Math.ceil((na+nc)/2);
			b = b - n*a; d = d - n*c; // try a column operation
			if (na == nc)             // didn't change the signs, so do row op
				if (-2*c <= a) {
					a = a+c; b = b+d;
				} else
					if (-c >= a) {
						c = a+c; d = b+d;
						if (c == 0.0) break;
					} else {
						n = a; a = -c; c = -(c+n);
						n = b; b = d; d = (d+n);
					}
		}
		return Math.floor((a*d-c*b)/e+0.5);
	}

	void makeSphereLD(int dead, TSphere sp) {
		// from the dead corner opp a live tetra, we make the sphere for pv.
		TSphere spL, spD; // live and dead spheres
		double L, D; // pv is in or on  live, inside dead: L >=0, D <0.
		double denom; // common denominator is spL*dead.v >0.

		sphereCnt++;
		spD = sph[TETRA(dead)]; // dead tetra
		spL = sph[TETRA(s[dead].opp)]; // live tetra
		denom = spdot(spL, s[dead].v);
		if (denom <= 0) {
			System.out.printf("denom %18.0f <= 0, which means that s[%d].v was in or on the live sphere!\n", denom, dead);
			cornerPrint4(s[dead].opp);
			cornerPrint4(dead);
			denom = 1;
		}
		D = spdot(spD, vert[pv]); D /= denom; // special case 0?
		L = spdot(spL, vert[pv]); D /= denom;
		sp.w = Math.floor((double)(L*spD.w  - D*spL.w) +0.5);
		sp.x = Math.floor((double)(L*spD.x  - D*spL.x) +0.5);
		sp.y = Math.floor((double)(L*spD.y  - D*spL.y) +0.5);
		sp.z = Math.floor((double)(L*spD.z  - D*spL.z) +0.5);
		sp.sq= Math.floor((double)(L*spD.sq - D*spL.sq)+0.5);
	}

	void makeSphereP(int dead, TSphere sp) {
		// from the dead corner opp a live tetra, we make the sphere for pv.
		int j;
		TSphere spL, spD; // live and dead spheres
		double L, D; // pv is outside or on live, inside dead: L >=0, D <0.
		double denom; // common denominator is spL*dead.v >0.

		sphereCnt++;
		spD = sph[TETRA(dead)]; // dead tetra
		spL = sph[TETRA(s[dead].opp)]; // live tetra
		D = spdot(spD, vert[pv]);
		L = spdot(spL, vert[pv]);
		ASSERT(L>=0 && D<0, "Unexpected L or D in makeSphereP");
		denom = spdot(spL, s[dead].v);
		ASSERT(denom > 0, "non-positive denom");
		sp.w = detdiv(L, spL.w, D, spD.w, denom);
		sp.x = detdiv(L, spL.x, D, spD.x, denom);
		sp.y = detdiv(L, spL.y, D, spD.y, denom);
		sp.z = detdiv(L, spL.z, D, spD.z, denom);
		sp.sq = detdiv(L, spL.sq, D, spD.sq, denom);

		ASSERT(spdot(sp,vert[pv]) == 0, "pv not on new sphere.");
		j = INDEX(dead);
		ASSERT(spdot(sp,s[dead+offset[1][j]].v) == 0, "s[dead+offset[1][j]].v not on new sphere.");
		ASSERT(spdot(sp,s[dead+offset[2][j]].v) == 0, "s[dead+offset[2][j]].v not on new sphere.");
		ASSERT(spdot(sp,s[dead+offset[3][j]].v) == 0, "s[dead+offset[3][j]].v not on new sphere.");
		ASSERT(Math.abs(sp.w)+Math.abs(sp.x)+Math.abs(sp.y)+Math.abs(sp.z)+Math.abs(sp.sq)>0, "New sphere is zero.");
	}


	void makeSphereV(TSphere sp, TPoint v0, TPoint v1, TPoint v2, TPoint pv) {
		double xy, xz, xs, yz, ys, zs, ts; // 2x2 minors
		// make sphere: only v0 or v1 may be infinte. 
		sphereCnt++;
		if(!infiniteV(v0)) subtractV(v0);
		if(!infiniteV(v1)) subtractV(v1);
		subtractV(v2);
		xy = ((v0).xv*(v1).yv-(v0).yv*(v1).xv);//xy = DET2(v0,v1,xv,yv);
		xz = ((v0).xv*(v1).zv-(v0).zv*(v1).xv);//xz = DET2(v0,v1,xv,zv);
		yz = ((v0).yv*(v1).zv-(v0).zv*(v1).yv);//yz = DET2(v0,v1,yv,zv);
		xs = ((v0).xv*(v1).sqv-(v0).sqv*(v1).xv);//xs = DET2(v0,v1,xv,sqv);
		ys = ((v0).yv*(v1).sqv-(v0).sqv*(v1).yv);//ys = DET2(v0,v1,yv,sqv);
		zs = ((v0).zv*(v1).sqv-(v0).sqv*(v1).zv);//zs = DET2(v0,v1,zv,sqv);

		sp.x  = -v2.yv*zs +v2.zv*ys -v2.sqv*yz;
		sp.y  =  v2.xv*zs -v2.zv*xs +v2.sqv*xz;
		sp.z  = -v2.xv*ys +v2.yv*xs -v2.sqv*xy;
		sp.sq =  v2.xv*yz -v2.yv*xz +v2.zv*xy;
		sp.w  = -pv.x*sp.x -pv.y*sp.y -pv.z*sp.z -pv.sq*sp.sq;
	}
	// when we initialize, this is what we fill in. 
	//const int initialopp[] = {11,5,15,21,25, 1,10,16,20,26, 6,0,17,22,27, 2,7,12,23,28, 8,3,13,18,29, 4,9,14,19,24};
	private static final int[][] initialopp = {
		{CORNER(1,1), CORNER(2,0), CORNER(3,1), CORNER(4,0)},
		{CORNER(2,1), CORNER(0,0), CORNER(3,0), CORNER(4,1)},
		{CORNER(0,1), CORNER(1,0), CORNER(3,2), CORNER(4,2)},
		{CORNER(1,2), CORNER(0,2), CORNER(2,2), CORNER(4,3)},
		{CORNER(0,3), CORNER(1,3), CORNER(2,3), CORNER(3,3)}
	};

	private static final int[][][] initialv = {
		{{1,2,3,4}, {2,0,3,4}, {0,1,3,4}, {1,0,2,4}, {0,1,2,3}},
		{{0,2,3,4}, {2,1,3,4}, {1,0,3,4}, {0,1,2,4}, {1,0,2,3}}
	};

	void initTetras()
	{
		for(int i=0;i<s.length;i++) s[i] = new TCorner();
		for(int i=0;i<sph.length;i++) sph[i] = new TSphere();
		
		int p,j;
		int last = -1; /* set up free list of tetrahedra */
		double d;

		freeTetra = MAXTETRA;
		do {
			freeTetra--;
			active[freeTetra]=0; // KILL(freeTetra);
			s[CORNER(freeTetra,0)].opp = last;
			last = freeTetra;
		} while (freeTetra > 5);

		active[4] = 2;
		makeSphereV(sph[4], vert[0], vert[1], vert[2], vert[3]);
		d = spdot(sph[4], vert[4]); // if d<0, then we need to swap

		if (d == 0.0) {
			System.err.printf("ERROR: Need first five vertices to be in general position"); System.exit(1);
		}

		if (d < 0) {
			sph[4].w = -sph[4].w; sph[4].x = -sph[4].x; sph[4].y = -sph[4].y; sph[4].z = -sph[4].z;
			sph[4].sq = -sph[4].sq;
		}

		for (p=0; p<5; p++) {
			for (j=0; j<4; j++) {       // pay attention to orientation when assigning vertices
				setCornerVC(CORNER(p,j), vert[initialv[d<0?1:0][p][j]], initialopp[p][j]);
			}
			makeSphereV(sph[p], s[CORNER(p,0)].v, s[CORNER(p,1)].v, s[CORNER(p,2)].v, s[CORNER(p,3)].v);
			active[p] = 1;
			// swap first two if d<0. 
			ASSERT(spdot(sph[p], vert[p + (d<0?1:0)*(p<2?1:0)*(1-2*p)]) > 0, "Somehow vertp is in or on sphere p in init.");
		}
		liveTetra = 4;
		maxTetra = 5;
	}

	void auditCornersAux(int SphereCheck) {
		int p, b,c, i,j,k, guard;
		double d;
		TSphere sp;
		TPoint vv;

		for (p = 0; p < maxTetra; p++) {
			if (DEAD(p)) continue; // don't audit tetras on free list
			b = CORNER(p,0);
			for (c=b; c < CORNER(p,4); c++) { // per corner checks
				i = s[c].opp; // check opposite
				if (s[i].opp != c) { System.err.printf("%%AUDIT: wrong opp.opp \n"); cornerPrint(c); cornerPrint(i); }
				if (s[c].v == s[i].v){
					//					printf("%%AUDIT: Same vertex  s[%d(%d,%d)].v %d == opp[%d(%d,%d)].v %d \n",
					//							c,TETRA(c), INDEX(c), s[c].v-vert, i, TETRA(i), INDEX(i), s[i].v-vert);
					System.out.printf("%%AUDIT: Same vertex  s[%d(%d,%d)].v == opp[%d(%d,%d)].v \n",
							c,TETRA(c), INDEX(c), i, TETRA(i), INDEX(i));
					cornerPrint4(c); cornerPrint4(i);
				}

				for (j = 0; j < 4; j++)
					if ((j != INDEX(c)) && (CORNERINOPP(j,c)<3))
						if (s[BASECORNER(c)+j].v != s[CORNERINOPP(j,c)].v) {
							//							printf("%%AUDIT:Bad auditCornerInOpp(%d,%d) = %d  since vertex %d != %d\n",
							//									j,c,CORNERINOPP(j,c), s[BASECORNER(c)+j].v-vert, s[CORNERINOPP(j,c)].v-vert);
							System.out.printf("%%AUDIT:Bad auditCornerInOpp(%d,%d) = %d  since vertices dont match\n",
									j,c,CORNERINOPP(j,c));
							cornerPrint4(BASECORNER(c)); cornerPrint4(BASECORNER(CORNERINOPP(j,c)));
							break;
						}

				for (j = 0; j<4; j++) {
					k = CORNERINOPP(j,c);
					if ((TETRA(k) != TETRA(s[c].opp)) || (s[b+j].v != s[k].v)) {
						if (TETRA(k) != TETRA(s[c].opp)) {
							System.out.printf("%%AUDIT: CORNERINOPP(%d,%d) ==>%d: Accessing [%d:%d][%d:%d][%d]\n",
									j, c, k-s[c].opp, c, INDEX(c), s[c].opp, INDEX(s[c].opp), j);
							k = s[c].opp;
						}
						else {
							if (b+j == c)
								if (k == s[c].opp) continue; // these vertices are supposed to differ; don't flag them
								else System.out.printf("%%AUDIT: CORNERINOPP(%d,%d) says %d(%d,%d) and %d(%d,%d) shouldn't happen\n",
										j, c, b+j, TETRA(b), j, k, TETRA(k), INDEX(k));
							else {
								System.out.printf("%%AUDIT: CORNERINOPP(%d,%d) says %d(%d,%d) and %d(%d,%d) should agree\n",
										j, c, b+j, TETRA(b), j, k, TETRA(k), INDEX(k));
								cornerPrint4(c);
								cornerPrint4(b+j);
								cornerPrint4(k);
								ASSERT(TETRA(s[c].opp) == TETRA(CORNERINOPP(i,c)), "CORNERINOP screws up tetras");

							}
						}
						cornerPrint4(b); cornerPrint4(BASECORNER(k)); break;
					}
				}

				// check sphere opposite corner c
				if (SphereCheck>0)        // check sphere opposite corner
				{
					k = s[c].opp;
					sp = sph[TETRA(k)];

					d = spdot(sp, s[c].v);
					if (d < 0) {
						//						printf("disp('AUDIT: corner %d v%d in or on sphere %d(%d) =%5.0f');\n",
						//								c, s[c].v-vert, TETRA(k), k, d);
						System.out.printf("disp('AUDIT: corner %d in or on sphere %d(%d) =%5.0f');\n",
								c, TETRA(k), k, d);
						cornerPrint4(k);
					}
				}
			}

			if (SphereCheck>0)  // check sphere sqs (orient dets) 
			{
				sp = sph[p]; // only spheres using pt at infty have sq==0; none have sq < 0.
				b = CORNER(p,0);

				if (sp.sq < 0 || (sp.sq ==0 && !infiniteV(s[b].v) && !infiniteV(s[b+1].v))) {
					System.out.printf("disp('AUDIT: sq<=0 in tetra %d(%d) =%5.0f'); \n", CORNER(p,0), p, sp.sq);
					cornerPrint4(b);
					System.out.printf("DetCheckH([");
					for (k = 0; k < 5; k++)
						System.out.printf(" %d %5.0f %5.0f %5.0f %5.0f; ",
								(infiniteV(s[b+k].v)?0:1), s[b+k].v.x, s[b+k].v.y, s[b+k].v.z, s[b+k].v.sq);
					System.out.printf("]);\n");
				} /**/

			}
		}
	}

	private void printVerts() {
		int i;
		for (i = 0; i<nvert; i++)
			System.out.printf("%4d: (%5.0f; %5.0f %5.0f %5.0f; %10f)\n",
					i, vert[i].p, vert[i].x, vert[i].y, vert[i].z, vert[i].sq);
	}

	private void printCorners() {
		int p;

		for (p = 0; p < maxTetra; p++)
			if (!DEAD(p))
				cornerPrint4(CORNER(p,0));
		auditCornersAux(1); // check spheres, too
	}

	public static void main(String[] args){
		Locale.setDefault(Locale.ENGLISH);
		new TESS().main();
	}
	void main() {
		int p; // tetra
		int b, c, newb; // corners
		int nc, ni, dead, jdead;
		int i, j, k, off; // indices
		int guard, rb;
		double I1, I2; // Insphere values for locate
		double d;
		//		TSphere s1, s2; // spheres for locate
		int s1, s2 = 0;
		int c1, c2; // corners for locate
		TPoint v0, v1, v2; // pointers to vertices

		int ii; // timing
		long tic = System.nanoTime(), toc;


		maxLocate = locateSideCnt = sphereCnt = startTetraCnt = freeTetraCnt= inSphereCnt=0;
		locfail = lochard = locsteps = rb = 0;
		//		  if (args.length > 1)
		//		    readHierarchy(args[0]);
		//		  else
		initVerts();
		for (ii = 0; ii<5; ii++) // print first few points
			System.out.printf("%4d: (%5.0f; %5.0f %5.0f %5.0f; %10f)\n", ii, vert[ii].p, vert[ii].x, vert[ii].y, vert[ii].z, vert[ii].sq);
		for (ii = nvert-5; ii<nvert; ii++) // print last few points
			System.out.printf("%4d: (%5.0f; %5.0f %5.0f %5.0f; %10f)\n", ii, vert[ii].p, vert[ii].x, vert[ii].y, vert[ii].z, vert[ii].sq);
		//#endif
		initTetras();
		for (pv = 5; pv < nvert; pv++) { // incrementally insert pv
			//      auditCorners(1); // audit, and check spheres, too
			//LOCATE: find some tetrahedron with sphere strictly containing pv
			guard = 2*(pv+4); // prevent infinite loops  
			//      rb = 0; /**/

			s1 = liveTetra; // live sphere to start search
			c1 = CORNER(liveTetra,0); // corner in liveTetra
			locateSideCnt++;
			I1 = InSpherev(sph[s1], vert[pv]); // Check if inside first one.
			if ( (I1 < 0) || (I1 == 0 && sph[s1].sq > 0) ){// found already
				p = liveTetra;

				k = 1;
				/*        printf("g:%5d, %2d ",rb,rb);/**/
			} else {
				//      rb = 0;
				while( (--guard)!=0 ) {
					c2 = s[c1].opp; s2 = TETRA(c2); // nhbr corner, sphere,
					locateSideCnt++;
					I2 = InSpherev(sph[s2], vert[pv]);
					if ( (I2 < 0)  || (I2 == 0 && sph[s2].sq > 0) ) break; // found one!
					d = sph[s2].sq * I1 - sph[s1].sq * I2;
					if (d < 0 || (d==0 && rb++!=0 && RANDBIT()==1)) // if on I1 side (flip coin for equals, to break out of infinite loops.)
						c1 = INCREMENT(c1);
					else {// on I2 side
						//        ASSERT(s2.sq * I1 > s1.sq * I2, "On plane in locate");
						c1 = INCREMENT(c2); s1 = s2; I1 = I2;
					}
				}
				ASSERT(guard!=0,"infinite loop in locate");
				if (guard!=0) { // locate the hard way
					for (p = 0; p < maxTetra; p++)
						if (!DEAD(p)) {
							lochard++;
							d = InSpherev(sph[p], vert[pv]);
							if ( (d < 0) || (d == 0 && sph[p].sq > 0))
								break;
						}
//					System.out.printf("AUDIT: Locate %d the hard way to find %d\n", pv+vert, p); /**/
//					System.out.printf("AUDIT: Locate <something> the hard way to find %d\n", p); /**/
					if (p >= maxTetra) {
						locfail++;
						System.out.printf("AUDIT: Still can't locate! <something> ( %5.0f %5.0f %5.0f %5.0f)\n",
								vert[pv].x,vert[pv].y,vert[pv].z,vert[pv].sq); /**/
						continue; // next pv; skip this one.
					}
				}
				else {
					p = s2;
					k = 2*(pv+4)-guard;
					//      printf("g:%5d, %2d ",k,rb);/**/
					if (k > maxLocate) {
						maxLocate = k;
						/*      printf("AUDIT: max locate %d\n", maxLocate); /**/
					}
				}
				k = 2*(pv+4)-guard;
				locsteps += k;
				/*  printf("g:%5d, %2d ",k,rb);/**/
			}
			// Tetrahedra containing pv are "dead", and are pushed onto kill stack.  
			// We use DFS with stack pst to find them and kill them
			// At live-dead boundary, we save dead tetras on stack nbhr,
			//  then make new tetras and hook in to live by setting the last opp pointer.
			//
			// Invariants/operations: Tetrahedron p is marked alive or dead on first visit.
			//    Corner c is pushed on stack when TETRA(s[c].opp) is marked dead. 
			//
			// On termination, stack nbhr contains dead corners with live neighbors 
			//    that have new tetras (so s[nhbr].opp != s[s[nhbr].opp].opp temporarily.)
			//    Stack kill contains old tetrahedra for final recycling.

			stkINIT(dfs);  // DFS stack holds corners opposite dead tetras
			stkINIT(idfs); // iDFS stack holds corners opposite infinite tetras with pv on bdry
			//    (special case: dead, but don't propagate)
			stkINIT(nbhr); // stack for dead corners with live nbhr tetras
			stkINIT(kill); // stack of dead tetras to recycle
			b = CORNER(p,0);
			PUSH(p, kill); KILL(p); // kill tetra initial p,
			PUSH(s[b++].opp, dfs); // stack neighbors
			PUSH(s[b++].opp, dfs);
			PUSH(s[b++].opp, dfs);
			PUSH(s[b  ].opp, dfs);

			while (!isEMPTY(dfs)) {
				c = POP(dfs); p = TETRA(c);
				/*    printf("::Popping %d with opp %d \n", c, s[c].opp);/**/
				ASSERT(DEAD(TETRA(s[c].opp)), "dfs stack element with non-dead neighbor");
				if (DEAD(p)) continue; // dead already
				d = InSpherev(sph[p], vert[pv]); // Is pv in, out, or on?
				if (d > 0) { // pv is outside, so c is live neighbor of dead opp tetra s[c].opp
					PUSH(s[c].opp, nbhr); // remember old corner, so we can hook tetra into mesh later
					STARTTETRA(); // make new tetrahedron  liveTetra
					newb = CORNER(liveTetra,3); // last corner of new tetra
					setCornerVCN(newb, vert[pv], c); // last corner is pv; also set opposite corner c. Do rest later.
				}
				else
					if (d < 0 || sph[p].sq > 0) { // kill and continue dfs if pv is strictly inside 
						KILL(p); PUSH(p, kill); // kill and stack tetra
						j = INDEX(c);
						PUSH(s[c+offset[1][j]].opp, dfs); // stack neighbors to check
						PUSH(s[c+offset[2][j]].opp, dfs);
						PUSH(s[c+offset[3][j]].opp, dfs);
					}
					else { // d==0 && sph[p] is infinite: handle two special cases
						if (sph[TETRA(s[c].opp)].sq == 0) { // if dead sphere is infinite, too
							PUSH(c, idfs); // then if c stays alive, we make tetra to it (flat, but infinite).
						} else { // dead sphere is finite; kill c and make tetras to neighbors, if they stay alive.
							KILL(p); PUSH(p, kill); // kill and stack tetra
							j = INDEX(c);
							PUSH(s[c+offset[1][j]].opp, idfs); // stack neighbors to check
							PUSH(s[c+offset[2][j]].opp, idfs);
							PUSH(s[c+offset[3][j]].opp, idfs);
						}
					}
			}

			while (!isEMPTY(idfs)) {
				c = POP(idfs); p = TETRA(c);
				/*    printf("::Popping %d with opp %d \n", c, s[c].opp);/**/
				ASSERT(DEAD(TETRA(s[c].opp)), "dfs stack element with non-dead neighbor");
				if (DEAD(p)) continue; // dead already
				ASSERT(DEAD(TETRA(s[c].opp)), "Live corner c should have dead neighbor");
				PUSH(s[c].opp, nbhr); // remember old corner, so we can hook tetra into mesh later
				STARTTETRA(); // make new tetrahedron liveTetra
				newb = CORNER(liveTetra,3); // last corner of new tetra
				setCornerVCN(newb, vert[pv], c); // last corner is pv; also set opposite corner c. Do rest later.
			}
			// Now, we have stack of dead neighbors of live tetras, and we've hooked new tetras to them.  
			while (!isEMPTY(nbhr)) {
				dead = POP(nbhr); jdead = INDEX(dead); //  dead tetra and index of dropped corner.
				/*    printf("--Popped %d(%d)\n", dead, jdead); /**/
				ASSERT(DEAD(TETRA(dead)), "corner on nbhr stack is not dead!?");
				newb = s[s[dead].opp].opp-3; // base of new tetra. 

				// makeSphereP(dead, sph+TETRA(newb)); // use either this or makeSphereV below

				dead -= jdead; // just use base of dead one.
				// new tetra has 0,1,2,3=pv; 
				// corresponding old indices before jdead is dropped: 
				//   drop[j][0],..,drop[j][3], (no corresp to pv)
				j = jdead;
				i = drop[jdead][0]; // old index of new corner 0;
				c = dead+i; // note i = INDEX(C);
				v0 = s[c].v; // copy vertex v0
				nc = s[c].opp; // go to neighbor
				// In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
				// To avoid index calculations, maintain i = INDEX(c), nc = s[c].opp, ni = INDEX(nc).
				while (DEAD(TETRA(nc))) {
					ni = INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
					j = ni; i = ni + off; c = nc + off; nc = s[c].opp; // fix new j, i, c, and try neighbor
				}
				nc = s[nc].opp; // go to new tetra
				ASSERT(s[nc].v == vert[pv], "Expected to find new tetra using pv after walking dead tetras. ");

				setCornerVC(newb, v0, nc-3+invdrop[i][j]); newb++;


				j = jdead;
				i = drop[jdead][1]; // old index of new corner 1;
				c = dead+i; // note i = INDEX(C);
				v1 = s[c].v; // copy vertex v1
				nc = s[c].opp; // go to neighbor
				// In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
				// To avoid index calculations, maintain i = INDEX(c), nc = s[c].opp, ni = INDEX(nc).
				while (DEAD(TETRA(nc))) {
					ni = INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
					j = ni; i = ni + off; c = nc + off; nc = s[c].opp; // fix new j, i, c, and try neighbor
				}
				nc = s[nc].opp; // go to new tetra
				ASSERT(s[nc].v == vert[pv], "Expected to find new tetra using pv after walking dead tetras. ");

				setCornerVC(newb, v1, nc-3+invdrop[i][j]); newb++;


				j = jdead;
				i = drop[jdead][2]; // old index of new corner 2;
				c = dead+i; // note i = INDEX(C);
				v2 = s[c].v; // copy vertex v2
				nc = s[c].opp; // go to neighbor
				// In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
				// To avoid index calculations, maintain i = INDEX(c), nc = s[c].opp, ni = INDEX(nc).
				while (DEAD(TETRA(nc))) {
					ni = INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
					j = ni; i = ni + off; c = nc + off; nc = s[c].opp; // fix new j, i, c, and try neighbor
				}
				nc = s[nc].opp; // go to new tetra
				ASSERT(s[nc].v == vert[pv], "Expected to find new tetra using pv after walking dead tetras. ");

				setCornerVC(newb, v2, nc-3+invdrop[i][j]); newb++;

				c = s[s[dead+jdead].opp].opp;
				ASSERT(v0==s[CORNERINOPP(0,c)].v, "v0 does not line up");
				ASSERT(v1==s[CORNERINOPP(1,c)].v, "v1 does not line up");
				ASSERT(v2==s[CORNERINOPP(2,c)].v, "v2 does not line up");

				makeSphereV(sph[TETRA(newb)], v0, v1, v2, vert[pv]); // use either this or makeSphereP above
			}
			while (!isEMPTY(kill)) {
				p = POP(kill);
				FREETETRA(p);
			}
			/*    printf("\n%% Done with setting spheres for %d\n", pv-vert); /**/
		}
		toc = System.nanoTime();
		System.out.printf("We performed:\n  %d\tinSphere tests\n %d\tplane tests\n %d\ttetrahedra created-\n %d\tfreed =\n %d\ttetrahedra\n %d\tsphere     equations computed.\n",
				inSphereCnt, locateSideCnt, startTetraCnt, freeTetraCnt, startTetraCnt-freeTetraCnt, sphereCnt);
		System.out.printf("   %d\t locations hard\n %d\t locations failed\n %d\t location tests\n  %d\t location randbits\n",
				lochard, locfail, locsteps, rb);
		System.out.printf("On vertex insertion: max locate steps=%d, max killed=%d, max inf=%d, max created=%d, maxdfs=%d\n",
				maxLocate, stackmax[kill], stackmax[idfs], stackmax[nbhr], stackmax[dfs]);
		System.out.printf("On vertex insertion: max locate steps=%d\n", maxLocate);
		System.out.printf("Auditing\n");
		System.out.flush();
		auditCornersAux(1);
		System.out.printf("Done auditing\n");
		System.out.flush();
		System.out.printf("Time: %.1f secs\n", ((double) (toc - tic)) / 1000000000);

	}


}
