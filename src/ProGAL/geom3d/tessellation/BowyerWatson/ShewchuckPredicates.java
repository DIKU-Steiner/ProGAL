package ProGAL.geom3d.tessellation.BowyerWatson;

public class ShewchuckPredicates {
	
	private static ShewchuckPredicates instance;
	public static ShewchuckPredicates getInstance(){
		if(instance==null) instance = new ShewchuckPredicates();
		return instance;
	}
	
	private ShewchuckPredicates(){
		exactinit();
	}
	
	private void exactinit(){
		double half;
		double check, lastcheck;
		int every_other;

		every_other = 1;
		half = 0.5;
		epsilon = 1.0;
		splitter = 1.0;
		check = 1.0;


		do {
			lastcheck = check;
			epsilon *= half;
			if (every_other!=0) {
				splitter *= 2.0;
			}
			every_other = every_other==0?1:0;
			check = 1.0 + epsilon;
		} while ((check != 1.0) && (check != lastcheck));
		splitter += 1.0;


		resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
		ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
		ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
		ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
		o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
		o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
		o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
		iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
		iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
		iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
		isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
		isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
		isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;
	}

	
	private double splitter;
	private double epsilon;

	private double resulterrbound;
	@SuppressWarnings("unused")
	private double ccwerrboundA, ccwerrboundB, ccwerrboundC;
	private double o3derrboundA, o3derrboundB, o3derrboundC;
	@SuppressWarnings("unused")
	private double iccerrboundA, iccerrboundB, iccerrboundC;
	private double isperrboundA, isperrboundB, isperrboundC;


	
	
	double insphere(double[] pa,double[]  pb,double[]  pc,double[]  pd,double[]  pe)
	{
		double aex, bex, cex, dex;
		double aey, bey, cey, dey;
		double aez, bez, cez, dez;
		double aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
		double aexcey, cexaey, bexdey, dexbey;
		double alift, blift, clift, dlift;
		double ab, bc, cd, da, ac, bd;
		double abc, bcd, cda, dab;
		double aezplus, bezplus, cezplus, dezplus;
		double aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
		double cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
		double aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
		double det;
		double permanent, errbound;

		aex = pa[0] - pe[0];
		bex = pb[0] - pe[0];
		cex = pc[0] - pe[0];
		dex = pd[0] - pe[0];
		aey = pa[1] - pe[1];
		bey = pb[1] - pe[1];
		cey = pc[1] - pe[1];
		dey = pd[1] - pe[1];
		aez = pa[2] - pe[2];
		bez = pb[2] - pe[2];
		cez = pc[2] - pe[2];
		dez = pd[2] - pe[2];

		aexbey = aex * bey;
		bexaey = bex * aey;
		ab = aexbey - bexaey;
		bexcey = bex * cey;
		cexbey = cex * bey;
		bc = bexcey - cexbey;
		cexdey = cex * dey;
		dexcey = dex * cey;
		cd = cexdey - dexcey;
		dexaey = dex * aey;
		aexdey = aex * dey;
		da = dexaey - aexdey;

		aexcey = aex * cey;
		cexaey = cex * aey;
		ac = aexcey - cexaey;
		bexdey = bex * dey;
		dexbey = dex * bey;
		bd = bexdey - dexbey;

		abc = aez * bc - bez * ac + cez * ab;
		bcd = bez * cd - cez * bd + dez * bc;
		cda = cez * da + dez * ac + aez * cd;
		dab = dez * ab + aez * bd + bez * da;

		alift = aex * aex + aey * aey + aez * aez;
		blift = bex * bex + bey * bey + bez * bez;
		clift = cex * cex + cey * cey + cez * cez;
		dlift = dex * dex + dey * dey + dez * dez;

		det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

		aezplus = ((aez) >= 0.0 ? (aez) : -(aez));
		bezplus = ((bez) >= 0.0 ? (bez) : -(bez));
		cezplus = ((cez) >= 0.0 ? (cez) : -(cez));
		dezplus = ((dez) >= 0.0 ? (dez) : -(dez));
		aexbeyplus = ((aexbey) >= 0.0 ? (aexbey) : -(aexbey));
		bexaeyplus = ((bexaey) >= 0.0 ? (bexaey) : -(bexaey));
		bexceyplus = ((bexcey) >= 0.0 ? (bexcey) : -(bexcey));
		cexbeyplus = ((cexbey) >= 0.0 ? (cexbey) : -(cexbey));
		cexdeyplus = ((cexdey) >= 0.0 ? (cexdey) : -(cexdey));
		dexceyplus = ((dexcey) >= 0.0 ? (dexcey) : -(dexcey));
		dexaeyplus = ((dexaey) >= 0.0 ? (dexaey) : -(dexaey));
		aexdeyplus = ((aexdey) >= 0.0 ? (aexdey) : -(aexdey));
		aexceyplus = ((aexcey) >= 0.0 ? (aexcey) : -(aexcey));
		cexaeyplus = ((cexaey) >= 0.0 ? (cexaey) : -(cexaey));
		bexdeyplus = ((bexdey) >= 0.0 ? (bexdey) : -(bexdey));
		dexbeyplus = ((dexbey) >= 0.0 ? (dexbey) : -(dexbey));
		permanent = ((cexdeyplus + dexceyplus) * bezplus
				+ (dexbeyplus + bexdeyplus) * cezplus
				+ (bexceyplus + cexbeyplus) * dezplus)
				* alift
				+ ((dexaeyplus + aexdeyplus) * cezplus
						+ (aexceyplus + cexaeyplus) * dezplus
						+ (cexdeyplus + dexceyplus) * aezplus)
						* blift
						+ ((aexbeyplus + bexaeyplus) * dezplus
								+ (bexdeyplus + dexbeyplus) * aezplus
								+ (dexaeyplus + aexdeyplus) * bezplus)
								* clift
								+ ((bexceyplus + cexbeyplus) * aezplus
										+ (cexaeyplus + aexceyplus) * bezplus
										+ (aexbeyplus + bexaeyplus) * cezplus)
										* dlift;
		errbound = isperrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}

		return insphereadapt(pa, pb, pc, pd, pe, permanent);
	}
	double insphereadapt(double[] pa,double[]  pb,double[]  pc,double[]  pd,double[]  pe,double permanent)

	{
		double aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
		double det, errbound;

		double aexbey1, bexaey1, bexcey1, cexbey1;
		double cexdey1, dexcey1, dexaey1, aexdey1;
		double aexcey1, cexaey1, bexdey1, dexbey1;
		double aexbey0, bexaey0, bexcey0, cexbey0;
		double cexdey0, dexcey0, dexaey0, aexdey0;
		double aexcey0, cexaey0, bexdey0, dexbey0;
		double[] ab = new double[4], bc = new double[4], cd = new double[4], da = new double[4], ac = new double[4], bd = new double[4];
		double ab3, bc3, cd3, da3, ac3, bd3;
		double abeps, bceps, cdeps, daeps, aceps, bdeps;
		double[] temp8a = new double[8], temp8b = new double[8], temp8c = new double[8], temp16 = new double[16], temp24 = new double[24], temp48 = new double[48];
		int temp8alen, temp8blen, temp8clen, temp16len, temp24len, temp48len;
		double[] xdet = new double[96], ydet = new double[96], zdet = new double[96], xydet = new double[192];
		int xlen, ylen, zlen, xylen;
		double[] adet = new double[288], bdet = new double[288], cdet = new double[288], ddet = new double[288];
		int alen, blen, clen, dlen;
		double[] abdet = new double[576], cddet = new double[576];
		int ablen, cdlen;
		double[] fin1 = new double[1152];
		int finlength;

		double aextail, bextail, cextail, dextail;
		double aeytail, beytail, ceytail, deytail;
		double aeztail, beztail, ceztail, deztail;

		double bvirt;
		double avirt, bround, around;
		double c;
		double abig;
		double ahi, alo, bhi, blo;
		double err1, err2, err3;
		double _i, _j;
		double _0;

		aex = (double) (pa[0] - pe[0]);
		bex = (double) (pb[0] - pe[0]);
		cex = (double) (pc[0] - pe[0]);
		dex = (double) (pd[0] - pe[0]);
		aey = (double) (pa[1] - pe[1]);
		bey = (double) (pb[1] - pe[1]);
		cey = (double) (pc[1] - pe[1]);
		dey = (double) (pd[1] - pe[1]);
		aez = (double) (pa[2] - pe[2]);
		bez = (double) (pb[2] - pe[2]);
		cez = (double) (pc[2] - pe[2]);
		dez = (double) (pd[2] - pe[2]);

		aexbey1 = (double) (aex * bey); c = (double) (splitter * aex); abig = (double) (c - aex); ahi = c - abig; alo = aex - ahi; c = (double) (splitter * bey); abig = (double) (c - bey); bhi = c - abig; blo = bey - bhi; err1 = aexbey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); aexbey0 = (alo * blo) - err3;
		bexaey1 = (double) (bex * aey); c = (double) (splitter * bex); abig = (double) (c - bex); ahi = c - abig; alo = bex - ahi; c = (double) (splitter * aey); abig = (double) (c - aey); bhi = c - abig; blo = aey - bhi; err1 = bexaey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bexaey0 = (alo * blo) - err3;
		_i = (double) (aexbey0 - bexaey0); bvirt = (double) (aexbey0 - _i); avirt = _i + bvirt; bround = bvirt - bexaey0; around = aexbey0 - avirt; ab[0] = around + bround; _j = (double) (aexbey1 + _i); bvirt = (double) (_j - aexbey1); avirt = _j - bvirt; bround = _i - bvirt; around = aexbey1 - avirt; _0 = around + bround; _i = (double) (_0 - bexaey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - bexaey1; around = _0 - avirt; ab[1] = around + bround; ab3 = (double) (_j + _i); bvirt = (double) (ab3 - _j); avirt = ab3 - bvirt; bround = _i - bvirt; around = _j - avirt; ab[2] = around + bround;
		ab[3] = ab3;

		bexcey1 = (double) (bex * cey); c = (double) (splitter * bex); abig = (double) (c - bex); ahi = c - abig; alo = bex - ahi; c = (double) (splitter * cey); abig = (double) (c - cey); bhi = c - abig; blo = cey - bhi; err1 = bexcey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bexcey0 = (alo * blo) - err3;
		cexbey1 = (double) (cex * bey); c = (double) (splitter * cex); abig = (double) (c - cex); ahi = c - abig; alo = cex - ahi; c = (double) (splitter * bey); abig = (double) (c - bey); bhi = c - abig; blo = bey - bhi; err1 = cexbey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cexbey0 = (alo * blo) - err3;
		_i = (double) (bexcey0 - cexbey0); bvirt = (double) (bexcey0 - _i); avirt = _i + bvirt; bround = bvirt - cexbey0; around = bexcey0 - avirt; bc[0] = around + bround; _j = (double) (bexcey1 + _i); bvirt = (double) (_j - bexcey1); avirt = _j - bvirt; bround = _i - bvirt; around = bexcey1 - avirt; _0 = around + bround; _i = (double) (_0 - cexbey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cexbey1; around = _0 - avirt; bc[1] = around + bround; bc3 = (double) (_j + _i); bvirt = (double) (bc3 - _j); avirt = bc3 - bvirt; bround = _i - bvirt; around = _j - avirt; bc[2] = around + bround;
		bc[3] = bc3;

		cexdey1 = (double) (cex * dey); c = (double) (splitter * cex); abig = (double) (c - cex); ahi = c - abig; alo = cex - ahi; c = (double) (splitter * dey); abig = (double) (c - dey); bhi = c - abig; blo = dey - bhi; err1 = cexdey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cexdey0 = (alo * blo) - err3;
		dexcey1 = (double) (dex * cey); c = (double) (splitter * dex); abig = (double) (c - dex); ahi = c - abig; alo = dex - ahi; c = (double) (splitter * cey); abig = (double) (c - cey); bhi = c - abig; blo = cey - bhi; err1 = dexcey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dexcey0 = (alo * blo) - err3;
		_i = (double) (cexdey0 - dexcey0); bvirt = (double) (cexdey0 - _i); avirt = _i + bvirt; bround = bvirt - dexcey0; around = cexdey0 - avirt; cd[0] = around + bround; _j = (double) (cexdey1 + _i); bvirt = (double) (_j - cexdey1); avirt = _j - bvirt; bround = _i - bvirt; around = cexdey1 - avirt; _0 = around + bround; _i = (double) (_0 - dexcey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - dexcey1; around = _0 - avirt; cd[1] = around + bround; cd3 = (double) (_j + _i); bvirt = (double) (cd3 - _j); avirt = cd3 - bvirt; bround = _i - bvirt; around = _j - avirt; cd[2] = around + bround;
		cd[3] = cd3;

		dexaey1 = (double) (dex * aey); c = (double) (splitter * dex); abig = (double) (c - dex); ahi = c - abig; alo = dex - ahi; c = (double) (splitter * aey); abig = (double) (c - aey); bhi = c - abig; blo = aey - bhi; err1 = dexaey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dexaey0 = (alo * blo) - err3;
		aexdey1 = (double) (aex * dey); c = (double) (splitter * aex); abig = (double) (c - aex); ahi = c - abig; alo = aex - ahi; c = (double) (splitter * dey); abig = (double) (c - dey); bhi = c - abig; blo = dey - bhi; err1 = aexdey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); aexdey0 = (alo * blo) - err3;
		_i = (double) (dexaey0 - aexdey0); bvirt = (double) (dexaey0 - _i); avirt = _i + bvirt; bround = bvirt - aexdey0; around = dexaey0 - avirt; da[0] = around + bround; _j = (double) (dexaey1 + _i); bvirt = (double) (_j - dexaey1); avirt = _j - bvirt; bround = _i - bvirt; around = dexaey1 - avirt; _0 = around + bround; _i = (double) (_0 - aexdey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - aexdey1; around = _0 - avirt; da[1] = around + bround; da3 = (double) (_j + _i); bvirt = (double) (da3 - _j); avirt = da3 - bvirt; bround = _i - bvirt; around = _j - avirt; da[2] = around + bround;
		da[3] = da3;

		aexcey1 = (double) (aex * cey); c = (double) (splitter * aex); abig = (double) (c - aex); ahi = c - abig; alo = aex - ahi; c = (double) (splitter * cey); abig = (double) (c - cey); bhi = c - abig; blo = cey - bhi; err1 = aexcey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); aexcey0 = (alo * blo) - err3;
		cexaey1 = (double) (cex * aey); c = (double) (splitter * cex); abig = (double) (c - cex); ahi = c - abig; alo = cex - ahi; c = (double) (splitter * aey); abig = (double) (c - aey); bhi = c - abig; blo = aey - bhi; err1 = cexaey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cexaey0 = (alo * blo) - err3;
		_i = (double) (aexcey0 - cexaey0); bvirt = (double) (aexcey0 - _i); avirt = _i + bvirt; bround = bvirt - cexaey0; around = aexcey0 - avirt; ac[0] = around + bround; _j = (double) (aexcey1 + _i); bvirt = (double) (_j - aexcey1); avirt = _j - bvirt; bround = _i - bvirt; around = aexcey1 - avirt; _0 = around + bround; _i = (double) (_0 - cexaey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cexaey1; around = _0 - avirt; ac[1] = around + bround; ac3 = (double) (_j + _i); bvirt = (double) (ac3 - _j); avirt = ac3 - bvirt; bround = _i - bvirt; around = _j - avirt; ac[2] = around + bround;
		ac[3] = ac3;

		bexdey1 = (double) (bex * dey); c = (double) (splitter * bex); abig = (double) (c - bex); ahi = c - abig; alo = bex - ahi; c = (double) (splitter * dey); abig = (double) (c - dey); bhi = c - abig; blo = dey - bhi; err1 = bexdey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bexdey0 = (alo * blo) - err3;
		dexbey1 = (double) (dex * bey); c = (double) (splitter * dex); abig = (double) (c - dex); ahi = c - abig; alo = dex - ahi; c = (double) (splitter * bey); abig = (double) (c - bey); bhi = c - abig; blo = bey - bhi; err1 = dexbey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dexbey0 = (alo * blo) - err3;
		_i = (double) (bexdey0 - dexbey0); bvirt = (double) (bexdey0 - _i); avirt = _i + bvirt; bround = bvirt - dexbey0; around = bexdey0 - avirt; bd[0] = around + bround; _j = (double) (bexdey1 + _i); bvirt = (double) (_j - bexdey1); avirt = _j - bvirt; bround = _i - bvirt; around = bexdey1 - avirt; _0 = around + bround; _i = (double) (_0 - dexbey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - dexbey1; around = _0 - avirt; bd[1] = around + bround; bd3 = (double) (_j + _i); bvirt = (double) (bd3 - _j); avirt = bd3 - bvirt; bround = _i - bvirt; around = _j - avirt; bd[2] = around + bround;
		bd[3] = bd3;

		temp8alen = scale_expansion_zeroelim(4, cd, bez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, -cez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, bc, dez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, aex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, -aex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, aey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, -aey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, aez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, -aez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		alen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, adet);

		temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, bex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, bex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, bey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, bey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, bez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, bez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		blen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, bdet);

		temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, cex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, -cex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, cey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, -cey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, cez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, -cez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		clen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cdet);

		temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b);
		temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
		temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, dex, temp48);
		xlen = scale_expansion_zeroelim(temp48len, temp48, dex, xdet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, dey, temp48);
		ylen = scale_expansion_zeroelim(temp48len, temp48, dey, ydet);
		temp48len = scale_expansion_zeroelim(temp24len, temp24, dez, temp48);
		zlen = scale_expansion_zeroelim(temp48len, temp48, dez, zdet);
		xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
		dlen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, ddet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, fin1);

		det = estimate(finlength, fin1);
		errbound = isperrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		bvirt = (double) (pa[0] - aex); avirt = aex + bvirt; bround = bvirt - pe[0]; around = pa[0] - avirt; aextail = around + bround;
		bvirt = (double) (pa[1] - aey); avirt = aey + bvirt; bround = bvirt - pe[1]; around = pa[1] - avirt; aeytail = around + bround;
		bvirt = (double) (pa[2] - aez); avirt = aez + bvirt; bround = bvirt - pe[2]; around = pa[2] - avirt; aeztail = around + bround;
		bvirt = (double) (pb[0] - bex); avirt = bex + bvirt; bround = bvirt - pe[0]; around = pb[0] - avirt; bextail = around + bround;
		bvirt = (double) (pb[1] - bey); avirt = bey + bvirt; bround = bvirt - pe[1]; around = pb[1] - avirt; beytail = around + bround;
		bvirt = (double) (pb[2] - bez); avirt = bez + bvirt; bround = bvirt - pe[2]; around = pb[2] - avirt; beztail = around + bround;
		bvirt = (double) (pc[0] - cex); avirt = cex + bvirt; bround = bvirt - pe[0]; around = pc[0] - avirt; cextail = around + bround;
		bvirt = (double) (pc[1] - cey); avirt = cey + bvirt; bround = bvirt - pe[1]; around = pc[1] - avirt; ceytail = around + bround;
		bvirt = (double) (pc[2] - cez); avirt = cez + bvirt; bround = bvirt - pe[2]; around = pc[2] - avirt; ceztail = around + bround;
		bvirt = (double) (pd[0] - dex); avirt = dex + bvirt; bround = bvirt - pe[0]; around = pd[0] - avirt; dextail = around + bround;
		bvirt = (double) (pd[1] - dey); avirt = dey + bvirt; bround = bvirt - pe[1]; around = pd[1] - avirt; deytail = around + bround;
		bvirt = (double) (pd[2] - dez); avirt = dez + bvirt; bround = bvirt - pe[2]; around = pd[2] - avirt; deztail = around + bround;
		if ((aextail == 0.0) && (aeytail == 0.0) && (aeztail == 0.0)
				&& (bextail == 0.0) && (beytail == 0.0) && (beztail == 0.0)
				&& (cextail == 0.0) && (ceytail == 0.0) && (ceztail == 0.0)
				&& (dextail == 0.0) && (deytail == 0.0) && (deztail == 0.0)) {
			return det;
		}

		errbound = isperrboundC * permanent + resulterrbound * ((det) >= 0.0 ? (det) : -(det));
		abeps = (aex * beytail + bey * aextail)
		- (aey * bextail + bex * aeytail);
		bceps = (bex * ceytail + cey * bextail)
		- (bey * cextail + cex * beytail);
		cdeps = (cex * deytail + dey * cextail)
		- (cey * dextail + dex * ceytail);
		daeps = (dex * aeytail + aey * dextail)
		- (dey * aextail + aex * deytail);
		aceps = (aex * ceytail + cey * aextail)
		- (aey * cextail + cex * aeytail);
		bdeps = (bex * deytail + dey * bextail)
		- (bey * dextail + dex * beytail);
		det += (((bex * bex + bey * bey + bez * bez)
				* ((cez * daeps + dez * aceps + aez * cdeps)
						+ (ceztail * da3 + deztail * ac3 + aeztail * cd3))
						+ (dex * dex + dey * dey + dez * dez)
						* ((aez * bceps - bez * aceps + cez * abeps)
								+ (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
								- ((aex * aex + aey * aey + aez * aez)
										* ((bez * cdeps - cez * bdeps + dez * bceps)
												+ (beztail * cd3 - ceztail * bd3 + deztail * bc3))
												+ (cex * cex + cey * cey + cez * cez)
												* ((dez * abeps + aez * bdeps + bez * daeps)
														+ (deztail * ab3 + aeztail * bd3 + beztail * da3))))
														+ 2.0 * (((bex * bextail + bey * beytail + bez * beztail)
																* (cez * da3 + dez * ac3 + aez * cd3)
																+ (dex * dextail + dey * deytail + dez * deztail)
																* (aez * bc3 - bez * ac3 + cez * ab3))
																- ((aex * aextail + aey * aeytail + aez * aeztail)
																		* (bez * cd3 - cez * bd3 + dez * bc3)
																		+ (cex * cextail + cey * ceytail + cez * ceztail)
																		* (dez * ab3 + aez * bd3 + bez * da3)));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		return insphereexact(pa, pb, pc, pd, pe);
	}
	
	double insphereexact(double[] pa,double[]  pb,double[]  pc,double[]  pd,double[]  pe)
	{
		double axby1, bxcy1, cxdy1, dxey1, exay1;
		double bxay1, cxby1, dxcy1, exdy1, axey1;
		double axcy1, bxdy1, cxey1, dxay1, exby1;
		double cxay1, dxby1, excy1, axdy1, bxey1;
		double axby0, bxcy0, cxdy0, dxey0, exay0;
		double bxay0, cxby0, dxcy0, exdy0, axey0;
		double axcy0, bxdy0, cxey0, dxay0, exby0;
		double cxay0, dxby0, excy0, axdy0, bxey0;
		double[] ab = new double[4], bc = new double[4], cd = new double[4], de = new double[4], ea = new double[4];
		double[] ac = new double[4], bd = new double[4], ce = new double[4], da = new double[4], eb = new double[4];
		double[] temp8a = new double[8], temp8b = new double[8], temp16 = new double[16];
		int temp8alen, temp8blen, temp16len;
		double[] abc = new double[24], bcd = new double[24], cde = new double[24], dea = new double[24], eab = new double[24];
		double[] abd = new double[24], bce = new double[24], cda = new double[24], deb = new double[24], eac = new double[24];
		int abclen, bcdlen, cdelen, dealen, eablen;
		int abdlen, bcelen, cdalen, deblen, eaclen;
		double[] temp48a = new double[48], temp48b = new double[48];
		int temp48alen, temp48blen;
		double[] abcd = new double[96], bcde = new double[96], cdea = new double[96], deab = new double[96], eabc = new double[96];
		int abcdlen, bcdelen, cdealen, deablen, eabclen;
		double[] temp192 = new double[192];
		double[] det384x = new double[384], det384y = new double[384], det384z = new double[384];
		int xlen, ylen, zlen;
		double[] detxy = new double[768];
		int xylen;
		double[] adet = new double[1152], bdet = new double[1152], cdet = new double[1152], ddet = new double[1152], edet = new double[1152];
		int alen, blen, clen, dlen, elen;
		double[] abdet = new double[2304], cddet = new double[2304], cdedet = new double[3456];
		int ablen, cdlen;
		double[] deter = new double[5760];
		int deterlen;
		int i;

		double bvirt;
		double avirt, bround, around;
		double c;
		double abig;
		double ahi, alo, bhi, blo;
		double err1, err2, err3;
		double _i, _j;
		double _0;

		axby1 = (double) (pa[0] * pb[1]); c = (double) (splitter * pa[0]); abig = (double) (c - pa[0]); ahi = c - abig; alo = pa[0] - ahi; c = (double) (splitter * pb[1]); abig = (double) (c - pb[1]); bhi = c - abig; blo = pb[1] - bhi; err1 = axby1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); axby0 = (alo * blo) - err3;
		bxay1 = (double) (pb[0] * pa[1]); c = (double) (splitter * pb[0]); abig = (double) (c - pb[0]); ahi = c - abig; alo = pb[0] - ahi; c = (double) (splitter * pa[1]); abig = (double) (c - pa[1]); bhi = c - abig; blo = pa[1] - bhi; err1 = bxay1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bxay0 = (alo * blo) - err3;
		_i = (double) (axby0 - bxay0); bvirt = (double) (axby0 - _i); avirt = _i + bvirt; bround = bvirt - bxay0; around = axby0 - avirt; ab[0] = around + bround; _j = (double) (axby1 + _i); bvirt = (double) (_j - axby1); avirt = _j - bvirt; bround = _i - bvirt; around = axby1 - avirt; _0 = around + bround; _i = (double) (_0 - bxay1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - bxay1; around = _0 - avirt; ab[1] = around + bround; ab[3] = (double) (_j + _i); bvirt = (double) (ab[3] - _j); avirt = ab[3] - bvirt; bround = _i - bvirt; around = _j - avirt; ab[2] = around + bround;

		bxcy1 = (double) (pb[0] * pc[1]); c = (double) (splitter * pb[0]); abig = (double) (c - pb[0]); ahi = c - abig; alo = pb[0] - ahi; c = (double) (splitter * pc[1]); abig = (double) (c - pc[1]); bhi = c - abig; blo = pc[1] - bhi; err1 = bxcy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bxcy0 = (alo * blo) - err3;
		cxby1 = (double) (pc[0] * pb[1]); c = (double) (splitter * pc[0]); abig = (double) (c - pc[0]); ahi = c - abig; alo = pc[0] - ahi; c = (double) (splitter * pb[1]); abig = (double) (c - pb[1]); bhi = c - abig; blo = pb[1] - bhi; err1 = cxby1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cxby0 = (alo * blo) - err3;
		_i = (double) (bxcy0 - cxby0); bvirt = (double) (bxcy0 - _i); avirt = _i + bvirt; bround = bvirt - cxby0; around = bxcy0 - avirt; bc[0] = around + bround; _j = (double) (bxcy1 + _i); bvirt = (double) (_j - bxcy1); avirt = _j - bvirt; bround = _i - bvirt; around = bxcy1 - avirt; _0 = around + bround; _i = (double) (_0 - cxby1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cxby1; around = _0 - avirt; bc[1] = around + bround; bc[3] = (double) (_j + _i); bvirt = (double) (bc[3] - _j); avirt = bc[3] - bvirt; bround = _i - bvirt; around = _j - avirt; bc[2] = around + bround;

		cxdy1 = (double) (pc[0] * pd[1]); c = (double) (splitter * pc[0]); abig = (double) (c - pc[0]); ahi = c - abig; alo = pc[0] - ahi; c = (double) (splitter * pd[1]); abig = (double) (c - pd[1]); bhi = c - abig; blo = pd[1] - bhi; err1 = cxdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cxdy0 = (alo * blo) - err3;
		dxcy1 = (double) (pd[0] * pc[1]); c = (double) (splitter * pd[0]); abig = (double) (c - pd[0]); ahi = c - abig; alo = pd[0] - ahi; c = (double) (splitter * pc[1]); abig = (double) (c - pc[1]); bhi = c - abig; blo = pc[1] - bhi; err1 = dxcy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dxcy0 = (alo * blo) - err3;
		_i = (double) (cxdy0 - dxcy0); bvirt = (double) (cxdy0 - _i); avirt = _i + bvirt; bround = bvirt - dxcy0; around = cxdy0 - avirt; cd[0] = around + bround; _j = (double) (cxdy1 + _i); bvirt = (double) (_j - cxdy1); avirt = _j - bvirt; bround = _i - bvirt; around = cxdy1 - avirt; _0 = around + bround; _i = (double) (_0 - dxcy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - dxcy1; around = _0 - avirt; cd[1] = around + bround; cd[3] = (double) (_j + _i); bvirt = (double) (cd[3] - _j); avirt = cd[3] - bvirt; bround = _i - bvirt; around = _j - avirt; cd[2] = around + bround;

		dxey1 = (double) (pd[0] * pe[1]); c = (double) (splitter * pd[0]); abig = (double) (c - pd[0]); ahi = c - abig; alo = pd[0] - ahi; c = (double) (splitter * pe[1]); abig = (double) (c - pe[1]); bhi = c - abig; blo = pe[1] - bhi; err1 = dxey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dxey0 = (alo * blo) - err3;
		exdy1 = (double) (pe[0] * pd[1]); c = (double) (splitter * pe[0]); abig = (double) (c - pe[0]); ahi = c - abig; alo = pe[0] - ahi; c = (double) (splitter * pd[1]); abig = (double) (c - pd[1]); bhi = c - abig; blo = pd[1] - bhi; err1 = exdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); exdy0 = (alo * blo) - err3;
		_i = (double) (dxey0 - exdy0); bvirt = (double) (dxey0 - _i); avirt = _i + bvirt; bround = bvirt - exdy0; around = dxey0 - avirt; de[0] = around + bround; _j = (double) (dxey1 + _i); bvirt = (double) (_j - dxey1); avirt = _j - bvirt; bround = _i - bvirt; around = dxey1 - avirt; _0 = around + bround; _i = (double) (_0 - exdy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - exdy1; around = _0 - avirt; de[1] = around + bround; de[3] = (double) (_j + _i); bvirt = (double) (de[3] - _j); avirt = de[3] - bvirt; bround = _i - bvirt; around = _j - avirt; de[2] = around + bround;

		exay1 = (double) (pe[0] * pa[1]); c = (double) (splitter * pe[0]); abig = (double) (c - pe[0]); ahi = c - abig; alo = pe[0] - ahi; c = (double) (splitter * pa[1]); abig = (double) (c - pa[1]); bhi = c - abig; blo = pa[1] - bhi; err1 = exay1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); exay0 = (alo * blo) - err3;
		axey1 = (double) (pa[0] * pe[1]); c = (double) (splitter * pa[0]); abig = (double) (c - pa[0]); ahi = c - abig; alo = pa[0] - ahi; c = (double) (splitter * pe[1]); abig = (double) (c - pe[1]); bhi = c - abig; blo = pe[1] - bhi; err1 = axey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); axey0 = (alo * blo) - err3;
		_i = (double) (exay0 - axey0); bvirt = (double) (exay0 - _i); avirt = _i + bvirt; bround = bvirt - axey0; around = exay0 - avirt; ea[0] = around + bround; _j = (double) (exay1 + _i); bvirt = (double) (_j - exay1); avirt = _j - bvirt; bround = _i - bvirt; around = exay1 - avirt; _0 = around + bround; _i = (double) (_0 - axey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - axey1; around = _0 - avirt; ea[1] = around + bround; ea[3] = (double) (_j + _i); bvirt = (double) (ea[3] - _j); avirt = ea[3] - bvirt; bround = _i - bvirt; around = _j - avirt; ea[2] = around + bround;

		axcy1 = (double) (pa[0] * pc[1]); c = (double) (splitter * pa[0]); abig = (double) (c - pa[0]); ahi = c - abig; alo = pa[0] - ahi; c = (double) (splitter * pc[1]); abig = (double) (c - pc[1]); bhi = c - abig; blo = pc[1] - bhi; err1 = axcy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); axcy0 = (alo * blo) - err3;
		cxay1 = (double) (pc[0] * pa[1]); c = (double) (splitter * pc[0]); abig = (double) (c - pc[0]); ahi = c - abig; alo = pc[0] - ahi; c = (double) (splitter * pa[1]); abig = (double) (c - pa[1]); bhi = c - abig; blo = pa[1] - bhi; err1 = cxay1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cxay0 = (alo * blo) - err3;
		_i = (double) (axcy0 - cxay0); bvirt = (double) (axcy0 - _i); avirt = _i + bvirt; bround = bvirt - cxay0; around = axcy0 - avirt; ac[0] = around + bround; _j = (double) (axcy1 + _i); bvirt = (double) (_j - axcy1); avirt = _j - bvirt; bround = _i - bvirt; around = axcy1 - avirt; _0 = around + bround; _i = (double) (_0 - cxay1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cxay1; around = _0 - avirt; ac[1] = around + bround; ac[3] = (double) (_j + _i); bvirt = (double) (ac[3] - _j); avirt = ac[3] - bvirt; bround = _i - bvirt; around = _j - avirt; ac[2] = around + bround;

		bxdy1 = (double) (pb[0] * pd[1]); c = (double) (splitter * pb[0]); abig = (double) (c - pb[0]); ahi = c - abig; alo = pb[0] - ahi; c = (double) (splitter * pd[1]); abig = (double) (c - pd[1]); bhi = c - abig; blo = pd[1] - bhi; err1 = bxdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bxdy0 = (alo * blo) - err3;
		dxby1 = (double) (pd[0] * pb[1]); c = (double) (splitter * pd[0]); abig = (double) (c - pd[0]); ahi = c - abig; alo = pd[0] - ahi; c = (double) (splitter * pb[1]); abig = (double) (c - pb[1]); bhi = c - abig; blo = pb[1] - bhi; err1 = dxby1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dxby0 = (alo * blo) - err3;
		_i = (double) (bxdy0 - dxby0); bvirt = (double) (bxdy0 - _i); avirt = _i + bvirt; bround = bvirt - dxby0; around = bxdy0 - avirt; bd[0] = around + bround; _j = (double) (bxdy1 + _i); bvirt = (double) (_j - bxdy1); avirt = _j - bvirt; bround = _i - bvirt; around = bxdy1 - avirt; _0 = around + bround; _i = (double) (_0 - dxby1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - dxby1; around = _0 - avirt; bd[1] = around + bround; bd[3] = (double) (_j + _i); bvirt = (double) (bd[3] - _j); avirt = bd[3] - bvirt; bround = _i - bvirt; around = _j - avirt; bd[2] = around + bround;

		cxey1 = (double) (pc[0] * pe[1]); c = (double) (splitter * pc[0]); abig = (double) (c - pc[0]); ahi = c - abig; alo = pc[0] - ahi; c = (double) (splitter * pe[1]); abig = (double) (c - pe[1]); bhi = c - abig; blo = pe[1] - bhi; err1 = cxey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cxey0 = (alo * blo) - err3;
		excy1 = (double) (pe[0] * pc[1]); c = (double) (splitter * pe[0]); abig = (double) (c - pe[0]); ahi = c - abig; alo = pe[0] - ahi; c = (double) (splitter * pc[1]); abig = (double) (c - pc[1]); bhi = c - abig; blo = pc[1] - bhi; err1 = excy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); excy0 = (alo * blo) - err3;
		_i = (double) (cxey0 - excy0); bvirt = (double) (cxey0 - _i); avirt = _i + bvirt; bround = bvirt - excy0; around = cxey0 - avirt; ce[0] = around + bround; _j = (double) (cxey1 + _i); bvirt = (double) (_j - cxey1); avirt = _j - bvirt; bround = _i - bvirt; around = cxey1 - avirt; _0 = around + bround; _i = (double) (_0 - excy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - excy1; around = _0 - avirt; ce[1] = around + bround; ce[3] = (double) (_j + _i); bvirt = (double) (ce[3] - _j); avirt = ce[3] - bvirt; bround = _i - bvirt; around = _j - avirt; ce[2] = around + bround;

		dxay1 = (double) (pd[0] * pa[1]); c = (double) (splitter * pd[0]); abig = (double) (c - pd[0]); ahi = c - abig; alo = pd[0] - ahi; c = (double) (splitter * pa[1]); abig = (double) (c - pa[1]); bhi = c - abig; blo = pa[1] - bhi; err1 = dxay1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); dxay0 = (alo * blo) - err3;
		axdy1 = (double) (pa[0] * pd[1]); c = (double) (splitter * pa[0]); abig = (double) (c - pa[0]); ahi = c - abig; alo = pa[0] - ahi; c = (double) (splitter * pd[1]); abig = (double) (c - pd[1]); bhi = c - abig; blo = pd[1] - bhi; err1 = axdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); axdy0 = (alo * blo) - err3;
		_i = (double) (dxay0 - axdy0); bvirt = (double) (dxay0 - _i); avirt = _i + bvirt; bround = bvirt - axdy0; around = dxay0 - avirt; da[0] = around + bround; _j = (double) (dxay1 + _i); bvirt = (double) (_j - dxay1); avirt = _j - bvirt; bround = _i - bvirt; around = dxay1 - avirt; _0 = around + bround; _i = (double) (_0 - axdy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - axdy1; around = _0 - avirt; da[1] = around + bround; da[3] = (double) (_j + _i); bvirt = (double) (da[3] - _j); avirt = da[3] - bvirt; bround = _i - bvirt; around = _j - avirt; da[2] = around + bround;

		exby1 = (double) (pe[0] * pb[1]); c = (double) (splitter * pe[0]); abig = (double) (c - pe[0]); ahi = c - abig; alo = pe[0] - ahi; c = (double) (splitter * pb[1]); abig = (double) (c - pb[1]); bhi = c - abig; blo = pb[1] - bhi; err1 = exby1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); exby0 = (alo * blo) - err3;
		bxey1 = (double) (pb[0] * pe[1]); c = (double) (splitter * pb[0]); abig = (double) (c - pb[0]); ahi = c - abig; alo = pb[0] - ahi; c = (double) (splitter * pe[1]); abig = (double) (c - pe[1]); bhi = c - abig; blo = pe[1] - bhi; err1 = bxey1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bxey0 = (alo * blo) - err3;
		_i = (double) (exby0 - bxey0); bvirt = (double) (exby0 - _i); avirt = _i + bvirt; bround = bvirt - bxey0; around = exby0 - avirt; eb[0] = around + bround; _j = (double) (exby1 + _i); bvirt = (double) (_j - exby1); avirt = _j - bvirt; bround = _i - bvirt; around = exby1 - avirt; _0 = around + bround; _i = (double) (_0 - bxey1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - bxey1; around = _0 - avirt; eb[1] = around + bround; eb[3] = (double) (_j + _i); bvirt = (double) (eb[3] - _j); avirt = eb[3] - bvirt; bround = _i - bvirt; around = _j - avirt; eb[2] = around + bround;

		temp8alen = scale_expansion_zeroelim(4, bc, pa[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, -pb[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, ab, pc[2], temp8a);
		abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				abc);

		temp8alen = scale_expansion_zeroelim(4, cd, pb[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, -pc[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, bc, pd[2], temp8a);
		bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				bcd);

		temp8alen = scale_expansion_zeroelim(4, de, pc[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ce, -pd[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, cd, pe[2], temp8a);
		cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				cde);

		temp8alen = scale_expansion_zeroelim(4, ea, pd[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, da, -pe[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, de, pa[2], temp8a);
		dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				dea);

		temp8alen = scale_expansion_zeroelim(4, ab, pe[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, eb, -pa[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, ea, pb[2], temp8a);
		eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				eab);

		temp8alen = scale_expansion_zeroelim(4, bd, pa[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, da, pb[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, ab, pd[2], temp8a);
		abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				abd);

		temp8alen = scale_expansion_zeroelim(4, ce, pb[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, eb, pc[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, bc, pe[2], temp8a);
		bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				bce);

		temp8alen = scale_expansion_zeroelim(4, da, pc[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ac, pd[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, cd, pa[2], temp8a);
		cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				cda);

		temp8alen = scale_expansion_zeroelim(4, eb, pd[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, bd, pe[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, de, pb[2], temp8a);
		deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				deb);

		temp8alen = scale_expansion_zeroelim(4, ac, pe[2], temp8a);
		temp8blen = scale_expansion_zeroelim(4, ce, pa[2], temp8b);
		temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
		temp8alen = scale_expansion_zeroelim(4, ea, pc[2], temp8a);
		eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				eac);

		temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, bcde);
		xlen = scale_expansion_zeroelim(bcdelen, bcde, pa[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pa[0], det384x);
		ylen = scale_expansion_zeroelim(bcdelen, bcde, pa[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pa[1], det384y);
		zlen = scale_expansion_zeroelim(bcdelen, bcde, pa[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pa[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, adet);

		temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, cdea);
		xlen = scale_expansion_zeroelim(cdealen, cdea, pb[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pb[0], det384x);
		ylen = scale_expansion_zeroelim(cdealen, cdea, pb[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pb[1], det384y);
		zlen = scale_expansion_zeroelim(cdealen, cdea, pb[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pb[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, bdet);

		temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, deab);
		xlen = scale_expansion_zeroelim(deablen, deab, pc[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pc[0], det384x);
		ylen = scale_expansion_zeroelim(deablen, deab, pc[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pc[1], det384y);
		zlen = scale_expansion_zeroelim(deablen, deab, pc[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pc[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cdet);

		temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, eabc);
		xlen = scale_expansion_zeroelim(eabclen, eabc, pd[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pd[0], det384x);
		ylen = scale_expansion_zeroelim(eabclen, eabc, pd[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pd[1], det384y);
		zlen = scale_expansion_zeroelim(eabclen, eabc, pd[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pd[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, ddet);

		temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
		temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
		for (i = 0; i < temp48blen; i++) {
			temp48b[i] = -temp48b[i];
		}
		abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, abcd);
		xlen = scale_expansion_zeroelim(abcdlen, abcd, pe[0], temp192);
		xlen = scale_expansion_zeroelim(xlen, temp192, pe[0], det384x);
		ylen = scale_expansion_zeroelim(abcdlen, abcd, pe[1], temp192);
		ylen = scale_expansion_zeroelim(ylen, temp192, pe[1], det384y);
		zlen = scale_expansion_zeroelim(abcdlen, abcd, pe[2], temp192);
		zlen = scale_expansion_zeroelim(zlen, temp192, pe[2], det384z);
		xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
		elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, edet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
		cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
		deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

		return deter[deterlen - 1];
	}
	
	
	

	double orient(double[] pa, double[] pb, double[] pc, double[] pd)	{
		double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
		double det;
		double permanent, errbound;

		adx = pa[0] - pd[0];
		bdx = pb[0] - pd[0];
		cdx = pc[0] - pd[0];
		ady = pa[1] - pd[1];
		bdy = pb[1] - pd[1];
		cdy = pc[1] - pd[1];
		adz = pa[2] - pd[2];
		bdz = pb[2] - pd[2];
		cdz = pc[2] - pd[2];

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;

		det = adz * (bdxcdy - cdxbdy)
				+ bdz * (cdxady - adxcdy)
				+ cdz * (adxbdy - bdxady);

		permanent = (((bdxcdy) >= 0.0 ? (bdxcdy) : -(bdxcdy)) + ((cdxbdy) >= 0.0 ? (cdxbdy) : -(cdxbdy))) * ((adz) >= 0.0 ? (adz) : -(adz))
				+ (((cdxady) >= 0.0 ? (cdxady) : -(cdxady)) + ((adxcdy) >= 0.0 ? (adxcdy) : -(adxcdy))) * ((bdz) >= 0.0 ? (bdz) : -(bdz))
				+ (((adxbdy) >= 0.0 ? (adxbdy) : -(adxbdy)) + ((bdxady) >= 0.0 ? (bdxady) : -(bdxady))) * ((cdz) >= 0.0 ? (cdz) : -(cdz));
		errbound = o3derrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}
		return orient3dadapt(pa, pb, pc, pd, permanent);
	}

	double orient3dadapt(double[] pa,double[] pb,double[] pc,double[] pd,double permanent)
	{
		double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		double det, errbound;

		double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
		double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
		double[] bc = new double[4], ca = new double[4], ab = new double[4];
		double bc3, ca3, ab3;
		double[] adet = new double[8], bdet = new double[8], cdet = new double[8];
		int alen, blen, clen;
		double[] abdet = new double[16];
		int ablen;
		double[] finnow, finother, finswap;
		double[] fin1 = new double[192], fin2 = new double[192];
		int finlength;

		double adxtail, bdxtail, cdxtail;
		double adytail, bdytail, cdytail;
		double adztail, bdztail, cdztail;
		double at_blarge, at_clarge;
		double bt_clarge, bt_alarge;
		double ct_alarge, ct_blarge;
		double[] at_b = new double[4], at_c = new double[4], bt_c = new double[4], bt_a = new double[4], ct_a = new double[4], ct_b = new double[4];
		int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
		double bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
		double adxt_cdy1, adxt_bdy1, bdxt_ady1;
		double bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
		double adxt_cdy0, adxt_bdy0, bdxt_ady0;
		double bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
		double adyt_cdx1, adyt_bdx1, bdyt_adx1;
		double bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
		double adyt_cdx0, adyt_bdx0, bdyt_adx0;
		double[] bct = new double[8], cat = new double[8], abt = new double[8];
		int bctlen, catlen, abtlen;
		double bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
		double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
		double bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
		double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
		double[] u = new double[4], v = new double[12], w = new double[16];
		double u3;
		int vlength, wlength;
		double negate;

		double bvirt;
		double avirt, bround, around;
		double c;
		double abig;
		double ahi, alo, bhi, blo;
		double err1, err2, err3;
		double _i, _j, _k;
		double _0;

		adx = (double) (pa[0] - pd[0]);
		bdx = (double) (pb[0] - pd[0]);
		cdx = (double) (pc[0] - pd[0]);
		ady = (double) (pa[1] - pd[1]);
		bdy = (double) (pb[1] - pd[1]);
		cdy = (double) (pc[1] - pd[1]);
		adz = (double) (pa[2] - pd[2]);
		bdz = (double) (pb[2] - pd[2]);
		cdz = (double) (pc[2] - pd[2]);

		bdxcdy1 = (double) (bdx * cdy); c = (double) (splitter * bdx); abig = (double) (c - bdx); ahi = c - abig; alo = bdx - ahi; c = (double) (splitter * cdy); abig = (double) (c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = bdxcdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxcdy0 = (alo * blo) - err3;
		cdxbdy1 = (double) (cdx * bdy); c = (double) (splitter * cdx); abig = (double) (c - cdx); ahi = c - abig; alo = cdx - ahi; c = (double) (splitter * bdy); abig = (double) (c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = cdxbdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxbdy0 = (alo * blo) - err3;
		_i = (double) (bdxcdy0 - cdxbdy0); bvirt = (double) (bdxcdy0 - _i); avirt = _i + bvirt; bround = bvirt - cdxbdy0; around = bdxcdy0 - avirt; bc[0] = around + bround; _j = (double) (bdxcdy1 + _i); bvirt = (double) (_j - bdxcdy1); avirt = _j - bvirt; bround = _i - bvirt; around = bdxcdy1 - avirt; _0 = around + bround; _i = (double) (_0 - cdxbdy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cdxbdy1; around = _0 - avirt; bc[1] = around + bround; bc3 = (double) (_j + _i); bvirt = (double) (bc3 - _j); avirt = bc3 - bvirt; bround = _i - bvirt; around = _j - avirt; bc[2] = around + bround;
		bc[3] = bc3;
		alen = scale_expansion_zeroelim(4, bc, adz, adet);

		cdxady1 = (double) (cdx * ady); c = (double) (splitter * cdx); abig = (double) (c - cdx); ahi = c - abig; alo = cdx - ahi; c = (double) (splitter * ady); abig = (double) (c - ady); bhi = c - abig; blo = ady - bhi; err1 = cdxady1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxady0 = (alo * blo) - err3;
		adxcdy1 = (double) (adx * cdy); c = (double) (splitter * adx); abig = (double) (c - adx); ahi = c - abig; alo = adx - ahi; c = (double) (splitter * cdy); abig = (double) (c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = adxcdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxcdy0 = (alo * blo) - err3;
		_i = (double) (cdxady0 - adxcdy0); bvirt = (double) (cdxady0 - _i); avirt = _i + bvirt; bround = bvirt - adxcdy0; around = cdxady0 - avirt; ca[0] = around + bround; _j = (double) (cdxady1 + _i); bvirt = (double) (_j - cdxady1); avirt = _j - bvirt; bround = _i - bvirt; around = cdxady1 - avirt; _0 = around + bround; _i = (double) (_0 - adxcdy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - adxcdy1; around = _0 - avirt; ca[1] = around + bround; ca3 = (double) (_j + _i); bvirt = (double) (ca3 - _j); avirt = ca3 - bvirt; bround = _i - bvirt; around = _j - avirt; ca[2] = around + bround;
		ca[3] = ca3;
		blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

		adxbdy1 = (double) (adx * bdy); c = (double) (splitter * adx); abig = (double) (c - adx); ahi = c - abig; alo = adx - ahi; c = (double) (splitter * bdy); abig = (double) (c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = adxbdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxbdy0 = (alo * blo) - err3;
		bdxady1 = (double) (bdx * ady); c = (double) (splitter * bdx); abig = (double) (c - bdx); ahi = c - abig; alo = bdx - ahi; c = (double) (splitter * ady); abig = (double) (c - ady); bhi = c - abig; blo = ady - bhi; err1 = bdxady1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxady0 = (alo * blo) - err3;
		_i = (double) (adxbdy0 - bdxady0); bvirt = (double) (adxbdy0 - _i); avirt = _i + bvirt; bround = bvirt - bdxady0; around = adxbdy0 - avirt; ab[0] = around + bround; _j = (double) (adxbdy1 + _i); bvirt = (double) (_j - adxbdy1); avirt = _j - bvirt; bround = _i - bvirt; around = adxbdy1 - avirt; _0 = around + bround; _i = (double) (_0 - bdxady1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - bdxady1; around = _0 - avirt; ab[1] = around + bround; ab3 = (double) (_j + _i); bvirt = (double) (ab3 - _j); avirt = ab3 - bvirt; bround = _i - bvirt; around = _j - avirt; ab[2] = around + bround;
		ab[3] = ab3;
		clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

		det = estimate(finlength, fin1);
		errbound = o3derrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		bvirt = (double) (pa[0] - adx); avirt = adx + bvirt; bround = bvirt - pd[0]; around = pa[0] - avirt; adxtail = around + bround;
		bvirt = (double) (pb[0] - bdx); avirt = bdx + bvirt; bround = bvirt - pd[0]; around = pb[0] - avirt; bdxtail = around + bround;
		bvirt = (double) (pc[0] - cdx); avirt = cdx + bvirt; bround = bvirt - pd[0]; around = pc[0] - avirt; cdxtail = around + bround;
		bvirt = (double) (pa[1] - ady); avirt = ady + bvirt; bround = bvirt - pd[1]; around = pa[1] - avirt; adytail = around + bround;
		bvirt = (double) (pb[1] - bdy); avirt = bdy + bvirt; bround = bvirt - pd[1]; around = pb[1] - avirt; bdytail = around + bround;
		bvirt = (double) (pc[1] - cdy); avirt = cdy + bvirt; bround = bvirt - pd[1]; around = pc[1] - avirt; cdytail = around + bround;
		bvirt = (double) (pa[2] - adz); avirt = adz + bvirt; bround = bvirt - pd[2]; around = pa[2] - avirt; adztail = around + bround;
		bvirt = (double) (pb[2] - bdz); avirt = bdz + bvirt; bround = bvirt - pd[2]; around = pb[2] - avirt; bdztail = around + bround;
		bvirt = (double) (pc[2] - cdz); avirt = cdz + bvirt; bround = bvirt - pd[2]; around = pc[2] - avirt; cdztail = around + bround;

		if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
				&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
				&& (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
			return det;
		}

		errbound = o3derrboundC * permanent + resulterrbound * ((det) >= 0.0 ? (det) : -(det));
		det += (adz * ((bdx * cdytail + cdy * bdxtail)
				- (bdy * cdxtail + cdx * bdytail))
				+ adztail * (bdx * cdy - bdy * cdx))
				+ (bdz * ((cdx * adytail + ady * cdxtail)
						- (cdy * adxtail + adx * cdytail))
						+ bdztail * (cdx * ady - cdy * adx))
						+ (cdz * ((adx * bdytail + bdy * adxtail)
								- (ady * bdxtail + bdx * adytail))
								+ cdztail * (adx * bdy - ady * bdx));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		finnow = fin1;
		finother = fin2;

		if (adxtail == 0.0) {
			if (adytail == 0.0) {
				at_b[0] = 0.0;
				at_blen = 1;
				at_c[0] = 0.0;
				at_clen = 1;
			} else {
				negate = -adytail;
				at_blarge = (double) (negate * bdx); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * bdx); abig = (double) (c - bdx); bhi = c - abig; blo = bdx - bhi; err1 = at_blarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); at_b[0] = (alo * blo) - err3;
				at_b[1] = at_blarge;
				at_blen = 2;
				at_clarge = (double) (adytail * cdx); c = (double) (splitter * adytail); abig = (double) (c - adytail); ahi = c - abig; alo = adytail - ahi; c = (double) (splitter * cdx); abig = (double) (c - cdx); bhi = c - abig; blo = cdx - bhi; err1 = at_clarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); at_c[0] = (alo * blo) - err3;
				at_c[1] = at_clarge;
				at_clen = 2;
			}
		} else {
			if (adytail == 0.0) {
				at_blarge = (double) (adxtail * bdy); c = (double) (splitter * adxtail); abig = (double) (c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (double) (splitter * bdy); abig = (double) (c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = at_blarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); at_b[0] = (alo * blo) - err3;
				at_b[1] = at_blarge;
				at_blen = 2;
				negate = -adxtail;
				at_clarge = (double) (negate * cdy); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * cdy); abig = (double) (c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = at_clarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); at_c[0] = (alo * blo) - err3;
				at_c[1] = at_clarge;
				at_clen = 2;
			} else {
				adxt_bdy1 = (double) (adxtail * bdy); c = (double) (splitter * adxtail); abig = (double) (c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (double) (splitter * bdy); abig = (double) (c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = adxt_bdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxt_bdy0 = (alo * blo) - err3;
				adyt_bdx1 = (double) (adytail * bdx); c = (double) (splitter * adytail); abig = (double) (c - adytail); ahi = c - abig; alo = adytail - ahi; c = (double) (splitter * bdx); abig = (double) (c - bdx); bhi = c - abig; blo = bdx - bhi; err1 = adyt_bdx1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adyt_bdx0 = (alo * blo) - err3;
				_i = (double) (adxt_bdy0 - adyt_bdx0); bvirt = (double) (adxt_bdy0 - _i); avirt = _i + bvirt; bround = bvirt - adyt_bdx0; around = adxt_bdy0 - avirt; at_b[0] = around + bround; _j = (double) (adxt_bdy1 + _i); bvirt = (double) (_j - adxt_bdy1); avirt = _j - bvirt; bround = _i - bvirt; around = adxt_bdy1 - avirt; _0 = around + bround; _i = (double) (_0 - adyt_bdx1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - adyt_bdx1; around = _0 - avirt; at_b[1] = around + bround; at_blarge = (double) (_j + _i); bvirt = (double) (at_blarge - _j); avirt = at_blarge - bvirt; bround = _i - bvirt; around = _j - avirt; at_b[2] = around + bround;

				at_b[3] = at_blarge;
				at_blen = 4;
				adyt_cdx1 = (double) (adytail * cdx); c = (double) (splitter * adytail); abig = (double) (c - adytail); ahi = c - abig; alo = adytail - ahi; c = (double) (splitter * cdx); abig = (double) (c - cdx); bhi = c - abig; blo = cdx - bhi; err1 = adyt_cdx1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adyt_cdx0 = (alo * blo) - err3;
				adxt_cdy1 = (double) (adxtail * cdy); c = (double) (splitter * adxtail); abig = (double) (c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (double) (splitter * cdy); abig = (double) (c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = adxt_cdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxt_cdy0 = (alo * blo) - err3;
				_i = (double) (adyt_cdx0 - adxt_cdy0); bvirt = (double) (adyt_cdx0 - _i); avirt = _i + bvirt; bround = bvirt - adxt_cdy0; around = adyt_cdx0 - avirt; at_c[0] = around + bround; _j = (double) (adyt_cdx1 + _i); bvirt = (double) (_j - adyt_cdx1); avirt = _j - bvirt; bround = _i - bvirt; around = adyt_cdx1 - avirt; _0 = around + bround; _i = (double) (_0 - adxt_cdy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - adxt_cdy1; around = _0 - avirt; at_c[1] = around + bround; at_clarge = (double) (_j + _i); bvirt = (double) (at_clarge - _j); avirt = at_clarge - bvirt; bround = _i - bvirt; around = _j - avirt; at_c[2] = around + bround;

				at_c[3] = at_clarge;
				at_clen = 4;
			}
		}
		if (bdxtail == 0.0) {
			if (bdytail == 0.0) {
				bt_c[0] = 0.0;
				bt_clen = 1;
				bt_a[0] = 0.0;
				bt_alen = 1;
			} else {
				negate = -bdytail;
				bt_clarge = (double) (negate * cdx); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * cdx); abig = (double) (c - cdx); bhi = c - abig; blo = cdx - bhi; err1 = bt_clarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bt_c[0] = (alo * blo) - err3;
				bt_c[1] = bt_clarge;
				bt_clen = 2;
				bt_alarge = (double) (bdytail * adx); c = (double) (splitter * bdytail); abig = (double) (c - bdytail); ahi = c - abig; alo = bdytail - ahi; c = (double) (splitter * adx); abig = (double) (c - adx); bhi = c - abig; blo = adx - bhi; err1 = bt_alarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bt_a[0] = (alo * blo) - err3;
				bt_a[1] = bt_alarge;
				bt_alen = 2;
			}
		} else {
			if (bdytail == 0.0) {
				bt_clarge = (double) (bdxtail * cdy); c = (double) (splitter * bdxtail); abig = (double) (c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (double) (splitter * cdy); abig = (double) (c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = bt_clarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bt_c[0] = (alo * blo) - err3;
				bt_c[1] = bt_clarge;
				bt_clen = 2;
				negate = -bdxtail;
				bt_alarge = (double) (negate * ady); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * ady); abig = (double) (c - ady); bhi = c - abig; blo = ady - bhi; err1 = bt_alarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bt_a[0] = (alo * blo) - err3;
				bt_a[1] = bt_alarge;
				bt_alen = 2;
			} else {
				bdxt_cdy1 = (double) (bdxtail * cdy); c = (double) (splitter * bdxtail); abig = (double) (c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (double) (splitter * cdy); abig = (double) (c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = bdxt_cdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxt_cdy0 = (alo * blo) - err3;
				bdyt_cdx1 = (double) (bdytail * cdx); c = (double) (splitter * bdytail); abig = (double) (c - bdytail); ahi = c - abig; alo = bdytail - ahi; c = (double) (splitter * cdx); abig = (double) (c - cdx); bhi = c - abig; blo = cdx - bhi; err1 = bdyt_cdx1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdyt_cdx0 = (alo * blo) - err3;
				_i = (double) (bdxt_cdy0 - bdyt_cdx0); bvirt = (double) (bdxt_cdy0 - _i); avirt = _i + bvirt; bround = bvirt - bdyt_cdx0; around = bdxt_cdy0 - avirt; bt_c[0] = around + bround; _j = (double) (bdxt_cdy1 + _i); bvirt = (double) (_j - bdxt_cdy1); avirt = _j - bvirt; bround = _i - bvirt; around = bdxt_cdy1 - avirt; _0 = around + bround; _i = (double) (_0 - bdyt_cdx1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - bdyt_cdx1; around = _0 - avirt; bt_c[1] = around + bround; bt_clarge = (double) (_j + _i); bvirt = (double) (bt_clarge - _j); avirt = bt_clarge - bvirt; bround = _i - bvirt; around = _j - avirt; bt_c[2] = around + bround;

				bt_c[3] = bt_clarge;
				bt_clen = 4;
				bdyt_adx1 = (double) (bdytail * adx); c = (double) (splitter * bdytail); abig = (double) (c - bdytail); ahi = c - abig; alo = bdytail - ahi; c = (double) (splitter * adx); abig = (double) (c - adx); bhi = c - abig; blo = adx - bhi; err1 = bdyt_adx1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdyt_adx0 = (alo * blo) - err3;
				bdxt_ady1 = (double) (bdxtail * ady); c = (double) (splitter * bdxtail); abig = (double) (c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (double) (splitter * ady); abig = (double) (c - ady); bhi = c - abig; blo = ady - bhi; err1 = bdxt_ady1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxt_ady0 = (alo * blo) - err3;
				_i = (double) (bdyt_adx0 - bdxt_ady0); bvirt = (double) (bdyt_adx0 - _i); avirt = _i + bvirt; bround = bvirt - bdxt_ady0; around = bdyt_adx0 - avirt; bt_a[0] = around + bround; _j = (double) (bdyt_adx1 + _i); bvirt = (double) (_j - bdyt_adx1); avirt = _j - bvirt; bround = _i - bvirt; around = bdyt_adx1 - avirt; _0 = around + bround; _i = (double) (_0 - bdxt_ady1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - bdxt_ady1; around = _0 - avirt; bt_a[1] = around + bround; bt_alarge = (double) (_j + _i); bvirt = (double) (bt_alarge - _j); avirt = bt_alarge - bvirt; bround = _i - bvirt; around = _j - avirt; bt_a[2] = around + bround;

				bt_a[3] = bt_alarge;
				bt_alen = 4;
			}
		}
		if (cdxtail == 0.0) {
			if (cdytail == 0.0) {
				ct_a[0] = 0.0;
				ct_alen = 1;
				ct_b[0] = 0.0;
				ct_blen = 1;
			} else {
				negate = -cdytail;
				ct_alarge = (double) (negate * adx); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * adx); abig = (double) (c - adx); bhi = c - abig; blo = adx - bhi; err1 = ct_alarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ct_a[0] = (alo * blo) - err3;
				ct_a[1] = ct_alarge;
				ct_alen = 2;
				ct_blarge = (double) (cdytail * bdx); c = (double) (splitter * cdytail); abig = (double) (c - cdytail); ahi = c - abig; alo = cdytail - ahi; c = (double) (splitter * bdx); abig = (double) (c - bdx); bhi = c - abig; blo = bdx - bhi; err1 = ct_blarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ct_b[0] = (alo * blo) - err3;
				ct_b[1] = ct_blarge;
				ct_blen = 2;
			}
		} else {
			if (cdytail == 0.0) {
				ct_alarge = (double) (cdxtail * ady); c = (double) (splitter * cdxtail); abig = (double) (c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (double) (splitter * ady); abig = (double) (c - ady); bhi = c - abig; blo = ady - bhi; err1 = ct_alarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ct_a[0] = (alo * blo) - err3;
				ct_a[1] = ct_alarge;
				ct_alen = 2;
				negate = -cdxtail;
				ct_blarge = (double) (negate * bdy); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * bdy); abig = (double) (c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = ct_blarge - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ct_b[0] = (alo * blo) - err3;
				ct_b[1] = ct_blarge;
				ct_blen = 2;
			} else {
				cdxt_ady1 = (double) (cdxtail * ady); c = (double) (splitter * cdxtail); abig = (double) (c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (double) (splitter * ady); abig = (double) (c - ady); bhi = c - abig; blo = ady - bhi; err1 = cdxt_ady1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxt_ady0 = (alo * blo) - err3;
				cdyt_adx1 = (double) (cdytail * adx); c = (double) (splitter * cdytail); abig = (double) (c - cdytail); ahi = c - abig; alo = cdytail - ahi; c = (double) (splitter * adx); abig = (double) (c - adx); bhi = c - abig; blo = adx - bhi; err1 = cdyt_adx1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdyt_adx0 = (alo * blo) - err3;
				_i = (double) (cdxt_ady0 - cdyt_adx0); bvirt = (double) (cdxt_ady0 - _i); avirt = _i + bvirt; bround = bvirt - cdyt_adx0; around = cdxt_ady0 - avirt; ct_a[0] = around + bround; _j = (double) (cdxt_ady1 + _i); bvirt = (double) (_j - cdxt_ady1); avirt = _j - bvirt; bround = _i - bvirt; around = cdxt_ady1 - avirt; _0 = around + bround; _i = (double) (_0 - cdyt_adx1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cdyt_adx1; around = _0 - avirt; ct_a[1] = around + bround; ct_alarge = (double) (_j + _i); bvirt = (double) (ct_alarge - _j); avirt = ct_alarge - bvirt; bround = _i - bvirt; around = _j - avirt; ct_a[2] = around + bround;

				ct_a[3] = ct_alarge;
				ct_alen = 4;
				cdyt_bdx1 = (double) (cdytail * bdx); c = (double) (splitter * cdytail); abig = (double) (c - cdytail); ahi = c - abig; alo = cdytail - ahi; c = (double) (splitter * bdx); abig = (double) (c - bdx); bhi = c - abig; blo = bdx - bhi; err1 = cdyt_bdx1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdyt_bdx0 = (alo * blo) - err3;
				cdxt_bdy1 = (double) (cdxtail * bdy); c = (double) (splitter * cdxtail); abig = (double) (c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (double) (splitter * bdy); abig = (double) (c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = cdxt_bdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxt_bdy0 = (alo * blo) - err3;
				_i = (double) (cdyt_bdx0 - cdxt_bdy0); bvirt = (double) (cdyt_bdx0 - _i); avirt = _i + bvirt; bround = bvirt - cdxt_bdy0; around = cdyt_bdx0 - avirt; ct_b[0] = around + bround; _j = (double) (cdyt_bdx1 + _i); bvirt = (double) (_j - cdyt_bdx1); avirt = _j - bvirt; bround = _i - bvirt; around = cdyt_bdx1 - avirt; _0 = around + bround; _i = (double) (_0 - cdxt_bdy1); bvirt = (double) (_0 - _i); avirt = _i + bvirt; bround = bvirt - cdxt_bdy1; around = _0 - avirt; ct_b[1] = around + bround; ct_blarge = (double) (_j + _i); bvirt = (double) (ct_blarge - _j); avirt = ct_blarge - bvirt; bround = _i - bvirt; around = _j - avirt; ct_b[2] = around + bround;

				ct_b[3] = ct_blarge;
				ct_blen = 4;
			}
		}

		bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
		wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
		finswap = finnow; finnow = finother; finother = finswap;

		catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
		wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
		finswap = finnow; finnow = finother; finother = finswap;

		abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
		wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
		finswap = finnow; finnow = finother; finother = finswap;

		if (adztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, bc, adztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
					finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
					finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
					finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}

		if (adxtail != 0.0) {
			if (bdytail != 0.0) {
				adxt_bdyt1 = (double) (adxtail * bdytail); c = (double) (splitter * adxtail); abig = (double) (c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (double) (splitter * bdytail); abig = (double) (c - bdytail); bhi = c - abig; blo = bdytail - bhi; err1 = adxt_bdyt1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxt_bdyt0 = (alo * blo) - err3;
				c = (double) (splitter * cdz); abig = (double) (c - cdz); bhi = c - abig; blo = cdz - bhi; _i = (double) (adxt_bdyt0 * cdz); c = (double) (splitter * adxt_bdyt0); abig = (double) (c - adxt_bdyt0); ahi = c - abig; alo = adxt_bdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (adxt_bdyt1 * cdz); c = (double) (splitter * adxt_bdyt1); abig = (double) (c - adxt_bdyt1); ahi = c - abig; alo = adxt_bdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (cdztail != 0.0) {
					c = (double) (splitter * cdztail); abig = (double) (c - cdztail); bhi = c - abig; blo = cdztail - bhi; _i = (double) (adxt_bdyt0 * cdztail); c = (double) (splitter * adxt_bdyt0); abig = (double) (c - adxt_bdyt0); ahi = c - abig; alo = adxt_bdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (adxt_bdyt1 * cdztail); c = (double) (splitter * adxt_bdyt1); abig = (double) (c - adxt_bdyt1); ahi = c - abig; alo = adxt_bdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
			if (cdytail != 0.0) {
				negate = -adxtail;
				adxt_cdyt1 = (double) (negate * cdytail); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * cdytail); abig = (double) (c - cdytail); bhi = c - abig; blo = cdytail - bhi; err1 = adxt_cdyt1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxt_cdyt0 = (alo * blo) - err3;
				c = (double) (splitter * bdz); abig = (double) (c - bdz); bhi = c - abig; blo = bdz - bhi; _i = (double) (adxt_cdyt0 * bdz); c = (double) (splitter * adxt_cdyt0); abig = (double) (c - adxt_cdyt0); ahi = c - abig; alo = adxt_cdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (adxt_cdyt1 * bdz); c = (double) (splitter * adxt_cdyt1); abig = (double) (c - adxt_cdyt1); ahi = c - abig; alo = adxt_cdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (bdztail != 0.0) {
					c = (double) (splitter * bdztail); abig = (double) (c - bdztail); bhi = c - abig; blo = bdztail - bhi; _i = (double) (adxt_cdyt0 * bdztail); c = (double) (splitter * adxt_cdyt0); abig = (double) (c - adxt_cdyt0); ahi = c - abig; alo = adxt_cdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (adxt_cdyt1 * bdztail); c = (double) (splitter * adxt_cdyt1); abig = (double) (c - adxt_cdyt1); ahi = c - abig; alo = adxt_cdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
		}
		if (bdxtail != 0.0) {
			if (cdytail != 0.0) {
				bdxt_cdyt1 = (double) (bdxtail * cdytail); c = (double) (splitter * bdxtail); abig = (double) (c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (double) (splitter * cdytail); abig = (double) (c - cdytail); bhi = c - abig; blo = cdytail - bhi; err1 = bdxt_cdyt1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxt_cdyt0 = (alo * blo) - err3;
				c = (double) (splitter * adz); abig = (double) (c - adz); bhi = c - abig; blo = adz - bhi; _i = (double) (bdxt_cdyt0 * adz); c = (double) (splitter * bdxt_cdyt0); abig = (double) (c - bdxt_cdyt0); ahi = c - abig; alo = bdxt_cdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (bdxt_cdyt1 * adz); c = (double) (splitter * bdxt_cdyt1); abig = (double) (c - bdxt_cdyt1); ahi = c - abig; alo = bdxt_cdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (adztail != 0.0) {
					c = (double) (splitter * adztail); abig = (double) (c - adztail); bhi = c - abig; blo = adztail - bhi; _i = (double) (bdxt_cdyt0 * adztail); c = (double) (splitter * bdxt_cdyt0); abig = (double) (c - bdxt_cdyt0); ahi = c - abig; alo = bdxt_cdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (bdxt_cdyt1 * adztail); c = (double) (splitter * bdxt_cdyt1); abig = (double) (c - bdxt_cdyt1); ahi = c - abig; alo = bdxt_cdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
			if (adytail != 0.0) {
				negate = -bdxtail;
				bdxt_adyt1 = (double) (negate * adytail); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * adytail); abig = (double) (c - adytail); bhi = c - abig; blo = adytail - bhi; err1 = bdxt_adyt1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxt_adyt0 = (alo * blo) - err3;
				c = (double) (splitter * cdz); abig = (double) (c - cdz); bhi = c - abig; blo = cdz - bhi; _i = (double) (bdxt_adyt0 * cdz); c = (double) (splitter * bdxt_adyt0); abig = (double) (c - bdxt_adyt0); ahi = c - abig; alo = bdxt_adyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (bdxt_adyt1 * cdz); c = (double) (splitter * bdxt_adyt1); abig = (double) (c - bdxt_adyt1); ahi = c - abig; alo = bdxt_adyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (cdztail != 0.0) {
					c = (double) (splitter * cdztail); abig = (double) (c - cdztail); bhi = c - abig; blo = cdztail - bhi; _i = (double) (bdxt_adyt0 * cdztail); c = (double) (splitter * bdxt_adyt0); abig = (double) (c - bdxt_adyt0); ahi = c - abig; alo = bdxt_adyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (bdxt_adyt1 * cdztail); c = (double) (splitter * bdxt_adyt1); abig = (double) (c - bdxt_adyt1); ahi = c - abig; alo = bdxt_adyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
		}
		if (cdxtail != 0.0) {
			if (adytail != 0.0) {
				cdxt_adyt1 = (double) (cdxtail * adytail); c = (double) (splitter * cdxtail); abig = (double) (c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (double) (splitter * adytail); abig = (double) (c - adytail); bhi = c - abig; blo = adytail - bhi; err1 = cdxt_adyt1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxt_adyt0 = (alo * blo) - err3;
				c = (double) (splitter * bdz); abig = (double) (c - bdz); bhi = c - abig; blo = bdz - bhi; _i = (double) (cdxt_adyt0 * bdz); c = (double) (splitter * cdxt_adyt0); abig = (double) (c - cdxt_adyt0); ahi = c - abig; alo = cdxt_adyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (cdxt_adyt1 * bdz); c = (double) (splitter * cdxt_adyt1); abig = (double) (c - cdxt_adyt1); ahi = c - abig; alo = cdxt_adyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (bdztail != 0.0) {
					c = (double) (splitter * bdztail); abig = (double) (c - bdztail); bhi = c - abig; blo = bdztail - bhi; _i = (double) (cdxt_adyt0 * bdztail); c = (double) (splitter * cdxt_adyt0); abig = (double) (c - cdxt_adyt0); ahi = c - abig; alo = cdxt_adyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (cdxt_adyt1 * bdztail); c = (double) (splitter * cdxt_adyt1); abig = (double) (c - cdxt_adyt1); ahi = c - abig; alo = cdxt_adyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
			if (bdytail != 0.0) {
				negate = -cdxtail;
				cdxt_bdyt1 = (double) (negate * bdytail); c = (double) (splitter * negate); abig = (double) (c - negate); ahi = c - abig; alo = negate - ahi; c = (double) (splitter * bdytail); abig = (double) (c - bdytail); bhi = c - abig; blo = bdytail - bhi; err1 = cdxt_bdyt1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxt_bdyt0 = (alo * blo) - err3;
				c = (double) (splitter * adz); abig = (double) (c - adz); bhi = c - abig; blo = adz - bhi; _i = (double) (cdxt_bdyt0 * adz); c = (double) (splitter * cdxt_bdyt0); abig = (double) (c - cdxt_bdyt0); ahi = c - abig; alo = cdxt_bdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (cdxt_bdyt1 * adz); c = (double) (splitter * cdxt_bdyt1); abig = (double) (c - cdxt_bdyt1); ahi = c - abig; alo = cdxt_bdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
				finswap = finnow; finnow = finother; finother = finswap;
				if (adztail != 0.0) {
					c = (double) (splitter * adztail); abig = (double) (c - adztail); bhi = c - abig; blo = adztail - bhi; _i = (double) (cdxt_bdyt0 * adztail); c = (double) (splitter * cdxt_bdyt0); abig = (double) (c - cdxt_bdyt0); ahi = c - abig; alo = cdxt_bdyt0 - ahi; err1 = _i - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); u[0] = (alo * blo) - err3; _j = (double) (cdxt_bdyt1 * adztail); c = (double) (splitter * cdxt_bdyt1); abig = (double) (c - cdxt_bdyt1); ahi = c - abig; alo = cdxt_bdyt1 - ahi; err1 = _j - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); _0 = (alo * blo) - err3; _k = (double) (_i + _0); bvirt = (double) (_k - _i); avirt = _k - bvirt; bround = _0 - bvirt; around = _i - avirt; u[1] = around + bround; u3 = (double) (_j + _k); bvirt = u3 - _j; u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
					finswap = finnow; finnow = finother; finother = finswap;
				}
			}
		}

		if (adztail != 0.0) {
			wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
					finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (bdztail != 0.0) {
			wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
					finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}
		if (cdztail != 0.0) {
			wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
					finother);
			finswap = finnow; finnow = finother; finother = finswap;
		}

		return finnow[finlength - 1];
	}
	private int fast_expansion_sum_zeroelim(int elen, double[] e, int flen, double[] f, double[] h)
	{
		double Q;
		double Qnew;
		double hh;
		double bvirt;
		double avirt, bround, around;
		int eindex, findex, hindex;
		double enow, fnow;

		enow = e[0];
		fnow = f[0];
		eindex = findex = 0;
		if ((fnow > enow) == (fnow > -enow)) {
			Q = enow;
			enow = e[++eindex];
		} else {
			Q = fnow;
			fnow = f[++findex];
		}
		hindex = 0;
		if ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Qnew = (double) (enow + Q); bvirt = Qnew - enow; hh = Q - bvirt;
				enow = e[++eindex];
			} else {
				Qnew = (double) (fnow + Q); bvirt = Qnew - fnow; hh = Q - bvirt;
				fnow = f[++findex];
			}
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
			while ((eindex < elen) && (findex < flen)) {
				if ((fnow > enow) == (fnow > -enow)) {
					Qnew = (double) (Q + enow); bvirt = (double) (Qnew - Q); avirt = Qnew - bvirt; bround = enow - bvirt; around = Q - avirt; hh = around + bround;
					enow = e[++eindex];
				} else {
					Qnew = (double) (Q + fnow); bvirt = (double) (Qnew - Q); avirt = Qnew - bvirt; bround = fnow - bvirt; around = Q - avirt; hh = around + bround;
					fnow = f[++findex];
				}
				Q = Qnew;
				if (hh != 0.0) {
					h[hindex++] = hh;
				}
			}
		}
		while (eindex < elen) {
			Qnew = (double) (Q + enow); bvirt = (double) (Qnew - Q); avirt = Qnew - bvirt; bround = enow - bvirt; around = Q - avirt; hh = around + bround;
			enow = e[++eindex];
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
		while (findex < flen) {
			Qnew = (double) (Q + fnow); bvirt = (double) (Qnew - Q); avirt = Qnew - bvirt; bround = fnow - bvirt; around = Q - avirt; hh = around + bround;
			fnow = f[++findex];
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	private int scale_expansion_zeroelim(int elen, double[] e, double b,double[] h)
	{
		double Q, sum;
		double hh;
		double product1;
		double product0;
		int eindex, hindex;
		double enow;
		double bvirt;
		double avirt, bround, around;
		double c;
		double abig;
		double ahi, alo, bhi, blo;
		double err1, err2, err3;

		c = (double) (splitter * b); abig = (double) (c - b); bhi = c - abig; blo = b - bhi;
		Q = (double) (e[0] * b); c = (double) (splitter * e[0]); abig = (double) (c - e[0]); ahi = c - abig; alo = e[0] - ahi; err1 = Q - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); hh = (alo * blo) - err3;
		hindex = 0;
		if (hh != 0) {
			h[hindex++] = hh;
		}
		for (eindex = 1; eindex < elen; eindex++) {
			enow = e[eindex];
			product1 = (double) (enow * b); c = (double) (splitter * enow); abig = (double) (c - enow); ahi = c - abig; alo = enow - ahi; err1 = product1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); product0 = (alo * blo) - err3;
			sum = (double) (Q + product0); bvirt = (double) (sum - Q); avirt = sum - bvirt; bround = product0 - bvirt; around = Q - avirt; hh = around + bround;
			if (hh != 0) {
				h[hindex++] = hh;
			}
			Q = (double) (product1 + sum); bvirt = Q - product1; hh = sum - bvirt;
			if (hh != 0) {
				h[hindex++] = hh;
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	private double estimate(int elen, double[] e)
	{
		double Q;
		int eindex;

		Q = e[0];
		for (eindex = 1; eindex < elen; eindex++) {
			Q += e[eindex];
		}
		return Q;
	}

}
