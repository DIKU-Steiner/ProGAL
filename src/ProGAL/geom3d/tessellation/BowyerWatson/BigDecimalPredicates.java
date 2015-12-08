package ProGAL.geom3d.tessellation.BowyerWatson;

import java.math.BigDecimal;

public class BigDecimalPredicates {
	
	private static BigDecimalPredicates instance;
	public static BigDecimalPredicates getInstance(){
		if(instance==null) instance = new BigDecimalPredicates();
		return instance;
	}
	
	private BigDecimalPredicates(){
	}
	
	BigDecimal insphere(double[] pa,double[]  pb,double[]  pc,double[]  pd,double[]  pe)
	{
		BigDecimal aex = new BigDecimal(pa[0] - pe[0]);
		BigDecimal bex = new BigDecimal(pb[0] - pe[0]);
		BigDecimal cex = new BigDecimal(pc[0] - pe[0]);
		BigDecimal dex = new BigDecimal(pd[0] - pe[0]);
		BigDecimal aey = new BigDecimal(pa[1] - pe[1]);
		BigDecimal bey = new BigDecimal(pb[1] - pe[1]);
		BigDecimal cey = new BigDecimal(pc[1] - pe[1]);
		BigDecimal dey = new BigDecimal(pd[1] - pe[1]);
		BigDecimal aez = new BigDecimal(pa[2] - pe[2]);
		BigDecimal bez = new BigDecimal(pb[2] - pe[2]);
		BigDecimal cez = new BigDecimal(pc[2] - pe[2]);
		BigDecimal dez = new BigDecimal(pd[2] - pe[2]);

		BigDecimal aexbey = aex.multiply(bey);
		BigDecimal bexaey = bex.multiply(aey);
		BigDecimal ab = aexbey.subtract(bexaey);
		BigDecimal bexcey = bex.multiply(cey);
		BigDecimal cexbey = cex.multiply(bey);
		BigDecimal bc = bexcey.subtract(cexbey);
		BigDecimal cexdey = cex.multiply(dey);
		BigDecimal dexcey = dex.multiply(cey);
		BigDecimal cd = cexdey.subtract(dexcey);
		BigDecimal dexaey = dex.multiply(aey);
		BigDecimal aexdey = aex.multiply(dey);
		BigDecimal da = dexaey.subtract(aexdey);

		BigDecimal aexcey = aex.multiply(cey);
		BigDecimal cexaey = cex.multiply(aey);
		BigDecimal ac = aexcey.subtract(cexaey);
		BigDecimal bexdey = bex.multiply(dey);
		BigDecimal dexbey = dex.multiply(bey);
		BigDecimal bd = bexdey.subtract(dexbey);

		BigDecimal abc = aez.multiply(bc).subtract(bez.multiply(ac)).add(cez.multiply(ab));
		BigDecimal bcd = bez.multiply(cd).subtract(cez.multiply(bd)).add(dez.multiply(bc));
		BigDecimal cda = cez.multiply(da).add(dez.multiply(ac)).add(aez.multiply(cd));
		BigDecimal dab = dez.multiply(ab).add(aez.multiply(bd)).add(bez.multiply(da));

		BigDecimal alift = aex.multiply(aex).add(aey.multiply(aey)).add(aez.multiply(aez));
		BigDecimal blift = bex.multiply(bex).add(bey.multiply(bey)).add(bez.multiply(bez));
		BigDecimal clift = cex.multiply(cex).add(cey.multiply(cey)).add(cez.multiply(cez));
		BigDecimal dlift = dex.multiply(dex).add(dey.multiply(dey)).add(dez.multiply(dez));

		BigDecimal det = (dlift.multiply(abc).subtract(clift.multiply(dab))).add(blift.multiply(cda).subtract(alift.multiply(bcd)));
		return det;
	}

	

	BigDecimal orient(double[] pa, double[] pb, double[] pc, double[] pd)	{
		BigDecimal adx = new BigDecimal(pa[0] - pd[0]);
		BigDecimal bdx = new BigDecimal(pb[0] - pd[0]);
		BigDecimal cdx = new BigDecimal(pc[0] - pd[0]);
		BigDecimal ady = new BigDecimal(pa[1] - pd[1]);
		BigDecimal bdy = new BigDecimal(pb[1] - pd[1]);
		BigDecimal cdy = new BigDecimal(pc[1] - pd[1]);
		BigDecimal adz = new BigDecimal(pa[2] - pd[2]);
		BigDecimal bdz = new BigDecimal(pb[2] - pd[2]);
		BigDecimal cdz = new BigDecimal(pc[2] - pd[2]);

		BigDecimal bdxcdy = bdx.multiply(cdy);
		BigDecimal cdxbdy = cdx.multiply(bdy);
		BigDecimal cdxady = cdx.multiply(ady);
		BigDecimal adxcdy = adx.multiply(cdy);
		BigDecimal adxbdy = adx.multiply(bdy);
		BigDecimal bdxady = bdx.multiply(ady);
		BigDecimal det = adz.multiply(bdxcdy.subtract(cdxbdy)).add(bdz.multiply(cdxady.subtract(adxcdy))).add(cdz.multiply(adxbdy.subtract(bdxady)));
		return det;
	}

}
