package ProGAL.geom3d.complex.alphaComplex;


public class TriangleAlphaProperties implements SimplexAlphaProperties {
	private Interval singular, regular, interior;
    private boolean attached;
    private boolean onConvexHull;
    private double inAlphaComplex;

    public TriangleAlphaProperties(Interval i, Interval r, Interval s, boolean a, double inAlpha){
    	this.singular = s;
    	this.regular = r;
    	this.interior = i;
    	this.attached = a;
    	this.inAlphaComplex = inAlpha;
    }

    public TriangleAlphaProperties(double muDown, double muUp,double rho, boolean ch, boolean a){
    	if(a)	this.singular = null;
    	else	this.singular = new Interval(rho,muDown);
    	if(ch)	this.regular = new Interval(muDown, Double.POSITIVE_INFINITY);
    	else	this.regular = new Interval(muDown, muUp);
    	if(ch)	this.interior = null;
    	else	this.interior = new Interval(muUp, Double.POSITIVE_INFINITY);
    	this.attached = a;
    	this.onConvexHull = ch;
    	if(a) 	this.inAlphaComplex = muDown;
    	else	this.inAlphaComplex = rho;
    }
    
    
    public boolean isAttached(){ 			return attached; 			}
    public boolean getOnConvexHull(){ 		return onConvexHull; 		}
    public double getInAlphaComplex(){ 		return inAlphaComplex; 		}

    public Interval getSingularInterval(){	return singular; 			}
    public Interval getRegularInterval(){	return regular; 			}
    public Interval getInteriorInterval(){	return interior; 			}
    
	public int getSimplexType() {	return 2;	}
}
