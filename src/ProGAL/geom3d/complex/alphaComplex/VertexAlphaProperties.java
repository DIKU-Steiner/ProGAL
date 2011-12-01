package ProGAL.geom3d.complex.alphaComplex;


public class VertexAlphaProperties implements SimplexAlphaProperties {
	private final Interval singular, regular, interior;
    private final boolean onConvexHull;
    private final double inAlphaComplex;

    public VertexAlphaProperties(Interval i, Interval r, Interval s, boolean ch, double inAlpha){
    	this.singular = s;
    	this.regular = r;
    	this.interior = i;
    	this.onConvexHull = ch;
    	this.inAlphaComplex = inAlpha;
    }
    

    public VertexAlphaProperties(double muDown, double muUp, boolean ch){
    	this.singular = new Interval(0,muDown);
    	if(ch){
    		this.regular = new Interval(muDown,Double.POSITIVE_INFINITY);
    		this.interior = null;
    	}else{
    		this.regular = new Interval(muDown,muUp);
    		this.interior = new Interval(muUp,Double.POSITIVE_INFINITY);
    	}
    	this.onConvexHull = ch;
    	this.inAlphaComplex = 0;
    }
    
    
    public boolean getOnConvexHull(){ 		return onConvexHull; 		}
    public double getInAlphaComplex(){ 		return inAlphaComplex; 		}
    public boolean isAttached(){			return false; 				}

    public Interval getSingularInterval(){	return singular; 			}
    public Interval getRegularInterval(){	return regular; 			}
    public Interval getInteriorInterval(){	return interior; 			}
    
	public int getSimplexType() {	return 0;	}
}
