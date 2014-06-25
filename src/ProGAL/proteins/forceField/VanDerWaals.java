package ProGAL.proteins.forceField;

//import java.awt.Color;
//import java.io.File;
//import java.io.IOException;
//import java.text.DecimalFormat;
//import java.util.ArrayList;
//import java.util.HashSet;
//
//import jxl.*;
//import jxl.write.*;
//import jxl.write.Number;
//
//import ProGAL.dataStructures.Graph;
//import ProGAL.geom3d.kineticDelaunay.KineticAlphaComplex;
//import ProGAL.geom3d.kineticDelaunay.KineticAlphaComplex.ProblemInstanceType;
//import ProGAL.geom3d.kineticDelaunay.Vertex;
//import ProGAL.geom3d.kineticDelaunay.Tet;
//import ProGAL.geom3d.viewer.J3DScene;
//import ProGAL.io.IOToolbox;
//import ProGAL.proteins.structure.AminoAcid;

public class VanDerWaals {
//
//	// van der Waals radii (rA) and van der Waals well depths (eA) from the AMBER force field
//	double rN = 1.8240;                 double eN = 0.1700;
//	double rC = 1.9080;                 double eC = 0.1094;
//	double rO = 1.6612;                 double eO = 0.2100;
//	double rP = 2.1000;                 double eP = 0.2000;
//	double rS = 2.0000;                 double eS = 0.2500;
//	double rH = 0.6000;                 double eH = 0.0157;
//	
////								      		  double rNN2 = rN*rN;		  
////	double rNC = Math.sqrt(rN * rC);		  double rNC2 = rNC*rNC;	
////	double rNO = Math.sqrt(rN * rO);		  double rNO2 = rNO*rNO;	
////	double rNP = Math.sqrt(rN * rP);		  double rNP2 = rNP*rNP;	
////	double rNS = Math.sqrt(rN * rS);          double rNS2 = rNS*rNS;
////	double rNH = Math.sqrt(rN * rH);		  double rNH2 = rNH*rNH;	
////									  		  double rCC2 = rC*rC;
////	double rCO = Math.sqrt(rC * rO);		  double rCO2 = rCO*rCO;	
////	double rCP = Math.sqrt(rC * rP);		  double rCP2 = rCP*rCP;	
////	double rCS = Math.sqrt(rC * rS);   		  double rCS2 = rCS*rCS;
////	double rCH = Math.sqrt(rC * rH);		  double rCH2 = rCH*rCH;
////											  double rOO2 = rO*rO;
////	double rOP = Math.sqrt(rO * rP);		  double rOP2 = rOP*rOP;
////	double rOS = Math.sqrt(rO * rS);   	 	  double rOS2 = rOS*rOS;
////	double rOH = Math.sqrt(rO * rH);		  double rOH2 = rOH*rOH;
////							     	  		  double rPP2 = rP*rP;	
////	double rPS = Math.sqrt(rP * rS);   	 	  double rPS2 = rPS*rPS;
////	double rPH = Math.sqrt(rP * rH);		  double rPH2 = rPH*rPH;
////								  	  		  double rSS2 = rS*rS;	
////	double rSH = Math.sqrt(rS * rH);	      double rSH2 = rSH*rSH;				
////											  double rHH2 = rH*rH;
//
//								      
//								      double rNN2 = rN*rN;		  
//	double rNC = (rN + rC)/2;		  double rNC2 = rNC*rNC;	
//	double rNO = (rN + rO)/2;		  double rNO2 = rNO*rNO;	
//	double rNP = (rN + rP)/2;		  double rNP2 = rNP*rNP;	
//	double rNS = (rN + rS)/2;         double rNS2 = rNS*rNS;
//	double rNH = (rN + rH)/2;		  double rNH2 = rNH*rNH;	
//									  double rCC2 = rC*rC;
//	double rCO = (rC + rO)/2;		  double rCO2 = rCO*rCO;	
//	double rCP = (rC + rP)/2;		  double rCP2 = rCP*rCP;	
//	double rCS = (rC + rS)/2;   	  double rCS2 = rCS*rCS;
//	double rCH = (rC + rH)/2;		  double rCH2 = rCH*rCH;
//									  double rOO2 = rO*rO;
//	double rOP = (rO + rP)/2;		  double rOP2 = rOP*rOP;
//	double rOS = (rO + rS)/2;   	  double rOS2 = rOS*rOS;
//	double rOH = (rO + rH)/2;		  double rOH2 = rOH*rOH;
//							     	  double rPP2 = rP*rP;	
//	double rPS = (rP + rS)/2;   	  double rPS2 = rPS*rPS;
//	double rPH = (rP + rH)/2;		  double rPH2 = rPH*rPH;
//								  	  double rSS2 = rS*rS;	
//	double rSH = (rS + rH)/2;	      double rSH2 = rSH*rSH;				
//								      double rHH2 = rH*rH;
//								      
//
//double eNC = Math.sqrt(eN*eC);	
//double eNO = Math.sqrt(eN*eO);	
//double eNP = Math.sqrt(eN*eP);	
//double eNS = Math.sqrt(eN*eS);	
//double eNH = Math.sqrt(eN*eH);	
//
//double eCO = Math.sqrt(eC*eO);	
//double eCP = Math.sqrt(eC*eP);	
//double eCS = Math.sqrt(eC*eS);	
//double eCH = Math.sqrt(eC*eH);	
//
//double eOP = Math.sqrt(eO*eP);	
//double eOS = Math.sqrt(eO*eS);	
//double eOH = Math.sqrt(eO*eH);	
//	        
//double ePS = Math.sqrt(eP*eS);	
//double ePH = Math.sqrt(eP*eH);	
//
//double eSH = Math.sqrt(eS*eH);	
//
//	
//	
//	
//	double[][] r = {{rNN2, rNC2, rNO2, rNP2, rNS2, rNH2},
//					{rNC2, rCC2, rCO2, rCP2, rCS2, rCH2},
//					{rNO2, rCO2, rOO2, rOP2, rOS2, rOH2},
//					{rNP2, rCP2, rOP2, rPP2, rPS2, rPH2},
//					{rNS2, rCS2, rOS2, rPS2, rSS2, rSH2},
//					{rNH2, rCH2, rOH2, rPH2, rSH2, rHH2}};
//
//	double[][] e = {{eN,  eNC, eNO, eNP, eNS, eNH},
//					{eNC, eC,  eCO, eCP, eCS, eCH},
//					{eNO, eCO, eO,  eOP, eOS, eOH},
//					{eNP, eCP, eOP, eP,  ePS, ePH},
//					{eNS, eCS, eOS, ePS, eS,  eSH},
//					{eNH, eCH, eOH, ePH, eSH, eH}};
//
//	Vertex u, v;
//	
//	private void setUpCharges(KineticAlphaComplex kDT) {
//		int i = 4; 
//		while (i < kDT.getNrVertices()) {
//			System.out.println(i + ": " + kDT.getVertex(i).atomName);
//			if (!kDT.getVertex(i).atomName.equals("N"))
//				System.out.println("STOP");
//			
//			char aminoAcid = kDT.getVertex(i).aaType; 
//			switch (aminoAcid) {
//				case 'A':  
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e =  0.0337; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e = -0.1825; // CB
//					kDT.getVertex(i++).e =  0.2719; // H bonded with N
//					kDT.getVertex(i++).e =  0.0823; // HA bonded with CA
//					kDT.getVertex(i++).e =  0.0603; // HB1 bonded with CB
//					kDT.getVertex(i++).e =  0.0603; // HB2 bonded with CB
//					kDT.getVertex(i++).e =  0.0603; // HB3 bonded with CB
//					break;
//				case 'R':
//					kDT.getVertex(i++).e = -0.3479; // N
//					kDT.getVertex(i++).e = -0.2637; // CA
//					kDT.getVertex(i++).e =  0.7341; // C
//					kDT.getVertex(i++).e = -0.5894; // O
//					kDT.getVertex(i++).e = -0.0007; // CB
//					kDT.getVertex(i++).e =  0.0390; // CG
//					kDT.getVertex(i++).e =  0.0486; // CD
//					kDT.getVertex(i++).e = -0.5295; // NE
//					kDT.getVertex(i++).e =  0.8076; // CZ
//					kDT.getVertex(i++).e = -0.8627; // NH1
//					kDT.getVertex(i++).e = -0.8627; // NH2
//					kDT.getVertex(i++).e =  0.2747; // H
//					kDT.getVertex(i++).e =  0.1560; // HA
//					kDT.getVertex(i++).e =  0.0327; // HB2
//					kDT.getVertex(i++).e =  0.0327; // HB3
//					kDT.getVertex(i++).e =  0.0285; // HG2
//					kDT.getVertex(i++).e =  0.0285; // HG3
//					kDT.getVertex(i++).e =  0.0687; // HD2
//					kDT.getVertex(i++).e =  0.0687; // HD3
//					kDT.getVertex(i++).e =  0.3456; // HE
//					kDT.getVertex(i++).e =  0.4478; // HH11
//					kDT.getVertex(i++).e =  0.4478; // HH12
//					kDT.getVertex(i++).e =  0.4478; // HH21
//					kDT.getVertex(i++).e =  0.4478; // HH22
//					break;
//				case 'N':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e =  0.0143; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e = -0.2041; // CB
//					kDT.getVertex(i++).e =  0.7130; // CG
//					kDT.getVertex(i++).e = -0.5931; // OD1
//					kDT.getVertex(i++).e = -0.9191; // ND2
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.1048; // HA
//					kDT.getVertex(i++).e =  0.0797; // HB2
//					kDT.getVertex(i++).e =  0.0797; // HB3
//					kDT.getVertex(i++).e =  0.4196; // HD21
//					kDT.getVertex(i++).e =  0.4196; // HD22
//					break;
//				case 'D': 
//					kDT.getVertex(i++).e = -0.5163; // N
//					kDT.getVertex(i++).e =  0.0381; // CA
//					kDT.getVertex(i++).e =  0.5366; // C
//					kDT.getVertex(i++).e = -0.5819; // O
//					kDT.getVertex(i++).e = -0.0303; // CB
//					kDT.getVertex(i++).e =  0.7994; // CG
//					kDT.getVertex(i++).e = -0.8014; // OD1
//					kDT.getVertex(i++).e = -0.8014; // OD2
//					kDT.getVertex(i++).e =  0.2936; // H
//					kDT.getVertex(i++).e =  0.0880; // HA
//					kDT.getVertex(i++).e = -0.0122; // HB2
//					kDT.getVertex(i++).e = -0.0122; // HB3
//					break;
//				case 'C':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e =  0.0213; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e = -0.1231; // CB
//					kDT.getVertex(i++).e = -0.3119; // SG
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.1124; // HA
//					kDT.getVertex(i++).e =  0.1112; // HB2
//					kDT.getVertex(i++).e =  0.1112; // HB3
////					kDT.getVertex(i++).e =  0.1933; // HG
//					break;
//				case 'Q':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0031; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e = -0.0036; // CB
//					kDT.getVertex(i++).e = -0.0645; // CG
//					kDT.getVertex(i++).e =  0.6951; // CD
//					kDT.getVertex(i++).e = -0.6086; // OE1
//					kDT.getVertex(i++).e = -0.9407; // NE2
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0850; // HA
//					kDT.getVertex(i++).e =  0.0171; // HB2
//					kDT.getVertex(i++).e =  0.0171; // HB3
//					kDT.getVertex(i++).e =  0.0352; // HG2
//					kDT.getVertex(i++).e =  0.0352; // HG3
//					kDT.getVertex(i++).e =  0.4251; // HE21
//					kDT.getVertex(i++).e =  0.4251; // HE22
//					break;
//				case 'E':
//					kDT.getVertex(i++).e = -0.5163; // N
//					kDT.getVertex(i++).e =  0.0397; // CA
//					kDT.getVertex(i++).e =  0.5366; // C
//					kDT.getVertex(i++).e = -0.5819; // O
//					kDT.getVertex(i++).e =  0.0560; // CB
//					kDT.getVertex(i++).e =  0.0136; // CG
//					kDT.getVertex(i++).e =  0.8054; // CD
//					kDT.getVertex(i++).e = -0.8188; // OE1
//					kDT.getVertex(i++).e = -0.8188; // OE2
//					kDT.getVertex(i++).e =  0.2936; // H
//					kDT.getVertex(i++).e =  0.1105; // HA
//					kDT.getVertex(i++).e = -0.0173; // HB2
//					kDT.getVertex(i++).e = -0.0173; // HB3
//					kDT.getVertex(i++).e = -0.0425; // HG2
//					kDT.getVertex(i++).e = -0.0425; // HG3
//					break;
//
//				case 'G':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0252; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e =  0.2719; // H1 bonded with N
//					kDT.getVertex(i++).e =  0.0698; // HA2 bonded with CA
//					kDT.getVertex(i++).e =  0.0698; // HA3 bonded with CA
//					break;
//				case 'H':
//					kDT.getVertex(i++).e = -0.3479; // N
//					kDT.getVertex(i++).e = -0.1354; // CA
//					kDT.getVertex(i++).e =  0.7341; // C
//					kDT.getVertex(i++).e = -0.5894; // O
//					kDT.getVertex(i++).e = -0.0414; // CB
//					kDT.getVertex(i++).e = -0.0012; // CG
//					kDT.getVertex(i++).e = -0.1513; // ND1
//					kDT.getVertex(i++).e = -0.1141; // CD2
//					kDT.getVertex(i++).e = -0.0170; // CE1
//					kDT.getVertex(i++).e = -0.1718; // NE2
//					kDT.getVertex(i++).e =  0.2747; // H
//					kDT.getVertex(i++).e =  0.1212; // HA
//					kDT.getVertex(i++).e =  0.0810; // HB2
//					kDT.getVertex(i++).e =  0.0810; // HB3
////					kDT.getVertex(i++).e =  0.3866; // HD1
//					kDT.getVertex(i++).e =  0.2317; // HD2
//					kDT.getVertex(i++).e =  0.2681; // HE1
////					kDT.getVertex(i++).e =  0.3911; // HE2
//					break;
//				case 'I':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0597; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e =  0.1303; // CB
//					kDT.getVertex(i++).e = -0.0430; // CG1
//					kDT.getVertex(i++).e = -0.3204; // CG2
//					kDT.getVertex(i++).e = -0.0660; // CD1
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0869; // HA
//					kDT.getVertex(i++).e =  0.0187; // HB
//					kDT.getVertex(i++).e =  0.0236; // HG12
//					kDT.getVertex(i++).e =  0.0236; // HG13
//					kDT.getVertex(i++).e =  0.0882; // HG21
//					kDT.getVertex(i++).e =  0.0882; // HG22
//					kDT.getVertex(i++).e =  0.0882; // HG23
//					kDT.getVertex(i++).e =  0.0186; // HD11
//					kDT.getVertex(i++).e =  0.0186; // HD12
//					kDT.getVertex(i++).e =  0.0186; // HD13
//					break;
//				case 'L':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0518; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5678; // O
//					kDT.getVertex(i++).e = -0.1102; // CB
//					kDT.getVertex(i++).e =  0.3531; // CG
//					kDT.getVertex(i++).e = -0.4121; // CD1
//					kDT.getVertex(i++).e = -0.4121; // CD2
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0922; // HA
//					kDT.getVertex(i++).e =  0.0457; // HB2
//					kDT.getVertex(i++).e =  0.0457; // HB3
//					kDT.getVertex(i++).e = -0.0361; // HG
//					kDT.getVertex(i++).e =  0.1000; // HD11
//					kDT.getVertex(i++).e =  0.1000; // HD12
//					kDT.getVertex(i++).e =  0.1000; // HD13
//					kDT.getVertex(i++).e =  0.1000; // HD21
//					kDT.getVertex(i++).e =  0.1000; // HD22
//					kDT.getVertex(i++).e =  0.1000; // HD23
//					break;
//				case 'K':
//					kDT.getVertex(i++).e = -0.3479; // N
//					kDT.getVertex(i++).e = -0.2400; // CA
//					kDT.getVertex(i++).e =  0.7341; // C
//					kDT.getVertex(i++).e = -0.3894; // O
//					kDT.getVertex(i++).e = -0.0094; // CB
//					kDT.getVertex(i++).e =  0.0187; // CG
//					kDT.getVertex(i++).e = -0.0479; // CD
//					kDT.getVertex(i++).e = -0.0143; // CE
//					kDT.getVertex(i++).e = -0.3854; // NZ
//					kDT.getVertex(i++).e =  0.2747; // H
//					kDT.getVertex(i++).e =  0.1426; // HA
//					kDT.getVertex(i++).e =  0.0362; // HB2
//					kDT.getVertex(i++).e =  0.0362; // HB3
//					kDT.getVertex(i++).e =  0.0103; // HG2
//					kDT.getVertex(i++).e =  0.0103; // HG3
//					kDT.getVertex(i++).e =  0.0621; // HD2
//					kDT.getVertex(i++).e =  0.0621; // HD3
//					kDT.getVertex(i++).e =  0.1135; // HE2
//					kDT.getVertex(i++).e =  0.1135; // HE3
//					kDT.getVertex(i++).e =  0.3400; // HZ1
//					kDT.getVertex(i++).e =  0.3400; // HZ2
//					kDT.getVertex(i++).e =  0.3400; // HZ3
//					break;
//				case 'M':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0237; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e =  0.0342; // CB
//					kDT.getVertex(i++).e =  0.0018; // CG
//					kDT.getVertex(i++).e = -0.2737; // SD
//					kDT.getVertex(i++).e = -0.0536; // CE
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0880; // HA
//					kDT.getVertex(i++).e =  0.0241; // HB2
//					kDT.getVertex(i++).e =  0.0241; // HB3
//					kDT.getVertex(i++).e =  0.0440; // HG2
//					kDT.getVertex(i++).e =  0.0440; // HG3
//					kDT.getVertex(i++).e =  0.0684; // HE1
//					kDT.getVertex(i++).e =  0.0684; // HE2
//					kDT.getVertex(i++).e =  0.0684; // HE3
//					break;
//				case 'F':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0024; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5879; // O
//					kDT.getVertex(i++).e = -0.0343; // CB
//					kDT.getVertex(i++).e =  0.0118; // CG
//					kDT.getVertex(i++).e = -0.1256; // CD1
//					kDT.getVertex(i++).e = -0.1256; // CD2
//					kDT.getVertex(i++).e = -0.1704; // CE1
//					kDT.getVertex(i++).e = -0.1704; // CE2
//					kDT.getVertex(i++).e = -0.1072; // CZ
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0978; // HA
//					kDT.getVertex(i++).e =  0.0295; // HB2
//					kDT.getVertex(i++).e =  0.0295; // HB3
//					kDT.getVertex(i++).e =  0.1330; // HD1
//					kDT.getVertex(i++).e =  0.1330; // HD2
//					kDT.getVertex(i++).e =  0.1430; // HE1
//					kDT.getVertex(i++).e =  0.1430; // HE2
//					kDT.getVertex(i++).e =  0.1297; // HZ
//					break;
//				case 'P':
//					kDT.getVertex(i++).e = -0.2548; // N
//					kDT.getVertex(i++).e = -0.0266; // CA
//					kDT.getVertex(i++).e =  0.5896; // C
//					kDT.getVertex(i++).e = -0.5748; // O
//					kDT.getVertex(i++).e = -0.0070; // CB
//					kDT.getVertex(i++).e =  0.0189; // CG
//					kDT.getVertex(i++).e =  0.0192; // CD
//					kDT.getVertex(i++).e =  0.0641; // HA
//					kDT.getVertex(i++).e =  0.0253; // HB2
//					kDT.getVertex(i++).e =  0.0253; // HB3
//					kDT.getVertex(i++).e =  0.0213; // HG2
//					kDT.getVertex(i++).e =  0.0213; // HG3
//					kDT.getVertex(i++).e =  0.0391; // HD2
//					kDT.getVertex(i++).e =  0.0391; // HD3
//					break;
//				case 'S':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0249; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e =  0.2117; // CB
//					kDT.getVertex(i++).e = -0.6546; // OG
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0843; // HA
//					kDT.getVertex(i++).e =  0.0352; // HB2
//					kDT.getVertex(i++).e =  0.0352; // HB3
//					kDT.getVertex(i++).e =  0.4275; // HG
//					break;
//				case 'T':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0389; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e =  0.3654; // CB
//					kDT.getVertex(i++).e = -0.6761; // OG1
//					kDT.getVertex(i++).e = -0.2438; // CG2
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.1007; // HA
//					kDT.getVertex(i++).e =  0.0043; // HB
//					kDT.getVertex(i++).e =  0.4102; // HG1
//					kDT.getVertex(i++).e =  0.0642; // HG21
//					kDT.getVertex(i++).e =  0.0642; // HG22
//					kDT.getVertex(i++).e =  0.0642; // HG23
//					break;
//				case 'W':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0275; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e = -0.0050; // CB
//					kDT.getVertex(i++).e = -0.1415; // CG
//					kDT.getVertex(i++).e = -0.1638; // CD1
//					kDT.getVertex(i++).e =  0.1234; // CD2
//					kDT.getVertex(i++).e = -0.3418; // NE1
//					kDT.getVertex(i++).e =  0.1380; // CE2
//					kDT.getVertex(i++).e = -0.2387; // CE3
//					kDT.getVertex(i++).e = -0.2601; // CZ2
//					kDT.getVertex(i++).e = -0.1972; // CZ3
//					kDT.getVertex(i++).e = -0.1134; // CH2
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.1123; // HA
//					kDT.getVertex(i++).e =  0.0339; // HB2
//					kDT.getVertex(i++).e =  0.0339; // HB3
//					kDT.getVertex(i++).e =  0.2062; // HD1
//					kDT.getVertex(i++).e =  0.3412; // HE1
//					kDT.getVertex(i++).e =  0.1700; // HE3
//					kDT.getVertex(i++).e =  0.1572; // HZ2
//					kDT.getVertex(i++).e =  0.1447; // HZ3
//					kDT.getVertex(i++).e =  0.1417; // HH2
//					break;
//				case 'Y':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0014; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e = -0.0152; // CB
//					kDT.getVertex(i++).e = -0.0011; // CG
//					kDT.getVertex(i++).e = -0.1906; // CD1
//					kDT.getVertex(i++).e = -0.1906; // CD2
//					kDT.getVertex(i++).e = -0.2341; // CE1
//					kDT.getVertex(i++).e = -0.2341; // CE2
//					kDT.getVertex(i++).e =  0.3226; // CZ
//					kDT.getVertex(i++).e = -0.5579; // OH
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0876; // HA
//					kDT.getVertex(i++).e =  0.0295; // HB2
//					kDT.getVertex(i++).e =  0.0295; // HB3
//					kDT.getVertex(i++).e =  0.1689; // HD1
//					kDT.getVertex(i++).e =  0.1689; // HD2
//					kDT.getVertex(i++).e =  0.1656; // HE1
//					kDT.getVertex(i++).e =  0.1656; // HE2
//					kDT.getVertex(i++).e =  0.3992; // HH
//					break;
//				case 'V':
//					kDT.getVertex(i++).e = -0.4157; // N
//					kDT.getVertex(i++).e = -0.0875; // CA
//					kDT.getVertex(i++).e =  0.5973; // C
//					kDT.getVertex(i++).e = -0.5679; // O
//					kDT.getVertex(i++).e =  0.2985; // CB
//					kDT.getVertex(i++).e = -0.3192; // CG1
//					kDT.getVertex(i++).e = -0.3192; // CG2
//					kDT.getVertex(i++).e =  0.2719; // H
//					kDT.getVertex(i++).e =  0.0969; // HA
//					kDT.getVertex(i++).e = -0.0297; // HB
//					kDT.getVertex(i++).e =  0.0791; // HG11
//					kDT.getVertex(i++).e =  0.0791; // HG12
//					kDT.getVertex(i++).e =  0.0791; // HG13
//					kDT.getVertex(i++).e =  0.0791; // HG21
//					kDT.getVertex(i++).e =  0.0791; // HD22
//					kDT.getVertex(i++).e =  0.0791; // HD23
//					break;
//				default: 
//					System.out.println("Unidentified amino acid: " + aminoAcid);
//					break;
//			}
//		}
//	}
//	
//	
//	public boolean areApart(int a, int b, char aaTypeA, char aaTypeB) {
//		switch (a) {
//			case 1: return b != 0;
//			case 2: if ((b >= 0) && (b <= 1)) return false;
//					switch (aaTypeB) {
//						case 'A': case 'S': return b != 5;
//						case 'R': case 'F': return b != 11;
//						case 'N': case 'D': case 'I': case 'L': case 'M': return b != 8;
//						case 'C': case 'T': case 'Y':  return b != 7;
//						case 'Q': case 'E': case 'K': return b != 9;
//						case 'G': return b != 4;
//						case 'H': return b != 10;
//						case 'P': return b != 6;
//						case 'W': return b != 12;
//						default: return true;
//					}
//			default: return true;
//			}			
//	}
//					
//
//	
//	public boolean areApart3(int a, int b, char aaType) {
//		switch (aaType) {
//			case 'A': // ALA
//				switch (a) {
//					case 0: return b==3||b>=7&&b<=9;
//					case 1: case 7: case 8: case 9: return false;
//					case 2: return b==5||b>=7&&b<=9;
//					case 3: case 5: case 6: return true;
//					case 4: return b==5;
//			}
//			break;			
//			case 'R': // ARG
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=10||b>=13;
//					case 1: return b>=6&&b<=10||b>=15;
//					case 2: return b>=5&&b<=11||b>=13;
//					case 3: case 11: case 12: return true;
//					case 4: return b>=7&&b<=11||b>=17;
//					case 5: return b>=8&&b<=12||b>=19;
//					case 6: return b>=9&&b<=14||b>=20;
//					case 7: return b>=11&&b<=16||b>=20;
//					case 8: return b>=11&&b<=18;
//					case 9: return b>=11&&b<=19||b>=22&&b<=23;
//					case 10: return b>=11&&b<=19||b>=20&&b<=21;
//					case 13: case 14: return b>=15;
//					case 15: case 16: return b>=17;
//					case 17: case 18: return b>=19;
//					case 19: return b>=20;
//					case 20: case 21: return b>=22;
//					case 22: case 23: return false;
//			}
//			break;
//			case 'N': // ASN
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=7||b>=10;
//					case 1: return b>=6&&b<=7||b>=12;
//					case 2: return b>=5&&b<=8||b>=10;
//					case 3: case 8: case 9: return true;
//					case 4: return b==8||b>=12;
//					case 5: return b>=8&&b<=9;
//					case 6: return b>=8&&b<=9||b>=10;
//					case 7: return b>=8&&b<=11;
//					case 10: case 11: return b >=12;
//					case 12: case 13: return false;
//				}
//			break;
//			case 'D': // ASP
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=7||b>=10;
//					case 1: return b>=6&&b<=7;
//					case 2: return b>=5&&b<=8||b>=10;
//					case 3: case 7: case 8: case 9: return true;
//					case 4: return b==8;
//					case 6: return b>=8;
//					case 10: case 11: return false;
//				}
//			break;
//			case 'C': // CYS
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=6||b>=9;
//					case 1: return b==6;
//					case 2: return b>=5&&b<=7||b>=9;
//					case 3: case 6: case 7: case 8: return true;
//					case 4: return b==7;
//					case 5: return b>=7&&b<=8;
//					case 9: case 10: return false;
//				}
//			break;
//			case 'Q': // GLN
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=8||b>=11;
//					case 1: return b>=6&&b<=8||b>=13;
//					case 2: return b>=5&&b<=9||b>=11;
//					case 3: case 9: case 10: return true;
//					case 4: return b>=7&&b<=9||b>=15;
//					case 5: return b>=9&&b<=10||b>=15;
//					case 6: return b>=9&&b<=12;
//					case 7: return b>=9;
//					case 8: return b>=9&&b<=14;
//					case 11: case 12: return b >=13;
//					case 13: case 14: return b>=15;
//					case 15: case 16: return false;
//				}
//			break;
//			case 'E': // GLU
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=8||b>=11;
//					case 1: return b>=6&&b<=8||b>=13;
//					case 2: return b>=5&&b<=9||b>=11;
//					case 3: case 9: case 10: return true;
//					case 4: return b>=7&&b<=9;
//					case 5: return b>=9&&b<=10;
//					case 6: return b>=9&&b<=12;
//					case 7: case 8: return b>=9;
//					case 11: case 12: return b >=13;
//					case 13: case 14: return false;
//				}
//				break;
//			case 'G': // GLY
//				switch (a) {
//					case 0: return b==3;
//					case 1: case 5: case 6: return false;
//					case 2: return b==4;
//					case 3: return true;
//			}
//			break;			
//			case 'H': // HIS
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=9||b>=12;
//					case 1: return b>=6&&b<=9||b>=14;
//					case 2: return b>=5&&b<=11||b>=12;
//					case 3: case 10: case 11: case 13: case 14: case 15: case 16: case 17: return true;
//					case 4: return b>=8&&b<=10||b>=14;
//					case 5: return b>=10&&b<=11||b>=16;
//					case 6: return b>=10&&b<=13||b==15||b==17;
//					case 7: return b>=10&&b<=14||b==16;
//					case 8: return b>=10&&b<=13||b==15;
//					case 9: return b>=10&&b<=14||b>=22&&b<=23;
//					case 12: return b>=14;
//				}
//			break;
//			case 'I': // ILE
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=7||b>=10;
//					case 1: return b==7||b>=11;
//					case 2: return b>=5&&b<=8||b>=10;
//					case 3: case 8: case 9: case 10: case 12: case 15: return true;
//					case 4: return b==8||b>=16;
//					case 5: return b>=8&&b<=9||b>=13&&b<=15;
//					case 6: return b>=7&&b<=9||b>=11&&b<=12||b>=16;
//					case 7: return b>=8&&b<=10||b>=13&&b<=15;
//					case 11: return b>=13;
//					case 13: case 14: return b>=16;
//					case 16: case 17: case 18: return false;
//				}
//			break;
//			case 'L': // LEU
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=7||b>=10;
//					case 1: return b>=6&&b<=7||b>=12;
//					case 2: return b>=5&&b<=8||b>=10;
//					case 3: case 8: case 9: case 11: case 12: case 15: return true;
//					case 4: return b==8||b>=13;
//					case 5: return b>=8&&b<=9||b>=16;
//					case 6: return b>=8&&b<=11||b>=16;
//					case 7: return b>=8&&b<=11||b>=13&&b<=15;
//					case 10: return b>=12;
//					case 13: case 14: return b>=16;
//					case 16: case 17: case 18: return false;
//			}
//			break;
//			case 'K': // LYS
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=8||b>=11;
//					case 1: return b>=6&&b<=8||b>=13;
//					case 2: return b>=5&&b<=9||b>=11;
//					case 3: case 9: case 10: return true;
//					case 4: return b>=7&&b<=9||b>=15;
//					case 5: return b>=8&&b<=10||b>=17;
//					case 6: return b>=9&&b<=12||b>=19;
//					case 7: return b>=9&&b<=14||b>=20;
//					case 8: return b>=9&&b<=16;
//					case 11: case 12: return b>=13;
//					case 13: case 14: return b>=15;
//					case 15: case 16: return b>=17;
//					case 17: case 18: return b>=19;
//					case 19: case 20: case 21: return false;
//			}
//			break;
//			case 'M': // MET
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=7||b>=10;
//					case 1: return b>=6&&b<=7||b>=12;
//					case 2: return b>=5&&b<=8||b>=12;
//					case 3: case 9: case 11: case 13: return true;
//					case 4: return b>=7&&b<=8||b>=14;
//					case 5: return b>=8&&b<=9||b>=14;
//					case 6: return b>=8&&b<=11;
//					case 7: return b>=8&&b<=13;
//					case 8: return b>=10;
//					case 10: return b>=12;
//					case 12: return b>=14;
//					case 14: case 15: case 16: return false;
//			}
//			break;
//			case 'F': // PHE
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=10||b>=13;
//					case 1: return b>=6&&b<=10||b>=15;
//					case 2: return b>=5&&b<=11||b>=13;
//					case 3: case 12: case 14: case 15: case 16: case 17: case 18: return true;
//					case 4: return b>=8&&b<=11||b>=15;
//					case 5: return b>=8&&b<=12||b>=17;
//					case 6: return b==9||b>=11&&b<=14||b==16||b>=18;
//					case 7: return b==8||b>=11&&b<=15||b==17||b==19;
//					case 8: return b>=11&&b<=14||b==16||b==18;
//					case 9: return b>=11&&b<=15||b==17;
//					case 10: return b>=10&&b<=16;
//					case 11: return b>=13;
//					case 13: return b>=15;
//					case 19: return false;
//			}
//			break;
//			case 'P': // PRO
//				switch (a) {
//					case 0: return b==3||b>=8&&b<=11;
//					case 1: return b>=10;
//					case 2: return b>=5&&b<=6||b>=8;
//					case 3: case 9: case 11: return true;
//					case 4: return b>=12;
//					case 5: return b==7;
//					case 6: return b>=7&&b<=9;
//					case 7: return b>=8&b<=13;
//					case 8: return b>=10&&b<=13;
//					case 10: return b >=12;
//					case 12: case 13: return false;
//			}
//			break;
//			case 'S': // SER
//				switch (a) {
//					case 0: return b==3||b>=8;
//					case 1: return b==5||b==10;
//					case 2: return b>=5&&b<=6||b>=7;
//					case 3: case 7: case 9: return true;
//					case 4: return b==6;
//					case 5: return b>=6&&b<=7;
//					case 6: return b==7||b==10;
//					case 8: return b==10;
//					case 10: return false;
//			}
//			break;			
//			case 'T': // THR
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=6||b>=9;
//					case 1: return b>=5&&b<=6||b>=11;
//					case 2: return b>=5&&b<=7||b>=9;
//					case 3: case 8: case 10: case 12: return true;
//					case 4: return b==7||b==13;
//					case 5: return b>=7||b<=8;
//					case 6: return b>=7&&b<=10;
//					case 7: return b>=9;
//					case 9: return b>=11;
//					case 11: return b==13;
//					case 13: return false;
//				}
//			break;
//			case 'W': // TRP
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=13||b>=16;
//					case 1: return b>=6&&b<=10||b>=18;
//					case 2: return b>=5&&b<=14||b>=16;
//					case 3: case 14: case 15: case 17: case 18: case 19: case 20: case 21: case 22: return true;
//					case 4: return b>=8&&b<=14||b>=18;
//					case 5: return b>=11&&b<=15||b>=19;
//					case 6: return b>=10&&b<=17||b>=20;
//					case 7: return b==13&&b<=19||b>=21;
//					case 8: return b==10||b>=12&&b<=17||b==22;
//					case 9: return b==12||b>=14&&b<=18||b==22;
//					case 10: return b==11||b>=14&&b<=19||b==21||b==23;
//					case 11: return b>=14&&b<=20||b==22;
//					case 12: return b>=14&&b<=19||b==21;
//					case 13: return b>=14&&b<=20;
//					case 16: return b>=18;
//					case 23: return false;
//			}
//			break;
//			case 'Y': // TYR
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=11||b>=13;
//					case 1: return b>=10;
//					case 2: return b>=5&&b<=7||b>=9;
//					case 3: case 7: case 8: case 9: case 12: return true;
//					case 4: return b==7;
//					case 5: return b>=7&&b<=8||b>=13&&b<=15;
//					case 6: return b>=7&&b<=8||b>=10&&b<=12;
//					case 10: return b>=13;
//					case 11: return b>=13;
//					case 13: case 14: case 15: return false;
//			}
//			break;
//			case 'V': // VAL
//				switch (a) {
//					case 0: return b==3||b>=5&&b<=6||b>=9;
//					case 1: return b>=9;
//					case 2: return b>=5&&b<=7||b>=9;
//					case 3: case 7: case 8: case 9: case 12: case 14: return true;
//					case 4: return b==7;
//					case 5: return b>=7&&b<=8||b>=13;
//					case 6: return b>=7&&b<=8||b>=10&&b<=12;
//					case 10: case 11: return b>=13;
//					case 13: return b>=15;
//					case 15: return false;
//				}
//			break;
//
//			
//		}
//		return false;
//	}
//
//	
//	
//	
//	
//	
//	public double getNonBondedEnergy(KineticAlphaComplex kDT) { return getNonBondedEnergy(kDT, 10000000.0, 10000000.0);}
//
//	static StringBuilder interactionLog = new StringBuilder();
//	
//	public double[] areToClose(KineticAlphaComplex kDT) {
//		double[] energy = new double[2];
//		energy[0] = 0.0;
//		energy[1] = 0.0;
//		double sqDist, rij2, rij6;
//		int n = kDT.getNrVertices();
//		int i = 4;
//		while (i < n+3) {
//			System.out.println(i);
//			u = kDT.getVertex(i);
//			int nrAtoms = AminoAcid.getSize(kDT.getVertex(i).aaType);
//			for  (int j = i; j < i + nrAtoms; j++) {
//				for (int k = j+1; k < i + nrAtoms; k++) {
//					if (!areApart3(j, k, u.aaType)) {
//						v = kDT.getVertex(k);
//						sqDist = u.distanceSquared(v);
//						rij2 = r[u.atomType][v.atomType]/sqDist;
//						rij6 = rij2*rij2*rij2;
//						energy[0] += e[u.atomType][v.atomType]*rij6*(rij6 - 2);
//						energy[1] += u.e*v.e/Math.sqrt(sqDist);
//
//					}
//				}
//			}
//			if (i + nrAtoms < n+3) {
//				if (i == 1460) 
//					System.out.println(i + " " + nrAtoms + " " + (n+3));
//				int nrAtomsNext = AminoAcid.getSize(kDT.getVertex(i + nrAtoms).aaType);
//				for  (int j = i; j < i + nrAtoms; j++) {					
//					for (int k = i+nrAtoms; k < i + nrAtoms+nrAtomsNext; k++) {
//						if (k == 1486) 
//							System.out.println("STOP");
//						v = kDT.getVertex(k);
//						if (!areApart(j, k, u.aaType, v.aaType)) {
//							sqDist = u.distanceSquared(v);
//							rij2 = r[u.atomType][v.atomType]/sqDist;
//							rij6 = rij2*rij2*rij2;
//							energy[0] += e[u.atomType][v.atomType]*rij6*(rij6 - 2);
//							energy[1] += u.e*v.e/Math.sqrt(sqDist);
//						}
//					}
//				}
//			}
//			i = i + nrAtoms;
//		}
//		return energy;
//	}
//	
//	/** returns non-bonded potential between all atom within specified cutoff distances */
//	public double[] getNBP(KineticAlphaComplex kDT, double cutoffVdW, double cutoffElS) {
//		double[] energy = new double[2];
//		double energyVdW = 0.0;
//		double energyElS = 0.0;
//		
//		double sqDist, rij2,rij6;
//		double cutoffVdW2 = cutoffVdW*cutoffVdW;
//		double cutoffElS2 = cutoffElS*cutoffElS;
//		int uAtomType;
//
//		int n = kDT.getNrVertices();
//		for (int i = 4; i < n+3; i++) {
//			u = kDT.getVertex(i);
//			uAtomType = u.atomType;
//			for (int j = i+1; j < n+4; j++) {
//				v = kDT.getVertex(j);
//				sqDist = u.distanceSquared(v);
//				if (sqDist < cutoffVdW2) {
//					rij2 = r[uAtomType][v.atomType]/sqDist;
//					rij6 = rij2*rij2*rij2;
//					energyVdW += e[uAtomType][v.atomType]*rij6*(rij6 - 2);
//				}
//				if (sqDist < cutoffElS2) energyElS += u.e*v.e/Math.sqrt(sqDist);
//			}
//		}		
//		System.out.println("energyVdW = " + energyVdW + ", energyElS = " + energyElS);
//		energy[0] = energyVdW;
//		energy[1] = 332*energyElS;
//		return energy;
//	}
//
//	/** return non-bonded potential between all atoms using alpha complexes */
//	public double[] getNBP(KineticAlphaComplex kDT, int depth) {
//		double[] energy = new double[2];
//		double energyVdW = 0.0;
//		double energyElS = 0.0;
//		int uAtomType;
//		double sqDist, rij2,rij6;
//		Graph ac = kDT.getGraph();
//		ArrayList<Integer> neighbours;
//
//		for (int i = 4; i < kDT.getNrVertices(); i++) {
//			u = kDT.getVertex(i);
//			uAtomType = u.atomType;
//			neighbours = ac.breadthFirst(i, depth);
//			for (int j : neighbours) {
//				if (i < j) {
//					v = kDT.getVertex(j);
//					sqDist = u.distanceSquared(v);
//					rij2 = r[uAtomType][v.atomType]/sqDist;
//					rij6 = rij2*rij2*rij2;
//					energyVdW += e[uAtomType][v.atomType]*rij6*(rij6 - 2);	
//					energyElS += u.e*v.e/Math.sqrt(sqDist);
//				}
//			}
//		}		
//		energy[0] = energyVdW;
//		energy[1] = 332*energyElS;
//		return energy;
//	}
//	
//	
//	public double getNonBondedEnergy(KineticAlphaComplex kDT, double cutoffVdW, double cutoffElS) {
//		double cutoffVdW2 = cutoffVdW*cutoffVdW;
//		double cutoffElS2 = cutoffElS*cutoffElS;
//		double sqDist;
//		
//		double energyVdW = 0.0;
//		double energyElS = 0.0;
//		if (kDT.getInstanceType() == ProblemInstanceType.pdb) {
//			int n = kDT.getNrVertices();
//			int uAtomType;
//			double rij2, rij6;
//			int i = 4;
//			while (i < n+4) {
//				int nrAtoms = AminoAcid.getSize(kDT.getVertex(i).aaType);
//				int k = 0;
//				while (k < nrAtoms) {
//					u = kDT.getVertex(i+k);
//					uAtomType = u.atomType;
//					// interactions between atom i+k and other atoms in the same amino acid
//					for (int j = k+1; j < nrAtoms; j++) {
//						if (areApart3(k, j, u.aaType)) {
//							v = kDT.getVertex(i + j);
//							sqDist = u.distanceSquared(v);
//							if (sqDist < cutoffVdW2) {
//								rij2 = r[uAtomType][v.atomType]/sqDist;
//								rij6 = rij2*rij2*rij2;
//								energyVdW += e[uAtomType][v.atomType]*rij6*(rij6 - 2);
//							}
//							if (sqDist < cutoffElS2) {
//								energyElS += u.e*v.e/Math.sqrt(sqDist);
////								if( u.e*v.e/Math.sqrt(sqDist)>0.01 )
////									IOToolbox.writeToFile(u.atomName+"_"+u.atomId+" "+v.atomName+"_"+v.atomId+" "+ u.e*v.e/Math.sqrt(sqDist)+'\n', "/Users/pawel/interactionLog.txt", true);
//							}
//						}
//					}
//					// interactions between atom i+k and atoms in the next amino acid
//					if (i + nrAtoms < n + 4) {
//						int nrAtomsNext = AminoAcid.getSize(kDT.getVertex(i+nrAtoms).aaType);
//						for (int j = 0; j < nrAtomsNext; j++) {
//							Vertex v = kDT.getVertex(i + nrAtoms + j);
//							if (areApart(k, j, u.aaType, v.aaType)) {
//								sqDist = u.distanceSquared(v);
//								if (sqDist < cutoffVdW2) {
//									rij2 = r[uAtomType][v.atomType]/sqDist;
//									rij6 = rij2*rij2*rij2;
//									energyVdW += e[uAtomType][v.atomType]*rij6*(rij6 - 2);
//								}
//								if (sqDist < cutoffElS2) {
//									energyElS += u.e*v.e/Math.sqrt(sqDist);
////									if( u.e*v.e/Math.sqrt(sqDist)>0.01 )
////										IOToolbox.writeToFile(u.atomName+"_"+u.atomId+" "+v.atomName+"_"+v.atomId+" "+ u.e*v.e/Math.sqrt(sqDist)+'\n', "/Users/pawel/interactionLog.txt", true);
//								}
//							}
//						}
//						// interactions between atom i+k and all other atoms  
//						int l = i + nrAtoms + nrAtomsNext;
//						while (l < n+4) {
//							nrAtomsNext = AminoAcid.getSize(kDT.getVertex(l).aaType);
//							for (int j = 0; j < nrAtomsNext; j++) {
//								Vertex v = kDT.getVertex(l+j);
//								sqDist = u.distanceSquared(v);
//								if (sqDist < cutoffVdW2) {
//									rij2 = r[uAtomType][v.atomType]/sqDist;
//									rij6 = rij2*rij2*rij2;
//									energyVdW += e[uAtomType][v.atomType]*rij6*(rij6 - 2);
//								}
//								if (sqDist < cutoffElS2){
//									energyElS += u.e*v.e/Math.sqrt(sqDist);
////									if( u.e*v.e/Math.sqrt(sqDist)>0.01 )
////										IOToolbox.writeToFile(u.atomName+"_"+u.atomId+" "+v.atomName+"_"+v.atomId+" "+ u.e*v.e/Math.sqrt(sqDist)+'\n', "/Users/pawel/interactionLog.txt", true);
//								}
//							}
//							l = l + nrAtomsNext;	
//						}
//					}
//					k++;
//				}
//				i = i + nrAtoms;
//			}
//		}
//		System.out.println(energyVdW + "   " + (332*energyElS));
//		return energyVdW + 332*energyElS;
//	}
///*
//	public double getNonBondedEnergy(KineticDelaunayTessellation kDT, Graph G) {
//		double sqDist;
//		
//		double energyVdW = 0.0;
//		double energyElS = 0.0;
//		if (kDT.getInstanceType() == ProblemInstanceType.pdb) {
//			int n = kDT.getNrVertices();
//			int uAtomType;
//			double rij2, rij6;
//			int i = 4;
//			while (i < n) {
//				int nrAtoms = AminoAcid.getSize(kDT.getVertex(i).aaType);
//				int k = 0;
//				while (k < nrAtoms) {
//					u = kDT.getVertex(i+k);
//					uAtomType = u.atomType;
//					for (int j = i + nrAtoms; j < n; j++) {
//						v = kDT.getVertex(j);
//						sqDist = u.distanceSquared(v);
//						if (sqDist < cutoffVdW2) {
//							rij2 = r[uAtomType][v.atomType]/sqDist;
//							rij6 = rij2*rij2*rij2;
////							if (e[uAtomType][v.atomType]*rij6*(rij6 - 2) > 10) 
////								System.out.println("jump");
//							energyVdW += e[uAtomType][v.atomType]*rij6*(rij6 - 2);
//						}
//						if (sqDist < cutoffElS2) energyElS += u.e*v.e/Math.sqrt(sqDist);
////						System.out.println((i+k) + "," + j + ": " + energyVdW + "   " + energyElS);
//					}
//					k++;
//				}
//				i = i + nrAtoms;
//			}
//			
//		}
//		System.out.println(energyVdW + "   " + (332*energyElS));
//		return energyVdW + 332*energyElS;
//	}
//*/
//	
//	public double getVdWPotentialDTVeryFast(KineticAlphaComplex kDT, int maxDepth) {
//		int n = kDT.getNrVertices();
//		double rmj2, rmj6;
//		int m, j, mType, jType;
//		Tet tet;
//		HashSet<Vertex> processed = new HashSet<Vertex>();
//		ArrayList<Vertex> reachedVertices = new ArrayList<Vertex>();
//		ArrayList<Vertex> tetVertices = new ArrayList<Vertex>();
//		for (int i = 0; i < 4; i++) {
//			kDT.getVertex(i).setDepth(-1);
//			processed.add(kDT.getVertex(i));
//		}
//		double energy = 0.0;
//		
//		for (int i = 4; i < n; i++) {
//			tet = kDT.getVertex(i).getTet();
//			tetVertices.clear();
//			for (int k = 0; k < 4; k++) if (!processed.contains(tet.getCorner(k))) tetVertices.add(tet.getCorner(k));
//			if (!tetVertices.isEmpty()) {
//				reachedVertices = tet.breadthFirstVertices(maxDepth);
//				for (Vertex v0 : tetVertices) {
//					processed.add(v0);
//					m = v0.getId();			
//					mType = (m-4)%3;
//					for (Vertex v : reachedVertices) {
//						j = v.getId();
//						if (m+4 <= j) {
//							jType = (j-4)%3;
//							rmj2 = r[mType][jType]/v0.distanceSquared(v);
//							rmj6 = rmj2*rmj2*rmj2;
//							energy += e[mType][jType]*rmj6*(rmj6 - 2);
//						}
//					}
//				}
//				for (Vertex v : reachedVertices) v.flag = false;
//			}
//		}
//		return energy;
//	}
//		
//	public double getVdWPotentialDTFast(KineticAlphaComplex kDT) {
//		double energyVdW = 0.0;
//		double energyElS = 0.0;
//				
//		int n = kDT.getNrVertices();
//		double rij2, rij6;
//		int iType, jType;
//		int i, j;
//		int v0Neighbors = 0;
//		int v1Neighbors = 0;
//		int v2Neighbors = 0;
//		
//		for (i = 4; i < n-4; i++) {
//			Vertex v = kDT.getVertex(i);		
//			if (kDT.testingScreen && (i == 4)) v.toScene(kDT.getScene(), 0.1, Color.blue);
//			ArrayList<Vertex> adjList = v.computeAdjVertices(kDT.getScene());									
//			
//			int indx1 = adjList.size();														
//			for (int k = 0; k < indx1; k++) {
//				for (Vertex v1 : adjList.get(k).computeAdjVertices(kDT.getScene())) 
//					if ((v1 != v) && !adjList.contains(v1)) adjList.add(v1); 
//			}
//			
//			int indx2 = adjList.size();
//			for (int k = indx1; k < indx2; k++) {
//				for (Vertex v2 : adjList.get(k).computeAdjVertices(kDT.getScene())) 
//					if ((v2 != v) && !adjList.contains(v2)) adjList.add(v2); 
//				
//			}
//			iType = (i-4)%3;
//			double sqDist;
//			for (Vertex b : adjList) {
//				j = b.getId() - 4;
//				if (i <= j) {
//					jType = j%3;
//					sqDist = v.distanceSquared(b);
//					rij2 = r[iType][jType]/v.distanceSquared(b);
//					rij6 = rij2*rij2*rij2;
//					energyVdW += e[iType][jType]*rij6*(rij6 - 2);
//					energyElS += u.e*v.e/Math.sqrt(sqDist);
//				}
//			}
//			for (Vertex u : adjList) u.flag = false;
//			for (Tet tet : kDT.getTetrahedra()) tet.setFlag(false);
//
//		}
//		System.out.println(" Averege number neighbors of level 1 vertices: " + v0Neighbors/(n-4));
//		System.out.println(" Averege number neighbors of level 2 vertices: " + v1Neighbors/(n-4));
//		System.out.println(" Averege number neighbors of level 3 vertices: " + v2Neighbors/(n-4));
//
//
//		
//		System.out.println(energyVdW + "   " + (332*energyElS));
//		return energyVdW + 332*energyElS;
//
//	}
//
//	/* test potential energy computations using a-complexes with a in [1.0, 6.0[ and depths d in [1, 5] */
//	private void testing(String name) throws IOException {
//		double[] kDTTime = new double[50];
//		double[][] acTime = new double[50][5];
//		double[][] energyACVdW = new double[50][5];
//		double[][] energyACElS = new double[50][5];
//		long[] cutoffTime = new long[13];
//		DecimalFormat df = new DecimalFormat("#.00");
//
//		VanDerWaals vDW = new VanDerWaals();
//		KineticAlphaComplex kDT = new KineticAlphaComplex(1.5, ProblemInstanceType.pdb);
//		vDW.setUpCharges(kDT);
//
//		double energyCorrection[] = areToClose(kDT);
//		System.out.println("Energy correction: VdW: " + df.format(energyCorrection[0]) + ", ElS: " + df.format(energyCorrection[1]));
//		long start = System.nanoTime();
//		double[] energyCutOff = vDW.getNBP(kDT, 100000000, 10000000);
//		cutoffTime[0] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff100 = vDW.getNBP(kDT, 100,100);
//		cutoffTime[1] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff50 = vDW.getNBP(kDT, 50,50);
//		cutoffTime[2] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff30 = vDW.getNBP(kDT, 30,30);
//		cutoffTime[3] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff20 = vDW.getNBP(kDT, 20,20);
//		cutoffTime[4] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff15 = vDW.getNBP(kDT, 15,15);
//		cutoffTime[5] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff12 = vDW.getNBP(kDT, 12,12);
//		cutoffTime[6] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff10 = vDW.getNBP(kDT, 10,10);
//		cutoffTime[7] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff8 = vDW.getNBP(kDT, 8,8);
//		cutoffTime[8] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff6 = vDW.getNBP(kDT, 6,6);
//		cutoffTime[9] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff4 = vDW.getNBP(kDT, 4,4);
//		cutoffTime[10] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff2 = vDW.getNBP(kDT, 2,2);
//		cutoffTime[11] = System.nanoTime() - start;
//		start = System.nanoTime();
//		double[] energyCutOff1 = vDW.getNBP(kDT, 1,1);
//		cutoffTime[12] = System.nanoTime() - start;
//		energyCutOff[0] -= energyCorrection[0];
//		energyCutOff[1] -= energyCorrection[1];
//
//		System.out.println("Energy cutoff VdW: " + df.format(energyCutOff[0]) + ", ElS: " + df.format(energyCutOff[1]));
//		
//		
//		for (int i = 0; i < 50; i++) {
//			System.out.println("alpha = " + (1.0+0.1*i));
//			start = System.nanoTime();
//			kDT = new KineticAlphaComplex(1.0+ 0.1*i, ProblemInstanceType.pdb);
//			vDW.setUpCharges(kDT);
//			long slut = System.nanoTime();
//			kDTTime[i] = (slut-start)/1000000.0;
//			
//			for (int depth = 1; depth < 6; depth++) {
//				System.out.print("      depth = " + depth + ": ");
//				start = System.nanoTime();
//				double energyAC[] = vDW.getNBP(kDT, depth);
//				slut = System.nanoTime();
//				acTime[i][depth-1] = (slut - start)/1000000.0;
//				energyACVdW[i][depth-1] = energyAC[0] - energyCorrection[0];
//				energyACElS[i][depth-1] = energyAC[1] - energyCorrection[1];
//
//				System.out.println("energyVdW = " + energyACVdW[i][depth-1] + ", energyElS = " + energyACElS[i][depth-1]);
//
//			}
//			
//			
//		}
//		
//		try {
//			String filename = "/Users/pawel/" + name + ".xls";
//			WritableWorkbook workbook = Workbook.createWorkbook(new File(filename));
//			WritableSheet sheet1 = workbook.createSheet("Sheet1", 0);
//			WritableSheet sheet2 = workbook.createSheet("Sheet2", 0);
//			WritableSheet sheet3 = workbook.createSheet("Sheet3", 0);
//			WritableSheet sheet4 = workbook.createSheet("Sheet4", 0);
//			
//			//adding labels
//			sheet1.addCell(new Label(0, 0, name + " VdW"));
//			sheet2.addCell(new Label(0, 0, name + " ElS"));
//			sheet3.addCell(new Label(0, 0, name + " time"));
//			sheet4.addCell(new Label(0, 0, name + " cutoff"));
//			
//			for (int i1 = 1; i1 < 6; i1++) {
//				sheet1.addCell(new Label(i1, 1, "d = " + i1));
//				sheet2.addCell(new Label(i1, 1, "d = " + i1));
//				sheet3.addCell(new Label(i1, 1, "d = " + i1));
//			}
//			sheet4.addCell(new Number(2, 1, 100000000));
//			sheet4.addCell(new Number(3, 1, 100));
//			sheet4.addCell(new Number(4, 1, 50));
//			sheet4.addCell(new Number(5, 1, 30));
//			sheet4.addCell(new Number(6, 1, 20));
//			sheet4.addCell(new Number(7, 1, 15));
//			sheet4.addCell(new Number(8, 1, 12));			
//			sheet4.addCell(new Number(9, 1, 10));			
//			sheet4.addCell(new Number(10, 1, 8));			
//			sheet4.addCell(new Number(11, 1, 6));
//			sheet4.addCell(new Number(12, 1, 4));			
//			sheet4.addCell(new Number(13, 1, 2));			
//			sheet4.addCell(new Number(14, 1, 1));			
//			sheet4.addCell(new Number(2, 3, energyCutOff[0]));
//			sheet4.addCell(new Number(2, 4, energyCutOff[1]));
//			sheet4.addCell(new Number(3, 3, energyCutOff100[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(3, 4, energyCutOff100[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(3, 5, Math.abs((energyCutOff100[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(3, 6, Math.abs((energyCutOff100[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(4, 3, energyCutOff50[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(4, 4, energyCutOff50[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(4, 5, Math.abs((energyCutOff50[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(4, 6, Math.abs((energyCutOff50[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(5, 3, energyCutOff30[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(5, 4, energyCutOff30[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(5, 5, Math.abs((energyCutOff30[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(5, 6, Math.abs((energyCutOff30[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(6, 3, energyCutOff20[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(6, 4, energyCutOff20[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(6, 5, Math.abs((energyCutOff20[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(6, 6, Math.abs((energyCutOff20[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(7, 3, energyCutOff15[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(7, 4, energyCutOff15[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(7, 5, Math.abs((energyCutOff15[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(7, 6, Math.abs((energyCutOff15[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(8, 3, energyCutOff12[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(8, 4, energyCutOff12[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(8, 5, Math.abs((energyCutOff12[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(8, 6, Math.abs((energyCutOff12[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(9, 3, energyCutOff10[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(9, 4, energyCutOff10[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(9, 5, Math.abs((energyCutOff10[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(9, 6, Math.abs((energyCutOff10[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(10, 3, energyCutOff8[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(10, 4, energyCutOff8[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(10, 5, Math.abs((energyCutOff8[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(10, 6, Math.abs((energyCutOff8[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(11, 3, energyCutOff6[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(11, 4, energyCutOff6[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(11, 5, Math.abs((energyCutOff6[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(11, 6, Math.abs((energyCutOff6[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(12, 3, energyCutOff4[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(12, 4, energyCutOff4[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(12, 5, Math.abs((energyCutOff4[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(12, 6, Math.abs((energyCutOff4[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(13, 3, energyCutOff2[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(13, 4, energyCutOff2[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(13, 5, Math.abs((energyCutOff2[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(13, 6, Math.abs((energyCutOff2[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//			sheet4.addCell(new Number(14, 3, energyCutOff1[0] - energyCorrection[0]));
//			sheet4.addCell(new Number(14, 4, energyCutOff1[1] - energyCorrection[1]));
//			sheet4.addCell(new Number(14, 5, Math.abs((energyCutOff1[0]-energyCorrection[0]-energyCutOff[0])/energyCutOff[0])));
//			sheet4.addCell(new Number(14, 6, Math.abs((energyCutOff1[1]-energyCorrection[1]-energyCutOff[1])/energyCutOff[1])));
//
//			for (int i = 2; i < 15; i++) sheet4.addCell(new Number(i, 7, cutoffTime[i-2]/1000000.0));
//
//			
//			
//			for (int i1 = 0; i1 < 50; i1++) {
//				sheet1.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				sheet2.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				sheet3.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				for (int j = 0; j < 5; j++) {
//					sheet1.addCell(new Number(j+1, i1+2, Math.abs(energyACVdW[i1][j] - energyCutOff[0])/Math.abs(energyCutOff[0])));
//					sheet2.addCell(new Number(j+1, i1+2, Math.abs(energyACElS[i1][j] - energyCutOff[1])/Math.abs(energyCutOff[1])));
//					sheet3.addCell(new Number(j+1, i1+2, acTime[i1][j]));
//				}
//			}
//			
//			workbook.write();
//			workbook.close();
//		} catch (WriteException e) {}
//		System.out.println("Done");
//
//	}
//
//	/* test the convergence of the potential energy computations using a-complexes with a in [1.0, 6.0[ and depths d in [1, 5] */
//	private void testingConvergence(String name, int factor) throws IOException {
//		double[] kDTTime = new double[50];
//		double[][] acTime = new double[50][5];
//		double[][] energyACVdW = new double[50][5];
//		double[][] energyACElS = new double[50][5];
//		DecimalFormat df = new DecimalFormat("#.00");
//
//		VanDerWaals vDW = new VanDerWaals();
//		KineticAlphaComplex kDT = new KineticAlphaComplex(1.5, ProblemInstanceType.pdb);
//		vDW.setUpCharges(kDT);
//
//		double energyCorrection[] = areToClose(kDT);
//		System.out.println("Energy correction: VdW: " + df.format(energyCorrection[0]) + ", ElS: " + df.format(energyCorrection[1]));
//		
//		double[] energyCutOff = vDW.getNBP(kDT, 100000000, 10000000);
//		energyCutOff[0] -= energyCorrection[0];
//		energyCutOff[1] -= energyCorrection[1];
//
//		System.out.println("Energy cutoff VdW: " + df.format(energyCutOff[0]) + ", ElS: " + df.format(energyCutOff[1]));
//		
//		
//		for (int i = 0; i < 10; i++) {
//			System.out.println("alpha = " + (1.0+0.1*i));
//			long start = System.nanoTime();
//			kDT = new KineticAlphaComplex(1.0+ 0.1*i, ProblemInstanceType.pdb);
//			vDW.setUpCharges(kDT);
//			long slut = System.nanoTime();
//			kDTTime[i] = (slut-start)/1000000.0;
//			
//			for (int depth = 1; depth < 6; depth++) {
//				System.out.print("      depth = " + factor*depth + ": ");
//				start = System.nanoTime();
//				double energyAC[] = vDW.getNBP(kDT, factor*depth);
//				slut = System.nanoTime();
//				acTime[i][depth-1] = (slut - start)/1000000.0;
//				energyACVdW[i][depth-1] = energyAC[0] - energyCorrection[0];
//				energyACElS[i][depth-1] = energyAC[1] - energyCorrection[1];
//
//				System.out.println("energyVdW = " + energyACVdW[i][depth-1] + ", energyElS = " + energyACElS[i][depth-1]);
//
//			}		
//		}
//		
//		try {
//			String filename = "/Users/pawel/" + name + ".conv";
//			WritableWorkbook workbook = Workbook.createWorkbook(new File(filename));
//			WritableSheet sheet1 = workbook.createSheet("Sheet1", 0);
//			WritableSheet sheet2 = workbook.createSheet("Sheet2", 0);
//			WritableSheet sheet3 = workbook.createSheet("Sheet3", 0);
//			
//			//adding labels
//			sheet1.addCell(new Label(0, 0, name + " VdW"));
//			sheet2.addCell(new Label(0, 0, name + " ElS"));
//			sheet3.addCell(new Label(0, 0, name + " time"));
//			
//			for (int i1 = 1; i1 < 6; i1++) {
//				sheet1.addCell(new Label(i1, 1, "d = " + factor*i1));
//				sheet2.addCell(new Label(i1, 1, "d = " + factor*i1));
//				sheet3.addCell(new Label(i1, 1, "d = " + factor*i1));
//			}		
//			
//			for (int i1 = 0; i1 < 10; i1++) {
//				sheet1.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				sheet2.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				sheet3.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				for (int j = 0; j < 5; j++) {
//					sheet1.addCell(new Number(j+1, i1+2, Math.abs(energyACVdW[i1][j] - energyCutOff[0])/Math.abs(energyCutOff[0])));
//					sheet2.addCell(new Number(j+1, i1+2, Math.abs(energyACElS[i1][j] - energyCutOff[1])/Math.abs(energyCutOff[1])));
//					sheet3.addCell(new Number(j+1, i1+2, acTime[i1][j]));
//				}
//			}
//			
//			workbook.write();
//			workbook.close();
//		} catch (WriteException e) {}
//		System.out.println("Done");
//
//	}
//
//	
//	/* tests potential energy computations using a-complexes for a in [1.3, 1.5[ and depths d in [1, 5]*/
//	private void testingDetail(String name) throws IOException {
//		double[] kDTTime = new double[50];
//		double[][] acTime = new double[50][5];
//		double[][] energyACVdW = new double[50][5];
//		double[][] energyACElS = new double[50][5];
//		DecimalFormat df = new DecimalFormat("#.00");
//
//		VanDerWaals vDW = new VanDerWaals();
//		KineticAlphaComplex kDT = new KineticAlphaComplex(1.5, ProblemInstanceType.pdb);
//		vDW.setUpCharges(kDT);
//
//		double energyCorrection[] = areToClose(kDT);
//		System.out.println("Energy correction: VdW: " + df.format(energyCorrection[0]) + ", ElS: " + df.format(energyCorrection[1]));
//		
//		double[] energyCutOff = vDW.getNBP(kDT, 100000000, 10000000);
//		energyCutOff[0] -= energyCorrection[0];
//		energyCutOff[1] -= energyCorrection[1];
//
//		System.out.println("Energy cutoff VdW: " + df.format(energyCutOff[0]) + ", ElS: " + df.format(energyCutOff[1]));
//		
//		
//		for (int i = 0; i < 50; i++) {
//			System.out.println("alpha = " + (1+0.01*i));
//			long start = System.nanoTime();
//			kDT = new KineticAlphaComplex(1+ 0.01*i, ProblemInstanceType.pdb);
//			vDW.setUpCharges(kDT);
//			long slut = System.nanoTime();
//			kDTTime[i] = (slut-start)/1000000.0;
//			
//			for (int depth = 1; depth < 6; depth++) {
//				System.out.print("      depth = " + depth + ": ");
//				start = System.nanoTime();
//				double energyAC[] = vDW.getNBP(kDT, depth);
//				slut = System.nanoTime();
//				acTime[i][depth-1] = (slut - start)/1000000.0;
//				energyACVdW[i][depth-1] = energyAC[0] - energyCorrection[0];
//				energyACElS[i][depth-1] = energyAC[1] - energyCorrection[1];
//
//				System.out.println("energyVdW = " + energyACVdW[i][depth-1] + ", energyElS = " + energyACElS[i][depth-1]);
//
//			}
//			
//			
//		}
//		
//		try {
//			String filename = "/Users/pawel/" + name + "Detail.xls";
//			WritableWorkbook workbook = Workbook.createWorkbook(new File(filename));
//			WritableSheet sheet1 = workbook.createSheet("Sheet1", 0);
//			WritableSheet sheet2 = workbook.createSheet("Sheet2", 0);
//			WritableSheet sheet3 = workbook.createSheet("Sheet3", 0);
//			
//			//adding labels
//			sheet1.addCell(new Label(0, 0, name + " VdW"));
//			sheet2.addCell(new Label(0, 0, name + " ElS"));
//			sheet3.addCell(new Label(0, 0, name + " time"));
//			
//			for (int i1 = 1; i1 < 6; i1++) {
//				sheet1.addCell(new Label(i1, 1, "d = " + i1));
//				sheet2.addCell(new Label(i1, 1, "d = " + i1));
//				sheet3.addCell(new Label(i1, 1, "d = " + i1));
//			}
//			
//			for (int i1 = 0; i1 < 50; i1++) {
//				sheet1.addCell(new Label(0, i1+2, df.format(1+0.01*i1)));
//				sheet2.addCell(new Label(0, i1+2, df.format(1+0.01*i1)));
//				sheet3.addCell(new Label(0, i1+2, df.format(1+0.01*i1)));
//				for (int j = 0; j < 5; j++) {
//					sheet1.addCell(new Number(j+1, i1+2, Math.abs(energyACVdW[i1][j] - energyCutOff[0])/Math.abs(energyCutOff[0])));
//					sheet2.addCell(new Number(j+1, i1+2, Math.abs(energyACElS[i1][j] - energyCutOff[1])/Math.abs(energyCutOff[1])));
//					sheet3.addCell(new Number(j+1, i1+2, acTime[i1][j]));
//				}
//			}
//			
//			workbook.write();
//			workbook.close();
//		} catch (WriteException e) {}
//		System.out.println("Done");
//
//	}
//	
//
//	private void testingKinetic(String name) throws IOException {
//		KineticAlphaComplex kDT;
//		ArrayList<Integer> qVertices = new ArrayList<Integer>();
//		ArrayList<Integer> nVertices = new ArrayList<Integer>();
//		double [][] average = new double[50][5];
//		int counter;
//		Tet nTet;
//		double[] kDTTime = new double[50];
//		for (int i = 0; i < 50; i++) {
//			System.out.println("alpha = " + (1+0.1*i));
//			long start = System.nanoTime();
//			kDT = new KineticAlphaComplex(1+ 0.1*i, ProblemInstanceType.pdb);
//			long slut = System.nanoTime();
//			kDTTime[i] = (slut-start)/1000000.0;
//			
//			for (int depth = 1; depth < 6; depth++) {
//				counter = 0;
//				System.out.print("      depth = " + depth + ": ");
//				start = System.nanoTime();
//				for (Tet tet : kDT.getTets()) {
//					if (!tet.isBig()) {
//						qVertices.clear();
//						for (int k = 0; k < 4; k++) qVertices.add(tet.getCorner(k).getId());
//						for (int k = 0; k < 4; k++) {
//							nTet = tet.getNeighbor(k);
//							int q = nTet.getCorner(nTet.apex(tet)).getId();
//							if (tet.getCorner(k).getId() < q) {
//								qVertices.add(q);
//								nVertices = kDT.getGraph().breadthFirst(qVertices, depth);
//								average[i][depth-1] += nVertices.size();
//								counter++;
//								qVertices.remove(4);
//							}
//						}
//					}
//				}
//				average[i][depth-1] /= counter; 
//			}
//			System.out.println();
//		}
//		DecimalFormat df = new DecimalFormat("#.00");
//		try {
//			String filename = "/Users/pawel/" + name + "neighborhoods.xls";
//			WritableWorkbook workbook = Workbook.createWorkbook(new File(filename));				
//			WritableSheet sheet1 = workbook.createSheet("Sheet1", 0);
//				
//			//adding labels
//			sheet1.addCell(new Label(0, 0, name + " average"));
//				
//			for (int i1 = 1; i1 < 6; i1++) {
//				sheet1.addCell(new Label(i1, 1, "d = " + i1));
//			}
//				
//			for (int i1 = 0; i1 < 50; i1++) {
//				sheet1.addCell(new Label(0, i1+2, df.format(1+0.1*i1)));
//				for (int j = 0; j < 5; j++) {						
//					sheet1.addCell(new Number(j+1, i1+2, average[i1][j]));
//				}
//			}
//				
//			workbook.write();
//			workbook.close();
//		} catch (WriteException e) {}
//		System.out.println("Done");
//	}
//	
//	public static void main(String[] args) throws IOException{
//		
//		
//		
//		VanDerWaals vDW = new VanDerWaals();
//		vDW.testing("sample_000009600000_0_302.425321");
//
//		
//		
///*		vDW.setUpCharges(kDT);
//		
//		int cutoffVdW = 10;
//		int cutoffElS = 40;
//		for (int i = 1; i < 3; i++) { 	
//			start = System.nanoTime();
//			double[] energyCutoff = vDW.getNBP(kDT, cutoffVdW, cutoffElS);
//			slut = System.nanoTime();
//			System.out.printf(" Simple cutoff for VdW: " + cutoffVdW + ", and for ElS: " + cutoffElS + ", time in miliseconds: %.2f\n", (slut - start)/1000000.0);
//			System.out.println("Energy = " + (energyCutoff[0]+energyCutoff[1]));
//		}
//		for (int i = 1; i < 6; i++) {
//			start = System.nanoTime();
//			double[] energyAC = vDW.getNBP(kDT, i);
//			slut = System.nanoTime();
//			System.out.printf(kDT.getAlpha() + "-complex for VdW, depth = " + i + ", time in miliseconds: %.2f\n", (slut - start)/1000000.0);
//			System.out.println("Energy = " + (energyAC[0]+energyAC[1]));
//		}
//		
//		J3DScene scene  = J3DScene.createJ3DSceneInFrame();
//		kDT.toScene(scene, kDT.getAlpha());
////		kDT.toSceneComplement(scene, 2.8, 10);
//		
//		start = System.nanoTime();
//		double energyFastDT = vDW.getVdWPotentialDTFast(kDT);
//		System.out.printf(" alpha complex (alpha = " + alpha + "), time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
//		System.out.println("EnergyDTFast = " + energyFastDT);
//*/		
///*		start = System.nanoTime();
//		double energyVeryFastDT = vDW.getVdWPotentialDTVeryFast(kDT, 20);
//		System.out.printf(" alpha complex (alpha = " + alpha + "), time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
//		System.out.println("EnergyDTVeryFast = " + energyVeryFastDT);
//*/
//	}
	
	
}
