package correlescence;


import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.plugin.filter.MaximumFinder;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class Temporal_ICS implements PlugIn {
	
	
	/** original stack **/
	ImagePlus imp;
	
	/** correlation calculation method 
	 * 0 = direct correlation
	 * 1 = correlation using FFT **/
	int nCalcMethod=1;
	
	/** correlation normalization method 
	 * 0 = overlap area
	 * 1 = full area **/
	int nNormMethod=1;
	
	/** object calculating correlation **/
	imCC1D x1D = new imCC1D();
	
	/** maximum delay in frames ICS**/
	int nMaxUserDelayICS;
	
	/** input stack size **/
	int nStackSize;
	
	/** image width in px**/
	int origW;
	
	/** image height in px**/
	int origH;
	
	/** index of reslicing **/
	int nindW;
	
	/** subtract projection **/
	boolean bSubtractProj;
	
	/** stack projection to subtract **/
	int nProjMethod;
	
	/** find first maximum **/
	boolean bFindFirstMax;
	
	/** estimate the period/frequency more precisely **/
	boolean bSubFramePeak;
	
	/** maximum finding tolerance **/
	double dMaxTol;
	
	/** whether or not show final TICS stack **/
	boolean bShowTICSStack;
	
	/** rotated stack containing TICS **/
	ImageStack ipICSRotated;
	
	/** rescale CC image? **/
	boolean bRescaleCC;
	
	/** new scaled CC width  **/
	int nScaledW = 20;
	/** new scaled CC height  **/
	int nScaledH = 20;
	
	/** whether to calculate frequency/phase map **/
	boolean bICSfreqphase;
	
	/** image processor containing dominating frequency **/
	ImageProcessor freqip = null;
	
	/** image processor containing dominating phase **/
	ImageProcessor phaseip = null; 
	
	//dialog stuff
	String[] projMethods = 
		{"Average Intensity", "Max Intensity", "Min Intensity", "Sum Slices", "Standard Deviation", "Median"}; 
	
	String [] sOutput = new String [] {
			"Frequency","Period"};
	String [] sOutputUnits = new String [] {
			"Frames","Stack's time units"};
	
	double dTimeUnits = 1.0;
	
	int nOutput;
	int nOutputUnits;

	
	@Override
	public void run(String arg) 
	{
		ImagePlus finalImp;
		ImagePlus finalImpMax=null;
		String sTitle;
		
		// TODO Auto-generated method stub
		imp = IJ.getImage();		

		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}		
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32 bit greyscale image required");
		    return;
		}
		
		nStackSize = imp.getStackSize();
		origW=imp.getWidth();
		origH=imp.getHeight();
		if(nStackSize ==1)
		{
		    IJ.error("Stack of at least 2 images is required, aborted.");
		    	return;
		}
		
		if(!xDialogTICS())
			return;
		
		sTitle=imp.getTitle();
		Calibration outCal = (Calibration) imp.getCalibration().clone();
		if(nOutputUnits == 1)
		{
			dTimeUnits = outCal.frameInterval;
			if(Math.abs(dTimeUnits)<0.0000000001)
			{
				IJ.showMessage("Stack's frame interval is equal to zero, proceeding with output in frame units.");
				IJ.log("Stack's frame interval is equal to zero, proceeding with output in frame units.");
				dTimeUnits = 1.0;
				nOutputUnits = 0;
			}
		}
		
		IJ.log(" --- Correlescence plugin version " + ConstantsCorrelescence.sVersion+ " --- ");
		IJ.log("Running Temporal ICS analysis.");
		IJ.log("Image title: \"" + sTitle + "\"");
		
		if(bSubtractProj)
		{
			subtractProj();
			IJ.log("Subtracted "+projMethods[nProjMethod]);
		}
		

		
		int nHalfStack=(int)Math.floor((double)nStackSize*0.5);
		//maximum delay in time
		if(nMaxUserDelayICS==0)
		{
			nMaxUserDelayICS = nHalfStack;
			IJ.log("Maximum delay value : "+String.format(" %d frames", nMaxUserDelayICS));
		}
		else
		{
			if(nMaxUserDelayICS> nHalfStack)
			{			
				nMaxUserDelayICS=nHalfStack;
				IJ.log("Maximum delay value is too high, reducing it to "+String.format(" %d frames", nMaxUserDelayICS));
			}
		}
		

		
		//function to calculate TICS,
		// returns "rotated stack", where time instead of
		// becomes x
		//stores it in ipICSRotated
		ICS_pixel();
		
		

		
		if(bRescaleCC)
		{
			IJ.showStatus("Resizing CC...");
			IJ.log("Output is resized to: " + Integer.toString(nScaledW) + "x" + Integer.toString(nScaledH) + " pixels");
			
			ImagePlus ipRescale = new ImagePlus("TICS_"+sTitle, resliceRotateStack(new ImagePlus("rotateback",ipICSRotated)));
			int nProgressTotal = ipRescale.getImageStackSize() + nScaledW;
			int nProgressCount = 0;	
			ImageStack isRescale = new ImageStack(nScaledW, nScaledH);
			ImageProcessor singleImage;
			for (int i=0;i<ipRescale.getImageStackSize();i++)
			{
				ipRescale.setPosition(i+1);
				singleImage = ipRescale.getProcessor();
				singleImage.setInterpolationMethod(ImageProcessor.BILINEAR);
				isRescale.addSlice(singleImage.resize(nScaledW, nScaledH, true));
				IJ.showProgress(nProgressCount, nProgressTotal);
				nProgressCount++;
			}
			ipRescale = new ImagePlus("rescaled", isRescale);
			
			ipICSRotated = new ImageStack(nMaxUserDelayICS+1,nScaledH);
			for (nindW=0;nindW<nScaledW;nindW++)
			{
				ipICSRotated.addSlice(getResliceRotatedSlice(ipRescale, nindW));	
				IJ.showProgress(nProgressCount, nProgressTotal);
				nProgressCount++;
			}
			
			IJ.showStatus("Resizing CC...done.");
			IJ.showProgress(3, 3);
			outCal.pixelWidth *= origW/nScaledW;
			outCal.pixelHeight *= origH/nScaledH;
			//new ImagePlus("rotated again", ipICSRotated).show();
		}
		
		//rotate stack so time becomes Z again
		if(bShowTICSStack)
		{
			finalImp = new ImagePlus("TICS_"+sTitle, resliceRotateStack(new ImagePlus("rotateback",ipICSRotated)));
			finalImp.setCalibration(outCal);
			finalImp.show();
		}
		
		//if we want to calculate maximum
		if(bFindFirstMax)
		{
			String sTitleMax = "TICS_";
			if(nOutput==0)
			{
				sTitleMax = sTitleMax+"frequency_";
				IJ.log("Reported TICS results represent frequency.");
			}
			else
			{
				sTitleMax = sTitleMax+"period_";
				IJ.log("Reported TICS results represent period.");
			}
			if(nOutputUnits==0)
			{
				IJ.log("Time units: frames.");
			}
			else
			{
				IJ.log("Frame interval "+Double.toString(dTimeUnits)+ " "+outCal.getTimeUnit());
			}
				
			finalImpMax = new ImagePlus(sTitleMax+sTitle,findMax(ipICSRotated,dMaxTol));
			finalImpMax.setCalibration(outCal);	
			finalImpMax.show();
			IJ.run(finalImpMax,"Fire","");
		}
		if(bICSfreqphase)
		{
			ImagePlus finalFFTfr;
			ImagePlus finalFFTphase;
			if(nOutput==0)
			{
				finalFFTfr = new ImagePlus("Dominant_FFT_frequency_"+sTitle, freqip);
			}
			else
			{
				finalFFTfr = new ImagePlus("Dominant_FFT_period_"+sTitle, freqip);
			}
			finalFFTfr.setCalibration(outCal);
			finalFFTfr.show();			
			IJ.run(finalFFTfr,"Fire","");
			finalFFTphase = new ImagePlus("Phase_"+sTitle, phaseip);
			finalFFTphase.setCalibration(outCal);
			finalFFTphase.show();
			IJ.run(finalFFTphase,"Spectrum","");
		}
		
	}
	
	/** subtracts a specified stack projection**/
	private void subtractProj() {
		ZProjector projector = new ZProjector(imp);
		projector.setMethod(nProjMethod);
		projector.doProjection();
		ImagePlus stackProjection = projector.getProjection();
		ImageCalculator ic = new ImageCalculator();
	    imp = ic.run("Subtract create stack", imp, stackProjection);
	}

	/** function to calculate TICS, returns "rotated stack", 
	 * where time instead of z becomes x.
	 * stores result in the ipICSRotated (!!) 
	 * in case 	bICSfreqphase is true, also calculates 
	 * dominant frequency and corresponding phase per pixel,
	 * stored in isFreqPhase **/

	public void ICS_pixel()
	{
		//ImageStack finalIS ;
		ImageProcessor ipRotatedSlice;
		//long startTime;
		//long detectionTime=0;

		//startTime = System.nanoTime();
		//int i;
		//float [][] fFreqPhase;
		
		if(bICSfreqphase)
		{ 
			freqip = new FloatProcessor(origW,origH);
			phaseip = new FloatProcessor(origW,origH);
		}
		
		IJ.showStatus("Calculating temporal correlation...");
		ipICSRotated = new ImageStack(nMaxUserDelayICS+1,origH);
		for (nindW=0;nindW<origW;nindW++)
		{
			ipRotatedSlice=getResliceRotatedSlice(imp, nindW);
			ipICSRotated.addSlice(xCorrSpace(ipRotatedSlice, nMaxUserDelayICS,nCalcMethod,nNormMethod));

			IJ.showProgress(nindW, origW-1);
		}
		
		IJ.showProgress(3, 3);
		IJ.showStatus("Calculating temporal correlation..done");
		//new ImagePlus("frequency",freqip).show();
		//new ImagePlus("phase",phaseip).show();
		//return finalIS;

		
		//detectionTime = System.nanoTime() - startTime;
		//IJ.log("Processing time: " + String.format("%.2f",((double)Math.abs(detectionTime))*0.000000001) + " s");
	}
	/** 
	 * Dialog displaying options of Temporal ICS correlation
	 * **/
	public boolean xDialogTICS() 
	{
		GenericDialog xTICSDial = new GenericDialog("Temporal image correlation spectroscopy");



		xTICSDial.addMessage("Preprocessing:");
		xTICSDial.addCheckbox("Subtract stack projection?", Prefs.get("Correlescence.ICSsubtract", false));
		xTICSDial.addChoice("Projection type:", projMethods, Prefs.get("Correlescence.ICSprojection", "Average Intensity"));
		xTICSDial.addMessage("~~~~~~~  TICS  ~~~~~~~");
		xTICSDial.addNumericField("Maximum period/delay (zero=half of duration)", Prefs.get("Correlescence.nMaxUserDelayICS", 0), 0, 4, " frames");
		//xTICSDial.addChoice("CC Normalization:", itemsN, Prefs.get("Correlescence.ICSnormmethod", "Normalize by overlap area"));
		xTICSDial.addCheckbox("Resize CC image?", Prefs.get("Correlescence.ICSrescale", false) );
		xTICSDial.addNumericField("Width resized:", Prefs.get("Correlescence.ICSrescaleW", 200), 0,5," ");
		xTICSDial.addNumericField("Height resized:", Prefs.get("Correlescence.ICSrescaleH", 200), 0,5," ");
		xTICSDial.addCheckbox("Estimate period/frequency?", Prefs.get("Correlescence.ICSfindfirstmax", false));
		xTICSDial.addNumericField("Tolerance of first CC maximum >", Prefs.get("Correlescence.dMaxTolICS", 0.25), 1,7," ");
		xTICSDial.addCheckbox("Subframe precision?", Prefs.get("Correlescence.ICSsubframepeak", true));
		xTICSDial.addCheckbox("Show TICS stack?", Prefs.get("Correlescence.ICSshowTICSstack", false));
		xTICSDial.addMessage("~~~~~~~  FFT  ~~~~~~~");
		xTICSDial.addCheckbox("Calculate dominant frequency/phase map?", Prefs.get("Correlescence.bICSfreqphase", false));
		xTICSDial.addMessage("~~~~~~~  Output  ~~~~~~~");
		xTICSDial.addChoice("Output image value:", sOutput, Prefs.get("Correlescence.ICSoutput", "Frequency"));
		xTICSDial.addChoice("Units of output:", sOutputUnits, Prefs.get("Correlescence.ICSoutputunits", "Frames"));
		xTICSDial.setResizable(false);
		xTICSDial.showDialog();
		if (xTICSDial.wasCanceled())
	        return false;

		bSubtractProj = xTICSDial.getNextBoolean(); 
		Prefs.set("Correlescence.ICSsubtract", bSubtractProj);
		
		nProjMethod= xTICSDial.getNextChoiceIndex();
		Prefs.set("Correlescence.ICSprojection", projMethods[nProjMethod]);
		
		nMaxUserDelayICS = (int)xTICSDial.getNextNumber();
		Prefs.set("Correlescence.nMaxUserDelayICS", nMaxUserDelayICS);
		
		bRescaleCC = xTICSDial.getNextBoolean(); 
		Prefs.set("Correlescence.ICSrescale", bRescaleCC);
		
		nScaledW=(int)xTICSDial.getNextNumber();
		Prefs.set("Correlescence.ICSrescaleW", nScaledW);
		
		nScaledH=(int)xTICSDial.getNextNumber();
		Prefs.set("Correlescence.ICSrescaleH", nScaledH);
		
		bFindFirstMax = xTICSDial.getNextBoolean(); 
		Prefs.set("Correlescence.ICSfindfirstmax", bFindFirstMax);
		
		dMaxTol=xTICSDial.getNextNumber();
		Prefs.set("Correlescence.dMaxTolICS", dMaxTol);
		
		bSubFramePeak = xTICSDial.getNextBoolean(); 
		Prefs.set("Correlescence.ICSsubframepeak", bSubFramePeak);
		
		bShowTICSStack = xTICSDial.getNextBoolean(); 
		Prefs.set("Correlescence.ICSshowTICSstack", bShowTICSStack);
		
		bICSfreqphase = xTICSDial.getNextBoolean(); 
		Prefs.set("Correlescence.bICSfreqphase", bICSfreqphase);
		
		nOutput = xTICSDial.getNextChoiceIndex();
		Prefs.set("Correlescence.ICSoutput", sOutput[nOutput]);
		
		nOutputUnits = xTICSDial.getNextChoiceIndex();
		Prefs.set("Correlescence.ICSoutputunits", sOutputUnits[nOutputUnits]);

		return true;
	}
	
	/** Function gets one slice from a resliced (rotated) stack at nWidthPos, 
	 * so time becomes x and height becomes y. 
	 * It does not verify if nWidthPos is inside the image
	 * **/
	public ImageProcessor getResliceRotatedSlice(ImagePlus imp_in, int nWidthPos)
	{
		ImageStack origStack;
		ImageProcessor ip,rotatedIP=null;
		int k;
		float[] line = null;
		int stackSize, inHeight; 
		
		origStack=imp_in.getStack();
		inHeight = imp_in.getHeight();
		
		stackSize = origStack.getSize();
		rotatedIP = new FloatProcessor(stackSize,inHeight);
		ip=origStack.getProcessor(1);

			for(k=0;k<stackSize;k++)
			{
				ip=origStack.getProcessor(k+1);
				line=FFTtools.getImageCol(ip,nWidthPos, inHeight);
				//putImageRow(ip2,k,inHeight,line);
				FFTtools.putImageCol(rotatedIP,k,inHeight,line);
			
			}

	
		
		return rotatedIP;
	}
	
	/** Function reslices (rotates) stack, so time becomes x and height becomes y) 
	 * **/
	public ImageStack resliceRotateStack(ImagePlus imp_in)
	{
		
		ImageStack rotatedIm;
		ImageStack origStack;
		ImageProcessor ip,ip2=null;
		int i,k;
		float[] line = null;
		int stackSize, inWidth, inHeight; 
		
		origStack=imp_in.getStack();
		inWidth = imp_in.getWidth();
		inHeight = imp_in.getHeight();
		
		stackSize = origStack.getSize();
		rotatedIm = new ImageStack(stackSize,inHeight);
		
		ip=origStack.getProcessor(1);
		for (i=0;i<inWidth;i++)
		{
			ip2=ip.createProcessor(stackSize,inHeight);
			for(k=0;k<stackSize;k++)
			{
				ip=origStack.getProcessor(k+1);
				line=FFTtools.getImageCol(ip,i, inHeight);
				//putImageRow(ip2,k,inHeight,line);
				FFTtools.putImageCol(ip2,k,inHeight,line);
			
			}
			rotatedIm.addSlice(ip2);
		}
		return rotatedIm;
	 
	}
	/** function to find first maximum with given tolerance 
	 * in the (rotated) stack of TICS **/
	public FloatProcessor findMax(ImageStack imp_in, double tolerance)
	{
		int stackSize, inWidth, inHeight; 
		int j,k;
		
		FloatProcessor ipFin;		
		
		ImageProcessor ip;
		
		double[] line = null;
		
		int [] maxpos;
		
		int maxInd;
		
		double finVal = 0.0;
		
		int arrlength;
		
		
		inWidth = imp_in.getWidth();
		inHeight = imp_in.getHeight();
		stackSize = imp_in.getSize();
		ipFin = new FloatProcessor(stackSize,inHeight);
		
		IJ.showStatus("Estimating frequency/period..");
		
		for (k=0;k<stackSize;k++)
		{
			ip = imp_in.getProcessor(k+1);
			IJ.showProgress(k, stackSize-1);
			
			for (j=0;j<inHeight;j++)
			{
				line = getImageRowD(ip,j,inWidth);
				
				maxpos = MaximumFinder.findMaxima(line, tolerance, true);
				
				Arrays.sort(maxpos);
								
				arrlength = maxpos.length;
				
				if(arrlength == 0)
				{
					maxInd = 0;					
				}
				else
				{
					maxInd = maxpos[0];
				}
				
				
				if(maxInd > 1 && maxInd <line.length-1 && bSubFramePeak)
				{
					//centroid estimate
					//according to https://iopscience.iop.org/article/10.3847/2515-5172/aae265
					// A Robust Method to Measure Centroids of Spectral Lines
					// Richard Teague  and Daniel Foreman-Mackey
					// DOI 10.3847/2515-5172/aae265
					finVal= maxInd - 0.5*(line[maxInd+1]-line[maxInd-1])/(line[maxInd+1]+line[maxInd-1]-2*line[maxInd]);					
				}
				else
				{
					finVal = maxInd;
				}
				//convert to time units (1 if frames)
				finVal *= dTimeUnits;
				
				//in case frequency needs to be reported
				if(nOutput == 0)
				{
					finVal = 1./finVal;
				}
				
				//special case
				if(maxInd == 0)
				{
					//frequency
					if(nOutput == 0)
					{
						finVal = 0.0;
					}
					//period
					else
					{
						finVal = Double.POSITIVE_INFINITY;
					}
				}
				
			
				//store result
				ipFin.putPixelValue(k, j, finVal);
				
			}
			
		}
		IJ.showProgress(3, 3);
		IJ.showStatus("Estimating frequency/period..done");
		return ipFin;
	}
	
	/** given input image it treats x coordinates as time,
	 * calculates FFT, finds a frequency with maximum amplitude
	 * and returns its value[0] and phase[1] in float array equal to image height**/
	public float [][] findFreqPhase(ImageProcessor ip)
	{
		int imW, imH;
		int x,f;
		imW=ip.getWidth();
		imH=ip.getHeight();
		float [][] fFreqPhase = new float [2][imH];
		float [] rowt;
		float [][] fft;
		float ampMax, amp;
		int indMax;
		int Nfin=Math.floorDiv(imW, 2)+1;
		
		for(x=0;x<imH;x++)
		{
			//get time sequence
			rowt=FFTtools.getImageRow(ip, x, imW);
			//remove average
			imCC1D.subtractAverage(rowt);
			//get FFT
			fft = GeneralFFT.transformR(rowt);
			//find index of frequency with max amlitude
			ampMax=Float.MIN_VALUE;
			indMax=0;
			for(f=0;f<Nfin+1;f++)
			{
				amp=(float)(Math.pow(fft[0][f],2)+Math.pow(fft[1][f],2));
				if(ampMax<amp)
				{
					ampMax=amp;
					indMax=f;
				}
			}
			if(indMax==0)
			{
				fFreqPhase[0][x]=0;
			}
			else
			{
				fFreqPhase[0][x]=(float)(indMax)/(float)imW;
			}
			fFreqPhase[1][x]=(float)Math.atan2(fft[1][indMax], fft[0][indMax]);
		}
		
		return fFreqPhase;
		
		
	}
	
	/** Function calculates temporal auto-correlation (in one direction, which is X) 
	 *  
	 *  
	 * **/
	
	public FloatProcessor xCorrSpace(ImageProcessor ipp1, int nMaxDelay,int nCalcMethod, int nNormMethod)
	{
		int imWidth, imHeight;
		FloatProcessor resultip;
		imWidth=ipp1.getWidth();
		imHeight=ipp1.getHeight();
		float [] fCurrCC;
		float [] row1;
		float [] rowx;
		imCC1D cc1D = new imCC1D();
		
		int  resWidth; //final width of the image, depending on the maximum time
		int i,t;
		
		int nBegin = (int)Math.round(imWidth*0.5-1);
		
		resWidth=nMaxDelay+1;
		
		resultip = new FloatProcessor(resWidth, imHeight);
		
		for(t=0;t<imHeight;t++)
		{
			//IJ.showProgress(dt, resHeight-1);
			row1=cc1D.getImageRow(ipp1, t, imWidth);
			imCC1D.subtractAverage(row1);

			//direct cross-correlation
			if(nCalcMethod==0)
			{
				fCurrCC = cc1D.calcCCRowsDirect(row1,row1);
			}
			//FFT cross-correlation
			else
			{
				rowx=cc1D.prepareForFFT(row1);
				fCurrCC = cc1D.calcACRowsFFT(rowx,bICSfreqphase);
				//fCurrCC = calcCCRowsFFT(rowx, rowx);
				fCurrCC = cc1D.trimFFTCC(fCurrCC, imWidth);
			}
			//normalize
			cc1D.normalizeCC(row1,row1,fCurrCC,nNormMethod);	
			
			//put in final
			for(i=0;i<resWidth;i++)
			//for(i=0;i<imWidth;i++)
			{
				resultip.setf(i, t,fCurrCC[nBegin+i]);
			}
			if(bICSfreqphase)
			{ 
				if(nOutput==0)
				{
					freqip.setf(nindW,t,(float) (cc1D.fFreqPhase[0]/dTimeUnits));
				}
				else
				{
					freqip.setf(nindW,t,(float) (dTimeUnits/cc1D.fFreqPhase[0]));
				}
				phaseip.setf(nindW,t,cc1D.fFreqPhase[1]);
			}
			
		}
		
		
		return resultip;
	}
	
	/** function returns float array corresponding to the image 
	 * @param ip image processor
	 * @param nRow row number
	 * @param nWidth image width
	 * **/
	public double [] getImageRowD(ImageProcessor ip, int nRow, int nWidth)
	{
		double [] data = new double [nWidth];

		for(int i=0;i<nWidth;i++)
			data[i]= (double)ip.getf(i,nRow);
		return data;
	}
	


}
