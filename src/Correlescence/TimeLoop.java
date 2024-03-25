package Correlescence;

import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.RankFilters;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class TimeLoop implements PlugIn {
	
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
	
	
	/** input stack size **/
	int nStackSize;
	
	/** image width in px**/
	int origW;
	
	/** image height in px**/
	int origH;
	
	/** index of reslicing **/
	int nindW;
	
	/** rotated stack containing TICS **/
	ImageStack ipICSRotated;
	
	/** rotated stack containing image **/
	ImageStack isInputRotated;
	
	/** maximum finding tolerance **/
	double dMaxTol;
	
	/** method to locate max **/
	int nMaxMethod;
	
	/** number of periods to average **/
	int nPeriodsN;
	
	/** image processor containing dominating frequency **/
	ImageProcessor freqip = null;
	
	/** image processor containing dominating phase **/
	ImageProcessor phaseip = null; 
	
	/** image processor containing max of CC **/
	ImageProcessor ipMaxCC = null;
	
	/** show a stack with max position **/
	boolean bShowMaxStackTL = false;
	
	/** rotated stack containing TICS **/
	ImageStack isMaxLocations;
	
	@Override
	public void run(String arg) {
		
		ImagePlus finalImpMax = null;
		ImagePlus finalImpMaxNew = null;
		String sTitle;
		
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
		
		if(!xDialogTimeLoop())
			return;
		
		sTitle=imp.getTitle();
		Calibration outCal = (Calibration) imp.getCalibration().clone();

		
		IJ.log(" --- Correlescence plugin version " + ConstantsCorrelescence.sVersion+ " --- ");
		IJ.log("Running Time Loop.");
		IJ.log("Image title: \"" + sTitle + "\"");
		

		//run TICS
		ICS_pixel();
		
		/*
		//filter CC stack
		ImagePlus ipRescale = new ImagePlus("TICS_"+sTitle, resliceRotateStack(new ImagePlus("rotateback",ipICSRotated)));
		
		ImageStack isRescale = new ImageStack(origW, origH);
		ImageProcessor singleImage;
		for (int i=0;i<ipRescale.getImageStackSize();i++)
		{
			ipRescale.setPosition(i+1);
			singleImage = ipRescale.getProcessor();
			singleImage.setInterpolationMethod(ImageProcessor.BILINEAR);
			singleImage.smooth();
			isRescale.addSlice(singleImage);
		}
		int nHalfStack=(int)Math.floor((double)nStackSize*0.5);
		ipICSRotated = new ImageStack(nHalfStack+1,origH);
		for (nindW=0;nindW<origW;nindW++)
		{
			ipICSRotated.addSlice(getResliceRotatedSlice(ipRescale, nindW));	

		}
		*/
		//finalImp = new ImagePlus("TICS_"+sTitle, Temporal_ICS.resliceRotateStack(new ImagePlus("rotateback",ipICSRotated)));
		//finalImp.setCalibration(outCal);
		//finalImp.show();
		
		FloatProcessor ipPeriod = null;
		ipMaxCC = new FloatProcessor(origW,origH);
		switch(nMaxMethod)
		{
			case 0:
				ipPeriod = findFirstMax(ipICSRotated,dMaxTol);
				break;
			case 1:
				ipPeriod = findMaxMax(ipICSRotated,dMaxTol);
				break;
			case 2:
				ipPeriod = findLastMax(ipICSRotated,dMaxTol);
				break;
				
		}
		//RankFilters filt = new RankFilters();
		//filt.rank(ipPeriod, 4.0, RankFilters.MEDIAN);
		//ipPeriod.smooth();
		new ImagePlus("MaxCC_"+sTitle, ipMaxCC).show();
		finalImpMax = new ImagePlus("Period_"+sTitle, ipPeriod);
		finalImpMax.setCalibration(outCal);	
		finalImpMax.show();
		IJ.run(finalImpMax,"Fire","");
		freqip = new FloatProcessor(origW,origH);
		phaseip = new FloatProcessor(origW,origH);
		
		FloatProcessor ipPeriodNew = estimateTimeLoop(isInputRotated, ipPeriod);
		finalImpMaxNew =  new ImagePlus("Period_NEW_"+sTitle, ipPeriodNew);
		IJ.run(finalImpMaxNew,"Fire","");
		finalImpMaxNew.show();
		
		new ImagePlus("looped_" + sTitle, Temporal_ICS.resliceRotateStack(new ImagePlus("looped",isInputRotated))).show();
		if(bShowMaxStackTL)
		{
			new ImagePlus("CC_max_pos_" + sTitle, Temporal_ICS.resliceRotateStack(new ImagePlus("looped",isMaxLocations))).show();
		}
		//new ImagePlus("rotated orig",isInputRotated).show();
	}
	
	public FloatProcessor estimateTimeLoop(ImageStack imp_in, FloatProcessor ipPeriod)
	{
		int stackSize, inWidth, inHeight; 
		int j,k,i;
		//int circ_index;
		
		FloatProcessor ipFinPeriod;	
		
		float fPeriodEstimate;
		int nPeriodEstimate;
		double dFinalPeriod = 0.0;
		int nFinalPeriod;
		double[] line = null;
		float [] firstIntPeriod = null;
		float [] currentIntPeriod = null;
		float [] firstIntPeriodPad = null;
		float [] currentIntPeriodPad = null;
		double [] averageIntPeriod = null;
		int nPerCycle;
		int [] dLags; 
		double [] dLagsPeriods;
		int nPeriodsMax;
		//int [] dLags = new int[nPeriodsN];
		//double [] dLagsPeriods = new double[nPeriodsN-1];
		ImageProcessor ip;
		int dCCZeroPos;
		int nStartIndex;
		
		float [] fCrossCorr;
		
		inWidth = imp_in.getWidth();
		inHeight = imp_in.getHeight();
		stackSize = imp_in.getSize();
		IJ.showStatus("TimeLooping..");
		
		ipFinPeriod = new FloatProcessor(stackSize,inHeight);
		
		for (k=0;k<stackSize;k++)
		{
			ip = imp_in.getProcessor(k+1);
			IJ.showProgress(k, stackSize-1);

			for (j=0;j<inHeight;j++)
			{
				fPeriodEstimate =ipPeriod.getPixelValue(k,j); 
				if(Float.isNaN(fPeriodEstimate)||fPeriodEstimate<1.0)
				{
					fPeriodEstimate = inWidth;
				}
				nPeriodsMax = Math.floorDiv(inWidth, Math.round(fPeriodEstimate));

				//if(Float.isNaN(fPeriodEstimate)||(fPeriodEstimate*(nPeriodsN+0.5)>inWidth))
				if(nPeriodsMax<3)
				//if(Float.isNaN(fPeriodEstimate)||(fPeriodEstimate*(nPeriodsMax+1.0)>inWidth))
				{
					
					for(i=0;i<inWidth;i++)
					{
						ip.putPixelValue(i, j, 0.0);
					}
					ipFinPeriod.putPixelValue(k, j, Double.NaN);
				}
				else
				{
					nPeriodsMax--;
					dLags = new int[nPeriodsMax];
					dLagsPeriods = new double[nPeriodsMax-1];
					line = Temporal_ICS.getImageRowD(ip,j,inWidth);
					nPeriodEstimate = Math.round(fPeriodEstimate);
					dCCZeroPos = (int)Math.round(nPeriodEstimate*0.5-1);
					firstIntPeriod = new float[nPeriodEstimate];
					currentIntPeriod = new float[nPeriodEstimate];
					//averageIntPeriod = new double[dPeriodEstimate];
					for(i=0;i<nPeriodEstimate;i++)
					{
						firstIntPeriod[i] = (float) line[i];
						//averageIntPeriod[i] = line[i];
					}
					
					imCC1D.subtractAverage(firstIntPeriod);
					firstIntPeriodPad = imCC1D.prepareForFFT(firstIntPeriod);
					
					//estimate lags for each period
					//nPeriodsMax
					for (nPerCycle = 1;nPerCycle<nPeriodsMax;nPerCycle++)
					//for (nPerCycle = 1;nPerCycle<nPeriodsN;nPerCycle++)
					{
						for(i=0;i<nPeriodEstimate;i++)
						{
							currentIntPeriod[i] = (float)line[i+nPeriodEstimate*nPerCycle];
						}
						imCC1D.subtractAverage(currentIntPeriod);
						currentIntPeriodPad = imCC1D.prepareForFFT(currentIntPeriod);
						
						//calculate cross-correlation with the first period
						fCrossCorr = imCC1D.calcCCRowsFFT(firstIntPeriodPad, currentIntPeriodPad);
						fCrossCorr = imCC1D.trimFFTCC(fCrossCorr, nPeriodEstimate);
						imCC1D.normalizeCC(firstIntPeriod,currentIntPeriod,fCrossCorr,0);
						//find max
						dLags[nPerCycle] = 0;
						for (i=1;i<fCrossCorr.length;i++)
						{
							if(fCrossCorr[i]>fCrossCorr[dLags[nPerCycle]])
							{
								dLags[nPerCycle]=i;
							}
						}
						dLags[nPerCycle] -=dCCZeroPos; 
						dLagsPeriods[nPerCycle-1]= (nPeriodEstimate*(nPerCycle)+dLags[nPerCycle])/((double)nPerCycle);
					}
					//get new period estimate
					dFinalPeriod = 0.0;
					for(i=0;i<dLagsPeriods.length;i++)
					{
						dFinalPeriod += dLagsPeriods[i];
					}
					
					dFinalPeriod/=(nPeriodsMax-1);
					//dFinalPeriod/=(nPeriodsN-1);
					ipFinPeriod.putPixelValue(k, j, dFinalPeriod);
					nFinalPeriod = (int) Math.ceil(dFinalPeriod);
					averageIntPeriod = new double[nFinalPeriod];
					
					for (nPerCycle = 0;nPerCycle<nPeriodsMax;nPerCycle++)
					//for (nPerCycle = 0;nPerCycle<nPeriodsN;nPerCycle++)
					{
						nStartIndex = nPeriodEstimate*nPerCycle+dLags[nPerCycle];
						for(i=0;i<nFinalPeriod;i++)
						{
							averageIntPeriod[i] +=line[nStartIndex+i];
						}
					}
					double [] xVals = new double[nFinalPeriod];
				 	for(i=0;i<nFinalPeriod;i++)
					{
						//averageIntPeriod[i] /= nPeriodsN;
						averageIntPeriod[i] /= nPeriodsMax;
						xVals[i]=i;
						
					}

					
					//put values
					for(i=0;i<inWidth;i++)
					{
						//if(i<nFinalPeriod)
						//{
						//	ip.putPixelValue(i, j, averageIntPeriod[i]);
						//}
						//else
						//{
							//ip.putPixelValue(i, j, 0.0);
						//}
						//circ_index = i-(int)Math.round(Math.floor(((float)i)/nFinalPeriod)*nFinalPeriod);
						//ip.putPixelValue(i, j, averageIntPeriod[circ_index]);
						
						ip.putPixelValue(i, j,LinearInterpolation.evalLinearInterp(xVals, averageIntPeriod, ((double)i)%dFinalPeriod));
					}
				}
			
				
			}
			
		}
		IJ.showProgress(3, 3);
		IJ.showStatus("TimeLooping..done");
		return ipFinPeriod;
	}
	
	public void ICS_pixel()
	{

		ImageProcessor ipRotatedSlice;
	
		freqip = new FloatProcessor(origW,origH);
		phaseip = new FloatProcessor(origW,origH);
		
		IJ.showStatus("Calculating Time Loop parameters...");
		int nHalfStack=(int)Math.floor((double)nStackSize*0.5);
		ipICSRotated = new ImageStack(nHalfStack+1,origH);
		isMaxLocations = new ImageStack(nHalfStack+1,origH);
		isInputRotated = new ImageStack(nStackSize,origH);
		for (nindW=0;nindW<origW;nindW++)
		{
			ipRotatedSlice = Temporal_ICS.getResliceRotatedSlice(imp, nindW);
			isInputRotated.addSlice(ipRotatedSlice);
			ipICSRotated.addSlice(xCorrSpace(ipRotatedSlice, nHalfStack,nCalcMethod,nNormMethod));
			isMaxLocations.addSlice(new FloatProcessor(nHalfStack+1,origH));

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
	
	
	/** function to find first maximum with given tolerance 
	 * in the (rotated) stack of TICS **/
	public FloatProcessor findFirstMax(ImageStack imp_in, double tolerance)
	{
		int stackSize, inWidth, inHeight; 
		int j,k;
		
		FloatProcessor ipFin;		
		
		ImageProcessor ip;
		
		double[] line = null;
		
		int [] maxpos;
		
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
				line = Temporal_ICS.getImageRowD(ip,j,inWidth);

				maxpos = MaximumFinder.findMaxima(line, tolerance, true);
				
				Arrays.sort(maxpos);
								
				arrlength = maxpos.length;
				
				
				if(arrlength == 0)
				{
					finVal = 0.0;
					//finVal = Double.NaN;					
				}
				else
				{

					finVal = maxpos[0];
				}
					
				//store result in the 
				ipFin.putPixelValue(k, j, finVal);
				
			}
			
		}
		IJ.showProgress(3, 3);
		IJ.showStatus("Estimating period..done");
		return ipFin;
	}
	
	/** function to find first maximum with given tolerance 
	 * in the (rotated) stack of TICS **/
	public FloatProcessor findLastMax(ImageStack imp_in, double tolerance)
	{
		int stackSize, inWidth, inHeight; 
		int j,k;
		
		FloatProcessor ipFin;		
		
		ImageProcessor ip;
		
		double[] line = null;
		
		int [] maxpos;
		
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
			
			for (j=00;j<inHeight;j++)
			{
				line = Temporal_ICS.getImageRowD(ip,j,inWidth);

				maxpos = MaximumFinder.findMaxima(line, tolerance, true);
				
				Arrays.sort(maxpos);
								
				arrlength = maxpos.length;
				
				
				if(arrlength == 0)
				{
					finVal = 0.0;					
				}
				else
				{

					finVal = maxpos[arrlength-1];
				}
					
				//store result in the 
				ipFin.putPixelValue(k, j, finVal);
				
			}
			
		}
		IJ.showProgress(3, 3);
		IJ.showStatus("Estimating period..done");
		return ipFin;
	}
	
	/** function to find first maximum with given tolerance 
	 * in the (rotated) stack of TICS **/
	public FloatProcessor findMaxMax(ImageStack imp_in, double tolerance)
	{
		int stackSize, inWidth, inHeight; 
		int j,k,i;
		
		FloatProcessor ipFin;		
		
		ImageProcessor ip;
		ImageProcessor ipMacLoc = null;
		
		double[] line = null;
		
		int [] maxpos;
		
		double [][] maxposInt;
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
			if(bShowMaxStackTL)
			{
				ipMacLoc = isMaxLocations.getProcessor(k+1);
			}
			IJ.showProgress(k, stackSize-1);
			
			for (j=0;j<inHeight;j++)
			{
				line = Temporal_ICS.getImageRowD(ip,j,inWidth);

				maxpos = MaximumFinder.findMaxima(line, tolerance, true);
				
				Arrays.sort(maxpos);
								
				arrlength = maxpos.length;
				
				
				if(arrlength == 0)
				{
					//finVal = Double.NaN;					
					finVal = 0.0;
				}
				else
				{
					maxposInt = new double [arrlength][2];
					for(i=0;i<arrlength;i++)
					{
						maxposInt[i][0]=maxpos[i];
						maxposInt[i][1]=line[maxpos[i]];
						if(bShowMaxStackTL)
						{
							ipMacLoc.putPixelValue(maxpos[i], j, 1.0);
							//ipMacLoc.putPixelValue(maxpos[i], j, line[maxpos[i]]);
						}
					}
					java.util.Arrays.sort(maxposInt, new java.util.Comparator<double[]>() {
					    public int compare(double[] a, double[] b) {
					        return (-1)*Double.compare(a[1], b[1]);
					    }
					});
					finVal = maxposInt[0][0];
					ipMaxCC.putPixelValue(k, j, maxposInt[0][1]);
					
					/*
					finVal = maxposInt[0][0];
					int indMin =0;
					boolean bFlag = true;
					while(bFlag)
					{
						finVal = maxposInt[indMin][0];
						//long rrrr = Math.floorDiv(inWidth, Math.round(finVal));
						if(Math.floorDiv(nStackSize, Math.round(finVal))<3)
						{
							indMin++;
							if(indMin==arrlength)
								bFlag=false;
						}
						else
						{
							bFlag=false;
						}
					}*/
					//if(Math.floorDiv(inWidth, Math.round(finVal))<3)
						//finVal =maxposInt[0][0];
					
				}
					
				//store result in the 
				ipFin.putPixelValue(k, j, finVal);
				
				
			}
			
		}
		IJ.showProgress(3, 3);
		IJ.showStatus("Estimating period..done");
		return ipFin;
	}
	/** function to find first maximum with given tolerance 
	 * in the (rotated) stack of TICS **/
	public FloatProcessor findMaxFFT(ImageStack imp_in, double tolerance)
	{
		int stackSize, inWidth, inHeight; 
		int j,k,i;
		
		FloatProcessor ipFin;		
		
		ImageProcessor ip;
		
		double[] line = null;
		
		int [] maxpos;
		
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
				line = Temporal_ICS.getImageRowD(ip,j,inWidth);

				maxpos = MaximumFinder.findMaxima(line, tolerance, true);
				
				Arrays.sort(maxpos);
							
				arrlength = maxpos.length;
				
				
				if(arrlength == 0)
				{
					finVal = Double.NaN;					
				}
				else
				{
					float fftPer = freqip.getPixelValue(k, j);
					finVal = maxpos[0];
					int nDiff = 1000*inWidth;
					for (i=0;i<arrlength;i++)
					{
						if(Math.abs(maxpos[i]-fftPer)<nDiff)
						{
							finVal=maxpos[i];
							nDiff=(int) Math.abs(maxpos[i]-fftPer);
						}
					}
					finVal=fftPer;
					
				}
					
				//store result in the 
				ipFin.putPixelValue(k, j, finVal);
				
			}
			
		}
		IJ.showProgress(3, 3);
		IJ.showStatus("Estimating period..done");
		return ipFin;
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
			//row1=smoothfl(row1);
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
				fCurrCC = cc1D.calcACRowsFFT(rowx,true);
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
			freqip.setf(nindW,t,(float) (1.0/cc1D.fFreqPhase[0]));
			phaseip.setf(nindW,t,cc1D.fFreqPhase[1]);
			
		}
		
		
		return resultip;
	}
	
	/** 
	 * Dialog displaying options of TimeLoop
	 * **/
	public boolean xDialogTimeLoop() 
	{
		GenericDialog xTimeLoopDial = new GenericDialog("Time Loop");
		String[] maxFindMethods = 
			{"First maximum", "Highest CC maximum", "Last Maximum"}; 

		xTimeLoopDial.addNumericField("Tolerance of first CC maximum >", Prefs.get("Correlescence.dMaxTolTL", 0.25), 1,7," ");
		xTimeLoopDial.addChoice("Max search:", maxFindMethods, Prefs.get("Correlescence.maxFindTL", "First maximum"));
		xTimeLoopDial.addCheckbox("Show Max stack?", Prefs.get("Correlescence.bShowMaxStackTL", false));
		xTimeLoopDial.addNumericField("Number of periods to average", Prefs.get("Correlescence.nPeriodsN",2), 2);

		xTimeLoopDial.setResizable(false);
		xTimeLoopDial.showDialog();
		if (xTimeLoopDial.wasCanceled())
	        return false;

		dMaxTol=xTimeLoopDial.getNextNumber();
		Prefs.set("Correlescence.dMaxTolTL", dMaxTol);
		
		nMaxMethod= xTimeLoopDial.getNextChoiceIndex();
		Prefs.set("Correlescence.maxFindTL", maxFindMethods[nMaxMethod]);
		bShowMaxStackTL = xTimeLoopDial.getNextBoolean(); 
		Prefs.set("Correlescence.bShowMaxStackTL", bShowMaxStackTL);
		nPeriodsN=(int)Math.round(xTimeLoopDial.getNextNumber());
		Prefs.set("Correlescence.nPeriodsN", nPeriodsN);

		return true;
	}
	
	/** Function reslices (rotates) stack, so time becomes x and height becomes y) 
	 * **/
	static public ImageStack resliceRotateStack(ImagePlus imp_in)
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
	
	/** Function gets one slice from a resliced (rotated) stack at nWidthPos, 
	 * so time becomes x and height becomes y. 
	 * It does not verify if nWidthPos is inside the image
	 * **/
	static public ImageProcessor getResliceRotatedSlice(ImagePlus imp_in, int nWidthPos)
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
	double[] smooth(double[] a) {
		int n = a.length;
		double [] out = new double[n];
		int i,j;
		//no point in averaging 2 points
		if(n<3)
			return a;
		
		
		//preserve end points
		out[0]=	a[0];
		out[n-1] = a[n-1];
		
		for (i=1; i<(n-1); i++)
		{
			out[i]=0.0f;
			for(j=(i-1);j<(i+2);j++)
			{
				out[i] += a[j];
			}
			out[i] /= 3.0f;
		}
		return out;
	}
	float[] smoothfl(float[] a) {
		int n = a.length;
		float [] out = new float[n];
		int i,j;
		//no point in averaging 2 points
		if(n<3)
			return a;
		
		
		//preserve end points
		out[0]=	a[0];
		out[n-1] = a[n-1];
		
		for (i=1; i<(n-1); i++)
		{
			out[i]=0.0f;
			for(j=(i-1);j<(i+2);j++)
			{
				out[i] += a[j];
			}
			out[i] /= 3.0f;
		}
		return out;
	}
}
