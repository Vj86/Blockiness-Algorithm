import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.util.Arrays;

public class UtilitiesBlock {

	/**
     * Takes a BufferedImage and convert this into a grayscale image
     * 
     *
     */
	static BufferedImage ConvertToGrayScale(BufferedImage inputImg)
	{
	    ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);
	    ColorConvertOp op = new ColorConvertOp(cs, null);
	    BufferedImage grayScaleImg = new BufferedImage(inputImg.getWidth(), inputImg.getHeight(), BufferedImage.TYPE_BYTE_INDEXED);
	    grayScaleImg = op.filter(inputImg, null);
	    return grayScaleImg;
	    
	}
	
	/**
     * Takes a BufferedImage and convert this into a two-dimensional double array
     * 
     * 
     */
	static double[][] ImageToMatrix(BufferedImage inputImg)
	{
		int height = inputImg.getHeight();
	    int width = inputImg.getWidth();
	    //int r=0, g=0, b=0;
	    double[][] matrixImage = new double[height][width];
	    for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {   
	    	   /*
	    	   r = inputImg.getRGB(i, j)>> 16 & 0xFF;
	    	   g = inputImg.getRGB(i, j)>> 8 & 0xFF;
	    	   b = inputImg.getRGB(i, j) & 0xFF;
	    	   */
	    	   //System.out.println( (inputImg.getRGB(j, i) & 0xFF) );
	    	   matrixImage[i][j] = ((inputImg.getRGB(j, i)>> 16) & 0xFF);
	       }
	    }
	    return matrixImage;
	}
	/**
     * Takes a two-dimensional double array and convert this into a BufferedImage
     *
     *
     */
	static BufferedImage MatrixToBufferedImage(double[][] inputMtrixImg, int height, int width )
	{
		BufferedImage output = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
	    for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	  output.setRGB(j, i, (int)inputMtrixImg[i][j]);
	       }
	    }
	    return output;
	}
	
	/**
     * Takes a two-dimensional double array and returns it's absolute value
     *
     *
     */
	static double[][] MatrixAbs(double[][] inputMatrix)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
	    double[][] output = new double[height][width];
		output = inputMatrix.clone();
		for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	  output[i][j] = Math.abs(output[i][j]);
	       }
	    }
		
		return output;
	}
	
	/**
     * Takes a two-dimensional double array and returns the maximu value
     * over row and colums
     *
     * 
     */
	static double GetMax(double[][] inputMatrix)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
		double maximum = 0.0;
		for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	  if(inputMatrix[i][j]>maximum)
	    		  maximum=inputMatrix[i][j];
	       }
	    }
		return maximum;
	}
	
	/**
     * Takes two two-dimensional double array and returns the maximum value
     * over row and colums
     *
     *
     */
	static double[][] GetMatrixMax(double[][] inputMatrix_A, double[][] inputMatrix_B)
	{
		int height = inputMatrix_A.length;
		int width = inputMatrix_A[0].length;
		double[][] output = new double[height][width];
		for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	   if(inputMatrix_A[i][j] > inputMatrix_B[i][j])
	    		   output[i][j] = inputMatrix_A[i][j];
	    	   else if(inputMatrix_A[i][j] <= inputMatrix_B[i][j])
	    		   output[i][j] = inputMatrix_B[i][j];
	       }
	    }
		return output;
	}
	
	/**
     * Takes a two-dimensional double array and divide each element by the 'div' value
     *
     
     */
	static double[][] MatrixDiv(double[][] inputMatrix, double divisor)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
	    double[][] output = new double[height][width];
		output = inputMatrix.clone();
		for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	  output[i][j] = output[i][j]/divisor;
	       }
	    }
		return output;
	}
	/**
     * Takes two two-dimensional double array and compute element-by-element product
     *
     
     */
	static double[][] MatrixRightProduct(double[][] inputMatrix_A, double[][] inputMatrix_B)
	{
		int height = inputMatrix_A.length;;
		int width = inputMatrix_A[0].length;
		
		//check Dimensions
		if(inputMatrix_A.length != inputMatrix_B.length)
		{
			System.out.println("Height of the two matrix is different. Get the smallest one.");
			if(inputMatrix_A.length < inputMatrix_B.length)
				height = inputMatrix_A.length;
			else
				height = inputMatrix_B.length;
		}
	    if(inputMatrix_A[0].length != inputMatrix_B[0].length)
	    {
	    	System.out.println("Width of the two matrix is different. Get the smallest one");
	    	if(inputMatrix_A[0].length < inputMatrix_B[0].length)
	    		width = inputMatrix_A[0].length;
	    	else
	    		width = inputMatrix_B[0].length;
	    }

		//compute the product element-by-element
	    double[][] output = new double[height][width];
		for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	  output[i][j] = ( (inputMatrix_A[i][j]*inputMatrix_A[i][j]) + (inputMatrix_B[i][j]*inputMatrix_B[i][j]) );
	       }
	    }
		
		return output;	
	}
	/**
     * Takes a two-dimensional double array and apply the square root element-by-element
     *
     * 
     */
	static double[][] MatrixSquareRoot(double[][] inputMatrix)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
		double[][] output = new double[height][width];
		for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	   output[i][j] = Math.sqrt(inputMatrix[i][j]);
	       }
	    }
		return output;
	}
	
	/**
     * Takes a two-dimensional double array and sort its colums in descending mode
     *
     *
     */	
	static double[][] SortMatrix(double[][] inputMatrix)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
		
		double[][] output = new double[height][width];
		double[] tmp = new double[width];
		
        for(int i=0; i<width; i++)
        {  
     	    tmp = inputMatrix[i].clone();
        	Arrays.sort(tmp);
            for(int j=width-1, k=0; j>=0; j--, k++)
        	   output[i][k] = tmp[j];  
        }
    	return output;
	 }

	/**
     * Takes a two-dimensional double array and sort its colums in descending mode
     *
     * 
     */	
	static double[] SortVector(double[][] inputMatrix)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
		if(height==1 && width==1)
		{
			double[] output = new double[1];
			output[0]=inputMatrix[0][0];
			return output;
		}
		else
		{
			double[] output = new double[height*width];
	        for(int i=0; i<height; i++)
	        {  
	            for(int j=0; j<width; j++)
	            {
	        	   output[j+(width*i)] = inputMatrix[i][j];
	            }
	        }
	        Arrays.sort(output);
	    	return output;
		}
	 }	
	
	
	/**
     * Takes a two-dimensional double array and return an array with the mean over its colums 
     *
     * 
     */	
	static double[] MeanMatrix(double[][] inputMatrix)
	{
		int width = inputMatrix[0].length;
		
		double[] output = new double[width];
        for(int i=0; i<width; i++)
        {  
        	output[i]= ArrayAvg(inputMatrix[i]);
        }
    	return output;
	 }

	/**
     * Takes a one-dimensional double array and return its mean 
     *
     * 
     */	
	static double ArrayAvg(double[] inputArray)
	{
		double avg = 0.0;
		if(inputArray.length == 1)
		{
			return (inputArray[0]);
		}
		else
		{
			for(int i=0; i<inputArray.length; i++)
	        {
				avg+=inputArray[i];
	        }
			if(inputArray.length>0)
				return (avg/inputArray.length);
			else
				return -1.0;	
		}

	}
	
	/**
     * Takes a two-dimensional double array and return a sub-matrix starting from index and 
     * of subMatrixSize size
     *
     *
     */
	static double[][] GetSubMatrix(double[][] inputMatrix, int x, int y, int subMatrixSize)
	{
		double[][] output = new double[subMatrixSize][subMatrixSize];
		
	    for(int i=x, h=0; h<subMatrixSize; i++, h++)
	    {
	    	for(int j=y, k=0; k<subMatrixSize; j++, k++)
	    	{
	    		output[h][k]=inputMatrix[i][j];
	    	}
	    }
		return output;
	}
	
	/**
     * Reset to zero the input matrix border
     *
     * 
     */
	static double[][] ResetBorderMatrix(double[][] inputMatrix)
	{
		int height = inputMatrix.length;
		int width = inputMatrix[0].length;
	    double[][] output = new double[height][width];
		output = inputMatrix.clone();		
		for(int i=0; i<height; i++)
        {
			//Sx(1:end,1)=0;
			//Row 1 to n; Column 1 = 0 
			output[i][0] = 0.0;
			//Sx(1:end,end)=0;
			//Row 1 to n; Column m = 0
			output[i][width-1] = 0.0;
        }
		for(int i=0; i<width; i++)
        {
			//Sx(1,1:end)=0;
			//Row 1; Column 1 to n = 0
			output[0][i] = 0.0;
			//Sx(end,1:end)=0;
			//Row n; Column 1 to m = 0
			output[height-1][i] = 0.0;
        }
		return output;
	}
	
	static double[][] ApplyZeroPadding(double[][] inputMatrix, int rowPad, int colPad)
	{
		int height = inputMatrix.length;
		int width =  inputMatrix[0].length;
		int startN, startM;
		//single padding row, add at end of matrix
		if(rowPad==1)
		{
			startN = 0;
		}
		//even padding row, half at start and half at end of matrix
		else if(rowPad%2==0)
		{
			startN = (rowPad/2);
		}
		//odd padding row, (int) rowPad/2 at start of matrix, (int) rowPad/2 + 1 at end
		else
		{
			startN =  (int) rowPad/2;
		}
		
		//-----column count-----
		
		//single padding row, add at end of matrix
		if(colPad==1)
		{
			startM = 0;
		}
		//even padding row, half at start and half at end of matrix
		else if(colPad%2==0)
		{
			startM = (colPad/2);
		}
		//odd padding row, (int) rowPad/2 at start of matrix, (int) rowPad/2 + 1 at end
		else
		{
			startM =  (int) colPad/2;
		}
		
		double[][] output = new double[height+rowPad][width+colPad];
	    int i = 0, j = 0, h = 0, k = 0;
		for( i=startN, h=0; h<height; i++, h++)
	    {
	    	for( j=startM, k=0; k<width; j++, k++)
	    	{
	    		output[i][j]=inputMatrix[h][k];
	    	}
	    }
		return output;
	}
	
	/**
     * Concatenate two arrays
     *
     *
     */
	static double[] ConcatenateArray(double[] a, double[] b)
	{
	   int aLen = a.length;
	   int bLen = b.length;
	   double[] c= new double[aLen+bLen];
	   System.arraycopy(a, 0, c, 0, aLen);
	   System.arraycopy(b, 0, c, aLen, bLen);
	   return c;
	}
	
	/**
     * Compute blockiness x-by-y over input sub-matrix for x directiom
     *
     *
     */
	double blockinessS(double[][] inputMatrix, int kernelSize)
	{
		double S = 0.0;  
		double[][] O3 = new double[kernelSize][kernelSize];
		//Kernel initializzation
		for(int i=0; i<kernelSize; i++)
	    {
	       for(int j=0; j<kernelSize; j++)
	       {  
	    	   O3[i][j] = 0.0;
	       }
	    }
		for(int i=1; i<kernelSize-1; i++)
	    {
	       for(int j=1; j<kernelSize-1; j++)
	       {  
	    	   O3[i][j] = 1.0;
	       }
	    }
		for(int i=2; i<kernelSize-2; i++)
	    {
	       for(int j=2; j<kernelSize-2; j++)
	       {  
	    	   O3[i][j] = 0.0;
	       }
	    }
		//product input submatrix * kernel
		for(int i=0; i<kernelSize; i++)
	    {
	       for(int j=0; j<kernelSize; j++)
	       {  
	    	   O3[i][j] = O3[i][j]*inputMatrix[i][j];
	       }
	    }
		O3 = MatrixAbs(O3);
		double sum = 0.0;
		//return array contains sum over column
		for(int i=0; i<kernelSize; i++)
	    {
			for(int j=0; j<kernelSize; j++)
		    {
				sum+=O3[j][i];
		    }
	    }
		S = sum/20.0;
		return S;
	}
	
	
	/**
     * Compute blockiness x-by-y over input sub-matrix for x directiom
     *
     * 
     */
	static double blockinessSx(double[][] inputMatrix, int kernelSize)
	{
		double Sx = 0.0;  
		double[][] O1 = new double[kernelSize][kernelSize];
		//Kernel initializzation
		for(int i=0; i<kernelSize; i++)
	    {
	       for(int j=0; j<kernelSize; j++)
	       {  
	    	   O1[i][j] = 0.0;
	       }
	    }
		for(int i=0; i<kernelSize; i++)
	    {
			O1[i][0] = 1.0;
			O1[i][kernelSize-1] = 1.0;
	    }
		//product input submatrix * kernel
		for(int i=0; i<kernelSize; i++)
	    {
	       for(int j=0; j<kernelSize; j++)
	       {  
	    	   O1[i][j] = O1[i][j]*inputMatrix[i][j];
	       }
	    }
		O1 = MatrixAbs(O1);
		double sum = 0.0;
		//return array contains sum over column
		for(int i=0; i<kernelSize; i++)
	    {
			for(int j=0; j<kernelSize; j++)
		    {
				sum+=O1[j][i];
		    }
	    }
		Sx = sum/16.0;
		//System.out.println("Sx " + Sx );
		return Sx;
	}
	
	
	/**
     * Compute blockiness x-by-y over input sub-matrix for y directiom
     *
     *
     */
	static double blockinessSy(double[][] inputMatrix, int kernelSize)
	{
		double Sy = 0.0;  
		double[][] O2 = new double[kernelSize][kernelSize];
		//Kernel initializzation
		for(int i=0; i<kernelSize; i++)
	    {
	       for(int j=0; j<kernelSize; j++)
	       {  
	    	   O2[i][j] = 0.0;
	       }
	    }
		for(int i=0; i<kernelSize; i++)
	    {
			O2[0][i] = 1.0;
			O2[kernelSize-1][i] = 1.0;
	    }
		//product input submatrix * kernel
		for(int i=0; i<kernelSize; i++)
	    {
	       for(int j=0; j<kernelSize; j++)
	       {  
	    	   O2[i][j] = O2[i][j]*inputMatrix[i][j];
	       }
	    }
		O2 = MatrixAbs(O2);
		double sum = 0.0;
		//return array contains sum over column
		for(int i=0; i<kernelSize; i++)
	    {
			for(int j=0; j<kernelSize; j++)
		    {
				sum+=O2[j][i];
		    }
	    }
		Sy = sum/16.0;
		return Sy;
	}
	
	/**
     * Calculate Sdiff
     *
     * 
     */
	static double[][] CalculateSdiff(double[][] inputMatrix_A, double[][] inputMatrix_B, double k)
	{
		
		int height = inputMatrix_A.length;
		int width = inputMatrix_A[0].length;
	    double[][] output = new double[height][width];
	    //abs(M.^k-S2.^k)./((abs(M.^k+S2.^k)+0.0000001)/2);
	    for(int i=0; i<height; i++)
	    {
	       for(int j=0; j<width; j++)
	       {  
	    	   output[i][j] = Math.abs( Math.pow(inputMatrix_A[i][j],k)-Math.pow(inputMatrix_B[i][j],k) ) / ( ( ( Math.abs( Math.pow(inputMatrix_A[i][j],k) + Math.pow(inputMatrix_B[i][j],k)) )+0.0000001) / 2.0);  
	       }
	    }
		return output;
	}
}
