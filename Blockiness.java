
import java.awt.image.BufferedImage;
import java.awt.image.ConvolveOp;
import java.awt.image.Kernel;
import javax.imageio.ImageIO;
import java.io.File;
import java.io.IOException;

import org.bytedeco.javacpp.BytePointer;
import org.bytedeco.javacpp.opencv_core.IplImage;
import org.bytedeco.javacpp.opencv_core.Mat;
import static org.bytedeco.javacpp.opencv_imgproc.CV_RGB2GRAY;
import static org.bytedeco.javacpp.opencv_imgproc.cvCvtColor;


public class Blockiness
{
	BufferedImage inputImage, grayScaleImg;
	IplImage inIPL;
	IplImage grey;	
	ConvolveOp conv; 
	String path;
	int w;
	int h;
	int blockSize = 0;  
	double max = 0.0;
	double blockiness = 0.0;
	double k_pow = 2.3;
	double[][] matrixGSImage;
	
	//---Calculation variables---
	//Matrix
	double[][] Dx;
	double[][] Dy; 
	double[][] D;
	double[][] Sx;
	double[][] Sy;
	double[][] S;
	double[][] M;
	double[][] Sdiff;
	double[] SdiffSort;
	
	//Convolution kernel

	float[] Mx_row = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
	float[] My_row = { -1, -2, -1, 0, 0, 0, 1, 2, 1 };
	
	/**
     
     * 
     * @param File Path to the image (jpeg or png).
     * @param blkSize representing the size of single block.
     */
	
	public Blockiness(String ImagePath, int blkSize)
	{
		//check blocksize value
		if(blkSize < 1 )
		{
			System.out.println("Error: Zero BlockSize not allowed!");
			inputImage = null;
			return;			
		}
		else
		{
			this.blockSize = blkSize;
		}
		//check path to image
		this.path = new String(ImagePath);
		
		if(!this.path.isEmpty())
		{
			File Fp = new File(path);
			if (!Fp.exists())
			{
			   System.out.println("Image File not found! " + path);
			   inputImage = null;
			   return;
			}
			else
			{
				try
				{
					inputImage = ImageIO.read(Fp);
				}
				catch (IOException e)
				{
					System.out.println("Unable to read input file " + path + "\nError: " + e.getMessage());
					inputImage = null;
					return;
				}
			}
		}
		else
		{
			System.out.println("Invalid Path!");
			inputImage = null;
			return;
		}
	}
	
	/**
    
     * 
     *  BufferedImage representing the image (jpeg or png).
     *  blkSize representing the size of single block.
     */
	public Blockiness(BufferedImage inputImage, int blkSize)
	{
		//Check input image 
		if(inputImage == null)
		{
			System.out.println("Null pointer to input image (BufferedImage). Error!");
			this.inputImage = null;
			return;
		}
		else
		{
			this.inputImage = inputImage;
		}
		//check block size
		if(blkSize < 1 )
		{
			System.out.println("Error: Zero BlockSize not allowed!");
			this.inputImage = null;
			return;			
		}
		else
		{
			this.blockSize = blkSize;
		}
		
		this.path = "BufferedImage_used";
	}
	
	/**
     * 
     * 
     *  IplImage representing the input image.
     *  blkSize representing the size of single block.
     */
	public Blockiness(IplImage inputImage, int blkSize)
	{
		//check block size
		if(blkSize < 1 )
		{
			System.out.println("Error: Zero BlockSize not allowed!");
			this.inputImage = null;
			return;			
		}
		else
		{
			this.blockSize = blkSize;
		}
		
		//check and convert IplImage to BufferedImage
		if(inputImage == null)
		{
			System.out.println("Null pointer to input image (IplImage). Error!");
			this.inputImage = null;
			return;
		}
		else
		{		
			this.inIPL = inputImage;
			this.path = "BufferedImage_used";
			int type = BufferedImage.TYPE_BYTE_GRAY;
			
	        if(this.inIPL.nChannels() > 1)
	        {
	        	grey = IplImage.create(this.inIPL.width(), this.inIPL.height(), this.inIPL.depth(), 1);
	            cvCvtColor(this.inIPL, grey, CV_RGB2GRAY);
	        }
	        
			//IplImage to BufferedImage Conversion
			this.inputImage = new BufferedImage(this.inIPL.width(), this.inIPL.height(), type); 
			Mat mat = new Mat(this.grey);
			BytePointer dt = mat.arrayData();
	        byte[] data = new byte[mat.rows()*mat.cols()*(int)(mat.elemSize())];
			dt.get(data);
	        this.inputImage.getRaster().setDataElements(0, 0, mat.cols(), mat.rows(), data);
		}
	}
	/**
     * Compute the Blockiness Factor
     */
	public double CalculateBlockinessImage() throws Exception
	{
		if(this.inputImage==null || blockSize == 0)
		{
			
			return -1.0;
		}	
		
		//Convert image into grayscale only for not IplImage input 
		if(this.inIPL == null && this.inputImage != null)
		{
			grayScaleImg = UtilitiesBlock.ConvertToGrayScale(inputImage);
		}
		else
		{
			grayScaleImg = inputImage;
		}

        h = grayScaleImg.getHeight();
        w = grayScaleImg.getWidth();
        
        //Convert Image into double 2D matrix
        matrixGSImage = UtilitiesBlock.ImageToMatrix(grayScaleImg);    

        //-------------------------------        
        conv = new ConvolveOp(new Kernel(3,3,Mx_row));
        BufferedImage tmpImage = conv.createCompatibleDestImage(grayScaleImg,null);
        conv.filter(grayScaleImg, tmpImage);
        //-------------------------------        
        
        Dx = UtilitiesBlock.ImageToMatrix(tmpImage);
        
        //Absolute Value of Convolved Matrix.
        Dx = UtilitiesBlock.MatrixAbs(Dx);
        
        //find maximum value over Dx
        max = UtilitiesBlock.GetMax(Dx);
        
        //divide Dx by the maximum value
        Dx = UtilitiesBlock.MatrixDiv(Dx, max);
        
        //Verify matrix size vs blockSize proportion and zero padding if necessary
        if( Dx.length%blockSize != 0 || Dx[0].length%blockSize != 0)
        {
        	int hPadding = blockSize - Dx.length%blockSize;
        	int wPadding = blockSize - Dx[0].length%blockSize;
        	//Apply zero padding to image
        	Dx = UtilitiesBlock.ApplyZeroPadding(Dx, hPadding, wPadding);
        }

        //Apply 2D Convolution using Canny Y-Filter
        //--------------------------
        conv = new ConvolveOp(new Kernel(3,3,My_row));
        tmpImage = conv.createCompatibleDestImage(grayScaleImg,null);
        conv.filter(grayScaleImg, tmpImage);
        //--------------------------
        
        Dy = UtilitiesBlock.ImageToMatrix(tmpImage);    
        
        //Absolute Value of Convolved Matrix
        Dy = UtilitiesBlock.MatrixAbs(Dy);
        
        //find maximum value over Dy
        max = UtilitiesBlock.GetMax(Dy);
        
        //divide Dy by the maximum value
        Dy = UtilitiesBlock.MatrixDiv(Dy, max);

        //Verify matrix size vs blockSize proportion and zero padding if necessary
        if( Dy.length%blockSize != 0 || Dy[0].length%blockSize != 0)
        {
        	int hPadding = blockSize - Dy.length%blockSize;
        	int wPadding = blockSize - Dy[0].length%blockSize;
        	//Apply zero padding to image
        	Dy = UtilitiesBlock.ApplyZeroPadding(Dy, hPadding, wPadding);
        }     
        
        //compute the product of Dx by Dy
        D = UtilitiesBlock.MatrixRightProduct(Dx, Dy);
		
        //compute the square root element by element
        D = UtilitiesBlock.MatrixSquareRoot(D);
        
        //compute the maximum value
        max = UtilitiesBlock.GetMax(D);
        
        //divide D by the maximum value
        D = UtilitiesBlock.MatrixDiv(D, max);
        //Verify matrix size vs blockSize proportion and zero padding  if necessary
        if( D.length%blockSize != 0 || D[0].length%blockSize != 0)
        {
        	int hPadding = blockSize - D.length%blockSize;
        	int wPadding = blockSize - D[0].length%blockSize;
        	//Apply zero padding to image
        	D = UtilitiesBlock.ApplyZeroPadding(D, hPadding, wPadding);
        }
        
        //new matrix image size
        h = D.length;
        w = D[0].length;
            
        //compute m-by-n blockbproc func x-direction
        Sx = new double[h/blockSize][w/blockSize];
        for(int i=0; i<h; i+=blockSize)
        {	
        	for(int j=0; j<w; j+=blockSize)
        	{
        		if(i==0 && j!=0)
        			Sx[0][j/blockSize]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dx,i,j,blockSize),blockSize);
        		else if(i!=0 && j==0)	
        			Sx[i/blockSize][0]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dx,i,j,blockSize),blockSize);
        		else if(i==0 && j==0)
        			Sx[0][0]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dx,i,j,blockSize),blockSize);
        		else
        			Sx[i/blockSize][j/blockSize]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dx,i,j,blockSize),blockSize);
        		
        	}
        }

        //compute m-by-n blockbproc func y-direction
        Sy = new double[h/blockSize][w/blockSize];
        for(int i=0; i<h; i+=blockSize)
        {	
        	for(int j=0; j<w; j+=blockSize)
        	{
        		if(i==0 && j!=0)
        			Sy[0][j/blockSize]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dy,i,j,blockSize),blockSize);
        		else if(i!=0 && j==0)	
        			Sy[i/blockSize][0]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dy,i,j,blockSize),blockSize);
        		else if(i==0 && j==0)
        			Sy[0][0]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dy,i,j,blockSize),blockSize);
        		else
        			Sy[i/blockSize][j/blockSize]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(Dy,i,j,blockSize),blockSize);
        	}
        }
        
        //compute m-by-n blockbproc func full matrix
        S = new double[h/blockSize][w/blockSize];
        for(int i=0; i<h; i+=blockSize)
        {	
        	for(int j=0; j<w; j+=blockSize)
        	{
        		if(i==0 && j!=0)
        			S[0][j/blockSize]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(D,i,j,blockSize),blockSize);
        		else if(i!=0 && j==0)	
        			S[i/blockSize][0]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(D,i,j,blockSize),blockSize);
        		else if(i==0 && j==0)
        			S[0][0]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(D,i,j,blockSize),blockSize);
        		else
        			S[i/blockSize][j/blockSize]= UtilitiesBlock.blockinessSx(UtilitiesBlock.GetSubMatrix(D,i,j,blockSize),blockSize);
        	}
        }
        //Matrix border Reset
        Sx = UtilitiesBlock.ResetBorderMatrix(Sx);
        Sy = UtilitiesBlock.ResetBorderMatrix(Sy);        
        S = UtilitiesBlock.ResetBorderMatrix(S);
        
        //Obtain Max from Sx and Sy
        M = UtilitiesBlock.GetMatrixMax(Sx,Sy);
        //compute Sdiff matrix
        Sdiff = UtilitiesBlock.CalculateSdiff(M, S, k_pow);

        //compute the maximum value
        max = UtilitiesBlock.GetMax(Sdiff);
        
        //Normalize Sdiff
        Sdiff = UtilitiesBlock.MatrixDiv(Sdiff, max);
        
        //get submatrix form Sdiff
        if( (Sdiff.length > 2) && (Sdiff[0].length > 2) )
        {
	        Sdiff = UtilitiesBlock.GetSubMatrix(Sdiff, 1, 1, Sdiff.length-2);
	        
	        //sort Sdiff over column
	        SdiffSort = new double[Sdiff.length*Sdiff[0].length];
	        SdiffSort = UtilitiesBlock.SortVector(Sdiff);
	        
	        //compute avg of sorted Sdiff
	        blockiness = UtilitiesBlock.ArrayAvg(SdiffSort);
	        //return blockiness value
	        return blockiness;
        }
        else
        {
        	System.out.println("Unable to calculate blockiness factor. Image size too small!");
            return -1.0;
        }

	}

}

