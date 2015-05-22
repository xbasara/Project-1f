#include <cmath>
#include <iostream>
#include <algorithm>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <math.h>

using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;

class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(void) {;};
    Matrix          CameraTransform(void) {;};
    Matrix          DeviceTransform(void) {;};
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

double linearInterpolateZ(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3) {
	double r;
	if(y1 != y2) {
		r = (y3 - y1)/(y2 - y1);
	} else {
		r = (x3 - x1)/(x2 - x1);
	}
	return z1 + r*(z2 - z1);
}

double * linearInterpolateColor(double x1, double y1, double color1[], double x2, double y2, double color2[], double x3, double y3) {
	double* returnColor = (double*) malloc(sizeof(double)*3);
	double r;
	if(y1 != y2) {
		r = (y3 - y1)/(y2 - y1);
	} else {
		r = (x3 - x1)/(x2 - x1);
	}
	for(int i = 0; i < 3; i++) {
		returnColor[i] = color1[i] + r*(color2[i] - color1[i]);
	}
	return returnColor;
}

double * linearInterpolateVector(double x1, double y1, double vector1[], double x2, double y2, double vector2[], double x3, double y3) {
	double* returnVector = (double*) malloc(sizeof(double)*3);
	double r;
	if(y1 != y2) {
		r = (y3 - y1)/(y2 - y1);
	} else {
		r = (x3 - x1)/(x2 - x1);
	}
	for(int i = 0; i < 3; i++) {
		returnVector[i] = vector1[i] + r*(vector2[i] - vector1[i]);
	}
	return returnVector;
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = +0.6;
         lightDir[1] = 0;
         lightDir[2] = +0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    }
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
	  double		 Z[3];
      double         colors[3][3];
	  double		 normals[3][3];

	  int leftBaseIndex;
	  int rightBaseIndex;
	  int offPointIndex;
	  
	  int topPointIndex;
	  int midPointIndex;
	  int bottomPointIndex;
	    
	  /*int id;
	  
	  void setId(int i) {
		id = i;
	  }
	  /int getId() {
		return id;
	  }*/
	  bool arrangePoints() {
		if(Y[0] == Y[1]) {
			offPointIndex = 2;
			if(X[0] < X[1]) {
				rightBaseIndex = 1;
				leftBaseIndex = 0;
			} else {
				rightBaseIndex = 0;
				leftBaseIndex = 1;
			}
			return true;
		} else if(Y[0] == Y[2]) {
			offPointIndex = 1;
			if(X[0] < X[2]) {
				rightBaseIndex = 2;
				leftBaseIndex = 0;
			} else {
				rightBaseIndex = 0;
				leftBaseIndex = 2;
			}
			return true;
		} else if(Y[1] == Y[2]) {
			offPointIndex = 0;
			if(X[1] < X[2]) {
				rightBaseIndex = 2;
				leftBaseIndex = 1;
			} else {
				rightBaseIndex = 1;
				leftBaseIndex = 2;
			}
			return true;
		} else {
			if(Y[0] < Y[1]) {
				if(Y[0] < Y[2]) {
					bottomPointIndex = 0;
					if(Y[1] < Y[2]) {
						midPointIndex = 1;
						topPointIndex = 2;
					} else {
						midPointIndex = 2;
						topPointIndex = 1;
					}
				} else {
					bottomPointIndex = 2;
					midPointIndex = 0;
					topPointIndex = 1;
				}
			} else {
				if(Y[1] < Y[2]) {
					bottomPointIndex = 1;
					if(Y[0] < Y[2]) {
						midPointIndex = 0;
						topPointIndex = 2;
					} else {
						midPointIndex = 2;
						topPointIndex = 0;
					}
				} else {
					bottomPointIndex = 2;
					midPointIndex = 1;
					topPointIndex = 0;
				}
			}			
			return false;
		}
	  }
	  
	  double getLeftIntercept(double val) {
		if(X[leftBaseIndex] == X[offPointIndex]) {
			return X[leftBaseIndex];
		}
		/* if(X[leftBaseIndex]/X[offPointIndex] < 1.01 && X[leftBaseIndex]/X[offPointIndex] > .99) {
			return X[rightBaseIndex];
		} */
		
		double leftM = (Y[leftBaseIndex] - Y[offPointIndex]) / (X[leftBaseIndex] - X[offPointIndex]);
		double leftB = Y[leftBaseIndex] - leftM * X[leftBaseIndex];

		return (val - leftB) / leftM;
	  }
	  
	  double getRightIntercept(double val) {
		if(X[rightBaseIndex] == X[offPointIndex]) {
			return X[rightBaseIndex];
		}
		/* if(X[rightBaseIndex]/X[offPointIndex] < 1.01 && X[rightBaseIndex]/X[offPointIndex] > .99) {
			return X[rightBaseIndex];
		} */ 
		double rightM = (Y[rightBaseIndex] - Y[offPointIndex]) / (X[rightBaseIndex] - X[offPointIndex]);
		double rightB = Y[rightBaseIndex] - rightM * X[rightBaseIndex];
		
		//cout << "val: " << val << ", m: " << rightM << ", b; " << rightB << ", ret: " << (val - rightB) / rightM << endl;
		return (val - rightB) / rightM;
	  }
	  
	  double getSplitX() {
		/*if(id == 748956) {
			cout << X[topPointIndex] << "," << Y[topPointIndex] << " " << X[bottomPointIndex] << "," << Y[bottomPointIndex] << endl;
		}*/
		if(X[topPointIndex] == X[bottomPointIndex]) {
			return X[topPointIndex];
		}
		/* if(X[topPointIndex]/X[bottomPointIndex] < 1.01 && X[topPointIndex]/X[bottomPointIndex] > .99) {
			return X[topPointIndex];
		} */
		double splitM = (Y[topPointIndex] - Y[bottomPointIndex]) / (X[topPointIndex] - X[bottomPointIndex]);
		double splitB = Y[topPointIndex] - splitM * X[topPointIndex];
		
		return (Y[midPointIndex] - splitB) / splitM;
	  }
	  
	  double getLeftBaseX() {
		return X[leftBaseIndex];
	  }
	  double getLeftBaseY() {
		return Y[leftBaseIndex];
	  }
	  double getLeftBaseZ() {
		return Z[leftBaseIndex];
	  }
	  double * getLeftBaseColor() {
		return colors[leftBaseIndex];
	  }
	  double * getLeftBaseVector() {
		return normals[leftBaseIndex];
	  }
	  double getRightBaseX() {
		return X[rightBaseIndex];
	  }
	  double getRightBaseY() {
		return Y[rightBaseIndex];
	  }
	  double getRightBaseZ() {
		return Z[rightBaseIndex];
	  }
	  double * getRightBaseColor() {
		return colors[rightBaseIndex];
	  }
	  double * getRightBaseVector() {
		return normals[rightBaseIndex];
	  }
	  double getOffX() {
		return X[offPointIndex];
	  }
	  double getOffY() {
		return Y[offPointIndex];
	  }
	  double getOffZ() {
		return Z[offPointIndex];
	  }
	  double * getOffColor() {
		return colors[offPointIndex];
	  }
	  double * getOffVector() {
		return normals[offPointIndex];
	  }
	  double getTopPointX() {
		return X[topPointIndex];
	  }
	  double getTopPointY() {
		return Y[topPointIndex];
	  }
	  double getTopPointZ() {
		return Z[topPointIndex];
	  }
	  double * getTopPointColor() {
		return colors[topPointIndex];
	  }
	  double * getTopPointVector() {
		return normals[topPointIndex];
	  }
	  double getMidPointX() {
		return X[midPointIndex];
	  }
	  double getMidPointY() {
		return Y[midPointIndex];
	  }
	  double getMidPointZ() {
		return Z[midPointIndex];
	  }
	  double * getMidPointColor() {
		return colors[midPointIndex];
	  }
	  double * getMidPointVector() {
		return normals[midPointIndex];
	  }
	  double getBottomPointX() {
		return X[bottomPointIndex];
	  }
	  double getBottomPointY() {
		return Y[bottomPointIndex];
	  }
	  double getBottomPointZ() {
		return Z[bottomPointIndex];
	  }
	  double * getBottomPointColor() {
		return colors[bottomPointIndex];
	  }
	  double * getBottomPointVector() {
		return normals[bottomPointIndex];
	  }
  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
	  double *depthBuffer;
      int width, height;
  
  
	  void setPixel(int x, int y, double color[3], double z, double shading) {
		if(x<0 || x > (width-1) || y<0 || y > (height-1) ) {
			return;
		}
		/*if(x == 0 && y == 0)
			cout << "(" << color[0] << "," << color[1] << "," << color[2] << ")" << endl;
		*/
		/*if(ceil441((255*color[0]))  > 255) {
		 cout << ceil441((255*color[0])) << endl;
		}*/
		if(z > depthBuffer[width*y + x]) {
			depthBuffer[width*y + x] = z;
			buffer[3*(width*y + x)] = (unsigned char)ceil441(255*min(1.0,color[0]*shading));
			buffer[3*(width*y + x) + 1] = (unsigned char)ceil441(255*min(1.0,color[1]*shading));
			buffer[3*(width*y + x) + 2] = (unsigned char)ceil441(255*min(1.0,color[2]*shading));
		}
	  }
	  
	  double calculateShading(LightingParameters & lp, double *view, double *vector) {
		//double vectorNorm = sqrt(pow(vector[0], 2) + pow(vector[1],2) + pow(vector[2],2));
		//cout << vectorNorm << endl;
		//cout << vector[0]*view[0] + vector[1]*view[1] + vector[2]*view[2] << ", " << vectorNorm << endl;
		
		double LN = vector[0]*lp.lightDir[0] + vector[1]*lp.lightDir[1] + vector[2]*lp.lightDir[2];
		double diffuse = LN; //*vectorNorm;
		diffuse = std::abs(diffuse);
		/*if(diffuse != 0 && diffuse != 1) {
			cout << diffuse << endl;
		}*/
		
		double r[3];
		for(int i = 0; i<3; i++) {
			r[i] = (2*LN*vector[i]) - lp.lightDir[i];
		}
		
		double specular = max(0.0, pow(view[0]*r[0]+view[1]*r[1]+view[2]*r[2], lp.alpha));
		
		return lp.Ka + (lp.Kd * diffuse) + (lp.Ks * specular);
	  }
	  
	  void paintTriangle(Triangle t, LightingParameters & lp, double view[3]) {
		
		if(t.getLeftBaseY() < t.getOffY()) {
			
			for(int j = (int)ceil441(t.getLeftBaseY()); j <= (int)floor441(t.getOffY()); j++) {
				double leftZ = linearInterpolateZ(t.getLeftBaseX(), t.getLeftBaseY(), t.getLeftBaseZ(), t.getOffX(), t.getOffY(), t.getOffZ(), t.getLeftIntercept(j), (double)j);
				double rightZ = linearInterpolateZ(t.getRightBaseX(), t.getRightBaseY(), t.getRightBaseZ(), t.getOffX(), t.getOffY(), t.getOffZ(), t.getRightIntercept(j), (double)j);
				double* leftColor = linearInterpolateColor(t.getLeftBaseX(), t.getLeftBaseY(), t.getLeftBaseColor(), t.getOffX(), t.getOffY(), t.getOffColor(), t.getLeftIntercept(j), (double)j);
				double* rightColor = linearInterpolateColor(t.getRightBaseX(), t.getRightBaseY(), t.getRightBaseColor(), t.getOffX(), t.getOffY(), t.getOffColor(), t.getRightIntercept(j), (double)j);
				double* leftVector = linearInterpolateVector(t.getLeftBaseX(), t.getLeftBaseY(), t.getLeftBaseVector(), t.getOffX(), t.getOffY(), t.getOffVector(), t.getLeftIntercept(j), (double)j);
				double* rightVector = linearInterpolateVector(t.getRightBaseX(), t.getRightBaseY(), t.getRightBaseVector(), t.getOffX(), t.getOffY(), t.getOffVector(), t.getRightIntercept(j), (double)j);
				for(int k = (int)ceil441( t.getLeftIntercept(j) ); k <= (int)floor441( t.getRightIntercept(j) ); k++) {
					double z = linearInterpolateZ(t.getLeftIntercept(j), (double)j, leftZ, t.getRightIntercept(j), (double)j, rightZ, (double)k, (double)j); 
					double* color = linearInterpolateColor(t.getLeftIntercept(j), (double)j, leftColor, t.getRightIntercept(j), (double)j, rightColor, (double)k, (double)j);
					double* vector = linearInterpolateVector(t.getLeftIntercept(j), (double)j, leftVector, t.getRightIntercept(j), (double)j, rightVector, (double)k, (double)j);
					double shading = calculateShading(lp, view, vector);
					//cout << shading << endl;
					//cerr << color[0] << "," << color[1] << "," << color[2] << endl;
					setPixel(k,j,color,z,shading);
				}
			}
		} else {
			
			for(int j = (int)ceil441(t.getOffY()); j <= (int)floor441(t.getLeftBaseY()); j++) {
				double leftZ = linearInterpolateZ(t.getOffX(), t.getOffY(), t.getOffZ(), t.getLeftBaseX(), t.getLeftBaseY(), t.getLeftBaseZ(), t.getLeftIntercept(j), (double)j);
				double rightZ = linearInterpolateZ(t.getOffX(), t.getOffY(), t.getOffZ(), t.getRightBaseX(), t.getRightBaseY(), t.getRightBaseZ(), t.getRightIntercept(j), (double)j);
				double* leftColor = linearInterpolateColor(t.getOffX(), t.getOffY(), t.getOffColor(), t.getLeftBaseX(), t.getLeftBaseY(), t.getLeftBaseColor(), t.getLeftIntercept(j), (double)j);
				double* rightColor = linearInterpolateColor(t.getOffX(), t.getOffY(), t.getOffColor(), t.getRightBaseX(), t.getRightBaseY(), t.getRightBaseColor(), t.getRightIntercept(j), (double)j);
				double* leftVector = linearInterpolateVector(t.getOffX(), t.getOffY(), t.getOffVector(), t.getLeftBaseX(), t.getLeftBaseY(), t.getLeftBaseVector(), t.getLeftIntercept(j), (double)j);
				double* rightVector = linearInterpolateVector(t.getOffX(), t.getOffY(), t.getOffVector(), t.getRightBaseX(), t.getRightBaseY(), t.getRightBaseVector(), t.getRightIntercept(j), (double)j);
				for(int k = (int)ceil441( t.getLeftIntercept(j) ); k <= (int)floor441( t.getRightIntercept(j) ); k++) {
					double z = linearInterpolateZ(t.getLeftIntercept(j), (double)j, leftZ, t.getRightIntercept(j), (double)j, rightZ, (double)k, (double)j); 
					double* color = linearInterpolateColor(t.getLeftIntercept(j), (double)j, leftColor, t.getRightIntercept(j), (double)j, rightColor, (double)k, (double)j);
					double* vector = linearInterpolateVector(t.getLeftIntercept(j), (double)j, leftVector, t.getRightIntercept(j), (double)j, rightVector, (double)k, (double)j);
					double shading = calculateShading(lp, view, vector);
					//cout << shading << endl;
					setPixel(k,j,color,z,shading);
				}
			}
		}
	  }
		
  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1f_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
/*
vtkDataSetWriter *writer = vtkDataSetWriter::New();
writer->SetInput(pd);
writer->SetFileName("hrc.vtk");
writer->Write();
 */

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double * crossProduct(double * vector1, double * vector2) {
	double * returnVector = (double*)malloc(sizeof(double)*3);
	returnVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
	returnVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
	returnVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
	return returnVector;
}

double dotProduct(double * vector1, double * vector2) {
	return vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];
}

void normalize(double * vector) {
	double norm = sqrt( pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2) );
	vector[0] = vector[0]/norm;
	vector[1] = vector[1]/norm;
	vector[2] = vector[2]/norm;
}

Matrix getTransformMatrix(Camera camera) {
	Matrix cameraMatrix;
	Matrix imageMatrix;
	Matrix deviceMatrix;
	
	/*
	cerr << "position[0] = " << camera.position[0] << endl;
	cerr << "position[1] = " << camera.position[1] << endl;
	cerr << "position[2] = " << camera.position[2] << endl;
	cerr << "focus[0] = " << camera.focus[0] << endl;
	cerr << "focus[1] = " << camera.focus[1] << endl;
	cerr << "focus[2] = " << camera.focus[2] << endl;
	cerr << "up[0] = " << camera.up[0] << endl;
	cerr << "up[1] = " << camera.up[1] << endl;
	cerr << "up[2] = " << camera.up[2] << endl;
	*/
	
	double * v1;
	double * v2;
	double * v3 = (double*)malloc(sizeof(double)*3);
	double * t = (double*)malloc(sizeof(double)*3);
	
	for(int i = 0; i < 3; i++) {
		v3[i] = camera.position[i] - camera.focus[i];
		t[i] = 0.0 - camera.position[i];
	}
	
	//cerr << "got here" <<  endl;
	v1 = crossProduct(camera.up, v3);
	normalize(v1);
	
	//cerr << v1[0] << "," << v1[1] << "," << v1[2] << endl;
	v2 = crossProduct(v3, v1);
	normalize(v2);
	normalize(v3);
	
	//cerr << v2[0] << "," << v2[1] << "," << v2[2] << endl;
	//cerr << v3[0] << "," << v3[1] << "," << v3[2] << endl;
	//cerr << t[0] << "," << t[1] << "," << t[2] << endl;
	
	for(int i = 0; i < 3; i++) {
		cameraMatrix.A[i][0] = v1[i];
		cameraMatrix.A[i][1] = v2[i];
		cameraMatrix.A[i][2] = v3[i];
		cameraMatrix.A[i][3] = 0.0;
	}
	
	cameraMatrix.A[3][0] = dotProduct(t, v1);
	cameraMatrix.A[3][1] = dotProduct(t, v2);
	cameraMatrix.A[3][2] = dotProduct(t, v3);
	cameraMatrix.A[3][3] = 1.0;
	// camera matrix formed
	
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			imageMatrix.A[i][j] = 0.0;
			deviceMatrix.A[i][j] = 0.0;
		}
	}
	
	//cerr << "image matrix\n";
	//imageMatrix.Print(cerr);
	//cerr << "device matrix\n";
	//deviceMatrix.Print(cerr);
	
	imageMatrix.A[0][0] = 1.0 / tan(camera.angle/2.0);
	imageMatrix.A[1][1] = imageMatrix.A[0][0];
	imageMatrix.A[2][2] = (camera.far+camera.near)/(camera.far-camera.near);
	imageMatrix.A[2][3] = -1.0;
	imageMatrix.A[3][2] = 2*camera.far*camera.near/(camera.far-camera.near);
	// image matrix formed
	
	deviceMatrix.A[0][0] = 500.0;
	deviceMatrix.A[1][1] = 500.0;
	deviceMatrix.A[2][2] = 1.0;
	deviceMatrix.A[3][0] = 500.0;
	deviceMatrix.A[3][1] = 500.0;
	deviceMatrix.A[3][3] = 1.0;
	/*
	cerr << "camera matrix\n";
	cameraMatrix.Print(cerr);
	cerr << "image matrix\n";
	imageMatrix.Print(cerr);
	cerr << "device matrix\n";
	deviceMatrix.Print(cerr);
	*/
	Matrix compose1 = Matrix::ComposeMatrices(cameraMatrix, imageMatrix);
	//compose1.Print(cerr);
	Matrix finalMatrix = Matrix::ComposeMatrices(compose1, deviceMatrix);
	//finalMatrix.Print(cerr);
	return finalMatrix;
}

int main()
{  
	std::vector<Triangle> triangles = GetTriangles();
	LightingParameters lp;
	for(int h = 0; h < 4; h++) {
		vtkImageData *image = NewImage(1000, 1000);
		unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
		double *depthBuffer = (double*)malloc(sizeof(double)*1000*1000);
		int npixels = 1000*1000;
		for (int i = 0 ; i < npixels*3 ; i++)
			buffer[i] = 0;
		
		for (int i = 0; i <  npixels; i++)
			depthBuffer[i] = -1;
		
		
		
		
		Screen screen;
		screen.buffer = buffer;
		screen.depthBuffer = depthBuffer;
		screen.width = 1000;
		screen.height = 1000;
	
		Camera camera = GetCamera(h*250, 1000);
		
		Matrix finalMatrix = getTransformMatrix(camera);
		//cerr << "final matrix" << endl;
		//finalMatrix.Print(cerr);
		double * view = (double*)malloc(sizeof(double)*3);
		for(int i = 0; i < 3; i++) {
			view[i] = camera.position[i] - camera.focus[i];
		}
		normalize(view);
		//cerr << view[0] << "," << view[1] << "," << view[2] << endl;
		
		for(int i=0; i< triangles.size(); i++) {
			/*triangles[i].setId(i);
			if(i == 748956) {
				bool flat = triangles[i].arrangePoints();
				double splitX = triangles[i].getSplitX();
				cout << "(" << triangles[i].X[0] << "," << triangles[i].Y[0] << ") (" << triangles[i].X[1] << "," << triangles[i].Y[1] << ") (" <<
					triangles[i].X[2] << "," << triangles[i].Y[2] <<  ") flat: " << flat << ", split X:" << splitX << endl; 
			}*/
			double vertexIn1[4];
			double vertexIn2[4];
			double vertexIn3[4];
			double vertexOut1[4];
			double vertexOut2[4];
			double vertexOut3[4];
			vertexIn1[0] = triangles[i].X[0];
			vertexIn1[1] = triangles[i].Y[0];
			vertexIn1[2] = triangles[i].Z[0];
			vertexIn1[3] = 1.0;
			vertexIn2[0] = triangles[i].X[1];
			vertexIn2[1] = triangles[i].Y[1];
			vertexIn2[2] = triangles[i].Z[1];
			vertexIn2[3] = 1.0;
			vertexIn3[0] = triangles[i].X[2];
			vertexIn3[1] = triangles[i].Y[2];
			vertexIn3[2] = triangles[i].Z[2];
			vertexIn3[3] = 1.0;
			finalMatrix.TransformPoint(vertexIn1, vertexOut1);
			finalMatrix.TransformPoint(vertexIn2, vertexOut2);
			finalMatrix.TransformPoint(vertexIn3, vertexOut3);
			
			/*if(i < 100) {
				cerr << endl;
				cerr << vertexIn1[0] << "," << vertexIn1[1] << "," << vertexIn1[2] << " "
					<< vertexOut1[0]/vertexOut1[3] << "," << vertexOut1[1]/vertexOut1[3] << "," << vertexOut1[2]/vertexOut1[3] << endl;
				cerr << vertexIn2[0] << "," << vertexIn2[1] << "," << vertexIn2[2] << " "
					<< vertexOut2[0]/vertexOut2[3] << "," << vertexOut2[1]/vertexOut2[3] << "," << vertexOut2[2]/vertexOut2[3] << endl;
				cerr << vertexIn3[0] << "," << vertexIn3[1] << "," << vertexIn3[2] << " "
					<< vertexOut3[0]/vertexOut3[3] << "," << vertexOut3[1]/vertexOut3[3] << "," << vertexOut3[2]/vertexOut3[3] << endl;
			}*/
			
			Triangle transformed;
			transformed.X[0] = vertexOut1[0]/vertexOut1[3];
			transformed.X[1] = vertexOut2[0]/vertexOut2[3];
			transformed.X[2] = vertexOut3[0]/vertexOut3[3];
			transformed.Y[0] = vertexOut1[1]/vertexOut1[3];
			transformed.Y[1] = vertexOut2[1]/vertexOut2[3];
			transformed.Y[2] = vertexOut3[1]/vertexOut3[3];
			transformed.Z[0] = vertexOut1[2]/vertexOut1[3];
			transformed.Z[1] = vertexOut2[2]/vertexOut2[3];
			transformed.Z[2] = vertexOut3[2]/vertexOut3[3];
			
			for(int j = 0; j < 3; j++) {
				for(int k = 0; k < 3; k++) {
					transformed.colors[j][k] = triangles[i].colors[j][k];
					transformed.normals[j][k] = triangles[i].normals[j][k];
				}
			}
			//if(i < 10) {	
			//	cerr << transformed.colors[0][0] << "," << transformed.colors[0][1] << "," << transformed.colors[0][2] << endl;
			//}
			
			if(transformed.arrangePoints() == true) {
				// cout << "(" << transformed.X[0] << "," << transformed.Y[0] << ") (" << transformed.X[1] << "," << transformed.Y[1] << ") (" <<
				//	transformed.X[3] << "," << transformed.Y[3] <<  ")" << endl; 
				screen.paintTriangle(transformed, lp, view);
				// cout << "flat triangle" << endl;
			} else {
				//if(i < 100) {
				//	cout << "(" << transformed.X[0] << "," << transformed.Y[0] << ") (" << transformed.X[1] << "," << transformed.Y[1] << ") (" <<
				//		transformed.X[3] << "," << transformed.Y[3] <<  ")" << endl;
				//}
				
				double splitX;
				splitX = transformed.getSplitX();
				double splitZ;
				splitZ = linearInterpolateZ(transformed.getBottomPointX(), transformed.getBottomPointY(), transformed.getBottomPointZ(),
												transformed.getTopPointX(), transformed.getTopPointY(), transformed.getTopPointZ(),
												splitX, transformed.getMidPointY());
												
				double* splitColor = linearInterpolateColor(transformed.getBottomPointX(), transformed.getBottomPointY(), transformed.getBottomPointColor(),
												transformed.getTopPointX(), transformed.getTopPointY(), transformed.getTopPointColor(),
												splitX, transformed.getMidPointY());
												
				double* splitVector = linearInterpolateVector(transformed.getBottomPointX(), transformed.getBottomPointY(), transformed.getBottomPointVector(),
												transformed.getTopPointX(), transformed.getTopPointY(), transformed.getTopPointVector(),
												splitX, transformed.getMidPointY());
				
				Triangle topTriangle;
				Triangle bottomTriangle;
				topTriangle.X[0] = transformed.getTopPointX();
				topTriangle.Y[0] = transformed.getTopPointY();
				topTriangle.Z[0] = transformed.getTopPointZ();
				topTriangle.colors[0][0] = transformed.getTopPointColor()[0];
				topTriangle.colors[0][1] = transformed.getTopPointColor()[1];
				topTriangle.colors[0][2] = transformed.getTopPointColor()[2];
				topTriangle.normals[0][0] = transformed.getTopPointVector()[0];
				topTriangle.normals[0][1] = transformed.getTopPointVector()[1];
				topTriangle.normals[0][2] = transformed.getTopPointVector()[2];
				
				
				topTriangle.X[1] = transformed.getMidPointX();
				topTriangle.Y[1] = transformed.getMidPointY();
				topTriangle.Z[1] = transformed.getMidPointZ();
				topTriangle.colors[1][0] = transformed.getMidPointColor()[0];
				topTriangle.colors[1][1] = transformed.getMidPointColor()[1];
				topTriangle.colors[1][2] = transformed.getMidPointColor()[2];
				topTriangle.normals[1][0] = transformed.getMidPointVector()[0];
				topTriangle.normals[1][1] = transformed.getMidPointVector()[1];
				topTriangle.normals[1][2] = transformed.getMidPointVector()[2];
				
				topTriangle.X[2] = splitX;
				topTriangle.Y[2] = transformed.getMidPointY();
				topTriangle.Z[2] = splitZ;
				topTriangle.colors[2][0] = splitColor[0];
				topTriangle.colors[2][1] = splitColor[1];
				topTriangle.colors[2][2] = splitColor[2];
				topTriangle.normals[2][0] = splitVector[0];
				topTriangle.normals[2][1] = splitVector[1];
				topTriangle.normals[2][2] = splitVector[2];
				
				//topTriangle.setId(i);
				
				bottomTriangle.X[0] = transformed.getBottomPointX();
				bottomTriangle.Y[0] = transformed.getBottomPointY();
				bottomTriangle.Z[0] = transformed.getBottomPointZ();
				bottomTriangle.colors[0][0] = transformed.getBottomPointColor()[0];
				bottomTriangle.colors[0][1] = transformed.getBottomPointColor()[1];
				bottomTriangle.colors[0][2] = transformed.getBottomPointColor()[2];
				bottomTriangle.normals[0][0] = transformed.getBottomPointVector()[0];
				bottomTriangle.normals[0][1] = transformed.getBottomPointVector()[1];
				bottomTriangle.normals[0][2] = transformed.getBottomPointVector()[2];
				
				bottomTriangle.X[1] = transformed.getMidPointX();
				bottomTriangle.Y[1] = transformed.getMidPointY();
				bottomTriangle.Z[1] = transformed.getMidPointZ();
				bottomTriangle.colors[1][0] = transformed.getMidPointColor()[0];
				bottomTriangle.colors[1][1] = transformed.getMidPointColor()[1];
				bottomTriangle.colors[1][2] = transformed.getMidPointColor()[2];
				bottomTriangle.normals[1][0] = transformed.getMidPointVector()[0];
				bottomTriangle.normals[1][1] = transformed.getMidPointVector()[1];
				bottomTriangle.normals[1][2] = transformed.getMidPointVector()[2];
				
				bottomTriangle.X[2] = splitX;
				bottomTriangle.Y[2] = transformed.getMidPointY();
				bottomTriangle.Z[2] = splitZ;
				bottomTriangle.colors[2][0] = splitColor[0];
				bottomTriangle.colors[2][1] = splitColor[1];
				bottomTriangle.colors[2][2] = splitColor[2];
				bottomTriangle.normals[2][0] = splitVector[0];
				bottomTriangle.normals[2][1] = splitVector[1];
				bottomTriangle.normals[2][2] = splitVector[2];
				
				//bottomTriangle.setId(-1*i);
				
				topTriangle.arrangePoints();
				bottomTriangle.arrangePoints();
				//if(i < 10) {
				//	cerr << endl;
				//	cerr << topTriangle.colors[0][0] << "," << topTriangle.colors[0][1] << "," << topTriangle.colors[0][2] << endl;
				//	cerr << bottomTriangle.colors[0][0] << "," << bottomTriangle.colors[0][1] << "," << bottomTriangle.colors[0][2] << endl;
				//}
				/*if(i == 748956) {
					cout << "(" << topTriangle.X[0] << "," << topTriangle.Y[0] << ") (" << topTriangle.X[1] << "," << topTriangle.Y[1] << ") (" <<
						topTriangle.X[2] << "," << topTriangle.Y[2] <<  ")" << endl; 
					cout << "(" << bottomTriangle.X[0] << "," << bottomTriangle.Y[0] << ") (" << bottomTriangle.X[1] << "," << bottomTriangle.Y[1] << ") (" <<
						bottomTriangle.X[2] << "," << bottomTriangle.Y[2] <<  ")" << endl;
				}*/
				screen.paintTriangle(topTriangle, lp, view);
				screen.paintTriangle(bottomTriangle, lp, view);
			}
			/*if(i == 758956) {
				cout << "finished triangle 0" << endl;
			}*/
		}
		//cout << ceil441(255*1) << endl;
		char picName[16];
		sprintf(picName, "Bubbles%d", h);
		WriteImage(image, picName);
		free(depthBuffer);
	}
}