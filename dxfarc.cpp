/**
 * \file main.cpp
 *\brief Reading 3D line and arc data from .dxf file
 *
 * it reads a .dxf file, determines start and end points of arcs in x, y, z and
 * write lines between those points as .dxf
 * the lines can be overlaid to original drawing for checking
 *
 * Basically it a check program for conversion from object coordinate system to world coordinate system
 *
 * \author rainer_jordan@<very, very warm>mail.com
 *
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date   August 2021
 */

#include <cstdlib>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>      // std::invalid_argument
#define M_PI 3.14159265

 /*!
 \def MINUS1
 \brief a -1 for size_t variable
	gives 0 when incremented
 */
#define MINUS1 0xffffffffffffffffULL

 /*!
 * \class _tube
 * \brief class handling one tube
 *
 */
class _line {
public: //just to avoid using structure 
	size_t PointIn;        ///< number of point at tube start (inlet)
	size_t PointOut;       ///< number of point at tube end (outlet)
	double Radius;
	std::string LayerInfo; ///< string containing the information of the layer
  // constructor

	_line() {
		PointIn = MINUS1; // number of point at tube start 
		PointOut = MINUS1; // number of point at tube end [1]
		Radius = 0.;
		LayerInfo = " ";
	}
};

/*!
* \class _point
* \brief data related to one point
*
*/
class _point {
public:
	double xCoord; ///< x Coordinate of point in mm
	double yCoord; ///< y Coordinate of point in mm
	double zCoord; ///< z Coordinate of point in mm
  //constructor

	_point() {
		xCoord = 0.; // x Coordinate of point in mm
		yCoord = 0.; // y Coordinate of point in mm
		zCoord = 0.; // z Coordinate of point in mm 
	}
};

using namespace std;
/*!
* \brief writes data for a line to .dxf file
*
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle (has to be unique number)
* \param layer layer name
* \param color color of line
* \param xStart x-coordinate of start point
* \param yStart y-coordinate of start point
* \param zStart z-coordinate of start point
* \param xEnd x-coordinate of end point
* \param yEnd y-coordinate of end point
* \param zEnd z-coordinate of end point
* \return int error code
*/
int DXFWriteLine(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
	double xStart, double yStart, double zStart, double xEnd, double yEnd, double zEnd) {

	outData << "  0\nLINE\n  5\n";
	outData << std::hex << ++handle << "\n100\nAcDbEntity\n  8\n" << layer << "\n 62\n";
	outData << std::dec << color;
	outData << "\n100\nAcDbLine\n 10\n" << xStart << "\n 20\n" << yStart << "\n 30\n" << zStart << "\n";
	outData << " 11\n" << xEnd << "\n 21\n" << yEnd << "\n 31\n" << zEnd << "\n";
	return 0;
}

/**
 * \brief write multiple line data to .dxf file
 *
 * \param outData file stream of .dxf file
 * \param Lines
 * \param Points
 * \return
 */
int DXFWriteShort(ostream& outData,
	vector<_line>& Lines,
	vector<_point>& Points
) {
	int error = 0;
	//using new handles 
	unsigned long long handle = 100uLL;
	unsigned short color = 6;
	outData << "999\nDXF2WC\n  0\nSECTION\n  2\nCLASSES\n  0\nENDSEC\n  0\nSECTION\n  2\nENTITIES\n";

	if (!Lines.empty()) {
		for (auto& iLine : Lines) {
			if (iLine.Radius > 1.) {
				error = DXFWriteLine(outData, handle, "0", color,
					Points[iLine.PointIn].xCoord,
					Points[iLine.PointIn].yCoord,
					Points[iLine.PointIn].zCoord,
					Points[iLine.PointOut].xCoord,
					Points[iLine.PointOut].yCoord,
					Points[iLine.PointOut].zCoord
				);
			}
		}
	}

	outData << "  0\nENDSEC\n  0\nSECTION\n  2\nOBJECTS\n  0\nENDSEC\n  0\nEOF\n";
	return error;
}

/**
 * \brief converts OCS data of arcs to WCS data.
 * 
 * \param OCSCenterX x-coordinate of center point in OCS
 * \param OCSCenterY y-coordinate of center point in OCS
 * \param OCSCenterZ z-coordinate of center point in OCS
 * \param Nx x-component of extrusion vector 
 * \param Ny y-component of extrusion vector
 * \param Nz z-component of extrusion vector
 * \param OCSRadius radius of arc
 * \param OCSAngle angle of point in OCS
 * \return 
 */
_point OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
	double Nx, double Ny, double Nz,
	double OCSRadius, double OCSAngle) {

	// Autocad help gives a limit of 1/64 for Nx and Ny to switch to different calculation
	// otherwise division by 0
	_point WCS;

	if (fabs(Nx) < 1. / 64. && fabs(Ny) < 1. / 64.) {
		double root = sqrt(1. - Ny * Ny);
		WCS.xCoord = OCSCenterX * Nz / root - OCSCenterY * Ny * Nx / root + OCSCenterZ * Nx;
		WCS.yCoord = OCSCenterY * root + OCSCenterZ * Ny;
		WCS.zCoord = -OCSCenterX * Nx / root - OCSCenterY * Ny * Nz / root + OCSCenterZ * Nz;
		if (OCSRadius < 1e-6) return WCS;
		double OCSx = cos(OCSAngle / 180. * M_PI) * OCSRadius;
		double OCSy = sin(OCSAngle / 180. * M_PI) * OCSRadius;
		WCS.xCoord += OCSx * Nz / root - OCSy * Ny * Nx / root;
		WCS.yCoord += OCSy * root;
		WCS.zCoord += -OCSx * Nx / root - OCSy * Ny * Nz / root;
	}
	else {
		double root = sqrt(1. - Nz * Nz);
		WCS.xCoord = -OCSCenterX * Ny / root - OCSCenterY * Nz * Nx / root + OCSCenterZ * Nx;
		WCS.yCoord = OCSCenterX * Nx / root - OCSCenterY * Nz * Ny / root + OCSCenterZ * Ny;
		WCS.zCoord = OCSCenterY * root + OCSCenterZ * Nz;
		if (OCSRadius < 1e-6) return WCS;
		double OCSx = cos(OCSAngle / 180. * M_PI) * OCSRadius;
		double OCSy = sin(OCSAngle / 180. * M_PI) * OCSRadius;
		WCS.xCoord += -OCSx * Ny / root - OCSy * Nz * Nx / root;
		WCS.yCoord += OCSx * Nx / root - OCSy * Nz * Ny / root;
		WCS.zCoord += OCSy * root;
	}
	return WCS;
}


/**
 * .
 *
 * \param mPt
 * \param px
 * \param py
 * \param pz
 * \param Points
 * \return
 */
size_t AddPointIfNotExist(size_t& mPt, double px, double py, double pz, vector<_point>& Points) {
	if (mPt != MINUS1) {
		for (size_t iPt = 0; iPt <= mPt; iPt++) {
			if (fabs(Points[iPt].xCoord - px) < 1. &&
				fabs(Points[iPt].yCoord - py) < 1. &&
				fabs(Points[iPt].zCoord - pz) < 1.) {
				return iPt;
			}
		}
	}
	// not found in Point vector -> add
	Points.push_back(_point());
	mPt++;
	Points[mPt].xCoord = px;
	Points[mPt].yCoord = py;
	Points[mPt].zCoord = pz;
	return mPt;
}
/**
 * .
 *
 * \param inData
 * \param NoLineDxf
 * \param mLn
 * \param error
 * \return
 */
double ReadCheckCoord(ifstream& inData, int NoLineDxf, size_t mLn, string error) {
	string row;
	string::size_type tailptr;
	double data;
	if (getline(inData, row).good()) {
		try {
			data = stod(row, &tailptr);
		}
		catch (const std::invalid_argument& ia) {
			cout << "Invalid " << error << "point of line " << mLn;
			cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file";
			cout << "\nPlease correct dxf-file" << endl;
			exit(1);
		}
	}
	else {
		cout << "Invalid " << error << "point of line " << mLn;
		cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
		cout << "\nPlease correct dxf-file" << endl;
		exit(1);
	}
	return data;
}

/**
 *
 * @param argc
 * @param argv can contain the project name, if empty user will be asked for file name
 * @return
 */
int main(int argc, char** argv) {
	/* Local variables */
	bool dxfError = false, entities = false;
	int NoLineDxf = 0, error;
	size_t mLn; // max index of Lines (max number of Lines - 1)
	size_t mPt; // max index of Points (number of points - 1)  
	string filename, readfile, writefile, row, LayerInf, token, Project;
	double p1x, p1y, p1z, p2x, p2y, p2z;
	double OCSCenterX, OCSCenterY, OCSCenterZ,
		OCSRadius, OCSStartAngle, OCSEndAngle,
		Nx, Ny, Nz;
	vector <_line> Lines;
	vector <_point> Points;
	vector <size_t> OrphanedPoints;
	_point center;
	ifstream inData;
	ofstream outData;
	ofstream prot;

	/**
	 * -# open the input file (.dxf file)
	 * -# read all lines until content is "EOF" (End Of File) or end of entities section
	 *    - there can be lines in some definition like blocks\n
	 *    those should not be processed\n
	 *    only the lines in ENTITIES section should be taken
	 * -# read data of new line or arc.
	 * -# no further checking performed
	 * -# check if p1 or p2 already exist\n
	 *		if not found in Point vector -> add
	 * -# write a new dxf-file with lines between arc start and end points\n
	 * 	The dxf file helps to identify the position of the arcs
	 */

	mPt = MINUS1;
	mLn = MINUS1;

	/*
	 * reading file name if not already passed as argument and
	 *  determine the name of input and output file
	 */
	if (argc > 1) {
		filename = argv[1];
	}
	else {
		cout << "\nFile name : ";
		cin >> filename;
	}
	readfile = filename + ".dxf";
	/*
	 * open the input file (.dxf file)
	 */
	inData.open(readfile.c_str());
	if (!inData.good()) {
		cout << "\n cannot open file " << readfile << endl;
		exit(1);
	}
	else {
		while (true) {
			if (getline(inData, row).good()) {
				NoLineDxf++;
				/*
				 * read all lines until content is "EOF" (End Of File)
				 */
				if ((entities && row == "ENDSEC") || row == "EOF") {
					break;
				}
				else if (row == "ENTITIES") {
					entities = true;
					/*
					 * there can be lines in some definition like blocks
					 * those should not be processed
					 * only the lines in ENTITIES section should be taken
					 */
				}
				else if (entities && row == "LINE") {
					//new tube 
					mLn++;
					Lines.push_back(_line());
					while (true) {
						NoLineDxf++;
						if (getline(inData, row).good()) {
							if (row == "  0") {
								break;
							}
							else if (row == "  8") {
								NoLineDxf++;
								if (getline(inData, LayerInf).good()) {
									Lines[mLn].LayerInfo = LayerInf;
								}
								else {
									cout << "Error reading layer name of line (tube) " << mLn;
									cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
								}
							}
							else if (row == " 10") {
								NoLineDxf++;
								p1x = ReadCheckCoord(inData, NoLineDxf, mLn, "x-coordinate of starting point");
							}
							else if (row == " 20") {
								NoLineDxf++;
								p1y = ReadCheckCoord(inData, NoLineDxf, mLn, "y-coordinate of starting point");
							}
							else if (row == " 30") {
								NoLineDxf++;
								p1z = ReadCheckCoord(inData, NoLineDxf, mLn, "z-coordinate of starting point");
							}
							else if (row == " 11") {
								NoLineDxf++;
								p2x = ReadCheckCoord(inData, NoLineDxf, mLn, "x-coordinate of end point");
							}
							else if (row == " 21") {
								NoLineDxf++;
								p2y = ReadCheckCoord(inData, NoLineDxf, mLn, "y-coordinate of end point");
							}
							else if (row == " 31") {
								NoLineDxf++;
								p2z = ReadCheckCoord(inData, NoLineDxf, mLn, "z-coordinate of end point");
							}
						}
						else {
							cout << "\nError reading line " << mLn;
							cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
							cout << "\nPlease correct dxf-file" << endl;
							exit(1);
						}
					}

					//set number of inlet and outlet point, check if they already exist
					Lines[mLn].PointIn = AddPointIfNotExist(mPt, p1x, p1y, p1z, Points);
					Lines[mLn].PointOut = AddPointIfNotExist(mPt, p2x, p2y, p2z, Points);

				} // end of this line
				else if (entities && row == "ARC") {
					//new arc
					mLn++;
					Lines.push_back(_line());
					// set standard values, will be overwritten if available
					Nx = 0.;
					Ny = 0.;
					Nz = 1.;
					while (true) {
						NoLineDxf++;
						if (getline(inData, row).good()) {
							if (row == "  0") {
								break;
							}
							else if (row == "  8") {
								NoLineDxf++;
								if (getline(inData, LayerInf).good()) {
									Lines[mLn].LayerInfo = LayerInf;
								}
								else {
									cout << "Error reading layer name of line (tube) " << mLn;
									cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
								}
							}
							else if (row == " 10") {
								NoLineDxf++;
								OCSCenterX = ReadCheckCoord(inData, NoLineDxf, mLn, "x-coordinate of OCS center point");
							}
							else if (row == " 20") {
								NoLineDxf++;
								OCSCenterY = ReadCheckCoord(inData, NoLineDxf, mLn, "y-coordinate of OCS center point");
							}
							else if (row == " 30") {
								NoLineDxf++;
								OCSCenterZ = ReadCheckCoord(inData, NoLineDxf, mLn, "z-coordinate of OCS center point");
							}
							else if (row == " 40") {
								NoLineDxf++;
								OCSRadius = ReadCheckCoord(inData, NoLineDxf, mLn, "Radius of arc");
							}
							else if (row == " 50") {
								NoLineDxf++;
								OCSStartAngle = ReadCheckCoord(inData, NoLineDxf, mLn, "Start angle of arc");
							}
							else if (row == " 51") {
								NoLineDxf++;
								OCSEndAngle = ReadCheckCoord(inData, NoLineDxf, mLn, "End angle of arc");
							}
							else if (row == "210") {
								NoLineDxf++;
								Nx = ReadCheckCoord(inData, NoLineDxf, mLn, "x-component of extrusion vector");
							}
							else if (row == "220") {
								NoLineDxf++;
								Ny = ReadCheckCoord(inData, NoLineDxf, mLn, "y-component of extrusion vector");
							}
							else if (row == "230") {
								NoLineDxf++;
								Nz = ReadCheckCoord(inData, NoLineDxf, mLn, "z-component of extrusion vector");
							}
						}
						else {
							cout << "\nError reading arc " << mLn;
							cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
							cout << "\nPlease correct dxf-file" << endl;
							exit(1);
						}
					}

					_point WCSPointIn = OCS2WCS(OCSCenterX, OCSCenterY, OCSCenterZ, Nx, Ny, Nz, OCSRadius, OCSStartAngle);
					_point WCSPointOut = OCS2WCS(OCSCenterX, OCSCenterY, OCSCenterZ, Nx, Ny, Nz, OCSRadius, OCSEndAngle);

					//set number of inlet and outlet point, check if they already exist
					Lines[mLn].PointIn = AddPointIfNotExist(mPt, WCSPointIn.xCoord, WCSPointIn.yCoord, WCSPointIn.zCoord, Points);
					Lines[mLn].PointOut = AddPointIfNotExist(mPt, WCSPointOut.xCoord, WCSPointOut.yCoord, WCSPointOut.zCoord, Points);
					Lines[mLn].Radius = OCSRadius;
				} // end of this tube (arc)
			} // valid line read
		} // while loop for line reading 
		inData.close();
	}//file opened correctly
	cout << "\n All " << mLn + 1 << " lines read" << endl;

	/* write the coordinates to dxf-file
	 * this dxf-file can be overlaid to the original  drawing
	 * the lines in this file are between start and end points of arcs */
	writefile = filename + "_Short.dxf";

	outData.open(writefile.c_str());

	if (!outData.good()) {
		cout << "\n cannot open file " << writefile << endl;
		exit(1);
	}
	else {
		error = DXFWriteShort(outData, Lines, Points);
		outData.close();
	}//file opened correctly


	return (EXIT_SUCCESS);
}

