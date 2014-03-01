///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef PIECEWISEBICUBIC_H
#define PIECEWISEBICUBIC_H

#include "useful.H"
#include "Field.H"

class PiecewiseBicubic : public Field {
private:
  int nx, ny;
  double x0, y0; // origin
  double lx, ly; // length along each dimension
  double dx, dy; // step along each dimension
  double x1, y1; // position of last node (x1=x0+lx-dx)
public:
  PiecewiseBicubic(const char* fileName, bool periodic0) {
    v0 = NULL;
    read(fileName, periodic0);
  }

  ~PiecewiseBicubic() {
  }

  PiecewiseBicubic(const PiecewiseBicubic& tab) {
    n = tab.n;
    nx = tab.nx;
    ny = tab.ny;
    x0 = tab.x0;
    y0 = tab.y0;
    lx = tab.lx;
    ly = tab.ly;
    dx = tab.dx;
    dy = tab.dy;
    x1 = tab.x1;
    y1 = tab.y1;

    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      err[i] = tab.err[i];
    }

    periodic = tab.periodic;
  }

  PiecewiseBicubic(double x00, double y00, double x10, double y10, int nx0, int ny0, bool periodic0) :
    Field(nx0*ny0,periodic0), nx(nx0), ny(ny0), x0(x00), y0(y00), lx(x10-x00), ly(y10-y00), dx(lx/nx), dy(ly/ny), x1(x10-dx), y1(y10-dy) {
  }

  void read(const String& fileName, bool periodic0) {
    periodic = periodic0;
    readDx(fileName.cs());
  }
  void readDx(const char* fileName) {
    // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      fprintf(stderr,"ERROR PiecewiseBicubic::readDx Couldn't open file %s.\n",fileName);
      exit(-1);
    }
    //printf("Reading dx file %s...\n", fileName);
    
    n = 0;
    nx = 0;
    ny = 0;
    int dummy;
    v0 = NULL;

    int size = 0;
    double x, y, z;
    char line[STRLEN];
    int p, nRead;
    int deltaCount = 0;
    Vector3 base[3];
    while (fgets(line, STRLEN, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
    
      if (String::isInt(line[0]) && size < n) {
	if (v0 == NULL) {
	  fprintf(stderr,"ERROR PiecewiseBicubic::readDx dx file `%s' does not specify size before having numerical entries.\n",fileName);
	  exit(-1);
	}

	// Read grid values.
	nRead = sscanf(line, "%lf %lf %lf", &x, &y, &z);
	switch(nRead) {
	case 1:
	  v0[size] = x;
	  size++;
	  if (n != size) {
	    fprintf(stderr,"ERROR PiecewiseBicubic::readDx Improperly formatted dx file %s.\n", fileName);
	    fprintf(stderr,"line `%s'\n", line);
	    exit(-1);
	  }
	  break;
	case 2:
	  v0[size] = x;
	  v0[size+1] = y;
	  size += 2;
	  if (n != size) {
	    fprintf(stderr,"ERROR PiecewiseBicubic::readDx Improperly formatted dx file %s.\n", fileName);
	    fprintf(stderr,"line `%s'\n", line);
	    exit(-1);
	  }
	  break;
	case 3:
	  if (n <= size+2) {
	    fprintf(stderr,"ERROR PiecewiseBicubic::readDx Too many lines in `%s'.\n", fileName);
	    fprintf(stderr,"line `%s'\n", line);
	    exit(-1);
	  }
	  v0[size] = x;
	  v0[size+1] = y;
	  v0[size+2] = z;
	  size += 3;
	  break;
	}
      } else if (len > 5) {
	// Read the grid parameters.
	char start[6];
	for (int i = 0; i < 5; i++) start[i] = line[i];
	start[5] = '\0';

	if(strcmp("origi", start) == 0) {
	  // Get an origin line.
	  p = String::firstSpace(line, STRLEN);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  x0 = x;
	  y0 = y;
	  //printf("Origin: %.12g %.12g %.12g\n", x, y, z);
	} else if(strcmp("delta", start) == 0) {
	  // Get a delta matrix line.
	  p = String::firstSpace(line, STRLEN);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  base[deltaCount] = Vector3(x, y, z);
	  //printf("Delta %d: %.12g %.12g %.12g\n", deltaCount, x, y, z);
	  if (deltaCount < 2) deltaCount = deltaCount + 1;
	} else if(strcmp("objec", start) == 0) {
	  //printf("%s", line);
	  // Get the system dimensions.
	  if (line[7] != '1') continue;
	  int read = sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &nx, &ny, &dummy);
	  //printf("Size: %d %d\n", nx, ny);
	  if (read == 3) {
	    n = nx*ny;
	    if (n == 0) {
	      fprintf(stderr,"ERROR PiecewiseBicubic::readDx Line `%s' of `%s' indicates that the grid has zero points.\n",line,fileName);
      	    }
	    v0 = new double[n];
	    err = new double[n];
	    // Zero the value and the error.
	    zero();
	  }
	}
      }
    }
    fclose(inp);

    dx = sqrt(base[0].x*base[0].x + base[0].y*base[0].y + base[0].z*base[0].z);
    dy = sqrt(base[1].x*base[1].x + base[1].y*base[1].y + base[1].z*base[1].z);
    x1 = x0 + (nx-1)*dx;
    y1 = y0 + (ny-1)*dy;
    lx = dx*nx;
    ly = dy*ny;

    if (n <= 0 || n != size) {
      fprintf(stderr,"ERROR PiecewiseBicubic::readDx Improperly formatted dx file %s.\n",fileName);
      fprintf(stderr,"declared size: %d, items: %d\n", n, size);
      fprintf(stderr,"first value: %10g, final value: %.10g\n", v0[0], v0[size-1]);
      exit(-1);
    }
  }  
  

  // Writes the grid as a file in the dx format.
  bool write(const String& fileName) const {
    Vector3 origin(x0, y0, 0.0);
    Matrix3 basis(dx, dy, 1.0);
    int nz = 1;

    // Open the file.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      printf("ERROR PiecewiseBicubic::write Couldn't open file %s.\n",fileName.cs());
      return false;
    }

    // Write the header.
    fprintf(out, "# %s\n", fileName.cs());
    fprintf(out, "object 1 class gridpositions counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "origin %.12g %.12g %.12g\n", origin.x, origin.y, origin.z);
    fprintf(out, "delta %.12g %.12g %.12g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "delta %.12g %.12g %.12g\n", basis.exy, basis.eyy, basis.ezy);
    fprintf(out, "delta %.12g %.12g %.12g\n", basis.exz, basis.eyz, basis.ezz);
    fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "object 3 class array type double rank 0 items %d data follows\n", n);
    
    // Write the data.
    int penultima = 3*(n/3);
    int mod = n - penultima;

    int i;
    for (i = 0; i < penultima; i+=3) {
      fprintf(out, "%.12g %.12g %.12g\n", v0[i], v0[i+1], v0[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.12g\n", v0[n-1]);
    } else if (mod == 2) {
      fprintf(out, "%.12g %.12g\n", v0[n-2], v0[n-1]);
    }
    fclose(out);
    return true;
  }

  // Writes the grid as a file in the dx format.
  bool writeArray(const String& filePrefix) const {
    String outFileV(filePrefix);
    outFileV.add(".v.dat");
    String outFileX(filePrefix);
    outFileX.add(".x.dat");
    String outFileY(filePrefix);
    outFileY.add(".y.dat");

    // Open the file.
    FILE* outV = fopen(outFileV.cs(),"w");
    FILE* outX = fopen(outFileX.cs(),"w");
    FILE* outY = fopen(outFileY.cs(),"w");
    if (outX == NULL || outY == NULL || outV == NULL ) {
      printf("ERROR PiecewiseBicubic::writeArray Couldn't open file %s.\n",outFileV.cs());
      return false;
    }

    // Write the positions and values.
    for (int ix = 0; ix < nx; ix++) {
      const double x = getPosX(ix);

      // Write the first item on the line.
      const double y0 = getPosY(0);
      const double v0 = get(ix*ny);
      fprintf(outX, "%.5g", x);
      fprintf(outY, "%.5g", y0);
      fprintf(outV, "%.5g", v0);

      // Write the rest of the line.
      for (int iy = 1; iy < ny; iy++) {
	const double y = getPosY(iy);
	const double v = get(iy + ix*ny);
	fprintf(outX, " %.5g", x);
	fprintf(outY, " %.5g", y);
	fprintf(outV, " %.5g", v);
      }
      fprintf(outX, "\n");
      fprintf(outY, "\n");
      fprintf(outV, "\n");
    }

    fclose(outV);
    fclose(outX);
    fclose(outY);
    return true;
  }

  PiecewiseBicubic& operator=(const PiecewiseBicubic& tab) {
    // Resize the arrays if necessary.
    if (n != tab.n) {
      delete[] v0;
      delete[] err;
      n = tab.n;
      v0 = new double[n];
      err = new double[n];
    }

    // Copy the values.
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      err[i] = tab.err[i];
    }
    dx = tab.dx;
    dy = tab.dy;
    nx = tab.nx;
    ny = tab.ny;
    x0 = tab.x0;
    y0 = tab.y0;
    x1 = tab.x1;
    y1 = tab.y1;
    lx = tab.lx;
    ly = tab.ly;
    periodic = tab.periodic;
    
    return (*this);
  }

  void scale(double s) {
    for (int i = 0; i < n; i++) v0[i] *= s;
  }
  void zero() {
    for (int i = 0; i < n; i++) {
      v0[i] = 0.0;
      err[i] = 0.0;
    }
  }

  // Get and set properties.
  double getX0() const { return x0; }
  double getY0() const { return y0; }
  double getX1() const { return x1; }
  double getY1() const { return y1; }
  double getDx() const { return dx; }
  double getDy() const { return dy; }
  int getNx() const { return nx; }
  int getNy() const { return ny; }

  double getBound0(int dir) const {
    switch(dir) {
    case 1:
      return x0;
    case 0:
    default:
      return y0;
    }
  }
  double getBound1(int dir) const {
    switch(dir) {
    case 1:
      return x1;
    case 0:
    default:
      return y1;
    }
  }

  // Things that must be implemented.
  bool set(int j, double v) {
    if (j < 0 || j >= n) return false;
    v0[j] = v;
    return true;
  }

  
  IndexList neighborsMinimal(int j) const {
    IndexList ret;
    int jy = j%ny;
    int jx = j/ny;
    
    if (periodic) {
      for (int ix = -2; ix <= 1; ix++) {
	int kx = wrap(jx+ix, nx);

	for (int iy = -2; iy <= 1; iy++) {
	  int ky = wrap(jy+iy, ny);
	  ret.add(ky + kx*ny);
	}
      } 

    } else {
      
      for (int ix = -2; ix <= 1; ix++) {
	int kx = jx+ix;
	if (kx < 0 || kx >= nx) continue;
	for (int iy = -2; iy <= 1; iy++) {
	  int ky = jy+iy;
	  if (ky < 0 || ky >= ny) continue;
	  ret.add(ky + kx*ny);
	}
      } 

    }
    
    return ret;
  }

  IndexList neighbors(int j) const {
    IndexList ret;
    int jy = j%ny;
    int jx = j/ny;
    
    if (periodic) {
      for (int ix = -2; ix <= 2; ix++) {
	int kx = wrap(jx+ix, nx);

	for (int iy = -2; iy <= 2; iy++) {
	  int ky = wrap(jy+iy, ny);
	  ret.add(ky + kx*ny);
	}
      } 

    } else {
      
      for (int ix = -2; ix <= 2; ix++) {
	int kx = jx+ix;
	if (kx < 0 || kx >= nx) continue;
	for (int iy = -2; iy <= 2; iy++) {
	  int ky = jy+iy;
	  if (ky < 0 || ky >= ny) continue;
	  ret.add(ky + kx*ny);
	}
      } 

    }
    
    return ret;
  }


  int nearestNode(const double* pos) const {
    int ix = wrap(int(floor((pos[0]-x0)/dx+0.5)), nx);
    int iy = wrap(int(floor((pos[1]-y0)/dy+0.5)), ny);
    return index(ix,iy);
  }
  // The nearest node with all coordinates less than those of pos.
  int getNode(const double* pos) const {
    int ix = wrap(int(floor((pos[0]-x0)/dx)), nx);
    int iy = wrap(int(floor((pos[1]-y0)/dy)), ny);
    return index(ix,iy);
  }
  inline int getNode(double x, double y) const {
    int ix = wrap(int(floor((x-x0)/dx)), nx);
    int iy = wrap(int(floor((y-y0)/dy)), ny);
    return index(ix,iy);
  }
  inline int getNodePos(double x, double y, double& hx, double& hy) const {
    int ix = wrap(int(floor((x-x0)/dx)), nx);
    int iy = wrap(int(floor((y-y0)/dy)), ny);
    // The position of the home node.
    hx = getPosX(ix);
    hy = getPosY(iy);
    return iy + ny*ix;
  }
  inline double getPosX(int ix) const { return ix*dx + x0; }
  inline double getPosY(int iy) const { return iy*dy + y0; }
  inline int index(int ix, int iy) const { return iy + ny*ix; }
  inline int indexX(int j) const { return j/ny; }
  inline int indexY(int j) const { return j%ny; }

  // Return the index of the previous node along a dimension.
  // Return -1 for none.
  int prevIndex(int j, int dim) const {
    // We don't check that j is valid.
    int ix = indexX(j);
    int iy = indexY(j);
    int ix0, iy0;
    
    switch(dim) {
    case 0:
      ix0 = ix-1;
      iy0 = iy;
      break;
    case 1:
      ix0 = ix;
      iy0 = iy-1;
      break;
    default:
      return -1;
    }

    if (periodic) return index(wrap(ix0,nx),wrap(iy0,ny));
    else {
      if (ix0 < 0) return -1;
      if (iy0 < 0) return -1;
      return index(ix0,iy0);
    }
  }
  // Return the index of the next node along a dimension.
  // Return -1 for none.
  int nextIndex(int j, int dim) const {
    // We don't check that j is valid.
    int ix = indexX(j);
    int iy = indexY(j);
    int ix1, iy1;

    switch(dim) {
    case 0:
      ix1 = ix+1;
      iy1 = iy;
      break;
    case 1:
      ix1 = ix;
      iy1 = iy+1;
      break;
    default:
      return -1;
    }

    if (periodic) return index(wrap(ix1,nx),wrap(iy1,ny));
    else {
      if (ix1 >= nx) return -1;
      if (iy1 >= ny) return -1;
      return index(ix1,iy1);
    }
  }
  // Grid spacing.
  double spacing(int j, int dim) const {
    switch(dim) {
    case 0:
      return dx;
    case 1:
      return dy;
    default:
      return 0.0;
    }
  }

  // Map. Don't forget to delete the returned value.
  // Map the values of 'srcVal' to the nodes of this.
  Field* map(const Field* srcVal) const {
    const PiecewiseBicubic* src = dynamic_cast<const PiecewiseBicubic*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR PiecewiseBicubic::map() requires fields of the same type.\n");
      exit(-1);
    }
    
    PiecewiseBicubic* ret = new PiecewiseBicubic(*this); 
    double r[2];
    for (int i = 0; i < n; i++) {
      r[0] = getPosX(indexX(i));
      r[1] = getPosY(indexY(i));
      double v = src->computeVal(r[0],r[1]);
      int near = src->nearestNode(r);
      ret->v0[i] = v;
      ret->err[i] = src->err[near];
    }
    return ret;
  }
  
  // Check that 'srcVal' is at least as big as '*this'.
  bool spannedBy(const Field* srcVal) const {
    const PiecewiseBicubic* src = dynamic_cast<const PiecewiseBicubic*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR PiecewiseBicubic::spannedBy() requires fields of the same type.\n");
      exit(-1);
    }

    if (x0 + dx < src->x0) return false;
    if (x1 - dx > src->x1) return false;
    if (y0 + dy < src->y0) return false;
    if (y1 - dy > src->y1) return false;
    return true;
  }

  bool writeIntegral(const String& fileName, double scale) const {
    return true;
  }
  // The log of the field.
  bool writeLog(const String& fileName, double scale) const {
    return true;
  }

  // Wraps floating point positions depending on the periodicity.
  inline double wrapX(double x) const {
    if (periodic) return x0 + wrapReal(x-x0,lx);
    if (x < x0) return x0;
    if (x > x1) return x1;
    return x;
  }
  inline double wrapY(double y) const {
    if (periodic) return y0 + wrapReal(y-y0,ly);
    if (y < y0) return y0;
    if (y > y1) return y1;
    return y;
  }

  double computeValZeroOrder(double x, double y) const {
    double g2[4][4];
    double wx, wy;
    double v = prepInterp(x, y, wx, wy, g2);

    return v;
  }

  // Do the interpolation.
  double computeVal(double x, double y) const {
    double g2[4][4];
    double wx, wy;
    prepInterp(x, y, wx, wy, g2);

    double wx2 = wx*wx;
    double wx3 = wx*wx2;
    double wy2 = wy*wy;
    double wy3 = wy*wy2;
    double a0, a1, a2, a3;

    // Mix along x.
    double g3[4];
    for (int iy = 0; iy < 4; iy++) {
      a3 = 0.5*(-g2[0][iy] + 3*g2[1][iy] - 3*g2[2][iy] + g2[3][iy]);
      a2 = 0.5*(2*g2[0][iy] - 5*g2[1][iy] + 4*g2[2][iy] - g2[3][iy]);
      a1 = 0.5*(-g2[0][iy] + g2[2][iy]);
      a0 = g2[1][iy];
   
      g3[iy] = a3*wx3 + a2*wx2 + a1*wx + a0;
    }

    // Mix along y.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return a3*wy3 + a2*wy2 + a1*wy + a0;
  }

  double computeGrad(int dir, double x, double y) const {
    switch(dir) {
    case 1:
	return computeGradY(x,y);
    case 0:
    default:
	return computeGradX(x,y);
    }
  }

  // Do the interpolation.
  double computeGradX(double x, double y) const {
    double g2[4][4];
    double wx, wy;
    prepInterp(x, y, wx, wy, g2);

    double wx2 = wx*wx;
    //double wx3 = wx*wx2;
    double wy2 = wy*wy;
    double wy3 = wy*wy2;
    double a0, a1, a2, a3;

    // Mix along x.
    double g3[4];
    for (int iy = 0; iy < 4; iy++) {
      a3 = 0.5*(-g2[0][iy] + 3*g2[1][iy] - 3*g2[2][iy] + g2[3][iy]);
      a2 = 0.5*(2*g2[0][iy] - 5*g2[1][iy] + 4*g2[2][iy] - g2[3][iy]);
      a1 = 0.5*(-g2[0][iy] + g2[2][iy]);
      a0 = g2[1][iy];
   
      //g3[iy] = a3*wx3 + a2*wx2 + a1*wx + a0;
      g3[iy] = 3.0*a3*wx2 + 2.0*a2*wx + a1;
    }

    // Mix along y.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return (a3*wy3 + a2*wy2 + a1*wy + a0)/dx;
  }

  // Do the interpolation.
  double computeGradY(double x, double y) const {
    double g2[4][4];
    double wx, wy;
    prepInterp(x, y, wx, wy, g2);

    double wx2 = wx*wx;
    double wx3 = wx*wx2;
    double wy2 = wy*wy;
    //double wy3 = wy*wy2;
    double a0, a1, a2, a3;

    // Mix along x.
    double g3[4];
    for (int iy = 0; iy < 4; iy++) {
      a3 = 0.5*(-g2[0][iy] + 3*g2[1][iy] - 3*g2[2][iy] + g2[3][iy]);
      a2 = 0.5*(2*g2[0][iy] - 5*g2[1][iy] + 4*g2[2][iy] - g2[3][iy]);
      a1 = 0.5*(-g2[0][iy] + g2[2][iy]);
      a0 = g2[1][iy];
   
      g3[iy] = a3*wx3 + a2*wx2 + a1*wx + a0;
    }

    // Mix along y.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    //return a3*wy3 + a2*wy2 + a1*wy + a0;
    return (3.0*a3*wy2 + 2.0*a2*wy + a1)/dy;
  }

private:
  inline double prepInterp(double x, double y, double& wx, double& wy, double val[4][4]) const {
    double rx = wrapX(x);
    double ry = wrapY(y);
    int indX = int(floor((rx-x0)/dx));
    int indY = int(floor((ry-y0)/dy));
    double homeX = indX*dx + x0;
    double homeY = indY*dy + y0;

    // The reduced coordinates.
    wx = (rx - homeX)/dx;
    wy = (ry - homeY)/dy;

    // Get the values at the neighbors.
    for (int ix = 0; ix < 4; ix++) {
      int jx = wrap(ix + indX - 1, nx);
      for (int iy = 0; iy < 4; iy++) {
	int jy = wrap(iy + indY - 1, ny);
	val[ix][iy] = v0[index(jx,jy)];
      }
    }

    return wx;
  }


};
#endif
