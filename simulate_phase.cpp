#include <iostream>
#include <fstream>
#include <cmath>

#define PI          3.141592653589793238462643383279502884L /* pi */

inline double sqr(double x) {return x*x;} //our own sqr(x) function

class Point 
{
    /*
     * The class Point is meant to ease using points in 3D cartesian 
     * coordinates, especially addition, substraction, reading and
     * writing points and vectors, it can be extended to include other
     * operations but that's all what I need now
     * */
public:
    Point(double X=0.0, double Y=0.0, double Z=0.0) :x{X}, y{Y}, z{Z} {}
    friend Point operator+(const Point &P1, const Point &P2);
    friend Point operator-(const Point &P1, const Point &P2);
    friend std::ostream& operator<<(std::ostream &out, const Point &point);
    friend std::istream& operator>>(std::istream &in, Point &point);
    double norm() { return sqrt(x*x + y*y + z*z);}
private:
    double x, y, z;
};

Point operator+(const Point &P1, const Point &P2)
{
    Point result;
    result.x = P1.x + P2.x;
    result.y = P1.y + P2.y;
    result.z = P1.z + P2.z;
    return result;
}

Point operator-(const Point &P1, const Point &P2)
{
    Point result;
    result.x = P1.x - P2.x;
    result.y = P1.y - P2.y;
    result.z = P1.z - P2.z;
    return result;
}

std::ostream& operator<<(std::ostream &out, const Point &point)
{
    out << "(" << point.x << " ," << point.y << " ," << point.z << ")";
    return out;
}

std::istream& operator>>(std::istream &in, Point &point)
{
    in >> point.x; in  >> point.y; in  >> point.z;
    return in;
}

struct Platform
{
    double gamma;
    double f_0;
    double *t;
    Point *pos;
    int npulses;
    int nsamples;
    double T_p;
};

long read_data(Platform *&plat, Point *&points, float *&amplitudes)
    /*
     * This function reads platform points and amplitudes from the files
     * platform.dat points.dat and amplitudes.dat files respectively and
     * put them in the input variables plat, points and amplitudes
     *
     * */
{
    using namespace std;

    //Allocating memory for platform
    //and reading data
    plat = new Platform;
    ifstream inf("platform.dat");
    inf >> plat->gamma;
    inf >> plat->f_0;
    inf >> plat->npulses;
    inf >> plat->nsamples;
    inf >> plat->T_p;
    plat->t = new double[plat->nsamples];
    plat->pos = new Point[plat->npulses];
    for (int i=0; i<(plat->nsamples); i++)
        inf >> plat->t[i];
    
    for (int i=0; i<plat->npulses; i++)
        inf >> plat->pos[i];

    //Allocating and reading points
    //assuming that first number is number of points read
    ifstream infP("points.dat");
    long npoints;
    infP >> npoints;
    points = new Point[npoints];
    for (int i=0; i<npoints; i++)
        infP >> points[i];

    //Allocating and reading amplitudes
    //number of amplitudes is same as npoints
    ifstream infA("amplitudes.dat");
    amplitudes = new float[npoints];
    for(int i=0; i<npoints; i++)
        infA >> amplitudes[i];
    
    return npoints;
}

void simulate_phase(const Platform* plat, Point* points, float *amplitudes, int npoints)
    /*
     * This function is meant to simulate SAR phase history data
     * converted to C++ from the RITSAR Python library
     * it is meant to be more efficient
     * */
{
    using namespace std;
    double c = 3.0e8;

    //initialize phs real and imaginary parts
    //to a [npulses,nsamples] array of 0
    float **phs_real, **phs_imag;
    phs_real = new float*[plat->npulses];
    phs_imag = new float*[plat->npulses];
    for (int i = 0; i<plat->npulses; i++) {
        phs_real[i] = new float[plat->nsamples]();
        phs_imag[i] = new float[plat->nsamples]();
    }

    //simulate pulses
    double *phase = new double[plat->nsamples];
    for (int i=0; i<plat->npulses; i++) {
        cout << "simulating pulse " << i+1 << std::endl;

        double R_0 = (plat->pos[i]).norm();

        for (int j = 0; j<npoints; j++){
            double R_t = (plat->pos[i] - points[j]).norm();
            double dr = R_t - R_0;
            for (int k=0; k<(plat->nsamples); k++){
                double tmp = (plat->t[k] - 2*dr/c)/plat->T_p;
                if(tmp<0.5 && tmp>-0.5){
                    phase[k] = PI*plat->gamma*sqr(2*dr/c) -
                        2*PI*(plat->f_0 + plat->gamma*plat->t[k])*2*dr/c;
                    phs_real[i][k] += amplitudes[j]*cos(phase[k]);
                    phs_imag[i][k] += amplitudes[j]*sin(phase[k]);
                }
            }
        }
    }
    delete[] phase;

    //write phases to files
    ofstream out_real("phase_real.dat");
    ofstream out_imag("phase_imag.dat");
    for (int i=0; i<plat->npulses;i++){
        out_real << phs_real[i][0];
        out_imag << phs_imag[i][0];
        for (int j=1;j<plat->nsamples; j++) {
            out_real << "," << phs_real[i][j];
            out_imag << "," << phs_imag[i][j];
        }
        out_real << ",";
        out_imag << ",";
    }

}
int main()
{
    using namespace std;

    Platform *plat;
    Point *points;
    float *amplitudes;

    long npoints;
    npoints = read_data(plat, points, amplitudes);
    
    simulate_phase(plat, points, amplitudes, npoints);
    
    return 0;
}
