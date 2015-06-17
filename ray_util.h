#pragma once
#ifndef _vector
#define _vector
#include <vector>
#endif

using namespace std;


class Vector {
public:
	Vector();
	Vector(double px, double py, double pz);
	void unit();
	double xc();
	double yc();
	double zc();
    Vector copy();
	Vector cross(Vector v);
	Vector add(Vector v);
	Vector sub(Vector v);
	double dot(Vector v);
private:
	double x;
	double y;
	double z;
};

class Vector4 {
public:
    Vector4();
    Vector4(double px, double py, double pz, double pw);
    double xc();
    double yc();
    double zc();
    double wc();
    double dot4(Vector4 v);
protected:
    double x;
    double y;
    double z;
    double w;
};

class matrix {
//0: translation, 1: scale, 2: rotation
//Warning: don't use the invmult on a matrix created by transposeInverse
public:
    matrix();
	matrix(double a, double b, double c, int mtype);
    matrix(double x1, double x2, double x3, double x4,
           double y1, double y2, double y3, double y4,
           double z1, double z2, double z3, double z4,
           double w1, double w2, double w3, double w4);
	Vector4 multiplyv(Vector4 v);
    Vector4 invmult(Vector4 v);
	void multiplym(matrix m);
    matrix transposeInverse();
    void printMatrix();
private:
    vector<Vector4> mtrx;
    vector<Vector4> inv;
};

class ray
{
public:
	ray(double xo, double yo, double zo, double dx, double dy, double dz);
    ray(Vector start, Vector d);
    ray();
	double xpos(double t);
	double ypos(double t);
	double zpos(double t);
    Vector start();
	Vector d();
private:
    Vector i;
	Vector v;
};

class material{
    public:
        material(double mka[], double mkd[], double mks[], double p, double mkr[]);
        material(){};
        material copy();
        double ka[3];
        double kd[3];
        double ks[3];
        double sp;
        double kr[3];
};

class shape
{
    public:
        virtual double intersect(ray r) {return 0.0;};
        virtual void transform(matrix m) {};
        virtual Vector getNormal(double xx, double yy, double zz){return Vector();};
        virtual void setMaterial(material* m);
        virtual int getID(){return id;};
        virtual void printPoints() {};
        virtual void printNormal(double xx, double yy, double zz) {};
        int id;
        material mat;
};

class sphere: public shape
{
    public:
        sphere(double xc, double yc, double zc, double r);
        Vector getNormal(double x, double y, double z);
        void transform(matrix m);
        double intersect(ray r);
        void printPoints();
        void printNormal(double xx, double yy, double zz);
        double x;
        double y;
        double z;
        double radius;
        matrix transformation;
};

class triangle: public shape
{
    public:
        triangle(double v_1 [3], double v_2 [3], double v_3 [3]);
        triangle(double v_1 [3], double n_1[3],double v_2 [3], double n_2[3],double v_3 [3],double n_3[3]);
        Vector getNormal(double x, double y, double z);
        double intersect(ray r);
        void transform(matrix m);
        void printPoints();
        void printNormal(double xx, double yy, double zz);
        vector<double> v1;
        vector<double> v2;
        vector<double> v3;
        vector<double> n1;
        vector<double> n2;
        vector<double> n3;
        bool interpolate;
};

double quadratic_solve(double a, double b, double c);

class light {
    public:
    vector<double> getIntensity(); 
    double rI;
	double gI;
	double bI;
};

class plight: public light {
    public:
        plight(double px,double py,double pz,double r, double g,double b,double f);
        plight(){};
        double x;
        double y;
        double z;
        double falloff;
};

class dlight: public light {
    public:
        dlight(double px,double py,double pz,double r, double g,double b);
        dlight(){};
        double x;
        double y;
        double z;

};

class alight: public light 
{
    public:
        alight(double ar,double ag,double ab);
        alight(){};

};

class camera {
    public:
        camera(double cxyz[], double cll[], double clr[], double cul[], double cur[]);
        camera(){};
        double* xyz;
        double* ll;
        double* lr;
        double* ul;
        double* ur;
};

class ray_container {
public:
    ray_container();
    ray_container(double r, double g, double b, ray ra);
    ray get_ray();
    vector<double> get_rgb();
    void set_rgb(double r, double g, double b);
private:
    vector<double> rgb;
    ray i;
};


