#include "ray_util.h"
#include <math.h>
#include <cmath>
#include "rowreduce.cpp"
#ifndef _vector
#define _vector
#include <vector>
#endif



#include <iostream>



using namespace std;
const double PI_rad = 3.141592653589793/180;
int num_shapes = 0;

//Vector class
Vector::Vector() {
	x = 0;
	y = 0;
	z = 0;
}

Vector::Vector(double px, double py, double pz) {
	x = px;
	y = py;
	z = pz;
}


void Vector::unit() {
	double length = pow(pow(x, 2) + pow(y, 2) + pow(z, 2), .5);
	x = x / length;
	y = y / length;
	z = z / length;
}


double Vector::xc() {
	return x;
}

double Vector::yc() {
	return y;
}

double Vector::zc() {
	return z;
}

Vector Vector::copy() {
    Vector retVal = Vector(x, y, z);
    return retVal;
}

Vector Vector::cross(Vector v) {
	return Vector(y*v.zc()-z*v.yc(), -(x*v.zc()-z*v.xc()), x*v.yc()-y*v.xc());
}

Vector Vector::add(Vector v) {
	return Vector(x+v.xc(), y+v.yc(), z+v.zc());
}

Vector Vector::sub(Vector v) {
	return Vector(x-v.xc(), y-v.yc(), z-v.zc());
}

double Vector::dot(Vector v) {
	return x * v.xc() + y * v.yc() + z * v.zc();
}

Vector4::Vector4() {
    x = 0;
    y = 0;
    z = 0;
    w = 0;
}

Vector4::Vector4(double px, double py, double pz, double pw) {
    x = px;
    y = py;
    z = pz;
    w = pw;
}

double Vector4::xc() {
    return x;
}

double Vector4::yc() {
    return y;
}

double Vector4::zc() {
    return z;
}

double Vector4::wc() {
    return w;
}

double Vector4::dot4(Vector4 v) {
    return x*v.xc() + y*v.yc() + z*v.zc() + w*v.wc();
}

matrix::matrix() {
	mtrx.clear();
	inv.clear();
    mtrx.push_back(Vector4(1.0, 0.0, 0.0, 0.0));
    mtrx.push_back(Vector4(0.0, 1.0, 0.0, 0.0));
    mtrx.push_back(Vector4(0.0, 0.0, 1.0, 0.0));
    mtrx.push_back(Vector4(0.0, 0.0, 0.0, 1.0));
    inv.push_back(Vector4(1.0, 0.0, 0.0, 0.0));
    inv.push_back(Vector4(0.0, 1.0, 0.0, 0.0));
    inv.push_back(Vector4(0.0, 0.0, 1.0, 0.0));
    inv.push_back(Vector4(0.0, 0.0, 0.0, 1.0));

}

matrix::matrix(double a, double b, double c, int mtype) {
	mtrx.clear();
	inv.clear();
    if (mtype == 0) { // translation 
        mtrx.push_back(Vector4(1.0, 0.0, 0.0, a));
        mtrx.push_back(Vector4(0.0, 1.0, 0.0, b));
        mtrx.push_back(Vector4(0.0, 0.0, 1.0, c));
        inv.push_back(Vector4(1.0, 0.0, 0.0, -a));
        inv.push_back(Vector4(0.0, 1.0, 0.0, -b));
        inv.push_back(Vector4(0.0, 0.0, 1.0, -c));
    } else if (mtype == 1) { // scale
        mtrx.push_back(Vector4(a, 0.0, 0.0, 0.0));
        mtrx.push_back(Vector4(0.0, b, 0.0, 0.0));
        mtrx.push_back(Vector4(0.0, 0.0, c, 0.0));
        inv.push_back(Vector4(1.0/a, 0.0, 0.0, 0.0));
        inv.push_back(Vector4(0.0, 1.0/b, 0.0, 0.0));
        inv.push_back(Vector4(0.0, 0.0, 1.0/c, 0.0));
    } else if (mtype == 2) { //rotation
        double theta = pow(pow(a, 2) + pow(b, 2) + pow(c, 2), .5)*PI_rad;
        Vector rotation = Vector(a, b, c);
        rotation.unit();
        double x = rotation.xc();
        double y = rotation.yc();
        double z = rotation.zc();
        mtrx.push_back(Vector4(pow(x, 2)+(pow(z, 2)+pow(y, 2))*cos(theta),
                         x*y-z*sin(theta)-x*y*cos(theta),
                         x*z+y*sin(theta)-x*z*cos(theta), 0.0));
        mtrx.push_back(Vector4(x*y+z*sin(theta)-x*y*cos(theta),
                        pow(y, 2)+(pow(x, 2)+pow(z, 2))*cos(theta),
                        y*z-x*sin(theta)-y*z*cos(theta), 0.0));
        mtrx.push_back(Vector4(x*z-y*sin(theta)-x*z*cos(theta),
                         y*z+x*sin(theta)-z*y*cos(theta),
                         pow(z, 2)+(pow(x, 2)+pow(y, 2))*cos(theta), 0.0));
        inv.push_back(Vector4(pow(x, 2)+(pow(z, 2)+pow(y, 2))*cos(-theta),
                         x*y-z*sin(-theta)-x*y*cos(-theta),
                         x*z+y*sin(-theta)-x*z*cos(-theta), 0.0));
        inv.push_back(Vector4(x*y+z*sin(-theta)-x*y*cos(-theta),
                        pow(y, 2)+(pow(x, 2)+pow(z, 2))*cos(-theta),
                        y*z-x*sin(-theta)-y*z*cos(-theta), 0.0));
        inv.push_back(Vector4(x*z-y*sin(-theta)-x*z*cos(-theta),
                         y*z+x*sin(-theta)-z*y*cos(-theta),
                         pow(z, 2)+(pow(x, 2)+pow(y, 2))*cos(-theta), 0.0));
        }
    mtrx.push_back(Vector4(0.0, 0.0, 0.0, 1.0));
    inv.push_back(Vector4(0.0, 0.0, 0.0, 1.0));
}

matrix::matrix(double x1, double x2, double x3, double x4,
               double y1, double y2, double y3, double y4,
               double z1, double z2, double z3, double z4,
               double w1, double w2, double w3, double w4) {
	mtrx.clear();
	inv.clear();
    mtrx.push_back(Vector4(x1, y1, z1, w1));
    mtrx.push_back(Vector4(x2, y2, z2, w2));
    mtrx.push_back(Vector4(x3, y3, z3, w3));
    mtrx.push_back(Vector4(x4, y4, z4, w4));

    for (int i=0; i<4; i++) {
        inv.push_back(Vector4());
    }
}

Vector4 matrix::multiplyv(Vector4 v) {
    return Vector4(mtrx[0].dot4(v), mtrx[1].dot4(v), mtrx[2].dot4(v), mtrx[3].dot4(v));
}

Vector4 matrix::invmult(Vector4 v) {
    return Vector4(inv[0].dot4(v), inv[1].dot4(v), inv[2].dot4(v), inv[3].dot4(v));
}

void matrix::multiplym(matrix m) {
    Vector4 a = Vector4(m.mtrx[0].xc(), m.mtrx[1].xc(), m.mtrx[2].xc(), m.mtrx[3].xc());
    Vector4 b = Vector4(m.mtrx[0].yc(), m.mtrx[1].yc(), m.mtrx[2].yc(), m.mtrx[3].yc());
    Vector4 c = Vector4(m.mtrx[0].zc(), m.mtrx[1].zc(), m.mtrx[2].zc(), m.mtrx[3].zc());
    Vector4 d = Vector4(m.mtrx[0].wc(), m.mtrx[1].wc(), m.mtrx[2].wc(), m.mtrx[3].wc());
    for (int i=0; i<4; i++) {
        mtrx[i] = Vector4(mtrx[i].dot4(a), mtrx[i].dot4(b), mtrx[i].dot4(c), mtrx[i].dot4(d));
    }
    a = Vector4(inv[0].xc(), inv[1].xc(), inv[2].xc(), inv[3].xc());
    b = Vector4(inv[0].yc(), inv[1].yc(), inv[2].yc(), inv[3].yc());
    c = Vector4(inv[0].zc(), inv[1].zc(), inv[2].zc(), inv[3].zc());
    d = Vector4(inv[0].wc(), inv[1].wc(), inv[2].wc(), inv[3].wc());
    for (int i=0; i<4; i++) {
        inv[i] = Vector4(m.inv[i].dot4(a), m.inv[i].dot4(b), m.inv[i].dot4(c), m.inv[i].dot4(d));
    }
}

matrix matrix::transposeInverse() {
    return matrix(inv[0].xc(), inv[0].yc(), inv[0].zc(), inv[0].wc(),
                  inv[1].xc(), inv[1].yc(), inv[1].zc(), inv[1].wc(),
                  inv[2].xc(), inv[2].yc(), inv[2].zc(), inv[2].wc(),
                  inv[3].xc(), inv[3].yc(), inv[3].zc(), inv[3].wc());
}

void matrix::printMatrix() {
	cout << "Printing matrix" << endl;
	for (int i=0; i<4; i++) {
		cout << mtrx[i].xc() << " " << mtrx[i].yc() << " " << mtrx[i].zc() << " " << mtrx[i].wc() << endl;
	}
	cout << endl;
	for (int i=0; i<4; i++) {
		cout << inv[i].xc() << " " << inv[i].yc() << " " << inv[i].zc() << " " << inv[i].wc() << endl;
	}
	cout << endl;
}

//ray class
ray::ray(double xo, double yo, double zo, double dx, double dy, double dz) {
	i = Vector(xo, yo, zo);
	v = Vector(dx, dy, dx);
}

ray::ray(Vector start, Vector d) {
    Vector o =  Vector(start.xc(), start.yc(), start.zc());
    Vector du = Vector(d.xc(), d.yc(), d.zc());
    i = o;
    v = du;
}

ray::ray() {
    i = Vector();
    v = i; 
}

double ray::xpos(double t) {
	return i.xc() + v.xc() * t;
}

double ray::ypos(double t) {
	return i.yc() + v.yc() * t;
}

double ray::zpos(double t) {
	return i.zc() + v.zc() * t;
}

Vector ray::start() {
    return i;
}

Vector ray::d() {
	return v;
}

//material class
material::material(double mka[], double mkd[], double mks[], double p, double mkr[]){
	for (int i=0; i<3; i++) {
		ka[i] = mka[i];
		kd[i] = mkd[i];
		ks[i] = mks[i];
		kr[i] = mkr[i];
	}
	sp = p;
}

material material::copy() {
	double ca[3];
	double cd[3];
	double cs[3];
	double csp = sp;
	double cr[3];
	for (int i=0; i<3; i++) {
		ca[i] = ka[i];
		cd[i] = kd[i];
		cs[i] = ks[i];
		cr[i] = kr[i];
	}
	return material(ca, cd, cs, csp, cr);
}

//shape super-class
void shape::setMaterial(material* m){
    mat = m->copy();
}

//sphere class
sphere::sphere(double xc, double yc, double zc, double r){
    x = xc;
    y = yc;
    z = zc;
    radius = r;
    transformation = matrix();
    id = num_shapes;
    num_shapes++;
}

double sphere::intersect(ray r){ 
    /*
    * returns the lowest positive t value or INFINITY 
    */
    Vector vi = r.start();
    Vector di = r.d();
    Vector4 start = Vector4(vi.xc(), vi.yc(), vi.zc(), 1.0);
    Vector4 d = Vector4(di.xc(), di.yc(), di.zc(), 0.0);
    start = transformation.invmult(start);
    d = transformation.invmult(d);
    ray r0 = ray(Vector(start.xc(), start.yc(), start.zc()), Vector(d.xc(), d.yc(), d.zc()));
    double a = pow(r0.d().xc(),2) + pow(r0.d().yc(),2) + pow(r0.d().zc(),2);
    double b = 2*r0.d().xc()*(r0.xpos(0) - x) + 2*r0.d().yc()*(r0.ypos(0) - y) + 2*(r0.d().zc()*(r0.zpos(0) - z));
    double c = pow(r0.xpos(0) -x,2) + pow(r0.ypos(0)-y,2) + pow(r0.zpos(0)-z,2) - pow(radius,2);
    return quadratic_solve(a,b,c);  
   
}

Vector sphere::getNormal(double xx, double yy, double zz) {
    Vector4 prep = Vector4(xx, yy, zz, 1.0);
    prep = transformation.invmult(prep);
    prep = Vector4(prep.xc()-x, prep.yc()-y, prep.zc()-z, 0.0);

    matrix transInv = transformation.transposeInverse();
    prep = transInv.multiplyv(prep);
    Vector result = Vector(prep.xc(), prep.yc(), prep.zc());
    result.unit();
    return result;
}

void sphere::transform(matrix m) {
    matrix copy_matrix = matrix();
    copy_matrix.multiplym(m);
    transformation = copy_matrix;
}

void sphere::printPoints() {
	cout << "Sphere center: " << x << ", " << y << ", " << z << endl;
}

void sphere::printNormal(double xx, double yy, double zz) {
	Vector v = getNormal(xx, yy, zz);
	cout << v.xc() << " " << v.yc() << " " << v.zc() << endl;
}

//triangle
triangle::triangle(double v_1[3], double v_2[3], double v_3[3]){
	v1.clear();
	v2.clear();
	v3.clear();
	n1.clear();
	n2.clear();
	n3.clear();
	for (int i=0; i<3; i++) {
		v1.push_back(v_1[i]);
		v2.push_back(v_2[i]);
		v3.push_back(v_3[i]);
	}
    Vector one_two = Vector(v2[0] -v1[0], v2[1] - v1[1], v2[2] - v1[2]);
    Vector one_three = Vector(v3[0] -v1[0], v3[1] - v1[1], v3[2] - v1[2]);
    Vector normalV = one_two.cross(one_three);
    normalV.unit();
    double normal[3] = {normalV.xc(), normalV.yc(), normalV.zc()};
    for (int i=0; i<3; i++) {
        n1.push_back(normal[i]);
        n2.push_back(normal[i]);
        n3.push_back(normal[i]);
    }
    interpolate = false;
    id = num_shapes;
    num_shapes++;
}



triangle::triangle(double v_1 [3], double n_1[3],double v_2 [3], double n_2[3],double v_3 [3],double n_3[3]){
	v1.clear();
	v2.clear();
	v3.clear();
	n1.clear();
	n2.clear();
	n3.clear();
	for (int i=0; i<3; i++) {
		v1.push_back(v_1[i]);
		v2.push_back(v_2[i]);
		v3.push_back(v_3[i]);
		n1.push_back(n_1[i]);
		n2.push_back(n_2[i]);
		n3.push_back(n_3[i]);
	}
    interpolate = true;
    id = num_shapes;
    num_shapes++;
}

void triangle::transform(matrix m){
    Vector4 one = Vector4(v1[0], v1[1], v1[2], 1.0);
    Vector4 two = Vector4(v2[0], v2[1], v2[2], 1.0);
    Vector4 thr = Vector4(v3[0], v3[1], v3[2], 1.0);
    Vector4 no1 = Vector4(n1[0], n1[1], n1[2], 0.0);
    Vector4 no2 = Vector4(n2[0], n2[1], n2[2], 0.0);
    Vector4 no3 = Vector4(n3[0], n3[1], n3[2], 0.0);
    one = m.multiplyv(one);
    two = m.multiplyv(two);
    thr = m.multiplyv(thr);
    matrix transinv = m.transposeInverse();
    no1 = transinv.multiplyv(no1);
    no2 = transinv.multiplyv(no2);
    no3 = transinv.multiplyv(no3);
    v1[0] = one.xc();
    v1[1] = one.yc();
    v1[2] = one.zc();
    v2[0] = two.xc();
    v2[1] = two.yc();
    v2[2] = two.zc();
    v3[0] = thr.xc();
    v3[1] = thr.yc();
    v3[2] = thr.zc();
    n1[0] = no1.xc();
    n1[1] = no1.yc();
    n1[2] = no1.zc();
    n2[0] = no2.xc();
    n2[1] = no2.yc();
    n2[2] = no2.zc();
    n3[0] = no3.xc();
    n3[1] = no3.yc();
    n3[2] = no3.zc();    
}

Vector triangle::getNormal(double x, double y, double z) {
    if(!interpolate){
        Vector return_val = Vector(n1[0],n1[1],n1[2]);
        return return_val;
    }
    double ret[3][4] = {{v1[0], v2[0], v3[0], x},
                        {v1[1], v2[1], v3[1], y},
                        {v1[2], v2[2], v3[2], z}};
    to_reduced_row_echelon_form(ret);
    double weight1 = ret[0][3];
    double weight2 = ret[1][3];
    double weight3 = ret[2][3];
    Vector one = Vector(n1[0]*weight1, n1[1]*weight1, n1[2]*weight1);
    Vector two = Vector(n2[0]*weight2, n2[1]*weight2, n2[2]*weight2);
    Vector thr = Vector(n3[0]*weight3, n3[1]*weight3, n3[2]*weight3);
    Vector result = (one.add(two)).add(thr);    
    result.unit();
    return result;
}

double triangle::intersect(ray r){
    /*
    * Möller–Trumbore intersection algorithm
    */
    double EPSILON = 0.000001; 
    Vector e1;  //Edge1
    Vector e2;  //Edge2
    //Vec3 P, Q, T;
    Vector P;
    Vector Q;
    Vector T;
    double det, inv_det, u, v;
    double t;
    
    //Find vectors for two edges sharing V1
    e1 = Vector(v2[0] - v1[0], v2[1]-v1[1], v2[2]-v1[2]);
    e2 = Vector(v3[0] - v1[0], v3[1]-v1[1], v3[2]-v1[2]);
    
    //Begin calculating determinant - also used to calculate u parameter
    P = r.d().cross(e2); 
    
    //if determinant is near zero, ray lies in plane of triangle
    det = e1.dot(P);
    if(det > -EPSILON && det < EPSILON) return INFINITY;
    inv_det = 1.0 / det;
 
    //calculate distance from V1 to ray origin
    T = Vector(r.xpos(0) - v1[0], r.ypos(0) - v1[1], r.zpos(0) - v1[2]);
    //Calculate u parameter and test bound
    u = T.dot(P) * inv_det;
    //The intersection lies outside of the triangle
    if(u < 0.f || u > 1.f) return INFINITY;
 
    //Prepare to test v parameter
    Q = T.cross(e1);
    //Calculate V parameter and test bound
    v = r.d().dot(Q) * inv_det;
    //The intersection lies outside of the triangle
    if(v < 0.f || u + v  > 1.f) return INFINITY;
 
    t = e2.dot(Q) * inv_det; 
    if(t > EPSILON) { //ray intersection
        return t;
    }
    return INFINITY;
}

void triangle::printPoints() {
	cout << "Vertex 1: " << v1[0] << ", " << v1[1] << ", " << v1[2] << endl;
	cout << "Vertex 2: " << v2[0] << ", " << v2[1] << ", " << v2[2] << endl;
	cout << "Vertex 3: " << v3[0] << ", " << v3[1] << ", " << v3[2] << endl;
}

void triangle::printNormal(double xx, double yy, double zz) {
	Vector v = getNormal(xx, yy, zz);
	cout << v.xc() << " " << v.yc() << " " << v.zc() << endl;

}

double quadratic_solve(double a, double b, double c){
    double values[2] = {INFINITY, INFINITY};
    double determinant = b*b - 4*a*c;
    if (determinant > 0) {  // two real values
        values[1] = (-b + pow(determinant,0.5)) / (2*a);
        values[0] = (-b - pow(determinant,0.5)) / (2*a);
    }
    else if(determinant == 0){ // one real value
        values[0] = (-b)/(2*a);
    }
    else if(determinant < 0){
        return INFINITY;
    }
    if(values[0] < 0){
        values[0] = INFINITY;
    }
    if(values[1] < 0){
        values[1] = INFINITY;
    }
    return fmin(values[0],values[1]);
}

vector<double> light::getIntensity(){
    vector<double> intensity;
    intensity.push_back(rI);
    intensity.push_back(gI);
    intensity.push_back(bI);
    return intensity;
}

plight::plight(double px,double py,double pz,double r, double g,double b,double f){
    x = px;
    y = py;
    z = pz;
    rI = r;
    gI = g;
    bI = b;
    falloff = f;
}

dlight::dlight(double px,double py,double pz,double r, double g,double b){
    x = px;
    y = py;
    z = pz;
    rI = r;
    gI = g;
    bI = b;
}

alight::alight(double ar, double ag, double ab){
    rI = ar;
    gI = ag;
    bI = ab;
}

camera::camera(double cxyz[], double cll[], double clr[], double cul[], double cur[]){
    xyz = cxyz;
    ll = cll;
    lr = clr;
    ul = cul;
    ur = cur;
}

ray_container::ray_container() {
    i = ray(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    double rgba[3] = {0.0, 0.0, 0.0};
    rgb.push_back(rgba[0]);
    rgb.push_back(rgba[1]);
    rgb.push_back(rgba[2]);
}

ray_container::ray_container(double r, double g, double b, ray ra) {
    i = ra;
    double rgba[3] = {r, g, b};
    rgb.push_back(rgba[0]);
    rgb.push_back(rgba[1]);
    rgb.push_back(rgba[2]);

}

ray ray_container::get_ray() {
    return i;
}

vector<double> ray_container::get_rgb() {
    return rgb;
}

void ray_container::set_rgb(double r, double g, double b) {
    double rgba[3] = {r, g, b};
    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;

}