#include <iostream>
#include <thread>
#include <fstream>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ray_util.h"
#include "lodepng.h"
#ifndef _vector_
#define _vector_
#include <vector>
#endif

using namespace std;
inline double sqr(float x) { return x*x; }

/* Scene components */
double EMPTY_ARRAY[3] = {0,0,0};
vector<shape*> shapes;
vector<plight> plights;
vector<dlight> dlights;
alight al;
int has_camera;
camera cam;

material ma = material(EMPTY_ARRAY,EMPTY_ARRAY,EMPTY_ARRAY,0,EMPTY_ARRAY);
material* m = &ma;
int last = 0;
string objFile;
vector<vector<double> > vertices;
vector<vector<double> > normals;
matrix trans = matrix();
matrix tri_trans = matrix();

vector<double> diffu(double lv[], double nv[], double rgb[],double kd[]) {
    double l_n = fmax(0.0, lv[0]*nv[0] + lv[1]*nv[1] + lv[2]*nv[2]);
    vector<double> diff;
    diff.push_back(kd[0]*rgb[0]*l_n);
    diff.push_back(kd[1]*rgb[1]*l_n);
    diff.push_back(kd[2]*rgb[2]*l_n);
    return diff;
}

vector<double> ambien(double rgb[], double ka[]) {
    vector<double> ambi;
    ambi.push_back(ka[0] * rgb[0]);
    ambi.push_back(ka[1] * rgb[1]); 
    ambi.push_back(ka[2] * rgb[2]);
    return ambi;
}

vector<double> specul(double lv[], double nv[], double rgb[], double ks[],double sp, double v[]) {
    double l_n = lv[0]*nv[0] + lv[1]*nv[1] + lv[2]*nv[2];              
    double new_normal[] = {2*l_n*nv[0], 2*l_n*nv[1], 2*l_n*nv[2]};
    double r[] = {new_normal[0] - lv[0], new_normal[1]-lv[1], new_normal[2]-lv[2]};
    double r_v = pow(fmax(0.0, r[0]*v[0] + r[1]*v[1] + r[2]*v[2]),sp);
    vector<double> spec;
    spec.push_back(ks[0]*rgb[0]*r_v);
    spec.push_back(ks[1]*rgb[1]*r_v);
    spec.push_back(ks[2]*rgb[2]*r_v);
    return spec;
}

vector<double> phong(ray in, double t, shape* s){
    vector<double> diffuse;
    vector<double> ambient;
    vector<double> specular;

    for(int i = 0; i < 3; i++){
        diffuse.push_back(0);
        ambient.push_back(0);
        specular.push_back(0);
    }

    if(t != INFINITY){
        double xyz[3] = {in.xpos(t), in.ypos(t), in.zpos(t)};
        shape* sph = s;
        /*
         * ------------------
         * POINT       LIGHTS
         * ------------------
         */
        for(std::vector<plight>::iterator it = plights.begin(); it != plights.end(); ++it) {
            plight cur = *it;
            Vector nurmal = (*sph).getNormal(xyz[0],xyz[1],xyz[2]);
            double normal[] = {nurmal.xc(),nurmal.yc(),nurmal.zc()};
            double plrgb[3] = {0,0,0};
            bool intersect = false;
            Vector start = Vector(cur.x,cur.y,cur.z);
            Vector direction = Vector(xyz[0],xyz[1],xyz[2]).sub(start);
            ray to_obj = ray(start,direction);
            double tuh = (*s).intersect(to_obj);
            for(std::vector<shape*>::iterator iter = shapes.begin(); iter != shapes.end(); ++iter){
                shape* shapeCheck = *iter;
                if(shapeCheck->intersect(to_obj) < tuh && shapeCheck->getID() != (*s).getID()) {
                    intersect = true;
                    break;
                }
            }
            if(!intersect){
                double dist = sqrt(sqr(cur.x-xyz[0]) + sqr(cur.y-xyz[1]) + sqr(cur.z - xyz[2]));
                if(cur.falloff == 0){
                    plrgb[0] = cur.rI;
                    plrgb[1] = cur.gI;
                    plrgb[2] = cur.bI;
                }
                else if(cur.falloff == 1){
                    plrgb[0] = cur.rI/dist;
                    plrgb[1] = cur.gI/dist;
                    plrgb[2] = cur.bI/dist;

                }
                else if(cur.falloff == 2){
                    plrgb[0] = cur.rI/(dist*dist);
                    plrgb[1] = cur.gI/(dist*dist);
                    plrgb[2] = cur.bI/(dist*dist); 
                }
            }
            /*
             * DIFFUSE DIFFUSIVE DIFFUSE
             */
            Vector alpha_d = Vector(cur.x, cur.y, cur.z);
            Vector beta_d = Vector(xyz[0], xyz[1], xyz[2]);
            alpha_d = alpha_d.sub(beta_d);
            alpha_d.unit();
            double plv[3] = {alpha_d.xc(), alpha_d.yc(), alpha_d.zc()}; // l vector
            vector<double> diff = diffu(plv, normal, plrgb,(*sph).mat.kd);
            for (int i=0; i<3; i++) {
                diffuse[i] = diffuse[i] + diff[i];
            }
            /*
             * SPECULATIVE SPECULATION
             */
            Vector v = Vector().sub(in.d()); 
            v.unit();
            double view[3] = {v.xc(),v.yc(),v.zc()};
            vector<double> spec = specul(plv, normal, plrgb, (*sph).mat.ks,(*sph).mat.sp, view);
            for (int i=0; i<3; i++) {
                specular[i] = specular[i] + spec[i];
            }
        }

        /*
         * ------------------
         * DIRECTIONAL LIGHTS
         * ------------------
         */
        for(std::vector<dlight>::iterator it = dlights.begin(); it != dlights.end(); ++it) {
            dlight cur = *it;
            Vector nurmal = (*sph).getNormal(xyz[0],xyz[1],xyz[2]);
            double normal[] = {nurmal.xc(),nurmal.yc(),nurmal.zc()};
            double dlrgb[3] = {0, 0, 0};
            bool intersect = false;

            ray to_light = ray(xyz[0],xyz[1],xyz[2],-cur.x, -cur.y, -cur.z);
            for(std::vector<shape*>::iterator iter = shapes.begin(); iter != shapes.end(); ++iter){
                shape* shapeCheck = *iter; 
                if(shapeCheck->intersect(to_light)!= INFINITY && shapeCheck->getID() != (*s).getID()){
                    intersect = true;
                    break;
                }
            }
            if(!intersect){
                dlrgb[0] = cur.rI;
                dlrgb[1] = cur.gI;
                dlrgb[2] = cur.bI;

            }

            /*
             * DIFFUSE DIFFUSIVE DIFFUSE
             */
            double normalizeFactor = sqrt(sqr(cur.x) + sqr(cur.y) + sqr(cur.z));
            double dl[3] = {-cur.x/normalizeFactor, -cur.y/normalizeFactor, -cur.z/normalizeFactor};
            vector<double> diff = diffu(dl, normal, dlrgb, (*sph).mat.kd);
            for (int i=0; i<3; i++) {
                diffuse[i] = diffuse[i] + diff[i];
            }

            /*
             * SPECULATIVE SPECULATION
             */

            Vector v = Vector().sub(in.d()); 
            v.unit();
            double view[3] = {v.xc(),v.yc(),v.zc()};
            vector<double> spec = specul(dl, normal, dlrgb, (*sph).mat.ks, (*sph).mat.sp, view);
            for (int i=0; i<3; i++) {
                specular[i] = specular[i] + spec[i];
            }
        }
        /*
         * ------------------
         * AMBIENT      LIGHT
         * ------------------
         */ 
        vector<double> intens = al.getIntensity();
        double al_rgb[3] = {intens[0], intens[1], intens[2]};
        vector<double> ambi = ambien(al_rgb, (*sph).mat.ka);
        for (int i=0; i<3; i++) {
            ambient[i] = ambient[i] + ambi[i];
        }
    }
    double r = fmin(diffuse[0] + ambient[0] + specular[0], 1.0f);
    double g = fmin(diffuse[1] + ambient[1] + specular[1], 1.0f);
    double b = fmin(diffuse[2] + ambient[2] + specular[2], 1.0f);
    vector<double> rgb;
    rgb.push_back(r);
    rgb.push_back(g);
    rgb.push_back(b);
    return rgb;

}
void start_ray(ray_container* r, double mkr[], int past_id, int max_depth) {
    max_depth--;
    ray incoming = r->get_ray();
    shape* current;
    double t = INFINITY;
    double t_old = t;
    for (int i=0; i<shapes.size(); i++) {
        shape* testing = (shape*)(shapes[i]);
        
        double test = testing->intersect(incoming);
        if (test < t && testing->getID() != past_id) {
            t = test;
            current = shapes[i];
        }
    }
    if (t == t_old) {
        return;
    }
    vector<double> rgb = phong(incoming, t, current);
    int negligible = 0;
    for (int i=0; i<3; i++) {
        if (mkr[i] <= .00001) {
            negligible = 1;
        }
    }
    if (negligible || max_depth <= 0) {
        r->set_rgb(rgb[0], rgb[1], rgb[2]);
        return;
    }
    Vector pnt = Vector(incoming.xpos(t), incoming.ypos(t), incoming.zpos(t));
    Vector copy_v = incoming.d().copy();
    copy_v.unit();
    Vector norm = current->getNormal(pnt.xc(), pnt.yc(), pnt.zc());
    Vector copy_n = norm.copy();
    copy_n.unit();
    double dotp = copy_v.dot(copy_n);
    Vector out = Vector(copy_v.xc()-2.0*dotp*copy_n.xc(),
                        copy_v.yc()-2.0*dotp*copy_n.yc(),
                        copy_v.zc()-2.0*dotp*copy_n.zc()); 
    ray outgoing = ray(pnt, out);
    ray_container* out_cont = new ray_container(0.0, 0.0, 0.0, outgoing);
    mkr[0] = mkr[0] * current->mat.kr[0];
    mkr[1] = mkr[1] * current->mat.kr[1];
    mkr[2] = mkr[2] * current->mat.kr[2];
    start_ray(out_cont, mkr, current->getID(), max_depth);
    vector<double> retRGB = out_cont->get_rgb();
    rgb[0] = rgb[0] + fmax(0,current->mat.kr[0]*retRGB[0]);
    rgb[1] = rgb[1] + fmax(0,current->mat.kr[1]*retRGB[1]);
    rgb[2] = rgb[2] + fmax(0,current->mat.kr[2]*retRGB[2]);
    delete out_cont;
    r->set_rgb(rgb[0], rgb[1], rgb[2]);
    return;
}

void begin_raytrace(string output) {
    cout << "Total number of shapes: " << shapes.size() << "." << endl;

    double increment = 0.001;

    int num_pixels = (int)sqr(1/increment);
    unsigned width = (unsigned)(1/increment);
    unsigned height = width;

    vector<ray_container*> rays;
    rays.resize(num_pixels);
    int limit = width;
    for(int vuck = 0; vuck < limit; vuck += 1){
        for(int uck = 0; uck < limit; uck += 1){
            double kr_array[3] = {1,1,1};
            double u = uck*increment + increment/2;
            double v = vuck*increment + increment/2;

            double p[3];
            p[0] = (1-u) *(v*cam.ll[0]+ (1-v)*cam.ul[0])+ u*(v*cam.lr[0]+ (1-v)*cam.ur[0]);
            p[1] = (1-u) *(v*cam.ll[1]+ (1-v)*cam.ul[1])+ u*(v*cam.lr[1]+ (1-v)*cam.ur[1]);
            p[2] = (1-u) *(v*cam.ll[2]+ (1-v)*cam.ul[2])+ u*(v*cam.lr[2]+ (1-v)*cam.ur[2]);
            Vector origin = Vector(cam.xyz[0],cam.xyz[1],cam.xyz[2]);

            double dir[3] = {p[0] - cam.xyz[0], p[1] - cam.xyz[1], p[2] - cam.xyz[2]};
            Vector direction = Vector(dir[0],dir[1],dir[2]);
            ray r = ray(origin,direction);
            
            int i = vuck*limit + uck;
            ray_container * a = new ray_container(0.0, 0.0, 0.0, r);
            rays[i] = a;
            start_ray(rays[i],kr_array,-1, 100);
        }
    }
    //Check if file already exists
    unsigned not_found = 1;
    string temp = output;
    int counter = 1;
    while (true) {
    	vector<unsigned char> check;
    	unsigned cwidth, cheight;
    	string outfile = temp + ".png";
    	not_found = lodepng::decode(check, cwidth, cheight, outfile);
    	if (not_found) {
    		output = temp;
    		cout << "Output file can be found at " << output << ".png" << endl;
    		break;
    	}
    	else {
    		temp = output + to_string(counter);
    		counter++;
    	}
    	if (counter > 999) {
    		cout << "Error. Too many existing files in output folder. Please empty." << endl;
    		exit(0);
    	}
    }
    vector<unsigned char> outgoing;
    outgoing.resize(width*height*4);
    for (int i=0; i<num_pixels; i++) {
        vector<double> rgb = rays[i]->get_rgb();
        outgoing[4 * i + 0] = (unsigned char)(rgb[0]*255);
        outgoing[4 * i + 1] = (unsigned char)(rgb[1]*255);
        outgoing[4 * i + 2] = (unsigned char)(rgb[2]*255);
        outgoing[4 * i + 3] = 255;
    }
    unsigned has_error = lodepng::encode(output + ".png", outgoing, width, height);
    if (has_error) {
        std::cout << "Unable to write to png file. Please review desired output '" << output << "' and usage." << std::endl;
    }
    return;
}

int main(int argc, char *argv[]) {
    const int MAX_CHARS_PER_LINE = 512;
    const int MAX_TOKENS_PER_LINE = 17;
    const char* const DELIMITER = " ";

    string out;
    string readFile;


    for (int i=1; i<argc; i++) {
        string arg;
        arg = string(argv[i]);
        if (arg=="input") {
            readFile = string(argv[i+1]);
            i++;
            continue;
        }
        if (arg=="write_img") {
            out = string(argv[i+1]);
            i++;
            continue;
        }
        else {
            std::cout << "Unrecognized argument. Please review usage." << std::endl;
            exit(0);
        }
    }
    ifstream myFile;
    myFile.open(readFile);
    if(readFile == ""){
        cout << "No input provided. Please review usage." << std::endl;
        exit(0);
    }
    else{
        if (!myFile.good()) {
            std::cout << "Could not open/find input file. Please review usage." << std::endl;
            exit(0);
        }
        else{
            cout << "Parsing SCN file: " << readFile << "." << endl;
            while (!myFile.eof()){
                char buf[MAX_CHARS_PER_LINE];
                myFile.getline(buf, MAX_CHARS_PER_LINE);

                const char* token[MAX_TOKENS_PER_LINE] = {}; 
                // parse the line
                token[0] = strtok(buf, DELIMITER); // first token
                if (token[0]){ 
                    int length = 0;
                    for (int n = 1; n < MAX_TOKENS_PER_LINE; n++) {
                        token[n] = strtok(0, DELIMITER); // subsequent tokens
                        length+=1;
                        if (!token[n]){
                            break;

                        }
                    }
                    string first = string(token[0]).c_str();
                    if(first == "cam"){
                        if (has_camera) {
                            cout << "Camera/Viewport already initialized (ignoring new camera)" << endl;
                            continue; 
                        }
                        if(length >= 16){
                            if(length > 16){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double xyz[3] = {atof(string(token[1]).c_str()),atof(string(token[2]).c_str()),atof(string(token[3]).c_str())};
                            double ll[3] = {atof(string(token[4]).c_str()),atof(string(token[5]).c_str()),atof(string(token[6]).c_str())};
                            double lr[3] = {atof(string(token[7]).c_str()),atof(string(token[8]).c_str()),atof(string(token[9]).c_str())};
                            double ul[3] = {atof(string(token[10]).c_str()),atof(string(token[11]).c_str()),atof(string(token[12]).c_str())};
                            double ur[3] = {atof(string(token[13]).c_str()),atof(string(token[14]).c_str()),atof(string(token[15]).c_str())};
                            cam = camera(xyz,ll,lr,ul,ur);
                            has_camera = 1;

                        }

                        else{
                            cout << "Parsing Error for cam" << endl;

                        }
                        continue;

                    }
                    else if(first == "sph"){
                        if(length >= 5){
                            if(length > 5){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double xc = atof(string(token[1]).c_str());
                            double yc = atof(string(token[2]).c_str());
                            double zc = atof(string(token[3]).c_str());
                            double r = atof(string(token[4]).c_str());
                            if(r <= 0){
                                cout << "Parsing Error for sph" << endl;
                                continue;
                            }
                            sphere* sph = new sphere(xc,yc,zc,r);
                            sph->setMaterial(m);
                            sph->transform(trans);
                            shapes.push_back(sph);
                            continue;
                        }
                        else{
                            cerr << "Parsing Error for sph" << endl;
                            continue;
                        }

                    }
                    else if(first == "tri"){
                        if(length >= 10){
                            if(length > 10){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double a[3] = {atof(string(token[1]).c_str()),atof(string(token[2]).c_str()),atof(string(token[3]).c_str())};
                            double b[3] = {atof(string(token[4]).c_str()),atof(string(token[5]).c_str()),atof(string(token[6]).c_str())};
                            double c[3] = {atof(string(token[7]).c_str()),atof(string(token[8]).c_str()),atof(string(token[9]).c_str())};
                            triangle* tri = new triangle(a,b,c);
                            tri->setMaterial(m);
                            tri->transform(trans);
                            shapes.push_back(tri);
                            continue;
                        }
                        else{
                            cout << "Parsing Error for tri" << endl;
                        }

                    }
                    else if(first == "obj"){
                        if(length >= 2){
                            if(length > 2){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            objFile =  string(token[1]).c_str();

                            ifstream myFile;
                            myFile.open(objFile);
                            if (!myFile.good()) {
                                std::cout << "Could not open/find obj file. Please review usage." << std::endl;
                            }
                            else{
                                cout << "Parsing OBJ file:" << objFile << "." << endl;
                                while (!myFile.eof()){
                                    char buf[MAX_CHARS_PER_LINE];
                                    myFile.getline(buf, MAX_CHARS_PER_LINE);

                                    const char* token[MAX_TOKENS_PER_LINE] = {}; 
                                    // parse the line
                                    token[0] = strtok(buf, DELIMITER); // first token
                                    int length = 0;
                                    if (token[0]){ 
                                        for (int n = 1; n < MAX_TOKENS_PER_LINE; n++) {
                                            token[n] = strtok(0, DELIMITER); // subsequent tokens
                                            length +=1;
                                            if (!token[n]){
                                                break;
                                            }
                                        }
                                        string first = string(token[0]).c_str();
                                        if(first == "v"){
                                            if(length >= 4){ 
                                                if(length > 4){
                                                    cout << "Parsing Error in v (ignoring extra arguments)" << endl;
                                                }
                                                vector<double> vert;
                                                double xc = atof(string(token[1]).c_str());
                                                vert.push_back(xc);
                                                double yc = atof(string(token[2]).c_str());
                                                vert.push_back(yc);
                                                double zc = atof(string(token[3]).c_str());
                                                vert.push_back(zc);
                                                vertices.push_back(vert);
                                            }
                                            else{
                                                cout << "Parsing Error for v" << endl;
                                            }
                                        }
                                        if(first == "vn"){
                                            if(length >= 4){
                                                if(length > 4){
                                                    cout << "Parsing Error in vn (ignoring extra arguments)" << endl;
                                                }
                                                vector<double> norm;
                                                double xc = atof(string(token[1]).c_str());
                                                norm.push_back(xc);
                                                double yc = atof(string(token[2]).c_str());
                                                norm.push_back(yc);
                                                double zc = atof(string(token[3]).c_str());
                                                norm.push_back(zc);
                                                normals.push_back(norm);
                                            }
                                            else{
                                                cout << "Parsing Error for vn" << endl;
                                            }
                                        }
                                        if(first == "f"){
                                            if(length >= 4){
                                                if(length > 4){
                                                    cout << "Parsing Error in f (ignoring extra arguments)" << endl;
                                                }
                                                try{
                                                    string first = string(token[1]).c_str();
                                                    string second = string(token[2]).c_str();
                                                    string third = string(token[3]).c_str();
                                                    bool found1 = false;
                                                    bool found2 = false;
                                                    bool found3 = false;
                                                    if (first.find("//") != std::string::npos) {
                                                        found1 = true;
                                                    }
                                                    if (second.find("//") != std::string::npos) {
                                                        found2 = true;
                                                    }
                                                    if (third.find("//") != std::string::npos) {
                                                        found3 = true;
                                                    }
                                                    if(!found1 || !found2 || !found3){
                                                        throw std::invalid_argument("Continuing on to other f attempt");
                                                    }
                                                    double* v1;
                                                    double* n1;
                                                    double* v2;
                                                    double* n2;
                                                    double* v3;
                                                    double* n3;

                                                    std::string delimiter = "//";
                                                    size_t pos = 0;
                                                    std::string token;
                                                    while ((pos = first.find(delimiter)) != std::string::npos) {
                                                        token = first.substr(0, pos);
                                                        v1 = &(vertices.at(atof(string(token).c_str()) -1)[0]);
                                                        first.erase(0, pos + delimiter.length());
                                                    }
                                                    n1 = &(normals.at(atof(string(first).c_str())-1)[0]);
                                                    while ((pos = second.find(delimiter)) != std::string::npos) {
                                                        token = second.substr(0, pos);
                                                        v2 = &(vertices.at(atof(string(token).c_str())-1)[0]);
                                                        second.erase(0, pos + delimiter.length());
                                                    }
                                                    n2 = &(normals.at(atof(string(second).c_str())-1)[0]);
                                                    while ((pos = third.find(delimiter)) != std::string::npos) {
                                                        token = third.substr(0, pos);
                                                        v3 = &(vertices.at(atof(string(token).c_str())-1)[0]);
                                                        third.erase(0, pos + delimiter.length());
                                                    }
                                                    n3 = &(normals.at(atof(string(third).c_str())-1)[0]);
                                                    triangle* tri = new triangle(v1,n1,v2,n2,v3,n3);
                                                    tri->setMaterial(m);
                                                    tri->transform(trans);
                                                    shapes.push_back(tri);

                                                }
                                                catch(...){
                                                    try{
                                                        double* v1;
                                                        double* v2;
                                                        double* v3;
                                                        v1 = &(vertices.at(atof(string(token[1]).c_str()) -1)[0]);
                                                        v2 = &(vertices.at(atof(string(token[2]).c_str()) -1)[0]);
                                                        v3 = &(vertices.at(atof(string(token[3]).c_str()) -1)[0]);

                                                        triangle* tri = new triangle(v1,v2,v3);
                                                        tri->setMaterial(m);
                                                        tri->transform(trans);

                                                        shapes.push_back(tri);
                                                        continue;

                                                    }
                                                    catch(...){
                                                        cout << "Parsing Error for f" << endl;
                                                    }


                                                }
                                            }
                                            else{
                                                cout << "Parsing Error for f: too few arguments." << endl;
                                            }
                                        }
                                    }
                                }
                            }
                            myFile.close();
                            continue;
                        }
                        else{
                            cout << "Parsing Error for obj" << endl;
                            continue;
                        }

                    }
                    else if(first == "ltp"){
                        if(length >= 7){
                            if(length > 8){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double x = atof(string(token[1]).c_str());
                            double y = atof(string(token[2]).c_str());
                            double z = atof(string(token[3]).c_str());
                            double r = atof(string(token[4]).c_str());
                            double g = atof(string(token[5]).c_str());
                            double b = atof(string(token[6]).c_str());
                            double f = 0;
                            if(length > 7){
                                f = atof(string(token[7]).c_str());
                            }
                            if(r > 1 || g > 1 || b > 1 || !(f == 0 || f == 1 || f ==2) || r < 0 || g < 0 || b < 0){
                                cout << "Parsing Error for ltp" << endl;
                                continue;
                            }
                            plight pl = plight(x,y,z,r,g,b,f);
                            plights.push_back(pl);
                            continue;
                        }
                        else{
                            cout << "Parsing Error for ltp" << endl;
                        }

                    }
                    else if(first == "ltd"){
                        if(length >= 7){
                            if(length > 7){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double x = atof(string(token[1]).c_str());
                            double y = atof(string(token[2]).c_str());
                            double z = atof(string(token[3]).c_str());
                            double r = atof(string(token[4]).c_str());
                            double g = atof(string(token[5]).c_str());
                            double b = atof(string(token[6]).c_str());
                            if(r > 1 || g > 1 || b > 1 || r < 0 || g < 0 || b < 0){
                                cout << "Parsing Error for ltd" << endl;
                                continue;
                            }
                            dlight dl = dlight(x,y,z,r,g,b);
                            dlights.push_back(dl);
                            continue;
                        }
                        else{
                            cout << "Parsing Error for ltd" << endl;
                        }

                    }
                    else if(first == "lta"){
                        if(length >= 4){
                            if(length > 4){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double r = atof(string(token[1]).c_str());
                            double g = atof(string(token[2]).c_str());
                            double b = atof(string(token[3]).c_str());
                            if(r > 1 || g > 1 || b > 1 || r < 0 || g < 0 || b < 0){
                                cout << "Parsing Error for ltp" << endl;
                                continue;
                            }
                            alight ali = alight(r,g,b);
                            al = ali;
                            continue;
                        }
                        else{
                            cout << "Parsing Error for lta" << endl;
                        }

                    }
                    else if(first == "mat"){
                        if(length >= 14){
                            if(length > 14){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }

                            double ka[3] = {atof(string(token[1]).c_str()),atof(string(token[2]).c_str()),atof(string(token[3]).c_str())};
                            double kd[3] = {atof(string(token[4]).c_str()),atof(string(token[5]).c_str()),atof(string(token[6]).c_str())};
                            double ks[3] = {atof(string(token[7]).c_str()),atof(string(token[8]).c_str()),atof(string(token[9]).c_str())};
                            double sp = atof(string(token[10]).c_str());
                            double kr[3] = {atof(string(token[11]).c_str()),atof(string(token[12]).c_str()),atof(string(token[13]).c_str())};
                            bool pass = true;
                            for (int j=0; j<3; j++) {
                                if (kd[j] > 1 || kd[j] < 0 || ka[j] > 1 || ka[j] < 0 || ks[j] > 1 || ks[j] < 0 || kr[j] > 1 || kr[j] < 0) {
                                    pass = false;
                                    std::cout << "K values not between 0 and 1" << std::endl;
                                }
                            }
                            if(pass){
                                m = new material(ka,kd,ks,sp,kr);
                            }
                            continue;
                        }
                        else{
                            cout << "Parsing Error for mat" << endl;
                        }

                    }
                    else if(first == "xft"){
                        if(length >= 4){
                            if(length > 4){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double x = atof(string(token[1]).c_str());
                            double y = atof(string(token[2]).c_str());
                            double z = atof(string(token[3]).c_str());
                            matrix mul = matrix(x,y,z,0);
                            trans.multiplym(mul);
                            continue;
                        }
                        else{
                            cout << "Parsing Error for xft" << endl;
                        }
                    }
                    else if(first == "xfr"){
                        if(length >= 4){
                            if(length > 4){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double x = atof(string(token[1]).c_str());
                            double y = atof(string(token[2]).c_str());
                            double z = atof(string(token[3]).c_str());
                            matrix mul = matrix(x,y,z,2);
                            trans.multiplym(mul);
                            continue;
                        }
                        else{
                            cout << "Parsing Error for xfr" << endl;
                        }
                    }
                    else if(first == "xfs"){
                        if(length >= 4){
                            if(length > 4){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }
                            double x = atof(string(token[1]).c_str());
                            double y = atof(string(token[2]).c_str());
                            double z = atof(string(token[3]).c_str());
                            if (x == 0 || y == 0 || z == 0) {
                                cout << "Scaling transformation must be non-zero." << endl;
                                continue;
                            }
                            matrix mul = matrix(x,y,z,1);
                            trans.multiplym(mul);
                            continue;
                        }
                        else{
                            cout << "Parsing Error for xfs" << endl;
                        }


                    }
                    else if(first == "xfz"){
                        if(length >= 1){
                            if(length > 1){
                                cout << "Parsing Error (ignoring extra arguments)" << endl;
                            }

                            trans = matrix();
                            continue;
                        }
                        else{
                            cout << "Parsing Error for xfz" << endl;
                        }
                    }
                }
            }
        }
    }
    if (!has_camera) {
        cout << "Did not provide camera/viewport. Please review usage." << endl;
        exit(0);
    }

    myFile.close();

    begin_raytrace(out);
    return 0;
}
