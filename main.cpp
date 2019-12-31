//Dingheng Zhang
//Project 3 test


#define GL_SILENCE_DEPRECATION
#ifdef WIN32
#include <windows.h>
#endif

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else //linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif


//other includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

#define PI acos(-1)


class VerticePoint{
public:
    GLfloat x,y,z;
    GLfloat r,g,b;
};

class FaceSpecs{
public:
    GLint p1,p2,p3;
    GLfloat s;
};

class wcPt3D {
public:
GLfloat x, y, z;
};

class LightIntense{
public:
    GLfloat r,g,b;
};

class ShowPoint{
public:
    GLfloat o1,o2;
    GLfloat r,g,b;
};

wcPt3D Pmin, Pmax;

string File_Address;

int num_poly;

int *num_vertices;
int *num_faces;
int **num_edgepoints;

VerticePoint **Vertices;
FaceSpecs **Faces;
VerticePoint ***EdgePoint;
ShowPoint ***AllPointXY;
ShowPoint ***AllPointYZ;
ShowPoint ***AllPointXZ;
int **num_allpoints_xy;
int **num_allpoints_yz;
int **num_allpoints_xz;

//initial value
int num_point_init = 3000;
int num_face_init = 5000;

//view window
int width=500;
int length=500;
float pixel_size;
int win_width=200;
int win_height=200;

//phong value
GLfloat IA,K_distance,IL;
LightIntense ka,ks;
wcPt3D front_P,light_P;

//phong data
LightIntense **phone_vertices;
int **faces_per_vertice;
GLfloat **specu_per_vertice;

//normal vector
wcPt3D **normal_vector;
wcPt3D **normal_vector_point;

//for all points
wcPt3D **face_min;
wcPt3D **face_max;
int **x_buffer_face;
int **y_buffer_face;
int **z_buffer_face;

//halftong
ShowPoint ***HalftongXY;
ShowPoint ***HalftongYZ;
ShowPoint ***HalftongXZ;

//operates
int phong_sign=0;
int halftone_sign=0;


void FMinMaxP();

int InitializeX(GLfloat num1);
int InitializeY(GLfloat num1);
int InitializeZ(GLfloat num1);

void ReadFile()
{
    cout<<"----------------------------------------------------------------"<<endl;
    cout<<"Please type in the file address. Input like 'cube_and_icosahedron.txt' "<<endl;
    cin>>File_Address;
    
    std::ifstream infile;
    infile.open(File_Address);
    
    infile>>num_poly;
    
    Vertices = new VerticePoint *[num_poly];
    Faces = new FaceSpecs *[num_poly];
    EdgePoint = new VerticePoint **[num_poly];
    num_vertices = new int[num_poly];
    num_faces = new int[num_poly];
    num_edgepoints = new int *[num_poly];
    normal_vector = new wcPt3D *[num_poly];
    normal_vector_point = new wcPt3D *[num_poly];
    x_buffer_face = new int *[num_poly];
    y_buffer_face = new int *[num_poly];
    z_buffer_face = new int *[num_poly];
    
    for (int i=0;i<num_poly;i++)
    {
        Vertices[i] = new VerticePoint[num_point_init];
        Faces[i] = new FaceSpecs[num_face_init];
        EdgePoint[i] = new VerticePoint *[num_face_init];
        for(int j=0;j<num_face_init;j++)
        {
            EdgePoint[i][j] = new VerticePoint[width*num_point_init];
        }
        num_edgepoints[i] = new int[num_face_init];
        normal_vector[i] = new wcPt3D[num_face_init];
        normal_vector_point[i] = new wcPt3D[num_point_init];
        x_buffer_face[i] = new int[num_face_init];
        y_buffer_face[i] = new int[num_face_init];
        z_buffer_face[i] = new int[num_face_init];
    }
    
    for(int i=0;i<num_poly;i++)
    {
        infile>>num_vertices[i];
        for(int j=0;j<num_vertices[i];j++)
        {
            infile>>Vertices[i][j].x;
            infile>>Vertices[i][j].y;
            infile>>Vertices[i][j].z;
        }
        for(int j=0;j<num_vertices[i];j++)
        {
            infile>>Vertices[i][j].r;
            Vertices[i][j].r = Vertices[i][j].r/255;
            infile>>Vertices[i][j].g;
            Vertices[i][j].g = Vertices[i][j].g/255;
            infile>>Vertices[i][j].b;
            Vertices[i][j].b = Vertices[i][j].b/255;
            //cout<<Vertices[i][j].r<<'\t'<<Vertices[i][j].g<<'\t'<<Vertices[i][j].b<<endl;
        }
        infile>>num_faces[i];
        for(int j=0;j<num_faces[i];j++)
        {
            infile>>Faces[i][j].p1;
            infile>>Faces[i][j].p2;
            infile>>Faces[i][j].p3;
        }
        for(int j=0;j<num_faces[i];j++)
            infile>>Faces[i][j].s;
    }
    FMinMaxP();
    for(int i =0;i<num_poly;i++)
        for(int j=0;j<num_vertices[i];j++)
        {
            Vertices[i][j].x = InitializeX(Vertices[i][j].x);
            Vertices[i][j].y = InitializeY(Vertices[i][j].y);
            Vertices[i][j].z = InitializeZ(Vertices[i][j].z);
        }
    cout<<"Finish Reading the File!"<<endl;
}

void FMinMaxP()
{
    Pmin.x = Pmin.y = Pmin.z = 999;
    Pmax.x = Pmax.y = Pmax.z = -999;
    
    for(int i =0;i<num_poly;i++)
        for(int j =0;j<num_vertices[i];j++)
        {

            if(Vertices[i][j].x<Pmin.x) Pmin.x = Vertices[i][j].x;
            if(Vertices[i][j].y<Pmin.y) Pmin.y = Vertices[i][j].y;
            if(Vertices[i][j].z<Pmin.z) Pmin.z = Vertices[i][j].z;
            if(Vertices[i][j].x>Pmax.x) Pmax.x = Vertices[i][j].x;
            if(Vertices[i][j].y>Pmax.y) Pmax.y = Vertices[i][j].y;
            if(Vertices[i][j].z>Pmax.z) Pmax.z = Vertices[i][j].z;
        }
}

void MinMaxFaceBuffer()
{
    int **x_buffer_sign;
    int **y_buffer_sign;
    int **z_buffer_sign;
    face_min = new wcPt3D *[num_poly];
    face_max = new wcPt3D *[num_poly];
    x_buffer_sign = new int*[num_poly];
    y_buffer_sign = new int*[num_poly];
    z_buffer_sign = new int*[num_poly];
    for(int i=0;i<num_poly;i++)
    {
        face_min[i] = new wcPt3D[num_face_init];
        face_max[i] = new wcPt3D[num_face_init];
        x_buffer_sign[i] = new int[num_face_init];
        y_buffer_sign[i] = new int[num_face_init];
        z_buffer_sign[i] = new int[num_face_init];
        for(int j=0;j<num_faces[i];j++)
        {
            x_buffer_face[i][j] = 0;
            y_buffer_face[i][j] = 0;
            z_buffer_face[i][j] = 0;
            x_buffer_sign[i][j] = 0;
            y_buffer_sign[i][j] = 0;
            z_buffer_sign[i][j] = 0;
            face_min[i][j].x = face_min[i][j].y = face_min[i][j].z = 999;
            face_max[i][j].x = face_max[i][j].y = face_max[i][j].z = -999;
        }
    }
    for(int i=0;i<num_poly;i++)
            for(int j=0;j<num_faces[i];j++)
            {
                if(face_min[i][j].x>Vertices[i][Faces[i][j].p1-1].x) face_min[i][j].x=Vertices[i][Faces[i][j].p1-1].x;
                if(face_min[i][j].y>Vertices[i][Faces[i][j].p1-1].y) face_min[i][j].y=Vertices[i][Faces[i][j].p1-1].y;
                if(face_min[i][j].z>Vertices[i][Faces[i][j].p1-1].z) face_min[i][j].z=Vertices[i][Faces[i][j].p1-1].z;
                if(face_max[i][j].x<Vertices[i][Faces[i][j].p1-1].x) face_max[i][j].x=Vertices[i][Faces[i][j].p1-1].x;
                if(face_max[i][j].y<Vertices[i][Faces[i][j].p1-1].y) face_max[i][j].y=Vertices[i][Faces[i][j].p1-1].y;
                if(face_max[i][j].z<Vertices[i][Faces[i][j].p1-1].z) face_max[i][j].z=Vertices[i][Faces[i][j].p1-1].z;
                
                if(face_min[i][j].x>Vertices[i][Faces[i][j].p2-1].x) face_min[i][j].x=Vertices[i][Faces[i][j].p2-1].x;
                if(face_min[i][j].y>Vertices[i][Faces[i][j].p2-1].y) face_min[i][j].y=Vertices[i][Faces[i][j].p2-1].y;
                if(face_min[i][j].z>Vertices[i][Faces[i][j].p2-1].z) face_min[i][j].z=Vertices[i][Faces[i][j].p2-1].z;
                if(face_max[i][j].x<Vertices[i][Faces[i][j].p2-1].x) face_max[i][j].x=Vertices[i][Faces[i][j].p2-1].x;
                if(face_max[i][j].y<Vertices[i][Faces[i][j].p2-1].y) face_max[i][j].y=Vertices[i][Faces[i][j].p2-1].y;
                if(face_max[i][j].z<Vertices[i][Faces[i][j].p2-1].z) face_max[i][j].z=Vertices[i][Faces[i][j].p2-1].z;
                
                if(face_min[i][j].x>Vertices[i][Faces[i][j].p3-1].x) face_min[i][j].x=Vertices[i][Faces[i][j].p3-1].x;
                if(face_min[i][j].y>Vertices[i][Faces[i][j].p3-1].y) face_min[i][j].y=Vertices[i][Faces[i][j].p3-1].y;
                if(face_min[i][j].z>Vertices[i][Faces[i][j].p3-1].z) face_min[i][j].z=Vertices[i][Faces[i][j].p3-1].z;
                if(face_max[i][j].x<Vertices[i][Faces[i][j].p3-1].x) face_max[i][j].x=Vertices[i][Faces[i][j].p3-1].x;
                if(face_max[i][j].y<Vertices[i][Faces[i][j].p3-1].y) face_max[i][j].y=Vertices[i][Faces[i][j].p3-1].y;
                if(face_max[i][j].z<Vertices[i][Faces[i][j].p3-1].z) face_max[i][j].z=Vertices[i][Faces[i][j].p3-1].z;
            
    //            cout<<"numnumnumnumnum"<<endl;
//                cout<<face_min[i][j].x<<'\t'<<face_min[i][j].y<<'\t'<<face_min[i][j].z<<endl;
//                cout<<face_max[i][j].x<<'\t'<<face_max[i][j].y<<'\t'<<face_max[i][j].z<<endl;
            }
    for(int i=0;i<num_poly;i++)
    {
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_faces[i];k++)
            {
                if(face_max[i][k].x>face_max[i][j].x) x_buffer_sign[i][k]++;
                if(face_max[i][k].y>face_max[i][j].y) y_buffer_sign[i][k]++;
                if(face_max[i][k].z>face_max[i][j].z) z_buffer_sign[i][k]++;
            }
        for(int j=0;j<num_faces[i];)
        {
            int sign = j;
            for(int k=0;k<num_faces[i];k++)
            {
                if(z_buffer_sign[i][k]==sign)
                {
                    z_buffer_face[i][j] = k;
//                    cout<<"j"<<'\t'<<j<<'\t'<<"k"<<'\t'<<k<<'\t'<<endl;
                    j++;
                }
            }
        }
        for(int j=0;j<num_faces[i];)
        {
                    int sign = j;
                    for(int k=0;k<num_faces[i];k++)
                    {
                        if(x_buffer_sign[i][k]==sign)
                        {
                            x_buffer_face[i][j] = k;
        //                    cout<<"j"<<'\t'<<j<<'\t'<<"k"<<'\t'<<k<<'\t'<<endl;
                            j++;
                        }
                    }
        }
        for(int j=0;j<num_faces[i];)
                {
                    int sign = j;
                    for(int k=0;k<num_faces[i];k++)
                    {
                        if(y_buffer_sign[i][k]==sign)
                        {
                            y_buffer_face[i][j] = k;
        //                    cout<<"j"<<'\t'<<j<<'\t'<<"k"<<'\t'<<k<<'\t'<<endl;
                            j++;
                        }
                    }
                }
        
    }
    for(int i=0;i<num_poly;i++)
    {
        delete [] x_buffer_sign[i];
        delete [] y_buffer_sign[i];
        delete [] z_buffer_sign[i];
    }
    delete [] x_buffer_sign;
    delete [] y_buffer_sign;
    delete [] z_buffer_sign;
}

int InitializeX(GLfloat num1)
{
    GLfloat px;
    px = Pmax.x-Pmin.x+0.2;
    num1 = (num1 - Pmin.x+0.1)/px*width;
    int num2;
    num2 = int(num1);
    return num2;
}
int InitializeY(GLfloat num1)
{
    GLfloat py;
    py = Pmax.y-Pmin.y+0.2;
    num1 = (num1 - Pmin.y+0.1)/py*width;
    int num2;
    num2 = int(num1);
    return num2;
}
int InitializeZ(GLfloat num1)
{
    GLfloat pz;
    pz = Pmax.z-Pmin.z+0.2;
    num1 = (num1 - Pmin.z+0.1)/pz*width;
    int num2;
    num2 = int(num1);
    return num2;
}


//float round_value(float v)
//{
//    v=v*100;
//    v = floor(v + 0.5);
//    v = v/100;
//    return v;
//}

void LineDDA(int p1,int p2,int i,int j)
{
    GLfloat dx,dy,dz,dmax;
    dx = Vertices[i][p2].x - Vertices[i][p1].x;
    dy = Vertices[i][p2].y - Vertices[i][p1].y;
    dz = Vertices[i][p2].z - Vertices[i][p1].z;
    if(fabs(dx)>fabs(dy))
        dmax = fabs(dx);
    else
        dmax = fabs(dy);
    if(dmax<fabs(dz))
        dmax = fabs(dz);
    GLfloat AInc,BInc,CInc;
    GLfloat a=Vertices[i][p1].x;
    GLfloat b=Vertices[i][p1].y;
    GLfloat c=Vertices[i][p1].z;
    AInc = dx/dmax;
    BInc = dy/dmax;
    CInc = dz/dmax;
//    EdgePoint[i][j][num_edgepoints[i][j]].x = round_value(Vertices[i][p1].x);
//    EdgePoint[i][j][num_edgepoints[i][j]].y = round_value(Vertices[i][p1].y);
//    EdgePoint[i][j][num_edgepoints[i][j]].z = round_value(Vertices[i][p1].z);
//    EdgePoint[i][j][num_edgepoints[i][j]].r = Vertices[i][p1].r;
//    EdgePoint[i][j][num_edgepoints[i][j]].g = Vertices[i][p1].g;
//    EdgePoint[i][j][num_edgepoints[i][j]].b = Vertices[i][p1].b;
//    num_edgepoints[i][j]++;
//    for (int k=1;k<dmax;k++)
//    {
//        GLfloat proport1,proport2;
//        a+=AInc;
//        b+=BInc;
//        c+=CInc;
//        EdgePoint[i][j][num_edgepoints[i][j]].x = round_value(a);
//        EdgePoint[i][j][num_edgepoints[i][j]].y = round_value(b);
//        EdgePoint[i][j][num_edgepoints[i][j]].z = round_value(c);
//        proport1 = (a-Vertices[i][p2].x)/(Vertices[i][p1].x-Vertices[i][p2].x);
//        proport2 = (Vertices[i][p1].x-a)/(Vertices[i][p1].x-Vertices[i][p2].x);
//        EdgePoint[i][j][num_edgepoints[i][j]].r = proport1*Vertices[i][p1].r + proport2*Vertices[i][p2].r;
//        EdgePoint[i][j][num_edgepoints[i][j]].g = proport1*Vertices[i][p1].g + proport2*Vertices[i][p2].g;
//        EdgePoint[i][j][num_edgepoints[i][j]].b = proport1*Vertices[i][p1].b + proport2*Vertices[i][p2].b;
//        num_edgepoints[i][j]++;
//    }
//    EdgePoint[i][j][num_edgepoints[i][j]].x = round_value(Vertices[i][p2].x);
//    EdgePoint[i][j][num_edgepoints[i][j]].y = round_value(Vertices[i][p2].y);
//    EdgePoint[i][j][num_edgepoints[i][j]].z = round_value(Vertices[i][p2].z);
//    EdgePoint[i][j][num_edgepoints[i][j]].r = Vertices[i][p2].r;
//    EdgePoint[i][j][num_edgepoints[i][j]].g = Vertices[i][p2].g;
//    EdgePoint[i][j][num_edgepoints[i][j]].b = Vertices[i][p2].b;
//    num_edgepoints[i][j]++;
    
    EdgePoint[i][j][num_edgepoints[i][j]].x = Vertices[i][p1].x;
    EdgePoint[i][j][num_edgepoints[i][j]].y = Vertices[i][p1].y;
    EdgePoint[i][j][num_edgepoints[i][j]].z = Vertices[i][p1].z;
    EdgePoint[i][j][num_edgepoints[i][j]].r = Vertices[i][p1].r;
    EdgePoint[i][j][num_edgepoints[i][j]].g = Vertices[i][p1].g;
    EdgePoint[i][j][num_edgepoints[i][j]].b = Vertices[i][p1].b;
//    cout<<num_edgepoints[i][j]<<endl;
//    cout<<EdgePoint[i][j][num_edgepoints[i][j]].x<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].y<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].z<<endl;
    num_edgepoints[i][j]++;
    for (int k=1;k<dmax;k++)
    {
        GLfloat proport1,proport2;
        a+=AInc;
        b+=BInc;
        c+=CInc;
//        a = round_value(a);
//        b = round_value(b);
//        c = round_value(c);
        EdgePoint[i][j][num_edgepoints[i][j]].x = a;
        EdgePoint[i][j][num_edgepoints[i][j]].y = b;
        EdgePoint[i][j][num_edgepoints[i][j]].z = c;
        if((Vertices[i][p1].x-Vertices[i][p2].x)!=0)
        {
            proport1 = (a-Vertices[i][p2].x)/(Vertices[i][p1].x-Vertices[i][p2].x);
            proport2 = (Vertices[i][p1].x-a)/(Vertices[i][p1].x-Vertices[i][p2].x);
        }
        else if ((Vertices[i][p1].y-Vertices[i][p2].y)!=0)
        {
            proport1 = (b-Vertices[i][p2].y)/(Vertices[i][p1].y-Vertices[i][p2].y);
            proport2 = (Vertices[i][p1].y-b)/(Vertices[i][p1].y-Vertices[i][p2].y);
        }
        else
        {
            proport1 = (c-Vertices[i][p2].z)/(Vertices[i][p1].z-Vertices[i][p2].z);
            proport2 = (Vertices[i][p1].z-c)/(Vertices[i][p1].z-Vertices[i][p2].z);
        }
//        cout<<"test"<<'\t'<<proport1<<'\t'<<proport2<<endl;
        EdgePoint[i][j][num_edgepoints[i][j]].r = proport1*Vertices[i][p1].r + proport2*Vertices[i][p2].r;
        EdgePoint[i][j][num_edgepoints[i][j]].g = proport1*Vertices[i][p1].g + proport2*Vertices[i][p2].g;
        EdgePoint[i][j][num_edgepoints[i][j]].b = proport1*Vertices[i][p1].b + proport2*Vertices[i][p2].b;
//        cout<<num_edgepoints[i][j]<<endl;
//        if(j==3)
//        {
//            cout<<proport1<<'\t'<<proport2<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].r<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].g<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].b<<endl;
//        }
        num_edgepoints[i][j]++;
    }
    EdgePoint[i][j][num_edgepoints[i][j]].x = Vertices[i][p2].x;
    EdgePoint[i][j][num_edgepoints[i][j]].y = Vertices[i][p2].y;
    EdgePoint[i][j][num_edgepoints[i][j]].z = Vertices[i][p2].z;
    EdgePoint[i][j][num_edgepoints[i][j]].r = Vertices[i][p2].r;
    EdgePoint[i][j][num_edgepoints[i][j]].g = Vertices[i][p2].g;
    EdgePoint[i][j][num_edgepoints[i][j]].b = Vertices[i][p2].b;
//    cout<<num_edgepoints[i][j]<<endl;
//    cout<<EdgePoint[i][j][num_edgepoints[i][j]].x<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].y<<'\t'<<EdgePoint[i][j][num_edgepoints[i][j]].z<<endl;
    num_edgepoints[i][j]++;
}

void Get_SidePoints()
{
    for(int i=0;i<num_poly;i++)
    {
        for(int j=0;j<num_faces[i];j++)
        {
            num_edgepoints[i][j] = 0;
//            if((Vertices[i][Faces[i][j].p1-1].r!=0)&&(Vertices[i][Faces[i][j].p1-1].g!=0)&&(Vertices[i][Faces[i][j].p1-1].b!=0))
//            {
                LineDDA(Faces[i][j].p1-1,Faces[i][j].p2-1,i,j);
                LineDDA(Faces[i][j].p2-1,Faces[i][j].p3-1,i,j);
                LineDDA(Faces[i][j].p3-1,Faces[i][j].p1-1,i,j);
//            }
          
        }
    }
}

void Get_AllPoints()
{
//    MinMaxFaceBuffer();
//    wcPt3D **face_min;
//    wcPt3D **face_max;
    AllPointXY = new ShowPoint **[num_poly];
    AllPointYZ = new ShowPoint **[num_poly];
    AllPointXZ = new ShowPoint **[num_poly];
//    face_min = new wcPt3D *[num_poly];
//    face_max = new wcPt3D *[num_poly];
    num_allpoints_xy = new int *[num_poly];
    num_allpoints_yz = new int *[num_poly];
    num_allpoints_xz = new int *[num_poly];
    for(int i=0;i<num_poly;i++)
    {
        AllPointXY[i] = new ShowPoint *[num_face_init];
        AllPointYZ[i] = new ShowPoint *[num_face_init];
        AllPointXZ[i] = new ShowPoint *[num_face_init];
//        face_min[i] = new wcPt3D[num_face_init];
//        face_max[i] = new wcPt3D[num_face_init];
        num_allpoints_xy[i] = new int[num_face_init];
        num_allpoints_yz[i] = new int[num_face_init];
        num_allpoints_xz[i] = new int[num_face_init];
        for(int j=0;j<num_faces[i];j++)
        {
            AllPointXY[i][j] = new ShowPoint[num_point_init*num_point_init];
            AllPointYZ[i][j] = new ShowPoint[num_point_init*num_point_init];
            AllPointXZ[i][j] = new ShowPoint[num_point_init*num_point_init];
//            face_min[i][j].x = face_min[i][j].y = face_min[i][j].z = 999;
//            face_max[i][j].x = face_max[i][j].y = face_max[i][j].z = -999;
        }
    }
    
//    for(int i=0;i<num_poly;i++)
//        for(int j=0;j<num_faces[i];j++)
//        {
//            for(int k=0;k<num_edgepoints[i][j];k++)
//            {
//                if(face_min[i][j].x>EdgePoint[i][j][k].x) face_min[i][j].x=EdgePoint[i][j][k].x;
//                if(face_min[i][j].y>EdgePoint[i][j][k].y) face_min[i][j].y=EdgePoint[i][j][k].y;
//                if(face_min[i][j].z>EdgePoint[i][j][k].z) face_min[i][j].z=EdgePoint[i][j][k].z;
//                if(face_max[i][j].x<EdgePoint[i][j][k].x) face_max[i][j].x=EdgePoint[i][j][k].x;
//                if(face_max[i][j].y<EdgePoint[i][j][k].y) face_max[i][j].y=EdgePoint[i][j][k].y;
//                if(face_max[i][j].z<EdgePoint[i][j][k].z) face_max[i][j].z=EdgePoint[i][j][k].z;
//            }
////            cout<<"numnumnumnumnum"<<endl;
////            cout<<face_min[i][j].x<<'\t'<<face_min[i][j].y<<'\t'<<face_min[i][j].z<<endl;
////            cout<<face_max[i][j].x<<'\t'<<face_max[i][j].y<<'\t'<<face_max[i][j].z<<endl;
//        }
    
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
        {
            num_allpoints_xy[i][j] = 0;
            num_allpoints_yz[i][j] = 0;
            num_allpoints_xz[i][j] = 0;
            GLfloat dx,dy,dz;
            dx = face_max[i][j].x - face_min[i][j].x;
            dy = face_max[i][j].y - face_min[i][j].y;
            dz = face_max[i][j].z - face_min[i][j].z;
            
//            cout<<"dy =   "<<dy<<endl;
            if(dy>1)
            {
                GLfloat min_sign=999;
                GLfloat max_sign=-999;
                int min_sign_num = 0;
                int max_sign_num = 0;
                GLfloat proport_1,proport_2;
                for(float k=1;k<dy;k=k+1)
                {
                    min_sign=999;
                    max_sign=-999;
                    min_sign_num = 0;
                    max_sign_num = 0;
                    for(int l=0;l<num_edgepoints[i][j];l++)
                    {
                        if((EdgePoint[i][j][l].y>=(face_min[i][j].y+k-1))&&(EdgePoint[i][j][l].y<(face_min[i][j].y+k)))
                        {
                            if(min_sign>EdgePoint[i][j][l].x)
                            {
                                min_sign = EdgePoint[i][j][l].x;
                                min_sign_num = l;
                                
                            }
                            if(max_sign<EdgePoint[i][j][l].x)
                            {
                                max_sign = EdgePoint[i][j][l].x;
                                max_sign_num = l;
                            }
                        }
                        
                    }
                    GLfloat dx_sign,step;
                    dx_sign = max_sign - min_sign;
                    step = dx_sign/dy/5;
                    if((min_sign!=999)&&(dx_sign>0))
                    {
                        
                        for(float l=step;l<(max_sign-min_sign);l=l+step)
                        {
                            AllPointXY[i][j][num_allpoints_xy[i][j]].o1 = min_sign+l;
                            AllPointXY[i][j][num_allpoints_xy[i][j]].o2 =face_min[i][j].y+k-0.5;
                            proport_1 = (max_sign-(min_sign+l))/(max_sign-min_sign);
                            proport_2 = ((min_sign+l)-min_sign)/(max_sign-min_sign);
                            AllPointXY[i][j][num_allpoints_xy[i][j]].r = proport_1 * EdgePoint[i][j][min_sign_num].r + proport_2 * EdgePoint[i][j][max_sign_num].r;
                            AllPointXY[i][j][num_allpoints_xy[i][j]].g = proport_1 * EdgePoint[i][j][min_sign_num].g + proport_2 * EdgePoint[i][j][max_sign_num].g;
                            AllPointXY[i][j][num_allpoints_xy[i][j]].b = proport_1 * EdgePoint[i][j][min_sign_num].b + proport_2 * EdgePoint[i][j][max_sign_num].b;
//                            if((AllPointXY[i][j][num_allpoints_xy[i][j]].r==0)&&(AllPointXY[i][j][num_allpoints_xy[i][j]].g==0)&&(AllPointXY[i][j][num_allpoints_xy[i][j]].b==0))
//                                cout<<"1"<<endl;
                            num_allpoints_xy[i][j]++;
                            
                        }
                        
                    }
                }
            }
            
            if(dz>1)
                        {
                            GLfloat min_sign=999;
                            GLfloat max_sign=-999;
                            int min_sign_num = 0;
                            int max_sign_num = 0;
                            GLfloat proport_1,proport_2;
                            for(float k=1;k<dz;k=k+1)
                            {
                                min_sign=999;
                                max_sign=-999;
                                min_sign_num = 0;
                                max_sign_num = 0;
                                for(int l=0;l<num_edgepoints[i][j];l++)
                                {
                                    if((EdgePoint[i][j][l].z>=(face_min[i][j].z+k-1))&&(EdgePoint[i][j][l].z<(face_min[i][j].z+k)))
                                    {
                                        if(min_sign>EdgePoint[i][j][l].y)
                                        {
                                            min_sign = EdgePoint[i][j][l].y;
                                            min_sign_num = l;
                                            
                                        }
                                        if(max_sign<EdgePoint[i][j][l].y)
                                        {
                                            max_sign = EdgePoint[i][j][l].y;
                                            max_sign_num = l;
                                            
                                        }
                                        
                                    }
                                    
                                }
                                GLfloat dy_sign,step;
                                dy_sign = max_sign - min_sign;
                                step = dy_sign/dz/5;
                                if((min_sign!=999)&&(dy_sign>0))
                                {
                                    
                                    for(float l=step;l<(max_sign-min_sign);l=l+step)
                                    {
                                        AllPointYZ[i][j][num_allpoints_yz[i][j]].o1 =min_sign+l;
                                        AllPointYZ[i][j][num_allpoints_yz[i][j]].o2 =face_min[i][j].z+k-0.5;
                                        proport_1 = (max_sign-(min_sign+l))/(max_sign-min_sign);
                                        proport_2 = ((min_sign+l)-min_sign)/(max_sign-min_sign);
                                        AllPointYZ[i][j][num_allpoints_yz[i][j]].r = proport_1 * EdgePoint[i][j][min_sign_num].r + proport_2 * EdgePoint[i][j][max_sign_num].r;
                                        AllPointYZ[i][j][num_allpoints_yz[i][j]].g = proport_1 * EdgePoint[i][j][min_sign_num].g + proport_2 * EdgePoint[i][j][max_sign_num].g;
                                        AllPointYZ[i][j][num_allpoints_yz[i][j]].b = proport_1 * EdgePoint[i][j][min_sign_num].b + proport_2 * EdgePoint[i][j][max_sign_num].b;
                                        //cout<< AllPointYZ[i][j][num_allpoints_yz[i][j]].r <<'\t'<< AllPointYZ[i][j][num_allpoints_yz[i][j]].g <<'\t'<< AllPointYZ[i][j][num_allpoints_yz[i][j]].b <<endl;
                                        num_allpoints_yz[i][j]++;
                                        
                                    }
                                    
                                }
                            }
                        }
            
            if(dz>1)
                        {
                            GLfloat min_sign=999;
                            GLfloat max_sign=-999;
                            int min_sign_num = 0;
                            int max_sign_num = 0;
                            GLfloat proport_1,proport_2;
                            for(float k=1;k<dz;k=k+1)
                            {
                                min_sign=999;
                                max_sign=-999;
                                min_sign_num = 0;
                                max_sign_num = 0;
                                for(int l=0;l<num_edgepoints[i][j];l++)
                                {
                                    if((EdgePoint[i][j][l].z>=(face_min[i][j].z+k-1))&&(EdgePoint[i][j][l].z<(face_min[i][j].z+k)))
                                    {
                                        if(min_sign>EdgePoint[i][j][l].x)
                                        {
                                            min_sign = EdgePoint[i][j][l].x;
                                            min_sign_num = l;
                                            
                                        }
                                        if(max_sign<EdgePoint[i][j][l].x)
                                        {
                                            max_sign = EdgePoint[i][j][l].x;
                                            max_sign_num = l;
                                            
                                        }
                                        
                                    }
                                    
                                }
                                GLfloat dx_sign,step;
                                dx_sign = max_sign - min_sign;
                                step = dx_sign/dz/5;
                                if((min_sign!=999)&&(dx_sign>0))
                                {
                                    
                                    for(float l=step;l<(max_sign-min_sign);l=l+step)
                                    {
                                        AllPointXZ[i][j][num_allpoints_xz[i][j]].o1 =min_sign+l;
                                        AllPointXZ[i][j][num_allpoints_xz[i][j]].o2 =face_min[i][j].z+k-0.5;
                                        proport_1 = (max_sign-(min_sign+l))/(max_sign-min_sign);
                                        proport_2 = ((min_sign+l)-min_sign)/(max_sign-min_sign);
                                        AllPointXZ[i][j][num_allpoints_xz[i][j]].r = proport_1 * EdgePoint[i][j][min_sign_num].r + proport_2 * EdgePoint[i][j][max_sign_num].r;
                                        AllPointXZ[i][j][num_allpoints_xz[i][j]].g = proport_1 * EdgePoint[i][j][min_sign_num].g + proport_2 * EdgePoint[i][j][max_sign_num].g;
                                        AllPointXZ[i][j][num_allpoints_xz[i][j]].b = proport_1 * EdgePoint[i][j][min_sign_num].b + proport_2 * EdgePoint[i][j][max_sign_num].b;
//                                        if(AllPointXZ[i][j][num_allpoints_xz[i][j]].r==0&&AllPointXZ[i][j][num_allpoints_xz[i][j]].g==0&&AllPointXZ[i][j][num_allpoints_xz[i][j]].b==0)
//                                            cout<<"warning"<<endl;
                                        num_allpoints_xz[i][j]++;
                                        
                                    }
                                    
                                }
                            }
                        }
            for(int k=0;k<num_edgepoints[i][j];k++)
            {
                AllPointXY[i][j][num_allpoints_xy[i][j]].o1 = EdgePoint[i][j][k].x;
                AllPointXY[i][j][num_allpoints_xy[i][j]].o2 = EdgePoint[i][j][k].y;
                AllPointXY[i][j][num_allpoints_xy[i][j]].r = EdgePoint[i][j][k].r;
                AllPointXY[i][j][num_allpoints_xy[i][j]].g = EdgePoint[i][j][k].g;
                AllPointXY[i][j][num_allpoints_xy[i][j]].b = EdgePoint[i][j][k].b;
                num_allpoints_xy[i][j]++;
                
                AllPointYZ[i][j][num_allpoints_yz[i][j]].o1 = EdgePoint[i][j][k].y;
                AllPointYZ[i][j][num_allpoints_yz[i][j]].o2 = EdgePoint[i][j][k].z;
                AllPointYZ[i][j][num_allpoints_yz[i][j]].r = EdgePoint[i][j][k].r;
                AllPointYZ[i][j][num_allpoints_yz[i][j]].g = EdgePoint[i][j][k].g;
                AllPointYZ[i][j][num_allpoints_yz[i][j]].b = EdgePoint[i][j][k].b;
                num_allpoints_yz[i][j]++;
                
                AllPointXZ[i][j][num_allpoints_xz[i][j]].o1 = EdgePoint[i][j][k].x;
                AllPointXZ[i][j][num_allpoints_xz[i][j]].o2 = EdgePoint[i][j][k].z;
                AllPointXZ[i][j][num_allpoints_xz[i][j]].r = EdgePoint[i][j][k].r;
                AllPointXZ[i][j][num_allpoints_xz[i][j]].g = EdgePoint[i][j][k].g;
                AllPointXZ[i][j][num_allpoints_xz[i][j]].b = EdgePoint[i][j][k].b;
                num_allpoints_xz[i][j]++;
            }
//            cout<<"num all side points on this face  "<<num_edgepoints[i][j]<<endl;
//            cout<<"num all points on this face  "<<num_allpoints[i][j]<<endl;
//            cout<<"num add points on this face  "<<(num_allpoints[i][j] - num_edgepoints[i][j])<<endl;
//            cout<<endl;
        }
    
//    for(int i=0;i<num_poly;i++)
//    {
//        delete [] face_min[i];
//        delete [] face_max[i];
//    }
//    delete [] face_min;
//    delete [] face_max;
    
    
}



void get_face_normal_vector()
{
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
        {
            wcPt3D l1,l2;
            l1.x = Vertices[i][Faces[i][j].p2-1].x - Vertices[i][Faces[i][j].p1-1].x;
            l1.y = Vertices[i][Faces[i][j].p2-1].y - Vertices[i][Faces[i][j].p1-1].y;
            l1.z = Vertices[i][Faces[i][j].p2-1].z - Vertices[i][Faces[i][j].p1-1].z;
            l2.x = Vertices[i][Faces[i][j].p3-1].x - Vertices[i][Faces[i][j].p1-1].x;
            l2.y = Vertices[i][Faces[i][j].p3-1].y - Vertices[i][Faces[i][j].p1-1].y;
            l2.z = Vertices[i][Faces[i][j].p3-1].z - Vertices[i][Faces[i][j].p1-1].z;
            normal_vector[i][j].x = l1.y*l2.z-l1.z*l2.y;
            normal_vector[i][j].y = l1.z*l2.x-l1.x*l2.z;
            normal_vector[i][j].z = l1.x*l2.y-l1.y*l2.x;
            float temp;
            temp = powf(normal_vector[i][j].x, 2.0) + powf(normal_vector[i][j].y, 2.0) + powf(normal_vector[i][j].z, 2.0);
            temp = sqrtf(temp);
            normal_vector[i][j].x = normal_vector[i][j].x / temp;
            normal_vector[i][j].y = normal_vector[i][j].y / temp;
            normal_vector[i][j].z = normal_vector[i][j].z / temp;
        }
}

void get_point_normal_vector()
{
    
    for(int i=0;i<num_poly;i++)
    {
        for(int j=0;j<num_vertices[i];j++)
        {
            normal_vector_point[i][j].x = normal_vector_point[i][j].y = normal_vector_point[i][j].z = 0;
        }
        for(int j=0;j<num_faces[i];j++)
        {
            normal_vector_point[i][Faces[i][j].p1-1].x = normal_vector_point[i][Faces[i][j].p1-1].x + normal_vector[i][j].x;
            normal_vector_point[i][Faces[i][j].p1-1].y = normal_vector_point[i][Faces[i][j].p1-1].y + normal_vector[i][j].y;
            normal_vector_point[i][Faces[i][j].p1-1].z = normal_vector_point[i][Faces[i][j].p1-1].z + normal_vector[i][j].z;
            normal_vector_point[i][Faces[i][j].p2-1].x = normal_vector_point[i][Faces[i][j].p2-1].x + normal_vector[i][j].x;
            normal_vector_point[i][Faces[i][j].p2-1].y = normal_vector_point[i][Faces[i][j].p2-1].y + normal_vector[i][j].y;
            normal_vector_point[i][Faces[i][j].p2-1].z = normal_vector_point[i][Faces[i][j].p2-1].z + normal_vector[i][j].z;
            normal_vector_point[i][Faces[i][j].p3-1].x = normal_vector_point[i][Faces[i][j].p3-1].x + normal_vector[i][j].x;
            normal_vector_point[i][Faces[i][j].p3-1].y = normal_vector_point[i][Faces[i][j].p3-1].y + normal_vector[i][j].y;
            normal_vector_point[i][Faces[i][j].p3-1].z = normal_vector_point[i][Faces[i][j].p3-1].z + normal_vector[i][j].z;
        }
        for(int j=0;j<num_vertices[i];j++)
        {
            float temp;
            temp = powf(normal_vector_point[i][j].x, 2) + powf(normal_vector_point[i][j].y, 2) + powf(normal_vector_point[i][j].z, 2);
            temp = sqrtf(temp);
//            cout<<"temp"<<'\t'<<temp<<endl;
            normal_vector_point[i][j].x = normal_vector_point[i][j].x / temp;
            normal_vector_point[i][j].y = normal_vector_point[i][j].y / temp;
            normal_vector_point[i][j].z = normal_vector_point[i][j].z / temp;
        }
    }
}

void get_specu_vertice()
{
    faces_per_vertice = new int *[num_poly];
    specu_per_vertice = new GLfloat *[num_poly];
    for(int i=0;i<num_poly;i++)
    {
        faces_per_vertice[i] = new int[num_point_init];
        specu_per_vertice[i] = new GLfloat[num_point_init];
        for(int j=0;j<num_vertices[i];j++)
        {
            specu_per_vertice[i][j] = 0;
            faces_per_vertice[i][j] = 0;
        }
    }
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
        {
            specu_per_vertice[i][Faces[i][j].p1-1] = specu_per_vertice[i][Faces[i][j].p1-1] + Faces[i][j].s;
            specu_per_vertice[i][Faces[i][j].p2-1] = specu_per_vertice[i][Faces[i][j].p2-1] + Faces[i][j].s;
            specu_per_vertice[i][Faces[i][j].p3-1] = specu_per_vertice[i][Faces[i][j].p3-1] + Faces[i][j].s;
            faces_per_vertice[i][Faces[i][j].p1-1]++;
            faces_per_vertice[i][Faces[i][j].p2-1]++;
            faces_per_vertice[i][Faces[i][j].p3-1]++;
        }
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_vertices[i];j++)
        {
            specu_per_vertice[i][j] = specu_per_vertice[i][j] / faces_per_vertice[i][j];
        }
}

void phong_model_vertices()
{
    get_face_normal_vector();
    get_point_normal_vector();
    get_specu_vertice();
    cout<<"Please type in the color of ambient light ka.r, ka.g, ka.b (Example: 1 1 1 for white))"<<endl;
    cin>>ka.r>>ka.g>>ka.b;
//    ka.r = ka.g = ka.b =1;
    cout<<"Please type in IA for intensity of ambient light (Example: 0.4)"<<endl;
    cin>>IA;
//    IA = 0.5;
    cout<<"Please type in IL for intensity of light source (Example: 4)"<<endl;
    cin>>IL;
//    IL = 1;
    cout<<"Please type in K for distance that light source to the face (Example:10)"<<endl;
    cin>>K_distance;
//    K_distance = 10;
    cout<<"Please type in f for the location of eye (Example: 0.5 0.5 0.5)"<<endl;
    cin>>front_P.x>>front_P.y>>front_P.z;
//    front_P.x = (Pmax.x-Pmin.x)/2 + Pmin.x;
//    front_P.y = (Pmax.y-Pmin.y)/2 + Pmin.y;
//    front_P.z = 0;
    front_P.x = InitializeX(front_P.x);
    front_P.y = InitializeY(front_P.y);
    front_P.z = InitializeZ(front_P.z);
    cout<<"Please type in light source point (Example: 1 1 1)"<<endl;
    cin>>light_P.x>>light_P.y>>light_P.z;
//    light_P.x = 1;
//    light_P.y = 1;
//    light_P.z = 0;
    light_P.x = InitializeX(light_P.x);
    light_P.y = InitializeY(light_P.y);
    light_P.z = InitializeZ(light_P.z);
    cout<<"Please type in the color of light source ks.r ks.g ks.b (Example: 1 1 1 for white)"<<endl;
    cin>>ks.r>>ks.g>>ks.b;
//    ks.r = ks.g = ks.b = 1;
    
    phone_vertices = new  LightIntense *[num_poly];
    for(int i=0;i<num_poly;i++)
        phone_vertices[i] = new LightIntense[num_point_init];
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_vertices[i];j++)
        {
            wcPt3D vector_v,vector_r,vector_l,temp_s1,temp_s2;
            GLfloat f_minus_p_plus_k,f_length,l_length,temp_s1_length,l_multi_n,r_multi_v,r_multi_v_pow_n,n_multi_v;
            vector_v.x = front_P.x - Vertices[i][j].x;
            vector_v.y = front_P.y - Vertices[i][j].y;
            vector_v.z = front_P.z - Vertices[i][j].z;
            f_length = sqrtf((powf(vector_v.x, 2) + powf(vector_v.y, 2) + powf(vector_v.z, 2)));
            f_minus_p_plus_k = f_length/width;
            f_minus_p_plus_k = f_minus_p_plus_k + K_distance;
            vector_v.x = vector_v.x / f_length;
            vector_v.y = vector_v.y / f_length;
            vector_v.z = vector_v.z / f_length;
            
            vector_l.x = light_P.x - Vertices[i][j].x;
            vector_l.y = light_P.y - Vertices[i][j].y;
            vector_l.z = light_P.z - Vertices[i][j].z;
            l_length = sqrtf((powf(vector_l.x, 2) + powf(vector_l.y, 2) + powf(vector_l.z, 2)));
            vector_l.x = vector_l.x / l_length;
            vector_l.y = vector_l.y / l_length;
            vector_l.z = vector_l.z / l_length;
            
            temp_s1_length = vector_l.x * normal_vector_point[i][j].x + vector_l.y * normal_vector_point[i][j].y + vector_l.z * normal_vector_point[i][j].z;
            temp_s1.x = normal_vector_point[i][j].x * temp_s1_length;
            temp_s1.y = normal_vector_point[i][j].y * temp_s1_length;
            temp_s1.z = normal_vector_point[i][j].z * temp_s1_length;
            temp_s2.x = vector_l.x - temp_s1.x;
            temp_s2.y = vector_l.y - temp_s1.y;
            temp_s2.z = vector_l.z - temp_s1.z;
            vector_r.x = vector_l.x - 2*temp_s2.x;
            vector_r.y = vector_l.y - 2*temp_s2.y;
            vector_r.z = vector_l.z - 2*temp_s2.z;
//            cout<<"test"<<'\t'<<(pow(vector_r.x, 2) + pow(vector_r.y, 2) + pow(vector_r.z, 2))<<endl;
            
            l_multi_n = temp_s1_length;
            r_multi_v = vector_r.x * vector_v.x + vector_r.y * vector_v.y+ vector_r.z * vector_v.z;
//            cout<<r_multi_v<<endl;
//            cout<<specu_per_vertice[i][j]<<'\t';
            r_multi_v_pow_n = powf(r_multi_v, specu_per_vertice[i][j]);
//            cout<<r_multi_v_pow_n<<endl;
            n_multi_v = normal_vector_point[i][j].x * vector_v.x + normal_vector_point[i][j].y * vector_v.y + normal_vector_point[i][j].z * vector_v.z;
            
//            if(n_multi_v>=0)
//            {
                if((n_multi_v>0)&&(l_multi_n>0))
                {
                    if(r_multi_v<0)
                    {
                        phone_vertices[i][j].r = ka.r * IA + ((IL/(f_minus_p_plus_k)) * (Vertices[i][j].r * l_multi_n));
                        phone_vertices[i][j].g = ka.g * IA + ((IL/(f_minus_p_plus_k)) * (Vertices[i][j].g * l_multi_n));
                        phone_vertices[i][j].b = ka.b * IA + ((IL/(f_minus_p_plus_k)) * (Vertices[i][j].b * l_multi_n));
//                        if((phone_vertices[i][j].r)<0||phone_vertices[i][j].g<0||phone_vertices[i][j].b<0)
//                            cout<<"l_multi_n"<<'\t'<<l_multi_n<<"r_multi_v"<<'\t'<<r_multi_v<<endl;
                    }
                    else
                    {
                        phone_vertices[i][j].r = ka.r * IA + ((IL/(f_minus_p_plus_k)) * ((Vertices[i][j].r * l_multi_n) + (ks.r * r_multi_v_pow_n)));
                        phone_vertices[i][j].g = ka.g * IA + ((IL/(f_minus_p_plus_k)) * ((Vertices[i][j].g * l_multi_n) + (ks.g * r_multi_v_pow_n)));
                        phone_vertices[i][j].b = ka.b * IA + ((IL/(f_minus_p_plus_k)) * ((Vertices[i][j].b * l_multi_n) + (ks.b * r_multi_v_pow_n)));
//                        if((phone_vertices[i][j].r)<0||phone_vertices[i][j].g<0||phone_vertices[i][j].b<0)
//                        cout<<"l_multi_n"<<'\t'<<l_multi_n<<"r_multi_v"<<'\t'<<r_multi_v<<endl;
                        
                    }
                    
                }
                else
                {
                    phone_vertices[i][j].r = ka.r * IA;
                    phone_vertices[i][j].g = ka.g * IA;
                    phone_vertices[i][j].b = ka.b * IA;
                    
                }
            
            
            
//            phone_vertices[i][j].r = ka.r * IA + ((IL/(f_minus_p_plus_k)) * ((Vertices[i][j].r * l_multi_n) + (ks.r * r_multi_v_pow_n)));
//            phone_vertices[i][j].g = ka.g * IA + ((IL/(f_minus_p_plus_k)) * ((Vertices[i][j].g * l_multi_n) + (ks.g * r_multi_v_pow_n)));
//            phone_vertices[i][j].b = ka.b * IA + ((IL/(f_minus_p_plus_k)) * ((Vertices[i][j].b * l_multi_n) + (ks.b * r_multi_v_pow_n)));
                Vertices[i][j].r = phone_vertices[i][j].r;
                Vertices[i][j].g = phone_vertices[i][j].g;
                Vertices[i][j].b = phone_vertices[i][j].b;
//            }
//            else
//            {
//                Vertices[i][j].r = 0;
//                Vertices[i][j].g = 0;
//                Vertices[i][j].b = 0;
//            }
//            cout<<Vertices[i][j].r<<'\t'<<Vertices[i][j].g<<'\t'<<Vertices[i][j].b<<endl;
        }
    
    
}

GLfloat normalize_for_halftone(GLfloat v)
{
    v = v/(3*width+2)*width;
    return v;
}

void halftoning()
{
    HalftongXY = new ShowPoint **[num_poly];
    HalftongYZ = new ShowPoint **[num_poly];
    HalftongXZ = new ShowPoint **[num_poly];
    for(int i=0;i<num_poly;i++)
    {
        HalftongXY[i] = new ShowPoint *[num_face_init];
        HalftongYZ[i] = new ShowPoint *[num_face_init];
        HalftongXZ[i] = new ShowPoint *[num_face_init];
        for(int j=0;j<num_faces[i];j++)
        {
            HalftongXY[i][j] = new ShowPoint[9*num_point_init*num_point_init];
            HalftongYZ[i][j] = new ShowPoint[9*num_point_init*num_point_init];
            HalftongXZ[i][j] = new ShowPoint[9*num_point_init*num_point_init];
        }
    }
    
    srand(rand());
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
        {
            
            //xy case
            for(int k=0;k<num_allpoints_xy[i][j];k++)
            {
                //cordinator change
                HalftongXY[i][j][9*k].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1);
                HalftongXY[i][j][9*k].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2);
                HalftongXY[i][j][9*k+1].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1+1);
                HalftongXY[i][j][9*k+1].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2);
                HalftongXY[i][j][9*k+2].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1+2);
                HalftongXY[i][j][9*k+2].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2);
                HalftongXY[i][j][9*k+3].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1);
                HalftongXY[i][j][9*k+3].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2+1);
                HalftongXY[i][j][9*k+4].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1+1);
                HalftongXY[i][j][9*k+4].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2+1);
                HalftongXY[i][j][9*k+5].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1+2);
                HalftongXY[i][j][9*k+5].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2+1);
                HalftongXY[i][j][9*k+6].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1);
                HalftongXY[i][j][9*k+6].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2+2);
                HalftongXY[i][j][9*k+7].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1+1);
                HalftongXY[i][j][9*k+7].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2+2);
                HalftongXY[i][j][9*k+8].o1 = normalize_for_halftone(3*AllPointXY[i][j][k].o1+2);
                HalftongXY[i][j][9*k+8].o2 = normalize_for_halftone(3*AllPointXY[i][j][k].o2+2);
                
                //intensity
                GLfloat max_intensity=-0.9;
                if(max_intensity<AllPointXY[i][j][k].r) max_intensity = AllPointXY[i][j][k].r;
                if(max_intensity<AllPointXY[i][j][k].g) max_intensity = AllPointXY[i][j][k].g;
                if(max_intensity<AllPointXY[i][j][k].b) max_intensity = AllPointXY[i][j][k].b;
                GLfloat intensity_number = 0;
                GLfloat intensity_number_integer = 0;
                GLfloat intensity_number_possibility = 0;
                intensity_number = max_intensity * 9;
                intensity_number_integer = floorf(intensity_number);
                intensity_number_possibility = intensity_number - intensity_number_integer;
                intensity_number_possibility = floorf(intensity_number_possibility*10);
                
                //get a draw order
                
                int a[9];
                int b[9];
                int c[9];
                for(int m=0;m<9;m++)
                {
                    
                    a[m] = (rand() % (8));
//                    if((j==0)&&(k<3))
//                    cout<<"a[m]"<<a[m]<<'\t';
                    b[m] = 0;
                    c[m] = 0;
                }
                for(int m=0;m<9;m++)
                    for(int n=0;n<9;n++)
                        if(a[n]<a[m]) b[n]++;
//                if((j==0)&&(k<3))
//                    cout<<endl<<"j"<<j<<endl;
                for(int m=0;m<9;)
                {
                    int sign = m;
                    for(int n=0;n<9;n++)
                    {
                        if(b[n]==sign)
                        {
                            c[m] = n;
//                            if((j==0)&&(k<3))
//                                cout<<"c[m]"<<'\t'<<c[m]<<'\t';
                            m++;
                        }
                    }
                }
                int rand_sign = (rand() % (9));
                if(halftone_sign == 1)
                {
                if(rand_sign<=(intensity_number_possibility-1))
                {
                    for(int m=0;m<intensity_number_integer+1;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 1;
                        HalftongXY[i][j][9*k+c[m]].g = 1;
                        HalftongXY[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_integer+1;m<9;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 0;
                        HalftongXY[i][j][9*k+c[m]].g = 0;
                        HalftongXY[i][j][9*k+c[m]].b = 0;
                    }
                }
                else
                {
                    for(int m=0;m<intensity_number_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 1;
                        HalftongXY[i][j][9*k+c[m]].g = 1;
                        HalftongXY[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_integer;m<9;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 0;
                        HalftongXY[i][j][9*k+c[m]].g = 0;
                        HalftongXY[i][j][9*k+c[m]].b = 0;
                    }
                }
                }
                
                if(halftone_sign==2)
                {
                
                if(rand_sign<=(intensity_number_possibility-1))
                {
                    GLfloat intensity_number_r = 0;
                    GLfloat intensity_number_g = 0;
                    GLfloat intensity_number_b = 0;
                    GLfloat intensity_number_r_integer = 0;
                    GLfloat intensity_number_g_integer = 0;
                    GLfloat intensity_number_b_integer = 0;
                    GLfloat intensity_number_r_possibility = 0;
                    GLfloat intensity_number_g_possibility = 0;
                    GLfloat intensity_number_b_possibility = 0;
                    intensity_number_r = AllPointXY[i][j][k].r/(AllPointXY[i][j][k].r+AllPointXY[i][j][k].g+AllPointXY[i][j][k].b) * (intensity_number_integer+1);
                    intensity_number_g = AllPointXY[i][j][k].g/(AllPointXY[i][j][k].r+AllPointXY[i][j][k].g+AllPointXY[i][j][k].b) * (intensity_number_integer+1);
                    intensity_number_b = AllPointXY[i][j][k].b/(AllPointXY[i][j][k].r+AllPointXY[i][j][k].g+AllPointXY[i][j][k].b) * (intensity_number_integer+1);
                    intensity_number_r_integer = floorf(intensity_number_r);
                    intensity_number_g_integer = floorf(intensity_number_g);
                    intensity_number_b_integer = floorf(intensity_number_b);
                    intensity_number_r_possibility = intensity_number_r - intensity_number_r_integer;
                    intensity_number_g_possibility = intensity_number_g - intensity_number_g_integer;
                    intensity_number_b_possibility = intensity_number_b - intensity_number_b_integer;
                    intensity_number_r_possibility = floorf(intensity_number_r_possibility*10);
                    intensity_number_g_possibility = floorf(intensity_number_g_possibility*10);
                    intensity_number_b_possibility = floorf(intensity_number_b_possibility*10);
                    int temp;
                    temp = intensity_number_r_possibility + intensity_number_g_possibility + intensity_number_b_possibility;
                    if(temp!=0)
                    {
                        intensity_number_r_possibility = intensity_number_r_possibility / temp;
                        intensity_number_g_possibility = intensity_number_g_possibility / temp;
                        intensity_number_b_possibility = intensity_number_b_possibility / temp;
                    }

                    for(int m=0;m<intensity_number_r_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 1;
                        HalftongXY[i][j][9*k+c[m]].g = 0;
                        HalftongXY[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer;m<intensity_number_r_integer + intensity_number_g_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 0;
                        HalftongXY[i][j][9*k+c[m]].g = 1;
                        HalftongXY[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer;m<intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 0;
                        HalftongXY[i][j][9*k+c[m]].g = 0;
                        HalftongXY[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m<9;m++)
                    {
                        float poss = 0;
                        poss = rand() / double(RAND_MAX);;
                        if(poss<intensity_number_r_possibility)
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 1;
                            HalftongXY[i][j][9*k+c[m]].g = 0;
                            HalftongXY[i][j][9*k+c[m]].b = 0;
                        }
                        else if(poss<intensity_number_r_possibility+intensity_number_g_possibility)
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 0;
                            HalftongXY[i][j][9*k+c[m]].g = 1;
                            HalftongXY[i][j][9*k+c[m]].b = 0;
                        }
                        else if(poss<intensity_number_r_possibility+intensity_number_g_possibility+intensity_number_b_possibility)
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 0;
                            HalftongXY[i][j][9*k+c[m]].g = 0;
                            HalftongXY[i][j][9*k+c[m]].b = 1;
                        }
                        else
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 0;
                            HalftongXY[i][j][9*k+c[m]].g = 0;
                            HalftongXY[i][j][9*k+c[m]].b = 0;
                        }
                    }

                }
                else
                {
                    GLfloat intensity_number_r = 0;
                    GLfloat intensity_number_g = 0;
                    GLfloat intensity_number_b = 0;
                    GLfloat intensity_number_r_integer = 0;
                    GLfloat intensity_number_g_integer = 0;
                    GLfloat intensity_number_b_integer = 0;
                    GLfloat intensity_number_r_possibility = 0;
                    GLfloat intensity_number_g_possibility = 0;
                    GLfloat intensity_number_b_possibility = 0;
                    intensity_number_r = AllPointXY[i][j][k].r/(AllPointXY[i][j][k].r+AllPointXY[i][j][k].g+AllPointXY[i][j][k].b) * (intensity_number_integer);
                    intensity_number_g = AllPointXY[i][j][k].g/(AllPointXY[i][j][k].r+AllPointXY[i][j][k].g+AllPointXY[i][j][k].b) * (intensity_number_integer);
                    intensity_number_b = AllPointXY[i][j][k].b/(AllPointXY[i][j][k].r+AllPointXY[i][j][k].g+AllPointXY[i][j][k].b) * (intensity_number_integer);
                    intensity_number_r_integer = floorf(intensity_number_r);
                    intensity_number_g_integer = floorf(intensity_number_g);
                    intensity_number_b_integer = floorf(intensity_number_b);
                    intensity_number_r_possibility = intensity_number_r - intensity_number_r_integer;
                    intensity_number_g_possibility = intensity_number_g - intensity_number_g_integer;
                    intensity_number_b_possibility = intensity_number_b - intensity_number_b_integer;
                    intensity_number_r_possibility = floorf(intensity_number_r_possibility*10);
                    intensity_number_g_possibility = floorf(intensity_number_g_possibility*10);
                    intensity_number_b_possibility = floorf(intensity_number_b_possibility*10);
                    int temp;
                    temp = intensity_number_r_possibility + intensity_number_g_possibility + intensity_number_b_possibility;
                    if(temp!=0)
                    {
                        intensity_number_r_possibility = intensity_number_r_possibility / temp;
                        intensity_number_g_possibility = intensity_number_g_possibility / temp;
                        intensity_number_b_possibility = intensity_number_b_possibility / temp;
                    }

                    for(int m=0;m<intensity_number_r_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 1;
                        HalftongXY[i][j][9*k+c[m]].g = 0;
                        HalftongXY[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer;m<intensity_number_r_integer + intensity_number_g_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 0;
                        HalftongXY[i][j][9*k+c[m]].g = 1;
                        HalftongXY[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer;m<intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m++)
                    {
                        HalftongXY[i][j][9*k+c[m]].r = 0;
                        HalftongXY[i][j][9*k+c[m]].g = 0;
                        HalftongXY[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m<9;m++)
                    {
                        float poss = 0;
                        poss = rand() / double(RAND_MAX);;
                        if(poss<intensity_number_r_possibility)
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 1;
                            HalftongXY[i][j][9*k+c[m]].g = 0;
                            HalftongXY[i][j][9*k+c[m]].b = 0;
                        }
                        else if(poss<intensity_number_r_possibility+intensity_number_g_possibility)
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 0;
                            HalftongXY[i][j][9*k+c[m]].g = 1;
                            HalftongXY[i][j][9*k+c[m]].b = 0;
                        }
                        else
                        {
                            HalftongXY[i][j][9*k+c[m]].r = 0;
                            HalftongXY[i][j][9*k+c[m]].g = 0;
                            HalftongXY[i][j][9*k+c[m]].b = 1;
                        }
                    }
                }
                }
                
            }
    
            
            //yz case
    for(int k=0;k<num_allpoints_yz[i][j];k++)
    {
        //cordinator change
        HalftongYZ[i][j][9*k].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1);
        HalftongYZ[i][j][9*k].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2);
        HalftongYZ[i][j][9*k+1].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1+1);
        HalftongYZ[i][j][9*k+1].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2);
        HalftongYZ[i][j][9*k+2].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1+2);
        HalftongYZ[i][j][9*k+2].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2);
        HalftongYZ[i][j][9*k+3].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1);
        HalftongYZ[i][j][9*k+3].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2+1);
        HalftongYZ[i][j][9*k+4].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1+1);
        HalftongYZ[i][j][9*k+4].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2+1);
        HalftongYZ[i][j][9*k+5].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1+2);
        HalftongYZ[i][j][9*k+5].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2+1);
        HalftongYZ[i][j][9*k+6].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1);
        HalftongYZ[i][j][9*k+6].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2+2);
        HalftongYZ[i][j][9*k+7].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1+1);
        HalftongYZ[i][j][9*k+7].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2+2);
        HalftongYZ[i][j][9*k+8].o1 = normalize_for_halftone(3*AllPointYZ[i][j][k].o1+2);
        HalftongYZ[i][j][9*k+8].o2 = normalize_for_halftone(3*AllPointYZ[i][j][k].o2+2);
        
        //intensity
        GLfloat max_intensity=-0.9;
        if(max_intensity<AllPointYZ[i][j][k].r) max_intensity = AllPointYZ[i][j][k].r;
        if(max_intensity<AllPointYZ[i][j][k].g) max_intensity = AllPointYZ[i][j][k].g;
        if(max_intensity<AllPointYZ[i][j][k].b) max_intensity = AllPointYZ[i][j][k].b;
        GLfloat intensity_number = 0;
        GLfloat intensity_number_integer = 0;
        GLfloat intensity_number_possibility = 0;
        intensity_number = max_intensity * 9;
        intensity_number_integer = floorf(intensity_number);
        intensity_number_possibility = intensity_number - intensity_number_integer;
        intensity_number_possibility = floorf(intensity_number_possibility*10);
        
        //get a draw order
        int a[9];
        int b[9];
        int c[9];
        for(int m=0;m<9;m++)
        {
            
            a[m] = (rand() % (8));
            b[m] = 0;
            c[m] = 0;
        }
        for(int m=0;m<9;m++)
            for(int n=0;n<9;n++)
                if(a[n]<a[m]) b[n]++;
        for(int m=0;m<9;)
        {
            int sign = m;
            for(int n=0;n<9;n++)
            {
                if(b[n]==sign)
                {
                    c[m] = n;
                    m++;
                }
            }
        }
        int rand_sign = (rand() % (9));
        if(halftone_sign == 1)
        {
        if(rand_sign<=(intensity_number_possibility-1))
        {
            for(int m=0;m<intensity_number_integer+1;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 1;
                HalftongYZ[i][j][9*k+c[m]].g = 1;
                HalftongYZ[i][j][9*k+c[m]].b = 1;
            }
            for(int m=intensity_number_integer+1;m<9;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 0;
                HalftongYZ[i][j][9*k+c[m]].g = 0;
                HalftongYZ[i][j][9*k+c[m]].b = 0;
            }
        }
        else
        {
            for(int m=0;m<intensity_number_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 1;
                HalftongYZ[i][j][9*k+c[m]].g = 1;
                HalftongYZ[i][j][9*k+c[m]].b = 1;
            }
            for(int m=intensity_number_integer;m<9;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 0;
                HalftongYZ[i][j][9*k+c[m]].g = 0;
                HalftongYZ[i][j][9*k+c[m]].b = 0;
            }
        }
        }
        
        if(halftone_sign==2)
        {
        
        if(rand_sign<=(intensity_number_possibility-1))
        {
            GLfloat intensity_number_r = 0;
            GLfloat intensity_number_g = 0;
            GLfloat intensity_number_b = 0;
            GLfloat intensity_number_r_integer = 0;
            GLfloat intensity_number_g_integer = 0;
            GLfloat intensity_number_b_integer = 0;
            GLfloat intensity_number_r_possibility = 0;
            GLfloat intensity_number_g_possibility = 0;
            GLfloat intensity_number_b_possibility = 0;
            intensity_number_r = AllPointYZ[i][j][k].r/(AllPointYZ[i][j][k].r+AllPointYZ[i][j][k].g+AllPointYZ[i][j][k].b) * (intensity_number_integer+1);
            intensity_number_g = AllPointYZ[i][j][k].g/(AllPointYZ[i][j][k].r+AllPointYZ[i][j][k].g+AllPointYZ[i][j][k].b) * (intensity_number_integer+1);
            intensity_number_b = AllPointYZ[i][j][k].b/(AllPointYZ[i][j][k].r+AllPointYZ[i][j][k].g+AllPointYZ[i][j][k].b) * (intensity_number_integer+1);
            intensity_number_r_integer = floorf(intensity_number_r);
            intensity_number_g_integer = floorf(intensity_number_g);
            intensity_number_b_integer = floorf(intensity_number_b);
            intensity_number_r_possibility = intensity_number_r - intensity_number_r_integer;
            intensity_number_g_possibility = intensity_number_g - intensity_number_g_integer;
            intensity_number_b_possibility = intensity_number_b - intensity_number_b_integer;
            intensity_number_r_possibility = floorf(intensity_number_r_possibility*10);
            intensity_number_g_possibility = floorf(intensity_number_g_possibility*10);
            intensity_number_b_possibility = floorf(intensity_number_b_possibility*10);
            int temp;
            temp = intensity_number_r_possibility + intensity_number_g_possibility + intensity_number_b_possibility;
            if(temp!=0)
            {
                intensity_number_r_possibility = intensity_number_r_possibility / temp;
                intensity_number_g_possibility = intensity_number_g_possibility / temp;
                intensity_number_b_possibility = intensity_number_b_possibility / temp;
            }

            for(int m=0;m<intensity_number_r_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 1;
                HalftongYZ[i][j][9*k+c[m]].g = 0;
                HalftongYZ[i][j][9*k+c[m]].b = 0;
            }
            for(int m=intensity_number_r_integer;m<intensity_number_r_integer + intensity_number_g_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 0;
                HalftongYZ[i][j][9*k+c[m]].g = 1;
                HalftongYZ[i][j][9*k+c[m]].b = 0;
            }
            for(int m=intensity_number_r_integer + intensity_number_g_integer;m<intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 0;
                HalftongYZ[i][j][9*k+c[m]].g = 0;
                HalftongYZ[i][j][9*k+c[m]].b = 1;
            }
            for(int m=intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m<9;m++)
            {
                float poss = 0;
                poss = rand() / double(RAND_MAX);;
                if(poss<intensity_number_r_possibility)
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 1;
                    HalftongYZ[i][j][9*k+c[m]].g = 0;
                    HalftongYZ[i][j][9*k+c[m]].b = 0;
                }
                else if(poss<intensity_number_r_possibility+intensity_number_g_possibility)
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 0;
                    HalftongYZ[i][j][9*k+c[m]].g = 1;
                    HalftongYZ[i][j][9*k+c[m]].b = 0;
                }
                else if(poss<intensity_number_r_possibility+intensity_number_g_possibility+intensity_number_b_possibility)
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 0;
                    HalftongYZ[i][j][9*k+c[m]].g = 0;
                    HalftongYZ[i][j][9*k+c[m]].b = 1;
                }
                else
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 0;
                    HalftongYZ[i][j][9*k+c[m]].g = 0;
                    HalftongYZ[i][j][9*k+c[m]].b = 0;
                }
            }

        }
        else
        {
            GLfloat intensity_number_r = 0;
            GLfloat intensity_number_g = 0;
            GLfloat intensity_number_b = 0;
            GLfloat intensity_number_r_integer = 0;
            GLfloat intensity_number_g_integer = 0;
            GLfloat intensity_number_b_integer = 0;
            GLfloat intensity_number_r_possibility = 0;
            GLfloat intensity_number_g_possibility = 0;
            GLfloat intensity_number_b_possibility = 0;
            intensity_number_r = AllPointYZ[i][j][k].r/(AllPointYZ[i][j][k].r+AllPointYZ[i][j][k].g+AllPointYZ[i][j][k].b) * (intensity_number_integer);
            intensity_number_g = AllPointYZ[i][j][k].g/(AllPointYZ[i][j][k].r+AllPointYZ[i][j][k].g+AllPointYZ[i][j][k].b) * (intensity_number_integer);
            intensity_number_b = AllPointYZ[i][j][k].b/(AllPointYZ[i][j][k].r+AllPointYZ[i][j][k].g+AllPointYZ[i][j][k].b) * (intensity_number_integer);
            intensity_number_r_integer = floorf(intensity_number_r);
            intensity_number_g_integer = floorf(intensity_number_g);
            intensity_number_b_integer = floorf(intensity_number_b);
            intensity_number_r_possibility = intensity_number_r - intensity_number_r_integer;
            intensity_number_g_possibility = intensity_number_g - intensity_number_g_integer;
            intensity_number_b_possibility = intensity_number_b - intensity_number_b_integer;
            intensity_number_r_possibility = floorf(intensity_number_r_possibility*10);
            intensity_number_g_possibility = floorf(intensity_number_g_possibility*10);
            intensity_number_b_possibility = floorf(intensity_number_b_possibility*10);
            int temp;
            temp = intensity_number_r_possibility + intensity_number_g_possibility + intensity_number_b_possibility;
            if(temp!=0)
            {
                intensity_number_r_possibility = intensity_number_r_possibility / temp;
                intensity_number_g_possibility = intensity_number_g_possibility / temp;
                intensity_number_b_possibility = intensity_number_b_possibility / temp;
            }

            for(int m=0;m<intensity_number_r_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 1;
                HalftongYZ[i][j][9*k+c[m]].g = 0;
                HalftongYZ[i][j][9*k+c[m]].b = 0;
            }
            for(int m=intensity_number_r_integer;m<intensity_number_r_integer + intensity_number_g_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 0;
                HalftongYZ[i][j][9*k+c[m]].g = 1;
                HalftongYZ[i][j][9*k+c[m]].b = 0;
            }
            for(int m=intensity_number_r_integer + intensity_number_g_integer;m<intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m++)
            {
                HalftongYZ[i][j][9*k+c[m]].r = 0;
                HalftongYZ[i][j][9*k+c[m]].g = 0;
                HalftongYZ[i][j][9*k+c[m]].b = 1;
            }
            for(int m=intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m<9;m++)
            {
                float poss = 0;
                poss = rand() / double(RAND_MAX);;
                if(poss<intensity_number_r_possibility)
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 1;
                    HalftongYZ[i][j][9*k+c[m]].g = 0;
                    HalftongYZ[i][j][9*k+c[m]].b = 0;
                }
                else if(poss<intensity_number_r_possibility+intensity_number_g_possibility)
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 0;
                    HalftongYZ[i][j][9*k+c[m]].g = 1;
                    HalftongYZ[i][j][9*k+c[m]].b = 0;
                }
                else
                {
                    HalftongYZ[i][j][9*k+c[m]].r = 0;
                    HalftongYZ[i][j][9*k+c[m]].g = 0;
                    HalftongYZ[i][j][9*k+c[m]].b = 1;
                }
            }
        }
        }
        
    }
            
            //xz
            for(int k=0;k<num_allpoints_xz[i][j];k++)
            {
                //cordinator change
                HalftongXZ[i][j][9*k].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1);
                HalftongXZ[i][j][9*k].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2);
                HalftongXZ[i][j][9*k+1].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1+1);
                HalftongXZ[i][j][9*k+1].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2);
                HalftongXZ[i][j][9*k+2].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1+2);
                HalftongXZ[i][j][9*k+2].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2);
                HalftongXZ[i][j][9*k+3].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1);
                HalftongXZ[i][j][9*k+3].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2+1);
                HalftongXZ[i][j][9*k+4].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1+1);
                HalftongXZ[i][j][9*k+4].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2+1);
                HalftongXZ[i][j][9*k+5].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1+2);
                HalftongXZ[i][j][9*k+5].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2+1);
                HalftongXZ[i][j][9*k+6].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1);
                HalftongXZ[i][j][9*k+6].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2+2);
                HalftongXZ[i][j][9*k+7].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1+1);
                HalftongXZ[i][j][9*k+7].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2+2);
                HalftongXZ[i][j][9*k+8].o1 = normalize_for_halftone(3*AllPointXZ[i][j][k].o1+2);
                HalftongXZ[i][j][9*k+8].o2 = normalize_for_halftone(3*AllPointXZ[i][j][k].o2+2);
                
                //intensity
                GLfloat max_intensity=-0.9;
                if(max_intensity<AllPointXZ[i][j][k].r) max_intensity = AllPointXZ[i][j][k].r;
                if(max_intensity<AllPointXZ[i][j][k].g) max_intensity = AllPointXZ[i][j][k].g;
                if(max_intensity<AllPointXZ[i][j][k].b) max_intensity = AllPointXZ[i][j][k].b;
                GLfloat intensity_number = 0;
                GLfloat intensity_number_integer = 0;
                GLfloat intensity_number_possibility = 0;
                intensity_number = max_intensity * 9;
                intensity_number_integer = floorf(intensity_number);
                intensity_number_possibility = intensity_number - intensity_number_integer;
                intensity_number_possibility = floorf(intensity_number_possibility*10);
                //if(intensity_number_integer==0)
                    //cout<<AllPointXZ[i][j][k].r<<'\t'<<AllPointXZ[i][j][k].g<<'\t'<<AllPointXZ[i][j][k].b<<'\t'<<endl;
                //get a draw order
                int a[9];
                int b[9];
                int c[9];
                for(int m=0;m<9;m++)
                {
                    
                    a[m] = (rand() % (8));
                    b[m] = 0;
                    c[m] = 0;
                }
                for(int m=0;m<9;m++)
                    for(int n=0;n<9;n++)
                        if(a[n]<a[m]) b[n]++;
                
                for(int m=0;m<9;)
                {
                    int sign = m;
                    for(int n=0;n<9;n++)
                    {
                        if(b[n]==sign)
                        {
                            c[m] = n;
                            
                            m++;
                        }
                    }
                }
                int rand_sign = (rand() % (9));
                if(halftone_sign == 1)
                {
                if(rand_sign<=(intensity_number_possibility-1))
                {
                    for(int m=0;m<intensity_number_integer+1;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 1;
                        HalftongXZ[i][j][9*k+c[m]].g = 1;
                        HalftongXZ[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_integer+1;m<9;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 0;
                        HalftongXZ[i][j][9*k+c[m]].g = 0;
                        HalftongXZ[i][j][9*k+c[m]].b = 0;
                    }
                }
                else
                {
                    for(int m=0;m<intensity_number_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 1;
                        HalftongXZ[i][j][9*k+c[m]].g = 1;
                        HalftongXZ[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_integer;m<9;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 0;
                        HalftongXZ[i][j][9*k+c[m]].g = 0;
                        HalftongXZ[i][j][9*k+c[m]].b = 0;
                    }
                }
                }
                
                if(halftone_sign==2)
                {
                
                if(rand_sign<=(intensity_number_possibility-1))
                {
                    GLfloat intensity_number_r = 0;
                    GLfloat intensity_number_g = 0;
                    GLfloat intensity_number_b = 0;
                    GLfloat intensity_number_r_integer = 0;
                    GLfloat intensity_number_g_integer = 0;
                    GLfloat intensity_number_b_integer = 0;
                    GLfloat intensity_number_r_possibility = 0;
                    GLfloat intensity_number_g_possibility = 0;
                    GLfloat intensity_number_b_possibility = 0;
                    intensity_number_r = AllPointXZ[i][j][k].r/(AllPointXZ[i][j][k].r+AllPointXZ[i][j][k].g+AllPointXZ[i][j][k].b) * (intensity_number_integer+1);
                    intensity_number_g = AllPointXZ[i][j][k].g/(AllPointXZ[i][j][k].r+AllPointXZ[i][j][k].g+AllPointXZ[i][j][k].b) * (intensity_number_integer+1);
                    intensity_number_b = AllPointXZ[i][j][k].b/(AllPointXZ[i][j][k].r+AllPointXZ[i][j][k].g+AllPointXZ[i][j][k].b) * (intensity_number_integer+1);
                    intensity_number_r_integer = floorf(intensity_number_r);
                    intensity_number_g_integer = floorf(intensity_number_g);
                    intensity_number_b_integer = floorf(intensity_number_b);
                    intensity_number_r_possibility = intensity_number_r - intensity_number_r_integer;
                    intensity_number_g_possibility = intensity_number_g - intensity_number_g_integer;
                    intensity_number_b_possibility = intensity_number_b - intensity_number_b_integer;
                    intensity_number_r_possibility = floorf(intensity_number_r_possibility*10);
                    intensity_number_g_possibility = floorf(intensity_number_g_possibility*10);
                    intensity_number_b_possibility = floorf(intensity_number_b_possibility*10);
                    int temp;
                    temp = intensity_number_r_possibility + intensity_number_g_possibility + intensity_number_b_possibility;
                    if(temp!=0)
                    {
                        intensity_number_r_possibility = intensity_number_r_possibility / temp;
                        intensity_number_g_possibility = intensity_number_g_possibility / temp;
                        intensity_number_b_possibility = intensity_number_b_possibility / temp;
                    }

                    for(int m=0;m<intensity_number_r_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 1;
                        HalftongXZ[i][j][9*k+c[m]].g = 0;
                        HalftongXZ[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer;m<intensity_number_r_integer + intensity_number_g_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 0;
                        HalftongXZ[i][j][9*k+c[m]].g = 1;
                        HalftongXZ[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer;m<intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 0;
                        HalftongXZ[i][j][9*k+c[m]].g = 0;
                        HalftongXZ[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m<9;m++)
                    {
                        float poss = 0;
                        poss = rand() / double(RAND_MAX);;
                        if(poss<intensity_number_r_possibility)
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 1;
                            HalftongXZ[i][j][9*k+c[m]].g = 0;
                            HalftongXZ[i][j][9*k+c[m]].b = 0;
                        }
                        else if(poss<intensity_number_r_possibility+intensity_number_g_possibility)
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 0;
                            HalftongXZ[i][j][9*k+c[m]].g = 1;
                            HalftongXZ[i][j][9*k+c[m]].b = 0;
                        }
                        else if(poss<intensity_number_r_possibility+intensity_number_g_possibility+intensity_number_b_possibility)
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 0;
                            HalftongXZ[i][j][9*k+c[m]].g = 0;
                            HalftongXZ[i][j][9*k+c[m]].b = 1;
                        }
                        else if(poss<1)
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 0;
                            HalftongXZ[i][j][9*k+c[m]].g = 0;
                            HalftongXZ[i][j][9*k+c[m]].b = 0;
                        }
                    }

                }
                else
                {
                    GLfloat intensity_number_r = 0;
                    GLfloat intensity_number_g = 0;
                    GLfloat intensity_number_b = 0;
                    GLfloat intensity_number_r_integer = 0;
                    GLfloat intensity_number_g_integer = 0;
                    GLfloat intensity_number_b_integer = 0;
                    GLfloat intensity_number_r_possibility = 0;
                    GLfloat intensity_number_g_possibility = 0;
                    GLfloat intensity_number_b_possibility = 0;
                    intensity_number_r = AllPointXZ[i][j][k].r/(AllPointXZ[i][j][k].r+AllPointXZ[i][j][k].g+AllPointXZ[i][j][k].b) * (intensity_number_integer);
                    intensity_number_g = AllPointXZ[i][j][k].g/(AllPointXZ[i][j][k].r+AllPointXZ[i][j][k].g+AllPointXZ[i][j][k].b) * (intensity_number_integer);
                    intensity_number_b = AllPointXZ[i][j][k].b/(AllPointXZ[i][j][k].r+AllPointXZ[i][j][k].g+AllPointXZ[i][j][k].b) * (intensity_number_integer);
                    intensity_number_r_integer = floorf(intensity_number_r);
                    intensity_number_g_integer = floorf(intensity_number_g);
                    intensity_number_b_integer = floorf(intensity_number_b);
                    intensity_number_r_possibility = intensity_number_r - intensity_number_r_integer;
                    intensity_number_g_possibility = intensity_number_g - intensity_number_g_integer;
                    intensity_number_b_possibility = intensity_number_b - intensity_number_b_integer;
                    intensity_number_r_possibility = floorf(intensity_number_r_possibility*10);
                    intensity_number_g_possibility = floorf(intensity_number_g_possibility*10);
                    intensity_number_b_possibility = floorf(intensity_number_b_possibility*10);
                    int temp;
                    temp = intensity_number_r_possibility + intensity_number_g_possibility + intensity_number_b_possibility;
                    if(temp!=0)
                    {
                        intensity_number_r_possibility = intensity_number_r_possibility / temp;
                        intensity_number_g_possibility = intensity_number_g_possibility / temp;
                        intensity_number_b_possibility = intensity_number_b_possibility / temp;
                    }

                    for(int m=0;m<intensity_number_r_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 1;
                        HalftongXZ[i][j][9*k+c[m]].g = 0;
                        HalftongXZ[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer;m<intensity_number_r_integer + intensity_number_g_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 0;
                        HalftongXZ[i][j][9*k+c[m]].g = 1;
                        HalftongXZ[i][j][9*k+c[m]].b = 0;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer;m<intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m++)
                    {
                        HalftongXZ[i][j][9*k+c[m]].r = 0;
                        HalftongXZ[i][j][9*k+c[m]].g = 0;
                        HalftongXZ[i][j][9*k+c[m]].b = 1;
                    }
                    for(int m=intensity_number_r_integer + intensity_number_g_integer + intensity_number_b_integer;m<9;m++)
                    {
                        float poss = 0;
                        poss = rand() / double(RAND_MAX);;
                        if(poss<intensity_number_r_possibility)
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 1;
                            HalftongXZ[i][j][9*k+c[m]].g = 0;
                            HalftongXZ[i][j][9*k+c[m]].b = 0;
                        }
                        else if(poss<intensity_number_r_possibility+intensity_number_g_possibility)
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 0;
                            HalftongXZ[i][j][9*k+c[m]].g = 1;
                            HalftongXZ[i][j][9*k+c[m]].b = 0;
                        }
                        else
                        {
                            HalftongXZ[i][j][9*k+c[m]].r = 0;
                            HalftongXZ[i][j][9*k+c[m]].g = 0;
                            HalftongXZ[i][j][9*k+c[m]].b = 1;
                        }
                    }
                }
                }
                
            }
        }
}



void DrawPoly()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, win_width, win_height);
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINES);
    glVertex2f(width/2, 0);
    glVertex2f(width/2, length);
    glVertex2f(0, length/2);
    glVertex2f(width, length/2);
    glEnd();

//// draw vertice
//    //xy, up left
//    glViewport(0, win_height, win_width, win_height);
//    glBegin(GL_POINTS);
//    for(int i=0;i<num_poly;i++)
//        for(int j=0;j<num_vertices[i];j++)
////            for(int k=0;k<num_edgepoints[i][j];k++)
//            {
//                glColor3f(Vertices[i][j].r, Vertices[i][j].g, Vertices[i][j].b);
//                glVertex2f(Vertices[i][j].x,Vertices[i][j].y);
//            }
//    glEnd();
//
//    //yz, bottom left
//    glViewport(0, 0, win_width, win_height);
//    glBegin(GL_POINTS);
//    for(int i=0;i<num_poly;i++)
//        for(int j=0;j<num_vertices[i];j++)
////            for(int k=0;k<num_edgepoints[i][j];k++)
//            {
//                glColor3f(Vertices[i][j].r, Vertices[i][j].g, Vertices[i][j].b);
//                glVertex2f(Vertices[i][j].y, Vertices[i][j].z);
//            }
//    glEnd();
//
//    //xz, up right
//    glViewport(win_width, win_height, win_width, win_height);
//    glBegin(GL_POINTS);
//    for(int i=0;i<num_poly;i++)
//        for(int j=0;j<num_vertices[i];j++)
////            for(int k=0;k<num_edgepoints[i][j];k++)
//            {
//                glColor3f(Vertices[i][j].r, Vertices[i][j].g, Vertices[i][j].b);
//                glVertex2f(Vertices[i][j].x, Vertices[i][j].z);
//            }
//    glEnd();

/*
//  draw side point
    //xy, up left
    glViewport(0, win_height/2, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_edgepoints[i][z_buffer_face[i][j]];k++)
            {
                glColor3f(EdgePoint[i][z_buffer_face[i][j]][k].r, EdgePoint[i][z_buffer_face[i][j]][k].g, EdgePoint[i][z_buffer_face[i][j]][k].b);
                glVertex2f(EdgePoint[i][z_buffer_face[i][j]][k].x,EdgePoint[i][z_buffer_face[i][j]][k].y);
            }
    glEnd();

    //yz, bottom left
    glViewport(0, 0, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_edgepoints[i][x_buffer_face[i][j]];k++)
            {
                glColor3f(EdgePoint[i][x_buffer_face[i][j]][k].r, EdgePoint[i][x_buffer_face[i][j]][k].g, EdgePoint[i][x_buffer_face[i][j]][k].b);
                glVertex2f(EdgePoint[i][x_buffer_face[i][j]][k].y, EdgePoint[i][x_buffer_face[i][j]][k].z);
            }
    glEnd();

    //xz, up right
    glViewport(win_width/2, win_height/2, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_edgepoints[i][y_buffer_face[i][j]];k++)
            {
                glColor3f(EdgePoint[i][y_buffer_face[i][j]][k].r, EdgePoint[i][y_buffer_face[i][j]][k].g, EdgePoint[i][y_buffer_face[i][j]][k].b);
                glVertex2f(EdgePoint[i][y_buffer_face[i][j]][k].x, EdgePoint[i][y_buffer_face[i][j]][k].z);
            }
    glEnd();
 
 
*/

       //draw all points
    if(halftone_sign == 0)
    {
    //xy, up left
    glViewport(0, win_height/2, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_allpoints_xy[i][z_buffer_face[i][j]];k++)
            {
                glColor3f(AllPointXY[i][z_buffer_face[i][j]][k].r, AllPointXY[i][z_buffer_face[i][j]][k].g, AllPointXY[i][z_buffer_face[i][j]][k].b);
                glVertex2f(AllPointXY[i][z_buffer_face[i][j]][k].o1,AllPointXY[i][z_buffer_face[i][j]][k].o2);
            }
    glEnd();

    glViewport(0, 0, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_allpoints_yz[i][x_buffer_face[i][j]];k++)
            {
                glColor3f(AllPointYZ[i][x_buffer_face[i][j]][k].r, AllPointYZ[i][x_buffer_face[i][j]][k].g, AllPointYZ[i][x_buffer_face[i][j]][k].b);
                glVertex2f(AllPointYZ[i][x_buffer_face[i][j]][k].o1,AllPointYZ[i][x_buffer_face[i][j]][k].o2);
            }
    glEnd();

    glViewport(win_width/2, win_height/2, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<num_allpoints_xz[i][y_buffer_face[i][j]];k++)
            {
                glColor3f(AllPointXZ[i][y_buffer_face[i][j]][k].r, AllPointXZ[i][y_buffer_face[i][j]][k].g, AllPointXZ[i][y_buffer_face[i][j]][k].b);
                glVertex2f(AllPointXZ[i][y_buffer_face[i][j]][k].o1,AllPointXZ[i][y_buffer_face[i][j]][k].o2);
            }
    glEnd();
    }
 
    

    if((halftone_sign==1)||(halftone_sign==2))
    {
    //xy, up left
    glViewport(0, win_height/2, win_width/2, win_height/2);
    glBegin(GL_POINTS);
    for(int i=0;i<num_poly;i++)
        for(int j=0;j<num_faces[i];j++)
            for(int k=0;k<9*num_allpoints_xy[i][z_buffer_face[i][j]]+8;k++)
            {
                glColor3f(HalftongXY[i][z_buffer_face[i][j]][k].r, HalftongXY[i][z_buffer_face[i][j]][k].g, HalftongXY[i][z_buffer_face[i][j]][k].b);
                glVertex2f(HalftongXY[i][z_buffer_face[i][j]][k].o1,HalftongXY[i][z_buffer_face[i][j]][k].o2);
            }
    glEnd();
        
        glViewport(0, 0, win_width/2, win_height/2);
        glBegin(GL_POINTS);
        for(int i=0;i<num_poly;i++)
            for(int j=0;j<num_faces[i];j++)
                for(int k=0;k<9*num_allpoints_yz[i][x_buffer_face[i][j]]+8;k++)
                {
                    glColor3f(HalftongYZ[i][x_buffer_face[i][j]][k].r, HalftongYZ[i][x_buffer_face[i][j]][k].g, HalftongYZ[i][x_buffer_face[i][j]][k].b);
                    glVertex2f(HalftongYZ[i][x_buffer_face[i][j]][k].o1,HalftongYZ[i][x_buffer_face[i][j]][k].o2);
                }
        glEnd();

        glViewport(win_width/2, win_height/2, win_width/2, win_height/2);
        glBegin(GL_POINTS);
        for(int i=0;i<num_poly;i++)
            for(int j=0;j<num_faces[i];j++)
                for(int k=0;k<9*num_allpoints_xz[i][y_buffer_face[i][j]]+8;k++)
                {
                    glColor3f(HalftongXZ[i][y_buffer_face[i][j]][k].r, HalftongXZ[i][y_buffer_face[i][j]][k].g, HalftongXZ[i][y_buffer_face[i][j]][k].b);
                    glVertex2f(HalftongXZ[i][y_buffer_face[i][j]][k].o1,HalftongXZ[i][y_buffer_face[i][j]][k].o2);
                }
        glEnd();
    }
 
 
    
    
    glFlush();
    //blits the current opengl framebuffer on the screen
    glutSwapBuffers();
    //checks for opengl errors
}

void operates()
{
    cout<<"----------------------------------------------------------------"<<endl;
    cout<<"Yes! Use phong model! Please press 1. No? I want to keep this color. Please 0."<<endl;
    cin>>phong_sign;
    if(phong_sign==1)
    {
        phong_model_vertices();
        Get_SidePoints();
        Get_AllPoints();
    }
    cout<<"----------------------------------------------------------------"<<endl;
    cout<<"Yes! Use black/white halftoning please press 1. For RGB halftoning, please press 2, No? Please press 0."<<endl;
    cin>>halftone_sign;
    if(halftone_sign==1||halftone_sign==2)
    {
        halftoning();
    }
    cout<<"Here is the result"<<endl;
    cout<<"Please press 'o' to make more operations!"<<endl;
}

void key(unsigned char ch, int x, int y)
{
    if(ch=='o')
    {
        operates();
        glutPostRedisplay();
    }
}

void reshape(int reshape_width, int reshape_height)
{
    /*set up projection matrix to define the view port*/
    //update the ne window width and height
    win_width = reshape_width;
    win_height = reshape_height;
    
    //creates a rendering area across the window
    glViewport(0,0,reshape_width,reshape_height);
    // up an orthogonal projection matrix so that
    // the pixel space is mapped to the grid space
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,width,0,length,-10,10);
    
    //clear the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    //set pixel size based on width, if the aspect ratio
    //changes this hack won't work as well
    pixel_size = reshape_width/(float)width;
    
    //set pixel size relative to the grid cell size
    glPointSize(pixel_size);
    //check for opengl errors
}

void Init()
{
  /* Set clear color to white */
  glClearColor(0,0,0,0);
  /* Set fill color to black */
  /* glViewport(0 , 0 , 640 , 480); */
  glMatrixMode(GL_PROJECTION);
  /* glLoadIdentity(); */
  gluOrtho2D(0, win_width , 0 , win_height);
}

int main(int argc, char **argv)
{
    cout<<"Welcome to Project3!"<<endl;
    cout<<"HAVE FUN!"<<endl;
    ReadFile();
    MinMaxFaceBuffer();
    Get_SidePoints();
    Get_AllPoints();
    cout<<"Here is the orignal object."<<endl;
    cout<<"Please press 'o' to make more operations like phong model and halftone!"<<endl;
    //operates();

    /* Initialise GLUT library */
    glutInit(&argc,argv);
    /* Set the initial display mode */
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    /* Set the initial window position and size */
    glutInitWindowPosition(100,100);
    glutInitWindowSize(win_width,win_height);
    glutCreateWindow("Project3");
    /* Initialize drawing colors */
    Init();
    /* Call the displaying function */
    glutDisplayFunc(DrawPoly);
    glutReshapeFunc(reshape);
    /* Keep displaying untill the program is closed */
    glutKeyboardFunc(key);
    glutMainLoop();
    
    
    for(int i =0;i<num_poly;i++)
    {
        for(int j=0;j<num_face_init;j++)
        {
            delete [] EdgePoint[i][j];
            delete [] AllPointXY[i][j];
            delete [] AllPointYZ[i][j];
            delete [] AllPointXZ[i][j];
            delete [] HalftongXY[i][j];
            delete [] HalftongXY[i][j];
            delete [] HalftongXY[i][j];
        }
        delete [] num_edgepoints[i];
        delete [] Vertices[i];
        delete [] Faces[i];
        delete [] EdgePoint[i];
        delete [] AllPointXY[i];
        delete [] AllPointYZ[i];
        delete [] AllPointXZ[i];
        delete [] HalftongXY[i];
        delete [] HalftongYZ[i];
        delete [] HalftongXZ[i];
        delete [] phone_vertices[i];
        delete [] faces_per_vertice[i];
        delete [] specu_per_vertice[i];
        delete [] normal_vector[i];
        delete [] normal_vector_point[i];
        delete [] num_allpoints_xy[i];
        delete [] num_allpoints_yz[i];
        delete [] num_allpoints_xz[i];
        delete [] face_min[i];
        delete [] face_max[i];
    }
    delete [] num_vertices;
    delete [] num_faces;
    delete [] num_edgepoints;
    delete [] Vertices;
    delete [] Faces;
    delete [] EdgePoint;
    delete [] AllPointXY;
    delete [] AllPointYZ;
    delete [] AllPointXZ;
    delete [] HalftongXY;
    delete [] HalftongYZ;
    delete [] HalftongXZ;
    delete [] phone_vertices;
    delete [] faces_per_vertice;
    delete [] specu_per_vertice;
    delete [] normal_vector;
    delete [] normal_vector_point;
    delete [] num_allpoints_xy;
    delete [] num_allpoints_yz;
    delete [] num_allpoints_xz;
    delete [] face_min;
    delete [] face_max;
    
}
