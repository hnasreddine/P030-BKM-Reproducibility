#include <quadmath.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
struct DD {double h,l; DD(double x=0):h(x),l(0){} DD(double x,double y):h(x),l(y){} };
inline DD operator+(DD a,DD b){double s=a.h+b.h;double v=s-a.h;double t=((b.h-v)+(a.h-(s-v)))+a.l+b.l;double z=s+t;return {z,t-(z-s)};}
inline DD operator-(DD a,DD b){double s=a.h-b.h;double v=s-a.h;double t=((-b.h-v)+(a.h-(s-v)))+a.l-b.l;double z=s+t;return {z,t-(z-s)};}
inline DD operator*(DD a,DD b){double p=a.h*b.h;double e=std::fma(a.h,b.h,-p)+a.h*b.l+a.l*b.h;double z=p+e;return {z,e-(z-p)};}
inline DD operator*(double a,DD b){return DD(a)*b;}
inline DD operator*(DD a,double b){return a*DD(b);}
inline DD ddparse(const char*s){__float128 q=strtoflt128(s,nullptr);double h=(double)q;return {h,(double)(q-(__float128)h)};}

#include <cstring>
int main(int argc,char**argv){
 FILE*f=fopen(argv[1],"r");
 char R[256], T[256],lam[256]; int N,K;
 fscanf(f,"%d %d %255s %255s %255s",&N,&K,R,T,lam);
 auto emit=[](const char* v){DD a=ddparse(v);printf("%s %a %a\n",v,a.h,a.l);};
 emit(lam);
 char token[256];
 for (int j=0;j<4*N+K+1;j++){if(fscanf(f,"%255s",token)!=1)return 8;emit(token);}
 fclose(f);
}
