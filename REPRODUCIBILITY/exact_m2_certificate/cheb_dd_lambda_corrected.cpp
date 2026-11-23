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
int main(int argc,char**argv){FILE*f=fopen(argv[1],"r");int N,K;double R,T;char lambda_token[256];fscanf(f,"%d %d %lf %lf %255s",&N,&K,&R,&T,lambda_token); DD lam=ddparse(lambda_token);
 std::vector<DD> x(N),a(N),b(N),d(N),J(K+1),ya(2*N),yb(2*N),yc(2*N);char str[256];
 for(int i=0;i<N;i++){for(DD*p:{&x[i],&a[i],&b[i],&d[i]}){fscanf(f,"%255s",str);*p=ddparse(str);}}
 for(int k=0;k<=K;k++){fscanf(f,"%255s",str);J[k]=ddparse(str);}fclose(f);
 DD lamR=lam*(1/R),twoR=DD(2.0/R),twolamR=DD(2.0)*lamR;
 FILE*out=fopen(argv[2],"w");
 for(int col=0;col<2;col++){
  std::fill(ya.begin(),ya.end(),DD());std::fill(yb.begin(),yb.end(),DD());
  DD pa0(col==0),pa1(col==1),pb0,pb1;
  DD re0=J[0]*pa0,re1=J[0]*pa1,im0,im1;
  for(int i=0;i<N;i++) {yb[2*i]=lamR*(a[i]*pa0+b[i]*pa1);yb[2*i+1]=lamR*(b[i]*pa0+d[i]*pa1);}
  for(int k=2;k<=K;k++){
   DD tmp0,tmp1;
   for(int i=0;i<N;i++){tmp0=tmp0+a[i]*yb[2*i]+b[i]*yb[2*i+1];tmp1=tmp1+b[i]*yb[2*i]+d[i]*yb[2*i+1];}
   DD pk0=twolamR*tmp0-pa0,pk1=twolamR*tmp1-pa1;
   for(int i=0;i<N;i++){
    yc[2*i]=twoR*(x[i]*yb[2*i]+lam*(a[i]*pb0+b[i]*pb1))-ya[2*i];
    yc[2*i+1]=twoR*(x[i]*yb[2*i+1]+lam*(b[i]*pb0+d[i]*pb1))-ya[2*i+1];
   }
   DD cf=DD(2)*J[k];
   switch(k%4){case 0:re0=re0+cf*pk0;re1=re1+cf*pk1;break;
    case 1:im0=im0-cf*pk0;im1=im1-cf*pk1;break;
    case 2:re0=re0-cf*pk0;re1=re1-cf*pk1;break;
    case 3:im0=im0+cf*pk0;im1=im1+cf*pk1;break;}
   pa0=pb0;pa1=pb1;pb0=pk0;pb1=pk1;ya.swap(yb);yb.swap(yc);
  }
  fprintf(out,"%.17g %.17g %.17g %.17g\n%.17g %.17g %.17g %.17g\n",re0.h,re0.l,im0.h,im0.l,re1.h,re1.l,im1.h,im1.l);fflush(out);
  fprintf(stderr,"dd col %d finished\n",col);
 }
 fclose(out);
}
