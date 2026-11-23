#include <quadmath.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <cfenv>
#include <cstdlib>
struct DD {double h,l; DD(double x=0):h(x),l(0){} DD(double x,double y):h(x),l(y){} };
inline DD operator+(DD a,DD b){double s=a.h+b.h;double v=s-a.h;double t=((b.h-v)+(a.h-(s-v)))+a.l+b.l;double z=s+t;return {z,t-(z-s)};}
inline DD operator-(DD a,DD b){double s=a.h-b.h;double v=s-a.h;double t=((-b.h-v)+(a.h-(s-v)))+a.l-b.l;double z=s+t;return {z,t-(z-s)};}
inline DD operator*(DD a,DD b){double p=a.h*b.h;double e=std::fma(a.h,b.h,-p)+a.h*b.l+a.l*b.h;double z=p+e;return {z,e-(z-p)};}
inline DD operator*(double a,DD b){return DD(a)*b;}
inline DD operator*(DD a,double b){return a*DD(b);}
inline DD ddparse(const char*s){__float128 q=strtoflt128(s,nullptr);double h=(double)q;return {h,(double)(q-(__float128)h)};}
// Local rounding ball radius for DD operators.
// Error-bound lemma under the stated IEEE754 assumptions, see PROOF_DD_AND_NUMERICAL_CERTIFICATE_2026-11-22.md: for normalized DD with
// |low|<=2^-53 |high|, absolute addition error <=100u²(|a|+|b|)+eta,
// multiplication error <=100u²|a||b|+eta, u=2^-53, eta=1e-300.
// This instrument checks the execution; its local roundoff lemma is proved separately.
static constexpr double U2=0x1.p-106;
static constexpr double CDD=100.0;
static constexpr double ETA=1e-300;
static long long bad_nonoverlap=0;
static double max_overlap_ratio=0;
static long long checked_dd=0;
static long long nonfinite=0;
static long long zero_high_nonzero_low=0;
inline double mg(DD v){return fabs((double)v.h)+fabs((double)v.l);}
inline void checkoverlap(DD v){
 checked_dd++;
 if (!std::isfinite(v.h)||!std::isfinite(v.l)) {nonfinite++; if(nonfinite<5) fprintf(stderr,"NONFINITE h %.17g l %.17g\n",v.h,v.l); return;}
 const double ah=std::fabs(v.h), al=std::fabs(v.l);
 if(ah==0){if(al!=0) {zero_high_nonzero_low++;if(zero_high_nonzero_low<5) fprintf(stderr,"ZERO_HIGH_WITH_LOW %.17g\n",v.l);}return;}
 const double ratio=al/ah;
 if(ratio>max_overlap_ratio)max_overlap_ratio=ratio;
 if(al>ah*0x1.0001p-53) {bad_nonoverlap++;if(bad_nonoverlap<5) fprintf(stderr,"BAD_DD %.17g %.17g ratio %.18g\n",v.h,v.l,ratio);}
}
struct E {DD v;double e;E(DD x=DD(),double y=0):v(x),e(y) {}};
inline E operator+(E a,E b){DD z=a.v+b.v;checkoverlap(z);return {z,a.e+b.e+CDD*U2*(mg(a.v)+mg(b.v))+ETA};}
inline E operator-(E a,E b){DD z=a.v-b.v;checkoverlap(z);return {z,a.e+b.e+CDD*U2*(mg(a.v)+mg(b.v))+ETA};}
inline E operator*(E a,E b){DD z=a.v*b.v;checkoverlap(z);return {z,mg(a.v)*b.e+mg(b.v)*a.e+a.e*b.e+CDD*U2*mg(a.v)*mg(b.v)+ETA};}

int main(int argc,char**argv){
if (!std::numeric_limits<double>::is_iec559 || std::fegetround()!=FE_TONEAREST) {fprintf(stderr,"IEEE_TONEAREST_CHECK_FAILED\n"); return 6;}
fprintf(stderr,"IEEE754_ROUND_NEAREST_PASS\n");
FILE*f=fopen(argv[1],"r");int N,K;double R,T;char lambda_token[256];fscanf(f,"%d %d %lf %lf %255s",&N,&K,&R,&T,lambda_token); DD lam=ddparse(lambda_token);
 std::vector<DD> x(N),a(N),b(N),d(N),J(K+1),ya(2*N),yb(2*N),yc(2*N);char str[256];
 for(int i=0;i<N;i++){for(DD*p:{&x[i],&a[i],&b[i],&d[i]}){fscanf(f,"%255s",str);*p=ddparse(str);}}
 for(int k=0;k<=K;k++){fscanf(f,"%255s",str);J[k]=ddparse(str);}fclose(f);
checkoverlap(lam);for (int i=0;i<N;i++){checkoverlap(x[i]);checkoverlap(a[i]);checkoverlap(b[i]);checkoverlap(d[i]);}for(int k=0;k<=K;k++)checkoverlap(J[k]);
 DD lamR=lam*(1/R),twoR=DD(2.0/R),twolamR=DD(2.0)*lamR;
 FILE*out=fopen(argv[2],"w");
 double sumR[2]={0,0},initR[2]={0,0},sumOut[2]={0,0},maxR[2]={0,0};
 for(int col=0;col<2;col++){
  std::fill(ya.begin(),ya.end(),DD());std::fill(yb.begin(),yb.end(),DD());
  DD pa0(col==0),pa1(col==1),pb0,pb1;
  E re0=E(J[0])*E(pa0),re1=E(J[0])*E(pa1),im0,im1;
  double initnormsq=0;
  for(int i=0;i<N;i++) {
   E z0=E(lamR)*(E(a[i])*E(pa0)+E(b[i])*E(pa1));
   E z1=E(lamR)*(E(b[i])*E(pa0)+E(d[i])*E(pa1));
   yb[2*i]=z0.v;yb[2*i+1]=z1.v;
   initnormsq+=z0.e*z0.e+z1.e*z1.e;
  }
  initR[col]=sqrt(initnormsq);
  for(int k=2;k<=K;k++){
   E tmp0,tmp1;
   for(int i=0;i<N;i++){
   tmp0=tmp0+E(a[i])*E(yb[2*i])+E(b[i])*E(yb[2*i+1]);
   tmp1=tmp1+E(b[i])*E(yb[2*i])+E(d[i])*E(yb[2*i+1]);
  }
   E zp0=E(twolamR)*tmp0-E(pa0),zp1=E(twolamR)*tmp1-E(pa1);
   DD pk0=zp0.v,pk1=zp1.v;
   double res2=zp0.e*zp0.e+zp1.e*zp1.e;
   for(int i=0;i<N;i++){
    E z0=E(twoR)*(E(x[i])*E(yb[2*i])+E(lam)*(E(a[i])*E(pb0)+E(b[i])*E(pb1)))-E(ya[2*i]);
    E z1=E(twoR)*(E(x[i])*E(yb[2*i+1])+E(lam)*(E(b[i])*E(pb0)+E(d[i])*E(pb1)))-E(ya[2*i+1]);
    yc[2*i]=z0.v;yc[2*i+1]=z1.v;
    res2+=z0.e*z0.e+z1.e*z1.e;
   }
   E cf=E(DD(2))*E(J[k]);
   double thisR=sqrt(res2);sumR[col]+=thisR;maxR[col]=std::max(maxR[col],thisR);
   if(k%4000==0)fprintf(stderr,"col %d k %d sumR %.4e curR %.4e\n",col,k,sumR[col],thisR);
   switch(k%4){case 0:re0=re0+cf*E(pk0);re1=re1+cf*E(pk1);break;
    case 1:im0=im0-cf*E(pk0);im1=im1-cf*E(pk1);break;
    case 2:re0=re0-cf*E(pk0);re1=re1-cf*E(pk1);break;
    case 3:im0=im0+cf*E(pk0);im1=im1+cf*E(pk1);break;}
   pa0=pb0;pa1=pb1;pb0=pk0;pb1=pk1;ya.swap(yb);yb.swap(yc);
  }
  fprintf(out,"%.17g %.17g %.17g %.17g\n%.17g %.17g %.17g %.17g\n",re0.v.h,re0.v.l,im0.v.h,im0.v.l,re1.v.h,re1.v.l,im1.v.h,im1.v.l);
  sumOut[col]=re0.e+im0.e+re1.e+im1.e;
  fprintf(stderr,"COL %d SUM_RES %.19e INIT_RES %.19e OUTPUT_SUM_ERR %.19e\n",col,sumR[col],initR[col],sumOut[col]);fflush(out);
  fprintf(stderr,"dd col %d finished\n",col);
 }
 fclose(out);
 fprintf(stderr,"GLOBAL_SUMRES %.19e INITRES %.19e OUTSUM %.19e MAXRES %.19e BAD_NONOVERLAP %lld MAX_OVERLAP %.19e\n",sumR[0]+sumR[1],initR[0]+initR[1],sumOut[0]+sumOut[1],std::max(maxR[0],maxR[1]),bad_nonoverlap,max_overlap_ratio);
 fprintf(stderr,"NONFINITE %lld ZERO_HIGH_NONZERO_LOW %lld CHECKED_DD %lld\n",nonfinite,zero_high_nonzero_low,checked_dd);
}
