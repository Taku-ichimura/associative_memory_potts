#include <iostream>
#include <string>
#include <random>
#include<fstream>
#include <cmath>
#include<omp.h>

using namespace std;

random_device rnd;//for nondeterministic random_num
//const int seed=rnd();
const int seed = 11;
mt19937 mt(seed);

//parameters
const int X=6;//width
const int Y=13;//height
const double J=0.0;//exchange interaction average
const double J_d = 1.0;//exchange disparsion
const int Ncol=5;//number of color+null
const int Nconne=4;//vanish  with  num of connect
const int Nmcstep=300;//num of MCMC step
const int Nms=2000;//num of MCMCstep with taking data
const int Ntrain=4000;//num of train
const int NT=10;
const double T_min=0.05;//temperature
const double T_max=1.05;
const double eta=0.005;
//

int field[X][Y];//playfield
int fldflag[X][Y];//for field serch
int depth[X];//height of col
int Nchain;//number of chain
int cou=0;
double T;
bool chainflg;

class BM{
  int Ns=1;
  int tumo[2];
  int whput[2][2];
  double dqs_ave[X*Y][Ncol]={};
  double mqs_ave[X*Y][Ncol]={};
  double dss_ave[(X*Y)*(X*Y-1)/2]={};
  double mss_ave[(X*Y)*(X*Y-1)/2]={};
  double Jij[(X*Y)*(X*Y-1)/2];
  double Hqi[X*Y][Ncol]={};
  string URL[10000];
  bool goFlag;

public:
  BM(){};
  void ini();
  void readData();
  void decode();
  void decoder(int Ns_i,int f[X][Y]);
  void dataave();
  void MCMC(double T);
  void modelave();
  void pupdate();
  void rndfill();
  void printfld();
  void outputips();
  void debug();
  void prn_d();
  void prn_m();
  void paramwrite();
  void readparam();
};

void BM::ini(){
  double tmp;
  double rl;
  double a,b;

  for(int i=0;i<X*Y;i++){
    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    //Hqi[i][0] = a;
    //Hqi[i][1] = b;


    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    //Hqi[i][2] = a;
    //Hqi[i][3] = b;


    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    //Hqi[i][4] = a;

  }

  for(int i=0;i<(X*Y)*(X*Y-1)/2;i++){
    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    Jij[i] = a;
  }

}

void BM::readData(){
  int i=0;
  ifstream ifs;
  ifs.open("train.d");
  while(!ifs.eof()){
    ifs>>URL[i];
    URL[i] = URL[i].erase(0,33);
    i++;
  }
  Ns = i;
}

void BM::decode(){
  string url="http://ips.karou.jp/simu/pe.html?";
  string urllist="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  string tmp;
  int f[X][Y];
  int a;
  int s=0;
  for(int i=0;i<Ns;i++){
    for(int j=Y-1;j>=0;j--){
      for(int k=0;k<X;k+=2){
        a = urllist.find(URL[i][s]);
        f[k+1][j] = a%8;
        f[k][j] = (a-a%8)/8;
        s++;
      }
    }
  }

}

void BM::decoder(int Ns_i,int f[X][Y]){
  string url="http://ips.karou.jp/simu/pe.html?";
  string urllist="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  string tmp;
  int a;
  int s=0;

  for(int j=Y-1;j>=0;j--){
    for(int k=0;k<X;k+=2){
      a = urllist.find(URL[Ns_i][s]);
      f[k+1][j] = a%8;
      f[k][j] = (a-a%8)/8;
      s++;
    }
  }

}

void BM::dataave(){
  int f[X][Y];
  int f2[X][Y];
  int f3[X][Y];
  int f4[X][Y];
  int s=0;


  for(int Ns_i=0;Ns_i<Ns;Ns_i++){
    decoder(Ns_i,f);
    //decoder(Ns_i,f2);
    //decoder(Ns_i,f3);
    //decoder(Ns_i,f4);

    // for(int i=0;i<X;i++){
    //   for(int j=0;j<Y;j++){
    //     if(f[i][j]!=0){
    //       f2[i][j]=(f2[i][j]-1+1)%(Ncol-1)+1;
    //       f3[i][j]=(f3[i][j]-1+2)%(Ncol-1)+1;
    //       f4[i][j]=(f4[i][j]-1+3)%(Ncol-1)+1;
    //     }
    //   }
    // }

    s=0;
    for(int i=0;i<X;i++){
       for(int j=0;j<Y;j++){
        if(f[i][j]==0){
          dqs_ave[i+X*j][0] +=1.0/Ns;
        }
        else{
          dqs_ave[i+X*j][1] +=1.0/Ns;
        }
          // dqs_ave[i+X*j][f[i][j]] +=1.0/Ns/4;
          //  dqs_ave[i+X*j][f2[i][j]] +=1.0/Ns/4;
          //  dqs_ave[i+X*j][f3[i][j]] +=1.0/Ns/4;
          //  dqs_ave[i+X*j][f4[i][j]] +=1.0/Ns/4;

          s=i+X*j;
          for(int l=0;l<s+1;l++){
            if(i+X*j!=l){
              if(f[i][j]==f[l%X][int(l/X)])dss_ave[s*(s-1)/2+l]+=1.0/Ns;
              //if(f2[i][j]==f2[l%X][int(l/X)])dss_ave[s*(s-1)/2+l]+=1.0/Ns/4;
              //if(f3[i][j]==f3[l%X][int(l/X)])dss_ave[s*(s-1)/2+l]+=1.0/Ns/4;
              //if(f4[i][j]==f4[l%X][int(l/X)])dss_ave[s*(s-1)/2+l]+=1.0/Ns/4;
            }
          }

        }

      }

  }

}

void BM::MCMC(double T){
  int s=0;
  int tmp,tmp2;
  double local;
  double dE,oldE,newE;
  double p;

  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
      local =0.0;

      s=i+X*j;
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(int l=0;l<X*Y;l++){

        if(s<l){
          if(field[i][j]==field[l%X][int(l/X)])local+=Jij[l*(l-1)/2+s];
        }
        else{
          if(s!=l && field[i][j]==field[l%X][int(l/X)])local+=Jij[s*(s-1)/2+l];
        }

      }

      oldE= 1.0*local+ Hqi[i+j*X][field[i][j]];

      tmp = (mt()%(Ncol));
      while(tmp==field[i][j])tmp = (mt()%(Ncol));//exclude now state
      tmp2 = field[i][j];
      field[i][j]=tmp;
      local=0;

      s=i+X*j;
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(int l=0;l<X*Y;l++){
        if(s<l){
          if(s!=l && field[i][j]==field[l%X][int(l/X)])local+=Jij[l*(l-1)/2+s];
        }
        else{
          if(s!=l && field[i][j]==field[l%X][int(l/X)])local+=Jij[s*(s-1)/2+l];
        }

      }

      newE= 1.0*local + Hqi[i+j*X][field[i][j]];

      p=1.0*mt()/mt.max();
      if(exp(-(newE-oldE)/T)>p){
      }
      else{
        field[i][j]=tmp2;//return state
      }


    }
  }
}


void BM::modelave(){
  int f[X][Y];
  int s=0;

  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){

      if(field[i][j]==0){
        mqs_ave[i+X*j][0] +=1.0/Nms;
      }
      else{
        mqs_ave[i+X*j][1] +=1.0/Nms;
      }
        //mqs_ave[i+X*j][field[i][j]] +=1.0/Nms;
        s=i+X*j;
        for(int l=0;l<s+1;l++){
          if(i+X*j!=l){
            if(field[i][j]==field[l%X][int(l/X)])mss_ave[s*(s-1)/2+l]+=1.0/Nms;
          }

        }
      }
    }

}


void BM::pupdate(){
  int s=0;
  string fname="Jijstep.d";
  ofstream outparamJij(fname,ios::app);
  string fname2="Hqistep.d";
  ofstream outparamHqi(fname2,ios::app);

  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
        Hqi[i+j*X][0]=Hqi[i+j*X][0]-eta*(dqs_ave[i+X*j][0]-mqs_ave[i+X*j][0]);
        Hqi[i+j*X][1]=Hqi[i+j*X][1]-eta*(dqs_ave[i+X*j][1]-mqs_ave[i+X*j][1]);
        // Hqi[i+j*X][2]=Hqi[i+j*X][2]-eta*(dqs_ave[i+X*j][2]-mqs_ave[i+X*j][2]);
        // Hqi[i+j*X][3]=Hqi[i+j*X][3]-eta*(dqs_ave[i+X*j][3]-mqs_ave[i+X*j][3]);
        // Hqi[i+j*X][4]=Hqi[i+j*X][4]-eta*(dqs_ave[i+X*j][4]-mqs_ave[i+X*j][4]);

        Hqi[i+j*X][2]=Hqi[i+j*X][1];
        Hqi[i+j*X][3]=Hqi[i+j*X][1];
        Hqi[i+j*X][4]=Hqi[i+j*X][1];

        mqs_ave[i+X*j][0]=0;
        mqs_ave[i+X*j][1]=0;
        mqs_ave[i+X*j][2]=0;
        mqs_ave[i+X*j][3]=0;
        mqs_ave[i+X*j][4]=0;
        s=i+X*j;
        for(int l=0;l<s+1;l++){
          if(i+X*j!=l){
            Jij[s*(s-1)/2+l]=Jij[s*(s-1)/2+l]-eta*(dss_ave[s*(s-1)/2+l]-mss_ave[s*(s-1)/2+l]);
            mss_ave[s*(s-1)/2+l]=0;
            //cout<<dss_ave[s*(s-1)/2+l]-mss_ave[s*(s-1)/2+l]<<endl;
          }

        }
      }
    }
    outparamJij<<cou<<" "<<Jij[0]<<" "<<Jij[1]<<" "<<Jij[2]<<endl;
    outparamHqi<<cou<<" "<<Hqi[0][0]<<" "<<Hqi[73][0]<<endl;
    cou++;
}


void BM::rndfill(){
  for(int i=0;i<X;i++)for(int j=0;j<Y;j++) field[i][j]=mt()%(Ncol);
}

void BM::printfld(){
  for(int i=Y-1;i>=0;i--){
    for(int j=0;j<X;j++){
      cout<<field[j][i];
    }
    cout<<""<<endl;
  }
  cout<<""<<endl;
}

void BM::outputips(){
  string url="http://ips.karou.jp/simu/pe.html?";
  string urllist="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  string tmp;
  for(int i=Y-1;i>=0;i--){
    for(int j=0;j<X;j+=2){
      tmp = urllist[field[j][i]*8+field[j+1][i]];
      url+=tmp;
    }
  }
  cout<<url<<endl;
}

void BM::debug(){
  double tmp;
  double rl;
  double a,b;

  for(int i=0;i<X*Y;i++){
    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    Hqi[i][0] = 0;
    Hqi[i][1] = 0;

    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    Hqi[i][2] = 0;
    Hqi[i][3] = 0;

    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    Hqi[i][4] = 0;

  }

  for(int i=0;i<(X*Y)*(X*Y-1)/2;i++){
    rl = sqrt(-2.0*log(1.0*mt()/mt.max()));
    tmp = M_PI*2*(1.0*mt()/mt.max());
    a = rl*cos(tmp);
    b = rl*sin(tmp);
    Jij[i] = -0.1;
  }

}

void BM::prn_d(){
  int s=0;
  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
        cout<<dqs_ave[i+X*j][1]<<endl;
        s=i+X*j;
        for(int l=0;l<s+1;l++){
          if(i+X*j!=l){
            cout<<dss_ave[s*(s-1)/2+l]<<s*(s-1)/2+l<<endl;
          }
        }

      }
    }

}

void BM::prn_m(){
  int s=0;
  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
        cout<<mqs_ave[i+X*j][1]<<endl;
        s=i+X*j;
        for(int l=0;l<s+1;l++){
          if(i+X*j!=l){
            cout<<mss_ave[s*(s-1)/2+l]<<s*(s-1)/2+l<<endl;
          }
        }

      }
    }

}

void BM::paramwrite(){
  string fname="Hqiparam.d";
  ofstream outparamHqi(fname);
  string fname2="Jijparam.d";
  ofstream outparamJij(fname2);
  int s=0;

  for(int i=0;i<X*Y;i++){
    outparamHqi<<i<<" "<<Hqi[i][0]<<" "<<Hqi[i][1]<<" "<<Hqi[i][2]<<" "<<Hqi[i][3]<<" "<<Hqi[i][4]<<endl;
  }

  // for(int i=0;i<(X*Y)*(X*Y-1)/2;i++){
  //   outparamJij<<i<<Jij[i]<<endl;
  // }

  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
        s=i+X*j;
        for(int l=0;l<s+1;l++){
          if(i+X*j!=l){
            outparamJij<<i<<" "<<j<<" "<<l%X<<" "<<l/X<<" "<<Jij[s*(s-1)/2+l]<<endl;
          }
        }

}
}
}

void BM::readparam(){
  ifstream paraJij;
  paraJij.open("Jijparam.d");
  ifstream paraHqi;
  paraHqi.open("Hqiparam.d");
  int s=0;
  int a,b,c,d;

  for(int i=0;i<X*Y;i++){
    paraJij>>a>>Hqi[i][0]>>Hqi[i][1]>>Hqi[i][2]>>Hqi[i][3];
  }

  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
        s=i+X*j;
        for(int l=0;l<s+1;l++){
          if(i+X*j!=l){
            paraJij>>a>>b>>c>>d>>Jij[s*(s-1)/2+l];
          }
        }
      }
    }

}

int main(){
  BM one ;
  double dT=(T_max-T_min)/NT;

  one.ini();
  one.readData();
  one.dataave();
  one.rndfill();
  //one.debug();

  for(int N_i=0;N_i<Ntrain;N_i++){
    T=T_max;

    for(int j=0;j<NT;j++){
      T=T_max-(j+1)*dT;
      for(int i=0;i<Nmcstep;i++){
        one.MCMC(T);
      }
    }

    for(int i=0;i<Nms;i++){
      one.MCMC(T);
      one.modelave();
    }

    one.pupdate();

  }
  one.paramwrite();

  one.rndfill();
  for(int j=0;j<NT;j++){
    T=T_max-(j+1)*dT;
    for(int i=0;i<Nmcstep*2;i++){
      one.MCMC(T);
    }
  }
  one.printfld();
  one.outputips();


  return 0;
}
