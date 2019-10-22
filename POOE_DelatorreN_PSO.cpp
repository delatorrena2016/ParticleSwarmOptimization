#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <thread>

using namespace std;

class           Particle{
private:
  vector<float> m_v, m_x, m_x_b;
  float         m_fitness, m_pbest;

public:

  int           m_id = 0;

  Particle() {}
  ~Particle() {}
  Particle(vector<float> rc) {
    m_v = m_x = m_x_b = rc;
    static int id = 1; m_id = id++;
  }

  bool          operator > (const Particle& par)
                const {return m_fitness > par.m_fitness;}

  vector<float> const &getV() const { return m_v; }
  vector<float>  &getX()  { return m_x; }
  vector<float> const &getXB() const { return m_x_b; }
  float         getFit() { return m_fitness; }
  float         getPbest() { return m_pbest; }

  void          setFit( float fit) { m_fitness = fit; }
  void          setPbest( float fit) { m_pbest = fit; }
  void          setX(const vector<float>& _x) { m_x = _x; }
  void          setV(const vector<float>& _x) { m_v = _x; }
  void          setXB(const vector<float>& _x) { m_x_b = _x; }
  void          log();
  void          logC();
};
void            SwarmAlgorithm(float (*function)(vector<float>& x),float& max,float& min,unsigned int& count,vector<Particle>& x_b_g,float& phi_1,
                float& phi_2,vector<float>& rand_coods,vector<Particle>& pop,unsigned int& M_pop,unsigned int& dimension,unsigned int N);
void            simMov(float& max,float& min,unsigned int& count,vector<Particle>& x_b_g,float& phi_1,float& phi_2,vector<float>& rand_coods,
                vector<Particle>& pop,unsigned int& M_pop,unsigned int& dimension);
vector<float>   randGenMS(unsigned int& dimension, float& max, float& min,vector<float>& rand_coods);
float           sphere(vector<float>& x);
float           Beale(vector<float>& x);
float           GoldPrice(vector<float>& x);
float           McCormick(vector<float>& x);
void            reevFit(float (*function)(vector<float>& x),unsigned int& M_pop,vector<Particle>& pop);
vector<float>   suma(const vector<float>& u,const vector<float>& v);
vector<float>   resta(const vector<float>& u,const vector<float>& v);
vector<float>   product(const vector<float>& u,const vector<float>& v);
vector<float>   inWeight(2,0.7);
void            sortBuble(unsigned int& M_pop,vector<Particle>& pop);
void            initSwarm(float (*function)(vector<float>& x),unsigned int& M_pop,unsigned int& dimension,
                float& min,float& max,vector<float>& rand_coods,vector<Particle>& pop);
float           zero = 0.0;
ofstream        fout;
unsigned int    c = thread::hardware_concurrency();

int             main() {
  float            min = -2, max = 2;
  unsigned int     dimension = 2, M_pop = 15, N=50, count=0;
  float            phi_1 = 2.0, phi_2 = 2.0;
  vector<Particle> pop, x_b_g;
  vector<float>    rand_coods;

  cout<<"Probable no. of cores: "<<c<<endl;//A system can support multiple threads per procesing unit
  float            (* function)(vector<float>& x) = McCormick;
  thread f1(SwarmAlgorithm,(*function),ref(max),ref(min),ref(count),ref(x_b_g),ref(phi_1),ref(phi_2),ref(rand_coods),ref(pop),ref(M_pop),ref(dimension),N);
                   function = sphere;
  thread f2(SwarmAlgorithm,(*function),ref(max),ref(min),ref(count),ref(x_b_g),ref(phi_1),ref(phi_2),ref(rand_coods),ref(pop),ref(M_pop),ref(dimension),N);
                   function = Beale;
  thread f3(SwarmAlgorithm,(*function),ref(max),ref(min),ref(count),ref(x_b_g),ref(phi_1),ref(phi_2),ref(rand_coods),ref(pop),ref(M_pop),ref(dimension),N);
                   function = GoldPrice;
  thread f4(SwarmAlgorithm,(*function),ref(max),ref(min),ref(count),ref(x_b_g),ref(phi_1),ref(phi_2),ref(rand_coods),ref(pop),ref(M_pop),ref(dimension),N);
  f1.join(); f2.join(); f3.join(); f4.join();
  //SwarmAlgorithm((*function),max,min,count,x_b_g,phi_1,phi_2,rand_coods,pop,M_pop,dimension,N);
  return 0;
}
void SwarmAlgorithm(float (*function)(vector<float>& x),float& max,float& min,unsigned int& count,vector<Particle>& x_b_g,float& phi_1,float& phi_2,vector<float>& rand_coods,vector<Particle>& pop,unsigned int& M_pop,unsigned int& dimension,unsigned int N) {
  initSwarm((*function),M_pop,dimension,min,max,rand_coods,pop);
  while (count < N) {
    sortBuble(M_pop,pop);
    for (size_t i = 0; i < M_pop; i++) { pop.at(i).log(); }
    simMov(max,min,count,x_b_g,phi_1,phi_2,rand_coods,pop,M_pop,dimension);
    reevFit((*function),M_pop,pop);
    x_b_g.at(count).logC();
    count++;
    fout.open("Convergence.csv",fstream::app); fout << count << endl; fout.close();
  }
}
void            simMov(float& max,float& min,unsigned int& count,vector<Particle>& x_b_g,float& phi_1,float& phi_2,vector<float>& rand_coods,vector<Particle>& pop,unsigned int& M_pop,unsigned int& dimension) {
  x_b_g.push_back(pop.at(0)); vector<float> p1, p2;
  for (size_t i = 0; i < M_pop; i++) {
    p1 = product(resta(pop.at(i).getXB(),pop.at(i).getX()),randGenMS(dimension,phi_1,zero,rand_coods));
    /*for (size_t i = 0; i <2; i++) {cout << rand_coods.at(i)<<endl;}*/ vector<float>().swap(rand_coods);
    p2 = product(resta(x_b_g.at(count).getX(),pop.at(i).getX()),randGenMS(dimension,phi_2,zero,rand_coods));
    /*for (size_t i = 0; i <2; i++) {cout << rand_coods.at(i)<<endl;}*/ vector<float>().swap(rand_coods);
    pop.at(i).setV(suma(suma(product(inWeight,pop.at(i).getV()),p1),p2));
    pop.at(i).setX(suma(pop.at(i).getX(),pop.at(i).getV()));
    //Position adjustments
    for (unsigned int j = 0; j < dimension; j++) {
      while (pop.at(i).getX().at(j) < min || pop.at(i).getX().at(j) > max){
        if (pop.at(i).getX().at(j) > max) { pop.at(i).getX().at(j) = (max - abs(pop.at(i).getX().at(j) - max)); }
        else if (pop.at(i).getX().at(j) < min) { pop.at(i).getX().at(j) = (min + abs(pop.at(i).getX().at(j) - (min))); }
      }
    }
  }
}
void            reevFit(float (*function)(vector<float>& x),unsigned int& M_pop,vector<Particle>& pop){
  for (size_t i = 0; i < M_pop; i++) {
    pop.at(i).setFit((*function)(pop.at(i).getX()));
    if (pop.at(i).getFit() < pop.at(i).getPbest()) {
      pop.at(i).setPbest(pop.at(i).getFit()); pop.at(i).setXB(pop.at(i).getX());
    }
  }
}
void            initSwarm(float (*function)(vector<float>& x),unsigned int& M_pop,unsigned int& dimension,float& min,float& max,vector<float>& rand_coods,vector<Particle>& pop){
  for (size_t i = 0; i < M_pop; i++) {
    randGenMS(dimension,max,min,rand_coods);
    pop.emplace_back(rand_coods); vector<float>().swap(rand_coods); //de-allocates vector's memory
    pop.at(i).setFit((*function)(pop.at(i).getX())); pop.at(i).setPbest(pop.at(i).getFit());
  }
  fout.open("Log.csv",fstream::app);
  fout << "Id" << "," << "Fitness" << "," << "Pbest" << ","
       << "X"  << "," << "Y"       << "," << "V1"    << ","
       << "V2" << "," << "XB"      << "," << "YB"    << endl;
  fout.close();
  fout.open("Convergence.csv",fstream::app);
  fout << "Id" << "," << "Fitness" << "," << "Pbest" << ","
       << "X"  << "," << "Y"       << "," << "V1"    << ","
       << "V2" << "," << "XB"      << "," << "YB"    << ","
       << "It" << endl;
  fout.close();
}
void            Particle::log() {
  ofstream fout;
  fout.open("Log.csv",fstream::app);
    fout << m_id << ","  << m_fitness << ","         << m_pbest << ","
         << m_x.at(0)    << ","       << m_x.at(1)   << ","
         << m_v.at(0)    << ","       << m_v.at(1)   << ","
         << m_x_b.at(0)  << ","       << m_x_b.at(1) << endl;
  fout.close();
}
void            Particle::logC() {
  ofstream fout;
  fout.open("Convergence.csv",fstream::app);
    fout << m_id << ","  << m_fitness << ","         << m_pbest << ","
         << m_x.at(0)    << ","       << m_x.at(1)   << ","
         << m_v.at(0)    << ","       << m_v.at(1)   << ","
         << m_x_b.at(0)  << ","       << m_x_b.at(1) << ",";
  fout.close();
}
void            sortBuble(unsigned int& M_pop,vector<Particle>& pop){
  bool swapped;
  Particle Temp;
  for (size_t i = 0; i < M_pop-1; i++) {
    swapped = false;
    for (size_t j = 0; j < M_pop-i-1; j++) {
      if (pop.at(j) > pop.at(j+1)) {
        Temp = pop.at(j); pop.at(j) = pop.at(j+1); pop.at(j+1) = Temp;
        swapped = true;
      }
    }
    if (swapped == false)
      break;
  }
}
float           sphere(vector<float>& x){
  float Temp = 0;
  for (size_t i = 0; i < x.size(); i++)
    Temp += pow(x.at(i), 2.0);

  return Temp;
}
float           Beale(vector<float>& x){
  float  func =
  pow(1.5 - x.at(0) + x.at(0)*x.at(1),2.0) +
  pow(2.25 - x.at(0) + x.at(0)*(pow(x.at(1),2.0)),2.0) +
  pow(2.625 - x.at(0) + x.at(0)*(pow(x.at(1),3.0)),2.0);
  return func;
}
float           GoldPrice(vector<float>& x){
  float  func =
    (1 + (pow(x.at(0) + x.at(1) + 1,2.0))*(19 - 14*x.at(0) + 3*pow(x.at(0),2.0)
   - 14*x.at(1) + 6*x.at(0)*x.at(1) + 3*pow(x.at(1),2.0)))   *   (30 + (pow(2*x.at(0)
   - 3*x.at(1),2.0))*(18 - 32*x.at(0) + 12*pow(x.at(0),2.0) + 48*x.at(1)
   - 36*x.at(0)*x.at(1) + 27*pow(x.at(1),2.0)));
  return func;
}
float           McCormick(vector<float>& x){
  float  func = sin(x.at(0) + x.at(1)) + pow(x.at(0) - x.at(1),2.0) - 1.5*x.at(0) + 2.5*x.at(1) + 1;
  return func;
}
vector<float>   suma(const vector<float>& u,const vector<float>& v){
  vector <float> sum;
  if (u.size() < v.size()){
    for (unsigned int i = 0; i < u.size() ; i++)
      sum.push_back( u.at(i) + v.at(i) );
    for (size_t i = u.size(); i < v.size(); i++)
      sum.push_back(v.at(i));
  }
  else if (u.size() > v.size()){
    for (unsigned int i = 0; i < v.size() ; i++)
      sum.push_back( u.at(i) + v.at(i) );
    for (size_t i = v.size(); i < u.size(); i++)
      sum.push_back(u.at(i));
  }
  else {
  for (unsigned int i = 0; i < u.size() ; i++)
    sum.push_back( u.at(i) + v.at(i) );
  }
  return sum;
}
vector<float>   resta(const vector<float>& u,const vector<float>& v){
  vector <float> subs;
  if (u.size() < v.size()){
    for (unsigned int i = 0; i < u.size() ; i++)
      subs.push_back( u.at(i) - v.at(i) );
    for (size_t i = u.size(); i < v.size(); i++)
      subs.push_back(-v.at(i));
  }
  else if (u.size() > v.size()){
    for (unsigned int i = 0; i < v.size() ; i++)
      subs.push_back( u.at(i) - v.at(i) );
    for (size_t i = v.size(); i < u.size(); i++)
      subs.push_back(-u.at(i));
  }
  else {
  for (unsigned int i = 0; i < u.size() ; i++)
    subs.push_back( u.at(i) - v.at(i) );
  }
  return subs;
}
vector<float>   product(const vector<float>& u,const vector<float>& v){
  vector<float> prods;
  for (size_t i = 0; i < u.size(); i++) { prods.push_back(u.at(i) * v.at(i)); }
  return prods;
}
vector<float>   randGenMS(unsigned int& dimension, float& max, float& min,vector<float>& rand_coods){
  static thread_local mt19937_64   mt(hash<thread::id>()(this_thread::get_id()));
  uniform_real_distribution<float> dist(min,max);
  for (unsigned int i = 0; i < dimension; i++) { rand_coods.push_back(dist(mt)); }
  return rand_coods;
}

/*
simMov cleaning vectors awkward behaviour
switch, headers, input pointer
*/
