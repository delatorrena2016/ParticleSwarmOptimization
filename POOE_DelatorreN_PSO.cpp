#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

class Particle{
private:
  vector<float> m_v, m_x, m_x_b;
  float         m_fitness = 1;

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
  vector<float> const &getX() const { return m_x; }
  vector<float> const &getXB() const { return m_x_b; }
  float         getFit() { return m_fitness; }

  void          setFit(float fit) { m_fitness = fit; }
  void          setX(const vector<float>& _x) { m_x = _x; }
  void          setV(const vector<float>& _x) { m_v = _x; }
  void          setXB(const vector<float>& _x) { m_x_b = _x; }
  void          log();
  void          evFit();
};

void            simMov(vector<Particle>& x_b_g,float& phi_1,float& phi_2,vector<float>& rand_coods,
                vector<Particle>& pop,unsigned int& M_pop,unsigned int& dimension);
vector<float>   randGen(unsigned int& dimension,float& max, float& min,vector<float>& rand_coods);
float           sphere(vector<float>& x);
vector<float>   suma(vector<float>&,vector<float>&);
vector<float>   resta(vector<float>&,vector<float>&);
vector<float>   product(vector<float>& u, vector<float>& v);
void            sortBuble(unsigned int& M_pop,vector<Particle>& pop);
void            initSwarm(unsigned int& M_pop,unsigned int& dimension,
                float& min,float& max,vector<float>& rand_coods,vector<Particle>& pop);
float           zero = 0.0;

int main() {
  float            min = -10, max = 10;
  unsigned int     dimension = 2, M_pop = 10, N=2, count=0;
  float            phi_1 = 2.2, phi_2 = 2.0;
  vector<Particle> pop, x_b_g, U;
  vector<float>    rand_coods;

  initSwarm(M_pop,dimension,min,max,rand_coods,pop);
  ofstream         fout;
  while (count < N) {
    for (size_t i = 0; i < M_pop; i++) { pop.at(i).evFit(); }
    fout.open("Bitacora.csv",fstream::app);
    fout << " id " << "," << " Fitness " << ","
         << " X " << "," << " V " << "," << " Previous best " << endl;
    fout.close();
    for (size_t i = 0; i < M_pop; i++) { pop.at(i).log(); }
    sortBuble(M_pop,pop);
    fout.open("Bitacora.csv",fstream::app);
    fout << "..............................................\n";
    fout.close();
    for (size_t i = 0; i < M_pop; i++) { pop.at(i).log(); }


    count++;
  }

  /*for (size_t i = 0; i < 2; i++) {
    for (size_t j = a; j < b; j++) {
      pop.at(j).log();
    }  a = 97; b = 100;
  }*/

  return 0;
}
void        simMov(vector<Particle>& x_b_g,float& phi_1,float& phi_2,vector<float>& rand_coods,vector<Particle>& pop,unsigned int& M_pop,unsigned int& dimension) {
  x_b_g.push_back(pop.at(M_pop-1));
  vector<float> _xbg = x_b_g.at(0).getX();
  for (size_t i = 0; i < M_pop; i++) {
    vector<float> _x = pop.at(i).getX();
    vector<float> _v = ;
    vector<float> _xb = pop.at(i).getXB();
    vector<float> _v1 = product(resta(_xbg,_x),randGen(dimension,phi_2,zero,rand_coods));
    vector<float> _v2 = product(resta(_xb,_x),randGen(dimension,phi_1,zero,rand_coods));
    vector<float> _v3 = suma(_v1,_v2);
    vector<float> _nv = suma(_v,_v3);
    vector<float> _nx = suma(_x,_v);
    pop.at(i).setX(_x); pop.at(i).setV(_v); pop.at(i).setXB(_xb);

    pop.at(i).setV(suma(pop.at(i).getV(),));
  }
}
void          Particle::evFit() { m_fitness = sphere(m_x); }
void          initSwarm(unsigned int& M_pop,unsigned int& dimension,float& min,float& max,vector<float>& rand_coods,vector<Particle>& pop){
  srand((unsigned)time(NULL)); //dont run for every random
  for (size_t i = 0; i < M_pop; i++) {//
    randGen(dimension,max,min,rand_coods);
    pop.emplace_back( rand_coods );
    rand_coods.clear();
  }
}
void          Particle::log() {
  ofstream fout;
  fout.open("Bitacora.csv",fstream::app);
    fout << m_id << "," << m_fitness << ","
         << m_x.at(0)    << "," << m_x.at(1)  << ","
         << m_v.at(0)    << "," << m_v.at(1)  << ","
         << m_x_b.at(0)   << "," << m_x_b.at(1) << endl;
  fout.close();
}
void          sortBuble(unsigned int& M_pop,vector<Particle>& pop){
  bool swapped;
  for (size_t i = 0; i < M_pop-1; i++) {
    swapped = false;
    for (size_t j = 0; j < M_pop-i-1; j++) {
      if (pop.at(j) > pop.at(j+1)) {
        Particle Temp = pop.at(j); pop.at(j) = pop.at(j+1); pop.at(j+1) = Temp;
        swapped = true;
      }
    }
    if (swapped == false)
      break;
  }
}
float         sphere(vector<float>& x){
  float Temp = 0;
  for (size_t i = 0; i < x.size(); i++)
    Temp += pow(x.at(i), 2.0);

  return Temp;
}
vector<float> suma(vector<float>& u, vector<float>& v){
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
vector<float> resta(vector<float>& u, vector<float>& v){
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
vector<float> product(vector<float>& u, vector<float>& v){
  vector<float> prods;
  for (size_t i = 0; i < u.size(); i++) { prods.push_back(u.at(i) * v.at(i)); }
  return prods;
}
vector<float> randGen(unsigned int& dimension, float& max, float& min,vector<float>& rand_coods) {
  for (unsigned int j = 0; j < dimension; j++) {
    rand_coods.push_back( min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (max - min))));
  }
  return rand_coods;
}
